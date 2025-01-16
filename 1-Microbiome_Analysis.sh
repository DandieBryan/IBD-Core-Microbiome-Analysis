#!/bin/bash -l 
#SBATCH --time=48:00:00
#SBATCH --mem=24G                
#SBATCH --cpus-per-task=1        
#SBATCH --nodes=1                
#SBATCH --ntasks=1  



##########################################################################################
### MADE BY BRYAN LE, NOVEMBER 2024

# This script is the first half of a pipeline I designed to download and analyze shotgun 
# metagenomic data. It retrieves data from 10 healthy individuals and 10 individuals with 
# inflammatory bowel disease (IBD), sourced from the European Nucleotide Archive under 
# study accession PRJEB7949. The pipeline processes the data by removing host reads with 
# Kraken2, performing quality control using SHI7, and initiating taxonomic analysis with 
# Kraken2.

# The taxon level is currently set as 'L7' (species).  To perform analysis at a different,
# find all instances of 'L7' and replace it with the desired level.

# As a side note, the pipeline had to be split into two parts because the second one
# requires you to input a rarefaction depth, but that can only be done after the first 
# part is complete.



##########################################################################################
### PART 0: DEFINING PATHWAYS AND DIRECTORIES
DEPENENDENCIES="INSERT PATHWAY TO DEPENDENCIES DIRECTORY HERE"
OUTPUT_DIR="INSERT OUTPUT DIRECTORY HERE"



##########################################################################################
### PART 1: DOWNLOADING METAGENOMIC SHOTGUN WGS FILES
job_ids_1=()

# Loop through all urls, submit batch job to download each FASTQ file, and collect job IDs
for i in {1..80}; do
    line=$(sed -n "${i}p" "${DEPENDENCIES}/url.txt")
    read url directory fastq <<< "$line"
    job_id=$(sbatch "${DEPENDENCIES}/download_files.sh" "$url" "$directory" "$fastq" "$OUTPUT_DIR" | awk '{print $NF}')
    job_ids_1+=("$job_id")
done

# Wait for downloads to finish
while true; do
    remaining_jobs=$(squeue -u $USER | awk 'NR>1 {print $1}' | grep -Ff <(printf "%s\n" "${job_ids_1[@]}") | wc -l)
    if [[ $remaining_jobs -eq 0 ]]; then
        break
    fi
    sleep 300
done



##########################################################################################
### PART 2: REMOVING HOST READS ###
job_ids_2=()

# Loop through R1 and R2 pairs for each sample and submit batch job to remove host reads
for R1 in ${OUTPUT_DIR}/raw_reads/*_R1.fastq.gz; do
    SAMPLE=$(basename $R1 _R1.fastq.gz)
    R2="${OUTPUT_DIR}/raw_reads/${SAMPLE}_R2.fastq.gz"
    job_id=$(sbatch ${DEPENDENCIES}/remove_host_reads.sh $R1 $R2 $SAMPLE $OUTPUT_DIR $DEPENDENCIES | awk '{print $NF}')
    job_ids_2+=($job_id)
done

# Wait for all jobs in the second batch to finish
while true; do
    remaining_jobs=$(squeue -u $USER | awk 'NR>1 {print $1}' | grep -Ff <(printf "%s\n" "${job_ids_2[@]}") | wc -l)
    if [[ $remaining_jobs -eq 0 ]]; then
        break
    fi
    sleep 300
done



##########################################################################################
### PART 3: QUALITY CONTROL AND PREPROCESSING ###

# Load modules
module load qiime/1.9.1_centos7
module load kraken
module load bowtie2
module load python/3.6.3

# Quality control for metagenomic shotgun sequencing data using default settings on Shi7
mkdir -p "${OUTPUT_DIR}/SHI7/Files"
python3 "${DEPENDENCIES}/shi7/shi7.py" -i "${OUTPUT_DIR}/host_removed/clean" -o "${OUTPUT_DIR}/SHI7/Files" --combine_fasta False

# Removing extra ".fna" in trimmed shotgun sequencing files
for f in "${OUTPUT_DIR}/SHI7/Files"/*.fa.fna; do
    mv "$f" "${OUTPUT_DIR}/SHI7/Files/$(basename "$f" .fa.fna).fna"
done

# Each sample was sequenced twice, lane 1 and lane 2, to increase accuracy.  We need to combine these FASTQ files for each sample for downstream analysis
mkdir -p "${OUTPUT_DIR}/SHI7/combined_Files"

for FILE1 in "${OUTPUT_DIR}/SHI7/Files"/*L001.fna; do
  BASENAME=$(basename "$FILE1" L001.fna)
  FILE2="${OUTPUT_DIR}/SHI7/Files/${BASENAME}L002.fna"

  cat "$FILE1" "$FILE2" > "${OUTPUT_DIR}/SHI7/combined_Files/${BASENAME}fna"
done



##########################################################################################
### PART 4: TAXONOMIC ANNOTATION

# Using Kraken2 for taxonomic classification
mkdir -p "${OUTPUT_DIR}/Kraken/Samples"
for f in "${OUTPUT_DIR}/SHI7/combined_Files"/*.fna; do
    echo "Processing $f"
    time kraken2 --db "${DEPENDENCIES}/minikraken2_v1_8GB" \
        --use-mpa-style \
        --output tmp \
        --report "${OUTPUT_DIR}/Kraken/Samples/$(basename $f .fna).txt" \
        --use-names $f
done

# Creating taxonomy table
mkdir -p  "${OUTPUT_DIR}/Kraken/taxon_tables"
wget https://raw.githubusercontent.com/danknights/mice5035/master/scripts/kraken2table.py -O kraken2table.py
python3 kraken2table.py "${OUTPUT_DIR}/Kraken/Samples"/*.txt "${OUTPUT_DIR}/Kraken/taxon_tables"

# Convert taxonomy table to biom format
for f in "${OUTPUT_DIR}/Kraken/taxon_tables"/*.txt; do
    echo "Converting $f to biom format"
    biom convert -i "$f" \
        --to-json \
        -o "$(dirname $f)/$(basename $f .txt).biom" \
        --process-obs-metadata taxonomy
done

# Make Stats Table
mkdir -p "${OUTPUT_DIR}/Kraken/Stats"
biom summarize-table \
    -i "${OUTPUT_DIR}/Kraken/taxon_tables/taxa_table_L7.biom" \
    -o "${OUTPUT_DIR}/Kraken/Stats/L7_stats.txt"
