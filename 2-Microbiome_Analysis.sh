#!/bin/bash -l 
#SBATCH --time=48:00:00
#SBATCH --mem=24G                
#SBATCH --cpus-per-task=1        
#SBATCH --nodes=1                
#SBATCH --ntasks=1  



##########################################################################################
### MADE BY BRYAN LE, NOVEMBER 2024

# The second half of this pipeline uses Qiime2 to perform rarefaction, OTU filtering, 
# alpha diversity analysis, and beta diversity analysis. Additionally, it uses MetaPhlAn 
# and HUMAnN for taxonomic and functional profiling. I also included additional steps for 
# data reformatting and cleanup for downstream analyses.

# CHANGE THE RAREFACTION DEPTH IN LINE 47 BEFORE RUNNING (currently set as -d 0000 for placeholder)
# CHANGE OTU THRESHOLD IN LINE 51 BEFORE RUNNING (currently set as -s 0000 for placeholder)

# The taxon level is currently set as 'L7' (species).  To perform analysis at a different,
# find all instances of 'L7' and replace it with the desired level.

# As a side note, the pipeline had to be split into two parts because the second one
# requires you to input a rarefaction depth, but that can only be done after the first 
# part is complete.



##########################################################################################
### PART 0: DEFINING PATHWAYS AND DIRECTORIES, AND LOADING MODULES

DEPENENDENCIES="INSERT PATHWAY TO DEPENDENCIES DIRECTORY HERE"
OUTPUT_DIR="INSERT OUTPUT DIRECTORY HERE"

module load qiime/1.9.1_centos7
module load kraken
module load bowtie2
module load python/3.6.3



##########################################################################################
### PART 5: ALPHA AND BETA DIVERSITY ANALYSIS
# Rarefy L7 data
mkdir -p "${OUTPUT_DIR}/Kraken/Rarefaction"
single_rarefaction.py -i "${OUTPUT_DIR}/Kraken/taxon_tables/taxa_table_L7.biom" -d 0000 -o "${OUTPUT_DIR}/Kraken/Rarefaction/taxa_table_L7_rarefied.biom" # Replace 1000 with your rarefaction depth

# Filter out low OTUs 
mkdir -p "${OUTPUT_DIR}/Kraken/Taxa"
filter_otus_from_otu_table.py -i "${OUTPUT_DIR}/Kraken/Rarefaction/taxa_table_L7_rarefied.biom" -o "${OUTPUT_DIR}/Kraken/Taxa/taxa_table_L7_final.biom" -s 0000

# Convert taxon biom file back to txt
biom convert -i "${OUTPUT_DIR}/Kraken/Taxa/taxa_table_L7_final.biom" -o "${OUTPUT_DIR}/Kraken/Taxa/taxa_table_L7_final.txt" --to-tsv --table-type "Taxon table"

# Calculate Shannon and Chao1 alpha diversity; Bray-Curtis and Jaccard beta diversity
mkdir -p "${OUTPUT_DIR}/Kraken/alpha"
mkdir -p "${OUTPUT_DIR}/Kraken/beta/L7_beta"
alpha_diversity.py -m "chao1,observed_otus,shannon" -i "${OUTPUT_DIR}/Kraken/Taxa/taxa_table_L7_final.biom" -o "${OUTPUT_DIR}/Kraken/alpha/L7_alpha-diversity.txt"
beta_diversity.py -i "${OUTPUT_DIR}/Kraken/Taxa/taxa_table_L7_final.biom" -o "${OUTPUT_DIR}/Kraken/beta/L7_beta" -m "bray_curtis,binary_jaccard"

# Run principal coordinate analysis on Bray-Curtis values
principal_coordinates.py -i "${OUTPUT_DIR}/Kraken/beta/L7_beta/bray_curtis_taxa_table_L7_final.txt" -o "${OUTPUT_DIR}/Kraken/beta/L7_beta/bray_curtis_taxa_table_L7_final_pc.txt"



##########################################################################################
### PART 6: FUNCTIONAL ANALYSIS

# Extract Samples names remaining after rarefaction
IFS=$'\t' read -r -a BASENAME < <(head -n 1 "taxa_table_L7.txt")
job_ids=()

# Submit job to run metaphlan and humann on each sample
for i in "${BASENAME[@]}"; do
    FILE="${OUTPUT_DIR}/SHI7/combined_Files/$i*"
    job_id=$(sbatch humann.sh "$FILE" "$i" "$OUTPUT_DIR" | awk '{print $NF}')
    job_ids+=("$job_id")
done

# Wait for all functional annotations to complete
while true; do
    remaining_jobs=$(squeue -u $USER | awk 'NR>1 {print $1}' | grep -Ff <(printf "%s\n" "${job_ids[@]}") | wc -l)
    if [[ $remaining_jobs -eq 0 ]]; then
        break
    fi
    sleep 1800
done

# Copy all annotation files to separate directory
mkdir -p "$OUTPUT_DIR/HUMAaN/All_Files"

for i in "$OUTPUT_DIR/HUMAaN/Samples"; do
  cp $i/* "$OUTPUT_DIR/HUMAaN/All_Files"

done

# Merge all samples
mkdir -p "$OUTPUT_DIR/HUMAaN/Merged_Files"
humann_join_tables --input "$OUTPUT_DIR/HUMAaN/All_Files" --output "$OUTPUT_DIR/HUMAaN/Merged_Files/merged_genefamilies.tsv" --file_name genefamilies
humann_join_tables --input "$OUTPUT_DIR/HUMAaN/All_Files" --output "$OUTPUT_DIR/HUMAaN/Merged_Files/merged_pathabundance.tsv" --file_name pathabundance
humann_join_tables --input "$OUTPUT_DIR/HUMAaN/All_Files"--output "$OUTPUT_DIR/HUMAaN/Merged_Files/merged_reactions.tsv" --file_name reactions

# Normalize output
humann_renorm_table --input "$OUTPUT_DIR/HUMAaN/Merged_Files/merged_genefamilies.tsv" --output "$OUTPUT_DIR/HUMAaN/Merged_Files/normalized_genefamilies.tsv" --units cpm
humann_renorm_table --input "$OUTPUT_DIR/HUMAaN/Merged_Files/merged_pathabundance.tsv" --output "$OUTPUT_DIR/HUMAaN/Merged_Files/normalized_pathabundance.tsv" --units cpm
humann_renorm_table --input "$OUTPUT_DIR/HUMAaN/Merged_Files/merged_reactions.tsv" --output "$OUTPUT_DIR/HUMAaN/Merged_Files/normalized_reactions.tsv" --units cpm

# Stratify Data
humann_split_stratified_table --input "$OUTPUT_DIR/HUMAaN/Merged_Files/normalized_genefamilies.tsv" --output "$OUTPUT_DIR/HUMAaN/Merged_Files/stratified_genefamilies.tsv"
humann_split_stratified_table --input "$OUTPUT_DIR/HUMAaN/Merged_Files/merged_pathabundance.tsv" --output "$OUTPUT_DIR/HUMAaN/Merged_Files/stratified_pathabundance"
humann_split_stratified_table --input "$OUTPUT_DIR/HUMAaN/Merged_Files/merged_reactions.tsv" --output "$OUTPUT_DIR/HUMAaN/Merged_Files/stratified_reactions"



##########################################################################################
### PART 7: REFORMAT TAXA FILE

# Read the header row and write it to the output file
header=$(head -n 1 "${OUTPUT_DIR}/Kraken/Taxa/taxa_table_L7_final.txt")
echo -e "$header" > "${OUTPUT_DIR}/Kraken/Taxa/taxa_table_L7_final_relative_abundance.txt"

# Calculate the relative abundance of all OTUs for each sample
awk '
NR == 1 { 
    # Read the header and store sample names
    for (i = 2; i <= NF; i++) col_sum[i] = 0
    next 
}
{
    # For each row, calculate column-wise sums
    for (i = 2; i <= NF; i++) {
        col_sum[i] += $i
    }
    rows[NR] = $0  # Store the row for later use
}
END {
    # Output relative abundances
    for (r = 2; r <= NR; r++) {
        split(rows[r], row_values)
        printf "%s", row_values[1]  # Print OTU name
        for (i = 2; i <= NF; i++) {
            printf "\t%.18f", row_values[i] / col_sum[i]
        }
        print ""
    }
}' "${OUTPUT_DIR}/Kraken/Taxa/taxa_table_L7_final.txt" >> "${OUTPUT_DIR}/Kraken/Taxa/taxa_table_L7_final_relative_abundance.txt"

# Transpose relative abundance data
awk 'BEGIN {FS="\t"; OFS="\t"} {
    for (i=1; i<=NF; i++) {
        a[i,NR]=$i
    }
    if (NF>max_fields) max_fields=NF
} 
END {
    for (i=1; i<=max_fields; i++) {
        for (j=1; j<=NR; j++) {
            printf "%s%s", a[i,j], (j==NR ? "\n" : OFS)
        }
    }
}' "${OUTPUT_DIR}/Kraken/Taxa/taxa_table_L7_final_relative_abundance.txt" > "${OUTPUT_DIR}/Kraken/Taxa/taxa_table_L7_final_transposed_relative_abundance.txt"


# Transpose absolute abundance data
awk 'BEGIN {FS="\t"; OFS="\t"} {
    for (i=1; i<=NF; i++) {
        a[i,NR]=$i
    }
    if (NF>max_fields) max_fields=NF
} 
END {
    for (i=1; i<=max_fields; i++) {
        for (j=1; j<=NR; j++) {
            printf "%s%s", a[i,j], (j==NR ? "\n" : OFS)
        }
    }
}' "${OUTPUT_DIR}/Kraken/Taxa/taxa_table_L7_final.txt" > "${OUTPUT_DIR}/Kraken/Taxa/taxa_table_L7_final_transposed_absolute_abundance.txt"