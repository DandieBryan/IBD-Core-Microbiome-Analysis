#!/bin/bash -l 
#SBATCH --time=48:00:00
#SBATCH --mem=24G                
#SBATCH --cpus-per-task=1        
#SBATCH --nodes=1                
#SBATCH --ntasks=1               

module load kraken

R1=$1
R2=$2
SAMPLE=$3
output=$4
dependencies=$5

# Create output directories
mkdir -p "${output}/host_removed/clean"
mkdir -p "${output}/host_removed/host"
mkdir -p "${output}/host_removed/report"


# Run Kraken2
kraken2 --db "${dependencies}/minikraken2_v1_8GB" \
    --paired "${R1}" "${R2}" \
    --threads 8 \
    --unclassified-out "${output}/host_removed/clean/${SAMPLE}_R#.fastq.gz" \
    --classified-out "${output}/host_removed/host/${SAMPLE}_R#.fastq.gz" \
    --output "${output}/host_removed/report/${SAMPLE}_kraken_report.txt"