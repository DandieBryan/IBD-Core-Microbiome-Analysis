#!/bin/bash -l
#SBATCH --time=06:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=8
#SBATCH --tmp=12G

url=$1
direct=$2
fastq=$3
output=$4

# Make output directory exists
mkdir -p "${output}/raw_reads"
cd "${output}/raw_reads"

lftp "$url" -e "cd \"$direct\"; get \"$fastq\"; bye"