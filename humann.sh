#!/bin/bash -l
#SBATCH --time=05:00:00
#SBATCH --mem=64G
#SBATCH --ntasks=8
#SBATCH --tmp=64G

### THIS SCRIPT GETS CALLED BY THE "QC_Preprocessing_and_Functional_Annotation.sh" FILE ###

# Input arguments
SAMPLE=${1}
BASENAME=${2}
OUTPUT_DIR=${3}

SITE_PACKAGES="INSERT PATHWAY TO SITE PACKAGES DIRECTORY"

# Generate MetaPhlAn taxonomic profile 
mkdir -p "${OUTPUT_DIR}/MetaPhlAn"

metaphlan "$SAMPLE" \
    --bowtie2db "${SITE_PACKAGES}/metaphlan" \
    --index mpa_vOct22_CHOCOPhlAnSGB_202403 \
    -t rel_ab_w_read_stats \
    -o "${OUTPUT_DIR}/MetaPhlAn/${BASENAME}_profile.tsv" \
    --input_type fasta

# Run HUMAnN using the MetaPhlAn taxonomic profile
mkdir -p "$OUTPUT_DIR/HUMAaN/Samples/$BASENAME"

humann --input "$SAMPLE" \
    --output "$OUTPUT_DIR/HUMAaN/Samples/$BASENAME" \
    --taxonomic-profile "${OUTPUT_DIR}/MetaPhlAn/${BASENAME}_profile.tsv" \
    --input-format fasta
