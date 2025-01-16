# IBD-Core-Microbiome-Analysis

This repository contains scripts and documentation for a pipeline designed to download, process, and analyze shotgun metagenomic data. The project focuses on investigating microbiome variations between healthy individuals and those with inflammatory bowel disease (IBD). The data used in this analysis are sourced from the European Nucleotide Archive under study accession [PRJEB7949].

## Overview

The pipeline is divided into three main components:

1. **Data Retrieval and Preprocessing**
   - Retrieves shotgun metagenomic data for 10 healthy individuals and 10 individuals with IBD.
   - Processes the data by:
     - **Host Read Removal**: Using Kraken2 to filter out host genomic reads.
     - **Quality Control**: Ensuring high-quality data using SHI7.
     - **Initial Taxonomic Analysis**: Performing taxonomic classification using Kraken2.

2. **Diversity Analysis, Taxonomic Analysis, and Functional Profiling**
   - Employs Qiime2 for:
     - Rarefaction and OTU filtering.
     - Alpha and beta diversity analyses to study species richness and community structure.
   - Utilizes MetaPhlAn and HUMAnN for:
     - Taxonomic profiling.
     - Functional profiling.
   - Includes steps for reformatting and cleaning data for seamless downstream analyses.

3. **Statistical Analysis**
   - An R script builds on the output of the preprocessing and analysis pipeline to:
     - Examine taxonomic differences between groups.
     - Compare alpha and beta diversity metrics.
     - Analyze functional differences using statistical methods.
   - The goal is to identify and characterize variations in the core microbiomes of healthy individuals and those with IBD, offering insights into microbiological distinctions associated with IBD.

## Prerequisites

### Software and Tools

- **Kraken2**: For host read removal and initial taxonomic analysis.
- **SHI7**: For quality control.
- **Qiime2**: For diversity and OTU analysis.
- **MetaPhlAn** and **HUMAnN**: For taxonomic and functional profiling.
- **R**: For statistical analysis (along with required libraries).

### Data

The pipeline uses publicly available shotgun metagenomic data from the European Nucleotide Archive under accession [PRJEB7949].
