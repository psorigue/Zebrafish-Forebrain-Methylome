#!/bin/bash

# Script to calculate methylation levels in predefined genomic regions from methylation call files of zebrafish forebrain methylome samples using bedtools. The coverage has been previously filtered to include only positions with coverage >=5 reads.
# Bedtools v2.31.1

# Define regions
REGIONS_NAME="cgi"
REGION_FILE="$HOME/Pol/Methylome/methylation_regions/regions/${REGIONS_NAME}.bed"

# Define paths
PATH_POSITIONS="$HOME/Pol/Methylome/Chaterjee/forebrain_datasets/meth_proportions_cpg"
PATH_OUT="$HOME/Pol/Methylome/Chaterjee/forebrain_datasets/methylation_regions/${REGIONS_NAME}"

# Create output directory
mkdir -p "${PATH_OUT}"

# Define samples
SAMPLES=(01 02 03 04 05 06)

for SAMPLE in "${SAMPLES[@]}" ; do

    POSITIONS="${PATH_POSITIONS}/${SAMPLE}_cpg.bed"
    OUT_NAME="${PATH_OUT}/${SAMPLE}_${REGIONS_NAME}.txt"
    
    # Output: chr, start (0-based), end (0-based), name of -a region (if applicable), mean methylation %, mean coverage, number of covered cpg positions in the region
    bedtools map \
      -a <(sort -k1,1 -k2,2n "${REGION_FILE}") \
      -b "${POSITIONS}" \
      -c 5,4,5 \
      -o mean,mean,count \
      > "${OUT_NAME}"
      
done

exit