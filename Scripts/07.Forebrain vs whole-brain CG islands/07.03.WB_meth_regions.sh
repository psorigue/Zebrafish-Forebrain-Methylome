#!/bin/bash


# Script to calculate methylation levels in predefined genomic regions from methylation call files of zebrafish whole-brain methylome samples using bedtools.
# Bedtools v2.31.1

# Define regions
REGIONS_NAME="cgi"
REGION_FILE="$HOME/Pol/Methylome/methylation_regions/regions/${REGIONS_NAME}.bed"

# Define paths
PATH_POSITIONS="$HOME/Pol/Methylome/Chaterjee/meth_calls/output"
PATH_OUT="$HOME/Pol/Methylome/Chaterjee/methylation_regions/${REGIONS_NAME}"

# Define samples
SAMPLES=(M1 M2 F1 F2)

# Create output directory
mkdir -p "${PATH_OUT}"

for SAMPLE in "${SAMPLES[@]}" ; do

    POSITIONS="${PATH_POSITIONS}/${SAMPLE}_cpg.bed"
    OUT_NAME="${PATH_OUT}/${SAMPLE}_${REGIONS_NAME}.txt"
    
    # Output: chr, start (0-based), end (0-based), name of -a region (if applicable), mean methylation %, mean coverage, number of covered cpg positions in the region. Only positions with coverage >=5 are considered.
    bedtools map \
      -a <(sort -k1,1 -k2,2n "${REGION_FILE}") \
      -b <(awk '$4 >= 5' "${POSITIONS}") \
      -c 5,4,5 \
      -o mean,mean,count \
      > "${OUT_NAME}"
      
done
  
exit