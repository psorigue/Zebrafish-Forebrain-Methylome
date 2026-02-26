#!/bin/bash


# Script to calculate methylation levels in predefined genomic regions from methylation call files of zebrafish whole-brain methylome samples using bedtools.
# Bedtools v2.31.1

# Define regions
REGION_FILE="$HOME/Data_methylation/methylation_regions/regions/cgi.bed"

# Define paths
PATH_POSITIONS="$HOME/Chaterjee/bismark_data_process/meth_calls/output"
PATH_OUT="$HOME/Chaterjee/methylation_cgi/whole-brain"

# Define samples
SAMPLES=(M1 M2 F1 F2)

# Create output directory
mkdir -p "${PATH_OUT}"

for SAMPLE in "${SAMPLES[@]}" ; do

    POSITIONS="${PATH_POSITIONS}/${SAMPLE}_cpg.bed"
    OUT_NAME="${PATH_OUT}/${SAMPLE}_cgi.txt"
    
    # Output: chr, start (0-based), end (0-based), name of -a region (if applicable), mean methylation %, mean coverage, number of covered cpg positions in the region. Only positions with coverage >=5 are considered.
    bedtools map \
      -a <(sort -k1,1 -k2,2n "${REGION_FILE}") \
      -b <(awk '$4 >= 5' "${POSITIONS}") \
      -c 5,4,5 \
      -o mean,mean,count \
      > "${OUT_NAME}"
      
done
  
exit