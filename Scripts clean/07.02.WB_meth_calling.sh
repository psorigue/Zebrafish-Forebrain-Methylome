#!/bin/bash

# Script to call methylation states from mapped bisulfite sequencing reads from zebrafish whole-brain methylome samples using Bismark.
# Bismark v0.25.1

# Set paths
PATH_BAMS="$HOME/fil/Methylome/Chaterjee/mapping"
PATH_METH="$HOME/fil/Methylome/Chaterjee/meth_calls"
PATH_OUT="$HOME/fil/Methylome/Chaterjee/meth_calls/output"

# Define sample names
SAMPLES=(M1 M2 F1 F2)

# Set number of threads
THREADS=8

# Create output directory
mkdir -p "${PATH_METH}"
cd "${PATH_METH}"

for SAMPLE in "${SAMPLES[@]}" ; do

    mkdir -p "${SAMPLE}"
    BAM_FILE="${PATH_BAMS}/${SAMPLE}/${SAMPLE}_bismark_bt2.bam"
    
    # 1. Methylation call
    bismark_methylation_extractor \
      --gzip \
      -o "${PATH_METH}/${SAMPLE}" \
      --comprehensive \
      --parallel "${THREADS}" \
      --bedGraph \
      --zero_based \
      "${BAM_FILE}"

    # 2. Arrange output to obtain: chr, pos (0-based), coverage, methylation percentage
    # 2.1. Decompress coverage
    gunzip -k "${PATH_METH}/${SAMPLE}/${SAMPLE}_bismark_bt2.bismark.cov.gz"
    
    # 2.2. Remove mitochondrial chromosome, reduce and compute new columns
    grep -v "NC_002333.2" "${PATH_METH}/${SAMPLE}/${SAMPLE}_bismark_bt2.bismark.cov" | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $3, $5+$6, $4}' > "${PATH_OUT}/${SAMPLE}_cpg.bed"
    
    # 2.3. Remove useless file
    rm -f "${PATH_METH}/${SAMPLE}/${SAMPLE}_bismark_bt2.bismark.cov"
    
      
done


exit