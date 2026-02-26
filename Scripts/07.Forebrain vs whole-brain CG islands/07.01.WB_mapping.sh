#!/bin/bash

# Script to map bisulfite sequencing reads from zebrafish whole-brain methylome samples to the reference genome using Bismark.
# Publication: Chaterjee et al., 2014
# Bismark v0.25.1

# Set path to fastq files
PATH_RUNS="${HOME}/Chaterjee/bismark_data_process/fastq_runs"

# Define sample names
SAMPLES=(M1 M2 F1 F2) # Sample names match those in the fastq file names

# Set number of threads
THREADS=8

mkdir -p "${HOME}/Chaterjee/bismark_data_process/mapping"
cd "${HOME}/Chaterjee/bismark_data_process/mapping"

for SAMPLE in "${SAMPLES[@]}" ; do

    FASTQ_FILE="${PATH_RUNS}/${SAMPLE}.fastq"
    
    mkdir -p "${SAMPLE}"
    
    # Bismark mapping
    # N 1 increases sensitivity, it makes it slower
    bismark \
      --genome_folder ../genome/ \
      --fastq "${FASTQ_FILE}" \
      -N 1 \
      --parallel "${THREADS}" \
      -p "${THREADS}" \
      --output_dir "${SAMPLE}"
  
done 


