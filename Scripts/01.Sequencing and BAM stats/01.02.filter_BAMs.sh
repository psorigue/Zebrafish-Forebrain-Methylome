#!/bin/bash

# This script filters BAM files based on mapping quality and read length, then sorts and indexes them. It also generates statistics on read counts before and after filtering, as well as the distribution of mapping qualities.
# Samtools v 1.21
# Modkit v 0.5.0

# Set input/output paths
PATH_IN="$HOME/Data/bams_genomics"
PATH_OUT="$HOME/Data/bams_filtered"
STATS_DIR="$HOME/QC/qc_stats"

# Define output files for statistics
COUNTS_OUT="${STATS_DIR}/read_counts.tsv"
MAPQ_OUT="${STATS_DIR}/mapq_distribution.tsv"

# Define parameters
THREADS=8 # Working threads
MAPQ_ACCEPTED=10 # Minimum mapping quality
MIN_LEN=200 # Minimum read length

# Define samples to process
SAMPLES=(01 02 03 04 05 06)

# Create output directories
mkdir -p "${PATH_OUT}" "${STATS_DIR}"

# Headers for output files
echo -e "sample\tbefore\tafter" > "${COUNTS_OUT}"
echo -e "sample\tmapq" > "${MAPQ_OUT}"

# Main loop over samples: filter, sort, index, and collect statistics
for SAM in "${SAMPLES[@]}"; do
    
    # Zero-pad the number to match file naming 
    SAM_PADDED=$(printf "%02d" "${SAM}")
    
    # Define input BAM file
    BAM_FILE="${PATH_IN}/barcode${SAM_PADDED}_aligned.bam"
    
    # Define sample name and output BAM file
    SAMPLE_NAME=$(basename "${BAM_FILE}" .bam)
    OUT_BAM="${PATH_OUT}/${SAMPLE_NAME}_filt.bam"

    # Count initial number of reads
    BEFORE=$(samtools view -c "${BAM_FILE}")

    # Collect MAPQ distribution (before filtering)
    samtools view "${BAM_FILE}" \
      | awk -v s="${SAMPLE_NAME}" '{print s "\t" $5}' \
      >> "${MAPQ_OUT}"

    # Filter, sort. '-F 2308' removes unmapped reads, secondary and supplementary aligned reads. '-q' keeps only reads with quality above ${MAPQ_ACCEPTED}. '-e' filters by read length.
    samtools view -@ "${THREADS}" -b \
        -q "${MAPQ_ACCEPTED}" \
        -F 2308 \
        -e "length(seq) >= ${MIN_LEN}" \
        "${BAM_FILE}" \
      | samtools sort -@ "${THREADS}" -T "${PATH_OUT}/${SAMPLE_NAME}_tmp" -o "${OUT_BAM}"
    
    # Index filtered BAM
    samtools index -@ "${THREADS}" "${OUT_BAM}"

    # Count final number of reads
    AFTER=$(samtools view -c "${OUT_BAM}")

    # Output read counts before and after filtering
    echo -e "${SAMPLE_NAME}\t${BEFORE}\t${AFTER}" >> "${COUNTS_OUT}"

done

exit
