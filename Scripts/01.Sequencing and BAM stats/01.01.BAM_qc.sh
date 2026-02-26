#!/bin/bash

# This script runs Qualimap to perform quality control on BAM files aligned to the zebrafish genome. It also computes read length from non-filtered BAM files.
# Qualimap v.2.3

# Set paths
PATH_BAMS="$HOME/Data/bams_genomics"
REF_GTF="$HOME/Ref_genome/GRCz12tu/genomic.gtf"
PATH_OUT="$HOME/QC/bam_qc/qualimap"
PATH_QC="$HOME/QC/qc_stats"

# Create output directories
mkdir -p "${PATH_OUT}" ; mkdir -p "${PATH_QC}"
cd "${PATH_OUT}"

# Loop over BAM files
for SAMPLE in {01..06} ; do

    # Define BAM file
    BAM_FILE="${PATH_BAMS}/barcode${SAMPLE}_aligned.bam"

    # Run Qualimap
    qualimap bamqc -bam "${BAM_FILE}" -gff "${REF_GTF}" -nw 400 -hm 3 --java-mem-size=8G

    # Extract read length information from non-filtered BAM files. It uses the 'rl' tag in the BAM file (ONT-specific).
    samtools view "${BAM_FILE}" \
    | awk '
        {
        for (i=12; i<=NF; i++) {
            if ($i ~ /^rl:i:/) {
            split($i,a,":")
            print a[3]
            }
        }
        }' > "${PATH_QC}/read_length_rep${SAMPLE}.txt"

done
exit