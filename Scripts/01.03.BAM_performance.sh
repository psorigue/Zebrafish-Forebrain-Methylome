#!/bin/bash

# This script evaluates the performance of BAM files by calculating read loss statistics based on various filtering criteria such as secondary/supplementary alignments, mapping quality, and read length.
# Samtools v 1.21

# Define paths
PATH_IN="$HOME/fil/Epi/Data_genomics/bams_genomics/Bam_Aligned_Cat"
QC_OUT="$HOME/fil/Epi/Data_genomics/qc_stats/read_loss_summary.tsv"

# Set number of threads
THREADS=2

# Create output directory if it doesn't exist
mkdir -p "$(dirname "${QC_OUT}")"

# Header
echo -e "sample\ttotal\tsecondary\tsupplementary\tMAPQ_lt10\tlen_lt200" > "${QC_OUT}"

# Quality Control filtering
for BAM_FILE in "${PATH_IN}"/*.bam; do

    # Define sample name
    SAMPLE_NAME=$(basename "${BAM_FILE}" .bam)
    echo "Processing ${SAMPLE_NAME}"

    # Calculate statistics
    TOTAL=$(samtools view -@ "${THREADS}" -c "${BAM_FILE}")
    SECONDARY=$(samtools view -@ "${THREADS}" -c -f 256 "${BAM_FILE}")
    SUPPLEMENTARY=$(samtools view -@ "${THREADS}" -c -f 2048 "${BAM_FILE}")

    # MAPQ < 10
    MAPQ_GE10=$(samtools view -@ "${THREADS}" -c -q "${MAPQ_ACCEPTED:-10}" "${BAM_FILE}")
    MAPQ_LT10=$(( TOTAL - MAPQ_GE10 ))

    # Read length < 200 bp (samtools ≥1.10)
    LEN_LT200=$(samtools view -@ "${THREADS}" -c -e "length(seq) < ${MIN_LEN:-200}" "${BAM_FILE}")

    # Append to output file
    echo -e "${SAMPLE_NAME}\t${TOTAL}\t${SECONDARY}\t${SUPPLEMENTARY}\t${MAPQ_LT10}\t${LEN_LT200}" >> "${QC_OUT}"
done

# Compute percentages
awk 'NR==1{
    print $0 "\tsecondary_pc\tsupplementary_pc\tMAPQ_lt10_pc\tlen_lt200_pc"
}
NR>1{
    printf "%s\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n",
    $1,$2,$3,$4,$5,$6,
    ($3/$2)*100,($4/$2)*100,($5/$2)*100,($6/$2)*100
}' "${qc_out}" > "${qc_out%.tsv}_pc.tsv"

exit
