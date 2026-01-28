#!/bin/bash

# This script evaluates the performance of BAM files by calculating read loss statistics based on various filtering criteria such as secondary/supplementary alignments, mapping quality, and read length.
# Samtools v 1.21

# Define paths
path_in=~/fil/Epi/Data_genomics/bams_genomics/Bam_Aligned_Cat/
qc_out=~/fil/Epi/Data_genomics/qc_stats/read_loss_summary.tsv

# Set number of threads
thr=2

# Create output directory if it doesn't exist
mkdir -p "$(dirname "${qc_out}")"

# Header
echo -e "sample\ttotal\tsecondary\tsupplementary\tMAPQ_lt10\tlen_lt200" > "${qc_out}"

# Qualtity Control filtering
for bam in "${path_in}"/*.bam; do
    sample=$(basename "${bam}" .bam)
    echo "Processing ${sample}"

    total=$(samtools view -@ "${thr}" -c "${bam}")
    secondary=$(samtools view -@ "${thr}" -c -f 256 "${bam}")
    supplementary=$(samtools view -@ "${thr}" -c -f 2048 "${bam}")

    # MAPQ < 10
    mapq_ge10=$(samtools view -@ "${thr}" -c -q 10 "${bam}")
    mapq_lt10=$(( total - mapq_ge10 ))

    # Read length < 200 bp (samtools ≥1.10)
    len_lt200=$(samtools view -@ "${thr}" -c -e 'length(seq) < 200' "${bam}")

    echo -e "${sample}\t${total}\t${secondary}\t${supplementary}\t${mapq_lt10}\t${len_lt200}" >> "${qc_out}"
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
