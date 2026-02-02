#!/bin/bash

# This script runs Qualimap to perform quality control on BAM files aligned to the zebrafish genome. It also computes read length from non-filtered BAM files.
# Qualimap v.2.3

# Set paths
path_bams=~/fil/Epi/Data_genomics/bams_genomics/Bam_Aligned_Cat/
ref_gtf=~/fil/Epi/Ref_genome/GRCz12tu/genome/genomic.gtf
path_out=~/fil/Methylome/QC/bam_qc/qualimap/
path_qc=~/fil/Methylome/QC/qc_stats/


mkdir -p $path_out ; mkdir -p $path_qc
cd $path_out

# Loop over BAM files
for num in {01..06} ; do

    # Define BAM file
    file="$path_bams"/barcode"$num"_aligned.bam

    # Run Qualimap
    qualimap bamqc -bam $file -gff $ref_gtf -nw 400 -hm 3 --java-mem-size=8G 

    # Extract read length information from non-filtered BAM files
    samtools view $file \
    | awk '
        {
        for (i=12; i<=NF; i++) {
            if ($i ~ /^rl:i:/) {
            split($i,a,":")
            print a[3]
            }
        }
        }' > "$path_qc"/read_length_rep"$num".txt  

done
exit