#!/bin/bash

# This script runs Qualimap to perform quality control on BAM files aligned to the zebrafish genome.
# Qualimap v.2.3

# Set paths
path_bams=~/fil/Epi/Data_genomics/bams_genomics/Bam_Aligned_Cat/
ref_gtf=~/fil/Epi/Ref_genome/GRCz12tu/genome/genomic.gtf
path_out=~/fil/Methylome/QC/bam_qc/qualimap/

cd $path_out

# Loop over BAM files
for num in {01..06} ; do

    # Define BAM file
    file="$path_bams"/barcode"$num"_aligned.bam

    # Run Qualimap
    qualimap bamqc -bam $file -gff $ref_gtf -nw 400 -hm 3 --java-mem-size=8G 

done
exit