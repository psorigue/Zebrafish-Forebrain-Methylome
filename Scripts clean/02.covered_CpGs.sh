#!/bin/bash

# Script to compute the number of CpG sites covered (depth >1) in each of the six zebrafish forebrain methylome samples.
# Bedtools v2.31.1

# Set variables
path_out=~/fil/Methylome/CpG_sites/coverage_CpG_samples/
cpg_bed=~/fil/Methylome/CpG_sites/CpG_genome/CpG_sites_genome.bed # File generated with 01.modkit_compute_CpGs.sh
path_bam=~/fil/Epi/Data_genomics/bams_filtered/ # Ready-to-use BAM files location

# Create output folders
mkdir -p "$path_out"
mkdir -p "$path_out/sites"
mkdir -p "$path_out/covered_sites"

cd "$path_out"

# Loop over samples
for sam in {1..6}; do 
    
    # Define BAM file name
    bam=barcode0"$sam"_aligned_filt.bam

    # Raw counts per CpG
    bedtools coverage -counts -a "$cpg_bed" -b "$path_bam/$bam" > sites/rep"$sam".counts

    # CpGs with coverage >1
    awk '$7>1 {print $1, $2, $3, $4, $5, $6}' sites/rep"$sam".counts > covered_sites/rep"$sam"_covered.bed

    # Simplified depth file for merging (2 columns: position and count)
    awk '{print $1":"$2"-"$3"\t"$7}' sites/rep"$sam".counts > sites/rep"$sam"_simplified.counts

done

# Merge all samples and compute "at least one" and "all samples"
paste sites/rep[1-6]_simplified.counts | \
awk '{ # awk function to compute number of CpGs covered in at least one and in all samples
    c=0 
    for(i=2;i<=12;i+=2) if($i>1) c++
    if(c>=1) any++
    if(c==6) all++
}
END{
    print "CpGs >1 in at least one sample:", any
    print "CpGs >1 in all samples:", all
}' > positions_covered_across_samples.txt


exit
