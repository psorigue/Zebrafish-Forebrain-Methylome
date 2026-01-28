#!/bin/bash

# Script to compute the number of CpG sites in the zebrafish genome and in each of the six zebrafish forebrain methylome samples.
# Modkit v 0.5.0
# Bedtools v2.31.1

# Set variables
ref_genome=~/fil/Epi/Ref_genome/GRCz12tu/genome/GCF_049306965.1_GRCz12tu_genomic.fna
path_out_samples=~/fil/Methylome/CpG_sites/coverage_CpG_samples/
path_bam=~/fil/Epi/Data_genomics/bams_filtered/
path_out_genome=~/fil/Methylome/CpG_sites/

chromosome_file=~/fil/Epi/Ref_genome/GRCz12tu/genome/chr_array.txt # Exclude MT chromosome
cpg_bed=~/fil/Methylome/CpG_sites/CpG_genome/CpG_sites_genome.bed # Genomic CpG sites

# Read chromosomes
chromosomes=( $( cat $chromosome_file ) )

# 1. Count CpG sites in genome. Exclude MT chromosome.
modkit motif bed "$ref_genome" "CG" 0 | grep -F '+' | grep -v "NC_002333.2" > "${path_out_genome}/CpG_sites_genome.bed"

# Count CpG sites per chromosome
for chr in "${chromosomes[@]}" ; do

    # Read number of occurrences of each chromosome
    chr_count=$( grep "$chr" "${path_out_genome}/CpG_sites_genome.bed" | wc -l )
    # Output chromosome and count
    echo -e "${chr}\t${chr_count}"
    
done > "${path_out_genome}/CpG_number_chromosomes.txt"


# 2. Samples CpG
# Create output folders
mkdir -p "$path_out_samples"
mkdir -p "$path_out_samples/sites"
mkdir -p "$path_out_samples/covered_sites"

cd "$path_out_samples"

# Loop over samples
for sam in {1..6}; do 
    bam=barcode0"$sam"_aligned_filt.bam

    # Raw counts per CpG
    bedtools coverage -counts -a "$cpg_bed" -b "$path_bam/$bam" > sites/rep"$sam".counts

    # CpGs with coverage >1
    awk '$7>1 {print $1, $2, $3, $4, $5, $6}' sites/rep"$sam".counts > covered_sites/rep"$sam"_covered.bed

    # Simplified depth file for merging
    awk '{print $1":"$2"-"$3"\t"$7}' sites/rep"$sam".counts > sites/rep"$sam"_simplified.counts
done

# Merge all samples and compute "at least one" and "all samples"
paste sites/rep[1-6]_simplified.counts | \
awk '{
    c=0
    for(i=2;i<=12;i+=2) if($i>1) c++
    if(c>=1) any++
    if(c==6 all++
}
END{
    print "CpGs >1 in at least one sample:", any
    print "CpGs >1 in all samples:", all
}' > positions_covered_across_samples.txt


exit

