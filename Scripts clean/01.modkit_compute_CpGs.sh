#!/bin/bash

# Script to compute the number of CpG sites in the zebrafish GRCz12tu genome and per chromosome.
# Modkit v 0.5.0

# Set variables
ref_genome=~/fil/Epi/Ref_genome/GRCz12tu/genome/GCF_049306965.1_GRCz12tu_genomic.fna
path_out=~/fil/Methylome/CpG_sites/
chromosome_file=~/fil/Epi/Ref_genome/GRCz12tu/genome/chr_array.txt # Excludes MT chromosome

# Define motif for modkit
motif="CG"

# Read chromosomes
chromosomes=( $( cat $chromosome_file ) )

# Count CpG sites in genome. Exclude MT chromosome.
modkit motif bed "$ref_genome" "$motif" 0 | grep -F '+' | grep -v "NC_002333.2" > "${path_out}/CpG_sites_genome.bed"

# Count CpG sites per chromosome
for chr in "${chromosomes[@]}" ; do

    # Read number of occurrences of each chromosome
    chr_count=$( grep "$chr" "${path_out}/CpG_sites_genome.bed" | wc -l )
    # Output chromosome and count
    echo -e "${chr}\t${chr_count}"
    
done > "${path_out}/CpG_number_chromosomes.txt"


exit