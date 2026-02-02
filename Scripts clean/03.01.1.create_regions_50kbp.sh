#!/bin/bash

# Create 50kb genomic regions for methylation level calculation. 
# Bedtools XXX

# Set paths
genome_index=~/fil//Epi/Ref_genome/GRCz12tu/genome/GCF_049306965.1_GRCz12tu_genomic.fna.fai
path_out=~/fil/Methylome/methylation_regions/regions/

mkdir -p $path_out

# Genome 50k bins excluding mitochondrial chromosome
bedtools makewindows -g <(cat $genome_index | cut -f1,2 | grep -v "NC_002333.2" ) -w 50000 > "$path_out"/genome_50kb_bins.bed # The output is 0-based


exit