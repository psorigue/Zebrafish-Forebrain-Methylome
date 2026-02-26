#!/bin/bash

# Create 50kb genomic regions for methylation level calculation. 
# Bedtools v2.31.1

# Set paths
GENOME_INDEX="$HOME/Ref_genome/GRCz12tu/GCF_049306965.1_GRCz12tu_genomic.fna.fai"
PATH_OUT="$HOME/methylation_regions/regions"

# Create output directory
mkdir -p "${PATH_OUT}"

# Genome 50k bins excluding mitochondrial chromosome
bedtools makewindows -g <(cat "${GENOME_INDEX}" | cut -f1,2 | grep -v "NC_002333.2" ) -w 50000 > "${PATH_OUT}/genome_50kb_bins.bed" # The output is 0-based


exit