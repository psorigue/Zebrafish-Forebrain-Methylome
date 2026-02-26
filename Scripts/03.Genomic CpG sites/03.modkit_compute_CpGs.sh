#!/bin/bash

# Script to compute the number of CpG sites in the zebrafish genome and in each of the six zebrafish forebrain methylome samples.
# Modkit v 0.5.0
# Bedtools v2.31.1

# Set variables
REF_GENOME="$HOME/Ref_genome/GRCz12tu/GCF_049306965.1_GRCz12tu_genomic.fna"
PATH_OUT_SAMPLES="$HOME/CpG_sites/coverage_CpG_samples"
PATH_BAM="$HOME/Data/bams_filtered"
PATH_OUT_GENOME="$HOME/CpG_sites"

CHROMOSOME_FILE="$HOME/Ref_genome/GRCz12tu/chr_array.txt" # Exclude MT chromosome
CPG_BED="$HOME/CpG_sites/CpG_genome/CpG_sites_genome.bed" # Genomic CpG sites

# Read chromosomes
CHROMOSOMES=( $( cat "${CHROMOSOME_FILE}" ) )

# 1. Count CpG sites in genome. Exclude MT chromosome.
modkit motif bed "${REF_GENOME}" "CG" 0 | grep -F '+' | grep -v "NC_002333.2" > "${PATH_OUT_GENOME}/CpG_sites_genome.bed"

# Count CpG sites per chromosome
for CHR in "${CHROMOSOMES[@]}" ; do

    # Read number of occurrences of each chromosome
    CHR_COUNT=$( grep "${CHR}" "${PATH_OUT_GENOME}/CpG_sites_genome.bed" | wc -l )
    # Output chromosome and count
    echo -e "${CHR}\t${CHR_COUNT}"
    
done > "${PATH_OUT_GENOME}/CpG_number_chromosomes.txt"


# 2. Samples CpG
# Create output folders
mkdir -p "${PATH_OUT_SAMPLES}"
mkdir -p "${PATH_OUT_SAMPLES}/sites"
mkdir -p "${PATH_OUT_SAMPLES}/covered_sites"

cd "${PATH_OUT_SAMPLES}"

# Loop over samples
for SAM in {1..6}; do 
    BAM_FILE="barcode0${SAM}_aligned_filt.bam"

    # Raw counts per CpG
    bedtools coverage -counts -a "${CPG_BED}" -b "${PATH_BAM}/${BAM_FILE}" > sites/rep"${SAM}".counts

    # CpGs with coverage >1
    awk '$7>1 {print $1, $2, $3, $4, $5, $6}' sites/rep"${SAM}".counts > covered_sites/rep"${SAM}"_covered.bed

    # Simplified depth file for merging
    awk '{print $1":"$2"-"$3"\t"$7}' sites/rep"${SAM}".counts > sites/rep"${SAM}"_simplified.counts
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

