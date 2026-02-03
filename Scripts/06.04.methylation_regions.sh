#!/bin/bash

# Script to calculate mean methylation levels in defined genomic regions (50kb bins, gene_bodies, promoters, CGIs) for CpG (5mC, 5hmC), and CH (5mC) contexts. Bedtools map function is used to compute mean methylation levels in each region, mean coverage and number of covered sites for that particular motif. Coverage filtering (> or =5) is already applied in previous scripts.
# Bedtools v2.31.1

# Define regions
REGIONS_NAME="genome_50kb_bins" # genome_50kb_bins, genes, promoters, cgi
REGION_FILE="$HOME/fil/Methylome/methylation_regions/regions/${REGIONS_NAME}.bed"

# Define paths
PATH_POSITIONS="$HOME/fil/Methylome/Data_methylation/datasets_proportions"
PATH_OUT="$HOME/fil/Methylome/methylation_regions/output/${REGIONS_NAME}"

# Create output directory
mkdir -p "${PATH_OUT}"

# For 5mC and 5hmC in CpG context
for mod in 5mC 5hmC ; do 

    mkdir -p "${PATH_OUT}/${mod}"
    
    for SAM in {1..6} ; do
    
        POSITIONS="${PATH_POSITIONS}/${mod}/${mod}_sample0${SAM}_mpct.txt"
        OUT_NAME="${PATH_OUT}/${mod}/meth_mean_${mod}_rep${SAM}_${REGIONS_NAME}.txt"

        # Output: chr, start (0-based), end (0-based), name of -a region (if applicable), mean methylation %, mean coverage, number of covered cpg positions in the region
        bedtools map \
          -a <(sort -k1,1 -k2,2n "${REGION_FILE}") \
          -b "${POSITIONS}" \
          -c 4,5,4 \
          -o mean,mean,count \
          > "${OUT_NAME}"
          
    done
      
done

# For 5mC in CH contexts
for motif in CA CC CT ; do

    for sam in {1..6} ; do
      
        mkdir -p "$path_out"/ch
    
        for st in pos neg ; do
    
            positions="$path_positions"/ch/ch_sample0"$sam"_"$motif"_"$st"_mpct.txt
            out_name="$path_out"/ch/meth_mean_"$motif"_"$st"_rep"$sam"_"$regions_name".txt

            # Output: chr, start (0-based), end (0-based), name of -a region (if applicable), mean methylation %, mean coverage, number of covered cpg positions in the region
            bedtools map \
              -a <(sort -k1,1 -k2,2n "$region_file") \
              -b "$positions" \
              -c 4,5,4 \
              -o mean,mean,count \
              > $out_name
              
        done
          
    done

done



