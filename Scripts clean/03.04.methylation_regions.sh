#!/bin/bash

# Script to calculate mean methylation levels in defined genomic regions (e.g., 50kb bins, genes, promoters, CGIs) for 5mC, 5hmC, 6mA, and CH contexts.
# Bedtools XXX

# Define regions
regions_name="genome_50kb_bins" # genome_50kb_bins, genes, promoters, cgi
region_file=~/fil/Methylome/methylation_regions/regions/"$regions_name".bed 

# Define paths
path_positions=~/fil/Methylome/Data_methylation/meth_proportions/
path_out=~/fil/Methylome/methylation_regions/"$regions_name"/


mkdir -p "$path_out"


# For 5mC and 5hmC
for mod in 5mC 5hmC ; do 

    mkdir "$path_out"/"$mod"
    
    for sam in {1..6} ; do
    
        positions="$path_positions"/"$mod"/"$mod"_sample0"$sam"_mpct.txt
        out_name="$path_out"/"$mod"/meth_mean_"$mod"_rep"$sam"_"$regions_name".txt
        
        
        # Coverage is already filtered in previous script to > or = 5 
    
        # Output: chr, start (0-based), end (0-based), name of -a region (if applicable), mean methylation %, mean coverage, number of covered cpg positions in the region
        bedtools map \
          -a <(sort -k1,1 -k2,2n "$region_file") \
          -b "$positions" \
          -c 4,5,4 \
          -o mean,mean,count \
          > $out_name
          
    done
      
done


# For 6mA
for sam in {1..6} ; do
  
    mkdir -p "$path_out"/6mA

    for st in pos neg ; do

        positions="$path_positions"/6mA/6mA_sample0"$sam"_mpct_"$st".txt
        out_name="$path_out"/6mA/meth_mean_6mA_"$st"_rep"$sam"_"$regions_name".txt
        
        
        # Coverage is already filtered in previous script to > or = 5 
    
        # Output: chr, start (0-based), end (0-based), name of -a region (if applicable), mean methylation %, mean coverage, number of covered cpg positions in the region
        bedtools map \
          -a <(sort -k1,1 -k2,2n "$region_file") \
          -b "$positions" \
          -c 4,5,4 \
          -o mean,mean,count \
          > $out_name
          
    done
      
done


# For CH
for motif in CA CC CT ; do

    for sam in {1..6} ; do
      
        mkdir -p "$path_out"/ch
    
        for st in pos neg ; do
    
            positions="$path_positions"/ch/ch_sample0"$sam"_"$motif"_"$st"_mpct.txt
            out_name="$path_out"/ch/meth_mean_"$motif"_"$st"_rep"$sam"_"$regions_name".txt
            
            
            # Coverage is already filtered in previous script to > or = 5 
        
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



