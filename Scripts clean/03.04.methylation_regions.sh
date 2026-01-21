#!/bin/bash




# STILL UNDER DEVELOPMENT







# Script to compute mean methylation levels across defined genomic regions.

regions_name="promoters"
region_file=                  # Make sure it's in bed-reading format


path_positions=~/fil/Methylome/meth_proportions/
path_out=~/fil/Methylome/methylation_levels/"$regions_name"/


mkdir -p "$path_out"


# For 5mC and 5hmC
for mod in (5mC 5hmC) ; do 

    mkdir "$path_out"/"$mod"
    
    for sam in {1..6} ; do
    
        positions="$path_positions"/"$mod"/"$mod"_sample0"$sam"_mpct.txt
        out_name="$path_out"/"$mod"/meth_mean_"$mod"_rep"$sam"_"$regions_name".txt
        
        bedtools map \
          -a $region_file \
          -b $positions \
          -c 4 \
          -o mean \
          > $out_name
          
    done
      
done


# For 6mA
for sam in {1..6} ; do
  
    mkdir -p "$path_out"/6mA

    for st in ("pos" "neg") ; do

        positions="$path_positions"/6mA/6mA_sample0"$sam"_mpct_"$st".txt
        out_name="$path_out"/6mA/meth_mean_6mA_"$st"_rep"$sam"_"$regions_name".txt
        
        bedtools map \
          -a $region_file \
          -b $positions \
          -c 4 \
          -o mean \
          > $out_name
    done
      
done
