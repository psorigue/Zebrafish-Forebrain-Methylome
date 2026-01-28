#!/bin/bash

# Script to compute BAM coverage statistics and modification summaries using samtools and modkit.
# Samtools v1.21
# Modkit v0.5.0

# Define paths
path_bam=~/fil/Epi/Data_genomics/bams_filtered/
path_stats=~/fil/Methylome/QC/bam_stats/

# Set number of threads
thr=2

mkdir -p "$path_stats" "$path_stats"/ind_samples/

# 1. Coverage for Individual samples
for num in {01..06} ; do
    
    bam="${path_bam}/barcode${num}_aligned_filt.bam"
    echo $bam
    # Per-chromosome coverage. Output: 
    samtools coverage $bam > "${path_stats}/ind_samples/barcode${num}.stats"
    
done


# 2. Coverage for all samples together
## Per-chromosome coverage values
samtools coverage "$path_bam"/barcode{01..06}_aligned_filt.bam > "${path_stats}/all_samples.stats"

# Genome-wide coverage values
samtools coverage "$path_bam"/barcode{01..06}_aligned_filt.bam \
    | awk '
    NR==1 { next }                          # skip header
    $1 == "NC_002333.2" { next }             # exclude mitochondrial DNA
    {
      len        += $3                       # contig length (start=1)
      covbases   += $5                       # bases =1X
      depth_sum  += $7 * $3                  # length-weighted depth
      baseq_sum  += $8 * $5                  # base-weighted baseQ
      mapq_sum   += $9 * $4                  # read-weighted mapQ
      reads      += $4                       # total reads
    }
    END {
      print "Fraction_genome_covered:", covbases / len
      print "Genome_coverage_%:",      100 * covbases / len
      print "Mean_depth:",             depth_sum / len
      print "Mean_baseQ:",             baseq_sum / covbases
      print "Mean_mapQ:",              mapq_sum / reads
    }' > "${path_stats}/genome_general_stats.txt"
    

# 3. Summary of the modifications:
for num in {01..06} ; do

    bam="${path_bam}/barcode${num}_aligned_filt.bam"

    modkit summary --filter-threshold C:0.8 --filter-threshold A:0.8 -n 100000 -t $thr $bam
    
done > "$path_stats"/modifications_summary_modkit.txt


exit