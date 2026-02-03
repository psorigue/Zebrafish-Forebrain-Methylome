#!/bin/bash

# Script to compute BAM coverage statistics and modification summaries using samtools and modkit.
# Samtools v1.21
# Modkit v0.5.0

# Define paths
PATH_BAM="$HOME/fil/Epi/Data_genomics/bams_filtered"
PATH_STATS="$HOME/fil/Methylome/QC/bam_stats"

# Set number of threads
THREADS=2

# Create output directories
mkdir -p "${PATH_STATS}" "${PATH_STATS}/ind_samples"

# 1. Coverage for Individual samples
for SAMPLE in {01..06} ; do
    
    BAM_FILE="${PATH_BAM}/barcode${SAMPLE}_aligned_filt.bam"
    echo "${BAM_FILE}"
    # Per-chromosome coverage. Output: 
    samtools coverage "${BAM_FILE}" > "${PATH_STATS}/ind_samples/barcode${SAMPLE}.stats"
    
done

# 2. Coverage for all samples together
# Per-chromosome coverage values
samtools coverage "${PATH_BAM}"/barcode{01..06}_aligned_filt.bam > "${PATH_STATS}/all_samples.stats"

# Genome-wide coverage values
samtools coverage "${PATH_BAM}"/barcode{01..06}_aligned_filt.bam \
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
    }' > "${PATH_STATS}/genome_general_stats.txt"
    

# 3. Summary of the modifications (built-in modkit function)
for SAMPLE in {01..06} ; do

    # Define input BAM file
    BAM_FILE="${PATH_BAM}/barcode${SAMPLE}_aligned_filt.bam"
    echo "sample ${SAMPLE}"

    # Compute modification summary
    modkit summary --filter-threshold C:0.8 --filter-threshold A:0.8 -n 100000 -t "${THREADS}" "${BAM_FILE}"
    
done > "${PATH_STATS}/modifications_summary_modkit.txt"


exit