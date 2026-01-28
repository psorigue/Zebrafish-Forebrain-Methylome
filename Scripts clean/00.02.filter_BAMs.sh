#!/bin/bash

# This script filters BAM files based on mapping quality and read length, then sorts and indexes them. It also generates statistics on read counts before and after filtering, as well as the distribution of mapping qualities.
# Samtools v 1.21
# Modkit v 0.5.0

# Set input/output paths
path_in=~/fil/Epi/Data_genomics/bams_genomics/Bam_Aligned_Cat/
path_out=~/fil/Epi/Data_genomics/bams_filtered/
stats_dir=~/fil/Methylome/QC/qc_stats/

# Define output files for statistics
counts_out="${stats_dir}/read_counts.tsv"
mapq_out="${stats_dir}/mapq_distribution.tsv"

# Define parameters
thr=8 # Working threads
mapq_accepted=10 # Minimum mapping quality
min_len=200 # Minimum read length

# Define samples to process
samples=(seq 1 6) 

mkdir -p "${path_out}" "${stats_dir}"

# Headers for output files
echo -e "sample\tbefore\tafter" > "${counts_out}"
echo -e "sample\tmapq" > "${mapq_out}"

for sam in "${samples[@]}"; do
    
    # Zero-pad the number to match file naming 
    sam_padded=$(printf "%02d" "$sam")
    
    # Define input BAM file
    bam="${path_in}/barcode${sam_padded}_aligned.bam"
    
    # Define sample name and output BAM file
    sample=$(basename "${bam}" .bam)
    outbam="${path_out}/${sample}_filt.bam"

    # Count initial number of reads
    before=$(samtools view -c "${bam}")

    # Collect MAPQ distribution (before filtering)
    samtools view "${bam}" \
      | awk -v s="${sample}" '{print s "\t" $5}' \
      >> "${mapq_out}"

    # Filter, sort. '-F 2308' removes unmapped reads, secondary and supplementary aligned reads. '-q' keeps only reads with quality above 10. '-e' filters by read length.
    samtools view -@ ${thr} -b \
        -q 10 \
        -F 2308 \
        -e 'length(seq) >= ${min_len}' \
        "${bam}" \
      | samtools sort -@ ${thr} -T "${path_out}/${sample}_tmp" -o "${outbam}"
    
    # Index filtered BAM
    samtools index -@ ${thr} "${outbam}"

    # Count final number of reads
    after=$(samtools view -c "${outbam}")

    # Output read counts before and after filtering
    echo -e "${sample}\t${before}\t${after}" >> "${counts_out}"

done

exit
