#!/bin/bash

# Script to compute the proportions of different modifications (m, h, a) and different motifs (CG, CA, CC, CT, A) in the zebrafish forebrain methylome samples using modkit's motif evaluation function.
# Modkit v0.5.0

# Set variables
PATH_PILEUP="$HOME/fil/Epi/Data_processed/pileups"
PATH_OUT="$HOME/fil/Methylome/methylation_proportions"
REF_GENOME="$HOME/fil/Epi/Ref_genome/GRCz12tu/genome/GCF_049306965.1_GRCz12tu_genomic.fna"
PATH_MOTIFS_TABLE="$HOME/fil/Methylome/methylation_proportions/motifs_CH.txt"
#PILEUP_CH="${PATH_PILEUP}/pileup_ch" Analyze one pilup at a time
PILEUP_ALL="${PATH_PILEUP}/pileup_all"


THREADS=6
SAMPLES=( 1 2 3 4 5 6 )

cd "${PATH_OUT}"

# Compute modification proportions
for SAMPLE in "${SAMPLES[@]}" ; do

    modkit motif evaluate \
        --known-motifs-table "${PATH_MOTIFS_TABLE}" \
        --out "mod_counts_fractions_${SAMPLE}.tsv" \
        --min-coverage 5 \
        --low-thresh 0.2 \
        --high-thresh 0.6 \
        -t "${THREADS}" \
        --log-filepath "mod_counts_fractions_rep${SAMPLE}.log" \
        --in-bedmethyl "${PILEUP_ALL}/pileup0${SAMPLE}.bed.gz" \
        --ref "${REF_GENOME}"

done

date
exit