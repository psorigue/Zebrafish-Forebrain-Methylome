#!/bin/bash

# The script performs three separate pileup analyses per sample BAM file:
# 1) A general pileup for all modified bases. From this pileup we extract 6mA data.
# 2) A pileup specifically for CpG sites, combining strands, to extract 5mC and 5hmC data.
# 3) A pileup for CH sites to extract 5mC data in non-CpG contexts.
# 4) After each pileup, it extracts relevant modification data and renames chromosomes based on a mapping file.
# Modkit version: 0.5.0

# Comuptaionally heavy script - recommended to run with multiple threads

# =========================
# Paths
# =========================
PATH_BAM="$HOME/Data/bams_filtered"

PATH_PILEUP_A="$HOME/Data/pileups/pileup_all"
PATH_PILEUP_CPG="$HOME/Data/pileups/pileup_cpg"
PATH_PILEUP_CH="$HOME/Data/pileups/pileup_ch"

PATH_OUT_5MC="$HOME/Data/datasets_by_mod/5mC"
PATH_OUT_5HMC="$HOME/Data/datasets_by_mod/5hmC"
PATH_OUT_6MA="$HOME/Data/datasets_by_mod/6mA"
PATH_OUT_CH="$HOME/Data/datasets_by_mod/ch"

CHR_MAPPING_FILE="$HOME/Data/chr_mapping.tsv"
REF_GENOME="$HOME/Ref_genome/GRCz12tu/GCF_049306965.1_GRCz12tu_genomic.fna"

THREADS=20


# =========================
# Setup
# =========================
mkdir -p "${PATH_PILEUP_A}" "${PATH_PILEUP_CPG}"
mkdir -p "${PATH_OUT_5MC}" "${PATH_OUT_5HMC}" "${PATH_OUT_6MA}"

# Export paths to import to awk
export PATH_OUT_5MC="${PATH_OUT_5MC}"
export PATH_OUT_5HMC="${PATH_OUT_5HMC}"
export PATH_OUT_6MA="${PATH_OUT_6MA}"

# =========================
# Main loop
# =========================
for SAMPLE in {01..06}; do
   
    BAM_FILE="${PATH_BAM}/barcode${SAMPLE}_aligned_filt.bam"

    echo "Processing sample ${SAMPLE}"

    # ==================================
    # 1) FIRST pileup: all modficiations
    # ==================================
    cd "${PATH_PILEUP_A}"

    # Perform pileup for all mods
    modkit pileup --threads "${THREADS}" "${BAM_FILE}" "pileup${SAMPLE}.bed" --log-filepath "pileup${SAMPLE}.log" --filter-threshold C:0.8 --filter-threshold A:0.8

    # Split modification 'a' and rename chromosomes
    awk -F'\t' -v OFS='\t' -v i="${SAMPLE}" '
    $4=="a" {
        print $1,$2,$3,$10,$12,$6 >> sprintf("%s/6mA_sample%s.txt", ENVIRON["PATH_OUT_6MA"], i)
    }
    ' "pileup${SAMPLE}.bed"

    # Compress and index
    bgzip --threads "${THREADS}" "pileup${SAMPLE}.bed"
    tabix "pileup${SAMPLE}.bed.gz"

    # ========================================
    # 2) SECOND pileup: CpG + combined strands
    # ========================================
    cd "${PATH_PILEUP_CPG}"

    # Perform pileup for CpG sites only, combining strands
    modkit pileup --threads "${THREADS}" "${BAM_FILE}" "pileup${SAMPLE}_cpg.bed" --log-filepath "pileup${SAMPLE}_cpg.log" --ref "${REF_GENOME}" --filter-threshold C:0.8 --filter-threshold A:0.8 --combine-strands --cpg

    # Extract m and h and rename chromosomes
    awk -F'\t' -v OFS='\t' -v i="${SAMPLE}" '
    $4=="m" {
        print $1,$2,$3,$10,$12,$6 >> sprintf("%s/5mC_sample%s.txt", ENVIRON["PATH_OUT_5MC"], i)
    }
    $4=="h" {
        print $1,$2,$3,$10,$12,$6 >> sprintf("%s/5hmC_sample%s.txt", ENVIRON["PATH_OUT_5HMC"], i)
    }
    ' "pileup${SAMPLE}_cpg.bed"

    # Compress and index
    bgzip --threads "${THREADS}" "pileup${SAMPLE}_cpg.bed"
    tabix "pileup${SAMPLE}_cpg.bed.gz"

    # =======================
    # 3) THIRD pileup: CH
    # =======================
    cd "${PATH_PILEUP_CH}"

    # Perform pileup for CH sites only
    modkit pileup --threads "${THREADS}" "${BAM_FILE}" "pileup${SAMPLE}_ch.bed" --ref "${REF_GENOME}" --log-filepath "pileup${SAMPLE}_ch.log" --filter-threshold C:0.8 --ignore h --motif CA 0 --motif CT 0 --motif CC 0

    # Extract m and rename chromosomes. In this case, the column $4 is always 'm' because we ignored 'h' in the pileup. This column contains the motif info (CA, CT, CC).
    awk -F'\t' -v OFS='\t' -v mapping="${CHR_MAPPING_FILE}" -v i="${SAMPLE}" '
    BEGIN {
        while ((getline < mapping) > 0) map[$1]=$2
        close(mapping)
    }
    NR==1 {
        out = sprintf("%s/ch_sample%s.txt", ENVIRON["PATH_OUT_CH"], i)
        system("rm -f " out)
    }
    $4 ~ /^m,/ {     
        chr=$1
        for (p in map) {
            while (index(chr,p)) {
                chr = substr(chr,1,index(chr,p)-1) map[p] substr(chr,index(chr,p)+length(p))
            }
        }
        print chr,$2,$3,$10,$12,$6,$4 >> out
    }
    ' "pileup${SAMPLE}_ch.bed"

    # Compress and index
    bgzip --threads "${THREADS}" "pileup${SAMPLE}_ch.bed"
    tabix "pileup${SAMPLE}_ch.bed.gz"
    
    echo "Finished sample ${SAMPLE}"
done


date
exit