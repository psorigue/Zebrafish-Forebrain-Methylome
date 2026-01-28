#!/bin/bash

#SBATCH --job-name=pileup
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pol.tomas@gimm.pt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=1gb
#SBATCH --time=14-00:00:00
#SBATCH --output=pileup.log

date
export SINGULARITY_BIND="/ifs"

# The script performs three separate pileup analyses per sample BAM file:
# 1) A general pileup for all modified bases to extract 6mA data.
# 2) A pileup specifically for CpG sites, combining strands, to extract 5mC and 5hmC data.
# 3) A pileup for CH sites to extract 5mC data in non-CpG contexts.
# Modkit version: 0.5.0

# =========================
# Paths
# =========================
path_bam=~/fil/Epi/Data_genomics/bams_filtered/

path_pileup_a=~/fil/Epi/Data_processed/pileups/pileup_all/
path_pileup_cpg=~/fil/Epi/Data_processed/pileups/pileup_cpg/
path_pileup_ch=~/fil/Epi/Data_processed/pileups/pileup_ch/

path_out_5mC=~/fil/Methylome/Data_methylation/datasets_by_mod/5mC/
path_out_5hmC=~/fil/Methylome/Data_methylation/datasets_by_mod/5hmC/
path_out_6mA=~/fil/Methylome/Data_methylation/datasets_by_mod/6mA/
path_out_ch=~/fil/Methylome/Data_methylation/datasets_by_mod/ch/

chr_mapping_file=~/fil/Epi/Data_processed/chr_mapping.tsv
ref_genome=~/fil/Epi/Ref_genome/GRCz12tu/genome/GCF_049306965.1_GRCz12tu_genomic.fna

thr=20


# =========================
# Setup
# =========================
mkdir -p "$path_pileup_a" "$path_pileup_cpg"
mkdir -p "$path_out_5mC" "$path_out_5hmC" "$path_out_6mA"

# Export paths to import to awk
export path_out_5mC="$path_out_5mC"
export path_out_5hmC="$path_out_5hmC"
export path_out_6mA="$path_out_6mA"

# =========================
# Main loop
# =========================
for num in {01..06}; do
   
    bam="${path_bam}/barcode${num}_aligned_filt.bam"

    echo "Processing sample ${num}"

    # =======================
    # 1) FIRST pileup: a mod
    # =======================
    cd "$path_pileup_a"

    # Perform pileup for all mods
    modkit pileup --threads "$thr" "$bam" "pileup${num}.bed" --log-filepath "pileup${num}.log" --filter-threshold C:0.8 --filter-threshold A:0.8

    # Split modification 'a' and rename chromosomes
    awk -F'\t' -v OFS='\t' -v i="$num" '
    $4=="a" {
        print $1,$2,$3,$10,$12,$6 >> sprintf("%s/6mA_sample%s.txt", ENVIRON["path_out_6mA"], i)
    }
    ' "pileup${num}.bed"

    # Compress and index
    bgzip --threads "$thr" "pileup${num}.bed"
    tabix "pileup${num}.bed.gz"

    # ========================================
    # 2) SECOND pileup: CpG + combined strands
    # ========================================
    cd "$path_pileup_cpg"

    # Perform pileup for CpG sites only, combining strands
    modkit pileup --threads "$thr" "$bam" "pileup${num}_cpg.bed" --log-filepath "pileup${num}_cpg.log" --ref $ref_genome --filter-threshold C:0.8 --filter-threshold A:0.8 --combine-strands --cpg

    # Extract m and h and rename chromosomes
    awk -F'\t' -v OFS='\t' -v i="$num" '
    $4=="m" {
        print $1,$2,$3,$10,$12,$6 >> sprintf("%s/5mC_sample%s.txt", ENVIRON["path_out_5mC"], i)
    }
    $4=="h" {
        print $1,$2,$3,$10,$12,$6 >> sprintf("%s/5hmC_sample%s.txt", ENVIRON["path_out_5hmC"], i)
    }
    ' "pileup${num}_cpg.bed"

    # Compress and index
    bgzip --threads "$thr" "pileup${num}_cpg.bed"
    tabix "pileup${num}_cpg.bed.gz"

    # =======================
    # 3) THIRD pileup: CH
    # =======================
    cd "$path_pileup_ch"

    # Perform pileup for CH sites only
    modkit pileup --threads "$thr" "$bam" "pileup${num}_ch.bed" --ref $ref_genome --log-filepath "pileup${num}_ch.log" --filter-threshold C:0.8 --ignore h --motif CA 0 --motif CT 0 --motif CC 0

    # Extract m and rename chromosomes. In this case, the column $4 is always 'm' because we ignored 'h' in the pileup. This column contains the motif info (CA, CT, CC).
    awk -F'\t' -v OFS='\t' -v mapping="$chr_mapping_file" -v i="$num" '
    BEGIN {
        while ((getline < mapping) > 0) map[$1]=$2
        close(mapping)
    }
    NR==1 {
        out = sprintf("%s/ch_sample%s.txt", ENVIRON["path_out_ch"], i)
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
    ' "pileup${num}_ch.bed"

    # Compress and index
    bgzip --threads "$thr" "pileup${num}_ch.bed"
    tabix "pileup${num}_ch.bed.gz"
    
    echo "Finished sample ${num}"
done


date
exit