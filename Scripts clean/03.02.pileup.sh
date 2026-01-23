#!/bin/bash

#SBATCH --job-name=cpg_pileup
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pol.tomas@gimm.pt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=1gb
#SBATCH --time=14-00:00:00
#SBATCH --output=cpg_pileup.log

date

export SINGULARITY_BIND="/ifs"




# STILL NEED TO CLEAN THIS SCRIPT



# I keep this script because I'll have to show it, but I use the pileups I did for Epi.

# To avoid running all this, run tmp.sbatch script, which will use the pileups already done and will create the dataset_by_mod files

# The only IMPORTANT differences are
# - Here I derive h mod to dataset_by_mod from the --cpg pileup
# - Here I don't rename chromosomes to keep matching modkit


# =========================
# Paths
# =========================
path_bam=~/fil/Epi/Data_genomics/bams_filtered/

path_pileup_a=~/fil/Epi/Data_processed/pileup_all/
path_pileup_cpg=~/fil/Epi/Data_processed/pileup_cpg/

path_out_5mC=~/fil/Methylome/Data_methylation/datasets_by_mod/5mC/
path_out_5hmC=~/fil/Methylome/Data_methylation/datasets_by_mod/5hmC/
path_out_6mA=~/fil/Methylome/Data_methylation/datasets_by_mod/6mA/

chr_mapping_file=~/fil/Epi/Data_processed/chr_mapping.tsv
ref_genome=~/fil/Epi/Ref_genome/GRCz12tu/genome/GCF_049306965.1_GRCz12tu_genomic.fna

thr=20

samples=(2 3 4 5)

# =========================
# Setup
# =========================
mkdir -p "$path_pileup_a" "$path_pileup_cpg"
mkdir -p "$path_out_5mC" "$path_out_5hmC" "$path_out_6mA"

export path_out_5mC="$path_out_5mC"
export path_out_5hmC="$path_out_5hmC"
export path_out_6mA="$path_out_6mA"

# =========================
# Main loop
# =========================
for number in "${samples[@]}"; do
    num=$(printf "%02d" "$number")
    bam="${path_bam}/barcode${num}_aligned_filt.bam"

    echo "Processing sample ${num}"

    # =======================
    # 1) FIRST pileup: a mod
    # =======================
    cd "$path_pileup_a"

    modkit pileup --threads "$thr" "$bam" "pileup${num}.bed" --log-filepath "pileup${num}.log" --filter-threshold C:0.8 --filter-threshold A:0.8

    # Split ONLY h + a. Rename chromosomes
    awk -F'\t' -v OFS='\t' -v i="$num" '
    $4=="a" {
        print $1,$2,$3,$10,$12,$6 >> sprintf("%s/6mA_sample%s.txt", ENVIRON["path_out_6mA"], i)
    }
    ' "pileup${num}.bed"

    bgzip --threads "$thr" "pileup${num}.bed"
    tabix "pileup${num}.bed.gz"

    # =========================
    # 2) SECOND pileup: CpG + combined strands
    # =========================
    cd "$path_pileup_cpg"

    modkit pileup --threads "$thr" "$bam" "pileup${num}_cpg.bed" --log-filepath "pileup${num}_cpg.log" --ref $ref_genome --filter-threshold C:0.8 --filter-threshold A:0.8 --combine-strands --cpg

    # Extract ONLY m from second pileup. Rename chromosomes
    awk -F'\t' -v OFS='\t' -v i="$num" '
    $4=="m" {
        print $1,$2,$3,$10,$12,$6 >> sprintf("%s/5mC_sample%s.txt", ENVIRON["path_out_5mC"], i)
    }
    $4=="h" {
        print $1,$2,$3,$10,$12,$6 >> sprintf("%s/5hmC_sample%s.txt", ENVIRON["path_out_5hmC"], i)
    }
    ' "pileup${num}_cpg.bed"

    bgzip --threads "$thr" "pileup${num}_cpg.bed"
    tabix "pileup${num}_cpg.bed.gz"

    echo "Finished sample ${num}"
done

cat /etc/os-release
date
exit