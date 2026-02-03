# Zebrafish-Forebrain-Methylome
Scripts used in Scientific Publication

## Scripts index

1. **01 — Preprocessing & QC**
   - [01.01.BAM_qc.sh](Scripts/01.01.BAM_qc.sh) — Run Qualimap on aligned BAMs; extract read-length statistics.
   - [01.02.filter_BAMs.sh](Scripts/01.02.filter_BAMs.sh) — Filter BAMs by MAPQ and read length; sort, index, and produce MAPQ/read-count stats.
   - [01.03.BAM_performance.sh](Scripts/01.03.BAM_performance.sh) — Compute read loss statistics (secondary/supplementary, MAPQ, length).
   - [01.04.bam_coverage_and_modifications.sh](Scripts/01.04.bam_coverage_and_modifications.sh) — Per-sample and genome-wide coverage; call modkit summaries.

2. **02 — Pileup & per-modification extraction**
   - [02.pileup.sh](Scripts/02.pileup.sh) — Run `modkit pileup` for 6mA, CpG, and CH contexts and generate per-modification+motif outputs.

3. **03 — CpG site counting**
   - [03.modkit_compute_CpGs.sh](Scripts/03.modkit_compute_CpGs.sh) — Compute CpG sites in the reference and per-sample coverage counts.

4. **04 — Motif proportions**
   - [04.modification_proportions.sh](Scripts/04.modification_proportions.sh) — Compute motif counts/fractions (CH and A motifs) using `modkit motif evaluate`.

5. **05 — Region definitions**
   - [05.01.create_regions_50kbp.sh](Scripts/05.01.create_regions_50kbp.sh) — Create 50 kb genome windows (exclude mitochondrial chr).
   - [05.02.create_regions_promoters.R](Scripts/05.02.create_regions_promoters.R) — Create promoter regions (R).

6. **06 — Methylation region creation**
   - [06.03.1.CpG_create_meth_regions.R](Scripts%20clean/06.03.1.CpG_create_meth_regions.R) — Create dataset for CpG motif, to be used in script 06.04 (R).
   - [06.03.2.A_create_meth_regions_6mA.R](Scripts%20clean/06.03.2.A_create_meth_regions_6mA.R) — Create dataset for A motif, with strand information, to be used in script 06.04 (R).
   - [06.03.3.CH_create_meth_regions.R](Scripts%20clean/06.03.3.CH_create_meth_regions.R) — Create dataset for nonCpG CH motif, with strand information, to be used in script 06.04 (R).
   - [06.03.4.CH_single_strand_create_meth_regions.R](Scripts%20clean/06.03.4.CH_single_strand_create_meth_regions.R) — Create dataset for nonCpG CH motif, without strand information, to be used in script 06.04 (R).
   - [06.04.methylation_regions.sh](Scripts%20clean/06.04.methylation_regions.sh) — Map methylation datasets onto regions using `bedtools map`.

7. **07 — Whole-brain & Forebrain comparisons**
   - [07.01.WB_mapping.sh](Scripts/07.01.WB_mapping.sh) — Bismark mapping of whole-brain samples.
   - [07.02.WB_meth_calling.sh](Scripts/07.02.WB_meth_calling.sh) — Run `bismark_methylation_extractor` and format outputs to per-sample CpG beds.
   - [07.03.WB_meth_regions.sh](Scripts/07.03.WB_meth_regions.sh) — Compute region-level methylation for whole-brain samples (bedtools map).
   - [07.04.prepare_FB_dataset.R](Scripts/07.04.prepare_FB_dataset.R) — (R) Prepare forebrain dataset for downstream comparison.
   - [07.05.FB_meth_regions.sh](Scripts/07.05.FB_meth_regions.sh) — Compute region-level methylation for forebrain samples.
   - [07.06.FB_vs_WB_cgi.R](Scripts/07.06.FB_vs_WB_cgi.R) — (R) Compare forebrain vs whole-brain methylation at CGI regions.

**Notes**
Along the files and scripts, the notation "5mC" and "5hmC" refer to CpG-associated cytosine methylation, while "CH" refer to nonCpG-associated cytosine methylation