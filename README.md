# Zebrafish-Forebrain-Methylome
Scripts used in Scientific Publication

---

## Project index — Scripts ("Scripts clean") 🔧
This repository contains analysis scripts used for the zebrafish forebrain methylome project. The main pipeline scripts are located in the `Scripts clean/` folder and are ordered by their numeric prefixes. Below is an index with brief descriptions and **direct links** to each script.

> Note: files use numbered prefixes to indicate execution order. Filenames may have been renumbered to maintain insertion of new steps.

---

### Core preprocessing & QC (01–01.04) ✅
- [01.01.BAM_qc.sh](Scripts%20clean/01.01.BAM_qc.sh) — Run Qualimap on aligned BAMs; extract read-length statistics. 🔍
- [01.02.filter_BAMs.sh](Scripts%20clean/01.02.filter_BAMs.sh) — Filter BAMs by MAPQ and read length; sort, index, and produce MAPQ/read-count stats. 🧹
- [01.03.BAM_performance.sh](Scripts%20clean/01.03.BAM_performance.sh) — Compute read loss statistics (secondary/supplementary, MAPQ, length). 📈
- [01.04.bam_coverage_and_modifications.sh](Scripts%20clean/01.04.bam_coverage_and_modifications.sh) — Per-sample and genome-wide coverage; call modkit summary for modification summaries. 🧬

### Pileup & motif analyses (02–04) 🧾
- [02.pileup.sh](Scripts%20clean/02.pileup.sh) — Run `modkit pileup` for 6mA, CpG, and CH contexts and generate per-modification outputs. ⚙️
- [03.modkit_compute_CpGs.sh](Scripts%20clean/03.modkit_compute_CpGs.sh) — Compute CpG sites in the reference and per-sample coverage counts. 🧮
- [04.modification_proportions.sh](Scripts%20clean/04.modification_proportions.sh) — Compute motif counts/fractions (CH motifs) using `modkit motif evaluate`. 📊

### Region definitions & methylation region creation (05–06) 🗺️
- [05.01.create_regions_50kbp.sh](Scripts%20clean/05.01.create_regions_50kbp.sh) — Create 50 kb genome windows (exclude MT). 🧩
- [05.02.create_regions_promoters.R](Scripts%20clean/05.02.create_regions_promoters.R) — Create promoter regions (R). 🧾
- [06.03.1.CpG_create_meth_regions.R](Scripts%20clean/06.03.1.CpG_create_meth_regions.R) — Create CpG methylation regions (R). 🧭
- [06.03.2.A_create_meth_regions_6mA.R](Scripts%20clean/06.03.2.A_create_meth_regions_6mA.R) — Create 6mA methylation regions (R). 🧭
- [06.03.3.CH_create_meth_regions.R](Scripts%20clean/06.03.3.CH_create_meth_regions.R) — Create CH methylation regions (R). 🧭
- [06.03.4.CH_single_strand_create_meth_regions.R](Scripts%20clean/06.03.4.CH_single_strand_create_meth_regions.R) — CH single-strand regions (R). 🧭
- [06.04.methylation_regions.sh](Scripts%20clean/06.04.methylation_regions.sh) — Map methylation datasets onto regions using `bedtools map`. 🗂️

### Whole-brain & Forebrain comparisons (07) 🔬
- [07.01.WB_mapping.sh](Scripts%20clean/07.01.WB_mapping.sh) — Bismark mapping of whole-brain samples. 🧩
- [07.02.WB_meth_calling.sh](Scripts%20clean/07.02.WB_meth_calling.sh) — Run `bismark_methylation_extractor` and format outputs to per-sample CpG beds. 🧾
- [07.03.WB_meth_regions.sh](Scripts%20clean/07.03.WB_meth_regions.sh) — Compute region-level methylation for whole-brain samples (bedtools map). 📐
- [07.04.prepare_FB_dataset.R](Scripts%20clean/07.04.prepare_FB_dataset.R) — (R) Prepare forebrain dataset for downstream comparison. 🔁
- [07.05.FB_meth_regions.sh](Scripts%20clean/07.05.FB_meth_regions.sh) — Compute region-level methylation for forebrain samples. 📐
- [07.06.FB_vs_WB_cgi.R](Scripts%20clean/07.06.FB_vs_WB_cgi.R) — (R) Compare forebrain vs whole-brain methylation at CGI regions. 🔬

---

## Conventions & notes ✍️
- Variable naming in shell scripts follows **UPPER_SNAKE_CASE** for constants and paths (e.g., `PATH_BAMS`, `THREADS`) and `SAMPLES`/`SAMPLE` for loops. This change was applied to improve reproducibility and readability.
- All file-path expansions are quoted (`"${VAR}"`) to avoid whitespace bugs.
- Scripts are intended to be run in the numeric order indicated by their filename prefixes; some scripts assume outputs from previous steps.

> Tip: run scripts stepwise on a single sample to validate paths and dependencies before running the full dataset.

---

If you'd like, I can:
- Add a **high-level workflow diagram** (MD) showing dependencies, or
- Update internal script references to reflect the recent filename renames.

If you want changes or different organization, tell me which additional details to include and I’ll update `README.md`. ✨