# Genome-wide DNA methylation profiling of the zebrafish forebrain
Pol Sorigue, Maeva Pinget, João Costa, Magda Teles, Rui F. Oliveira

DOI pending 

## Data availability
The raw sequencing data for this study have been deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession number **PRJEB108899**.
Intermediate files, including methylation matrices (output of script *02.pileup.sh*) and per-region methylation (output of script *06.04.methylation_regions*) have been deposited at ArrayExpress under accession **E-MTAB-16780**.
Date of data public release: March 26th, 2026.

## Annotations Used:
- Zebrafish genome assembly: **GRCz12tu** (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_049306965.1/)
- Zebrafish predicted CpG Islands annotation, downloaded Oct 22nd, 2025 (https://genome.ucsc.edu/cgi-bin/hgTables)


## Scripts index

1. **01 — Preprocessing & QC**
   - [01.01.BAM_qc.sh](<Scripts/01.Sequencing and BAM stats/01.01.BAM_qc.sh>) — Run Qualimap on aligned BAMs; extract read‑length statistics. (Tools: Qualimap v2.3; samtools v1.21 and awk from system)
   - [01.02.filter_BAMs.sh](<Scripts/01.Sequencing and BAM stats/01.02.filter_BAMs.sh>) — Filter BAMs by MAPQ and read length; sort, index, and produce MAPQ/read‑count stats. (Tools: samtools v1.21; modkit v0.5.0)
   - [01.03.BAM_performance.sh](<Scripts/01.Sequencing and BAM stats/01.03.BAM_performance.sh>) — Compute read loss statistics (secondary/supplementary, MAPQ, length). (Tools: samtools v1.21)
   - [01.04.bam_coverage_and_modifications.sh](<Scripts/01.Sequencing and BAM stats/01.04.bam_coverage_and_modifications.sh>) — Per‑sample and genome‑wide coverage; call modkit summaries. (Tools: samtools v1.21; modkit v0.5.0)

2. **02 — Pileup & per‑modification extraction**
   - [02.pileup.sh](<Scripts/02.pileup.sh>) — Run `modkit pileup` for 6mA, CpG, and CH contexts and generate per‑modification+motif outputs. (Tools: modkit v0.5.0)

3. **03 — CpG site counting**
   - [03.modkit_compute_CpGs.sh](<Scripts/03.modkit_compute_CpGs.sh>) — Compute CpG sites in the reference and per‑sample coverage values. (Tools: modkit v0.5.0; bedtools v2.31.1)

4. **04 — Motif proportions**
   - [04.modification_proportions.sh](<Scripts/04.modification_proportions.sh>) — Compute motif counts/fractions (CH and A motifs) using `modkit motif evaluate`. (Tools: modkit v0.5.0)

5. **05 — Region definitions**
   - [05.01.create_regions_50kbp.sh](<Scripts/05.Create regions files/05.01.create_regions_50kbp.sh>) — Create 50 kb genome windows (exclude mitochondrial chr). (Tools: bedtools v2.31.1)
   - [05.02.create_regions_promoters.R](<Scripts/05.Create regions files/05.02.create_regions_promoters.R>) — Create promoter regions (R). (R v4.5.2; packages GenomicFeatures v1.62, dplyr v1.1.4)

6. **06 — Methylation regions & strand symmetry**
   - [06.01.CpG_create_meth_regions.R](<Scripts/06.Methylation region and strand symmetry/06.01.CpG_create_meth_regions.R>) — Build CpG methylation proportion datasets (5mC/5hmC). (R v4.5.2; package data.table v1.18)
   - [06.02.CH_create_meth_regions.R](<Scripts/06.Methylation region and strand symmetry/06.02.CH_create_meth_regions.R>) — Build stranded CH methylation datasets. (R v4.5.2; package data.table v1.18)
   - [06.03.CH_single_strand_create_meth_regions.R](<Scripts/06.Methylation region and strand symmetry/06.03.CH_single_strand_create_meth_regions.R>) — Build CH datasets without strand information. (R v4.5.2; package data.table v1.18)
   - [06.04.methylation_regions.sh](<Scripts/06.Methylation region and strand symmetry/06.04.methylation_regions.sh>) — Map methylation tables onto genomic regions with `bedtools map`. (Tools: bedtools v2.31.1)
   - [06.05.strand_symmetry.R](<Scripts/06.Methylation region and strand symmetry/06.05.strand_symmetry.R>) — Analyze strand symmetry statistics across motifs and genomic partitions (R). (R v4.5.2; packages dplyr v1.1.4, tidyr v0.0.6, ggplot2 v4.0.1)

7. **07 — Whole‑brain & forebrain comparisons**
   - [07.01.WB_mapping.sh](<Scripts/07.Forebrain vs whole-brain CG islands/07.01.WB_mapping.sh>) — Bismark mapping of whole‑brain samples. (Tools: Bismark v0.25.1)
   - [07.02.WB_meth_calling.sh](<Scripts/07.Forebrain vs whole-brain CG islands/07.02.WB_meth_calling.sh>) — Run `bismark_methylation_extractor` and format outputs to per‑sample CpG beds. (Tools: Bismark v0.25.1)
   - [07.03.WB_meth_regions.sh](<Scripts/07.Forebrain vs whole-brain CG islands/07.03.WB_meth_regions.sh>) — Compute region‑level methylation for whole‑brain samples. (Tools: bedtools v2.31.1)
   - [07.04.prepare_FB_dataset.R](<Scripts/07.Forebrain vs whole-brain CG islands/07.04.prepare_FB_dataset.R>) — (R) Prepare forebrain dataset for downstream comparison. (R v4.5.2; packages data.table v1.18, dplyr v1.1.4)
   - [07.05.FB_meth_regions.sh](<Scripts/07.Forebrain vs whole-brain CG islands/07.05.FB_meth_regions.sh>) — Compute region‑level methylation for forebrain samples. (Tools: bedtools v2.31.1)
   - [07.06.FB_vs_WB_cgi.R](<Scripts/07.Forebrain vs whole-brain CG islands/07.06.FB_vs_WB_cgi.R>) — (R) Compare forebrain vs whole‑brain methylation at CpG island regions. (R v4.5.2; packages dplyr v1.1.4, tidyr v0.0.6, ggplot2 v4.0.1)
   - [07.07.overlap_cgi_genes.R](<Scripts/07.Forebrain vs whole-brain CG islands/07.07.overlap_cgi_genes.R>) — (R) Find genes overlapping CpG islands. (R v4.5.2; package GenomicRanges v1.62.1)

**Notes**

> Along the files and scripts, the notation “5mC” and “5hmC” refer to CpG-associated cytosine methylation, while “CH” refers to non‑CpG cytosine methylation.


### Directory structure assumed by the scripts

Large files are not included in this repository. They can be obtained from:
- Raw sequencing data: ENA (accession PRJEB108899) (see *Data availability* section)
- Intermediate files: ArrayExpress (accession E-MTAB-16780) (see *Data availability* section)
- Online Annotations (see *Annotations Used* section)

All scripts assume the following directory structure under `$HOME`.  
After downloading the raw data, the pipeline can be run sequentially using the provided scripts.

```bash
$HOME/
├── Data/
│   ├── bams_genomics/
│   ├── bams_filtered/
│   ├── pileups/
│   │   ├── pileup_all/
│   │   ├── pileup_cpg/
│   │   └── pileup_ch/
│   └── datasets_by_mod/
│       ├── 5mC/
│       ├── 5hmC/
│       ├── 6mA/
│       └── ch/                # also used for CH single‑strand
│
├── Data_methylation/           
│   ├── datasets_proportions/
│   │   ├── 5mC/
│   │   ├── 5hmC/
│   │   ├── 6mA/
│   │   ├── ch/
│   │   └── ch_ss/            # single‑strand output
│   │
│   └── methylation_regions/
│       ├── output/
│       │  ├── cgi
│       │  ├── genes
│       │  ├── genome_50kb_bins
│       │  └── promoters
│       ├── regions/			# region coordinates templates
│       └── stats/
│
├── Ref_genome/
│   └── GRCz12tu/
│       ├── genomic.gtf
│       ├── GCF_049306965.1_GRCz12tu_genomic.fna
│       ├── GCF_049306965.1_GRCz12tu_genomic.fna.fai   # genome index
│       └── chr_array.txt       
│
├── QC/
│   ├── bam_qc/
│   │   └── qualimap/
│   └── qc_stats/
│
├── methylation_proportions/
│   └── motifs_CH.txt
│
├── CpG_sites/
│   ├── coverage_CpG_samples/
│   └── CpG_genome/
│       └── CpG_sites_genome.bed
│
└── Chaterjee/
    ├── bismark_data_process/
    │   ├── fastq_runs/
    │   ├── mapping/
    │   └── meth_calls/
    │       └── output/
    ├── forebrain_cpg_sites
    └── methylation_cgi/
        ├── forebrain/
        └── whole-brain/
```
