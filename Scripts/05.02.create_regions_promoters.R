# This script creates promoters regions (2kb upstream TSS) in BED format for later methylation level calculation. 
# It uses the GTF annotation file to create a TxDb object and extract promoter regions.

library(GenomicFeatures) # version 1.62
library(dplyr) # version 1.1.4

# Define output file and path
out_file <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/regions/promoters.bed"
path_genome <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Epi/Ref_genome/GRCz12tu/"

# Load GTF annotation file and create TxDb object
file_annot <- paste0(path_genome, "genome/genomic.gtf")
txdb <- txdbmaker::makeTxDbFromGFF(file_annot)
# Create transcript -> gene link
txdf <- AnnotationDbi::select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")

# Gene annotation info, downloaded from NCBI
gene_ranges <- read.csv(paste0(path_genome, "danrer_gene_info.tsv"),
                        sep = "\t")

# Make promoter dataset
# 1. Define regions into a GRanges object
promoters <- promoters(txdb, upstream=2000, downstream=0)
# 2. Join Gene ID linking it from transcript
df_prom <- merge(
  promoters,
  txdf,
  by.x = "tx_name",   
  by.y = "TXNAME",    
  all.x = TRUE        # keep all rows from promoters
)

# 3. Remove mitochondrial chromosome
df_prom_clean <- df_prom[!df_prom$seqnames == "NC_002333.2",]

# Covert to 0-based
df_prom_clean$start <- df_prom_clean$start - 1

# Add dummy column to match bed format
df_prom_clean$dum <- 0

# Reduce columns
df_red <- df_prom_clean[,c("seqnames", "start", "end", "GENEID", "dum", "strand")]

# Delete promoters that comprise the same region, strand and gene (transcripts with the same TSS)
df_red_unique <- df_red %>%
  distinct(seqnames, start, end, strand, GENEID, .keep_all = TRUE)

# Collapse GENEIDs that are the same in the rest of variables
df_collapsed <- df_red_unique %>%
  group_by(seqnames, start, end, strand, dum) %>%
  summarise(
    GENEID = paste(unique(GENEID), collapse = ";"),
    .groups = "drop"
  )

# Reorder columns
df_collapsed <- df_collapsed[,c("seqnames", "start", "end", "GENEID", "dum", "strand")]

# Write file
write.table(df_collapsed, file = out_file, sep = "\t", col.names = F, row.names = F, quote = F)
