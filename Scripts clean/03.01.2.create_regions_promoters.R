# Create promoters regions (2kb upstream TSS) in BED format for later methylation level calculation. For CG islands, we annotated manually using UCSC genome browser.

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

# Reduce columns
df_red <- df_prom_clean[,c("seqnames", "start", "end", "tx_name")]

# Write file
write.table(df_red, file = out_file, sep = "\t", col.names = F, row.names = F, quote = F)