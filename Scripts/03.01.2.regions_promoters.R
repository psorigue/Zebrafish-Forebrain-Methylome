library(GenomicFeatures)
library(dplyr)

# build a TxDb from GTF
file_annot <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Epi/Ref_genome/GRCz12tu/genome/genomic.gtf"
txdb <- txdbmaker::makeTxDbFromGFF(file_annot)
# Create transcript -> gene link
txdf <- AnnotationDbi::select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")

# Gene annotation info
gene_ranges <- read.csv("//files1.igc.gulbenkian.pt/folders/ANB/Pol/Epi/Ref_genome/GRCz12tu/danrer_gene_info.tsv",
                        sep = "\t")


# Make promoter dataset
# 1. Define regions into a GRanges object
promoters <- promoters(txdb, upstream=2000, downstream=0)
# 2. Join Gene ID linking it from transcript
df_prom <- merge(
  promoters,
  txdf,
  by.x = "tx_name",   # column in your BED dataframe
  by.y = "TXNAME",    # column in txdf
  all.x = TRUE        # keep all rows from df_red
)
# 3. Remove mitochondrial chromosome
df_prom_clean <- df_prom[!df_prom$seqnames == "NC_002333.2",]

# Covert to 0-based
df_prom_clean$start <- df_prom_clean$start - 1

# Reduce columns
df_red <- df_prom_clean[,c("seqnames", "start", "end", "tx_name")]

# Write file
file_name <- "//files1.igc.gulbenkian.pt/folders/ANB/Pol/Methylome/methylation_regions/regions/promoters.bed"
write.table(df_red, file = file_name, sep = "\t", col.names = F, row.names = F, quote = F)
