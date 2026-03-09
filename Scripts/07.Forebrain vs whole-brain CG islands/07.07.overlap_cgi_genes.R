# This script identifies genes that overlap with CpG islands (CGIs) in the zebrafish genome. It reads gene regions and CGI coordinates, converts them to GenomicRanges objects, finds overlaps, and outputs a list of genes that overlap with CGIs.

library(GenomicRanges) # version 1.62.1

# 1. Load input files
home <- path.expand("~")
file_gene_regions <- paste0(home, "/methylation_regions/regions/genes_plus_promoters.bed")
cgi_file <- paste0(home, "/methylation_regions/regions/cgi.bed")

# 2. Read data 
regions_ds <- read.table(file_gene_regions, header = F, col.names = c("chr", "start", "end", "gene_id", "dum", "strand"), sep = "\t", stringsAsFactors = FALSE)
cgi_ds <- read.table(cgi_file, header = F, col.names = c("chr", "start", "end", "name"), sep = "\t", stringsAsFactors = FALSE)

# Reduce dataset
cgi_ds <- cgi_ds[,c("chr", "start", "end", "name")]


# 3. Convert to GenomicRanges objects 
genes_gr <- GRanges(
  seqnames = regions_ds$chr,
  ranges = IRanges(start = regions_ds$start + 1, end = regions_ds$end), # +1 because GRanges interprets 1-based
  strand = regions_ds$strand,
  gene_id = regions_ds$gene_id
)

cgi_gr <- GRanges(
  seqnames = cgi_ds$chr,
  ranges = IRanges(start = cgi_ds$start + 1, end = cgi_ds$end), # +1 because GRanges interprets 1-based
  id = cgi_ds$name
)

# 4. Find overlaps 
hits <- findOverlaps(cgi_gr, genes_gr)

# Extract overlapping ranges
cgi_hits   <- cgi_gr[queryHits(hits)]
gene_hits  <- genes_gr[subjectHits(hits)]


overlap_df <- data.frame(
  cgi_chr   = as.character(seqnames(cgi_hits)),
  cgi_start = start(cgi_hits) - 1,  # back to BED (0-based)
  cgi_end   = end(cgi_hits),
  cgi_id    = mcols(cgi_hits)$id,
  
  gene_chr   = as.character(seqnames(gene_hits)),
  gene_start = start(gene_hits) - 1,
  gene_end   = end(gene_hits),
  gene_id    = mcols(gene_hits)$gene_id,
  strand     = as.character(strand(gene_hits)),
  
  stringsAsFactors = FALSE
)

out_file <- paste0(home, "/Chaterjee/overlap_genes_cgi.txt")
write.table(overlap_df, out_file, sep = "\t", quote = F, col.names = T, row.names = F)

# Reduced dataset
cgi_gene_collapsed <- aggregate(
  gene_id ~ cgi_id,
  data = overlap_df,
  FUN = function(x) paste(unique(x), collapse = ";")
)
out_file <- paste0(home, "/Chaterjee/overlap_genes_cgi_reduced.txt")
write.table(cgi_gene_collapsed, out_file, sep = "\t", quote = F, col.names = T, row.names = F)
