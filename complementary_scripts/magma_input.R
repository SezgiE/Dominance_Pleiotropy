library(tidyverse)
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install("biomaRt")
library(biomaRt)

rm(list = ls(all=TRUE)) # clean memory


# SNP positions b37
snps_pleio <- read.delim("/Users/sezgi/Documents/dominance_pleiotropy/loci_level/coloc_results/coloc_snps.tsv")

snps_pleio <- snps_pleio %>%
  separate(variant, into = c("CHR", "BP", "REF", "ALT"), sep = ":", remove = FALSE, convert = TRUE) %>%
  select(variant, CHR, BP)

write.table(snps_pleio, "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/magma_v1/snp_loc.bim", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# SNP positions b38
snps_pleio_b38 <- read.table("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/b37_to_b38/pleio_snps_b38.bed", stringsAsFactors = FALSE)

snps_pleio_b38$V1 <- gsub("chr", "", snps_pleio_b38$V1)
magma_snps <- snps_pleio_b38[, c("V4", "V1", "V3")]
magma_snps <- na.omit(magma_snps)

write.table(magma_snps, "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/magma_v1/snp_loc38.txt", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# ---------------------  GENELOC_FILE from ensembl ----------------------------
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve every gene 
all_genes_raw <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype", "chromosome_name", "start_position", "end_position", "strand"),
  mart = ensembl
)

write.table(all_genes_raw, file = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/magma_v1/ensembl_genes_raw38.tsv", 
            sep = "\t", row.names = F, quote = F)

# Format for MAGMA requirements
magma_loc <- all_genes_raw[all_genes_raw$chromosome_name %in% as.character(1:22), ]
magma_loc$chromosome_name <- as.numeric(magma_loc$chromosome_name)

# Recode Strand to + and -
magma_loc$strand <- ifelse(magma_loc$strand == 1, "+", "-")

# order: ID, Chrom, Start, Stop, Strand
magma_loc <- magma_loc[, c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand")]

# Export without header or quotes
write.table(magma_loc, file = "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/magma_v1/genes_loc38.loc", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

