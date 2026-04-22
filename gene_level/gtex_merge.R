library(arrow)
library(tidyverse)

rm(list = ls(all=TRUE)) # clean memory


snps_eqtl <- read_parquet("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_Analysis_v11_eQTL/Adipose_Subcutaneous.v11.eQTLs.signif_pairs.parquet")
genes_eqtl <- read.table(gzfile("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_Analysis_v11_eQTL/Adipose_Subcutaneous.v11.eGenes.txt.gz"), header=TRUE, sep="\t")
susie_eqtl <- read_parquet("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_v11_Susie/Adipose_Subcutaneous.v11.eQTLs.SuSiE_summary.parquet")
snps_pleio <- read.table("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/b37_to_b38/pleio_snps_b38.bed")

# Adjust SNP_b38 list
colnames(snps_pleio) <- c("chr", "start", "pos", "variant_id_b37")
snps_pleio <- snps_pleio %>%
  dplyr::select(chr, pos, variant_id_b37)

# Adjust Susie eqtl
susie_eqtl <-susie_eqtl %>% 
  separate(variant_id, into = c("chr", "pos", "ref", "alt", "build"), sep = "_")

# Check case-sensitivity and whitespace
snps_pleio$chr <- tolower(trimws(as.character(snps_pleio$chr)))
snps_pleio$pos <- trimws(as.character(snps_pleio$pos))

susie_eqtl$chr <- tolower(trimws(as.character(susie_eqtl$chr)))
susie_eqtl$pos <- trimws(as.character(susie_eqtl$pos))

# Overlapping causal variants
snps_overlapped <- susie_eqtl %>% 
  inner_join(snps_pleio, by = c("chr", "pos"))
colnames(snps_overlapped)[1] <- "gene_id"

# Obtain information from gene data
genes_eqtl_filt <- genes_eqtl %>%
  dplyr::select(gene_id, gene_start, gene_end, strand, beta_shape1, beta_shape2, qval)
  
snps_overlapped <- snps_overlapped %>%
  inner_join(genes_eqtl_filt, by = c("gene_id"))

# Filter based on Gene qvals
results_df <- snps_overlapped %>%
  filter(qval <= 0.05)


asd <- read.table("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/merged_gtex_susie_snps.tsv", sep="\t", header = T)

length(unique(asd$variant_id_b37))




