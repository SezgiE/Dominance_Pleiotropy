library(arrow)
library(tidyverse)

rm(list = ls(all=TRUE)) # clean memory


snps_eqtl <- read_parquet("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_Analysis_v11_eQTL/Adipose_Subcutaneous.v11.eQTLs.signif_pairs.parquet")
genes_eqtl <- read.table(gzfile("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_Analysis_v11_eQTL/Adipose_Subcutaneous.v11.eGenes.txt.gz"), header=TRUE, sep="\t")
susie_eqtl <- read_parquet("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/GTEx_v11_Susie/Adipose_Subcutaneous.v11.eQTLs.SuSiE_summary.parquet")
snps_pleio <- read.table("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/b37_to_b38/pleio_snps_b38.bed")

pip_vals <- read_parquet("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/merged_pip_values.parquet")

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


asd <- read.table("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/res_all_pleio/gtex_susie_pleio_snps.tsv", sep="\t", header = T)

table(asd$biotype)


x_tissue <- susie_eqtl %>%
  group_by(phenotype_id) %>%
  summarise(indep_sig = n_distinct(cs_id, na.rm = TRUE))

sum(x_tissue$indep_sig)

x_tissue <- susie_eqtl %>%
  group_by(biotype, phenotype_id) %>%
  summarise(unique_cs = n_distinct(cs_id), .groups = "drop") %>%
  group_by(biotype) %>%
  summarise(total_cs = sum(unique_cs))

x2_tissue <- susie_eqtl %>%
  filter(phenotype_id == "ENSG00000292994.2")

length(unique(x_tissue$variant_id_b37))


susie_genes <- unique(susie_eqtl$phenotype_id)
susie_snps <- unique(susie_eqtl$variant_id)

all_snps <- unique(snps_eqtl$variant_id)


sig_genes <- genes_eqtl %>%
  filter(qval <= 0.05)

all_genes <- unique(sig_genes$gene_id)


all(susie_genes %in% all_genes)
diff <- setdiff(all_genes, susie_genes)


"""--- Resampling Test Results ---
  Parent Population Size: 12699515
Subset Size: 1891
Parent Metric (Median): 0.1709
Observed Subset Metric (Median): 0.1906
P-value (2-tailed): 0.3020
0.302"""