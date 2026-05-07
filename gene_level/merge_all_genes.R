library(data.table)
library(tidyverse)
library(ggrepel)
library(viridis)
library(scales)
library(jsonlite)
library(openxlsx)
library(biomaRt)

rm(list = ls(all.names = TRUE))

# MAGMA Reader
magma_reader <- function(filepath) {
  
  raw_lines <- readLines(filepath)
  data_lines <- raw_lines[!grepl("^#", raw_lines)]
  

  long_df <- tibble(line = data_lines) %>%
    mutate(split_data = str_split(line, "\\s+")) %>%
    mutate(
      gene_id = map_chr(split_data, 1),
      region  = map_chr(split_data, 2),
      variant = map(split_data, ~ .x[-(1:2)])
    ) %>%
    dplyr::select(gene_id, variant) %>%
    unnest(variant) %>%
    filter(variant != "")
  
  return(long_df)
  
}


# SNPs
coloc_snps <- fread("/Users/sezgi/Documents/dominance_pleiotropy/loci_level/coloc_results/coloc_snp_info.tsv", select = c("phen_code", "phen_name",
                                                                                                                          "category", "variant", "rsID"))


# GTEx
gtex_genes <- fread("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/gtex/gtex_res/gtex_susie_pleio_snps.tsv", select = c("tissue_name", "ENSG",
                                                                                                                                "gene_name", "variant_id_b37", "biotype",
                                                                                                                                "qval", "cs_id", "cs_size", "pip")) 

                                                                                                                                
# MAGMA
magma_genes <- magma_reader("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/magma/magma_v1/magma_pleio_mapping.genes.annot")
manual_genes <- magma_reader("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/magma/magma_v1/manuel_pleio_mapping.genes.annot")
manual_genes_dist <- fread("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/magma/magma_v1/manuel_pleio_mapping.genes_log.txt",
                           select = c("SNP", "MIN_DISTANCE"))
pos_mapped_genes <- rbind(magma_genes, manual_genes)

pos_mapped_genes <- merge(pos_mapped_genes, manual_genes_dist, by.x = "variant", by.y = "SNP", all.x = TRUE, all.y = TRUE)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_annotations <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'),
  filters = 'ensembl_gene_id',
  values = unique(pos_mapped_genes$gene_id),
  mart = ensembl
)

pos_mapped_genes <- merge(
  pos_mapped_genes, 
  gene_annotations, 
  by.x = "gene_id", 
  by.y = "ensembl_gene_id", 
  all.x = TRUE
)


# Merging
df_merged <- merge(coloc_snps, pos_mapped_genes, by.x = "variant", by.y = "variant", all.x = TRUE, all.y = TRUE, allow.cartesian = TRUE)

df_merged <- df_merged %>%
  rename(gene_id_pos = gene_id,
         gene_name_pos = external_gene_name,
         gene_biotype_pos = gene_biotype,
         nearest_gene_dist = MIN_DISTANCE)


df_merged_all <- merge(df_merged, gtex_genes, by.x = "variant", by.y = "variant_id_b37", all.x = TRUE, all.y = TRUE, allow.cartesian = TRUE)

df_merged_all <- df_merged_all %>%
  rename(gene_id_eqtl = ENSG,
         gene_name_eqtl = gene_name,
         tissue_name_eqtl = tissue_name,
         gene_biotype_eqtl = biotype,
         gene_qval_eqtl = qval,
         cs_id_eqtl = cs_id, 
         cs_size_eqtl = cs_size, 
         pip_eqtl = pip) %>%
  dplyr::select(variant, rsID, phen_code, phen_name, category, gene_id_pos, gene_name_pos, gene_biotype_pos, nearest_gene_dist,
                gene_id_eqtl, gene_name_eqtl, gene_biotype_eqtl, tissue_name_eqtl, gene_qval_eqtl, cs_id_eqtl,
                cs_size_eqtl, pip_eqtl)



fwrite(df_merged_all, "/Users/sezgi/Documents/dominance_pleiotropy/gene_level/genes_all/genes_all.tsv", sep = "\t")
