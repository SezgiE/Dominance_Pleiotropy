library(data.table)
library(tidyverse)
library(ggrepel)
library(viridis)
library(scales)
library(jsonlite)
library(openxlsx)

rm(list = setdiff(ls(all.names = TRUE), c("all_sig_df", "clumped_df")))

#  Load the data 
merged_path <- "/Users/sezgi/Downloads/FUMA_job727114/snps.txt"
coloc_merged <- fread(merged_path)


fuma <- fread("/Users/sezgi/Documents/dominance_pleiotropy/gene_level/fuma_input/fuma_input.txt")
unqvar <- unique(fuma$rsID)

missing_var <- fuma %>%
  filter(!(rsID %in% coloc_merged$rsID))
write.table(missing_var, "/Users/sezgi/Downloads/FUMA_job727114/missing_snps.txt", sep = "\t", quote = FALSE, row.names = FALSE)




snp_info <- coloc_merged %>%
  distinct(locus, .keep_all = TRUE)


canonical <- coloc_merged %>%
  filter(CANONICAL == "YES")

n_canon <- length(unique(canonical$variant))


non_coding <- coloc_merged %>%
  filter(Consequence == "intergenic_variant")

n_noncoding <- length(non_coding$variant)

variants <- c(unique(non_coding$variant), unique(canonical$variant))

missing_var <- coloc_merged %>%
  filter(!(variant %in% variants))

coloc_merged2 <- coloc_merged %>%
  select(-Beta, -SE)
write.table(coloc_merged2, "/Users/sezgi/Documents///fuma_input/fuma_input1.txt", sep = "\t", quote = FALSE, row.names = FALSE)
