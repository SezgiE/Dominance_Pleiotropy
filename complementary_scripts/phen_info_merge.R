library(data.table)
library(tidyverse)
library(ggrepel)
library(viridis)
library(scales)
library(openxlsx)

rm(list = setdiff(ls(all.names = TRUE), c("all_sig_df")))

#  Load the data 
all_sig_path <- "/Users/sezgi/Documents/dominance_pleiotropy/SNP_level/significant_SNPs/all_sig_SNPs.tsv.gz"
all_sig_df <- fread(all_sig_path, select = c("variant", "rsid", "chr", "pos", "add_sig_total", "dom_sig_total",
                                             "over_dom_sig_total", "sig_dom_traits", "sig_over_dom_traits"),
                    colClasses = c(
                      variant = "character",
                      rsid    = "character",
                      chr     = "numeric",
                      pos     = "numeric",
                      add_sig_total = "numeric",
                      dom_sig_total = "numeric"
                    ))

dom_pleio <- all_sig_df %>%
  filter(dom_sig_total > 1)


unique_traits <- dom_pleio %>%
  mutate(sig_dom_traits = str_replace_all(sig_dom_traits, "[()\\s]", "")) %>% # Remove parentheses and spaces
  separate_rows(sig_dom_traits, sep = ",") %>%                               # Expand into individual rows
  filter(sig_dom_traits != "") %>%                                           # Remove empty strings if any
  distinct(sig_dom_traits) %>%
  rename("phenotype_code" = sig_dom_traits)


phen_info_Neale <- fread("/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/all_phenotypes_info.tsv")
phen_domains <- fread("/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/supp_table1.tsv", select = c("phenotype_code", "category"))

phen_info_filtered <- unique_traits %>%
  inner_join(phen_info_Neale, by = c("phenotype_code" = "phenotype"))

phen_info_filtered <- merge(phen_info_filtered, phen_domains, by.x = "phenotype_code", by.y = "phenotype_code", all.x = TRUE, all.y = FALSE)


write.xlsx(phen_info_filtered, file = "/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/phen_dict.xlsx", firstRow = TRUE, columnWidths = "auto")
