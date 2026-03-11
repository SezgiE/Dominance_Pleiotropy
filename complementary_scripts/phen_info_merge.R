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



unique_traits <- all_sig_df %>%
  mutate(sig_dom_traits = str_replace_all(sig_dom_traits, "[()\\s]", "")) %>% # Remove parentheses and spaces
  separate_rows(sig_dom_traits, sep = ",") %>%                               # Expand into individual rows
  filter(sig_dom_traits != "") %>%                                           # Remove empty strings if any
  distinct(sig_dom_traits) %>%
  rename("phen_code" = sig_dom_traits)


phen_info_Neale <- fread("/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/all_phenotypes_info.tsv.gz")
phen_info_UKB <- fread("/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/UKB_info.tsv", select = c("field_id", "main_category"))
phen_domains <- fread("/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/UKB_domains.tsv", select = c("category_id", "title"))

phen_info_filtered <- unique_traits %>%
  inner_join(phen_info_Neale, by = c("phen_code" = "phenotype"))

phen_info_filtered <- phen_info_filtered %>%
  mutate(field_id = str_extract(phen_code, "^[^_]+")) %>%
  mutate(field_id = as.character(field_id)) %>% 
  left_join(
    phen_info_UKB %>% mutate(field_id = as.character(field_id)), 
    by = "field_id"
  ) 

phen_info <- phen_info_filtered %>%
  left_join(
    phen_domains, by=c("main_category" = "category_id")
  ) %>%
  rename("category_code"= main_category,
         "category"= title)


write.xlsx(phen_info, file = "/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/phen_dict.xlsx", firstRow = TRUE, columnWidths = "auto")
