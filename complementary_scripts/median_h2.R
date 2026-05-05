library(data.table)
library(tidyverse)
library(ggrepel)
library(viridis)
library(scales)
library(jsonlite)
library(openxlsx)

rm(list = setdiff(ls(all.names = TRUE), c("all_sig_df", "clumped_df")))


df_h2 <- fread("/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/add_h2_UKB.tsv", select = c("phenotype", "h2_observed"))
df_traits <- fread("/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/supp_table1.tsv")

df_merged <- merge(df_traits, df_h2, by.x = "phenotype_code", by.y = "phenotype", all.x = TRUE, all.y = FALSE)

median(df_merged$h2_observed)


d_h2 <- fread("/Users/sezgi/Documents/dominance_heritability/dominance_h2_results.csv")
mean(d_h2$Dominance_h2)
