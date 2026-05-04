library(data.table)
library(tidyverse)
rm(list =ls(all.names = TRUE))

add <- readxl::read_xlsx("/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/a_sumStats.xlsx")
dom <- readxl::read_xlsx("/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/d_sumStats.xlsx")

phen_info <- fread("/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/all_phenotypes_info.tsv", sep = "\t", header = TRUE)
phen_info_UKB <- fread("/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/UKB_info.tsv", select = c("field_id", "main_category"))
phen_domains <- fread("/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/UKB_domains.tsv", select = c("category_id", "title"))


add <- add %>%
  select(phenotype_code, wget) %>%
  rename(additive_wget = wget)

dom <- dom %>%
  rename(phenotype_description = description)

merged_df <- merge(dom, add, by = "phenotype_code")

merged_df <- merge(merged_df, phen_info, by.x = "phenotype_code", by.y = "phenotype", all = FALSE)

merged_df <- merged_df %>%
  mutate(phen_id = sub("_.*", "", phenotype_code))

merged_df <- merge(merged_df, phen_info_UKB, by.x = "phen_id", by.y = "field_id", all.x = TRUE, all.y = FALSE)
merged_df <- merge(merged_df, phen_domains, by.x = "main_category", by.y = "category_id", all.x = TRUE, all.y = FALSE)


merged_df <- merged_df %>%
  mutate(
    clean_code = str_trim(toupper(phenotype_code)),
    first_letter = str_sub(clean_code, 1, 1),
    
    Clinical_Category = case_when(
      str_detect(clean_code, "^[IVX]+_") ~ paste0("ICD10 - ", phenotype_description),
      str_starts(clean_code, "ASTHMA|BRONCHITIS|COPD|PULMONARY") ~ "ICD10 - Diseases of the respiratory system",
      str_starts(clean_code, "CARDIAC") ~ "ICD10 - Diseases of the circulatory system",
      str_starts(clean_code, "KNEE|COX_ARTH") ~ "ICD10 - Diseases of the musculoskeletal system and connective tissue",
      str_starts(clean_code, "AB1_") ~ "ICD10 - Infectious and parasitic diseases",
      
      first_letter %in% c("A", "B") ~ "ICD10 - Infectious and parasitic diseases",
      first_letter == "C" | (first_letter == "D" & str_sub(clean_code, 2, 2) %in% as.character(0:4)) ~ "ICD10 - Neoplasms",
      (first_letter == "D" & str_sub(clean_code, 2, 2) %in% as.character(4:9)) ~ "ICD10 - Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism",
      first_letter == "E" ~ "ICD10 - Endocrine, nutritional and metabolic diseases",
      first_letter == "F" ~ "ICD10 - Mental and behavioural disorders",
      first_letter == "G" ~ "ICD10 - Diseases of the nervous system",
      (first_letter == "H" & str_sub(clean_code, 2, 2) %in% as.character(0:5)) ~ "ICD10 - Diseases of the eye and adnexa",
      (first_letter == "H" & str_sub(clean_code, 2, 2) %in% as.character(6:9)) ~ "ICD10 - Diseases of the ear and mastoid process",
      first_letter == "I" ~ "ICD10 - Diseases of the circulatory system",
      first_letter == "J" ~ "ICD10 - Diseases of the respiratory system",
      first_letter == "K" ~ "ICD10 - Diseases of the digestive system",
      first_letter == "L" ~ "ICD10 - Diseases of the skin and subcutaneous tissue",
      first_letter == "M" ~ "ICD10 - Diseases of the musculoskeletal system and connective tissue",
      first_letter == "N" ~ "ICD10 - Diseases of the genitourinary system",
      first_letter == "O" ~ "ICD10 - Pregnancy, childbirth and the puerperium",
      first_letter == "P" ~ "ICD10 - Certain conditions originating in the perinatal period",
      first_letter == "Q" ~ "ICD10 - Congenital malformations, deformations and chromosomal abnormalities",
      first_letter == "R" ~ "ICD10 - Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified",
      first_letter %in% c("S", "T") ~ "ICD10 - Injury, poisoning and certain other consequences of external causes",
      first_letter == "Z" ~ "ICD10 - Factors influencing health status and contact with health services",
      
      TRUE ~ "Uncategorized"
    )
  ) %>%
  select(-clean_code, -first_letter)

merged_df <- merged_df %>%
  mutate(title = if_else(source == "icd10" | source == "finngen", Clinical_Category, title))


final_df <- merged_df %>%
  select(phenotype_code, phenotype_description, title, main_category, sex, variable_type, source, n_non_missing, n_missing, n_controls, n_cases, wget, additive_wget) %>%
  rename(category = title)

length(unique(final_df$category))
category_counts <- as.data.frame(table(final_df$category))

write.table(final_df, "/Users/sezgi/Documents/dominance_pleiotropy/UKB_sumstats_Neale/supp_table1.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
