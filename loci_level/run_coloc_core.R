#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(susieR))
suppressPackageStartupMessages(library(coloc))
rm(list =ls(all.names = TRUE))

# Estimate phenotypic variance
compute_yty <- function(beta, se, p, R, n, k) {
  beta_s <- beta * sqrt(2 * p * (1 - p))
  se_s <- se * sqrt(2 * p * (1 - p))
  
  XjtXj <- (n - 1) * diag(R)
  yty <- beta_s**2 * XjtXj + se_s**2 * XjtXj * (n - k)
  return(median(yty))
}


# Run SuSiE for coloc
run_susie_pipeline <- function(data_file, matrix_file, n, k_covariates, trait_type, prop_cases = NULL) {
  
  # 1. Load Data
  df <- fread(data_file, data.table = FALSE)
  R_matrix <- as.matrix(fread(matrix_file, header = FALSE))
  dimnames(R_matrix) <- NULL

  # 2. Format LD Matrix
  rownames(R_matrix) <- df$variant
  colnames(R_matrix) <- df$variant
  
  # 3. Compute Trait Variance Estimation
  yty <- compute_yty(beta = df$dominance_beta, 
                     se = df$dominance_se, 
                     p = df$minor_AF, 
                     R = R_matrix, 
                     n = n, 
                     k = k_covariates)
  
  var_y <- yty / (n - 1)
  sd_y <- sqrt(var_y)
  
  # 4. Construct coloc dataset
  dataset <- list(
    beta = df$dominance_beta,
    varbeta = df$dominance_se^2,
    N = n,                  
    MAF = df$minor_AF,
    type = trait_type,              
    snp = df$variant,
    position = df$pos,
    LD = R_matrix       
  )
  
  # Trait type "cc" for case-control, "quant" for quantitative
  if (trait_type == "cc") {
    dataset$s <- prop_cases
  } else if (trait_type == "quant") {
    dataset$sdY <- sd_y
  }
  
  # 5. Run SuSiE
  susie_res <- runsusie(d = dataset,
                        var_y = var_y, 
                        L = 10,
                        scaled_prior_variance = 0.1,
                        estimate_residual_variance = TRUE,
                        estimate_prior_variance = TRUE,
                        standardize = TRUE,
                        check_input = FALSE)
  
  return(susie_res)
}


# --- THE EXECUTION FLOW ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript run_coloc_core.R <stats_file1.tsv> <ld_matrix1.csv> <n_samples1> <k_samples1>
       <stats_file2.tsv> <ld_matrix2.csv> <n_samples2> <k_samples2> <out_file.tsv>")
}

# # Trait 1
data_file1 <- args[1]
ld_file1 <- args[2]
n_samples1 <- as.numeric(args[3])
k1_covariates <- as.numeric(args[4])

# Trait 2
data_file2 <- args[5]
ld_file2 <- args[6]
n_samples2 <- as.numeric(args[7])
k2_covariates <- as.numeric(args[8])

# tmp output
out_file <- args[9]

# # Run Trait 1
susie_trait1 <- run_susie_pipeline(
  matrix_file = ld_file1,
  data_file = data_file1,
  n = n_samples1,
  k_covariates = k1_covariates
)

# # Run Trait 2
susie_trait2 <- run_susie_pipeline(
  matrix_file = ld_file2,
  data_file = data_file2,
  n = n_samples2,
  k_covariates = k2_covariates
)


# susie_trait1 <- run_susie_pipeline(
#   data_file = "/Users/sezgi/Documents/dominance_pleiotropy/loci_level/susie_results/susie_raw_files/1747_1_16:89395438:90122562_data.csv",
#   matrix_file = '/Users/sezgi/Documents/dominance_pleiotropy/loci_level/susie_results/susie_raw_files/1747_1_16:89395438:90122562_matrix.csv',
#   n = 360000,
#   k_covariates = 13,
#   trait_type ="quant"
# )

# # Run Trait 2
# susie_trait2 <- run_susie_pipeline(
#   data_file = '/Users/sezgi/Documents/dominance_pleiotropy/loci_level/susie_results/susie_raw_files/1747_1_16:89395438:90122562_data.csv',
#   matrix_file = '/Users/sezgi/Documents/dominance_pleiotropy/loci_level/susie_results/susie_raw_files/1747_1_16:89395438:90122562_matrix.csv',
#   n = 360000,
#   k_covariates = 13,
#   trait_type ="quant"
# )


susie_coloc_res <- coloc.susie(susie_trait1, susie_trait2)

# Rows in the summary table have PP.H4 >= 0.80
sig_rows <- which(susie_coloc_res$summary$PP.H4.abf >= 0.80)

if (length(sig_rows) == 0) {
  
  cat("No pleiotropic signals met the >= 0.80 threshold at this locus.\n")
  
} else {
  
  res_df <- as.data.frame(susie_coloc_res$results)
  cols_to_keep <- c("snp", paste0("SNP.PP.H4.row", sig_rows))
  filtered_res_df <- res_df[, cols_to_keep, drop = FALSE]
  cs1_indices <- susie_coloc_res$summary$idx1[sig_rows]
  cs2_indices <- susie_coloc_res$summary$idx2[sig_rows]
  new_col_names <- paste0("cs", cs1_indices, "_cs", cs2_indices)
  colnames(filtered_res_df) <- c("snp", new_col_names)
  
  head(filtered_res_df)
  
  # Write the output
  #fwrite(filtered_res_df, out_file, sep = "\t", quote = FALSE)
}



