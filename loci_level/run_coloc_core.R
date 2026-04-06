#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(susieR))
suppressPackageStartupMessages(library(coloc))
rm(list =ls(all.names = TRUE))
options(warn = -1)

# Estimate phenotypic variance
compute_yty <- function(beta, se, p, R, n, k) {
  beta_s <- beta * sqrt(2 * p * (1 - p))
  se_s <- se * sqrt(2 * p * (1 - p))
  
  XjtXj <- (n - 1) * diag(R)
  yty <- beta_s**2 * XjtXj + se_s**2 * XjtXj * (n - k)
  return(median(yty))
}


# Run SuSiE for coloc
run_susie_pipeline <- function(data_file, matrix_file, n, k_covariates, prop_cases, trait_type) {
  
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


# Format and save significant colocalization results
format_coloc_results <- function(susie_coloc_res, out_file, h4_threshold = 0.50) {
  
  # Find significant rows
  sig_rows <- which(susie_coloc_res$summary$PP.H4.abf >= h4_threshold)
  
  if (length(sig_rows) == 0) {
    cat(sprintf("No pleiotropic signals met the >= %.2f threshold at this locus. Skipping...\n", h4_threshold))
    return(FALSE) # Signal to the main script that nothing was found
  } 
  
  # Extract base results
  res_df <- as.data.frame(susie_coloc_res$results)
  
  cs1_indices <- susie_coloc_res$summary$idx1[sig_rows]
  cs2_indices <- susie_coloc_res$summary$idx2[sig_rows]
  
  # Initialize an empty list to store the long-format blocks
  long_format_list <- list()
  
  # Loop through each significant signal and create a stacked block
  for (i in seq_along(sig_rows)) {
    
    # Identify the specific column holding the PIP for this signal
    row_num <- sig_rows[i]
    pip_col_name <- paste0("SNP.PP.H4.row", row_num)

    # Check if the expected PIP column exists; if not, try alternative naming conventions 
    # as coloc drops the .row1 suffix if one credible set for both traits
    if (!(pip_col_name %in% names(res_df))) {
      if ("SNP.PP.H4" %in% names(res_df)) {
        pip_col_name <- "SNP.PP.H4" # Use the suffix-less name
      } else {
        # Absolute fallback: find whatever column contains "PP.H4"
        pip_col_name <- grep("PP\\.H4", names(res_df), value = TRUE)[1]
      }
    }
    
    # Construct the signal label
    cs_label <- paste0("cs", cs1_indices[i], "_cs", cs2_indices[i])
    
    # Extract the global H4 for this specific signal
    h4_val <- susie_coloc_res$summary$PP.H4.abf[row_num]
    
    # Create the temporary long dataframe for this specific signal
    temp_df <- data.frame(
      variant = res_df$snp,
      cs = cs_label,
      cs_H4 = h4_val,
      PIP = res_df[[pip_col_name]]
    )
    
    # Add to our list
    long_format_list[[i]] <- temp_df
  }
  
  final_long_df <- do.call(rbind, long_format_list)
  final_long_df <- final_long_df[final_long_df$PIP > 0, ]
  final_long_df <- final_long_df[order(final_long_df$cs, -final_long_df$PIP), ]
  fwrite(final_long_df, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  return(TRUE) # Signal success
}


# ---------------------- THE EXECUTION FLOW -------------------------#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 13) {
  stop("Usage: Rscript run_coloc_core.R <stats_file1.tsv> <ld_matrix1.csv> <n_samples1> <k_samples1>
       <stats_file2.tsv> <ld_matrix2.csv> <n_samples2> <k_samples2> <out_file.tsv>")
}

# ---------------------- Read Data -------------------------#
# # Trait 1
data_file1 <- args[1]
ld_file1 <- args[2]
n_samples1 <- as.numeric(args[3])
k1_covariates <- as.numeric(args[4])
prop_cases1 <- as.numeric(args[5])
type1 <- args[6]

# Trait 2
data_file2 <- args[7]
ld_file2 <- args[8]
n_samples2 <- as.numeric(args[9])
k2_covariates <- as.numeric(args[10])
prop_cases2 <- as.numeric(args[11])
type2 <- args[12]

# tmp output
out_file <- args[13]


# ---------------------- Run Trait 1 -------------------------#
susie_trait1 <- tryCatch({
  run_susie_pipeline(
    data_file = data_file1,
    matrix_file = ld_file1,
    n = n_samples1,
    k_covariates = k1_covariates,
    prop_cases = prop_cases1,
    trait_type = type1
  )
}, error = function(e) {
  cat("\n[SUSIE ERROR] Trait 1 fine-mapping failed.\n", file=stderr())
  cat("R Error:", conditionMessage(e), "\n", file=stderr())
  return(NULL)
})


# ---------------------- Run Trait 2 -------------------------#
susie_trait2 <- tryCatch({
  run_susie_pipeline(
    matrix_file = ld_file2,
    data_file = data_file2,
    n = n_samples2,
    k_covariates = k2_covariates,
    prop_cases = prop_cases2,
    trait_type = type2
  )
}, error = function(e) {
  cat("\n[SUSIE ERROR] Trait 2 fine-mapping failed.\n", file=stderr())
  cat("R Error:", conditionMessage(e), "\n", file=stderr())
  return(NULL)
})


if (is.null(susie_trait1) || is.null(susie_trait2)) {
  quit(save = "no", status = 1)
}


# ----------------------- Run coloc --------------------------#
susie_coloc_res <- tryCatch({
  
  coloc.susie(susie_trait1, susie_trait2)
  
}, error = function(e) {
  
  # Print the error to the console for your logs, but do NOT halt execution
  cat("\n[COLOC ERROR] Coloc failed for this pair (likely zero overlapping SNPs).\n", file=stderr())
  cat("R Error:", conditionMessage(e), "\n", file=stderr())
  
  # Return NULL to signal the failure to the rest of the script
  return(NULL)
  
})


if (is.null(susie_coloc_res)) {
  cat("Skipping...\n")
  quit(save = "no", status = 1)
}


# ---------------------- Save Results -------------------------#
found_signals <- format_coloc_results(susie_coloc_res, out_file, h4_threshold = 0.50)

# If no signals met the threshold, exit cleanly with status 0
if (!found_signals) {
  cat("No pleiotropic signals met the threshold value at this locus. Skipping...\n")
  quit(save = "no", status = 0)
}


