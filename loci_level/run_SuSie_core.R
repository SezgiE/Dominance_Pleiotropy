#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(susieR))

# Implemented based on github.com/mkanai/finemapping-pipeline/blob/master/R/run_susieR.R
# --- THE FUNCTIONS  ---
compute_yty <- function(beta, se, p, R, n, k) {
  beta_s <- beta * sqrt(2 * p * (1 - p))
  se_s <- se * sqrt(2 * p * (1 - p))
  
  XjtXj <- (n - 1) * diag(R)
  yty <- beta_s**2 * XjtXj + se_s**2 * XjtXj * (n - k)
  return(median(yty))
}


summarize.susie.cs <- function(object, orig_vars, R, low_purity_threshold = 0.5) {
  if (is.null(object$sets)) return(list(vars = orig_vars, cs = NULL))
  
  variables <- data.frame(cbind(1:length(object$pip), object$pip, -1, NA, NA, NA, NA))
  colnames(variables) <- c("variable", "variable_prob", "cs", "cs_specific_prob", "low_purity", "lead_r2", "lambda")
  added_vars <- c()
  
  if (object$null_index > 0) variables <- variables[-object$null_index, ]
  
  if (!is.null(object$sets$cs)) {
    for (i in 1:length(object$sets$cs)) {
      if (any(object$sets$cs[[i]] %in% added_vars)) next
      added_vars <- append(added_vars, object$sets$cs[[i]])
      
      in_cs_idx <- which(variables$variable %in% object$sets$cs[[i]])
      variables$cs[in_cs_idx] <- object$sets$cs_index[[i]]
      variables[in_cs_idx, "cs_specific_prob"] <- object$alpha[object$sets$cs_index[[i]], object$sets$cs[[i]]]
      variables$low_purity[in_cs_idx] <- object$sets$purity$min.abs.corr[i] < low_purity_threshold
      
      lead_pip_idx <- in_cs_idx[which.max(variables$variable_prob[in_cs_idx])]
      variables$lead_r2[in_cs_idx] <- R[lead_pip_idx, in_cs_idx]^2
    }
  }
  return(list(vars = variables))
}


# --- THE EXECUTION FLOW ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript run_susie_cor.R <stats_file.tsv> <ld_matrix.csv> <n_samples> <out_file.tsv>")
}

stats_file <- args[1]
ld_file <- args[2]
n_samples <- as.numeric(args[3])
out_file <- args[4]
k_covariates <- 25

# Load data
df <- fread(stats_file, data.table = FALSE)
R_matrix <- as.matrix(fread(ld_file, header = FALSE))
dimnames(R_matrix) <- NULL

# Calculate phenotypic variance estimate
yty <- compute_yty(beta = df$dominance_beta, se = df$dominance_se, p = df$minor_AF, R = R_matrix, n = n_samples, k = k_covariates)
var_y <- yty / (n_samples - 1)

# Run SuSiE 
fitted_bhat <- tryCatch({
  susie_rss(
    bhat = df$dominance_beta,
    shat = df$dominance_se,
    R = R_matrix,
    n = n_samples,
    var_y = var_y,
    L = 10,
    scaled_prior_variance = 0.1,
    estimate_residual_variance = TRUE,
    estimate_prior_variance = TRUE,
    standardize = TRUE,
    check_input = FALSE
  )
}, error = function(e) {
  message("SuSiE execution failed: ", e$message)
  q(status = 1)
})


# Estimate the lambda scalar for the locus
lambda_val <- estimate_s_rss(df$dom_z_score, R_matrix, n = n_samples)

# Extract data 
cs_summary <- summarize.susie.cs(fitted_bhat, df, R_matrix)

# Merge back into the original dataframe
df$PIP <- cs_summary$vars$variable_prob
df$CS <- cs_summary$vars$cs
df$CS_prob <- cs_summary$vars$cs_specific_prob
df$low_purity <- as.logical(cs_summary$vars$low_purity)
df$lead_r2 <- cs_summary$vars$lead_r2
df$lambda <- lambda_val

# Calculate posterior mean and standard deviation
df$post_mean <- susie_get_posterior_mean(fitted_bhat)
df$post_sd <- susie_get_posterior_sd(fitted_bhat)

# Clean up CS assignments (convert -1 to NA or 0)
df$CS[df$CS == -1] <- 0

# Save the final output
fwrite(df, out_file, sep = "\t", quote = FALSE)


