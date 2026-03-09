library(tidyverse)
# Population Parameters
set.seed(42)
N <- 1000000       # 1 million individuals
p1 <- 0.2          # MAF for SNP 1 (20%)
p2 <- 0.3          # MAF for SNP 2 (30%)
D <- 0.04          # Linkage Disequilibrium coeff (D) between SNP 1 and SNP 2 

#  Haplotype Frequencies
# formula: freq = expected frequency +/- D
freq_11 <- (p1 * p2) + D            # Both minor alleles
freq_10 <- (p1 * (1 - p2)) - D      # SNP 1 minor, SNP 2 major
freq_01 <- ((1 - p1) * p2) - D      # SNP 1 major, SNP 2 minor
freq_00 <- ((1 - p1) * (1 - p2)) + D  # Both major alleles

hap_probs <- c(freq_11, freq_10, freq_01, freq_00)

# Simulating Independent Chromosomes
simulate_chromosome <- function(N, probs) {
  # Randomly draw a haplotype based on the calculated frequencies
  draws <- sample(1:4, size = N, replace = TRUE, prob = probs)
  snp1_raw <- ifelse(draws %in% c(1, 2), 1, 0)
  snp2_raw <- ifelse(draws %in% c(1, 3), 1, 0)
  
  return(data.frame(snp1_raw, snp2_raw))
}

# Generating the Diploid Population
# mom and dad completely independent (HWE Assumption)
maternal <- simulate_chromosome(N, hap_probs)
paternal <- simulate_chromosome(N, hap_probs)

# Mean-Centering the Alleles
m1 <- maternal$snp1_raw - p1
m2 <- maternal$snp2_raw - p2

f1 <- paternal$snp1_raw - p1
f2 <- paternal$snp2_raw - p2

# Building the Additive and Dominance Genotypes
A1 <- m1 + f1
A2 <- m2 + f2

# Dominance = Interaction (product) of the chromosomes
d1 <- 2 * m1 * f1
d2 <- 2 * m2 * f2


# dataframe with all the required variables
geno_df <- data.frame(A1, A2, d1, d2) %>%
  rename(SNP1_centered = A1,
         SNP2_centered = A2,
         SNP1_dom_centered = d1,
         SNP2_dom_centered = d2)%>%
  mutate(SNP1 = SNP1_centered + 2*p1,
         SNP2 = SNP2_centered + 2*p2,
         SNP1_dom = case_when(  # hard-coding SNP1 dominance based on 0, 2*MAF, 4*MAF-2
           SNP1 == 0 ~ 0,
           SNP1 == 1 ~ 2*p1,
           SNP1 == 2 ~ 4*p1 - 2),
         SNP2_dom = case_when( # hard-coding SNP2 dominance based on 0, 2*MAF, 4*MAF-2
           SNP2 == 0 ~ 0,
           SNP2 == 1 ~ 2*p2,
           SNP2 == 2 ~ 4*p2 - 2)
  )


# Calculate and Compare Correlations
rD_hard_coded <- cor(geno_df$SNP1_dom, geno_df$SNP2_dom)

r_D_centered <- cor(geno_df$SNP1_dom_centered, geno_df$SNP2_dom_centered)
r_A_centered <- cor(geno_df$SNP1_centered, geno_df$SNP2_centered)
r_A_squared <- r_A_centered^2



# the Results
cat("==========================================\n")
cat("SIMULATION RESULTS\n")
cat("==========================================\n")
cat("Additive Correlation (r_A): ", round(r_A_centered, 6), "\n")
cat("Dominance Correlation (r_D): ", round(r_D_centered, 6), "\n")
cat("Dominance Correlation (r_D) coded as 0, 2*MAF, 4*MAF-2: ", round(rD_hard_coded, 6), "\n")
cat("Squared Additive Correlation (r_A^2): ", round(r_A_squared, 6), "\n")
cat("Difference: r_A^2 - r_D                          ", round(abs(r_A_squared - r_D), 6), "\n")
