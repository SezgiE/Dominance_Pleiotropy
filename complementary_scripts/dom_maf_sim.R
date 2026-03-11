# Define Adjustable Parameters
N   <- 35000  # Sample size
maf <- 0.45      # Minor Allele Frequency (p)
a   <- 0      # Additive effect
d   <- 0.0001      # Dominance deviation
h2  <- 0.05     # Narrow-sense heritability (to scale environmental noise)


p <- maf
q <- 1 - p

# Probabilities under Hardy-Weinberg Equilibrium
prob_0 <- q^2      # major homozygote (aa)
prob_1 <- 2 * p * q  # heterozygote (Aa)
prob_2 <- p^2      # minor homozygote (AA)

#  Simulate Genotypes (0, 1, 2)
# 0 = aa, 1 = Aa, 2 = AA
g <- sample(c(0, 1, 2), size = N, replace = TRUE, prob = c(prob_0, prob_1, prob_2))

#  Construct the Regression Columns
# Additive column (standard 0, 1, 2)
A_vec <- g

# Dominance column (Cockerham orthogonal coding: 0, 2p, 4p-2)
D_vec <- ifelse(g == 0, 0, 
                ifelse(g == 1, 2 * p, 4 * p - 2))

# Simulate the Phenotype (y)
# aa = 0, AA = 2a. The additive midpoint is 'a'.
# The dominance deviation 'd' pulls the heterozygote away from 'a'.
y_genetic <- ifelse(g == 0, 0,
                    ifelse(g == 1, a + d, 2 * a))

# Add environmental noise
var_g <- var(y_genetic)
var_e <- var_g * (1 - h2) / h2
y <- y_genetic + rnorm(N, mean = 0, sd = sqrt(var_e))

# Fit the Regression Model
add_model <- lm(y ~ A_vec)
dom_model <- lm(y ~ D_vec)


# 7. Results Comparison
cat("==========================================\n")
cat("SIMULATION PARAMETERS: a =", a, "| d =", d, "| p =", p, "\n")
cat("==========================================\n")
print(summary(add_model)$coefficients)
print(summary(dom_model)$coefficients)
cat("==========================================\n")
cat("Theoretical Expected Alpha (a + d(q-p)): ", a + d * (q - p), "\n")
cat("Theoretical Expected Beta_D (d):         ", d, "\n")