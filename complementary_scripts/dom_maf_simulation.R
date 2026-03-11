library(data.table)
library(tidyverse)

# 1. Set Parameters
n <- 1e6               # 1 Million individuals
p <- 0.3               # Minor allele frequency
q <- 1 - p
a <- 1.5               # Additive effect
d <- 1               # Dominance effect (d > a indicates overdominance)
intercept <- 10
heritability <- 0.3


genotypes <- rbinom(n, 2, p)

dt <- data.table(
  genotype = genotypes,
  x_a = genotypes - 1,               
  x_d = as.numeric(genotypes == 1)  
)

genetic_var <- var(dt$x_a * a + dt$x_d * d)
error_var <- genetic_var * (1 - heritability) / heritability

dt[, phenotype := intercept + (x_a * a) + (x_d * d) + rnorm(n, 0, sqrt(error_var))]

df <- dt %>%
  select(genotype, phenotype) %>%
  mutate(geno_dom_coded = case_when(
    genotype == 0 ~ 0,
    genotype == 1 ~ 2*p,
    genotype == 2 ~ 4*p - 2),
    geno_add_coded = genotype - (2 * p))


cov(df$genotype, df$geno_dom_coded)
cov(dt$x_a, dt$x_d)



model <- lm(phenotype ~ x_d, data = dt)
summary(model)


