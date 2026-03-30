library(susieR)
rm(list =ls(all.names = TRUE))
set.seed(1)

mtrx <- read.csv('/Users/sezgi/Documents/dominance_pleiotropy/loci_level/susie_results/7:37985280:38151544_matrix.csv', header = F)

data <- read.csv("/Users/sezgi/Documents/dominance_pleiotropy/loci_level/susie_results/7:37985280:38151544_data.csv", sep = '\t')

z_scr <- data$add_z_score
mtrx <- as.matrix(mtrx)
mtrx2 <- (mtrx + t(mtrx)) / 2.0
dimnames(mtrx2) <- NULL
isSymmetric(mtrx2)

#mtrx_blend <- (0.999 * mtrx2) + (0.001 * diag(nrow(mtrx2)))
#diag(mtrx_blend) <- 1.000

fitted_rss3 <- susie_rss(z = z_scr, 
                         R = mtrx_blend, 
                         n = 361194,                          
                         L = 10,
                         estimate_prior_variance = FALSE,
                         max_iter = 500
                         )

fitted_rss3$sets


attr(mtrx2, "eigen") = eigen(mtrx2, symmetric = TRUE)
susie_plot(z_scr, y = "z", b=b)

lambda = estimate_s_rss(z_scr, mtrx2, n=361194)
lambda
condz_in = kriging_rss(z_scr, mtrx2, n=361194)
condz_in$plot



ev <- eigen(mtrx2, symmetric = TRUE, only.values = TRUE)$values

# 1. Check the minimum eigenvalue
min_ev <- min(ev)
print(paste("Minimum Eigenvalue:", min_ev))

# 2. Count how many are negative
neg_count <- sum(ev < 0)
print(paste("Number of negative eigenvalues:", neg_count))

# 3. Proportion of 'bad' variance
# If this ratio is very small (e.g., < 0.001), it's usually safe noise
bad_variance <- sum(abs(ev[ev < 0])) / sum(abs(ev))
print(paste("Proportion of negative variance:", bad_variance))





# Get indices of the top 5 values in row 155, sorted descending
top_5_indices <- order(mtrx[155, ], decreasing = TRUE)[1:15]

# Get the actual values
top_5_values <- mtrx[155, top_5_indices]
top_5_values

data(N3finemapping)
attach(N3finemapping)
n = nrow(X)

dim(Y)

Rin = cor(N3finemapping$X)


b <- true_coef[,1]
plot(b, pch=16, ylab='effect size')

which(b != 0)

ss = compute_suff_stat(N2finemapping$X, N2finemapping$Y[,1])


sumstats <- univariate_regression(X, Y[,1])
z_scores <- sumstats$betahat / sumstats$sebetahat
susie_plot(z_scores, y = "z", b=b)

R <- cor(X)


fitted_rss3 <- susie_rss(z = z_scr, R=mtrx, L = 10)


fitted_rss3$pip
fitted_rss3$sets

listss <- susie_get_cs(fitted_rss3)
listss$cs$L1
