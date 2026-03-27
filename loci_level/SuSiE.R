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

mtrx_blend <- (0.999 * mtrx2) + (0.001 * diag(nrow(mtrx2)))
diag(mtrx_blend) <- 1.000

fitted_rss3 <- susie_rss(z = z_scr, 
                         R = mtrx_blend, 
                         n = 361194,                          
                         L = 10,
                         estimate_prior_variance = FALSE,
                         max_iter = 500
                         )

fitted_rss3$sets


# Get indices of the top 5 values in row 155, sorted descending
top_5_indices <- order(mtrx[155, ], decreasing = TRUE)[1:15]

# Get the actual values
top_5_values <- mtrx[155, top_5_indices]
top_5_values

data(N3finemapping)
attach(N3finemapping)
n = nrow(X)

dim(Y)

b <- true_coef[,1]
plot(b, pch=16, ylab='effect size')

which(b != 0)


sumstats <- univariate_regression(X, Y[,1])
z_scores <- sumstats$betahat / sumstats$sebetahat
susie_plot(z_scores, y = "z", b=b)

R <- cor(X)


fitted_rss3 <- susie_rss(z = z_scr, R=mtrx, L = 10)


fitted_rss3$pip
fitted_rss3$sets

listss <- susie_get_cs(fitted_rss3)
listss$cs$L1
