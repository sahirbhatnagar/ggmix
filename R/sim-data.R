
# # genotypes
# N.MAF <- 1000
# MAF <- runif(N.MAF, 0.05, 0.25)
# x <- as.matrix(sapply(MAF, fun <- function(x) rbinom(600, 2, x)))
# dim(x)
# x[1:5, 1:5]
# # link to a 600 x 600 theoretical kinship matrix https://www.dropbox.com/s/bk8zwqtzgvvm89z/kin1.Rdata?dl=0
# # load kinship matrix from dropbox link
# load("~/Dropbox/PhD/Year4/penfam/data/kin1.Rdata")
# Phi <- 2 * kin1
# Phi[1:10, 1:10]
# # random effect
# s.g = 6 # variance of the polygenic random effect
# P = mvrnorm(1, rep(0, 600), s.g*Phi)
# s.e = 4 # residual-effect variance
# b <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, rep(0,980))
# #
# X <- mvrnorm(600, rep(1,1000), diag(1,1000,1000))
# dim(X)
# X[1:5,1:5]
# Y <- 3 + X%*%b + P + rnorm(600,0,s.e)
#
#
# diag(s.g*Phi)



library(magrittr)
library(MASS)
load("~/Dropbox/PhD/Year4/penfam/data/kin1.Rdata")
Phi <- 2 * kin1
Phi[1:5,1:5]
Phi_eigen <- eigen(Phi)
Lambda <- Phi_eigen$values
any(Lambda < 1e-3)
all(Lambda > 0)
rcond(Phi)
kappa(Phi)


# simulation parameters
eta <- 0.35
sigma2 <- 2
p <- 1000
b0 <- 3
b <- c(runif(10, 0.8,1.2), rep(0,980), runif(10, -1.2, -0.8))

n <- nrow(Phi)

# polygenic random effect
P <- mvrnorm(1, mu = rep(0, n), Sigma = eta * sigma2 * Phi)

# environment random effect
E <- mvrnorm(1, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))

X <- mvrnorm(n, rep(1,p), diag(p))

Y <- b0 + X %*% b + P + E




