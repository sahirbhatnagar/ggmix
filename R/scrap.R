# options(warnPartialMatchArgs = FALSE, warnPartialMatchDollar = TRUE, warnPartialMatchAttr = TRUE)
library(magrittr)
library(glmnet)
library(MASS)

rm(list=ls())
source("~/git_repositories/penfam/R/fitting.R")
source("~/git_repositories/penfam/R/functions.R")
source("~/git_repositories/penfam/R/methods.R")
source("~/git_repositories/penfam/R/plot.R")
source("~/git_repositories/penfam/R/sim-data.R")

x <- X
y <- Y
phi <- Phi
lambda_min_ratio <- 0.001
nlambda <- 100
#convergence criterion
epsilon <- 1e-7
maxit <- 1e6
an = log(log(600)) * log(600)
lambda <- 0.10
#======================================

np <- dim(x)
n <- np[[1]]
p <- np[[2]]

# add column of 1s to x for intercept
x <- cbind(rep(1, n), x)
# x[1:5, 1:5]

phi_eigen <- eigen(phi)
# this is a N_T x N_T matrix
U <- phi_eigen$vectors

# dim(U)
# vector of length N_T
Lambda <- phi_eigen$values

utx <- crossprod(U, x)
uty <- crossprod(U, y)

# get sequence of tuning parameters
lamb <- lambda_sequence(x = utx, y = uty, phi = phi, lambda_min_ratio = lambda_min_ratio)

res <- penfam(x = X, y = Y, phi = Phi, lambda_min_ratio = 0.001,
              nlambda = 100,
              an = log(log(n)) * log(n))


#' Check of KKT
#' @param x should be U^T X, where U is the matrix of eigenvectors and X
#'   contains the first column of ones for the intercept. x should be a mtrix of
#'   dimension n x (p+1)
#' @param beta should include intercept as well. A vector of length p+1.
#' @param tol.beta Tolerance for determining if a coefficient is zero
#' @param tol.kkt Tolerance for determining if an entry of the subgradient is
#'   zero

kkt_check <- function(eta, sigma2, beta, eigenvalues, x, y, nt,
                      lambda, tol.beta = 1e-5, tol.kkt = 0.1){

  # x = utx
  # y = uty
  # tol.beta = 1e-5 ; tol.kkt = 0.1
  # eigenvalues = Lambda
  # nt = dim(x)[[1]]
  # p = dim(x)[[2]]-1
  # eta = lamb$eta
  # sigma2 = lamb$sigma2
  # beta = c(lamb$beta0, rep(0, p))
  # =======================

  di <- 1 + eta * (eigenvalues - 1)
  wi <- (1 / sigma2) * diag(1 / di)

  # KKT for beta0
  kkt_beta0 <- grr_beta0(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt)

  # KKT for eta
  kkt_eta <- grr_eta(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt)

  # KKT for sigma2
  kkt_sigma2 <- grr_sigma2(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt)

  # KKT for beta
  g0 <- crossprod(x[,-1, drop = F], wi) %*% (y - x %*% beta)

  g <- g0 - lambda * sign(beta)

  gg <- g0/lambda

  # which of the betas are non-zero, subject to the tolerance level for beta
  oo <- abs(bhatt) > tol.beta

  # cat(
  #   c(
  #     max(abs(g[oo])) > tol.kkt,
  #     min(gg[!oo]) < -1 - tol.kkt,
  #     max(gg[!oo]) > 1 + tol.kkt
  #   ),
  #   fill = T)

  kkt_beta <- c(max(abs(g[oo])) > tol.kkt,
                min(gg[!oo]) < -1 - tol.kkt,
                max(gg[!oo]) > 1 + tol.kkt)

  return(list(kkt_beta0 = kkt_beta0, kkt_eta = kkt_eta,
              kkt_sigma2 = kkt_sigma2, kkt_beta = kkt_beta))


}

















# Simulate Data -----------------------------------------------------------


rm(list = ls())
source("~/git_repositories/penfam/R/functions.R")
source("~/git_repositories/penfam/R/sim-data.R")



# lambda <- lambda_sequence(x = X, y = Y, phi = Phi, lambda_min_ratio = 1e-4)

res <- penfam(x = X, y = Y, phi = Phi, lambda_min_ratio = 0.001)
res$lambda_min
coef(res, s=res$lambda_min)
predict(res, type = "nonzero", s = res$lambda_min)
coef(res)

plot(res, xvar = "norm")
plot(res, xvar = "lambda")
plot(res, xvar = "dev")

plot(res, type = "BIC", sign.lambda = 1)



plot(coef(res, s = res$lambda_min), pch = 19, col = "red")
points(seq_along(c(b0,b, eta, sigma2)), c(b0,b, eta, sigma2), pch = 19, col = "blue")
legend("bottomleft",
       legend = c("Estimated", "Truth"),
       col = c("red","blue"),
       pch = c(19, 19),
       bg = "gray90")





# lambda sequence ---------------------------------------------------------

rm(list = ls())
library(microbenchmark)
source("~/git_repositories/penfam/R/functions.R")
source("~/git_repositories/penfam/R/sim-data.R")









# grplasso ----------------------------------------------------------------

pacman::p_load(grplasso)
View(grplasso:::grplasso.default)


# glmnet ------------------------------------------------------------------

fit <- glmnet::cv.glmnet(x = x, y = y)
plot(fit)
coef(fit)[(coef(fit)==0) , ]

beta_hats <- as.matrix(coef(fit, s = "lambda.1se"))
s0 <- beta_hats[beta_hats!=0,, drop=F]

# true positive rate
sum(paste0("V",1:20) %in% rownames(s0))/20

# true negative rate
sum(paste0("V",21:1000) %ni% rownames(s0))/980



# lmmlasso ----------------------------------------------------------------

data(classroomStudy)
head(classroomStudy)
fit1 <- lmmlasso(x=classroomStudy$X,y=classroomStudy$y,z=classroomStudy$Z,
                 grp=classroomStudy$grp,lambda=15,pdMat="pdIdent")
summary(fit1)
plot(fit1)


Xmat <- model.matrix(as.formula(paste("~", paste(colnames(x), collapse = "+"))), data = as.data.frame(x))
Xmat[1:5,1:5]
x[1:5,1:5]
fit1 <- lmmlasso(x = Xmat, y = y, z = as.matrix(rep(1, 600)), grp = 1:600, lambda = 29, pdMat = "pdIdent")
summary(fit1)
plot(fit1)

beta_hats <- as.matrix(fit1$coefficients)
s0 <- beta_hats[beta_hats!=0,, drop=F]

# true positive rate
sum(paste0("V",1:20) %in% rownames(s0))/20

# true negative rate
sum(paste0("V",21:1000) %ni% rownames(s0))/980




# quadform ----------------------------------------------------------------

pacman::p_load(emulator)
pacman::p_load(microbenchmark)

jj <- matrix(rnorm(100*1000),ncol = 1000)
M <- crossprod(jj,jj)
M.lower <- t(chol(M))
x <- matrix(rnorm(8),4,2)

jj.1 <- t(x) %*% M %*% x
jj.2 <- quad.form(M,x)
jj.3 <- quad.form(M.lower,x,chol=TRUE)
print(jj.1)
print(jj.2)
print(jj.3)




emulator::quad.form
cprod()
