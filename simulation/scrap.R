
# Trying low rank ---------------------------------------------------------

rm(list = ls())
pacman::p_load(gaston)
pacman::p_load(glmnet)
pacman::p_load(magrittr)
pacman::p_load(snpStats)
pacman::p_load_gh('StoreyLab/popkin')
pacman::p_load_gh('StoreyLab/bnpsd')
pacman::p_load(MASS)
devtools::load_all()
data("admixed")
data("karim")
karim$b %>% plot
Phi <- 2 * karim$kin1
P <- mvrnorm(1, rep(0,600), karim$s.g * Phi)
X <- karim$G
colnames(X) <- paste0("V",seq(ncol(X)))
ncausal <- 10
beta_mean <- 1
causal <- sample(colnames(X), ncausal, replace = FALSE)
not_causal <- setdiff(colnames(X), causal)
beta <- rep(0, length = ncol(X))
beta[which(colnames(X) %in% causal)] <- runif(n = length(causal), beta_mean - 0.2, beta_mean + 0.2)
plot(beta)
# beta[which(colnames(x) %in% causal)] <- rnorm(n = length(causal))
mu <- as.numeric(X %*% beta)
# mu <- as.numeric(X %*% karim$b)
# eta <- karim$s.g / (karim$s.e + karim$s.g)
# P <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = eta * sigma2 * kin)
# E <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))
# # y <- mu + sigma * matrix(rnorm(nsim * n), n, nsim)
# # y <- b0 + mu + t(P) + t(E)
# y <- MASS::mvrnorm(1, mu = mu, Sigma = eta * sigma2 * kin + (1 - eta) * sigma2 * diag(n))

y <- mu + P + rnorm(600, 0, karim$s.e)
hist(y)
phi_eigen <- eigen(Phi)
Phi[1:5,1:5]
# popkin::plotPopkin(karim$kin1)
phi_eigen$vectors[,2]
U_kinship <- phi_eigen$vectors
Lambda <- phi_eigen$values
any(Lambda < 1e-5)
plot(Lambda)
any(Lambda == 0)
PC <- sweep(phi_eigen$vectors, 2, sqrt(phi_eigen$values), "*")
plot(PC[,1],PC[,2])
X_lasso <- cbind(X,PC[,1:10])
dim(X_lasso)

source("simulation/model_functions.R")
dat <- make_INDmixed_model_not_simulator(n = 1000, p = 10000, ncausal = 100, k = 5, s = 0.5, Fst = 0.1,
                                         b0 = 1, beta_mean = 1,
                                         eta = 0.10, sigma2 = 4)
dat <- make_ADmixed_model_not_simulator(n = 1000,
                                        p_test = 2000,
                                        p_kinship = 5000,
                                        geography = "circ",
                                        percent_causal = 0.01,
                                        percent_overlap = "100",
                                        k = 5, s = 0.5, Fst = 0.1,
                                        b0 = 0, beta_mean = 1,
                                        eta = 0.1, sigma2 = 1)
phi_eigen <- eigen(dat$kin)
dat$kin[1:5,1:5]
popkin::plotPopkin(dat$kin)
U_kinship <- phi_eigen$vectors
Lambda <- phi_eigen$values
length(Lambda) # length n
dim(U_kinship) # n x n
any(Lambda < 1e-5)
which(Lambda < 1e-5)
Lambda[which(Lambda < 1e-5)] <- 1e-05
plot(Lambda)
any(Lambda == 0)
dev.off()
hist(dat$y)

correct_sparsity <- function(causal, not_causal, active, p){
  correct_nonzeros <- sum(active %in% causal)
  correct_zeros <- length(setdiff(not_causal, active))
  #correct sparsity
  (1 / p) * (correct_nonzeros + correct_zeros)
}

# ggmix -------------------------------------------------------------------

devtools::load_all()
# res <- lowrank(x = X, y = y,  d = Lambda, u = U_kinship)
# this is for karim data
res <- gic.penfam(x = X, y = y,  d = Lambda, u = U_kinship, an = log(length(y)))
res <- ggmix(x = admixed$x, y = admixed$y, kinship = admixed$kin,
             n_nonzero_eigenvalues = 10, estimation = "low")

# for make_INDmixed_model_not_simulator data and make_ADmixed_model_not_simulator data
res <- gic.penfam(x = dat$x, y = dat$y,  d = Lambda, u = U_kinship, an = log(length(dat$y)))
dev.off()
plot(res)
res$penfam.fit$result
(nonzero = res$penfam.fit$coef[,res$lambda.min.name,drop = F][nonzeroCoef(res$penfam.fit$coef[,res$lambda.min.name,drop = F]),,drop = F])
nonzero_names = setdiff(rownames(nonzero), c("beta0","eta","sigma2"))
length(intersect(nonzero_names, causal))/length(causal)
length(intersect(nonzero_names, dat$causal))/length(dat$causal)
length(nonzero_names)
res$penfam.fit$sigma2[,res$lambda.min.name]
res$penfam.fit$eta[,res$lambda.min.name]

correct_sparsity(causal = dat$causal, not_causal = dat$not_causal,
                 active = nonzero_names, p = ncol(dat$x))

l2norm(dat$x %*% res$penfam.fit$coef[colnames(dat$x),res$lambda.min.name,drop = F] -
         dat$x %*% matrix(dat$beta))
l2norm(res$penfam.fit$coef[colnames(dat$x),res$lambda.min.name,drop = F] -
         matrix(dat$beta))
plot(res$penfam.fit$coef[colnames(dat$x),res$lambda.min.name,drop = F])
sum(res$penfam.fit$coef[dat$causal,res$lambda.min.name,drop = F] > 0) / sum(dat$beta>0)
plot(res$penfam.fit$coef[dat$causal,res$lambda.min.name,drop = F], pch = 19, col = "blue")
points(dat$beta[which(dat$causal==colnames(dat$x))], pch = 19, col = "red")
###########$%$%#$%^#$%# Make sure that the first lambda sets everything to 0. its not
# curently doing this
# now fixed (june 14,2018)
# res$penfam.fit$coef[,1][which(res$penfam.fit$coef[,1] != 0)]

# lasso -------------------------------------------------------------------

# for karim data
fitglmnet <- cv.glmnet(x = X_lasso, y = y, penalty.factor = c(rep(1, 1000), rep(0, 10)))
plot(fitglmnet)
(nonzlasso <- setdiff(rownames(coef(fitglmnet, s = "lambda.min")[nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop=F]),c("(Intercept)","")))
(tprlasso <- length(intersect(nonzlasso, causal))/length(causal))
length(nonzlasso)

# for make_INDmixed_model_not_simulator data and make_ADmixed_model_not_simulator data
fitglmnet2 <- glmnet::cv.glmnet(x = dat$x_lasso, y = dat$y, standardize = T, alpha = 1, intercept = T,
                                penalty.factor = c(rep(1, 10000), rep(0, 10)))
plot(fitglmnet2)
# yhat2 = predict(fitglmnet2, newx = dat$x_lasso, s = "lambda.min")
# as.numeric(sqrt(crossprod(dat$y - yhat2)))
(nonzlasso <- setdiff(rownames(coef(fitglmnet2, s = "lambda.min")[nonzeroCoef(coef(fitglmnet2, s = "lambda.min")),,drop=F]),c("(Intercept)","")))
(tprlasso <- length(intersect(nonzlasso, dat$causal))/length(dat$causal))
length(nonzlasso)

correct_sparsity(causal = dat$causal, not_causal = dat$not_causal,
                 active = nonzlasso, p = ncol(dat$x))
l2norm(dat$x %*% coef(fitglmnet2, s = "lambda.min")[colnames(dat$x),,drop = F] -
         dat$x %*% matrix(dat$beta))

l2norm(coef(fitglmnet2, s = "lambda.min")[colnames(dat$x),,drop = F] -
         matrix(dat$beta))
plot(coef(fitglmnet2, s = "lambda.min")[colnames(dat$x),,drop = F])
plot(coef(fitglmnet2, s = "lambda.min")[dat$causal,,drop = F])
# two-step ----------------------------------------------------------------

#for karim data
pheno_dat <- data.frame(Y = y, id = paste0("ID",1:length(y)))
x1 <- cbind(rep(1, nrow(X)))
fit <- gaston::lmm.aireml(y, x1, K = Phi)
gaston_resid <- y - (fit$BLUP_omega + fit$BLUP_beta)
hist(gaston_resid)
twostep <- glmnet::cv.glmnet(x = X, y = gaston_resid, standardize = T, alpha = 1, intercept = T)
plot(twostep)
(nonz2step <- setdiff(rownames(coef(twostep, s = "lambda.min")[nonzeroCoef(coef(twostep, s = "lambda.min")),,drop=F]),c("(Intercept)")))
(tpr2step <- length(intersect(nonz2step, causal))/length(causal))
length(nonz2step)

# for make_INDmixed_model_not_simulator data and make_ADmixed_model_not_simulator data
pheno_dat <- data.frame(Y = dat$y, id = paste0("ID",1:length(dat$y)))
x1 <- cbind(rep(1, nrow(dat$x)))
fit <- gaston::lmm.aireml(dat$y, x1, K = dat$kin)
gaston_resid <- dat$y - (fit$BLUP_omega + fit$BLUP_beta)
hist(gaston_resid)
fitglmnet <- glmnet::cv.glmnet(x = dat$x, y = gaston_resid, standardize = T, alpha = 1, intercept = T)
plot(fitglmnet)
(nonz2step <- setdiff(rownames(coef(fitglmnet, s = "lambda.min")[nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop = F]),c("(Intercept)")))
(tpr2step <- length(intersect(nonz2step, dat$causal))/length(dat$causal))
length(nonz2step)

correct_sparsity(causal = dat$causal, not_causal = dat$not_causal,
                 active = nonz2step, p = ncol(dat$x))
l2norm(dat$x %*% coef(fitglmnet, s = "lambda.min")[colnames(dat$x),,drop = F] -
         dat$x %*% matrix(dat$beta))
l2norm(coef(fitglmnet, s = "lambda.min")[colnames(dat$x),,drop = F] -
         matrix(dat$beta))
plot(coef(fitglmnet, s = "lambda.min")[colnames(dat$x),,drop = F])
plot(coef(fitglmnet, s = "lambda.min")[dat$causal,,drop = F])



# dat <- make_mixed_model_not_simulator(b0 = 1, eta = 0.3, sigma2 = 2, type = "causal_400", related = TRUE)
# w_svd <- svd(karim$G)
# U <- w_svd$u
# dim(U)
# length(w_svd)
# This part shows the equivalence between the eigenvectors and eigenvalues
# from the kinship matrix vs. the SNP matrix used to construct the kinship
# U_w=U
# plot(U_w[,1],U_w[,2])
# phi_eigen <- eigen(phi)
# U_kinship <- phi_eigen$vectors
# plot(U_kinship[,1], U_w[,1])
# abline(a=0, b=-1)
# dim(U)
# vector of length N_T
# Lambda <- phi_eigen$values
# plot(Lambda, w_svd$d^2 / (p-1))
# abline(a=0,b=1,col="red")
# all.equal(Lambda, w_svd$d^2 / (p-1))
# we dived by p-1 because thats how the matrix was standardized
# Lambda <- w_svd$d^2 / (ncol(dat$x)-1)
# Lambda <- w_svd$d^2
# if (any(Lambda<1e-5)) Lambda[Lambda<1e-5] <- 1e-5


bic <- function(eta, sigma2, beta, eigenvalues, x, y, nt, c, df_lambda) {

  -2 * log_lik(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt) + c * df_lambda

}




bic(eta = res$eta, sigma2 = res$sigma2)

lasso <- cv.glmnet(x = dat$x, y = dat$y, alpha = 1)
fitted(lasso)
plot(lasso)

coef(lasso, s="lambda.min")[nonzeroCoef(coef(lasso, s="lambda.min")),,drop=F]
grep("rs", rownames(coef(lasso)[nonzeroCoef(coef(lasso)),,drop=F]), value = T) %in% dat$causal
grep("rs", rownames(coef(res)[nonzeroCoef(coef(res)),,drop=F]), value = T) %in% dat$causal

coef(lasso)[nonzeroCoef(coef(lasso)),,drop=F]
coef(res)[nonzeroCoef(coef(res)),,drop=F]

dat$beta[dat$beta!=0]

plot(coef(lass0)[-1], dat$beta)
abline(a=0, b=1)
all(colnames(dat$x)==rownames(coef(lass0)[-1]))

coef(res, s = res$lambda_min) %>% head
plot(coef(resenet, s = resenet$lambda_min)[c(-1,-4002,-4003),,drop=T], dat$beta)
abline(a=0, b=1)
all(colnames(dat$x)==rownames(coef(lass0)[-1]))













# this needs work, as lambda.min is different
resenet <- lowrank(x = dat$x, y = dat$y, w = dat$w, alpha = 1)

dim(dat$x)

par(mfrow=c(3,2))
plot(res)
# c("coef","BIC", "QQranef","QQresid", "predicted", "Tukey-Anscombe")
plot(res, type = "BIC")
plot(res, type = "QQranef")
plot(res, type = "QQresid")
plot(res, type = "predicted")
plot(res, type = "Tukey")

plot(resenet, type = "BIC")
predict(res, type = "nonzero", s = resenet$lambda_min)

grep("rs", rownames(predict(res, type = "nonzero", s = res$lambda_min)), value = T) %in% dat$causal
grep("rs", rownames(predict(resenet, type = "nonzero", s = resenet$lambda_min)), value = T) %in% dat$causal

plot(resenet, type = "QQranef")
plot(resenet, type = "QQresid")
plot(resenet, type = "predicted")
plot(resenet, type = "Tukey")


plot(res, type = "")
dev.off()
predict(res, type = "nonzero", s = res$lambda_min)
dat$causal
res$eta
res$sigma2

lass0 <- cv.glmnet(x = dat$x, y = dat$y, alpha = 0.5)
plot(lass0)

coef(lass0)[nonzeroCoef(coef(lass0)),,drop=F]
grep("rs", rownames(coef(lass0)[nonzeroCoef(coef(lass0)),,drop=F]), value = T) %in% dat$causal

dat$beta[dat$beta!=0]

plot(coef(lass0)[-1], dat$beta)
abline(a=0, b=1)
all(colnames(dat$x)==rownames(coef(lass0)[-1]))

coef(res, s = res$lambda_min) %>% head
plot(coef(resenet, s = resenet$lambda_min)[c(-1,-4002,-4003),,drop=T], dat$beta)
abline(a=0, b=1)
all(colnames(dat$x)==rownames(coef(lass0)[-1]))


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


# checking lowrank --------------------------------------------------------

pacman::p_load(gaston)
pacman::p_load(glmnet)
pacman::p_load(magrittr)
pacman::p_load(snpStats)
pacman::p_load_gh('StoreyLab/popkin')
pacman::p_load_gh('StoreyLab/bnpsd')
pacman::p_load(MASS)
pacman::p_load(RSpectra)
devtools::load_all()
k=5;n=10000;p_kinship=3000;s = 0.5; Fst = 0.1;
FF <- 1:k # subpopulation FST vector, up to a scalar
obj <- bnpsd::q1d(n = n, k = k, s = s, F = FF, Fst = Fst)
Q <- obj$Q
FF <- obj$F
out <- bnpsd::rbnpsd(Q, FF, p_kinship)
Xall <- t(out$X) # genotypes are columns, rows are subjects
dim(Xall)

# On the WTCCC data, use of fewer than 200 eigenvectors yielded univariate
# P values comparable to those obtained from many thousands of eigenvectors. (Lippert et al. 2009)

ir <- RSpectra::svds(Xall, k = 200)
ir$d %>% hist
summary(ir$d)
dim(ir$u)
PC <- sweep(ir$u, 2, STATS = ir$d, FUN = "*")
dim(PC)
plot(PC[,1],PC[,2], pch=19, col = rep(RColorBrewer::brewer.pal(5,"Paired"), each = n/k))
svdX <- svd(Xall,nu = 5)
svdX$u %>% dim
svdX$d %>% length()
svdX$v %>% dim
U <- svdX$u
dim(U)
U1 <- U[,1:40]
U2 <- U[,41:100]

set.seed(1)
A <- matrix(rnorm(5000*5000), 5000)
t1 <- proc.time()
L <- irlba(A, 5)
print(proc.time() - t1)

# checking lower rank approx ----------------------------------------------

library(magrittr)
options(digits = 4, scipen = 999)
X <- matrix(rnorm(n=100*40), ncol = 40)
X <- matrix(rnorm(n=100*40), ncol = 100)
dim(X)
Y <- rnorm(100)
Y <- rnorm(40)
xtx <- tcrossprod(scale(X))
xtx <- tcrossprod(X)
dim(xtx)
eigX <- eigen(xtx)

U <- eigX$vectors
eigX$values %>% length()
dim(U)
crossprod(U)[1:5,1:5]
tcrossprod(U)[1:5,1:5]

svdX <- svd(X)
dim(X)
svdX$u %>% dim
svdX$d %>% length()

plot(svdX$d^2,eigX$values)
all.equal(svdX$d^2,eigX$values)
abline(a=0,b=1)
plot(U[,1],svdX$u[,1])

svdX$v %>% dim
U <- svdX$u
dim(U)
U1 <- U[,1:40]
U2 <- U[,41:100]

(t(U2) %*% U2)[1:5,1:5]
(t(U2) %*% U2) %>% dim
(U2 %*% t(U2))[1:5,1:5]
crossprod(U2)[1:5,1:5]

round(tcrossprod(U2)[1:5,1:5] + tcrossprod(U1)[1:5,1:5], 2)

svdX$v %>% dim


crossprod(U1,Y) %>% dim


U1 %*% crossprod(U1,Y) %>% dim

dim(U1)


Matrix::rankMatrix(X)



eta = 0.6
I_nk = diag(5)
I_nk - eta * I_nk

all.equal((1 - eta) * I_nk, I_nk - eta * I_nk)


all.equal(solve((1 - eta) * I_nk), 1/(1 - eta) * I_nk)

