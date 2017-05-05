# options(warnPartialMatchArgs = FALSE, warnPartialMatchDollar = TRUE, warnPartialMatchAttr = TRUE)
library(magrittr)
library(bit64)
library(data.table)
library(glmnet)
library(MASS)
# library(lmmlasso)

# Simulate Data -----------------------------------------------------------


rm(list=ls())
source("~/git_repositories/penfam/R/functions.R")

# genotypes
N.MAF <- 1000
MAF <- runif(N.MAF,0.05,0.25)
x <- as.matrix(sapply(MAF,fun<-function(x) rbinom(600,2,x)))
dim(x)
x[1:5,1:5]
# link to a 600 x 600 theoretical kinship matrix https://www.dropbox.com/s/bk8zwqtzgvvm89z/kin1.Rdata?dl=0
# load kinship matrix from dropbox link
load("~/Dropbox/PhD/Year4/penfam/data/kin1.Rdata")
Phi <- 2 * kin1
Phi[1:10,1:10]
# random effect
s.g = 2 # variance of the polygenic random effect
P = mvrnorm(1, rep(0,600), s.g*Phi)
s.e = 1 # residual-effect variance
b <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, rep(0,980))
#
X <- mvrnorm(600, rep(1,1000), diag(1,1000,1000))
dim(X)
X[1:5,1:5]
Y <- 3 + X%*%b + P + rnorm(600,0,s.e)




penfam <- function(x, y, phi, lambda) {

  x = X
  y = Y
  phi = Phi
  lambda = 0.5
  #======================================

  np <- dim(x)
  n <- np[[1]]
  p <- np[[2]]

  # add column of 1s to x for intercept
  x <- cbind(rep(1,n),x)
  x[1:5,1:5]

  phi_eigen <- eigen(phi)
  U <- phi_eigen$vectors
  dim(U)
  Lambda <- phi_eigen$values

  utx <- crossprod(U,x)
  uty <- crossprod(U,y)

  # initial values for beta
  beta_init_fit <- glmnet(x = utx,
                          y = uty,
                          family = "gaussian",
                          penalty.factor = c(0, rep(1, p)),
                          standardize = FALSE,
                          intercept = FALSE,
                          lambda = lambda)

  plot(beta_init_fit)
  coef(beta_init_fit)[nonzeroCoef(coef(beta_init_fit)),,drop = F]

  # remove intercept since first row is intercept
  beta_init <- coef(beta_init_fit)[-1,,drop=F]
  head(beta_init)

  # initial value for eta
  eta_init = 0.5

  # initial value for sigma^2
  sigma2_init = 0.5

  # initial observation weights
  wi_init <- (1/sigma2_init) * (1/ (1 + eta_init * (Lambda-1)))
  length(wi_init)
  plot(wi_init)
  # are all weights positive?
  all(wi_init>0)

  #iteration counter
  k = 0

  #convergence criterion
  epsilon = 1e-7

  # to enter while loop
  converged = FALSE

  while(!converged) {

    # fit beta
    beta_next_fit <- glmnet(x = utx,
                            y = uty,
                            family = "gaussian",
                            weights = wi_init,
                            penalty.factor = c(0, rep(1, p)),
                            standardize = FALSE,
                            intercept = FALSE,
                            lambda = lambda)
    coef(beta_next_fit)[nonzeroCoef(coef(beta_next_fit)),,drop = F]

    beta_next <- coef(beta_next_fit)[-1, , drop = F]
    plot(beta_next)



    # fit eta

    optim(
      method = "L-BFGS-B"
    )




    k <- k + 1

  }


  }



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
