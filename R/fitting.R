# options(warnPartialMatchArgs = FALSE, warnPartialMatchDollar = TRUE, warnPartialMatchAttr = TRUE)
library(magrittr)
library(bit64)
library(data.table)
library(glmnet)
library(MASS)
# library(lmmlasso)

# Simulate Data -----------------------------------------------------------


rm(list = ls())
source("~/git_repositories/penfam/R/functions.R")
source("~/git_repositories/penfam/R/sim-data.R")




# penfam <- function(x, y, phi, lambda) {

  x <- X
  y <- Y
  phi <- Phi
  lambda <- 0.10
  #======================================

  np <- dim(x)
  n <- np[[1]]
  p <- np[[2]]

  # add column of 1s to x for intercept
  x <- cbind(rep(1, n), x)
  x[1:5, 1:5]

  phi_eigen <- eigen(phi)
  U <- phi_eigen$vectors
  dim(U)
  Lambda <- phi_eigen$values

  utx <- crossprod(U, x)
  uty <- crossprod(U, y)

  # initial values for beta
  beta_init_fit <- glmnet(x = utx,
                          y = uty,
                          family = "gaussian",
                          penalty.factor = c(0, rep(1, p)),
                          standardize = FALSE,
                          intercept = FALSE,
                          lambda = lambda)

  # plot(beta_init_fit)
  coef(beta_init_fit)[nonzeroCoef(coef(beta_init_fit)), , drop = F]

  # remove intercept since V1 is intercept
  beta_init <- coef(beta_init_fit)[-1, , drop = F]
  head(beta_init)

  # initial value for eta
  eta_init <- .01

  # closed form solution value for sigma^2
  sigma2_init <- (1 / n) * sum(((uty - utx %*% beta_init) ^ 2) / (1 + eta_init * (Lambda - 1)))

  #iteration counter
  k <- 0

  #convergence criterion
  epsilon <- 1e-7

  # to enter while loop
  converged <- FALSE

  while (!converged) {

    Theta_init <- c(drop(beta_init), eta_init, sigma2_init)

    # observation weights
    wi <- (1 / sigma2_init) * (1 / (1 + eta_init * (Lambda - 1)))
    length(wi)
    # plot(wi)
    # are all weights positive?
    all(wi > 0)

    # fit beta
    beta_next_fit <- glmnet(x = utx,
                            y = uty,
                            family = "gaussian",
                            weights = wi,
                            penalty.factor = c(0, rep(1, p)),
                            standardize = FALSE,
                            intercept = FALSE,
                            lambda = lambda)
    coef(beta_next_fit)[nonzeroCoef(coef(beta_next_fit)),, drop = F]

    beta_next <- coef(beta_next_fit)[-1,, drop = F]
    # plot(beta_next)

    # fit eta
    eta_next <- optim(par = eta_init,
                      fn = fr_eta,
                      # gr = grr_eta,
                      method = "L-BFGS-B",
                      control = list(fnscale = 1),
                      lower = 1e-4,
                      upper = 1 - 1e-4,
                      sigma2 = sigma2_init,
                      beta = beta_next,
                      eigenvalues = Lambda,
                      x = utx,
                      y = uty,
                      nt = n)$par

    # fit sigma (closed form)
    sigma2_next <- (1 / n) * sum(((uty - utx %*% beta_next) ^ 2) / (1 + eta_next * (Lambda - 1)))

    k <- k + 1

    Theta_next <- c(drop(beta_next), eta_next, sigma2_next)

    converged <- crossprod(Theta_next - Theta_init) < epsilon

    beta_init <- beta_next
    eta_init <- eta_next
    sigma2_init <- sigma2_next

    message(sprintf("l2 norm of Theta: %f \n log-lik: %f", crossprod(Theta_next - Theta_init),
                    log_lik(eta = eta_next, sigma2 = sigma2_next, beta = beta_next, eigenvalues = Lambda,x = utx, y = uty, nt = n)))

  }


  sigma2_next
  eta_next
  plot(Theta_next[1:1001])

  }

beta_init %>% head
utx %*% beta_init






# lambda sequence ---------------------------------------------------------

rm(list = ls())
library(microbenchmark)
source("~/git_repositories/penfam/R/functions.R")
source("~/git_repositories/penfam/R/sim-data.R")
# lambda_sequence <- function(x, y, phi, weights = NULL,
#                             lambda.factor = ifelse(nobs < nvars, 0.01, 1e-06),
#                             nlambda = 100, scale_x = F, center_y = F) {

x <- X
y <- drop(Y)
phi <- Phi
#======================================

n <- length(y)

# column of 1s for intercept
x0 <- cbind(rep(1, n))

phi_eigen <- eigen(phi)
U <- phi_eigen$vectors
dim(U)
Lambda <- phi_eigen$values

utx0 <- crossprod(U, x0)
utx <- crossprod(U, x)
uty <- crossprod(U, y)

# initial value for eta
eta_init <- .1

# weights
di <- 1 + eta_init * (Lambda - 1)
di_inverse <- diag(1 / di)

# initial value for beta0
beta0_init <- drop((t(utx0) %*% di_inverse %*% uty) / (t(utx0) %*% di_inverse %*% utx0))

# closed form for sigma^2
sigma2_init <- (1 / n) * sum ((uty - beta0_init * utx0) ^ 2 / di)

# sum version is faster
# mb <- microbenchmark(
#   mat = (1 / n) * t(uty - beta0_init * utx0) %*% di_inverse %*% (uty - beta0_init * utx0),
#   sum = (1 / n) * sum ((uty - beta0_init * utx0)^2 / (1 + eta_init * (Lambda - 1))),
#   times = 1000)
# ggplot2::autoplot(mb)

#iteration counter
k <- 0

#convergence criterion
epsilon <- 1e-7

# to enter while loop
converged <- FALSE

while (!converged) {

  Theta_init <- c(beta0_init, eta_init, sigma2_init)

  # fit eta
  eta_next <- optim(par = eta_init,
                    fn = fr_eta,
                    # gr = grr_eta,
                    method = "L-BFGS-B",
                    control = list(fnscale = 1),
                    lower = 1e-5,
                    upper = 1 - 1e-5,
                    sigma2 = sigma2_init,
                    beta = beta0_init,
                    eigenvalues = Lambda,
                    x = utx0,
                    y = uty,
                    nt = n)$par

  # weights
  di <- 1 + eta_next * (Lambda - 1)
  di_inverse <- diag(1 / di)

  # next value for beta0
  beta0_next <- drop((t(utx0) %*% di_inverse %*% uty) / (t(utx0) %*% di_inverse %*% utx0))

  # closed form for sigma^2
  sigma2_next <- (1 / n) * sum ((uty - beta0_next * utx0) ^ 2 / di)

  k <- k + 1

  Theta_next <- c(beta0_next, eta_next, sigma2_next)

  converged <- crossprod(Theta_next - Theta_init) < epsilon

  beta0_init <- beta0_next
  eta_init <- eta_next
  sigma2_init <- sigma2_next

  message(sprintf("l2 norm of Theta: %f \n log-lik: %f", crossprod(Theta_next - Theta_init),
                  log_lik(eta = eta_next, sigma2 = sigma2_next, beta = beta0_next, eigenvalues = Lambda,x = utx0, y = uty, nt = n)))

}

eta_next
sigma2_next
wi <- (1 / sigma2_next) * (1 + eta_next * (Lambda - 1))

lambda.max <- max(abs(colSums(as.vector(uty) * utx / wi))) / nrow(utx)


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
