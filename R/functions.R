"%ni%" <- Negate("%in%")

# function used to optimize eta
fr_eta <- function(eta, sigma2, beta, eigenvalues, x, y, nt) {

  # this is based on the negative log-lik

  # eta = 0.5
  # sigma2 = 1
  # beta = beta_next
  # eigenvalues = Lambda
  # x = utx
  # y = uty
  # nt = n
  # ============

  di <- 1 + eta * (eigenvalues - 1)

  (nt / 2) * log(2 * pi) +
    (nt / 2) * log(sigma2) +
    0.5 * sum(log(di)) +
    (1 / (2 * sigma2)) * sum((y - x %*% beta) ^ 2 / di)

}


#' @param x should be U^T X, where U is the matrix of eigenvectors and X
#'   contains the first column of ones for the intercept. Used for gradient of eta
#'   (currently being passed to optim)
grr_eta <- function(eta, sigma2, beta, eigenvalues, x, y, nt) {

  di <- 1 + eta * (eigenvalues - 1)

  (1 / 2) * sum(((eigenvalues - 1) / di) * (1 - (((y - x %*% beta) ^ 2) / (sigma2 * di))))
}

# gradient of sigma2 (used for KKT check)
grr_sigma2 <- function(eta, sigma2, beta, eigenvalues, x, y, nt) {

  di <- 1 + eta * (eigenvalues - 1)

  sigma2 - (1 / nt) * sum((((y - x %*% beta) ^ 2) / di))

}


grr_beta0 <- function(eta, sigma2, beta, eigenvalues, x, y, nt) {

  di <- 1 + eta * (eigenvalues - 1)
  wi <- (1 / sigma2) * diag(1 / di)

  as.numeric(crossprod(x[,1, drop = FALSE], wi) %*% (y - x %*% beta))

}

#' Check of KKT for a given value of Lambda
#' @param x should be U^T X, where U is the matrix of eigenvectors and X
#'   contains the first column of ones for the intercept. x should be a mtrix of
#'   dimension n x (p+1)
#' @param beta should include intercept as well. A 1 column matrix of dimension
#'   (p+1) x 1.
#' @param tol.kkt Tolerance for determining if an entry of the subgradient is
#'   zero

kkt_check <- function(eta, sigma2, beta, eigenvalues, x, y, nt,
                      lambda, tol.kkt = 1e-9){

  # eta = eta_next; sigma2 = sigma2_next;
  # beta = beta_next;
  # # beta = as.matrix(c(beta0_next,rep(0,p)))
  # eigenvalues = Lambda;
  # x = utx;
  # y = uty; nt = n;
  # lambda = lambda; tol.kkt = 0.000001
  # =======================

  di <- 1 + eta * (eigenvalues - 1)
  wi <- (1 / sigma2) * (1 / di)


  # scale the weights to sum to nvars
  wi_scaled <- as.vector(wi) / sum(as.vector(wi)) * nt
  # wi_mat <- diag(wi_scaled)

  # KKT for beta0
  # grr_beta0(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt)
  # kkt_beta0 <- abs(grr_beta0(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt)) < tol.kkt

  kkt_beta0 <- grr_beta0(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt)

  # KKT for eta
  # grr_eta(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt)
  # kkt_eta <- grr_eta(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt) < tol.kkt
  kkt_eta <- grr_eta(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt)

  # KKT for sigma2
  # grr_sigma2(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt)
  # kkt_sigma2 <- grr_sigma2(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt) < tol.kkt
  kkt_sigma2 <- grr_sigma2(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt)

  # KKT for beta
  # g0 <- (1 / nt) * crossprod(x[,-1, drop = F], wi) %*% (y - x %*% beta) / (colSums(sweep(x[,-1, drop = F]^2, MARGIN = 1, wi_vec, '*')))

  # KKT for beta
  g0 <- (1 / sum(wi_scaled)) * crossprod(x[,-1, drop = F] * wi_scaled, (y - x %*% beta ))

  # this gives same result as g0
  # g1 <- colSums((1 / nt) * sweep(sweep(x[,-1], MARGIN = 1, wi_vec, '*'), MARGIN = 1, drop((y - x %*% beta)),'*'))
  # plot(drop(g0),g1)
  # abline(a=0,b=1)

  g <- g0 - lambda * sign(beta[-1])

  # this is for when beta=0 and should be between -1 and 1
  gg <- g0 / lambda

  # which of the betas are non-zero
  oo <- abs(beta[-1]) > 0

  # if all betas are 0 then set to TRUE, else abs(g[oo]) will give error since 'oo' is all FALSE
  # kkt_beta_nonzero <- if (all(!oo)) TRUE else max(abs(g[oo])) < tol.kkt
  # kkt_beta_subgr1 <- min(gg[!oo]) > -1
  # kkt_beta_subgr2 <- max(gg[!oo]) < 1

  kkt_beta_nonzero <- if (all(!oo)) 0 else sum(abs(g[oo]) > tol.kkt)
  kkt_beta_subgr <- sum(abs(gg[!oo]) > 1)
  if (sum(abs(g[oo]) > tol.kkt) > 0) plot(abs(g[oo]))

  return(c(kkt_beta0 = kkt_beta0,
           kkt_eta = kkt_eta,
           kkt_sigma2 = kkt_sigma2,
           kkt_beta_nonzero = kkt_beta_nonzero,
           kkt_beta_subgr = kkt_beta_subgr))


}

# log likelihood (used to calculate BIC)
log_lik <- function(eta, sigma2, beta, eigenvalues, x, y, nt) {

  # this returns the log-likelihood

  # eta = 0.5
  # sigma2 = 1
  # beta = beta_next
  # eigenvalues = Lambda
  # x = utx
  # y = uty
  # nt = n
  # ============

  kernel <- 1 + eta * (eigenvalues - 1)

  -1 * (
     (nt / 2) * log(2 * pi) +
       (nt / 2) * log(sigma2) +
       0.5 * sum(log(kernel)) +
      (1 / (2 * sigma2)) * sum((y - x %*% beta) ^ 2 / kernel)
  )

}


#' Calculate Sequence of Tuning Parameters
#'
#' @description Function to calculate the sequence of tuning parameters based on
#'   the design matrix \code{x} and the response variable {y}. This is used in
#'   the \code{\link{shim_once}} function to calculate the tuning parameters
#'   applied to the main effects
#'
#' @inheritParams uni_fun
#' @param weights Separate penalty factors can be applied to each coefficient.
#'   This is a number that multiplies lambda to allow differential shrinkage,
#'   and can be used to apply adaptive LASSO. Can be 0 for some variables, which
#'   implies no shrinkage, and that variable is always included in the model.
#'   Default is 1 for all variables (and implicitly infinity for variables
#'   listed in exclude). Note: the penalty factors are internally rescaled to
#'   sum to nvars, and the lambda sequence will reflect this change.
#' @param lambda.factor The factor for getting the minimal lambda in lambda
#'   sequence, where \code{min(lambda) = lambda.factor * max(lambda).
#'   max(lambda)} is the smallest value of lambda for which all coefficients are
#'   zero. The default depends on the relationship between \code{N} (the number
#'   of rows in the matrix of predictors) and \code{p} (the number of
#'   predictors). If \code{N > p}, the default is \code{1e-6}, close to zero. If
#'   \code{N < p}, the default is \code{0.01}. A very small value of
#'   lambda.factor will lead to a saturated fit.
#' @param nlambda the number of lambda values - default is 100.
#' @param scale_x should the columns of x be scaled - default is FALSE
#' @param center_y should y be mean centered - default is FALSE
#' @return numeric vector of length \code{q}
#' @details The maximum lambda is calculated using the following inequality:
#'   \deqn{(N*w_j)^-1 | \sum x_ij y_i | \le \lambda_max}
#'
#'   The minimum lambda is given by lambda.factor*lambda_max. The sequence of
#'   nlambda values are decreasing from lambda_max to lambda_min on the log
#'   scale.
#'
#'   The penalty factors are internally rescaled to sum to the number of
#'   predictor variables in glmnet. Therefore, to get the correct sequence of
#'   lambdas when there are weights, this function first rescales the weights
#'   and then calclated the sequence of lambdas.
#'
#'   This formula is taken from section 2.5 of the \code{glmnet} paper in the
#'   Journal of Statistical Software (see references for details)
#'
#' @author
#' Sahir Bhatnagar
#'
#' Maintainer: Sahir Bhatnagar \email{sahir.bhatnagar@@mail.mcgill.ca}
#'
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008)
#'   \emph{Regularization Paths for Generalized Linear Models via Coordinate
#'   Descent}, \url{http://www.stanford.edu/~hastie/Papers/glmnet.pdf}
#'   \emph{Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010}
#'   \url{http://www.jstatsoft.org/v33/i01/}
#'
#'   Yang, Y., & Zou, H. (2015). A fast unified algorithm for solving
#'   group-lasso penalize learning problems. \emph{Statistics and Computing},
#'   25(6), 1129-1141.
#'   \url{http://www.math.mcgill.ca/yyang/resources/papers/gglasso.pdf}
#'
#'
#' @examples
#' # number of observations
#' n <- 100
#'
#' # number of predictors
#' p <- 5
#'
#' # environment variable
#' e <- sample(c(0,1), n, replace = T)
#'
#' # main effects
#' x <- cbind(matrix(rnorm(n*p), ncol = p), e)
#'
#' # need to label columns
#' dimnames(x)[[2]] <- c(paste0("x",1:p), "e")
#'
#' # design matrix without intercept
#' X <- model.matrix(~(x1+x2+x3+x4+x5)*e-1, data = as.data.frame(x))
#'
#' # response
#' Y <- X %*% rbinom(ncol(X), 1, 0.2) + 3*rnorm(n)
#'
#' lambda_sequence(X,Y)
#' @export

lambda_sequence <- function(x, y, eigenvalues, weights = NULL,
                            # lambda_min_ratio = ifelse(n < p, 0.01, 0.001),
                            lambda_min_ratio,
                            thresh_penfam = 1e-14,
                            tol.kkt = 1e-9,
                            eta_init = 0.5,
                            nlambda = 100, scale_x = F, center_y = F) {

  # rm(list=ls())
  # source("~/git_repositories/penfam/R/fitting.R")
  # source("~/git_repositories/penfam/R/functions.R")
  # source("~/git_repositories/penfam/R/methods.R")
  # source("~/git_repositories/penfam/R/plot.R")
  # source("~/git_repositories/penfam/R/sim-data.R")
  # x <- X
  # y <- drop(Y)
  # phi <- Phi
  # eta_init <- .5
  # nlambda = 100
  # lambda_min_ratio = 0.001
  # # column of 1s for intercept
  # x <- cbind(1, x)
  # dim(x); x[1:5,1:5]
  # phi_eigen <- eigen(phi)
  # U <- phi_eigen$vectors
  # dim(U)
  # Lambda <- phi_eigen$values
  #
  # # utx0 <- crossprod(U, x0)
  # # to be consistent with the paper and all other functions
  # # lets keep the intercept column in utx, and then define
  # # utx0 as the first column of utx, instead of keeping them separate
  # x <- crossprod(U, x)
  # y <- crossprod(U, y)
  # eigenvalues <- Lambda
  # # utx0 <- utx[, 1, drop = F]
  #convergence criterion
  # epsilon <- 1e-7

  #======================================

  np <- dim(x)
  n <- np[[1]]

  # assuming the first column is the intercept, so we subtract 1
  p <- np[[2]] - 1

  # weights
  di <- 1 + eta_init * (eigenvalues - 1)
  # di_inverse <- diag(1 / di)

  # initial value for beta0
  beta0_init <- (sum(x[, 1] * y / di)) / (sum(x[,1] ^ 2 / di))

  # this includes all other betas which are 0 by definition
  beta_init <- as.matrix(c(beta0_init, rep(0,p)))
  # sum is faster
  # microbenchmark::microbenchmark(
  #   mat = drop((t(utx0) %*% di_inverse %*% uty) / (t(utx0) %*% di_inverse %*% utx0)),
  #     sum = (sum(utx0 * uty / di)) / (sum(utx0 ^ 2 / di)), times = 1000
  # )

  # closed form for sigma^2
  sigma2_init <- (1 / n) * sum ((y - x %*% beta_init) ^ 2 / di)

  # sum version is faster
  # mb <- microbenchmark(
  #   mat = (1 / n) * t(uty - beta0_init * utx0) %*% di_inverse %*% (uty - beta0_init * utx0),
  #   sum = (1 / n) * sum ((uty - beta0_init * utx0)^2 / (1 + eta_init * (eigenvalues - 1))),
  #   times = 1000)
  # ggplot2::autoplot(mb)

  #iteration counter
  k <- 0

  # to enter while loop
  converged <- FALSE

  while (!converged) {

    Theta_init <- c(beta_init, eta_init, sigma2_init)

    # fit eta
    eta_next <- optim(par = eta_init,
                      fn = fr_eta,
                      gr = grr_eta,
                      method = "L-BFGS-B",
                      control = list(fnscale = 1),
                      lower = .01,
                      upper = .99,
                      sigma2 = sigma2_init,
                      beta = beta_init,
                      eigenvalues = eigenvalues,
                      x = x,
                      y = y,
                      nt = n)$par

    # weights
    di <- 1 + eta_next * (eigenvalues - 1)

    # di_inverse <- diag(1 / di)

    # next value for beta0
    beta0_next <- (sum(x[,1] * y / di)) / (sum(x[,1] ^ 2 / di))

    beta_next <- as.matrix(c(beta0_next, rep(0, p)))

    # closed form for sigma^2
    sigma2_next <- (1 / n) * sum ((y - x %*% beta_next) ^ 2 / di)

    k <- k + 1

    Theta_next <- c(beta_next, eta_next, sigma2_next)

    converged <- crossprod(Theta_next - Theta_init) < thresh_penfam
    # converged <- max(abs(Theta_next - Theta_init) / (1 + abs(Theta_next))) < epsilon

    message(sprintf("l2 norm squared of Theta_k+1 - Theta_k: %f \n log-lik: %f",
                    crossprod(Theta_next - Theta_init),
                    log_lik(eta = eta_next, sigma2 = sigma2_next, beta = beta_next,
                            eigenvalues = eigenvalues,
                            x = x, y = y, nt = n)))

    beta0_init <- beta0_next
    beta_init <- beta_next
    eta_init <- eta_next
    sigma2_init <- sigma2_next

  }

  # eta_next
  # sigma2_next
  di <- 1 + eta_next * (eigenvalues - 1)
  wi <- (1 / sigma2_next) * (1 / di)
  if (any(wi < 0)) stop("weights are negative")

  # scale the weights to sum to nvars
  wi_scaled <- as.vector(wi) / sum(as.vector(wi)) * n

  # wi_scaled <- as.vector(wi) * n

  # lambda.max <- max(abs(colSums((wi * x[,-1]) * drop(y - x %*% beta_next))))

  lambda.max <- max(abs(colSums(((1 / sum(wi_scaled)) * (wi_scaled * x[,-1]) * drop(y - x %*% beta_next)))))

  # (x[,-1, drop = F]) %>% dim
  # a <- colSums(x[,-1, drop = F]^2 * wi)
  # b <- colSums(sweep(x[,-1, drop = F]^2, MARGIN = 1, wi, '*'))
  # all(a == b)


  kkt <- kkt_check(eta = eta_next, sigma2 = sigma2_next, beta = beta_next,
                   eigenvalues = eigenvalues, x = x, y = y, nt = n,
                   lambda = lambda.max, tol.kkt = tol.kkt)
  # message(kkt)
  out <- list(sequence = rev(exp(seq(log(lambda_min_ratio * lambda.max), log(lambda.max), length.out = nlambda))),
              eta = eta_next, sigma2 = sigma2_next, beta0 = beta0_next, kkt = kkt)

}

bic <- function(eta, sigma2, beta, eigenvalues, x, y, nt, c, df_lambda) {

  -2 * log_lik(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt) + c * df_lambda

}






