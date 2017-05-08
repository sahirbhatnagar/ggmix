"%ni%" <- Negate("%in%")

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

  kernel <- 1 + eta * (eigenvalues - 1)

  (nt / 2) * log(2 * pi) +
    (nt / 2) * log(sigma2) +
    0.5 * sum(log(kernel)) +
    (1 / (2 * sigma2)) * sum((y - x %*% beta) ^ 2 / kernel)

}

grr_eta <- function(eta, sigma2, beta, eigenvalues, x, y, nt) {

  kernel <- 1 + eta * (eigenvalues - 1)

  (1 / 2) * sum(((eigenvalues - 1) / kernel) * (((y - x %*% beta) ^ 2) / (sigma2 * kernel) - 1))
}


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

# lambda_sequence <- function(x, y, weights = NULL,
#                             lambda.factor = ifelse(nobs < nvars, 0.01, 1e-06),
#                             nlambda = 100, scale_x = F, center_y = F) {
#
#   # when scaling, first you center then you standardize
#   if (any(as.vector(weights) < 0)) stop("Weights must be positive")
#   np <- dim(x)
#   nobs <- as.integer(np[1])
#   nvars <- as.integer(np[2])
#
#   if (!is.null(weights) & length(as.vector(weights)) < nvars)
#     stop("You must provide weights for every column of x")
#
#   # scale the weights to sum to nvars
#   w <- if (is.null(weights)) rep(1, nvars) else as.vector(weights) / sum(as.vector(weights)) * nvars
#
#   sx <- if (scale_x) apply(x,2, function(i) scale(i, center = TRUE, scale = mysd(i))) else x
#   sy <- if (center_y) as.vector(scale(y, center = T, scale = F)) else as.vector(y)
#   lambda.max <- max(abs(colSums(sy * sx) / w)) / nrow(sx)
#
#   rev(exp(seq(log(lambda.factor * lambda.max), log(lambda.max), length.out = nlambda)))
# }





