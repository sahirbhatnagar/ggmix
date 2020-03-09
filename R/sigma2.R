#' Estimation of Sigma2 for Linear Mixed Model with Lasso Penalty
#'
#' @description \code{sigma2lasso} estimates the value of the sigma2 for the
#'   linear mixed model with lasso penalty
#'
#' @seealso \code{\link{ggmix}}
#' @param ggmix_object A ggmix_object object of class \code{lowrank} or
#'   \code{fullrank}
#' @param beta current estimate of the beta parameter including the intercept.
#'   this should be of length p+1, where p is the number of variables.
#' @param eta current estimate of the eta parameter
#' @param n number of observations
#' @param ... Extra parameters. Currently ignored.
#' @return A decreasing sequence of tuning parameters
#' @note There is a closed form solution for sigma^2, given beta and eta. This
#'   function isn't meant to be called directly by the user.
sigma2lasso <- function(ggmix_object, ...) UseMethod("sigma2lasso")

#' @rdname sigma2lasso
sigma2lasso.default <- function(ggmix_object, ...) {
  stop(strwrap("This function should be used with a ggmix object of class
               lowrank or fullrank"))
}

#' @rdname sigma2lasso
sigma2lasso.fullrank <- function(ggmix_object,
                                 ...,
                                 n, beta, eta) {
  x <- ggmix_object[["x"]]
  y <- ggmix_object[["y"]]
  eigenvalues <- ggmix_object[["D"]]

  (1 / n) * sum(((y - x %*% beta)^2) / (1 + eta * (eigenvalues - 1)))
}


# still need to create sigma2lasso.lowrank, sigma2gglasso.fullrank, sigma2gglasso.lowrank
