#' Estimation of Log-likelihood for Linear Mixed Model with Lasso Penalty
#'
#' @description \code{sigma2lasso} estimates the value of the sigma2 for the
#'   linear mixed model with lasso penalty
#'
#' @seealso \code{\link{ggmix}}
#' @param ggmix_object A ggmix_object object of class \code{lowrank} or
#'   \code{fullrank}
#' @inheritParams ggmix
#' @param beta current estimate of the beta parameter including the intercept.
#'   this should be of length p+1, where p is the number of variables.
#' @param eta current estimate of the eta parameter
#' @param sigma2 current estimate of the sigma2 parameter
#' @param nt total number of observations
#' @param ... Extra parameters. Currently ignored.
#' @return A decreasing sequence of tuning parameters
#' @note This function isn't meant to be called directly by the user.
#' @seealso \code{\link{ggmix_data_object}}
logliklasso <- function(ggmix_object, ...) UseMethod("logliklasso")

#' @rdname logliklasso
logliklasso.default <- function(ggmix_object, ...) {
  stop(strwrap("This function should be used with a ggmix object of class
               lowrank or fullrank"))
}

#' @rdname logliklasso
logliklasso.fullrank <- function(ggmix_object,
                                 ...,
                                 eta, sigma2, beta, nt, x, y) {

  # this returns the log-likelihood

  # this allows the user to input the x and y if they want else
  # the default is from the ggmix_object
  # this flexibility is used to calculate the saturated and intercept loglik
  # in the lmmlasso function
  if (missing(x)) {
    x <- ggmix_object[["x"]]
  }

  if (missing(y)) {
    y <- ggmix_object[["y"]]
  }

  eigenvalues <- ggmix_object[["D"]]

  kernel <- 1 + eta * (eigenvalues - 1)

  -1 * (
    (nt / 2) * log(2 * pi) +
      (nt / 2) * log(sigma2) +
      0.5 * sum(log(kernel)) +
      (1 / (2 * sigma2)) * sum((y - x %*% beta)^2 / kernel)
  )
}


# still need to create logliklasso.lowrank, loglikgglasso.fullrank, loglikgglasso.lowrank
