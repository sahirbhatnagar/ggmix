#' Generalised Information Criterion
#'
#' @description Calculates the generalised information criterion for each value
#'   of the tuning parameter lambda
#'
#' @seealso \code{\link{ggmix}}
#' @param ggmix_fit An object of class \code{ggmix_fit} which is outputted by
#'   the \code{\link{ggmix}} function
#' @param ... other parameters. currently ignored.
#' @details the generalised information criterion used for gaussian response is
#'   given by \deqn{-2 * loglikelihood(\hat{\Theta}) + an * df} where df is the
#'   number of non-zero estimated parameters, including variance components
#' @references Fan Y, Tang CY. Tuning parameter selection in high dimensional
#'   penalized likelihood. Journal of the Royal Statistical Society: Series B
#'   (Statistical Methodology). 2013 Jun 1;75(3):531-52.
#'
#'   Nishii R. Asymptotic properties of criteria for selection of variables in
#'   multiple regression. The Annals of Statistics. 1984;12(2):758-65.
#' @return an object with S3 class \code{"ggmix_gic"}, \code{"ggmix_fit"},
#'   \code{"*"} and \code{"**"} where \code{"*"} is "lasso" or "gglasso" and
#'   \code{"**"} is fullrank or lowrank. Results are provided for converged
#'   values of lambda only. \describe{\item{ggmix_fit}{the ggmix_fit
#'   object}\item{lambda}{the sequence of converged tuning parameters}
#'   \item{nzero}{the number of non-zero estimated coefficients including the 2
#'   variance parameters which are not penalized and therefore always
#'   included}\item{gic}{gic value. a numeric vector with length equal to
#'   \code{length(lambda)}} \item{lambda.min.name}{a character corresponding to
#'   the name of the tuning parameter lambda which minimizes the
#'   gic}\item{lambda.min}{the value of lambda which minimizes the gic}}
#' @export
gic <- function(ggmix_fit, ...) UseMethod("gic")

#' @rdname gic
#' @export
gic.default <- function(ggmix_fit, ...) {
  stop(strwrap("This function should be used with an object of class
               ggmix_fit"))
}


#' @param an numeric, the penalty per parameter to be used; the default is an =
#'   log(log(n))*log(p) where n is the number of subjects and p is the number of
#'   parameters
#' @rdname gic
#' @export
gic.ggmix_fit <- function(ggmix_fit,
                          ...,
                          an = log(log(n)) * log(p)) {
  n <- ggmix_fit[["n_design"]]
  p <- ggmix_fit[["p_design"]]
  df <- ggmix_fit$result[, "Df"]
  model_loglik <- ggmix_fit$result[, "loglik"]

  model_gic <- -2 * model_loglik + an * df

  out <- list(
    ggmix_fit = ggmix_fit, # used in predict.ggmix_gic function
    lambda = ggmix_fit[["lambda"]],
    nzero = df,
    gic = model_gic,
    lambda.min.name = names(which.min(model_gic)),
    lambda.min = ggmix_fit$result[names(which.min(model_gic)), "Lambda"]
  )
  obj <- c(out)
  class(obj) <- c("ggmix_gic", attr(ggmix_fit, "class"))
  return(obj)
}
