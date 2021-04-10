#' Extract Random Effects
#'
#' @description Generic function for extracting the random effects. This is the
#'   same generic (and same name) defined in the nlme and lme4 package.
#'
#' @return a numeric vector of length equal to the number of observations of
#'   subject-specific random effects
#'
#' @seealso \code{\link{gic}}
#' @param object any fitted model object from which random effects estimates can
#'   be extracted. Currently supports "ggmix_gic" objects outputted by the
#'   \code{\link{gic}} function
#' @param ... other parameters. currently ignored
#' @details For objects of class "ggmix_gic", this function returns the
#'   subject-specific random effect value for the model which minimizes the GIC
#'   using the maximum a posteriori principle
#' @examples
#' data("admixed")
#' fit <- ggmix(x = admixed$xtrain, y = admixed$ytrain,
#'             kinship = admixed$kin_train)
#' gicfit <- gic(fit)
#' # random effect at selected value of lambda
#' plot(ggmix::ranef(gicfit))
#' # random effects at specific values of lambda
#' head(ggmix::ranef(gicfit, s = c(0.1,0.2)))
#' @export
ranef <- function(object, ...) {
  UseMethod("ranef")
}

#' @rdname ranef
#' @export
random.effects <- function(object, ...) {
  UseMethod("ranef")
}

#' @rdname ranef
#' @method random.effects default
#' @export
random.effects.default <- function(object, ...) {
  stop(strwrap("This function should be used with an object of class
               ggmix_gic"))
}

#' @rdname ranef
#' @method ranef default
#' @export
ranef.default <- function(object, ...) {
  stop(strwrap("This function should be used with an object of class
               ggmix_gic"))
}


#' @rdname ranef
#' @method ranef ggmix_gic
#' @inheritParams predict.ggmix_gic
#' @export
ranef.ggmix_gic <- function(object, s = "lambda.min", ...) {

  if (is.numeric(s)) {
    lambda <- s
  } else
    if (is.character(s)) {
      s <- match.arg(s)
      lambda <- object[[s]]
    }
  else {
    stop("Invalid form for s")
  }

  if (inherits(object, "lassofullrank")) {

    nall <- stats::coef(object, s = lambda, type = "all")

    if (length(lambda) == 1) {
      # browser()
      eta <- nall["eta", 1]
      beta <- nall[object[["ggmix_fit"]][["cov_names"]], 1, drop = FALSE]

      return(bi_lassofullrank(
        eta = eta,
        beta = beta,
        eigenvalues = object[["ggmix_fit"]][["ggmix_object"]][["D"]],
        eigenvectors = object[["ggmix_fit"]][["ggmix_object"]][["U"]],
        x = object[["ggmix_fit"]][["ggmix_object"]][["x"]],
        y = object[["ggmix_fit"]][["ggmix_object"]][["y"]])
      )
    } else {

      nall <- stats::coef(object, s = lambda, type = "all")
      bis <- lapply(seq_along(s), function(i) {
        eta <- nall["eta", i]
        beta <- nall[object[["ggmix_fit"]][["cov_names"]], i, drop = FALSE]
        bi_lassofullrank(
          eta = eta,
          beta = beta,
          eigenvalues = object[["ggmix_fit"]][["ggmix_object"]][["D"]],
          eigenvectors = object[["ggmix_fit"]][["ggmix_object"]][["U"]],
          x = object[["ggmix_fit"]][["ggmix_object"]][["x"]], # these are the transformed x
          y = object[["ggmix_fit"]][["ggmix_object"]][["y"]]) # these are the transformed y
      })

      bisall <- do.call(cbind, bis)
      dimnames(bisall) <- list(rownames(object[["ggmix_fit"]][["ggmix_object"]][["x"]]), paste(seq(along = s)))
      return(bisall)
    }
  } else {
    stop(strwrap("ranef currently only implemented for lasso full rank"))
  }
}


bi_lassofullrank <- function(eta, beta, eigenvalues, eigenvectors, x, y) {
  D_inv <- diag(1 / eigenvalues)
  p1 <- solve((diag(length(y)) + (1 / eta) * (eigenvectors %*% D_inv %*% t(eigenvectors))))
  as.vector( p1 %*% (y - x %*% beta))
}


# this is for future observations used by predict.ggmix_fit when type="individual"
# this contains the random part only, i.e. the 2nd part of eq 35 in manuscript section 3.7
bi_future_lassofullrank <- function(eta, beta, eigenvalues, eigenvectors, x, y, covariance) {
  di <- 1 + eta * (eigenvalues - 1)
  D_tilde_inv <- diag(1 / di)
  (eta) * as.vector(covariance %*% eigenvectors %*% D_tilde_inv %*% (y - x %*% beta))
}
