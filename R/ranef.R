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

    nall <- coef(object, s = lambda, type = "all")

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

      nall <- coef(object, s = lambda, type = "all")
      bis <- lapply(seq_along(s), function(i) {
        eta <- nall["eta", i]
        beta <- nall[object[["ggmix_fit"]][["cov_names"]], i, drop = FALSE]
        bi_lassofullrank(
          eta = eta,
          beta = beta,
          eigenvalues = object[["ggmix_fit"]][["ggmix_object"]][["D"]],
          eigenvectors = object[["ggmix_fit"]][["ggmix_object"]][["U"]],
          x = object[["ggmix_fit"]][["ggmix_object"]][["x"]],
          y = object[["ggmix_fit"]][["ggmix_object"]][["y"]])
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
  di <- 1 + eta * (eigenvalues - 1)
  D_tilde_inv <- diag(1 / di)
  as.vector(eigenvectors %*% diag(1 / (1 / di + 1 / (eta * eigenvalues))) %*%
              t(eigenvectors) %*% eigenvectors %*% D_tilde_inv %*% (y - x %*% beta))
}
