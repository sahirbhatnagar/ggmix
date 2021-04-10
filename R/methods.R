#' Print Method for Objects of Class \code{ggmix_fit}
#'
#' @description print method for objects of class \code{ggmix_fit}
#' @param x an object of class objects of class \code{ggmix_fit}
#' @param digits significant digits in printout. Default: \code{max(3,
#'   getOption("digits") - 3)}
#' @param ... other arguments passed to \code{print}
#' @seealso \code{\link{ggmix}}
#' @return The call that produced the object \code{x} is printed, followed by a
#'   three-column matrix with columns \code{Df}, \code{\%Dev}, and and
#'   \code{Lambda}. The \code{Df} columns correspond to the number of nonzero
#'   coefficients including variance components. \code{\%dev} is the percent
#'   deviance explained (relative to the null deviance). \code{Lambda} is the
#'   sequence of converged tuning parameters.
#' @export
print.ggmix_fit <- function(x, ..., digits = max(3, getOption("digits") - 3)) {

  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(
    Df = x$result[, "Df"],
    `%Dev` = signif(x$result[, "%Dev"], digits),
    Lambda = signif(x$result[, "Lambda"], digits)
  ))
}

#' @export
#' @rdname print.ggmix_fit
print.ggmix_gic <- function(x, ..., digits = max(3, getOption("digits") - 3)) {
  xx <- x$ggmix_fit
  cat("\nCall: ", deparse(xx$call), "\n\n")
  print(cbind(
    Df = xx$result[, "Df"],
    `%Dev` = signif(xx$result[, "%Dev"], digits),
    Lambda = signif(xx$result[, "Lambda"], digits),
    GIC = signif(x$gic, digits)
  ))
}


#' Make predictions from a \code{ggmix_fit} object
#'
#' @description Similar to other predict methods, this functions predicts fitted
#'   values, coefficients and more from a fitted \code{ggmix_fit} object.
#' @param object Fitted \code{ggmix_fit} model object from the
#'   \code{\link{ggmix}} function
#' @param newx matrix of values for \code{x} at which predictions are to be
#'   made. Do not include the intercept. Must be a matrix. This argument is not
#'   used for \code{type = c("coefficients","nonzero","all")}. This matrix must
#'   have the same number of columns originally supplied to the
#'   \code{\link{ggmix}} fitting function.
#' @param s Value(s) of the penalty parameter \code{lambda} at which predictions
#'   are required. Default is the entire sequence used to create the model.
#' @param type Type of prediction required. Type \code{"link"} gives the fitted
#'   values \eqn{X \beta}. Type \code{"response"} is equivalent to type
#'   \code{"link"}. Type \code{"coefficients"} computes the coefficients at the
#'   requested values for \code{s} and returns the regression coefficients only,
#'   including the intercept. Type \code{"all"} returns both the regression
#'   coefficients and variance components at the requested value of \code{s}.
#'   Type \code{"nonzero"} returns a 1 column matrix of the the nonzero fixed
#'   effects, as well as variance components for each value of \code{s}. If more
#'   than one \code{s} is provided, then \code{"nonzero"} will return a list of
#'   1 column matrices. Default: "link"
#' @param covariance covariance between test and training individuals. if there
#'   are q testing individuals and N-q training individuals, then this
#'   covariance matrix is q x (N-q)
#' @return The object returned depends on type.
#' @method predict ggmix_fit
#' @details \code{s} is the new vector at which predictions are requested. If
#'   \code{s} is not in the lambda sequence used for fitting the model, the
#'   predict function will use linear interpolation to make predictions. The new
#'   values are interpolated using a fraction of predicted values from both left
#'   and right lambda indices. \code{coef(...)} is equivalent to
#'   \code{predict(ggmix_fit, type="coefficients",...)}. To get individual level
#'   predictions at each value of lambda, you must provide the lambda sequence
#'   to the s argument. You can pass either a ggmix_fit or ggmix_gic object. See
#'   examples for more details.
#' @examples
#' data("admixed")
#' fitlmm <- ggmix(x = admixed$xtrain, y = admixed$ytrain,
#'                 kinship = admixed$kin_train,
#'                 estimation = "full")
#' bicGGMIX <- gic(fitlmm,
#'                 an = log(length(admixed$ytrain)))
#' plot(bicGGMIX)
#' coef(bicGGMIX, s = "lambda.min")
#' yhat_test <- predict(bicGGMIX, s="lambda.min",
#'                      newx = admixed$xtest, type = "individual",
#'                      covariance = admixed$kin_test_train)
#' cor(yhat_test, admixed$ytest)
#' yhat_test_population <- predict(bicGGMIX, s="lambda.min",
#'                                 newx = admixed$xtest,
#'                                 type = "response")
#' @export
predict.ggmix_fit <- function(object, newx, s = NULL,
                              type = c(
                                "link", "response", "coefficients",
                                "all", "nonzero", "individual"), covariance, ...) {

  # if you use predict on a gic_fit object, then by the time the function gets here
  # s has been converted to a numeric

  type <- match.arg(type)

  if (missing(newx)) {
    if (!match(type, c("coefficients", "nonzero", "all"), FALSE)) {
      stop("You need to supply a value for 'newx' when type is link or response or individual")
    }
  }

  if (missing(covariance)) {
    if (type == "individual") {
      stop("You need to supply a value for 'covariance' when type is individual")
    }
  }

  # a0 <- t(as.matrix(object[["b0"]]))
  # rownames(a0) <- "(Intercept)"
  # nbeta <- rbind(a0, object[["beta"]])
  nall <- object[["coef"]] #first row is intercept,then betas, last two rows are eta and sigma2

  if (!is.null(s)) {
    vnames <- dimnames(nall)[[1]]
    dimnames(nall) <- list(NULL, NULL)
    lambda <- object[["lambda"]]
    lamlist <- lambda.interp(lambda, s)
    if (length(s) == 1) {
      nall <- nall[, lamlist$left, drop = FALSE] * lamlist$frac +
        nall[, lamlist$right, drop = FALSE] * (1 -
                                                 lamlist$frac)
    } else {
      nall <- nall[, lamlist$left, drop = FALSE] %*%
        diag(lamlist$frac) + nall[, lamlist$right,
                                  drop = FALSE
                                  ] %*% diag(1 - lamlist$frac)
    }
    dimnames(nall) <- list(vnames, paste(seq(along = s)))
  }

  if (type == "all") return(nall)

  nbeta <- nall[object[["cov_names"]], , drop = FALSE]

  if (type == "coefficients") return(nbeta)

  if (type == "nonzero") {
    nall.mat <- as.matrix(nall)
    if (length(s) == 1) {
      return(nall.mat[nonzeroCoef(nall.mat, bystep = TRUE)[[1]], , drop = FALSE])
    } else {
      nzs <- nonzeroCoef(nall.mat, bystep = TRUE)
      return(lapply(seq_along(nzs), function(i) nall.mat[nzs[[i]], i, drop = FALSE]))
    }
  }

  if (type == "link" | type == "response") {
    nfit <- as.matrix(cbind(1, newx) %*% nbeta) # this will result in a n x nlambda matrix!!!!!
    # The user must not input the first column as a intercept
    # once the rotation is done on the Xs and Y, we use them for fitting the function
    # but after that we dont use the rotated Xs or Y anymore. We use the original Xs and Ys for
    # prediction, residuals, ect.
    return(nfit)
  }


  if (type == "individual") {


    if (inherits(object, "lassofullrank")) {

      if (length(s) == 1) {
        # browser()
        eta <- nall["eta", 1]
        beta <- nall[object[["cov_names"]], 1, drop = FALSE]
        nfit <- as.matrix(cbind(1, newx) %*% nbeta)


        # see ranef.R
        return(
          as.vector(nfit) +
          bi_future_lassofullrank(
          eta = eta,
          beta = beta,
          eigenvalues = object[["ggmix_object"]][["D"]],
          eigenvectors = object[["ggmix_object"]][["U"]],
          x = object[["ggmix_object"]][["x"]],
          y = object[["ggmix_object"]][["y"]],
          covariance = covariance)
        )
      } else {

        nfit <- as.matrix(cbind(1, newx) %*% nbeta)

        bis <- lapply(seq_along(s), function(i) {
          # browser()
          eta <- nall["eta", i]
          sigma2 <- nall["sigma2", i]
          beta <- nall[object[["cov_names"]], i, drop = FALSE]

          nfit[, i] +
          bi_future_lassofullrank(
            eta = eta,
            beta = beta,
            eigenvalues = object[["ggmix_object"]][["D"]],
            eigenvectors = object[["ggmix_object"]][["U"]],
            x = object[["ggmix_object"]][["x"]], # these are the transformed x
            y = object[["ggmix_object"]][["y"]],
            covariance = covariance) # these are the transformed y
        })

        bisall <- do.call(cbind, bis)
        dimnames(bisall) <- list(rownames(object[["ggmix_object"]][["x"]]), paste(seq(along = s)))
        return(bisall)
      }
    } else {
      stop(strwrap("predict with type='individual' currently only implemented for lasso full rank"))
    }


    # The user must not input the first column as a intercept
    # once the rotation is done on the Xs and Y, we use them for fitting the function
    # but after that we dont use the rotated Xs or Y anymore. We use the original Xs and Ys for
    # prediction, residuals, ect.

  }


}

#' @rdname predict.ggmix_fit
#' @param ... additional arguments to pass to predict function
#' @export
coef.ggmix_fit <- function(object, s = NULL, type, ...) {

  if (missing(type)) {
  stats::predict(object, s = s, type = "coefficients", ...)
  } else {
    stats::predict(object, s = s, type = type, ...)
  }
}




#' @title Make predictions from a \code{ggmix_gic} object
#' @description This function makes predictions from a \code{ggmix_gic} object,
#'   using the stored "ggmix_fit" object, and the optimal value chosen for
#'   lambda using the gic.
#' @param object fitted \code{ggmix_gic} object
#' @inheritParams predict.ggmix_fit
#' @param s Value(s) of the penalty parameter \code{lambda} at which predictions
#'   are required. Default is the value \code{s="lambda.min"} can be used. If
#'   \code{s} is numeric, it is taken as the value(s) of \code{lambda} to be
#'   used.
#' @param ... other arguments passed to \code{\link{predict.ggmix_fit}}
#' @return The object returned depends the ... argument which is passed on to
#'   the predict method for \code{ggmix_fit} objects.
#' @details This function makes it easier to use the results of gic chosen model
#'   to make a prediction.
#' @rdname predict.ggmix_gic
#' @seealso \code{\link{predict.ggmix_fit}}
#' @export
predict.ggmix_gic <- function(object, newx, s = c("lambda.min"), ...) {
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
  stats::predict(object[["ggmix_fit"]], newx = newx, s = lambda, ...)
}


#' @inheritParams predict.ggmix_gic
#' @rdname predict.ggmix_gic
#' @export
coef.ggmix_gic <- function(object, s = c("lambda.min"), type, ...) {
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

  if (missing(type)) {
  stats::coef(object[["ggmix_fit"]], s = lambda, type = "coefficients", ...)
  } else {
    stats::coef(object[["ggmix_fit"]], s = lambda, type = type, ...)
  }
}






# not used under this line ------------------------------------------------


#' @param s lamda at which to predict the random effects. current option is only
#'   "lambda.min"
#'
#' @details For objects of class "gic.ggmix", this function returns the
#'   subject-specific random effect value for the model which minimizes the GIC
#'   using the maximum a posteriori principle
#'
#' @method ranef gic.ggmix
#' @rdname ranef
# ranef.gic.ggmix <- function(object, s = "lambda.min", ...) {
#
#   # object = res
#   # s = "lambda.min"
#   # ==================
#
#   if (s == "lambda.min") {
#     ranef(object = object$ggmix.fit, s = object$lambda.min, ...)
#   } else if (is.numeric(s)) {
#
#   }
# }



#' @param s index of tuning parameter. Must be a character and an element of
#'   "s1","s2",...."s100", where "s100" is the index of the last pair of tuning
#'   parameters. Default is \code{NULL}
#' @details For objects of class "ggmix", this function returns the
#'   subject-specific random effect value for the model which minimizes the GIC
#'   using the maximum a posteriori principle
#'
#' @method ranef ggmix
#' @rdname ranef
# ranef.ggmix <- function(object, new.x, new.u, new.d, s = NULL,
#                         type = c("fitted", "predicted")) {
#
#   # object = res$ggmix.fit
#   # s = c(0.5, 0.3, 0.1)
#   # type = "link"
#   # new.x = dat$x[,1:500]
#   # new.u = U
#   # new.d = Lambda
#   # type = "fitted"
#   # s = "lambda.min"
#   # ==================
#
#   type <- match.arg(type)
#
#   if (any(missing(new.x), missing(new.u), missing(new.d))) {
#     if (!match(type, c("fitted"), FALSE)) {
#       stop("You need to supply a value for 'new.x', 'new.u' and 'new.d'")
#     }
#   }
#
#   a0 <- t(as.matrix(object$b0))
#   eta <- as.matrix(object$eta)
#   sigma2 <- as.matrix(object$sigma2)
#   rownames(a0) <- "(Intercept)"
#   rownames(eta) <- "eta"
#   rownames(sigma2) <- "sigma2"
#   nbeta <- rbind(a0, object$beta, eta, sigma2)
#
#   if (!is.null(s)) {
#     vnames <- dimnames(nbeta)[[1]]
#     dimnames(nbeta) <- list(NULL, NULL)
#     lambda <- object$lambda
#     lamlist <- lambda.interp(lambda, s)
#     if (length(s) == 1) {
#       nbeta <- nbeta[, lamlist$left, drop = FALSE] * lamlist$frac +
#         nbeta[, lamlist$right, drop = FALSE] * (1 -
#           lamlist$frac)
#     } else {
#       nbeta <- nbeta[, lamlist$left, drop = FALSE] %*%
#         diag(lamlist$frac) + nbeta[, lamlist$right,
#           drop = FALSE
#         ] %*% diag(1 - lamlist$frac)
#     }
#     dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
#   }
#
#   if (type == "fitted") {
#     bis <- lapply(seq_along(s), function(i) {
#       eta_next <- nbeta["eta", i]
#       beta_next <- nbeta[c("(Intercept)", object$cov_names[-1]), i, drop = F]
#       bi(
#         eta = eta_next, beta = beta_next, eigenvalues = object$eigenvalues,
#         eigenvectors = object$u, x = object$utx, y = object$uty
#       )
#     })
#
#     bisall <- do.call(cbind, bis)
#     dim(bisall)
#     dimnames(bisall) <- list(rownames(object$x), paste(seq(along = s)))
#     return(bisall)
#   }
# }








#' @method random.effects gic.ggmix
#' @rdname ranef
# random.effects.gic.ggmix <- function(object, s = "lambda.min") {
#
#   # object = res
#   # s = "lambda.min"
#   # ==================
#
#   U <- object$ggmix.fit[["u"]]
#   estimates <- coef(object, s = s)
#   eta_next <- estimates["eta", ]
#   beta_next <- estimates[c("(Intercept)", object$ggmix.fit$cov_names[-1]), , drop = F]
#   eigenvalues <- object$ggmix.fit$eigenvalues
#
#   di <- 1 + eta_next * (eigenvalues - 1)
#   D_tilde_inv <- diag(1 / di)
#   bi <- as.vector(U %*% diag(1 / (1 / di + 1 / (eta_next * eigenvalues))) %*% t(U) %*%
#     U %*% D_tilde_inv %*% (object$ggmix.fit$uty - object$ggmix.fit$utx %*%
#       beta_next))
#   bi
# }
