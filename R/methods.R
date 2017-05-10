#' Print Method for penfam function
#'
#' @description print method for penfam function
#' @export

print.penfam <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(Df = x$x[,"Df"],
              `%Dev` = signif(x$x[,"%Dev"], digits),
              Lambda = signif(x$x[,"Lambda"], digits),
              BIC = signif(x$x[,"BIC"], digits)))
}



#' Get coefficients from a "penfam" object
#'
#' @rdname predict.penfam
#' @export

coef.penfam <- function(object, s = NULL,  ...) {
  predict(object, s = s, type = "coefficients", ...)
}




#' Make predictions from a penfam object
#'
#' @description this function only works for tuning parameter values defined by
#'   the shim_multiple_faster function. The interpolation feature is not working
#'   yet
#' @param s index of tuning parameter. Must be a character and an element of
#'   "s1","s2",...."s100", where "s100" is the index of the last pair of tuning
#'   parameters. Default is \code{NULL}
#' @param object Fitted shim model object
#' @param type Type of prediction required. Type "link" gives the fitted values.
#'   Type "response" for "gaussian" type "response" is equivalent to type
#'   "link". Type "coefficients" computes the coefficients at the requested
#'   values for s. Type "nonzero" returns a list of the indices of the nonzero
#'   coefficients for each value of s.
#' @export

predict.penfam <- function(object, newx, s = NULL,
                           type = c("link", "response", "coefficients",
                                    "nonzero")) {

  # object = res
  # s = "s88"
  # ==================

  type = match.arg(type)

  if (missing(newx)) {
    if (!match(type, c("coefficients", "nonzero"), FALSE))
      stop("You need to supply a value for 'newx'")
  }

  a0 = t(as.matrix(object$b0))
  rownames(a0) = "(Intercept)"
  # this includes tuning parameters pairs that didnt converge
  nbeta = rbind(a0, object$beta, object$eta, object$sigma2)
  # nbeta@Dimnames <- list(X = c("(Intercept)", object$cov_names, "eta", "sigma2"),
  #                        Y = paste0("s",seq_len(object$nlambda)))

  # dimnames(nbeta)[[1]] <- list(c("(Intercept)", object$cov_names, "eta", "sigma2"))
  # dimnames(nbeta)[[2]] <- paste0("s",seq_len(object$nlambda))

  # this is the default returned by coef.shim i.e. any object of class shim
  # it will return all tuning parameters (including those that didnt converge)
  if (type == "coefficients" && is.null(s)) {
    return(nbeta)
  }

  if (type == "coefficients" && !is.null(s)) {
    return(nbeta[ , s, drop = F])
  }

  if (type == "nonzero" && is.null(s)) {
    return(glmnet::nonzeroCoef(nbeta[, , drop = FALSE], bystep = TRUE))
  }

  if (type == "nonzero" && !is.null(s)) {
    if (length(s) > 1) stop("type=nonzero is only valid for a single value for s")
    return(nbeta[glmnet::nonzeroCoef(nbeta[, s , drop = FALSE], bystep = TRUE)[[1]], s, drop = F])
  }

  if (inherits(newx, "sparseMatrix")) {
    newx = as(newx, "dgCMatrix")
  }

  # this is used by the cv_lspath function to calculate predicted values
  # which will subsequently be used for calculating MSE for each fold
  if (type == "link") {

    nfit = as.matrix(cbind2(1, newx) %*% nbeta)

    return(nfit)
  }

}




plot.penfam <- function(x, type=c("coef","BIC"), xvar=c("norm","lambda","dev"),
                        label=FALSE, sign.lambda = 1, ...){
  xvar=match.arg(xvar)
  type=match.arg(type)

  if (type == "coef") {
    plotCoef(x$beta, lambda=drop(x$x[,"Lambda"]),
             df=drop(x$x[,"Df"]), dev=drop(x$x[,"%Dev"]),
             label=label, xvar=xvar,...)
  }

  if (type == "BIC"){
    plotBIC(object = x$x, sign.lambda=sign.lambda, lambda.min = x$lambda_min, ...)
  }
}

