#' Print Method for penfam function
#'
#' @description print method for penfam function
#' @export

print.penfam <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(Df = x$result[,"Df"],
              `Deviance` = signif(x$result[,"Deviance"], digits),
              Lambda = signif(x$result[,"Lambda"], digits),
              BIC = signif(x$result[,"BIC"], digits)))
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




plot.penfam <- function(x, type = c("coef","BIC", "QQranef","QQresid", "predicted", "Tukey-Anscombe"),
                        xvar=c("norm","lambda","dev"), s = x$lambda_min,
                        label=FALSE, sign.lambda = 1, ...){
  xvar <- match.arg(xvar)
  type <- match.arg(type, several.ok = TRUE)

  if (any(type == "coef")) {
    plotCoef(x$beta, lambda=drop(x$result[,"Lambda"]),
             df=drop(x$result[,"Df"]), dev=drop(x$result[,"Deviance"]),
             label=label, xvar=xvar,...)
  }

  if (any(type == "BIC")){
    plotBIC(object = x$result, sign.lambda=sign.lambda, lambda.min = x$lambda_min, ...)
  }


  if (any(type == "QQranef")){
    if (s %ni% rownames(x$result)) stop("value for s not in lambda sequence")
    qqnorm(x$randomeff[, s], main = sprintf("QQ-Plot of the random effects at lambda = %.2f", x$result[s,"Lambda"]))
    qqline(x$randomeff[, s], col = "red")
  }

  if (any(type == "QQresid")){
    if (s %ni% rownames(x$result)) stop("value for s not in lambda sequence")
    qqnorm(x$residuals[, s], main = sprintf("QQ-Plot of the residuals at lambda = %.2f", x$result[s,"Lambda"]))
    qqline(x$residuals[, s], col = "red")
  }

  if (any(type == "predicted")){
    if (s %ni% rownames(x$result)) stop("value for s not in lambda sequence")
    plot(x$predicted[, s], drop(x$y), xlab = "predicted response (XB + b)", ylab = "observed response",
         main = "Observed vs. Predicted response")
    abline(a = 0, b = 1, col = "red")
  }

  if (any(type == "predicted")){
    if (s %ni% rownames(x$result)) stop("value for s not in lambda sequence")
    plot(x$predicted[, s], drop(x$y), xlab = "predicted response (XB + b)", ylab = "observed response",
         main = "Observed vs. Predicted response")
    abline(a = 0, b = 1, col = "red")
  }

  if (any(type == "Tukey-Anscombe")){
    plot(x$fitted[, s], x$residuals[, s], main = "Tukey-Anscombe Plot",
         xlab = "fitted values (XB)", ylab = "residuals")
    abline(h = 0, col = "red")
  }


}

