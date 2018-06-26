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
  print(cbind(Df = x$result[,"Df"],
              `%Dev` = signif(x$result[,"%Dev"], digits),
              Lambda = signif(x$result[,"Lambda"], digits)))
}



#' Get coefficients from a "ggmix" object
#'
#' @param object object of class ggmix
#' @param s lambda value
#' @param which a character of either "beta", "eta", or "sigma2" if the user
#'   wants only th
#' @param ... additional arguments to pass to predict function
#'
#' @rdname predict.ggmix
#' @export
coef.ggmix <- function(object, s = NULL, which = NULL, ...) {
  if (is.null(which)) predict(object, s = s, type = "coefficients", ...) else
    predict(object, s = s, type = which, ...)
}

#' @rdname predict.ggmix
coef.gic.ggmix <- function(object, s = "lambda.min", ...) {
  if (is.numeric(s))
    lambda = s
  else if (is.character(s)) {
    s = match.arg(s)
    lambda = object[[s]]
  }
  else stop("Invalid form for s")
  coef(object$ggmix.fit, s = lambda, ...)
}


#' Make predictions from a ggmix object
#'
#' @description this function only works for tuning parameter values defined by
#'   the shim_multiple_faster function. The interpolation feature is not working
#'   yet
#' @param s index of tuning parameter. Must be a character and an element of
#'   "s1","s2",...."s100", where "s100" is the index of the last pair of tuning
#'   parameters. Default is \code{NULL}
#' @param object Fitted ggmix model object
#' @param type Type of prediction required. Type "link" gives the fitted values
#'   given by \deqn{X\beta + b} where b is a subject-specific random effect. For
#'   "gaussian" type "response" is equivalent to type "link". Type
#'   "coefficients" computes the coefficients at the requested values for s.
#'   Type "nonzero" returns a list of the indices of the nonzero coefficients
#'   for each value of s.
#' @export
predict.ggmix_fit <- function(object, new.x, new.u, new.d, s = NULL,
                              type = c("link", "response", "coefficients","ranef",
                                       "nonzero", "beta", "eta", "sigma2")) {

  # object = res$ggmix.fit
  # s = c(0.1,0.2,0.3)
  # s=NULL
  # type = "link"
  # newx = dat$x[,1:500]
  # new.u = U
  # new.d = Lambda
  # ==================


  # need to think about how to project newx... cant use eigenvectors of training sample beacuse thats for a sample i
  # probably need to recalculate eigenvectors for newx. will not re-calculate. just use original x

  type = match.arg(type)

  if (any(missing(new.x), missing(new.u), missing(new.d))) {
    if (!match(type, c("coefficients", "nonzero"), FALSE))
      stop("You need to supply a value for 'new.x', 'new.u' and 'new.d'")
  }

  a0 = t(as.matrix(object$b0))
  rownames(a0) = "(Intercept)"
  nbeta = rbind2(a0, object$beta)


  if (!is.null(s)) {
    vnames <- dimnames(nbeta)[[1]]
    dimnames(nbeta) <- list(NULL, NULL)
    lambda <- object$lambda
    lamlist <- glmnet::lambda.interp(lambda, s)
    if (length(s) == 1) {
      nbeta = nbeta[, lamlist$left, drop = FALSE] * lamlist$frac +
        nbeta[, lamlist$right, drop = FALSE] * (1 -
                                                  lamlist$frac)
    } else {
      nbeta = nbeta[, lamlist$left, drop = FALSE] %*%
        diag(lamlist$frac) + nbeta[, lamlist$right,
                                   drop = FALSE] %*% diag(1 - lamlist$frac)
    }
    dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
  }


  if (type == "nonzero")
    return(glmnet::nonzeroCoef(nbeta[-1, , drop = FALSE], bystep = TRUE))

  if (inherits(newx, "sparseMatrix"))
    newx = as(newx, "dgCMatrix")


  # this is used by the cv_lspath function to calculate predicted values
  # which will subsequently be used for calculating MSE for each fold
  if (type == "link") {

    nfit = as.matrix(cbind2(1, newx) %*% nbeta) # this will result in a n x nlambda matrix!!!!!
    ranef.ggmix(object)
    return(nfit)
  }




  if (type == "coefficients" && is.null(s)) {
    return(nbeta)
  }

  if (type == "coefficients" && !is.null(s)) {
    return(nbeta[ , s, drop = F])
  }

  if (type == "beta" && is.null(s)) {
    return(object$beta)
  }

  if (type == "beta" && !is.null(s)) {
    return(object$beta[ , s, drop = F])
  }


  if (type == "varparams" && !is.null(s)) {
    return(object$beta[ , s, drop = F])
  }


  if (type == "response" && is.null(s)) {
    return(object$predicted)
  }

  if (type == "response" && !is.null(s)) {
    return(object$predicted[, s, drop = F])
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


}


#' Make predictions from a "gic.ggmix" object.
#'
#' @param object a fitted \code{\link{gic.ggmix}} object
#' @param newx matrix of new values for x at which predictions are to be made
#' @param s value(s) of the penalty parameter lambda at which predictions are
#'   required. Default is the value s="lambda.min" stored on the
#'   \code{\link{gic.ggmix}} object. If s is numeric, it is taken as the
#'   value(s) of lambda to be used.
#' @param ... not used. Other arguments to predict
#'
#' @description This function makes predictions from a ggmix model with GIC,
#'   using the stored "ggmix.fit" object, and the optimal value chosen for lambda
#'   based on the minimum GIC.
#' @details This function makes it easier to use the results of cross-validation
#'   to make a prediction
#' @method predict gic.ggmix
#' @export
predict.gic.ggmix <- function(object, newx, s = c("lambda.1se",
                                                   "lambda.min"), ...) {
  if (is.numeric(s))
    lambda <- s else if (is.character(s)) {
      s <- match.arg(s)
      lambda <- object[[s]]
    } else stop("Invalid form for s")
  predict(object$ggmix.fit, newx, s = lambda, ...)
}



plot.ggmix_fit <- function(x,
                           type = c("coef","QQranef","QQresid", "predicted", "Tukey-Anscombe"),
                           xvar=c("norm","lambda","dev"), s = x$lambda_min,
                           label=FALSE, sign.lambda = 1, ...){
  xvar <- match.arg(xvar)
  type <- match.arg(type, several.ok = TRUE)

  if (any(type == "coef")) {
    plotCoef(x$beta, lambda=drop(x$result[,"Lambda"]),
             df=drop(x$result[,"Df"]), dev=drop(x$result[,"Deviance"]),
             label=label, xvar=xvar,...)
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
         main = sprintf("Observed vs. Predicted response
                        R2 = %g", cor(x$predicted[, s], drop(x$y))))
    abline(a = 0, b = 1, col = "red")
  }


  if (any(type == "Tukey-Anscombe")){
    plot(x$fitted[, s], x$residuals[, s], main = "Tukey-Anscombe Plot",
         xlab = "fitted values (XB)", ylab = "residuals")
    abline(h = 0, col = "red")
  }


}




#' @param s lamda at which to predict the random effects. current option is only
#'   "lambda.min"
#'
#' @details For objects of class "gic.ggmix", this function returns the
#'   subject-specific random effect value for the model which minimizes the GIC
#'   using the maximum a posteriori principle
#'
#' @method ranef gic.ggmix
#' @rdname ranef
#' @export
ranef.gic.ggmix <- function(object, s = "lambda.min", ...) {

  # object = res
  # s = "lambda.min"
  #==================

  if (s == "lambda.min") {
  ranef(object = object$ggmix.fit, s = object$lambda.min, ...)
  } else if(is.numeric(s)) {

  }

}



#' @param s index of tuning parameter. Must be a character and an element of
#'   "s1","s2",...."s100", where "s100" is the index of the last pair of tuning
#'   parameters. Default is \code{NULL}
#' @details For objects of class "ggmix", this function returns the
#'   subject-specific random effect value for the model which minimizes the GIC
#'   using the maximum a posteriori principle
#'
#' @method ranef ggmix
#' @rdname ranef
#' @export
ranef.ggmix <- function(object, new.x, new.u, new.d, s = NULL,
                         type = c("fitted", "predicted")) {

  # object = res$ggmix.fit
  # s = c(0.5, 0.3, 0.1)
  # type = "link"
  # new.x = dat$x[,1:500]
  # new.u = U
  # new.d = Lambda
  # type = "fitted"
  # s = "lambda.min"
  #==================

  type = match.arg(type)

  if (any(missing(new.x), missing(new.u), missing(new.d))) {
    if (!match(type, c("fitted"), FALSE))
      stop("You need to supply a value for 'new.x', 'new.u' and 'new.d'")
  }

  a0 <- t(as.matrix(object$b0))
  eta <- as.matrix(object$eta)
  sigma2 <- as.matrix(object$sigma2)
  rownames(a0) <- "(Intercept)"
  rownames(eta) <- "eta"
  rownames(sigma2) <- "sigma2"
  nbeta = rbind(a0, object$beta, eta, sigma2)

  if (!is.null(s)) {
    vnames <- dimnames(nbeta)[[1]]
    dimnames(nbeta) <- list(NULL, NULL)
    lambda <- object$lambda
    lamlist <- glmnet::lambda.interp(lambda, s)
    if (length(s) == 1) {
      nbeta = nbeta[, lamlist$left, drop = FALSE] * lamlist$frac +
        nbeta[, lamlist$right, drop = FALSE] * (1 -
                                                  lamlist$frac)
    } else {
      nbeta = nbeta[, lamlist$left, drop = FALSE] %*%
        diag(lamlist$frac) + nbeta[, lamlist$right,
                                   drop = FALSE] %*% diag(1 - lamlist$frac)
    }
    dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
  }

  if (type == "fitted") {

    bis <- lapply(seq_along(s), function(i) {
      eta_next <- nbeta["eta",i]
      beta_next <- nbeta[c("(Intercept)", object$cov_names[-1]),i,drop=F]
      bi(eta = eta_next, beta = beta_next, eigenvalues = object$eigenvalues,
         eigenvectors = object$u, x = object$utx, y = object$uty)
    })

    bisall <- do.call(cbind, bis)
    dim(bisall)
    dimnames(bisall) <- list(rownames(object$x), paste(seq(along = s)))
    return(bisall)

  }


}




bi <- function(eta, beta, eigenvalues, eigenvectors, x, y){
  di <- 1 + eta * (eigenvalues - 1)
  D_tilde_inv <- diag(1 / di)
  as.vector(eigenvectors %*% diag(1 / (1/di + 1/(eta * eigenvalues))) %*%
              t(eigenvectors) %*% eigenvectors %*% D_tilde_inv %*% (y - x %*% beta))
}



#' @method random.effects gic.ggmix
#' @rdname ranef
#' @export
random.effects.gic.ggmix <- function(object, s = "lambda.min") {

  # object = res
  # s = "lambda.min"
  #==================

  U <- object$ggmix.fit[["u"]]
  estimates <- coef(object, s = s)
  eta_next <- estimates["eta",]
  beta_next <- estimates[c("(Intercept)",object$ggmix.fit$cov_names[-1]),,drop=F]
  eigenvalues <- object$ggmix.fit$eigenvalues

  di <- 1 + eta_next * (eigenvalues - 1)
  D_tilde_inv <- diag(1 / di)
  bi <- as.vector(U %*% diag(1 / (1/di + 1/(eta_next*eigenvalues))) %*% t(U) %*%
                    U %*% D_tilde_inv %*% (object$ggmix.fit$uty - object$ggmix.fit$utx %*%
                                             beta_next))
  bi

}
