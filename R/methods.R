#' Print Method for penfam function
#'
#' @description print method for penfam function
#' @export

print.penfam <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(Df = x$result[,"Df"],
              `%Dev` = signif(x$result[,"%Dev"], digits),
              Lambda = signif(x$result[,"Lambda"], digits)))
}



#' Get coefficients from a "penfam" object
#'
#' @param object object of class penfam
#' @param s lambda value
#' @param which a character of either "beta", "eta", or "sigma2" if the user wants only th
#' @param ... additional arguments to pass to predict function
#'
#' @rdname predict.penfam
#' @export

coef.penfam <- function(object, s = NULL, which = NULL, ...) {
  if(is.null(which)) predict(object, s = s, type = "coefficients", ...) else
    predict(object, s = s, type = which, ...)
}

#' @rdname predict.penfam
coef.gic.penfam <- function (object, s = "lambda.min", ...) {
  if (is.numeric(s))
    lambda = s
  else if (is.character(s)) {
    s = match.arg(s)
    lambda = object[[s]]
  }
  else stop("Invalid form for s")
  coef(object$penfam.fit, s = lambda, ...)
}







#' Make predictions from a penfam object
#'
#' @description this function only works for tuning parameter values defined by
#'   the shim_multiple_faster function. The interpolation feature is not working
#'   yet
#' @param s index of tuning parameter. Must be a character and an element of
#'   "s1","s2",...."s100", where "s100" is the index of the last pair of tuning
#'   parameters. Default is \code{NULL}
#' @param object Fitted penfam model object
#' @param type Type of prediction required. Type "link" gives the fitted values
#'   given by \deqn{X\beta + b} where b is a subject-specific random effect. For
#'   "gaussian" type "response" is equivalent to type "link". Type
#'   "coefficients" computes the coefficients at the requested values for s.
#'   Type "nonzero" returns a list of the indices of the nonzero coefficients
#'   for each value of s.
#' @export

predict.penfam <- function(object, new.x, new.u, new.d, s = NULL,
                           type = c("link", "response", "coefficients","ranef",
                                    "nonzero", "beta", "eta", "sigma2")) {

  object = res$penfam.fit
  s = c(0.1,0.2,0.3)
  s=NULL
  type = "link"
  newx = dat$x[,1:500]
  new.u = U
  new.d = Lambda
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
    ranef.penfam(object)
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


#' Make predictions from a "gic.penfam" object.
#'
#' @param object a fitted \code{\link{gic.penfam}} object
#' @param newx matrix of new values for x at which predictions are to be made
#' @param s value(s) of the penalty parameter lambda at which predictions are
#'   required. Default is the value s="lambda.min" stored on the
#'   \code{\link{gic.penfam}} object. If s is numeric, it is taken as the
#'   value(s) of lambda to be used.
#' @param ... not used. Other arguments to predict
#'
#' @description This function makes predictions from a penfam model with GIC,
#'   using the stored "penfam.fit" object, and the optimal value chosen for lambda
#'   based on the minimum GIC.
#' @details This function makes it easier to use the results of cross-validation
#'   to make a prediction
#' @method predict gic.penfam
#' @export
predict.gic.penfam <- function(object, newx, s = c("lambda.1se",
                                                   "lambda.min"), ...) {
  if (is.numeric(s))
    lambda <- s else if (is.character(s)) {
      s <- match.arg(s)
      lambda <- object[[s]]
    } else stop("Invalid form for s")
  predict(object$penfam.fit, newx, s = lambda, ...)
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

#' Plot the Generalised Information Criteria curve produced by gic.penfam
#'
#' @param x fitted "gic.penfam" object
#' @param sign.lambda Either plot against log(lambda) (default) or its negative
#'   if sign.lambda=-1
#' @param ... Other graphical parameters to plot
#' @details A plot is produced, and nothing is returned.
#' @seealso \code{\link{penfam}} and \code{\link{gic.penfam}}
#' @description Plots the Generalised Information Criteria curve, as a function
#'   of the lambda values used
#' @export
plot.gic.penfam <- function(x, sign.lambda = 1, ...) {

  # x <- res
  # sign.lambda = 1
  # lambda.min = res$lambda_min
  # ===============

  bicobj = x
  xlab="log(Lambda)"
  if(sign.lambda<0) xlab=paste("-",xlab,sep="")
  plot.args=list(x=sign.lambda*log(bicobj$lambda),
                 y=bicobj$bic,
                 ylim=range(bicobj$bic),
                 xlab=xlab,
                 ylab="BIC", type="n")
  new.args=list(...)
  if (length(new.args)) plot.args[names(new.args)]=new.args
  do.call("plot",plot.args)
  points(sign.lambda*log(bicobj$lambda),
         bicobj$bic,pch=20,col="red")
  axis(side=3,at=sign.lambda*log(bicobj$lambda),
       labels=paste(bicobj$nzero), tick=FALSE, line=0)
  abline(v=sign.lambda*log(bicobj$lambda.min.value),lty=3)
  # abline(v=sign.lambda*log(.1605),lty=3)
  # abline(v=sign.lambda*log(cvobj$lambda.1se),lty=3)
  # invisible()
}


#' @param s lamda at which to predict the random effects. current option is only
#'   "lambda.min"
#'
#' @details For objects of class "gic.penfam", this function returns the
#'   subject-specific random effect value for the model which minimizes the GIC
#'   using the maximum a posteriori principle
#'
#' @method ranef gic.penfam
#' @rdname ranef
#' @export
ranef.gic.penfam <- function(object, s = "lambda.min", ...) {

  # object = res
  # s = "lambda.min"
  #==================

  if (s == "lambda.min") {
  ranef(object = object$penfam.fit, s = object$lambda.min, ...)
  } else if(is.numeric(s)) {

  }

}



#' @param s index of tuning parameter. Must be a character and an element of
#'   "s1","s2",...."s100", where "s100" is the index of the last pair of tuning
#'   parameters. Default is \code{NULL}
#' @details For objects of class "penfam", this function returns the
#'   subject-specific random effect value for the model which minimizes the GIC
#'   using the maximum a posteriori principle
#'
#' @method ranef penfam
#' @rdname ranef
#' @export
ranef.penfam <- function(object, new.x, new.u, new.d, s = NULL,
                         type = c("fitted", "predicted")) {

  # object = res$penfam.fit
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
  as.vector(eigenvectors %*% diag(1 / (1/di + 1/(eta * eigenvalues))) %*% t(eigenvectors) %*% eigenvectors %*% D_tilde_inv %*% (y - x %*% beta))
}



#' @method random.effects gic.penfam
#' @rdname ranef
#' @export
random.effects.gic.penfam <- function(object, s = "lambda.min") {

  # object = res
  # s = "lambda.min"
  #==================

  U <- object$penfam.fit[["u"]]
  estimates <- coef(object, s = s)
  eta_next <- estimates["eta",]
  beta_next <- estimates[c("(Intercept)",object$penfam.fit$cov_names[-1]),,drop=F]
  eigenvalues <- object$penfam.fit$eigenvalues

  di <- 1 + eta_next * (eigenvalues - 1)
  D_tilde_inv <- diag(1 / di)
  bi <- as.vector(U %*% diag(1 / (1/di + 1/(eta_next*eigenvalues))) %*% t(U) %*% U %*% D_tilde_inv %*% (object$penfam.fit$uty - object$penfam.fit$utx %*% beta_next))
  bi

}
