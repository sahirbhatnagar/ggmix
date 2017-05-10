"%ni%" <- Negate("%in%")

fr_eta <- function(eta, sigma2, beta, eigenvalues, x, y, nt) {

  # this is based on the negative log-lik

  # eta = 0.5
  # sigma2 = 1
  # beta = beta_next
  # eigenvalues = Lambda
  # x = utx
  # y = uty
  # nt = n
  # ============

  kernel <- 1 + eta * (eigenvalues - 1)

  (nt / 2) * log(2 * pi) +
    (nt / 2) * log(sigma2) +
    0.5 * sum(log(kernel)) +
    (1 / (2 * sigma2)) * sum((y - x %*% beta) ^ 2 / kernel)

}

grr_eta <- function(eta, sigma2, beta, eigenvalues, x, y, nt) {

  kernel <- 1 + eta * (eigenvalues - 1)

  (1 / 2) * sum(((eigenvalues - 1) / kernel) * (((y - x %*% beta) ^ 2) / (sigma2 * kernel) - 1))
}


log_lik <- function(eta, sigma2, beta, eigenvalues, x, y, nt) {

  # this returns the log-likelihood

  # eta = 0.5
  # sigma2 = 1
  # beta = beta_next
  # eigenvalues = Lambda
  # x = utx
  # y = uty
  # nt = n
  # ============

  kernel <- 1 + eta * (eigenvalues - 1)

  -1 * (
    (nt / 2) * log(2 * pi) +
      (nt / 2) * log(sigma2) +
      0.5 * sum(log(kernel)) +
      (1 / (2 * sigma2)) * sum((y - x %*% beta) ^ 2 / kernel)
  )

}





#' Print Method for shim function
#'
#' @description print method for shim function
#' @export

print.penfam <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(Df = x$x[,"Df"],
              `%Dev` = signif(x$x[,"%Dev"], digits),
              Lambda = signif(x$x[,"Lambda"], digits),
              BIC = signif(x$x[,"BIC"], digits)))
}



#' Get coefficients from a "shim" object
#'
#' @rdname predict.shim
#' @export

coef.penfam <- function(object, s = NULL,  ...) {
  predict(object, s = s, type = "coefficients", ...)
}



#' Make predictions from a shim object
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



plotCoef <- function(beta,norm,lambda,df,dev,label=FALSE,
                     xvar=c("norm","lambda","dev"),
                     xlab=iname, ylab="Coefficients", ...){
  # as(x,"CsparseMatrix")
  # beta = as(res$beta,"sparseMatrix")
  # beta = res$beta
  # lambda = drop(res$x[,"Lambda"])
  # df = drop(res$x[,"Df"])
  # dev = drop(res$x[,"%Dev"])
  # xvar = "norm"
  # ===========================
  ##beta should be in "dgCMatrix" format
  ### bystep = FALSE means which variables were ever nonzero
  ### bystep = TRUE means which variables are nonzero for each step
  which=nonzeroCoef(beta, bystep = FALSE)
  nwhich=length(which)
  switch(nwhich+1,#we add one to make switch work
         "0"={
           warning("No plot produced since all coefficients zero")
           return()
         },
         "1"=warning("1 or less nonzero coefficients; glmnet plot is not meaningful")
  )
  beta=as.matrix(beta[which,,drop=FALSE])
  xvar=match.arg(xvar)
  switch(xvar,
         "norm"={
           index=if(missing(norm))apply(abs(beta),2,sum)else norm
           # index=apply(abs(beta),2,sum)
           iname="L1 Norm"
           approx.f=1
         },
         "lambda"={
           index=log(lambda)
           iname="Log Lambda"
           approx.f=0
         },
         "dev"= {
           index=dev
           iname="Fraction Deviance Explained"
           approx.f=1
         }
  )
  dotlist=list(...)
  type=dotlist$type
  if(is.null(type))
    matplot(index,t(beta),lty=1,xlab=xlab,ylab=ylab,type="l",...)
  else matplot(index,t(beta),lty=1,xlab=xlab,ylab=ylab,...)
  atdf=pretty(index)
  ### compute df by interpolating to df at next smaller lambda
  ### thanks to Yunyang Qian
  prettydf=approx(x=index,y=df,xout=atdf,rule=2,method="constant",f=approx.f)$y
  # prettydf=ceiling(approx(x=index,y=df,xout=atdf,rule=2)$y)
  axis(3,at=atdf,labels=prettydf,tcl=NA)
  if(label){
    nnz=length(which)
    xpos=max(index)
    pos=4
    if(xvar=="lambda"){
      xpos=min(index)
      pos=2
    }
    xpos=rep(xpos,nnz)
    ypos=beta[,ncol(beta)]
    text(xpos,ypos,paste(which),cex=.5,pos=pos)
  }

}


plotBIC <- function(object, sign.lambda, lambda.min, ...) {

  # object = res$x
  # sign.lambda = 1
  # lambda.min = res$lambda_min
  # ===============

  xlab="log(Lambda)"
  lambda_min <- drop(object[lambda.min,"Lambda"])
  if(sign.lambda<0) xlab=paste("-",xlab,sep="")
  plot.args=list(x=sign.lambda*log(drop(object[,"Lambda"])),
                 y=drop(object[,"BIC"]),
                 ylim=range(drop(object[,"BIC"])),
                 xlab=xlab,
                 ylab="BIC", type="n")
  new.args=list(...)
  if (length(new.args)) plot.args[names(new.args)]=new.args
  do.call("plot",plot.args)
  points(sign.lambda*log(drop(object[,"Lambda"])),
         drop(object[,"BIC"]),pch=20,col="red")
  axis(side=3,at=sign.lambda*log(drop(object[,"Lambda"])),
       labels=paste(drop(object[,"Df"])), tick=FALSE, line=0)
  abline(v=sign.lambda*log(lambda_min),lty=3)
  # abline(v=sign.lambda*log(.1605),lty=3)
  # abline(v=sign.lambda*log(cvobj$lambda.1se),lty=3)
  # invisible()
}


#' Calculate Sequence of Tuning Parameters
#'
#' @description Function to calculate the sequence of tuning parameters based on
#'   the design matrix \code{x} and the response variable {y}. This is used in
#'   the \code{\link{shim_once}} function to calculate the tuning parameters
#'   applied to the main effects
#'
#' @inheritParams uni_fun
#' @param weights Separate penalty factors can be applied to each coefficient.
#'   This is a number that multiplies lambda to allow differential shrinkage,
#'   and can be used to apply adaptive LASSO. Can be 0 for some variables, which
#'   implies no shrinkage, and that variable is always included in the model.
#'   Default is 1 for all variables (and implicitly infinity for variables
#'   listed in exclude). Note: the penalty factors are internally rescaled to
#'   sum to nvars, and the lambda sequence will reflect this change.
#' @param lambda.factor The factor for getting the minimal lambda in lambda
#'   sequence, where \code{min(lambda) = lambda.factor * max(lambda).
#'   max(lambda)} is the smallest value of lambda for which all coefficients are
#'   zero. The default depends on the relationship between \code{N} (the number
#'   of rows in the matrix of predictors) and \code{p} (the number of
#'   predictors). If \code{N > p}, the default is \code{1e-6}, close to zero. If
#'   \code{N < p}, the default is \code{0.01}. A very small value of
#'   lambda.factor will lead to a saturated fit.
#' @param nlambda the number of lambda values - default is 100.
#' @param scale_x should the columns of x be scaled - default is FALSE
#' @param center_y should y be mean centered - default is FALSE
#' @return numeric vector of length \code{q}
#' @details The maximum lambda is calculated using the following inequality:
#'   \deqn{(N*w_j)^-1 | \sum x_ij y_i | \le \lambda_max}
#'
#'   The minimum lambda is given by lambda.factor*lambda_max. The sequence of
#'   nlambda values are decreasing from lambda_max to lambda_min on the log
#'   scale.
#'
#'   The penalty factors are internally rescaled to sum to the number of
#'   predictor variables in glmnet. Therefore, to get the correct sequence of
#'   lambdas when there are weights, this function first rescales the weights
#'   and then calclated the sequence of lambdas.
#'
#'   This formula is taken from section 2.5 of the \code{glmnet} paper in the
#'   Journal of Statistical Software (see references for details)
#'
#' @author
#' Sahir Bhatnagar
#'
#' Maintainer: Sahir Bhatnagar \email{sahir.bhatnagar@@mail.mcgill.ca}
#'
#' @references Friedman, J., Hastie, T. and Tibshirani, R. (2008)
#'   \emph{Regularization Paths for Generalized Linear Models via Coordinate
#'   Descent}, \url{http://www.stanford.edu/~hastie/Papers/glmnet.pdf}
#'   \emph{Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010}
#'   \url{http://www.jstatsoft.org/v33/i01/}
#'
#'   Yang, Y., & Zou, H. (2015). A fast unified algorithm for solving
#'   group-lasso penalize learning problems. \emph{Statistics and Computing},
#'   25(6), 1129-1141.
#'   \url{http://www.math.mcgill.ca/yyang/resources/papers/gglasso.pdf}
#'
#'
#' @examples
#' # number of observations
#' n <- 100
#'
#' # number of predictors
#' p <- 5
#'
#' # environment variable
#' e <- sample(c(0,1), n, replace = T)
#'
#' # main effects
#' x <- cbind(matrix(rnorm(n*p), ncol = p), e)
#'
#' # need to label columns
#' dimnames(x)[[2]] <- c(paste0("x",1:p), "e")
#'
#' # design matrix without intercept
#' X <- model.matrix(~(x1+x2+x3+x4+x5)*e-1, data = as.data.frame(x))
#'
#' # response
#' Y <- X %*% rbinom(ncol(X), 1, 0.2) + 3*rnorm(n)
#'
#' lambda_sequence(X,Y)
#' @export

# lambda_sequence <- function(x, y, weights = NULL,
#                             lambda.factor = ifelse(nobs < nvars, 0.01, 1e-06),
#                             nlambda = 100, scale_x = F, center_y = F) {
#
#   # when scaling, first you center then you standardize
#   if (any(as.vector(weights) < 0)) stop("Weights must be positive")
#   np <- dim(x)
#   nobs <- as.integer(np[1])
#   nvars <- as.integer(np[2])
#
#   if (!is.null(weights) & length(as.vector(weights)) < nvars)
#     stop("You must provide weights for every column of x")
#
#   # scale the weights to sum to nvars
#   w <- if (is.null(weights)) rep(1, nvars) else as.vector(weights) / sum(as.vector(weights)) * nvars
#
#   sx <- if (scale_x) apply(x,2, function(i) scale(i, center = TRUE, scale = mysd(i))) else x
#   sy <- if (center_y) as.vector(scale(y, center = T, scale = F)) else as.vector(y)
#   lambda.max <- max(abs(colSums(sy * sx) / w)) / nrow(sx)
#
#   rev(exp(seq(log(lambda.factor * lambda.max), log(lambda.max), length.out = nlambda)))
# }





