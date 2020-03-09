#' Plot the Generalised Information Criteria curve produced by \code{gic}
#'
#' @description Plots the Generalised Information Criteria curve, as a function
#'   of the lambda values used
#' @param x fitted linear mixed model object of class \code{ggmix_gic} from the
#'   \code{\link{gic}} function
#' @param sign.lambda Either plot against log(lambda) (default) or its negative
#'   if sign.lambda=-1
#' @param lambda.min the value of lambda which minimizes the gic
#' @param type \code{gic} returns a plot of the GIC vs. log(lambda).
#'   \code{QQranef} return a qqplot of the random effects. \code{QQresid}
#'   returns a qqplot of the residuals which is \eqn{y - X\beta - b_i} where b_i
#'   is the subject specific random effect. \code{predicted} returns a plot of
#'   the predicted response (\eqn{X \beta} + b_i) vs. the observed response,
#'   where b_i is the subject specific random effect. \code{Tukey-Anscombe}
#'   returns a plot of the residuals vs. fitted values (\eqn{X \beta})
#' @param s Value of the penalty parameter \code{lambda} at which predictions
#'   are required. Default is the value \code{s="lambda.min"}. If \code{s} is
#'   numeric, it is taken as the value of \code{lambda} to be used. Must be a
#'   single value of the penalty parameter \code{lambda} at which coefficients
#'   will be extracted via the \code{coef} method for objects of class
#'   \code{ggmix_gic}. If more than one is supplied, only the first one will be
#'   used.
#' @param newy the response variable that was provided to \code{ggmix}. this is
#'   only required for \code{type="QQresis"}, \code{type="Tukey-Anscombe"} and
#'   \code{type="predicted"}
#' @param newx matrix of values for \code{x} at which predictions are to be
#'   made. Do not include the intercept. this is only required for
#'   \code{type="QQresis"}, \code{type="Tukey-Anscombe"} and
#'   \code{type="predicted"}
#' @param ... Other graphical parameters to plot
#' @return plot depends on the type selected
#' @details A plot is produced, and nothing is returned.
#' @seealso \code{\link{gic}}
#' @examples
#' data("admixed")
#' fit <- ggmix(x = admixed$xtrain,
#'              y = admixed$ytrain,
#'              kinship = admixed$kin_train)
#' hdbic <- gic(fit)
#'
#' # plot solution path
#' plot(fit)
#'
#' # plot HDBIC curve as a function of lambda
#' plot(hdbic)
#' @export
plot.ggmix_gic <- function(x, ..., sign.lambda = 1,
                           type = c("gic", "QQranef", "QQresid", "predicted", "Tukey-Anscombe"),
                           s = "lambda.min", newy, newx) {

  type <- match.arg(type, several.ok = FALSE)

  if (length(s) > 1) {
    s <- s[[1]]
    warning("More than 1 s value provided. Only first element will be used for the estimated coefficients.")
  }

  if (is.numeric(s)) {
    lambda <- s
  } else
    if (is.character(s)) {
      s <- match.arg(s)
      lambda <- x[[s]]
    }
  else {
    stop("Invalid form for s")
  }

  if (type == "gic") {
  plotGIC(
    x = x,
    sign.lambda = sign.lambda,
    lambda.min = lambda, ...
  )
  }

  if (type == "QQranef") {
    stats::qqnorm(ranef(x, s = lambda), main = sprintf("QQ-Plot of the random effects at lambda = %.2f", lambda))
    stats::qqline(ranef(x, s = lambda), col = "red")
  }

  if (type == "QQresid") {
    if (missing(newy) | missing(newx))
      stop("newy and newx must be provided when type='QQresid'")

    resids <- newy -
      stats::predict(x, s = lambda, newx = newx) -
      ranef(x, s = lambda)

    stats::qqnorm(resids, main = sprintf("QQ-Plot of the residuals at lambda = %.2f", lambda))
    stats::qqline(resids, col = "red")
  }

  if (type == "predicted") {
    if (missing(newy) | missing(newx))
      stop("newy and newx must be provided when type='QQresid'")

    preds <- stats::predict(x, s = lambda, newx = newx) +
      ranef(x, s = lambda)

    graphics::plot(preds, drop(newy),
                   xlab = "predicted response (XB + b_i)", ylab = "observed response",
                   main = strwrap(sprintf("Observed vs. Predicted response\n
                                          corr(observed,predicted)^2 = %g", stats::cor(preds, drop(newy))^2))
    )
    graphics::abline(a = 0, b = 1, col = "red")
  }

  if (type == "Tukey-Anscombe") {
    if (missing(newy) | missing(newx))
      stop("newy and newx must be provided when type='QQresid'")

    resids <- newy -
      stats::predict(x, s = lambda, newx = newx) -
      ranef(x, s = lambda)
    fitted <- stats::predict(x, s = lambda, newx = newx)
    graphics::plot(fitted, resids,
                   main = "Tukey-Anscombe Plot",
                   xlab = "fitted values (XB)", ylab = "residuals"
    )
    graphics::abline(h = 0, col = "red")
  }
}

#' @rdname plot.ggmix_gic
plotGIC <- function(x, sign.lambda, lambda.min, ...) {
  object <- x
  xlab <- "log(Lambda)"

  if (sign.lambda < 0) xlab <- paste("-", xlab, sep = "")

  plot.args <- list(
    x = sign.lambda * log(drop(object[["lambda"]])),
    y = drop(object[["gic"]]),
    ylim = range(drop(object[["gic"]])),
    xlab = xlab,
    ylab = "GIC", type = "n"
  )

  new.args <- list(...)

  if (length(new.args)) plot.args[names(new.args)] <- new.args

  do.call("plot", plot.args)

  graphics::points(sign.lambda * log(drop(object[["lambda"]])),
    drop(object[["gic"]]),
    pch = 20, col = "red"
  )

  graphics::axis(
    side = 3, at = sign.lambda * log(drop(object[["lambda"]])),
    labels = paste(drop(object[["nzero"]])), tick = FALSE, line = 0
  )

  graphics::abline(v = sign.lambda * log(lambda.min), lty = 3)
}







#' @title Plot Method for \code{ggmix_fit} object
#' @description Produces a coefficient profile plot of the coefficient paths for
#'   a fitted \code{ggmix_fit} object.
#' @param x a \code{ggmix_fit} object
#' @param xvar What is on the X-axis. "norm" plots against the L1-norm of the
#'   coefficients, "lambda" against the log-lambda sequence, and "dev" against
#'   the percent deviance explained.
#' @param label If TRUE, label the curves with variable sequence numbers.
#' @param sign.lambda Either plot against log(lambda) (default) or its negative
#'   if sign.lambda=-1
#' @param ... other graphical parameters passed to \code{plot}
#' @param beta fixed effects estimates
#' @param norm l1 norm of fixed effect estimates. if missing, (default) this
#'   function will calculate it
#' @param lambda sequence of tuning parameters
#' @param df number of non-zero fixed + random effects
#' @param dev percent deviance
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @details A coefficient profile plot is produced
#' @return A plot is produced and nothing is returned
#' @export
plot.ggmix_fit <- function(x,...,
                           xvar = c("norm", "lambda", "dev"),
                           label = FALSE, sign.lambda = 1) {
  xvar <- match.arg(xvar)

  plotCoef(x[["beta"]],
           lambda = drop(x[["result"]][, "Lambda"]),
           df = drop(x[["result"]][,"Df"]), dev = drop(x[["result"]][,"%Dev"]),
           label = label, xvar = xvar, ...)
}



#' @rdname plot.ggmix_fit
plotCoef <- function(beta, norm, lambda, df, dev, label = FALSE,
                     xvar = c("norm", "lambda", "dev"),
                     xlab = iname, ylab = "Coefficients", ...) {

  ## beta should be in "dgCMatrix" format
  ### bystep = FALSE means which variables were ever nonzero
  ### bystep = TRUE means which variables are nonzero for each step
  which <- nonzeroCoef(beta, bystep = FALSE)
  nwhich <- length(which)
  switch(nwhich + 1, # we add one to make switch work
    "0" = {
      warning("No plot produced since all coefficients zero")
      return()
    },
    "1" = warning("1 or less nonzero coefficients; glmnet plot is not meaningful")
  )
  beta <- as.matrix(beta[which, , drop = FALSE])
  xvar <- match.arg(xvar)
  switch(xvar,
    "norm" = {
      index <- if (missing(norm)) apply(abs(beta), 2, sum) else norm
      # index=apply(abs(beta),2,sum)
      iname <- "L1 Norm"
      approx.f <- 1
    },
    "lambda" = {
      index <- log(lambda)
      iname <- "Log Lambda"
      approx.f <- 0
    },
    "dev" = {
      index <- dev
      iname <- "Fraction Deviance Explained"
      approx.f <- 1
    }
  )
  dotlist <- list(...)
  type <- dotlist$type
  if (is.null(type)) {
    graphics::matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, type = "l", ...)
  } else {
    graphics::matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, ...)
  }
  atdf <- pretty(index)
  ### compute df by interpolating to df at next smaller lambda
  ### thanks to Yunyang Qian
  prettydf <- stats::approx(x = index, y = df, xout = atdf, rule = 2, method = "constant", f = approx.f)$y
  # prettydf=ceiling(approx(x=index,y=df,xout=atdf,rule=2)$y)
  graphics::axis(3, at = atdf, labels = prettydf, tcl = NA)
  if (label) {
    nnz <- length(which)
    xpos <- max(index)
    pos <- 4
    if (xvar == "lambda") {
      xpos <- min(index)
      pos <- 2
    }
    xpos <- rep(xpos, nnz)
    ypos <- beta[, ncol(beta)]
    graphics::text(xpos, ypos, paste(which), cex = .5, pos = pos)
  }
}



