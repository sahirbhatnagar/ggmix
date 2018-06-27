#' Plot the Generalised Information Criteria curve produced by \code{gic}
#'
#' @description Plots the Generalised Information Criteria curve, as a function
#'   of the lambda values used
#' @param x fitted linear mixed model object of class \code{ggmix_gic} from the
#'   \code{\link{gic}} function
#' @param sign.lambda Either plot against log(lambda) (default) or its negative
#'   if sign.lambda=-1
#' @param lambda.min the value of lambda which minimizes the gic
#' @param ... Other graphical parameters to plot
#' @details A plot is produced, and nothing is returned.
#' @seealso \code{\link{ggmix}} and \code{\link{gic}}
#' @export
plot.ggmix_gic <- function(x, sign.lambda = 1, ...) {
  plotGIC(
    x = x,
    sign.lambda = sign.lambda,
    lambda.min = x$lambda.min, ...
  )
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





plot.ggmix_fit <- function(x,
                           type = c("coef", "QQranef", "QQresid", "predicted", "Tukey-Anscombe"),
                           xvar = c("norm", "lambda", "dev"), s = x$lambda_min,
                           label = FALSE, sign.lambda = 1, ...) {
  xvar <- match.arg(xvar)
  type <- match.arg(type, several.ok = TRUE)

  if (any(type == "coef")) {
    plotCoef(x$beta,
             lambda = drop(x$result[, "Lambda"]),
             df = drop(x$result[, "Df"]), dev = drop(x$result[, "Deviance"]),
             label = label, xvar = xvar, ...
    )
  }

  if (any(type == "QQranef")) {
    if (s %ni% rownames(x$result)) stop("value for s not in lambda sequence")
    stats::qqnorm(x$randomeff[, s], main = sprintf("QQ-Plot of the random effects at lambda = %.2f", x$result[s, "Lambda"]))
    stats::qqline(x$randomeff[, s], col = "red")
  }

  if (any(type == "QQresid")) {
    if (s %ni% rownames(x$result)) stop("value for s not in lambda sequence")
    stats::qqnorm(x$residuals[, s], main = sprintf("QQ-Plot of the residuals at lambda = %.2f", x$result[s, "Lambda"]))
    stats::qqline(x$residuals[, s], col = "red")
  }

  if (any(type == "predicted")) {
    if (s %ni% rownames(x$result)) stop("value for s not in lambda sequence")
    graphics::plot(x$predicted[, s], drop(x$y),
         xlab = "predicted response (XB + b)", ylab = "observed response",
         main = sprintf("Observed vs. Predicted response
                        R2 = %g", stats::cor(x$predicted[, s], drop(x$y)))
    )
    graphics::abline(a = 0, b = 1, col = "red")
  }


  if (any(type == "Tukey-Anscombe")) {
    graphics::plot(x$fitted[, s], x$residuals[, s],
         main = "Tukey-Anscombe Plot",
         xlab = "fitted values (XB)", ylab = "residuals"
    )
    graphics::abline(h = 0, col = "red")
  }
}



plotCoef <- function(beta, norm, lambda, df, dev, label = FALSE,
                     xvar = c("norm", "lambda", "dev"),
                     xlab = iname, ylab = "Coefficients", ...) {
  # as(x,"CsparseMatrix")
  # beta = as(res$beta,"sparseMatrix")
  # beta = res$beta
  # lambda = drop(res$x[,"Lambda"])
  # df = drop(res$x[,"Df"])
  # dev = drop(res$x[,"%Dev"])
  # xvar = "norm"
  # ===========================
  ## beta should be in "dgCMatrix" format
  ### bystep = FALSE means which variables were ever nonzero
  ### bystep = TRUE means which variables are nonzero for each step
  which <- glmnet::nonzeroCoef(beta, bystep = FALSE)
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
