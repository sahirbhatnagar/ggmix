rm(list=ls())
library(sail)
library(doMC)
registerDoMC(cores = 8)

## ---- toy-example ----
gendataIntro <- function (n, p, corr = 0, E = truncnorm::rtruncnorm(n, a = -1, b = 1),
                          betaE = 2, SNR = 2, hierarchy = c("strong", "weak", "none"),
                          nonlinear = TRUE, interactions = TRUE, causal,
                          not_causal) {
  if (!requireNamespace("truncnorm", quietly = TRUE)) {
    stop("Package \"truncnorm\" needed for this function to simulate data. Please install it.",
         call. = FALSE)
  }

  E <- rbinom(n, 1, prob = 0.5)
  hierarchy <- match.arg(hierarchy)

  W <- replicate(n = p, truncnorm::rtruncnorm(n, a = 0, b = 1))
  U <- truncnorm::rtruncnorm(n, a = 0, b = 1)
  V <- truncnorm::rtruncnorm(n, a = 0, b = 1)

  X1 <- (W[, 1] + corr * U)/(1 + corr)
  X2 <- (W[, 2] + corr * U)/(1 + corr)
  X <- (W[, 3:p] + corr * V)/(1 + corr)
  Xall <- cbind(X1, X2, X)
  colnames(Xall) <- paste0("X", seq_len(p))

  f1 <- function(x) -3 * x
  f2 <- function(x) 2 * (2 * x - 1)^3

  f1.inter <- function(x, e) e * f1(x)
  f2.inter <- function(x, e) e * f2(x)

  error <- stats::rnorm(n)

  Y.star <- f1(X1) + f2(X2) + betaE * E + 1.5 * f2.inter(X2, E)

  k <- sqrt(stats::var(Y.star)/(SNR * stats::var(error)))
  Y <- Y.star + as.vector(k) * error
  return(list(x = Xall, y = Y, e = E, Y.star = Y.star,
              f1 = f1(X1), f2 = f2(X2), betaE = betaE,
              f1.f = f1, f2.f = f2,
              f1.inter = f1.inter, f2.inter = f2.inter,
              X1 = X1, X2 = X2))
}

set.seed(54321)
DT <- gendataIntro(n = 100, p = 20, corr = 0, SNR = 2, betaE = 1.75)

f.basis <- function(i) splines::bs(i, degree = 3)
cvfit <- cv.sail(x = DT$x, y = DT$y, e = DT$e,
                 center.e = FALSE,
                 verbose = 1,
                 nlambda = 100,
                 basis = f.basis, nfolds = 10, parallel = TRUE)

## ---- toy-solution-path ----

plotSailCoefIntro <- function(cvfit, coefs, lambda, group, df, dev, vnames, environ,
                              alpha = 1, legend.loc, label = TRUE, log.l = TRUE,
                              norm = FALSE, ...) {

  if (alpha < 0 | alpha > 1) {
    warning("alpha must be in the range [0,1]. Setting alpha = 1")
    alpha <- 1
  }

  lambda.min.index <- which(cvfit[["lambda.min"]]==cvfit$lambda)
  lambda.1se.index <- which(cvfit[["lambda.1se"]]==cvfit$lambda)
  nzind <- cvfit$sail.fit$group[abs(cvfit$sail.fit$beta[,lambda.1se.index])>0]
  zind <- cvfit$sail.fit$group[abs(cvfit$sail.fit$beta[,lambda.1se.index])==0]

  if (norm) {
    # not implemented for now
  } else {
    if (length(dim(coefs)) == 3) {
      beta <- matrix(coefs[, -1, , drop = FALSE], ncol = dim(coefs)[3])
    } else {
      beta <- coefs
    }
    penalized <- which(group != 0)
    nonzero <- which(apply(abs(beta), 1, sum) != 0)
    ind <- intersect(penalized, nonzero)
    Y <- as.matrix(beta[ind, , drop = FALSE])
    g <- as.numeric(as.factor(group[ind]))
  }
  p <- nrow(Y)
  l <- lambda
  n.g <- max(g)
  if (log.l) {
    l <- log(l)
    index <- l
    approx.f <- 0
    xlab <- expression(log(lambda))
  } else {
    xlab <- expression(lambda)
    index <- lambda
    approx.f <- 0
  }

  ylims <- if (!missing(environ)) range(Y, environ) else range(Y)
  plot.args <- list(
    x = l, y = 1:length(l), ylim = ylims,
    xlab = xlab, ylab = "", type = "n",
    xlim = c(rev(range(l))[1], rev(range(l))[2]),
    cex.lab = 1.05,
    cex.axis = 1.05,
    cex = 1.05,
    family = "serif"
  )
  new.args <- list(...)
  if (length(new.args)) {
    new.plot.args <- new.args[names(new.args) %in% c(
      names(par()),
      names(formals(plot.default))
    )]
    plot.args[names(new.plot.args)] <- new.plot.args
  }
  do.call("plot", plot.args)
  if (plot.args$ylab == "") {
    ylab <- if (norm) {
      expression("||" * hat(theta) * "||")
    } else {
      expression(hat(theta))
    }
    mtext(ylab, 2, 3.5, las = 1, adj = 0, cex = 2*0.5)
  }
  abline(h = 0, lwd = 0.8, col = "gray")
  cols <- hcl(
    h = seq(15, 375, len = max(4, n.g + 1)), l = 60,
    c = 150, alpha = alpha
  )
  cols <- if (n.g == 2) cols[c(1, 3)] else cols[1:n.g]
  line.args <- list(
    col = cols, lwd = 1 + 2 * exp(-p / 20),
    lty = 1, pch = ""
  )
  if (length(new.args)) {
    line.args[names(new.args)] <- new.args
  }

  line.args$x <- l
  line.args$y <- t(Y)

  newcols <- rep("#D3D3D3", n.g)
  newcols[unique(nzind)] <- sail:::cbbPalette[c(4,7)]
  line.args$col <- newcols[g]

  line.args$lty <- rep(line.args$lty, length.out = max(g))
  line.args$lty <- line.args$lty[g]

  line.args$lwd <- rep(line.args$lwd, length.out = max(g))
  line.args$lwd[unique(nzind)] <- 2
  line.args$lwd <- line.args$lwd[g]

  do.call("matlines", line.args)
  # browser()
  line.args$y <- t(Y)[,c("X1_1","X1_2","X1_3","X2_1","X2_2","X2_3")]
  line.args$col <- c(rep(sail:::cbbPalette[c(4)], 3),
                     rep(sail:::cbbPalette[c(7)], 3))
  line.args$lwd <- 2
  do.call("matlines", line.args)

  if (!missing(environ)) lines(l, environ, lwd = 2)

  if (!missing(legend.loc)) {
    legend.args <- list(
      col = cols, lwd = line.args$lwd,
      lty = line.args$lty, legend = vnames
    )
    if (length(new.args)) {
      new.legend.args <- new.args[names(new.args) %in%
                                    names(formals(legend))]
      legend.args[names(new.legend.args)] <- new.legend.args
    }
    legend.args$x <- legend.loc
    do.call("legend", legend.args)
  }
  # browser()
  if (label) {
    ypos <- Y[c("X1_1","X1_2","X1_3","X2_1","X2_2","X2_3"), ncol(Y)]
    ypos["X2_3"] <- 3.1
    ypos["X2_1"] <- 2.5
    ypos["X2_2"] <- -0.15
    ypos["X1_2"] <- -0.68
    ypos["X1_1"] <- -1.2
    # mtext(text = names(ypos), side = 4, at = ypos, las = 1)
    mtext(text = c(TeX("$X1_1$"),TeX("$X1_2$"),TeX("$X1_3$"),
                   TeX("$X2_1$"),TeX("$X2_2$"),TeX("$X2_3$")), side = 4, at = ypos, las = 1, cex = 0.8)
    mtext(text = "E", side = 4, at = environ[length(environ)]-.1, las = 1, cex = 0.8)
  }

  abline(v = log(cvfit$lambda.1se), lty=2)
  mtext(text = latex2exp::TeX("$\\lambda_{1SE}$"), side = 3, at = log(cvfit$lambda.1se)+.1, line = 0, cex=0.9)
  abline(v = log(cvfit$lambda.min), lty=2)
  mtext(text = latex2exp::TeX("$\\lambda_{min}$"), side = 3, at = log(cvfit$lambda.min)-.1, line = 0, cex=0.9)

}

plotSailCoefInterIntro <- function(cvfit, coefs, lambda, group, df, dev, vnames, environ,
                              alpha = 1, legend.loc, label = TRUE, log.l = TRUE,
                              norm = FALSE, ...) {

  # browser()
  if (alpha < 0 | alpha > 1) {
    warning("alpha must be in the range [0,1]. Setting alpha = 1")
    alpha <- 1
  }

  lambda.min.index <- which(cvfit[["lambda.min"]]==cvfit$lambda)
  lambda.1se.index <- which(cvfit[["lambda.1se"]]==cvfit$lambda)
  nzind <- cvfit$sail.fit$group[abs(cvfit$sail.fit$alpha[,lambda.1se.index])>0]
  zind <- cvfit$sail.fit$group[abs(cvfit$sail.fit$alpha[,lambda.1se.index])==0]

  if (norm) { # not implemented for now
  } else {
    if (length(dim(coefs)) == 3) {
      beta <- matrix(coefs[, -1, , drop = FALSE], ncol = dim(coefs)[3])
    } else {
      beta <- coefs
    }
    penalized <- which(group != 0)
    nonzero <- which(apply(abs(beta), 1, sum) != 0)
    ind <- intersect(penalized, nonzero)
    Y <- as.matrix(beta[ind, , drop = FALSE])
    g <- as.numeric(as.factor(group[ind]))
  }
  p <- nrow(Y)
  l <- lambda
  n.g <- max(g)
  if (log.l) {
    l <- log(l)
    index <- l
    approx.f <- 0
    xlab <- expression(log(lambda))
  } else {
    xlab <- expression(lambda)
    index <- lambda
    approx.f <- 0
  }

  ylims <- if (!missing(environ)) range(Y, environ) else range(Y)
  plot.args <- list(
    x = l, y = 1:length(l), ylim = ylims,
    xlab = xlab, ylab = "", type = "n",
    xlim = c(rev(range(l))[1], rev(range(l))[2]),
    cex.lab = 1.05,
    cex.axis = 1.05,
    cex = 1.05,
    family = "serif"
  )
  new.args <- list(...)
  if (length(new.args)) {
    new.plot.args <- new.args[names(new.args) %in% c(
      names(par()),
      names(formals(plot.default))
    )]
    plot.args[names(new.plot.args)] <- new.plot.args
  }
  do.call("plot", plot.args)
  if (plot.args$ylab == "") {
    ylab <- if (norm) {
      expression("||" * hat(theta) * "||")
    } else {
      expression(hat(theta))
    }
    mtext(ylab, 2, 3.5, las = 1, adj = 0, cex = 2*.5)
  }
  abline(h = 0, lwd = 0.8, col = "gray")
  cols <- hcl(
    h = seq(15, 375, len = max(4, n.g + 1)), l = 60,
    c = 150, alpha = alpha
  )
  cols <- if (n.g == 2) cols[c(1, 3)] else cols[1:n.g]
  line.args <- list(
    col = cols, lwd = 1 + 2 * exp(-p / 20),
    lty = 1, pch = ""
  )
  if (length(new.args)) {
    line.args[names(new.args)] <- new.args
  }

  line.args$x <- l
  line.args$y <- t(Y)
# browser()
  newcols <- rep("#D3D3D3", n.g)
  newcols[unique(nzind)] <- sail:::cbbPalette[c(7)]
  line.args$col <- newcols[g]

  line.args$lty <- rep(line.args$lty, length.out = max(g))
  line.args$lty <- line.args$lty[g]

  line.args$lwd <- rep(line.args$lwd, length.out = max(g))
  line.args$lwd[unique(nzind)] <- 2
  line.args$lwd <- line.args$lwd[g]

  do.call("matlines", line.args)
  line.args$y <- t(Y)[,c("X2_1:E","X2_2:E","X2_3:E")]
  line.args$col <- rep(sail:::cbbPalette[c(7)], 3)
  line.args$lwd <- 2
  do.call("matlines", line.args)

  if (!missing(environ)) lines(l, environ, lwd = 2)

  if (!missing(legend.loc)) {
    legend.args <- list(
      col = cols, lwd = line.args$lwd,
      lty = line.args$lty, legend = vnames
    )
    if (length(new.args)) {
      new.legend.args <- new.args[names(new.args) %in%
                                    names(formals(legend))]
      legend.args[names(new.legend.args)] <- new.legend.args
    }
    legend.args$x <- legend.loc
    do.call("legend", legend.args)
  }

  if (label) {
    ypos <- Y[c("X2_1:E","X2_2:E","X2_3:E"), ncol(Y)]
    ypos["X2_1:E"] <- ypos["X2_1:E"]-0.7
    mtext(text = c(TeX("$E\\cdot X2_1$"),TeX("$E\\cdot X2_2$"),TeX("$E\\cdot X2_3$")), side = 4, at = ypos, las = 1, cex = 0.8)
  }

  abline(v = log(cvfit$lambda.1se), lty=2)
  abline(v = log(cvfit$lambda.min), lty=2)

}

# c(bottom, left, top, right)
# dev.off()
trim <- 1:100
x <- cvfit$sail.fit
par(mfrow=c(2,1), tcl=-0.5, family="serif", omi=c(0.2,0.2,0,0))
par(mar=c(0,4,2,3.2))
plotSailCoefIntro(
  cvfit = cvfit,
  coefs = x$beta[,trim],
  environ = x$bE[trim],
  lambda = x$lambda[trim],
  df = (x$dfbeta + x$dfenviron)[trim],
  group = x$group,
  dev = x$dev.ratio[trim],
  vnames = x$vnames,
  ylim = c(-3.5,4),
  ylab = "Main effects", xaxt="n")
par(mar=c(4,4,0,3.2))
plotSailCoefInterIntro(
  cvfit = cvfit,
  coefs = x$alpha[,trim],
  lambda = x$lambda[trim],
  df = x$dfalpha[trim],
  group = x$group,
  dev = x$dev.ratio[trim],
  vnames = x$vnames,
  ylab = "Interactions",
  xlab = TeX("$\\log(\\lambda)$"))


## ---- toy-effects ----

# c(bottom, left, top, right)
# dev.off()
par(mfrow=c(1,2), tcl=-0.5, family="serif", omi=c(0.2,0.2,0,0))
x1 <- DT$X1
# truth
lin_pred <- DT$f1
# estimated
lin_pred_est <- as.vector(cvfit$sail.fit$design[,c("X1_1", "X1_2", "X1_3")] %*% coef(cvfit, s="lambda.min")[c("X1_1", "X1_2", "X1_3"),,drop=F])
min.length.top <- range(c(lin_pred, lin_pred_est))[1] ; max.length.top <- range(c(lin_pred, lin_pred_est))[2]
par(mar=c(4,4,0,0))
plot(x1, lin_pred,
     pch = 19,
     ylab = latex2exp::TeX("$f(X_1)$"),
     xlab = latex2exp::TeX("$X_1$"),
     col = sail:::cbbPalette[c(4,7)],
     bty="n",
     xaxt="n",
     type = "n",
     cex.lab = 1,
     cex.axis = 1,
     cex = 1,
     main = "",
     cex.main = 1,
     ylim = c(min.length.top, max.length.top+1))
axis(1, labels = T, cex.axis = 1)
lines(x1[order(x1)],lin_pred[order(x1)], lwd = 2)
lines(x1[order(x1)],lin_pred_est[order(x1)], lwd = 2, lty = 2)
legend(x=0.2,y=1.3,c("Truth", "Estimated"),
       cex = 1, bty = "n", lty = c(1,2), lwd = 1)
rug(x1)


x2 <- DT$X2
e <- DT$e

# truth
lin_pred <- #DT$betaE * DT$e +
  1.5 * DT$f2.inter(DT$X2, DT$e)

# estimated
lin_pred_est <- #coef(cvfit, s="lambda.min")["E",] * DT$e +
  as.vector(cvfit$sail.fit$design[,c("X2_1:E", "X2_2:E", "X2_3:E")] %*% coef(cvfit, s="lambda.min")[c("X2_1:E", "X2_2:E", "X2_3:E"),,drop=F])

unexposed_index <- which(e==0)
exposed_index <- which(e==1)

e0 <- lin_pred[unexposed_index]
e1 <- lin_pred[exposed_index]
x2e0 <- x2[unexposed_index]
x2e1 <- x2[exposed_index]

e0_est <- lin_pred_est[unexposed_index]
e1_est <- lin_pred_est[exposed_index]

min.length.top <- range(c(lin_pred, lin_pred_est))[1] ; max.length.top <- range(c(lin_pred, lin_pred_est))[2]

par(mar=c(4,4,0,0))
plot(x2, lin_pred,
     pch = 19,
     ylab = latex2exp::TeX("$f(X_2) \\cdot E$"),
     xlab = latex2exp::TeX("$X_2$"),
     col = sail:::cbbPalette[c(4,7)],
     bty="n",
     xaxt="n",
     type = "n",
     cex.lab = 1,
     cex.axis = 1,
     cex = 1,
     main = "",
     cex.main = 1,
     ylim = c(min.length.top, max.length.top))
axis(1, labels = T, cex.axis = 1)
lines(x2e0[order(x2e0)],e0[order(x2e0)], col = sail:::cbbPalette[c(4)], lwd = 2)
lines(x2e1[order(x2e1)],e1[order(x2e1)], col = sail:::cbbPalette[c(7)], lwd = 2)
lines(x2e1[order(x2e1)],e1_est[order(x2e1)], col = sail:::cbbPalette[c(7)], lwd = 2, lty = 2)
legend(x=0.05,y = 2.95, c("E=0", "E=1"),
       col = sail:::cbbPalette[c(4,7)], pch = 19, cex = 1, bty = "n")
rug(x2)


## ---- toy-plots ----

plot(cvfit)
predict(cvfit, s="lambda.min", type = "nonzero")
predict(cvfit, s="lambda.1se", type = "nonzero")
plot(cvfit$sail.fit)


plot(cvfit$sail.fit)#, xlim = c(1,-3), ylim = c(-5,5))
abline(v = log(cvfit[["lambda.1se"]]))
abline(v = log(cvfit[["lambda.min"]]))

plotMain(cvfit$sail.fit, x = DT$x, xvar = "X1",
         s = cvfit$lambda.1se, f.truth = DT$f1.f)
plotMain(cvfit$sail.fit, x = DT$x, xvar = "X2",
         s = cvfit$lambda.1se, f.truth = DT$f2.f)
plotInter(cvfit$sail.fit, x = DT$x, xvar = "X4",
          f.truth = DT$f4.inter,
          s = cvfit$lambda.min,
          title_z = "Estimated")

which(cvfit$lambda.min==cvfit$lambda)
which(cvfit$lambda.1se==cvfit$lambda)
cbind(cvfit$sail.fit$beta[,40,drop=F], cvfit$sail.fit$beta[,30,drop=F])
cbind(cvfit$sail.fit$alpha[,40,drop=F], cvfit$sail.fit$alpha[,30,drop=F])
