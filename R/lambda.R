#' Estimation of Lambda Sequence for Linear Mixed Model with Lasso Penalty
#'
#' @description \code{lambdalasso} estimates a decreasing sequence of tuning
#'   parameters
#'
#' @seealso \code{\link{ggmix}}
#' @param ggmix_object A ggmix_object object of class \code{lowrank} or
#'   \code{fullrank}
#' @inheritParams ggmix
#' @param scale_x should the columns of x be scaled - default is FALSE
#' @param center_y should y be mean centered - default is FALSE.
#' @param tol.kkt KKT tolerance. Currently ignored
#' @param ... Extra parameters. Currently ignored.
#' @note This function isn't meant to be called directly by the user.
#' @return A decreasing sequence of tuning parameters
lambdalasso <- function(ggmix_object, ...) UseMethod("lambdalasso")

#' @rdname lambdalasso
lambdalasso.default <- function(ggmix_object, ...) {
  stop(strwrap("This function should be used with a ggmix object of class
               lowrank or fullrank"))
}

#' @rdname lambdalasso
lambdalasso.fullrank <- function(ggmix_object,
                                 ...,
                                 penalty.factor,
                                 lambda_min_ratio,
                                 epsilon = 1e-14,
                                 tol.kkt = 1e-9,
                                 eta_init = 0.5,
                                 nlambda = 100, scale_x = F, center_y = F) {
  utx <- ggmix_object[["x"]]
  uty <- ggmix_object[["y"]]
  eigenvalues <- ggmix_object[["D"]]

  np <- dim(utx)
  n <- np[[1]]

  # assuming the first column is the intercept, so we subtract 1
  p <- np[[2]] - 1

  # weights
  di <- 1 + eta_init * (eigenvalues - 1)

  # initial value for beta0
  beta0_init <- (sum(utx[, 1] * uty / di)) / (sum(utx[, 1]^2 / di))

  # this includes all other betas which are 0 by definition
  beta_init <- as.matrix(c(beta0_init, rep(0, p)))

  # sum is faster
  # microbenchmark::microbenchmark(
  #   mat = drop((t(utx0) %*% di_inverse %*% uty) / (t(utx0) %*% di_inverse %*% utx0)),
  #     sum = (sum(utx0 * uty / di)) / (sum(utx0 ^ 2 / di)), times = 1000
  # )

  # closed form for sigma^2
  sigma2_init <- (1 / n) * sum((uty - utx %*% beta_init)^2 / di)

  # sum version is faster
  # mb <- microbenchmark(
  #   mat = (1 / n) * t(uty - beta0_init * utx0) %*% di_inverse %*% (uty - beta0_init * utx0),
  #   sum = (1 / n) * sum ((uty - beta0_init * utx0)^2 / (1 + eta_init * (eigenvalues - 1))),
  #   times = 1000)
  # ggplot2::autoplot(mb)

  # iteration counter
  k <- 0

  # to enter while loop
  converged <- FALSE

  while (!converged) {
    Theta_init <- c(beta_init, eta_init, sigma2_init)

    # fit eta
    eta_next <- stats::optim(
      par = eta_init,
      fn = fn_eta_lasso_fullrank,
      gr = gr_eta_lasso_fullrank,
      method = "L-BFGS-B",
      control = list(fnscale = 1),
      lower = .01,
      upper = .99,
      sigma2 = sigma2_init,
      beta = beta_init,
      eigenvalues = eigenvalues,
      x = utx,
      y = uty,
      nt = n
    )$par

    # weights
    di <- 1 + eta_next * (eigenvalues - 1)

    # next value for beta0
    beta0_next <- (sum(utx[, 1] * uty / di)) / (sum(utx[, 1]^2 / di))

    beta_next <- as.matrix(c(beta0_next, rep(0, p)))

    # closed form for sigma^2
    sigma2_next <- (1 / n) * sum((uty - utx %*% beta_next)^2 / di)

    k <- k + 1

    Theta_next <- c(beta_next, eta_next, sigma2_next)

    converged <- crossprod(Theta_next - Theta_init) < epsilon

    # message(sprintf("l2 norm squared of Theta_k+1 - Theta_k: %f \n log-lik: %f",
    #                 crossprod(Theta_next - Theta_init),
    #                 log_lik(eta = eta_next, sigma2 = sigma2_next, beta = beta_next,
    #                         eigenvalues = eigenvalues,
    #                         x = utx, y = uty, nt = n)))

    beta0_init <- beta0_next
    beta_init <- beta_next
    eta_init <- eta_next
    sigma2_init <- sigma2_next
  }

  # eta_next
  # sigma2_next
  di <- 1 + eta_next * (eigenvalues - 1)
  wi <- (1 / sigma2_next) * (1 / di)
  if (any(wi < 0)) stop("weights are negative")

  # scale the weights to sum to nvars
  # wi_scaled <- as.vector(wi) / sum(as.vector(wi)) * n

  # wi_scaled <- as.vector(wi) * n

  # lambda.max <- max(abs(colSums((wi * utx[,-1]) * drop(uty - utx %*% beta_next))))

  # this gives the same answer (see paper for details)
  # we divide by sum(wi) here and not in glmnet because the sequence is determined
  # on the log scale
  # lambda.max <- max(abs(colSums(((1 / sum(wi_scaled)) * (wi_scaled * utx[,-1]) * drop(uty - utx %*% beta_next)))))
  # browser()
  lambdas <- (1 / penalty.factor) *
    abs(colSums(((1 / sum(wi)) * (wi * utx[, -1]) *
      drop(uty - utx %*% beta_next))))
  # need to check for Inf, in case some penalty factors are 0
  lambda.max <- max(lambdas[lambdas != Inf])
  # lambda.max <- lambda.max * sum(wi)
  # (utx[,-1, drop = F]) %>% dim
  # a <- colSums(utx[,-1, drop = F]^2 * wi)
  # b <- colSums(sweep(utx[,-1, drop = F]^2, MARGIN = 1, wi, '*'))
  # all(a == b)

  # kkt <- kkt_check(eta = eta_next, sigma2 = sigma2_next, beta = beta_next,
  #                  eigenvalues = eigenvalues, x = utx, y = uty, nt = n,
  #                  lambda = lambda.max, tol.kkt = tol.kkt)
  # message(kkt)
  out <- list(
    sequence = rev(exp(seq(log(lambda_min_ratio * lambda.max),
      log(lambda.max),
      length.out = nlambda
    ))),
    eta = eta_next, sigma2 = sigma2_next, beta0 = beta0_next
  )

  return(out)
}


# still need to create lambdalasso.lowrank, lambdagglasso.fullrank, lambdagglasso.lowrank
