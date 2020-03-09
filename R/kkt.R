#' Check of KKT Conditions for Linear Mixed Model
#'
#' @description This function checks the KKT conditions
#' @param x rotated x. Should be U^T X, where U is the matrix of eigenvectors
#'   and X contains the first column of ones for the intercept. x should be a
#'   mtrix of dimension n x (p+1). These are outputted by the constructor
#'   functions. See \code{\link{ggmix_data_object}} for details
#' @param y rotated y. Should be U^T Y, where U is the matrix of eigenvectors
#'   and Y is the response.
#' @param eigenvalues non-zero eigenvalues of the kinship matrix, or the square
#'   of the singular values of the matrix used to construct the kinship matrix
#' @inheritParams ggmix
#' @inheritParams logliklasso
#' @param tol.kkt Tolerance for determining if an entry of the subgradient is
#'   zero
#' @rdname kkt_check
#' @note \code{grr_sigma2} and \code{grr_beta0} are functions for the gradient
#'   of sigma2 and beta0, respectively
kkt_check <- function(eta, sigma2, beta, eigenvalues, x, y, nt,
                                    lambda, tol.kkt = 1e-3) {

  # eta = eta_next; sigma2 = sigma2_next;
  # beta = beta_next;
  # # beta = as.matrix(c(beta0_next,rep(0,p)))
  # eigenvalues = Lambda;
  # x = utx;
  # y = uty; nt = n;
  # lambda = lambda; tol.kkt = 0.000001
  # =======================

  di <- 1 + eta * (eigenvalues - 1)
  wi <- (1 / sigma2) * (1 / di)


  # scale the weights to sum to nvars
  wi_scaled <- as.vector(wi) / sum(as.vector(wi)) * nt
  # wi_mat <- diag(wi_scaled)

  # KKT for beta0
  # grr_beta0(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt)
  # kkt_beta0 <- abs(grr_beta0(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt)) < tol.kkt

  kkt_beta0 <- grr_beta0(
    eta = eta, sigma2 = sigma2, beta = beta,
    eigenvalues = eigenvalues, x = x, y = y, nt = nt
  )

  # KKT for eta
  # gr_eta_lasso_fullrank(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt)
  # kkt_eta <- gr_eta_lasso_fullrank(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt) < tol.kkt
  kkt_eta <- gr_eta_lasso_fullrank(
    eta = eta, sigma2 = sigma2, beta = beta,
    eigenvalues = eigenvalues, x = x, y = y, nt = nt
  )

  # KKT for sigma2
  # grr_sigma2(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt)
  # kkt_sigma2 <- grr_sigma2(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt) < tol.kkt
  kkt_sigma2 <- grr_sigma2(
    eta = eta, sigma2 = sigma2, beta = beta,
    eigenvalues = eigenvalues, x = x, y = y, nt = nt
  )

  # KKT for beta
  # g0 <- (1 / nt) * crossprod(x[,-1, drop = F], wi) %*% (y - x %*% beta) / (colSums(sweep(x[,-1, drop = F]^2, MARGIN = 1, wi_vec, '*')))

  # KKT for beta
  # g0 <- (1 / sum(wi_scaled)) * crossprod(x[,-1, drop = F] * wi_scaled, (y - x %*% beta ))
  g0 <- (1 / sum(wi)) * crossprod(x[, -1, drop = F] * wi, (y - x %*% beta))

  # this gives same result as g0
  # g1 <- colSums((1 / nt) * sweep(sweep(x[,-1], MARGIN = 1, wi_vec, '*'), MARGIN = 1, drop((y - x %*% beta)),'*'))
  # plot(drop(g0),g1)
  # abline(a=0,b=1)

  g <- g0 - lambda * sign(beta[-1])

  # this is for when beta=0 and should be between -1 and 1
  gg <- g0 / lambda

  # which of the betas are non-zero
  oo <- abs(beta[-1]) > 0

  # if all betas are 0 then set to TRUE, else abs(g[oo]) will give error since 'oo' is all FALSE
  # kkt_beta_nonzero <- if (all(!oo)) TRUE else max(abs(g[oo])) < tol.kkt
  # kkt_beta_subgr1 <- min(gg[!oo]) > -1
  # kkt_beta_subgr2 <- max(gg[!oo]) < 1

  kkt_beta_nonzero <- if (all(!oo)) 0 else sum(abs(g[oo]) > tol.kkt)
  kkt_beta_subgr <- sum(abs(gg[!oo]) > 1)
  # if (sum(abs(g[oo]) > tol.kkt) > 0) plot(abs(g[oo]))

  return(c(
    kkt_beta0 = kkt_beta0,
    kkt_eta = kkt_eta,
    kkt_sigma2 = kkt_sigma2,
    kkt_beta_nonzero = kkt_beta_nonzero,
    kkt_beta_subgr = kkt_beta_subgr
  ))
}

#' @rdname kkt_check
grr_sigma2 <- function(eta, sigma2, beta, eigenvalues, x, y, nt) {
  di <- 1 + eta * (eigenvalues - 1)

  sigma2 - (1 / nt) * sum((((y - x %*% beta)^2) / di))
}

#' @rdname kkt_check
grr_beta0 <- function(eta, sigma2, beta, eigenvalues, x, y, nt) {
  di <- 1 + eta * (eigenvalues - 1)
  wi <- (1 / sigma2) * diag(1 / di)

  as.numeric(crossprod(x[, 1, drop = FALSE], wi) %*% (y - x %*% beta))
}
