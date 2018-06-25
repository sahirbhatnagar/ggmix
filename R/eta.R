#' Gradient of eta parameter
#'
#' @description Used for gradient of eta. Currently being passed to optim and
#'   used in \code{kkt_check}
#' @param x should be U^T X, where U is the matrix of eigenvectors and X
#'   contains the first column of ones for the intercept.
#' @inheritParams logliklasso
#' @rdname kkt_check
grr_eta <- function(eta, sigma2, beta, eigenvalues, x, y, nt) {

  di <- 1 + eta * (eigenvalues - 1)

  (1 / 2) * sum(((eigenvalues - 1) / di) * (1 - (((y - x %*% beta) ^ 2) / (sigma2 * di))))
}



#' Function used to optimize eta
#' @rdname kkt_check
fr_eta <- function(eta, sigma2, beta, eigenvalues, x, y, nt) {

  # this is based on the negative log-lik

  di <- 1 + eta * (eigenvalues - 1)

  (nt / 2) * log(2 * pi) +
    (nt / 2) * log(sigma2) +
    0.5 * sum(log(di)) +
    (1 / (2 * sigma2)) * sum((y - x %*% beta) ^ 2 / di)

}
