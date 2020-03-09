#' Constructor functions for the different ggmix objects
#'
#' \code{new_fullrank_kinship}, \code{new_fullrank_K}, \code{new_fullrank_UD},
#' \code{new_lowrank_kinship}, \code{new_lowrank_K} and \code{new_lowrank_UD}
#' create the ggmix objects from the provided data that are necessary to fit the
#' penalized linear mixed model according to the user's parameters.
#'
#' @seealso \code{\link{ggmix}}
#' @inheritParams ggmix
#' @param n_zero_eigenvalues the number of desired or specified zero eigenvalues.
#'   This is only needed when \code{estimation="lowrank"}, and is calculated
#'   internally by the \code{\link{ggmix}} function. It is equal to the number
#'   of observations minus \code{n_nonzero_eigenvalues}
#'
#' @return A ggmix object, of the class that corresponds to the estimation
#'   method. These objects are lists that contain the data necessary for
#'   computation. These functions are not meant to be called directly by the
#'   user
#' @name ggmix_data_object
NULL

#' @rdname ggmix_data_object
new_fullrank_kinship <- function(x, y, kinship) {
  phi_eigen <- eigen(kinship)
  U_kinship <- phi_eigen$vectors
  Lambda <- phi_eigen$values
  if (any(Lambda < 1e-5)) {
    Lambda[which(Lambda < 1e-5)] <- 1e-05
  }

  x <- cbind("(Intercept)" = 1, x)
  utx <- crossprod(U_kinship, x)
  uty <- crossprod(U_kinship, y)

  structure(list(
    x = utx,
    y = uty,
    U = U_kinship,
    D = Lambda
  ),
  class = c("fullrank")
  )
}

#' @rdname ggmix_data_object
new_fullrank_K <- function(x, y, K) {
  svdX <- svd(K)
  U_K <- svdX$u
  Lambda <- svdX$d^2
  if (any(Lambda < 1e-5)) {
    Lambda[which(Lambda < 1e-5)] <- 1e-05
  }

  x <- cbind("(Intercept)" = 1, x)
  utx <- crossprod(U_K, x)
  uty <- crossprod(U_K, y)

  structure(list(
    x = utx,
    y = uty,
    U = U_K,
    D = Lambda
  ),
  class = c("fullrank")
  )
}

#' @rdname ggmix_data_object
new_fullrank_UD <- function(x, y, U, D) {
  x <- cbind("(Intercept)" = 1, x)
  utx <- crossprod(U, x)
  uty <- crossprod(U, y)

  structure(list(
    x = utx,
    y = uty,
    U = U,
    D = D
  ),
  class = c("fullrank")
  )
}


#' @rdname ggmix_data_object
new_lowrank_kinship <- function(x, y, kinship,
                                n_nonzero_eigenvalues,
                                n_zero_eigenvalues) {
  phi_eigen <- RSpectra::eigs(kinship, k = n_nonzero_eigenvalues)
  U_kinship <- phi_eigen$vectors
  Lambda <- phi_eigen$values
  if (any(Lambda < 1e-5)) {
    Lambda[which(Lambda < 1e-5)] <- 1e-05
  }

  # for lowrank gglasso, and lasso, we use the original X and Y
  # because of the W matrix
  # we should calculate the W matrix in the lasso and gglaso methods
  # since it depends on sigma and eta

  structure(list(
    x = x, y = y,
    U = U_kinship,
    D = Lambda,
    n_nonzero_eigenvalues = n_nonzero_eigenvalues,
    n_zero_eigenvalues = n_zero_eigenvalues
  ),
  class = c("lowrank")
  )
}

#' @rdname ggmix_data_object
new_lowrank_K <- function(x, y, K,
                          n_nonzero_eigenvalues,
                          n_zero_eigenvalues) {
  svdX <- RSpectra::svds(K, k = n_nonzero_eigenvalues)
  U_K <- svdX$u
  Lambda <- svdX$d^2
  if (any(Lambda < 1e-5)) {
    Lambda[which(Lambda < 1e-5)] <- 1e-05
  }

  structure(list(
    x = x, y = y,
    U = U_K,
    D = Lambda,
    n_nonzero_eigenvalues = n_nonzero_eigenvalues,
    n_zero_eigenvalues = n_zero_eigenvalues
  ),
  class = c("lowrank")
  )
}

#' @rdname ggmix_data_object
new_lowrank_UD <- function(x, y, U, D,
                           n_nonzero_eigenvalues,
                           n_zero_eigenvalues) {
  structure(list(
    x = x, y = y,
    U = U,
    D = D,
    n_nonzero_eigenvalues = n_nonzero_eigenvalues,
    n_zero_eigenvalues = n_zero_eigenvalues
  ),
  class = c("lowrank")
  )
}
