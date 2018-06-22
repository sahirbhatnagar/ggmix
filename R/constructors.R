#' Constructor functions for the different ggmix objects
#'
#' \code{new_fullrank_kinship}, \code{new_fullrank_K}, \code{new_fullrank_UD},
#' \code{new_lowrank_kinship}, \code{new_lowrank_K} and \code{new_lowrank_UD}
#' create the ggmix objects from the provided data that are necessary to fit the
#' penalized linear mixed model according to the user's parameters.
#'
#' @seealso \code{\link{ggmix}}
#' @inheritParams ggmix
#' @param n_zero_eigenvalues the number of desired or specifed zero eigenvalues.
#'   This is only needed when \code{estimation="lowrank"}, and is calculated
#'   internally by the \code{\link{ggmix}} function. It is equal to the number
#'   of observations minus \code{n_nonzero_eigenvalues}
#'
#' @return A ggmix object, of the class that corresponds to the estimation
#'   method. These objects are lists that contain the data necessary for
#'   computation.
#' @name ggmix_data_object
NULL

#' @rdname ggmix_data_object
#' @export
new_fullrank_kinship <- function(kinship) {

  phi_eigen <- eigen(kinship)
  U_kinship <- phi_eigen$vectors
  Lambda <- phi_eigen$values
  if (any(Lambda < 1e-5))
    Lambda[which(Lambda < 1e-5)] <- 1e-05

  structure(list(U = U_kinship,
                 D = Lambda),
            class = c("fullrank"))
}


#' @rdname ggmix_data_object
#' @export
new_fullrank_K <- function(K) {

  svdX <- svd(K)
  U_K <- svdX$u
  Lambda <- svdX$d^2
  if (any(Lambda < 1e-5))
    Lambda[which(Lambda < 1e-5)] <- 1e-05

  structure(list(U = U_K,
                 D = Lambda),
            class = c("fullrank"))
}

#' @rdname ggmix_data_object
#' @export
new_fullrank_UD <- function(U, D) {

  structure(list(U = U,
                 D = D),
            class = c("fullrank"))
}



#' @rdname ggmix_data_object
#' @export
new_lowrank_kinship <- function(kinship,
                                n_nonzero_eigenvalues,
                                n_zero_eigenvalues) {

  phi_eigen <- RSpectra::eigs(kinship, k = n_nonzero_eigenvalues)
  U_kinship <- phi_eigen$vectors
  Lambda <- phi_eigen$values
  if (any(Lambda < 1e-5))
    Lambda[which(Lambda < 1e-5)] <- 1e-05

  structure(list(U = U_kinship,
                 D = Lambda,
                 n_nonzero_eigenvalues = n_nonzero_eigenvalues,
                 n_zero_eigenvalues = n_zero_eigenvalues),
            class = c("lowrank"))
}


#' @rdname ggmix_data_object
#' @export
new_lowrank_K <- function(K,
                          n_nonzero_eigenvalues,
                          n_zero_eigenvalues) {

  svdX <- RSpectra::svds(K, k = n_nonzero_eigenvalues)
  U_K <- svdX$u
  Lambda <- svdX$d^2
  if (any(Lambda < 1e-5))
    Lambda[which(Lambda < 1e-5)] <- 1e-05

  structure(list(U = U_K,
                 D = Lambda,
                 n_nonzero_eigenvalues = n_nonzero_eigenvalues,
                 n_zero_eigenvalues = n_zero_eigenvalues),
            class = c("lowrank"))
}


#' @rdname ggmix_data_object
#' @export
new_lowrank_UD <- function(U, D,
                           n_nonzero_eigenvalues,
                           n_zero_eigenvalues) {

  structure(list(U = U,
                 D = D,
                 n_nonzero_eigenvalues = n_nonzero_eigenvalues,
                 n_zero_eigenvalues = n_zero_eigenvalues),
            class = c("lowrank"))
}

