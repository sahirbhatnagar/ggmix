#' Constructor functions for the different ggmix objects
#'
#' \code{PcevClassical}, \code{PcevBlock} and \code{PcevSingular} create the
#' pcev objects from the provided data that are necessary to compute the PCEV
#' according to the user's parameters.
#'
#' @seealso \code{\link{estimatePcev}}, \code{\link{computePCEV}}
#' @param response A matrix of response variables.
#' @param covariate A matrix or a data frame of covariates.
#' @param confounder A matrix or data frame of confounders
#' @return A pcev object, of the class that corresponds to the estimation
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
new_fullrank_Kmat <- function(Kmat) {

  svdX <- svd(Kmat)
  U_Kmat <- svdX$u
  Lambda <- svdX$d^2
  if (any(Lambda < 1e-5))
    Lambda[which(Lambda < 1e-5)] <- 1e-05

  structure(list(U = U_Kmat,
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
new_lowrank_kinship <- function(kinship, nnzeigen, nzeigen) {

  phi_eigen <- RSpectra::eigs(kinship, k = nnzeigen)
  U_kinship <- phi_eigen$vectors
  Lambda <- phi_eigen$values
  if (any(Lambda < 1e-5))
    Lambda[which(Lambda < 1e-5)] <- 1e-05

  structure(list(U = U_kinship,
                 D = Lambda,
                 nnzeigen = nnzeigen,
                 nzeigen = nzeigen),
            class = c("lowrank"))
}


#' @rdname ggmix_data_object
#' @export
new_lowrank_Kmat <- function(Kmat, nnzeigen, nzeigen) {

  svdX <- RSpectra::svds(Kmat, k = nnzeigen)
  U_Kmat <- svdX$u
  Lambda <- svdX$d^2
  if (any(Lambda < 1e-5))
    Lambda[which(Lambda < 1e-5)] <- 1e-05

  structure(list(U = U_Kmat,
                 D = Lambda,
                 nnzeigen = nnzeigen,
                 nzeigen = nzeigen),
            class = c("lowrank"))
}


#' @rdname ggmix_data_object
#' @export
new_lowrank_UD <- function(U, D, nnzeigen, nzeigen) {

  structure(list(U = U,
                 D = D,
                 nnzeigen = nnzeigen,
                 nzeigen = nzeigen),
            class = c("lowrank"))
}

