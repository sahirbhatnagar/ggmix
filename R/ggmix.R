#' @title Fit Linear Mixed Model with Lasso or Group Lasso Penalty
#' @description Main function to fit the linear mixed model with lasso or group
#'   lasso penalty for a sequence of tuning parameters. This is a penalized
#'   regression method that accounts for population structure using either the
#'   kinship matrix or the factored realized relationship matrix
#' @export
ggmix <- function(x, y,
                  U, D,
                  kinship,
                  K,
                  n_nonzero_eigenvalues,
                  estimation = c("full", "low"),
                  penalty = c("lasso", "gglasso"),
                  group,
                  penalty.factor,
                  lambda = NULL,
                  lambda_min_ratio  = ifelse(n < p, 0.01, 0.0001),
                  nlambda = 100,
                  eta_init = 0.5,
                  maxit = 100,
                  fdev = 1e-5,
                  standardize = FALSE,
                  alpha = 1, # elastic net mixing param. 1 is lasso, 0 is ridge
                  thresh_glmnet = 1e-8, # this is for glmnet
                  epsilon = 1e-4 # this is for ggmix
                  ) {

  this.call <- match.call()

  # Check input -------------------------------------------------------------

  estimation <- tryCatch(match.arg(estimation),
                         error = function(c) {
                           stop(strwrap("Estimation method should be
                                        \"full_rank\" or \"low_rank\""),
                                call. = FALSE)
                         })

  penalty <- tryCatch(match.arg(penalty),
                        error = function(c) {
                          stop(strwrap("Inference method should be \"lasso\" or
                                       \"group_lasso\""),
                               call. = FALSE)
                        })

  if (!is.matrix(x))
    stop("x has to be a matrix")
  if (any(is.na(x)))
    stop("Missing values in x not allowed.")
  if (any(is.na(y)))
    stop("Missing values in y not allowed.")
  y <- drop(y)
  if (!is.numeric(y))
    stop("y must be numeric. Factors are not allowed.")

  if ((!missing(U) & missing(D)) | (!missing(D) & missing(U)))
    stop("both U and D must be specified.")

  if (estimation == "low" & missing(rank))
    stop("the rank argument must be specified with low rank estimation")

  if (penalty == "gglasso" & missing(group))
    stop(strwrap("group cannot be missing when using the group lasso
                 penalty"))

  np_design <- dim(x)
  if (is.null(np_design) | (np_design[2] <= 1))
    stop("x should be a matrix with 2 or more columns")

  n_design <- np_design[[1]]
  p_design <- np_design[[2]]

  is_kinship <- !missing(kinship)
  is_UD <- !missing(U) & !missing(D)
  is_K <- !missing(K)

  if (all(is_kinship, is_UD, is_K)) {
    warning(strwrap("kinship, U, D and K arguments have all been specified.
                    Only one of either kinship, U and D, or K need to be
                    specified. Ignoring U, D and K. Only using kinship in the
                    analysis"))
    is_UD <- FALSE
    is_K <- FALSE
  }

  if (all(is_kinship, is_UD, !is_K) | all(is_kinship, !is_UD, is_K)) {
    warning(strwrap("more than one of kinship, U, D or K arguments have
                    been specified. Only one of either kinship, U and D,
                    or K need to be specified. Only using kinship in the
                    analysis."))
    is_UD <- FALSE
    is_K <- FALSE
  }

  if (!is_kinship) {
    if (all(is_UD, is_K)) {
      warning(strwrap("U, D and K arguments have been specified. Only one of
                    either U and D, or K need to be specified. Only using
                    U and D in the analysis and ignoring K."))
      is_K <- FALSE
    }
  }

  if (is_kinship) {

    if (!is.matrix(kinship))
      stop("kinship has to be a matrix")
    np_kinship <- dim(kinship)
    n_kinship <- np_kinship[[1]]
    p_kinship <- np_kinship[[2]]

    if (n_kinship != p_kinship)
      stop("kinship should be a square matrix")
    if (n_kinship != n_design)
      stop(strwrap("number of rows in kinship matrix must equal the
                   number of rows in x matrix"))

    corr_type <- "kinship"
  }

  if (is_K) {
    if (!is.matrix(K))
      stop("K has to be a matrix")
    np_K <- dim(K)
    n_K <- np_K[[1]]
    p_K <- np_K[[2]]

    if (n_K != n_design)
      stop(strwrap("number of rows in K matrix must equal the
                   number of rows in x matrix"))

    corr_type <- "Kmat"
  }

  if (is_UD) {
    if (!is.matrix(U))
      stop("U has to be a matrix")
    np_U <- dim(U)
    n_U <- np_U[[1]]
    p_U <- np_U[[2]]

    D <- drop(D)
    if (!is.numeric(D)) stop(strwrap("D must be numeric"))
    p_D <- length(D)

    if (n_U != n_design)
      stop(strwrap("number of rows in U matrix must equal the
                   number of rows in x matrix"))
    if (p_U != p_D)
      stop(strwrap("Length of D should equal the number of columns of U."))

    corr_type <- "UD"
  }

  vnames <- colnames(x)
  if (is.null(vnames)) {
    vnames <- paste("V", seq(p), sep = "")
    colnames(x) <- vnames
  }

  if (!missing(penalty.factor)) {
    if (length(penalty.factor) != p_design)
      stop(strwrap("length of penalty.factor should be equal to number
                   of columns in x."))
  }


  # Create ggmix objects ----------------------------------------------------

  if (estimation == "full") {
    ggmix_data_object <- switch(corr_type,
                           kinship = new_fullrank_kinship(kinship = kinship),
                           Kmat = new_fullrank_Kmat(Kmat = K),
                           UD = new_fullrank_UD(U = U, D = D)
    )
  }

  if (estimation == "low") {

    if (missing(n_nonzero_eigenvalues) | n_nonzero_eigenvalues <= 0)
      stop(strwrap("n_nonzero_eigenvalues can't be missing, equal to 0 or negative
                 when low rank estimation is specified (estimation=\"low\").
                 Use estimation=\"full\" if you want all eigenvalues to be
                 computed."))

    if (!requireNamespace("RSpectra", quietly = TRUE)) {
      stop(strwrap("Package \"RSpectra\" needed when low rank estimation is
                   specified (estimation=\"low\"). Please install it."),
           call. = FALSE
      )
    }

    n_zero_eigenvalues <- n_design - n_nonzero_eigenvalues

    ggmix_data_object <- switch(corr_type,
                                kinship = new_lowrank_kinship(kinship = kinship,
                                                              nnzeigen = n_nonzero_eigenvalues,
                                                              nzeigen = n_zero_eigenvalues),
                                Kmat = new_lowrank_Kmat(K = K,
                                                        nnzeigen = n_nonzero_eigenvalues,
                                                        nzeigen = n_zero_eigenvalues),
                                UD = new_lowrank_UD(U = U, D = D,
                                                    nnzeigen = n_nonzero_eigenvalues,
                                                    nzeigen = n_zero_eigenvalues)
    )
  }






  }
