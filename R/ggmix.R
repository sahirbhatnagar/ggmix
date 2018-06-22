#' @title Fit Linear Mixed Model with Lasso or Group Lasso Regularization
#' @description Main function to fit the linear mixed model with lasso or group
#'   lasso penalty for a sequence of tuning parameters. This is a penalized
#'   regression method that accounts for population structure using either the
#'   kinship matrix or the factored realized relationship matrix
#' @param x input matrix, of dimension n x p; where n is the number of
#'   observations and p are the number of predictors.
#' @param y response variable. must be a quantitative variable
#' @param D eigenvalues
#' @param U left singular vectors corresponding to the eigenvalues provided in
#'   the \code{D} argument
#' @param kinship positive definite kinship matrix
#' @param K the matrix of SNPs used to determine the kinship matrix
#' @param n_nonzero_eigenvalues the number of nonzero eigenvalues. This argument
#'   is only used when \code{estimation="low"} and either \code{kinship} or
#'   \code{K} is provided. This argument will limit the function to finding the
#'   \code{n_nonzero_eigenvalues} largest eigenvalues.
#' @param estimation type of estimation
#' @param penalty type of regularization penalty. if \code{penalty="gglasso"}
#'   then the \code{group} argument must also be specified
#' @param group a vector of consecutive integers describing the grouping of the
#'   coefficients
#' @param penalty.factor Separate penalty factors can be applied to each
#'   coefficient. This is a number that multiplies lambda to allow differential
#'   shrinkage. Can be 0 for some variables, which implies no shrinkage, and
#'   that variable is always included in the model. Default is 1 for all
#'   variables
#' @param lambda A user supplied lambda sequence (this is the tuning parameter).
#'   Typical usage is to have the program compute its own lambda sequence based
#'   on nlambda and lambda.min.ratio. Supplying a value of lambda overrides
#'   this. WARNING: use with care. Do not supply a single value for lambda (for
#'   predictions after CV use predict() instead). Supply instead a decreasing
#'   sequence of lambda values. glmnet relies on its warms starts for speed, and
#'   its often faster to fit a whole path than compute a single fit.
#' @param lambda_min_ratio Smallest value for lambda, as a fraction of
#'   lambda.max, the (data derived) entry value (i.e. the smallest value for
#'   which all coefficients are zero). The default depends on the sample size
#'   nobs relative to the number of variables nvars. If nobs > nvars, the
#'   default is 0.0001, close to zero. If nobs < nvars, the default is 0.01. A
#'   very small value of lambda.min.ratio will lead to a saturated fit in the
#'   nobs < nvars case.
#' @param nlambda the number of lambda values - default is 100.
#' @param eta_init initial value for the eta parameter, with \eqn{0 < \eta < 1}
#'   used in determining lambda.max and starting value for fitting algorithm.
#' @param maxit Maximum number of passes over the data for all lambda values;
#'   default is 10^2.
#' @param fdev Fractional deviance change theshold. If change in deviance
#'   between adjacent lambdas is less than fdev, the solution path stops early.
#'   factory default = 1.0e-5
#' @param alpha The elasticnet mixing parameter, with \eqn{0 \leq \alpha \leq
#'   1}. alpha=1 is the lasso penalty, and alpha=0 the ridge penalty.
#' @param thresh_glmnet Convergence threshold for coordinate descent for
#'   updating beta parameters. Each inner coordinate-descent loop continues
#'   until the maximum change in the objective after any coefficient update is
#'   less than thresh times the null deviance. Defaults value is 1E-7
#' @param standardize Logical flag for x variable standardization, prior to
#'   fitting the model sequence. The coefficients are always returned on the
#'   original scale. Default is standardize=FALSE. If variables are in the same
#'   units already, you might not wish to standardize.
#' @param epsilon Convergence threshold for block relaxation of the entire
#'   parameter vector \eqn{\Theta = ( \beta, \eta, \sigma^2 )}. The algorithm
#'   converges when \deqn{crossprod(\Theta_{j+1} - \Theta_{j}) < \epsilon}.
#'   Defaults value is 1E-7
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

    corr_type <- "K"
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
                           Kmat = new_fullrank_K(K = K),
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
                                  n_nonzero_eigenvalues = n_nonzero_eigenvalues,
                                  n_zero_eigenvalues = n_zero_eigenvalues),
                                K = new_lowrank_K(K = K,
                                  n_nonzero_eigenvalues = n_nonzero_eigenvalues,
                                  n_zero_eigenvalues = n_zero_eigenvalues),
                                UD = new_lowrank_UD(U = U, D = D,
                                  n_nonzero_eigenvalues = n_nonzero_eigenvalues,
                                  n_zero_eigenvalues = n_zero_eigenvalues)
    )
  }

  browser()





  }
