#' @title Fit Linear Mixed Model with Lasso or Group Lasso Regularization
#' @description Main function to fit the linear mixed model with lasso or group
#'   lasso penalty for a sequence of tuning parameters. This is a penalized
#'   regression method that accounts for population structure using either the
#'   kinship matrix or the factored realized relationship matrix
#' @param x input matrix, of dimension n x p; where n is the number of
#'   observations and p are the number of predictors.
#' @param y response variable. must be a quantitative variable
#' @param D non-zero eigenvalues. This option is provided to the user should
#'   they decide or need to calculate the eigen decomposition of the kinship
#'   matrix or the singular value decomposition of the matrix of SNPs used to
#'   calculate the kinship outside of this function. This may occur, if for
#'   example, it is easier (e.g. because of memory issues, it's easier to
#'   calculate in plink). This should correspond to the non-zero eigenvalues
#'   only. Note that if you are doing an \code{svd} on the matrix of SNPs used
#'   to calculate the kinship matrix, then you must provide the square of the
#'   singular values so that they correspond to the eigenvalues of the kinship
#'   matrix. If you want to use the low rank estimation algorithm, you must
#'   provide the truncated eigenvalues and eigenvectors to the \code{D} and
#'   \code{U} arguments, respectively. If you want \code{ggmix} to truncate the
#'   eigenvectors and eigenvalues for low rank estimation, then provide either
#'   \code{K} or \code{kinship} instead and specify
#'   \code{n_nonzero_eigenvalues}.
#' @param U left singular vectors corresponding to the non-zero eigenvalues
#'   provided in the \code{D} argument.
#' @param kinship positive definite kinship matrix
#' @param K the matrix of SNPs used to determine the kinship matrix
#' @param n_nonzero_eigenvalues the number of nonzero eigenvalues. This argument
#'   is only used when \code{estimation="low"} and either \code{kinship} or
#'   \code{K} is provided. This argument will limit the function to finding the
#'   \code{n_nonzero_eigenvalues} largest eigenvalues. If \code{U} and \code{D}
#'   have been provided, then \code{n_nonzero_eigenvalues} defaults to the
#'   length of \code{D}.
#' @param n_zero_eigenvalues Currently not being used. Represents the number of
#'   zero eigenvalues. This argument must be specified when \code{U} and
#'   \code{D} are specified and \code{estimation="low"}. This is required for
#'   low rank estimation because the number of zero eigenvalues and their
#'   corresponding eigenvalues appears in the likelihood. In general this would
#'   be the rank of the matrix used to calculate the eigen or singular value
#'   decomposition. When \code{kinship} is provided and \code{estimation="low"}
#'   the default value will be \code{nrow(kinship) - n_nonzero_eigenvalues}.
#'   When \code{K} is provided and \code{estimation="low"}, the default value is
#'   \code{rank(K) - n_nonzero_eigenvalues}
#' @param estimation type of estimation. Currently only \code{type="full"} has
#'   been implemented.
#' @param penalty type of regularization penalty. Currently, only
#'   penalty="lasso" has been implemented.
#' @param group a vector of consecutive integers describing the grouping of the
#'   coefficients. Currently not implemented, but will be used when
#'   penalty="gglasso" is implemented.
#' @param dfmax limit the maximum number of variables in the model. Useful for
#'   very large \code{p} (the total number of predictors in the design matrix),
#'   if a partial path is desired. Default is the number of columns in the
#'   design matrix + 2 (for the variance components)
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
#' @param fdev Fractional deviance change threshold. If change in deviance
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
#' @param verbose display progress. Can be either 0,1 or 2. 0 will not display
#'   any progress, 2 will display very detailed progress and 1 is somewhere in
#'   between. Default: 0.
#'
#' @examples
#' data(admixed)
#' fitlmm <- ggmix(x = admixed$xtrain, y = admixed$ytrain,
#'                 kinship = admixed$kin_train,
#'                 estimation = "full")
#' gicfit <- gic(fitlmm)
#' coef(gicfit, type = "nonzero")
#' predict(gicfit, newx = admixed$xtest)[1:5,,drop=FALSE]
#' plot(gicfit)
#' plot(fitlmm)
#'
#' @references Bhatnagar, Sahir R, Yang, Yi, Lu, Tianyuan, Schurr, Erwin,
#'   Loredo-Osti, JC, Forest, Marie, Oualkacha, Karim, and Greenwood, Celia MT.
#'   (2020) \emph{Simultaneous SNP selection and adjustment for population
#'   structure in high dimensional prediction models}
#'   \url{https://doi.org/10.1101/408484}
#'
#'   Friedman, J., Hastie, T. and Tibshirani, R. (2008) \emph{Regularization
#'   Paths for Generalized Linear Models via Coordinate Descent},
#'   \url{http://www.stanford.edu/~hastie/Papers/glmnet.pdf} \emph{Journal of
#'   Statistical Software, Vol. 33(1), 1-22 Feb 2010}
#'   \url{http://www.jstatsoft.org/v33/i01/}
#'
#'   Yang, Y., & Zou, H. (2015). A fast unified algorithm for solving
#'   group-lasso penalize learning problems. \emph{Statistics and Computing},
#'   25(6), 1129-1141.
#'   \url{http://www.math.mcgill.ca/yyang/resources/papers/gglasso.pdf}
#' @export
ggmix <- function(x, y,
                  U, D,
                  kinship,
                  K,
                  n_nonzero_eigenvalues,
                  n_zero_eigenvalues,
                  estimation = c("full"),
                  penalty = c("lasso"),
                  group,
                  penalty.factor = rep(1, p_design),
                  lambda = NULL,
                  lambda_min_ratio = ifelse(n_design < p_design, 0.01, 0.0001),
                  nlambda = 100,
                  eta_init = 0.5,
                  maxit = 100,
                  fdev = 1e-20,
                  standardize = FALSE,
                  alpha = 1, # elastic net mixing param. 1 is lasso, 0 is ridge
                  thresh_glmnet = 1e-8, # this is for glmnet
                  epsilon = 1e-4, # this is for ggmix
                  dfmax = p_design + 2,
                  verbose = 0) {
  this.call <- match.call()

  # Check input -------------------------------------------------------------

  estimation <- tryCatch(match.arg(estimation),
    error = function(c) {
      stop(strwrap("Estimation method should be
                   \"full\""),
        call. = FALSE
      )
    }
  )

  penalty <- tryCatch(match.arg(penalty),
    error = function(c) {
      stop(strwrap("Penalty type should be \"lasso\""),
        call. = FALSE
      )
    }
  )

  if (!is.matrix(x)) {
    stop("x has to be a matrix")
  }
  if (any(is.na(x))) {
    stop("Missing values in x not allowed.")
  }
  if (any(is.na(y))) {
    stop("Missing values in y not allowed.")
  }
  y <- drop(y)
  if (!is.numeric(y)) {
    stop("y must be numeric. Factors are not allowed.")
  }

  if ((!missing(U) & missing(D)) | (!missing(D) & missing(U))) {
    stop("both U and D must be specified.")
  }

  if (penalty == "gglasso" & missing(group)) {
    stop(strwrap("group cannot be missing when using the group lasso
                 penalty"))
  }

  np_design <- dim(x)
  if (is.null(np_design) | (np_design[2] <= 1)) {
    stop("x should be a matrix with 2 or more columns")
  }

  # note that p_design doesn't contain the intercept
  # whereas the x in the ggmix_object will have the intercept
  n_design <- np_design[[1]]
  p_design <- np_design[[2]]

  dfmax <- as.double(dfmax)

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

  if (!missing(n_nonzero_eigenvalues)) {
    n_nonzero_eigenvalues <- as.double(n_nonzero_eigenvalues)
  }

  if (!missing(n_zero_eigenvalues)) {
    n_zero_eigenvalues <- as.double(n_zero_eigenvalues)
  }


  # check kinship input -----------------------------------------------------

  if (is_kinship) {
    if (!is.matrix(kinship)) {
      stop("kinship has to be a matrix")
    }
    np_kinship <- dim(kinship)
    n_kinship <- np_kinship[[1]]
    p_kinship <- np_kinship[[2]]

    if (n_kinship != p_kinship) {
      stop("kinship should be a square matrix")
    }
    if (n_kinship != n_design) {
      stop(strwrap("number of rows in kinship matrix must equal the
                   number of rows in x matrix"))
    }

    if (estimation == "low") {
      if (missing(n_nonzero_eigenvalues)) {
        stop(strwrap("when kinship is specified and estimation=\"low\",
                     n_nonzero_eigenvalues must be specified."))
      }
    }

    if (missing(n_zero_eigenvalues) & estimation == "low") {
      n_zero_eigenvalues <- nrow(kinship) - n_nonzero_eigenvalues
      warning(sprintf(
        "n_zero_eigenvalues missing. setting to %g",
        n_zero_eigenvalues
      ))
    }

    corr_type <- "kinship"
  }


  # check K matrix input ----------------------------------------------------

  if (is_K) {
    if (!is.matrix(K)) {
      stop("K has to be a matrix")
    }
    np_K <- dim(K)
    n_K <- np_K[[1]]
    p_K <- np_K[[2]]

    if (n_K != n_design) {
      stop(strwrap("number of rows in K matrix must equal the
                   number of rows in x matrix"))
    }

    if (estimation == "low") {
      if (missing(n_nonzero_eigenvalues)) {
        stop(strwrap("when K is specified and estimation=\"low\",
                     n_nonzero_eigenvalues must be specified."))
      }
    }

    if (missing(n_zero_eigenvalues) & estimation == "low") {
      n_zero_eigenvalues <- min(n_K, p_K) - n_nonzero_eigenvalues
      warning(sprintf(
        "n_zero_eigenvalues missing. setting to %g",
        n_zero_eigenvalues
      ))
    }

    corr_type <- "K"
  }


  # check U and D input -----------------------------------------------------

  if (is_UD) {
    if (!is.matrix(U)) {
      stop("U has to be a matrix")
    }
    if (missing(n_zero_eigenvalues) & estimation == "low") {
      stop(strwrap("n_zero_eigenvalues must be specified when U and D have
                   been provided and estimation=\"low\". In general this would
                   be the rank of the matrix used to calculate the eigen or
                   singular value decomposition."))
    }
    np_U <- dim(U)
    n_U <- np_U[[1]]
    p_U <- np_U[[2]]

    D <- drop(D)
    if (!is.numeric(D)) stop(strwrap("D must be numeric"))
    p_D <- length(D)

    if (n_U != n_design) {
      stop(strwrap("number of rows in U matrix must equal the
                   number of rows in x matrix"))
    }
    if (p_U != p_D) {
      stop(strwrap("Length of D should equal the number of columns of U."))
    }

    if (missing(n_nonzero_eigenvalues)) {
      n_nonzero_eigenvalues <- p_D
    }

    corr_type <- "UD"
  }

  vnames <- colnames(x)
  if (is.null(vnames)) {
    vnames <- paste("V", seq(p_design), sep = "")
    colnames(x) <- vnames
  }

  if (!missing(penalty.factor)) {
    penalty.factor <- as.double(penalty.factor)
    if (length(penalty.factor) != p_design) {
      stop(strwrap("length of penalty.factor should be equal to number
                   of columns in x."))
    }
  }


  # Create ggmix objects ----------------------------------------------------

  if (estimation == "full") {
    ggmix_data_object <- switch(corr_type,
      kinship = new_fullrank_kinship(x = x, y = y, kinship = kinship),
      K = new_fullrank_K(x = x, y = y, K = K),
      UD = new_fullrank_UD(x = x, y = y, U = U, D = D)
    )
  }

  if (estimation == "low") {
    if (missing(n_nonzero_eigenvalues) | n_nonzero_eigenvalues <= 0) {
      stop(strwrap("n_nonzero_eigenvalues can't be missing, equal to 0 or negative
                 when low rank estimation is specified (estimation=\"low\").
                 Use estimation=\"full\" if you want all eigenvalues to be
                 computed."))
    }

    if (!requireNamespace("RSpectra", quietly = TRUE)) {
      stop(strwrap("Package \"RSpectra\" needed when low rank estimation is
                   specified (estimation=\"low\"). Please install it."),
        call. = FALSE
      )
    }

    ggmix_data_object <- switch(corr_type,
      kinship = new_lowrank_kinship(
        x = x, y = y,
        kinship = kinship,
        n_nonzero_eigenvalues = n_nonzero_eigenvalues,
        n_zero_eigenvalues = n_zero_eigenvalues
      ),
      K = new_lowrank_K(
        x = x, y = y, K = K,
        n_nonzero_eigenvalues = n_nonzero_eigenvalues,
        n_zero_eigenvalues = n_zero_eigenvalues
      ),
      UD = new_lowrank_UD(
        x = x, y = y, U = U, D = D,
        n_nonzero_eigenvalues = n_nonzero_eigenvalues,
        n_zero_eigenvalues = n_zero_eigenvalues
      )
    )
  }


  # fit linear mixed model --------------------------------------------------

  # browser()
  if (penalty == "lasso") {
    fit <- lmmlasso(ggmix_data_object,
      penalty.factor = penalty.factor,
      lambda = lambda,
      lambda_min_ratio = lambda_min_ratio,
      nlambda = nlambda,
      n_design = n_design,
      p_design = p_design,
      eta_init = eta_init,
      maxit = maxit,
      fdev = fdev,
      standardize = standardize,
      alpha = alpha, # elastic net mixing param. 1 is lasso, 0 is ridge
      thresh_glmnet = thresh_glmnet, # this is for glmnet
      epsilon = epsilon,
      dfmax = dfmax,
      verbose = verbose
    )
  } else if (penalty == "gglasso") {
    # not yet implemented
  }

  fit$call <- this.call

  return(fit)
}
