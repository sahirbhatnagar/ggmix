#' Simulation Scenarion from Bhatnagar et al. (2018+) ggmix paper
#'
#' @description Function that generates data of the different simulation studies
#'   presented in the accompanying paper. This function requires the
#'   \code{popkin} and \code{bnpsd} package to be installed.
#' @param n number of observations to simulate
#' @param p_test number of variables in X_test, i.e., the design matrix
#' @param p_kinship number of variable in X_kinship, i.e., matrix used to
#'   calculate kinship
#' @param k number of intermediate subpopulations.
#' @param s the desired bias coefficient, which specifies sigma indirectly.
#'   Required if sigma is missing
#' @param Fst The desired final FST of the admixed individuals. Required if
#'   sigma is missing
#' @param b0 the true intercept parameter
#' @param eta the true eta parameter, which has to be \code{0 < eta < 1}
#' @param sigma2 the true sigma2 parameter
#' @param nPC number of principal components to include in the design matrix
#'   used for regression adjustment for population structure via principal
#'   components. This matrix is used as the input in a standard lasso regression
#'   routine, where there are no random effects.
#' @param geography the type of geography for simulation the kinship matrix.
#'   "ind" is independent populations where every individuals is actually
#'   unadmixed, "1d" is a 1D geography and "circ" is circular geography.
#'   Default: "ind". See the functions in the \code{bnpsd} for details on how
#'   this data is actually generated.
#' @param percent_causal percentage of \code{p_test} that is causal. must be
#'   \eqn{0 \leq percent_causal \leq 1}. The true regression coefficients are
#'   generated from a standard normal distribution.
#' @param percent_overlap this represents the percentage of causal SNPs that
#'   will also be included in the calculation of the kinship matrix
#' @details The kinship is estimated using the \code{popkin} function from the
#'   \code{popkin} package. This function will multiple that kinship matrix by 2
#'   to give the expected covariance matrix which is subsequently used in the
#'   linear mixed models
#' @return A list with the following elements \describe{\item{y}{simulated
#'   response vector} \item{x}{simulated design matrix} \item{causal}{character
#'   vector of the names of the causal SNPs} \item{beta}{the vector of true
#'   regression coefficients} \item{kin}{2 times the estimated kinship}
#'   \item{Xkinship}{the matrix of SNPs used to estimate the kinship matrix}
#'   \item{not_causal}{character vector of the non-causal SNPs}
#'   \item{causal_positive}{character vector of the causal SNPs with positive
#'   regression coefficient} \item{causal_negative}{character vector of the
#'   causal SNPs with negative regression coefficient}\item{x_lasso}{the design
#'   matrix which also includes \code{nPC} principal components} }
#' @seealso \code{\link[bnpsd]{q1d}},\code{\link[bnpsd]{qis}},
#'   \code{\link[bnpsd]{q1dc}}, \code{\link[bnpsd]{rbnpsd}}
gen_structured_model <- function(n, p_test, p_kinship, k, s, Fst, b0, nPC = 10,
                                 eta, sigma2, geography = c("ind", "1d", "circ"),
                                 percent_causal, percent_overlap) {

  # p_test:
  # p_kinship:
  # k:	N
  # s:
  # F: The length-k vector of inbreeding coefficients (or FST's) of the intermediate subpopulations,
  # up to a scaling factor (which cancels out in calculations). Required if sigma is missing
  # Fst: T
  # browser()
  # define population structure

  if (!requireNamespace("bnpsd", quietly = TRUE)) {
    stop(strwrap("Package \"bnpsd\" needed to simulate data.
                 Please install it."),
      call. = FALSE
    )
  }

  if (!requireNamespace("popkin", quietly = TRUE)) {
    stop(strwrap("Package \"popkin\" needed to simulate data.
                 Please install it."),
      call. = FALSE
    )
  }


  # FF <- 1:k # subpopulation FST vector, up to a scalar
  # s <- 0.5 # desired bias coefficient
  # Fst <- 0.1 # desired FST for the admixed individuals
  geography <- match.arg(geography)
  if (geography == "1d") {
    FF <- 1:k # subpopulation FST vector, up to a scalar
    obj <- bnpsd::admix_prop_1d_linear(n_ind = n,
                                       k_subpops = k,
                                       bias_coeff = s,
                                       coanc_subpops = FF,
                                       fst = Fst)
    # Q <- obj$Q
    # FF <- obj$F
    admix_proportions <- obj$admix_proportions
    # rescaled inbreeding vector for intermediate subpopulations
    inbr_subpops <- obj$coanc_subpops

    # get pop structure parameters of the admixed individuals
    coancestry <- coanc_admix(admix_proportions, inbr_subpops)
    kinship <- coanc_to_kinship(coancestry)

  } else if (geography == "ind") {
    ngroup <- n / k # equal sized groups
    # here’s the labels (for simplicity, list all individuals of S1 first, then S2, then S3)
    labs <- rep(paste0("S", 1:k), each = ngroup)
    # data dimensions infered from labs:
    length(labs) # number of individuals "n"
    # desired admixture matrix ("is" stands for "Independent Subpopulations")
    # number of subpopulations "k_subpops"
    k_subpops <- length(unique(labs))

    # desired admixture matrix
    admix_proportions <- admix_prop_indep_subpops(labs)

    # subpopulation FST vector, unnormalized so far
    inbr_subpops <- 1 : k_subpops
    # normalized to have the desired FST
    # NOTE fst is a function in the `popkin` package
    inbr_subpops <- inbr_subpops / popkin::fst(inbr_subpops) * Fst
    # verify FST for the intermediate subpopulations
    # fst(inbr_subpops)
    #> [1] 0.2

    # get coancestry of the admixed individuals
    coancestry <- coanc_admix(admix_proportions, inbr_subpops)
    # before getting FST for individuals, weigh then inversely proportional to subpop sizes
    weights <- popkin::weights_subpops(labs) # function from `popkin` package

    kinship <- coanc_to_kinship(coancestry)

  } else if (geography == "circ") {
    FF <- 1:k # subpopulation FST vector, up to a scalar
    # obj <- bnpsd::admix_prop_1d_circular(n_ind = n, k_subpops = k, s = s, F = FF, Fst = Fst)
    # Q <- obj$Q
    # FF <- obj$F

    # admixture proportions from *circular* 1D geography
    obj <- admix_prop_1d_circular(
      n_ind = n,
      k_subpops = k,
      bias_coeff = s,
      coanc_subpops = FF,
      fst = Fst
    )
    admix_proportions <- obj$admix_proportions
    inbr_subpops <- obj$coanc_subpops

    # get pop structure parameters of the admixed individuals
    coancestry <- coanc_admix(admix_proportions, inbr_subpops)

    kinship <- coanc_to_kinship(coancestry)
  }

  ncausal <- p_test * percent_causal
  # browser()
  if (percent_overlap == "100") {
    total_snps_to_simulate <- p_test + p_kinship - ncausal
    # this contains all SNPs (X_{Testing}:X_{kinship})
    out <- bnpsd::rbnpsd(Q, FF, total_snps_to_simulate)
    Xall <- t(out$X) # genotypes are columns, rows are subjects
    cnames <- paste0("X", 1:total_snps_to_simulate)
    colnames(Xall) <- cnames
    rownames(Xall) <- paste0("id", 1:n)
    Xall[1:5, 1:5]
    dim(Xall)
    subpops <- ceiling((1:n) / n * k)
    table(subpops) # got k=10 subpops with 100 individuals each

    # Snps used for kinship
    snps_kinships <- sample(cnames, p_kinship, replace = FALSE)
    length(snps_kinships)

    # all causal snps are in kinship matrix
    causal <- sample(snps_kinships, ncausal, replace = FALSE)

    snps_test <- c(setdiff(cnames, snps_kinships), causal)
    # length(snps_test)
    # setdiff(cnames, snps_kinships) %>% length()
    not_causal <- setdiff(snps_test, causal)

    Xkinship <- Xall[, snps_kinships]
    Xtest <- Xall[, snps_test]

    # now estimate kinship using popkin
    # PhiHat <- popkin::popkin(X, subpops, lociOnCols = TRUE)
    PhiHat <- popkin::popkin(Xkinship, lociOnCols = TRUE)
  } else if (percent_overlap == "0") {
    total_snps_to_simulate <- p_test + p_kinship
    # this contains all SNPs (X_{Testing}:X_{kinship})
    out <- bnpsd::rbnpsd(Q, FF, total_snps_to_simulate)
    Xall <- t(out$X) # genotypes are columns, rows are subjects
    cnames <- paste0("X", 1:total_snps_to_simulate)
    colnames(Xall) <- cnames
    rownames(Xall) <- paste0("id", 1:n)
    subpops <- ceiling((1:n) / n * k)
    table(subpops) # got k=10 subpops with 100 individuals each

    # Snps used for kinship
    snps_kinships <- sample(cnames, p_kinship, replace = FALSE)
    length(snps_kinships)

    snps_test <- setdiff(cnames, snps_kinships)
    # length(snps_test)
    # setdiff(cnames, snps_kinships) %>% length()
    causal <- sample(snps_test, ncausal, replace = FALSE)
    not_causal <- setdiff(snps_test, causal)

    Xkinship <- Xall[, snps_kinships]
    Xtest <- Xall[, snps_test]

    # now estimate kinship using popkin
    # PhiHat <- popkin::popkin(X, subpops, lociOnCols = TRUE)
    PhiHat <- popkin::popkin(Xkinship, lociOnCols = TRUE)
  }

  kin <- 2 * PhiHat
  eiK <- eigen(kin)
  if (any(eiK$values < 1e-5)) {
    eiK$values[ eiK$values < 1e-5 ] <- 1e-5
  }
  PC <- sweep(eiK$vectors, 2, sqrt(eiK$values), "*")
  # plot(eiK$values)
  # plot(PC[,1],PC[,2], pch = 19, col = rep(RColorBrewer::brewer.pal(5,"Paired"), each = 200))

  np <- dim(Xtest)
  n <- np[[1]]
  p <- np[[2]]

  x_lasso <- cbind(Xtest, PC[, 1:10])

  beta <- rep(0, length = p)
  beta[which(colnames(Xtest) %in% causal)] <- stats::rnorm(n = length(causal))
  causal_positive <- colnames(Xtest)[which(beta > 0)]
  causal_negative <- colnames(Xtest)[which(beta < 0)]

  mu <- as.numeric(Xtest %*% beta)

  P <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = eta * sigma2 * kin)
  E <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))
  y <- b0 + mu + P + E

  return(list(
    y = y, x = Xtest, causal = causal, beta = beta, kin = kin,
    Xkinship = Xkinship,
    not_causal = not_causal, causal_positive = causal_positive,
    causal_negative = causal_negative,
    x_lasso = x_lasso
  ))
}


l2norm <- function(x) sqrt(sum(x^2))


"%ni%" <- Negate("%in%")



