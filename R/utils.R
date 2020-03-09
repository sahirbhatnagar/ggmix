#' Simulation Scenario from Bhatnagar et al. (2018+) ggmix paper
#'
#' @description Function that generates data of the different simulation studies
#'   presented in the accompanying paper. This function requires the
#'   \code{popkin} and \code{bnpsd} package to be installed.
#' @param n number of observations to simulate
#' @param p_design number of variables in X_test, i.e., the design matrix
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
#' @param percent_causal percentage of \code{p_design} that is causal. must be
#'   \eqn{0 \leq percent_causal \leq 1}. The true regression coefficients are
#'   generated from a standard normal distribution.
#' @param percent_overlap this represents the percentage of causal SNPs that
#'   will also be included in the calculation of the kinship matrix
#' @param train_tune_test the proportion of sample size used for training tuning
#'   parameter selection and testing. default is 60/20/20 split
#' @details The kinship is estimated using the \code{popkin} function from the
#'   \code{popkin} package. This function will multiple that kinship matrix by 2
#'   to give the expected covariance matrix which is subsequently used in the
#'   linear mixed models
#' @return A list with the following elements \describe{\item{ytrain}{simulated
#'   response vector for training set} \item{ytune}{simulated response vector
#'   for tuning parameter selection set} \item{ytest}{simulated response vector
#'   for test set} \item{xtrain}{simulated design matrix for training
#'   set}\item{xtune}{simulated design matrix for tuning parameter selection
#'   set}\item{xtest}{simulated design matrix for testing set}
#'   \item{xtrain_lasso}{simulated design matrix for training set for lasso
#'   model. This is the same as xtrain, but also includes the nPC principal
#'   components} \item{xtune_lasso}{simulated design matrix for tuning parameter
#'   selection set for lasso model. This is the same as xtune, but also includes
#'   the nPC principal components}\item{xtest}{simulated design matrix for
#'   testing set for lasso model. This is the same as xtest, but also includes
#'   the nPC principal components} \item{causal}{character vector of the names
#'   of the causal SNPs} \item{beta}{the vector of true regression coefficients}
#'   \item{kin_train}{2 times the estimated kinship for the training set
#'   individuals} \item{kin_tune_train}{The covariance matrix between the tuning
#'   set and the training set individuals} \item{kin_test_train}{The covariance
#'   matrix between the test set and training set individuals}
#'   \item{Xkinship}{the matrix of SNPs used to estimate the kinship matrix}
#'   \item{not_causal}{character vector of the non-causal SNPs} \item{PC}{the
#'   principal components for population structure adjustment} }
#' @seealso \code{\link[bnpsd]{admix_prop_1d_linear}}
#' @export
#' @examples
#' admixed <- gen_structured_model(n = 100,
#'                                 p_design = 50,
#'                                 p_kinship = 5e2,
#'                                 geography = "1d",
#'                                 percent_causal = 0.10,
#'                                 percent_overlap = "100",
#'                                 k = 5, s = 0.5, Fst = 0.1,
#'                                 b0 = 0, nPC = 10,
#'                                 eta = 0.1, sigma2 = 1,
#'                                 train_tune_test = c(0.8, 0.1, 0.1))
#' names(admixed)
gen_structured_model <- function(n, p_design, p_kinship, k, s, Fst, b0, nPC = 10,
                                 eta, sigma2, geography = c("ind", "1d", "circ"),
                                 percent_causal, percent_overlap, train_tune_test = c(0.6, 0.2, 0.2)) {

  # p_design:
  # p_kinship:
  # k:	N
  # s:
  # F: The length-k vector of inbreeding coefficients (or FST's) of the intermediate subpopulations,
  # up to a scaling factor (which cancels out in calculations). Required if sigma is missing
  # Fst: T
  # browser()
  # define population structure

  if(sum(train_tune_test) != 1) stop("Training/tune/test split must be equal to 1")

  if (!requireNamespace("bnpsd", quietly = TRUE)) {
    stop(strwrap("Package \"bnpsd\" needed to simulate data.
                 Please install it."),
      call. = FALSE
    )
  }

  # if (!requireNamespace("simtrait", quietly = TRUE)) {
  #   stop(strwrap("Package \"simtrait\" needed to simulate data.
  #                Please install it."),
  #        call. = FALSE
  #   )
  # }

  if (!requireNamespace("popkin", quietly = TRUE)) {
    stop(strwrap("Package \"popkin\" needed to simulate data.
                 Please install it."),
      call. = FALSE
    )
  }


  ########################
  ### Generate Kinship ###
  ########################

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
    coancestry <- bnpsd::coanc_admix(admix_proportions, inbr_subpops)
    kinship <- bnpsd::coanc_to_kinship(coancestry)

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
    admix_proportions <- bnpsd::admix_prop_indep_subpops(labs)

    # subpopulation FST vector, unnormalized so far
    inbr_subpops <- 1 : k_subpops
    # normalized to have the desired FST
    # NOTE fst is a function in the `popkin` package
    inbr_subpops <- inbr_subpops / popkin::fst(inbr_subpops) * Fst
    # verify FST for the intermediate subpopulations
    # fst(inbr_subpops)
    #> [1] 0.2

    # get coancestry of the admixed individuals
    coancestry <- bnpsd::coanc_admix(admix_proportions, inbr_subpops)
    # before getting FST for individuals, weigh then inversely proportional to subpop sizes
    weights <- popkin::weights_subpops(labs) # function from `popkin` package

    kinship <- bnpsd::coanc_to_kinship(coancestry)

  } else if (geography == "circ") {
    FF <- 1:k # subpopulation FST vector, up to a scalar
    # obj <- bnpsd::admix_prop_1d_circular(n_ind = n, k_subpops = k, s = s, F = FF, Fst = Fst)
    # Q <- obj$Q
    # FF <- obj$F

    # admixture proportions from *circular* 1D geography
    obj <- bnpsd::admix_prop_1d_circular(
      n_ind = n,
      k_subpops = k,
      bias_coeff = s,
      coanc_subpops = FF,
      fst = Fst
    )
    admix_proportions <- obj$admix_proportions
    inbr_subpops <- obj$coanc_subpops

    # get pop structure parameters of the admixed individuals
    coancestry <- bnpsd::coanc_admix(admix_proportions, inbr_subpops)

    kinship <- bnpsd::coanc_to_kinship(coancestry)
  }


  #####################
  ### GENERATE SNPS ###
  #####################

  ncausal <- p_design * percent_causal
  # browser()
  if (percent_overlap == "100") {

    total_snps_to_simulate <- p_design + p_kinship - ncausal

    # this contains all SNPs (X_{Design}:X_{kinship})
    # out <- bnpsd::rbnpsd(Q = Q, F = FF, m = total_snps_to_simulate)

    # draw all random Allele Freqs (AFs) and genotypes
    # reuse the previous inbr_subpops, admix_proportions
    out <- bnpsd::draw_all_admix(
      admix_proportions = admix_proportions,
      inbr_subpops = inbr_subpops,
      m_loci = total_snps_to_simulate,
      # NOTE by default p_subpops and p_ind are not returned, but here we will ask for them
      want_p_subpops = TRUE,
      want_p_ind = TRUE
    )

    Xall <- t(out$X) # genotypes are columns, rows are subjects
    cnames <- paste0("X", 1:total_snps_to_simulate)
    colnames(Xall) <- cnames
    rownames(Xall) <- paste0("id", 1:n)
    Xall[1:5,1:5]
    dim(Xall)
    subpops <- ceiling( (1:n)/n*k )
    table(subpops) # got k=10 subpops with 100 individuals each

    ###################
    ### CAUSAL LOCI ###
    ###################
    # browser()
    # Snps used for kinship
    snps_kinships <- sample(cnames, p_kinship, replace = FALSE)

    # all causal snps are in kinship matrix
    if (percent_causal != 0 ) {
      # browser()
      # compute marginal allele frequencies
      p_anc_hat <- colMeans(Xall[,snps_kinships,drop=FALSE], na.rm = TRUE)/2

      # select random SNPs! this performs the magic...
      # also runs additional checks
      causal_indexes <- select_loci(maf = p_anc_hat, m_causal = ncausal, maf_cut = 0.05)

      # draw random SNP coefficients for selected loci
      # causal_coeffs <- stats::rnorm(ncausal, 0, 1)

      causal <- snps_kinships[causal_indexes]
      snps_design <- c(setdiff(cnames, snps_kinships), causal)
      not_causal <- setdiff(snps_design, causal)

      # subset data to consider causal loci only
      p_anc_hat <- p_anc_hat[causal_indexes]

    } else if (percent_causal == 0) {
      causal <- ""
      snps_design <- setdiff(cnames, snps_kinships)
      not_causal <- snps_design
    }

    Xkinship <- Xall[,snps_kinships]
    Xdesign <- Xall[,snps_design]

    # Xdesign_causal <- Xall[,causal_indexes,drop=F] # the subset of causal data (keep as a matrix even if m_causal == 1)

    # now estimate kinship using popkin
    PhiHat <- popkin::popkin(Xkinship, subpops = subpops, loci_on_cols = TRUE)
    mean_kinship <- mean(PhiHat)
    # PhiHat <- popkin::popkin(Xkinship, lociOnCols = TRUE)

  } else if (percent_overlap == "0") {

    total_snps_to_simulate <- p_design + p_kinship

    # this contains all SNPs (X_{Testing}:X_{kinship})
    # out <- bnpsd::rbnpsd(Q = Q, F = FF, m = total_snps_to_simulate)

    # draw all random Allele Freqs (AFs) and genotypes
    # reuse the previous inbr_subpops, admix_proportions
    out <- bnpsd::draw_all_admix(
      admix_proportions = admix_proportions,
      inbr_subpops = inbr_subpops,
      m_loci = total_snps_to_simulate,
      # NOTE by default p_subpops and p_ind are not returned, but here we will ask for them
      want_p_subpops = TRUE,
      want_p_ind = TRUE
    )


    Xall <- t(out$X) # genotypes are columns, rows are subjects
    cnames <- paste0("X", 1:total_snps_to_simulate)
    colnames(Xall) <- cnames
    rownames(Xall) <- paste0("id", 1:n)
    Xall[1:5,1:5]
    dim(Xall)
    subpops <- ceiling( (1:n)/n*k )
    table(subpops) # got k=10 subpops with 100 individuals each


    # Snps used for kinship
    snps_kinships <- sample(cnames, p_kinship, replace = FALSE)
    length(snps_kinships)

    snps_design <- setdiff(cnames, snps_kinships)
    # length(snps_design)
    # setdiff(cnames, snps_kinships) %>% length()
    if (percent_causal !=0) {

      # compute marginal allele frequencies
      p_anc_hat <- colMeans(Xall[,snps_design, drop=FALSE], na.rm = TRUE)/2

      # select random SNPs! this performs the magic...
      # also runs additional checks
      causal_indexes <- select_loci(maf = p_anc_hat, m_causal = ncausal, maf_cut = 0.05)

      # draw random SNP coefficients for selected loci
      # causal_coeffs <- stats::rnorm(ncausal, 0, 1)

      causal <- snps_design[causal_indexes]
      # causal <- sample(snps_design, ncausal, replace = FALSE)
    } else if (percent_causal == 0) {
      causal <- ""
    }

    not_causal <- setdiff(snps_design, causal)

    Xkinship <- Xall[,snps_kinships]
    Xdesign <- Xall[,snps_design]

    # now estimate kinship using popkin
    PhiHat <- popkin::popkin(Xkinship, subpops = subpops, loci_on_cols = TRUE)
    # PhiHat <- popkin::popkin(Xkinship, lociOnCols = TRUE)

  }

  np <- dim(Xdesign)
  n <- np[[1]]
  p <- np[[2]]

  beta <- rep(0, length = p)
  if (percent_causal != 0) {
    beta[which(colnames(Xdesign) %in% causal)] <- stats::rnorm(n = length(causal))
  }

  mu <- as.numeric(Xdesign %*% beta)

  # y <- trait
  kin <- 2 * PhiHat

  tt <- eta * sigma2 * kin
  if (!all(eigen(tt)$values > 0)) {
    message("eta * sigma2 * kin not PD, using Matrix::nearPD")
    tt <- Matrix::nearPD(tt)$mat
  }

  P <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = tt)
  E <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))
  # y <- mu + sigma * matrix(rnorm(nsim * n), n, nsim)
  # y <- b0 + mu + t(P) + t(E)
  # y <- MASS::mvrnorm(1, mu = mu, Sigma = eta * sigma2 * kin + (1 - eta) * sigma2 * diag(n))
  y <- b0 + mu + P + E

  # partition the data into train/tune/test
  spec = c(train = train_tune_test[1], tune = train_tune_test[2], test = train_tune_test[3])

  g = sample(cut(
    seq(nrow(Xdesign)),
    nrow(Xdesign)*cumsum(c(0,spec)),
    labels = names(spec)
  ))
  # g %>% table

  train_ind <- which(g == "train")
  tune_ind <- which(g == "tune")
  test_ind <- which(g == "test")
  # res = split(admixed$x, g)

  xtrain <- Xdesign[train_ind,,drop=FALSE]
  xtune <- Xdesign[tune_ind,,drop=FALSE]
  xtest <- Xdesign[test_ind,,drop=FALSE]

  ytrain <- y[train_ind]
  ytune <- y[tune_ind]
  ytest <- y[test_ind]

  kin_train <- kin[train_ind,train_ind]
  kin_tune_train <- kin[tune_ind,train_ind]
  kin_test_train <- kin[test_ind,train_ind]

  # xtrain <- Xdesign[ind,,drop=FALSE]
  # xtune <- Xdesign[-ind,,drop=FALSE]
  PC_all <- stats::prcomp(Xdesign)
  PC <- stats::prcomp(xtrain)
  xtrain_lasso <- cbind(xtrain, PC$x[,1:nPC])
  xtune_pc <- stats::predict(PC, newdata = xtune)
  xtune_lasso <- cbind(xtune, xtune_pc[,1:nPC])
  xtest_pc <- stats::predict(PC, newdata = xtest)
  xtest_lasso <- cbind(xtest, xtest_pc[,1:nPC])

  mu_train <- mu[train_ind]


  return(list(ytrain = ytrain,
              ytune = ytune,
              ytest = ytest,

              xtrain = xtrain,
              xtune = xtune,
              xtest = xtest,

              xtrain_lasso = xtrain_lasso,
              xtune_lasso = xtune_lasso,
              xtest_lasso = xtest_lasso,

              Xkinship = Xall[train_ind,snps_kinships, drop = F],
              kin_train = kin_train,
              kin_tune_train = kin_tune_train,
              kin_test_train = kin_test_train, # covaraince between train and test

              mu_train = mu_train, # Xbeta for training set

              causal = causal,
              beta = beta,

              not_causal = not_causal,

              # used in manuscript to generate simulation study figures
              kinship = kinship,
              coancestry = coancestry,
              PC = PC_all$x[,1:nPC],
              subpops = subpops
  ))
}


l2norm <- function(x) sqrt(sum(x^2))


"%ni%" <- Negate("%in%")

# internal function
# taken verbatim from glmnet package.
# it used to be exported by glmnet, but is no longer exported.
lambda.interp <- function (lambda, s) {
  if (length(lambda) == 1) {
    nums = length(s)
    left = rep(1, nums)
    right = left
    sfrac = rep(1, nums)
  }
  else {
    s[s > max(lambda)] = max(lambda)
    s[s < min(lambda)] = min(lambda)
    k = length(lambda)
    sfrac <- (lambda[1] - s)/(lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda)/(lambda[1] - lambda[k])
    coord <- stats::approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac = (sfrac - lambda[right])/(lambda[left] - lambda[right])
    sfrac[left == right] = 1
    sfrac[abs(lambda[left] - lambda[right]) < .Machine$double.eps] = 1
  }
  list(left = left, right = right, frac = sfrac)
}



# internal function
# taken verbatim from glmnet package.
# it used to be exported by glmnet, but is no longer exported.
nonzeroCoef <- function (beta, bystep = FALSE) {
  nr = nrow(beta)
  if (nr == 1) {
    if (bystep)
      apply(beta, 2, function(x) if (abs(x) > 0)
        1
        else NULL)
    else {
      if (any(abs(beta) > 0))
        1
      else NULL
    }
  }
  else {
    beta = abs(beta) > 0
    which = seq(nr)
    ones = rep(1, ncol(beta))
    nz = as.vector((beta %*% ones) > 0)
    which = which[nz]
    if (bystep) {
      if (length(which) > 0) {
        beta = as.matrix(beta[which, , drop = FALSE])
        nzel = function(x, which) if (any(x))
          which[x]
        else NULL
        which = apply(beta, 2, nzel, which)
        if (!is.list(which))
          which = data.frame(which)
        which
      }
      else {
        dn = dimnames(beta)[[2]]
        which = vector("list", length(dn))
        names(which) = dn
        which
      }
    }
    else which
  }
}



# internal function
# taken verbatim from https://github.com/OchoaLab/simtrait/blob/master/R/select_loci.R
select_loci <- function(maf, m_causal, maf_cut = 0.05) {
  # check for missing parameters
  if (missing(maf))
    stop('marginal allele frequency vector `maf` is required!')
  if (missing(m_causal))
    stop('the number of causal loci `m_causal` is required!')

  # data dimensions
  m <- length(maf)
  # other checks
  if (m_causal > m)
    stop('the number of causal loci cannot be larger than the total number of loci (', m_causal, ' > ', m, ')')

  # select random loci!
  # we might not want to pick extremely rare alleles, so set MAF thresholds
  i <- which(maf_cut <= maf & maf <= 1 - maf_cut) # candidate locus indexes
  sample(i, m_causal) # these are the chosen locus indeces!
}
