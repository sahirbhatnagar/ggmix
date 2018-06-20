gen_structured_model <- function(n, p_test, p_kinship, k, s, Fst, b0, beta_mean,
                                 eta, sigma2, geography = c("ind", "1d","circ"),
                                 percent_causal, percent_overlap) {

  # p_test: number of variables in X_test, i.e., the design matrix
  # p_kinship: number of variable in X_kinship, i.e., matrix used to calculate kinship
  # k:	Number of intermediate subpopulations
  # s: The desired bias coefficient, which specifies σ indirectly. Required if sigma is missing
  # F: The length-k vector of inbreeding coefficients (or FST's) of the intermediate subpopulations,
  # up to a scaling factor (which cancels out in calculations). Required if sigma is missing
  # Fst: The desired final FST of the admixed individuals. Required if sigma is missing
  # browser()
  # define population structure

  FF <- 1:k # subpopulation FST vector, up to a scalar
  # s <- 0.5 # desired bias coefficient
  # Fst <- 0.1 # desired FST for the admixed individuals
  geography <- match.arg(geography)
  if (geography == "1d") {
    obj <- bnpsd::q1d(n = n, k = k, s = s, F = FF, Fst = Fst)
    Q <- obj$Q
    FF <- obj$F
  } else if (geography == "ind") {
    ngroup <- n / k # equal sized groups
    # here’s the labels (for simplicity, list all individuals of S1 first, then S2, then S3)
    labs <- rep(paste0("S",1:k), each = ngroup)
    # data dimensions infered from labs:
    length(labs) # number of individuals "n"
    # desired admixture matrix ("is" stands for "Independent Subpopulations")
    Q <- bnpsd::qis(labs)
    FF <- 1:k # subpopulation FST vector, unnormalized so far
    FF <- FF/popkin::fst(FF)*Fst # normalized to have the desired Fst
  } else if (geography == "circ") {
    obj <- bnpsd::q1dc(n = n, k = k, s = s, F = FF, Fst = Fst)
    Q <- obj$Q
    FF <- obj$F
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
    Xall[1:5,1:5]
    dim(Xall)
    subpops <- ceiling( (1:n)/n*k )
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

    Xkinship <- Xall[,snps_kinships]
    Xtest <- Xall[,snps_test]

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
    subpops <- ceiling( (1:n)/n*k )
    table(subpops) # got k=10 subpops with 100 individuals each

    # Snps used for kinship
    snps_kinships <- sample(cnames, p_kinship, replace = FALSE)
    length(snps_kinships)

    snps_test <- setdiff(cnames, snps_kinships)
    # length(snps_test)
    # setdiff(cnames, snps_kinships) %>% length()
    causal <- sample(snps_test, ncausal, replace = FALSE)
    not_causal <- setdiff(snps_test, causal)

    Xkinship <- Xall[,snps_kinships]
    Xtest <- Xall[,snps_test]

    # now estimate kinship using popkin
    # PhiHat <- popkin::popkin(X, subpops, lociOnCols = TRUE)
    PhiHat <- popkin::popkin(Xkinship, lociOnCols = TRUE)
  }

  kin <- 2 * PhiHat
  eiK <- eigen(kin)
  if (any(eiK$values < 1e-5)) { eiK$values[ eiK$values < 1e-5 ] <- 1e-5 }
  PC <- sweep(eiK$vectors, 2, sqrt(eiK$values), "*")
  # plot(eiK$values)
  # plot(PC[,1],PC[,2], pch = 19, col = rep(RColorBrewer::brewer.pal(5,"Paired"), each = 200))

  np <- dim(Xtest)
  n <- np[[1]]
  p <- np[[2]]

  x_lasso <- cbind(Xtest,PC[,1:10])

  beta <- rep(0, length = p)
  beta[which(colnames(Xtest) %in% causal)] <- rnorm(n = length(causal))
  mu <- as.numeric(Xtest %*% beta)

  P <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = eta * sigma2 * kin)
  E <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))
  y <- b0 + mu + P + E

  return(list(y = y, x = Xtest, causal = causal, beta = beta, kin = kin,
              not_causal = not_causal,
              x_lasso = x_lasso))
}
