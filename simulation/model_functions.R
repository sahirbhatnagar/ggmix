## @knitr models


# USE THIS ONE MAY 7, 2019
make_ADmixed_model <- function(n, p_design, p_kinship, k, s, Fst, b0, beta_mean,
                               eta, sigma2, geography = c("ind", "1d","circ"),
                               percent_causal, percent_overlap) {

  # p_design: number of variables in X_design, i.e., the design matrix
  # p_kinship: number of variable in X_kinship, i.e., matrix used to calculate kinship
  # k:	Number of intermediate subpopulations
  # s: The desired bias coefficient, which specifies σ indirectly. Required if sigma is missing
  # F: The length-k vector of inbreeding coefficients (or FST's) of the intermediate subpopulations,
  # up to a scaling factor (which cancels out in calculations). Required if sigma is missing
  # Fst: The desired final FST of the admixed individuals. Required if sigma is missing
  # browser()
  # define population structure


  # s <- 0.5 # desired bias coefficient
  # Fst <- 0.1 # desired FST for the admixed individuals
  geography <- match.arg(geography)


  new_model(name = "ggmix_05_07_2019",
            label = sprintf("percent_causal = %s, percent_overlap = %s, eta = %s,
                            sigma2 = %s, geography = %s, p_design = %s, p_kinship = %s, beta_mean = %s",
                            percent_causal,
                            percent_overlap,
                            eta,
                            sigma2,
                            geography,
                            p_design,
                            p_kinship,
                            beta_mean),
            params = list(n = n,
                          p_design = p_design,
                          p_kinship = p_kinship,
                          k = k,
                          s = s,
                          Fst = Fst,
                          b0 = b0,
                          beta_mean = beta_mean,
                          eta = eta,
                          sigma2 = sigma2,
                          geography = geography,
                          percent_causal = percent_causal,
                          percent_overlap = percent_overlap),
            simulate = function(n,
                                p_design,
                                p_kinship,
                                k,
                                s,
                                Fst,
                                b0,
                                beta_mean,
                                eta,
                                sigma2,
                                geography,
                                percent_causal,
                                percent_overlap,
                                nsim) {

              models <- list()

              for (i in seq(nsim)) {

                ########################
                ### Generate Kinship ###
                ########################

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

                  n1 <- n2 <- n3 <- n4 <- n5 <- 200*2

                  # here’s the labels (for simplicity, list all individuals of S1 first, then S2, then S3)
                  labs <- c( rep.int("S1", n1), rep.int("S2", n2), rep.int("S3", n3),
                             rep.int("S4", n4), rep.int("S5", n5))
                  # data dimensions infered from labs:
                  length(labs) # number of individuals "n"

                  # train
                  # desired admixture matrix ("is" stands for "Independent Subpopulations")
                  # Q <- bnpsd::admix_prop_indep_subpops(labs)
                  # FF <- 1:k # subpopulation FST vector, unnormalized so far
                  # FF <- FF/popkin::fst(FF)*Fst # normalized to have the desired Fst


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

                #####################
                ### GENERATE SNPS ###
                #####################


                ncausal <- p_design * percent_causal

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
                    causal_indexes <- simtrait:::select_loci(maf = p_anc_hat, m_causal = ncausal, maf_cut = 0.05)

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
                  out <- draw_all_admix(
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
                    causal_indexes <- simtrait:::select_loci(maf = p_anc_hat, m_causal = ncausal, maf_cut = 0.05)

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
                  # beta[which(colnames(Xdesign) %in% causal)] <- runif(n = length(causal), beta_mean - 0.3, beta_mean + 0.3)
                  # beta[which(colnames(Xdesign) %in% causal)] <- c(2,4,3,3,1)
                  beta[which(colnames(Xdesign) %in% causal)] <- rnorm(n = length(causal))
                }

                #############
                ### SCALE ###
                #############

                # causal_coeffs <- beta[which(colnames(Xdesign) %in% causal)]
                #
                # # to scale causal_coeffs to give correct heritability, we need to estimate the pq = p(1-p) vector
                # # calculate pq = p_anc * (1 - p_anc) in one of two ways
                # # indirect, infer from genotypes and mean kinship
                # # recall E[ p_anc_hat * (1 - p_anc_hat) ] = pq * (1 - mean_kinship), so we solve for pq:
                # pq <- p_anc_hat * (1 - p_anc_hat) / (1 - mean_kinship)
                #
                # # the initial genetic variance is
                # sigma_0 <- sqrt( 2 * sum( pq * causal_coeffs^2 ) )
                #
                # sigma_sq = sigma2
                # herit = eta
                # # adjust causal_coeffss so final variance is sigma_sq*herit as desired!
                # # causal_coeffs <- causal_coeffs * sqrt( sigma_sq * herit ) / sigma_0 # scale by standard deviations
                #
                # G <- drop( Xdesign_causal %*% causal_coeffs ) # this is a vector
                # # NOTE by construction:
                # # Cov(G) = 2 * herit * kinship
                #
                # # works very well assuming causal_coeffs and p_anc are uncorrelated!
                # muXB <- 2 * sum( causal_coeffs ) * mean( p_anc_hat )
                #
                # # in all cases:
                # # - remove the mean from the genotypes (muXB)
                # # - add the desired mean
                # mu = 0
                # # G <- G - muXB + mu
                #
                # # draw noise
                # E <- stats::rnorm(n, 0, (1 - herit) * sigma_sq ) # noise has mean zero but variance ((1-herit) * sigma_sq)
                # # NOTE by construction:
                # # Cov(E) = (1-herit) * sigma_sq * I
                #
                # # lastly, here's the trait:
                # trait <- G + E

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

                # ind <- caret::createDataPartition(y, p = 0.5, list = FALSE)[,1]

                # partition the data into train/test/validate
                # here i used train/test/validate, but I should have used train/tune/test like I did in
                # the gen_structured_model function in utils.R
                spec = c(train = .6, test = .2, validate = .2)

                g = sample(cut(
                  seq(nrow(Xdesign)),
                  nrow(Xdesign)*cumsum(c(0,spec)),
                  labels = names(spec)
                ))
                # g %>% table

                train_ind <- which(g == "train")
                validate_ind <- which(g == "validate")
                test_ind <- which(g == "test")
                # res = split(admixed$x, g)

                xtrain <- Xdesign[train_ind,,drop=FALSE]
                xtest <- Xdesign[test_ind,,drop=FALSE]
                xvalidate <- Xdesign[validate_ind,,drop=FALSE]

                ytrain <- y[train_ind]
                ytest <- y[test_ind]
                yvalidate <- y[validate_ind]

                kin_train <- kin[train_ind,train_ind]
                kin_test_train <- kin[test_ind,train_ind]
                kin_validate_train <- kin[validate_ind,train_ind]

                # xtrain <- Xdesign[ind,,drop=FALSE]
                # xtest <- Xdesign[-ind,,drop=FALSE]

                PC <- prcomp(xtrain)
                xtrain_lasso <- cbind(xtrain, PC$x[,1:10])
                xtest_pc <- predict(PC, newdata = xtest)
                xtest_lasso <- cbind(xtest, xtest_pc[,1:10])
                xvalidate_pc <- predict(PC, newdata = xvalidate)
                xvalidate_lasso <- cbind(xvalidate, xvalidate_pc[,1:10])

                # kin_train <- kin[ind,ind]
                # kin_test_train <- kin[-ind,ind]

                # ytrain <- y[ind]
                # ytest <- y[-ind]

                mu_train <- mu[train_ind]

                models[[i]] <- list(ytrain = ytrain,
                                    ytest = ytest,
                                    yvalidate = yvalidate,

                                    xtrain = xtrain,
                                    xtest = xtest,
                                    xvalidate = xvalidate,

                                    xtrain_lasso = xtrain_lasso,
                                    xtest_lasso = xtest_lasso,
                                    xvalidate_lasso = xvalidate_lasso,

                                    kin_train = kin_train,
                                    kin_test_train = kin_test_train, # covaraince between train and test
                                    kin_validate_train = kin_validate_train,

                                    mu_train = mu_train, # Xbeta for training set

                                    causal = causal,
                                    beta = beta,

                                    not_causal = not_causal
                                    )
              }
              return(models)
            })

}



## @knitr models-not-used

make_ADmixed_model_simtrait <- function(n, p_design, p_kinship, k, s, Fst, b0, beta_mean,
                               eta, sigma2, geography = c("ind", "1d","circ"),
                               percent_causal, percent_overlap) {

  # p_design: number of variables in X_design, i.e., the design matrix
  # p_kinship: number of variable in X_kinship, i.e., matrix used to calculate kinship
  # k:	Number of intermediate subpopulations
  # s: The desired bias coefficient, which specifies σ indirectly. Required if sigma is missing
  # F: The length-k vector of inbreeding coefficients (or FST's) of the intermediate subpopulations,
  # up to a scaling factor (which cancels out in calculations). Required if sigma is missing
  # Fst: The desired final FST of the admixed individuals. Required if sigma is missing
  # browser()
  # define population structure


  # s <- 0.5 # desired bias coefficient
  # Fst <- 0.1 # desired FST for the admixed individuals
  geography <- match.arg(geography)


  new_model(name = "ggmix_05_07_2019",
            label = sprintf("percent_causal = %s, percent_overlap = %s, eta = %s,
                            sigma2 = %s, geography = %s, p_design = %s, p_kinship = %s, beta_mean = %s",
                            percent_causal,
                            percent_overlap,
                            eta,
                            sigma2,
                            geography,
                            p_design,
                            p_kinship,
                            beta_mean),
            params = list(n = n,
                          p_design = p_design,
                          p_kinship = p_kinship,
                          k = k,
                          s = s,
                          Fst = Fst,
                          b0 = b0,
                          beta_mean = beta_mean,
                          eta = eta,
                          sigma2 = sigma2,
                          geography = geography,
                          percent_causal = percent_causal,
                          percent_overlap = percent_overlap),
            simulate = function(n,
                                p_design,
                                p_kinship,
                                k,
                                s,
                                Fst,
                                b0,
                                beta_mean,
                                eta,
                                sigma2,
                                geography,
                                percent_causal,
                                percent_overlap,
                                nsim) {

              models <- list()

              for (i in seq(nsim)) {

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

                  n1 <- n2 <- n3 <- n4 <- n5 <- 200*2

                  # here’s the labels (for simplicity, list all individuals of S1 first, then S2, then S3)
                  labs <- c( rep.int("S1", n1), rep.int("S2", n2), rep.int("S3", n3),
                             rep.int("S4", n4), rep.int("S5", n5))
                  # data dimensions infered from labs:
                  length(labs) # number of individuals "n"

                  # train
                  # desired admixture matrix ("is" stands for "Independent Subpopulations")
                  # Q <- bnpsd::admix_prop_indep_subpops(labs)
                  # FF <- 1:k # subpopulation FST vector, unnormalized so far
                  # FF <- FF/popkin::fst(FF)*Fst # normalized to have the desired Fst


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

                total_snps_to_simulate <- p_design

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
                subpops <- ceiling( (1:n)/n*k )

                # ancestral allele frequencies
                p_anc <- out$p_anc

                np <- dim(Xall)
                n <- np[[1]]
                p <- np[[2]]

                # parameters of simulation
                m_causal <- p_design * percent_causal
                herit <- eta

                # create simulated trait and associated data

                # version 1: known p_anc (prefered, only applicable to simulated data)
                obj <- sim_trait(X = out$X, m_causal = m_causal, herit = herit, p_anc = p_anc)
                # version 2: known kinship (more broadly applicable but fewer guarantees)
                # obj <- sim_trait(X = X, m_causal = m_causal, herit = herit, kinship = kinship)

                # outputs in both versions:
                # trait vector
                y <- obj$y
                # randomly-picked causal locus index
                causal <- cnames[obj$i]
                not_causal <- setdiff(cnames, causal)


                # locus effect size vector
                beta <- obj$beta

                # theoretical covariance of the simulated traits
                # kin <- cov_trait(kinship = kinship, herit = herit)
                kin <- kinship
                # partition the data into train/test/validate
                spec = c(train = .6, test = .2, validate = .2)

                g = sample(cut(
                  seq(nrow(Xall)),
                  nrow(Xall)*cumsum(c(0,spec)),
                  labels = names(spec)
                ))
                # g %>% table

                Xdesign <- Xall

                train_ind <- which(g == "train")
                validate_ind <- which(g == "validate")
                test_ind <- which(g == "test")
                # res = split(admixed$x, g)

                xtrain <- Xdesign[train_ind,,drop=FALSE]
                xtest <- Xdesign[test_ind,,drop=FALSE]
                xvalidate <- Xdesign[validate_ind,,drop=FALSE]

                ytrain <- y[train_ind]
                ytest <- y[test_ind]
                yvalidate <- y[validate_ind]

                kin_train <- kin[train_ind,train_ind]
                kin_test_train <- kin[test_ind,train_ind]
                kin_validate_train <- kin[validate_ind,train_ind]

                # xtrain <- Xdesign[ind,,drop=FALSE]
                # xtest <- Xdesign[-ind,,drop=FALSE]

                PC <- prcomp(xtrain)
                xtrain_lasso <- cbind(xtrain, PC$x[,1:10])
                xtest_pc <- predict(PC, newdata = xtest)
                xtest_lasso <- cbind(xtest, xtest_pc[,1:10])
                xvalidate_pc <- predict(PC, newdata = xvalidate)
                xvalidate_lasso <- cbind(xvalidate, xvalidate_pc[,1:10])

                # kin_train <- kin[ind,ind]
                # kin_test_train <- kin[-ind,ind]

                # ytrain <- y[ind]
                # ytest <- y[-ind]

                mu_train <- (xtrain[,obj$i,drop = FALSE] %*% beta)[train_ind]

                models[[i]] <- list(ytrain = ytrain,
                                    ytest = ytest,
                                    yvalidate = yvalidate,

                                    xtrain = xtrain,
                                    xtest = xtest,
                                    xvalidate = xvalidate,

                                    xtrain_lasso = xtrain_lasso,
                                    xtest_lasso = xtest_lasso,
                                    xvalidate_lasso = xvalidate_lasso,

                                    kin_train = kin_train,
                                    kin_test_train = kin_test_train, # covaraince between train and test
                                    kin_validate_train = kin_validate_train,

                                    mu_train = mu_train, # Xbeta for training set

                                    causal = causal,
                                    beta = beta,

                                    not_causal = not_causal
                )
              }
              return(models)
            })

}

make_ADmixed_model_not_sim <- function(n, p_test, p_kinship, k, s, Fst, b0, beta_mean,
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
    n_train_k <- n / k
    n1 <- n_train_k; n2 <- n_train_k; n3 <- n_train_k; n4 <- n_train_k; n5 <- n_train_k

    # here’s the labels (for simplicity, list all individuals of S1 first, then S2, then S3)
    labs <- c( rep.int("S1", n1), rep.int("S2", n2), rep.int("S3", n3),
               rep.int("S4", n4), rep.int("S5", n5))
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
    if (percent_causal != 0 ) {
      causal <- sample(snps_kinships, ncausal, replace = FALSE)
      snps_design <- c(setdiff(cnames, snps_kinships), causal)
      not_causal <- setdiff(snps_design, causal)
    } else if (percent_causal == 0) {
      causal <- ""
      snps_design <- setdiff(cnames, snps_kinships)
      not_causal <- snps_design
    }

    # length(snps_design)
    # setdiff(cnames, snps_kinships) %>% length()

    # browser()
    Xkinship <- Xall[,snps_kinships]
    Xtest <- Xall[,snps_design]

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
      causal <- sample(snps_design, ncausal, replace = FALSE)
    } else if (percent_causal == 0) {
      causal <- ""
    }

    not_causal <- setdiff(snps_design, causal)

    Xkinship <- Xall[,snps_kinships]
    Xtest <- Xall[,snps_design]

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
  x_lasso[1:5,1:5]
  # kin <- snpStats::xxt(dat$genotypes)/p

  beta <- rep(0, length = p)
  if (percent_causal != 0) {
    beta[which(colnames(Xtest) %in% causal)] <- runif(n = length(causal), beta_mean - 0.3, beta_mean + 0.3)
  }
  # beta[which(colnames(Xtest) %in% causal)] <- rnorm(n = length(causal))
  mu <- as.numeric(Xtest %*% beta)

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

  models <- list(y = y, Xtest = Xtest, causal = causal,
                 beta = beta, kin = kin,
                 mu = mu,
                 Xkinship = Xkinship,
                 not_causal = not_causal,
                 x_lasso = x_lasso)

  return(models)


}


make_mixed_model_SSC <- function(b0, beta_mean, eta, sigma2, percent_causal, percent_overlap) {

  if (percent_causal == 1) {

    file_paths <- switch(percent_overlap,
                         `0` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_8k.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_8k.rel.id",
                                #X_phi is bed files used to make kinship
                                X_Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_other_8k",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_1k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_10")
                         },
                         `100` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_7990_causal_10.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_7990_causal_10.rel.id",
                                #X_phi is bed files used to make kinship
                                X_Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_other_7990_causal_10",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_1k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_10")
                         })

  } else if (percent_causal == 2.5) {

    file_paths <- switch(percent_overlap,
                         `0` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_4k.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_4k.rel.id",
                                #X_phi is bed files used to make kinship
                                X_Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_other_4k",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_100")
                         },
                         `100` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_7900_causal_100.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_7900_causal_100.rel.id",
                                #X_phi is bed files used to make kinship
                                X_Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_other_7900_causal_100",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_100")
                         })

  } else if (percent_causal == 10) {

    file_paths <- switch(percent_overlap,
                         `0` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_4k.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_4k.rel.id",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_400")
                         },
                         `50` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_3800_causal_200.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_3800_causal_200.rel.id",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_400")
                         },
                         `100` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_3600_causal_400.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_3600_causal_400.rel.id",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_400")
                         })

  } else if (percent_causal == 50) {

    file_paths <- switch(percent_overlap,
                         `0` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_4k.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_4k.rel.id",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_2000")
                         },
                         `50` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_3000_causal_1000.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_3000_causal_1000.rel.id",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_2000")
                         },
                         `100` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_2000_causal_2000.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_2000_causal_2000.rel.id",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_2000")
                         })
  }

  # gaston kinship
  Phi <- gaston::read.bed.matrix(file_paths$X_Phi)
  kin <- gaston::GRM(Phi, autosome.only = FALSE)
  kin[1:5,1:5]
  all(complete.cases(kin))
  eiK <- eigen(kin)
  # all(rownames(as.matrix(x))==rownames(kin))
  # deal with a small negative eigen value
  any(eiK$values < 0)
  eiK$values[ eiK$values < 0 ] <- 0
  PC <- sweep(eiK$vectors, 2, sqrt(eiK$values), "*")

  dat <- gaston::read.bed.matrix(file_paths$bedfile)
  gaston::standardize(dat) <- "p"
  X <- as.matrix(dat)
  any(is.na(X))
  X[is.na(X)] <- 0
  any(is.na(X))
  X[1:5,1:5]

  # need to re-order
  all(rownames(X)==rownames(kin))
  all(rownames(X) %in% rownames(kin))
  X <- X[match(rownames(kin), rownames(X)),]
  all(rownames(X)==rownames(kin))

  np <- dim(X)
  n <- np[[1]]
  p <- np[[2]]

  all(rownames(X)==rownames(kin))

  x_lasso <- cbind(X,PC[,1:10])
  x_lasso[1:5,1:5]

  causal <- data.table::fread(file_paths$causal_list, header = FALSE)$V1
  not_causal <- setdiff(colnames(X), causal)
  beta <- rep(0, length = p)
  beta[which(colnames(X) %in% causal)] <- runif(n = length(causal), beta_mean - 0.1, beta_mean + 0.1)
  mu <- as.numeric(X %*% beta)

  new_model(name = "ggmixSSCv3", label = sprintf("percent_causal = %s, percent_overlap = %s, eta = %s, sigma = %s, beta_mean = %s",
                                                 percent_causal, percent_overlap, eta, sigma2, beta_mean),
            params = list(mu = mu, n = n, x = X, x_lasso = x_lasso,
                          beta = beta, percent_causal = percent_causal,
                          percent_overlap = percent_overlap,
                          not_causal = not_causal,
                          kin = kin, b0 = b0, sigma2 = sigma2, eta = eta, causal = causal, beta_mean = beta_mean),
            simulate = function(mu, sigma2, eta, kin, n, nsim) {
              P <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = eta * sigma2 * kin)
              E <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))
              y <- b0 + mu + t(P) + t(E)
              return(split(y, col(y))) # make each col its own list element
            })
}




make_sparse_linear_model <- function(b0, eta, sigma2, type) {
  file_paths <- switch(type,
                       causal_400 = {
                         list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_3600_causal_400.rel",
                              Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_3600_causal_400.rel.id",
                              bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                              causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_400")
                       },
                       causal_2k = {
                         list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_2000_causal_2000.rel",
                              Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_2000_causal_2000.rel.id",
                              bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                              causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_2000")
                       },
                       causal_4k = {
                         list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_causal_4000.rel",
                              Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_causal_4000.rel.id",
                              bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                              causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_4000")
                       })

  kin <- as.matrix(read.table(file_paths$Phi))
  kin_names <- data.table::fread(file_paths$Phi_names)
  dimnames(kin)[[1]] <- kin_names$V1
  dimnames(kin)[[2]] <- kin_names$V1
  dat <- snpStats::read.plink(file_paths$bedfile)
  x <- as(dat$genotypes, "numeric")
  np <- dim(x)
  n <- np[[1]]
  p <- np[[2]]
  causal <- data.table::fread(file_paths$causal_list, header = FALSE)$V1
  not_causal <- setdiff(colnames(x), causal)
  beta <- rep(0, length = p)
  beta[which(colnames(x) %in% causal)] <- rnorm(n = length(causal))
  mu <- as.numeric(x %*% beta)


  # kin <- as.matrix(read.table(file_paths$Phi))[1:200,1:200]
  # kin_names <- data.table::fread(file_paths$Phi_names)
  # dimnames(kin)[[1]] <- kin_names$V1[1:200]
  # dimnames(kin)[[2]] <- kin_names$V1[1:200]
  # dat <- snpStats::read.plink(file_paths$bedfile)
  # x <- as(dat$genotypes, "numeric")[1:200,1:100]
  # np <- dim(x)
  # n <- np[[1]]
  # p <- np[[2]]
  # causal <- colnames(x)[1:10]
  # not_causal <- setdiff(colnames(x), causal)
  # beta <- rep(0, length = p)
  # beta[which(colnames(x) %in% causal)] <- rnorm(n = length(causal))
  # mu <- as.numeric(x %*% beta)


  new_model(name = "penfam", label = sprintf("type = %s, eta = %s", type, eta),
            params = list(mu = mu, n = n, x = x, beta = beta, type = type, not_causal = not_causal,
                          kin = kin, b0 = b0, sigma2 = sigma2, eta = eta, causal = causal),
            simulate = function(mu, sigma2, eta, kin, n, nsim) {
              P <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = eta * sigma2 * kin)
              E <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))
              # y <- mu + sigma * matrix(rnorm(nsim * n), n, nsim)
              y <- b0 + mu + t(P) + t(E)
              return(split(y, col(y))) # make each col its own list element
            })
}


make_mixed_model <- function(b0, eta, sigma2, type) {
  file_paths <- switch(type,
                       causal_400 = {
                         list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_3600_causal_400.rel",
                              Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_3600_causal_400.rel.id",
                              bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                              causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_400")
                       },
                       causal_2k = {
                         list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_2000_causal_2000.rel",
                              Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_2000_causal_2000.rel.id",
                              bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                              causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_2000")
                       },
                       causal_4k = {
                         list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_causal_4000.rel",
                              Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_causal_4000.rel.id",
                              bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                              causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_4000")
                       })

  kin <- as.matrix(read.table(file_paths$Phi))
  kin_names <- data.table::fread(file_paths$Phi_names)
  dimnames(kin)[[1]] <- kin_names$V1
  dimnames(kin)[[2]] <- kin_names$V1
  isPD <- all(eigen(kin)$values > 0)
  how_many_neg_eigenvalues <- sum(eigen(kin)$values <= 0)
  if (!isPD) {
    kinPD <- as(nearPD(kin)$mat,"matrix")
    dimnames(kinPD)[[1]] <- kin_names$V1
    dimnames(kinPD)[[2]] <- kin_names$V1
    kin <- kinPD
  }
  dat <- snpStats::read.plink(file_paths$bedfile)
  x <- as(dat$genotypes, "numeric")
  np <- dim(x)
  n <- np[[1]]
  p <- np[[2]]

  # kin <- snpStats::xxt(dat$genotypes)/p

  causal <- data.table::fread(file_paths$causal_list, header = FALSE)$V1
  not_causal <- setdiff(colnames(x), causal)
  beta <- rep(0, length = p)
  beta[which(colnames(x) %in% causal)] <- runif(n = length(causal), 0.1, 0.3)
  mu <- as.numeric(x %*% beta)

  new_model(name = "ggmix", label = sprintf("type = %s, eta = %s", type, eta),
            params = list(mu = mu, n = n, x = x, beta = beta, type = type, not_causal = not_causal,
                          kin = kin, b0 = b0, sigma2 = sigma2, eta = eta, causal = causal),
            simulate = function(mu, sigma2, eta, kin, n, nsim) {
              P <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = eta * sigma2 * kin)
              E <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))
              # y <- mu + sigma * matrix(rnorm(nsim * n), n, nsim)
              y <- b0 + mu + t(P) + t(E)
              return(split(y, col(y))) # make each col its own list element
            })
}





make_mixed_model_not_simulator <- function(b0,betamean,eta, sigma2,  percent_causal, percent_overlap) {

  if (percent_causal == 1) {

    file_paths <- switch(percent_overlap,
                         `0` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_8k.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_8k.rel.id",
                                #X_phi is bed files used to make kinship
                                X_Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_other_8k",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_1k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_10")
                         },
                         `100` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_7990_causal_10.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_7990_causal_10.rel.id",
                                #X_phi is bed files used to make kinship
                                X_Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_other_7990_causal_10",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_1k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_10")
                         })

  } else if (percent_causal == 2.5) {

    file_paths <- switch(percent_overlap,
                         `0` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_4k.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_4k.rel.id",
                                #X_phi is bed files used to make kinship
                                X_Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_other_4k",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_100")
                         },
                         `100` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_7900_causal_100.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_7900_causal_100.rel.id",
                                #X_phi is bed files used to make kinship
                                X_Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_other_7900_causal_100",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_100")
                         })

  } else if (percent_causal == 10) {

    file_paths <- switch(percent_overlap,
                         `0` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_4k.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_4k.rel.id",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_400")
                         },
                         `50` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_3800_causal_200.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_3800_causal_200.rel.id",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_400")
                         },
                         `100` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_3600_causal_400.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_3600_causal_400.rel.id",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_400")
                         })

  } else if (percent_causal == 50) {

    file_paths <- switch(percent_overlap,
                         `0` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_4k.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_4k.rel.id",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_2000")
                         },
                         `50` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_3000_causal_1000.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_3000_causal_1000.rel.id",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_2000")
                         },
                         `100` = {
                           list(Phi = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_2000_causal_2000.rel",
                                Phi_names = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/kinship/Kinship_other_2000_causal_2000.rel.id",
                                bedfile = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_test_4k",
                                causal_list = "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/snplist/causal_SNP_from_X_test_2000")
                         })


  }
  # browser()

  # this uses plink kinship matrix
  # kin <- as.matrix(read.table(file_paths$Phi))
  # kin_names <- data.table::fread(file_paths$Phi_names)
  # dimnames(kin)[[1]] <- kin_names$V1
  # dimnames(kin)[[2]] <- kin_names$V1
  # kin[1:5,1:5]
  # isPD <- all(eigen(kin)$values > 0)
  # how_many_neg_eigenvalues <- sum(eigen(kin)$values <= 0)
  # if (!isPD) {
  #   kinPD <- as(nearPD(kin)$mat,"matrix")
  #   dimnames(kinPD)[[1]] <- kin_names$V1
  #   dimnames(kinPD)[[2]] <- kin_names$V1
  #   kin <- kinPD
  # }

  # gaston kinship
  Phi <- gaston::read.bed.matrix(file_paths$X_Phi)
  kin <- gaston::GRM(Phi, autosome.only = FALSE)
  kin[1:5,1:5]
  all(complete.cases(kin))
  eiK <- eigen(kin)
  # all(rownames(as.matrix(x))==rownames(kin))
  # deal with a small negative eigen value
  any(eiK$values < 0)
  eiK$values[ eiK$values < 0 ] <- 0
  PC <- sweep(eiK$vectors, 2, sqrt(eiK$values), "*")

  # dat <- snpStats::read.plink(file_paths$bedfile)
  # x <- as(dat$genotypes, "numeric")
  # np <- dim(x)
  # n <- np[[1]]
  # p <- np[[2]]

  # need to re-order
  # all(rownames(x)==rownames(kin))
  # all(rownames(x) %in% rownames(kin))
  # x <- x[match(rownames(kin), rownames(x)),]
  # all(rownames(x)==rownames(kin))
  #
  #
  # # pcaCars <- princomp(kin, cor = TRUE)
  # eiK <- eigen(kin)
  # all(rownames(as.matrix(x))==rownames(kin))
  # # deal with a small negative eigen value
  # any(eiK$values < 0)
  # eiK$values[ eiK$values < 0 ] <- 0
  # PC <- sweep(eiK$vectors, 2, sqrt(eiK$values), "*")

  dat <- gaston::read.bed.matrix(file_paths$bedfile)
  gaston::standardize(dat) <- "p"
  X <- as.matrix(dat)
  any(is.na(X))
  X[is.na(X)] <- 0
  any(is.na(X))
  X[1:5,1:5]

  # need to re-order
  all(rownames(X)==rownames(kin))
  all(rownames(X) %in% rownames(kin))
  X <- X[match(rownames(kin), rownames(X)),]
  all(rownames(X)==rownames(kin))

  np <- dim(X)
  n <- np[[1]]
  p <- np[[2]]
  #
  # evv <- eigen(kin, symmetric=TRUE)
  # pcs <- evv$vectors[,1:10]
  #
  # plot(pcs[,3], pcaCars$loadings[,3])

  # # view objects stored in pcaCars
  # names(pcaCars)
  # # proportion of variance explained
  # summary(pcaCars)
  # # scree plot
  # plot(pcaCars, type = "l")
  # plot(eiK$values)
  all(rownames(X)==rownames(kin))


  x_lasso <- cbind(X,PC[,1:10])
  x_lasso[1:5,1:5]
  # kin <- snpStats::xxt(dat$genotypes)/p

  causal <- data.table::fread(file_paths$causal_list, header = FALSE)$V1
  not_causal <- setdiff(colnames(X), causal)
  beta <- rep(0, length = p)
  beta[which(colnames(X) %in% causal)] <- runif(n = length(causal), 0.9, 1.1)
  # beta[which(colnames(x) %in% causal)] <- rnorm(n = length(causal))
  mu <- as.numeric(X %*% beta)

  # sum(complete.cases(X))
  # dim(X)
  # length(beta)
  # kin <- as.matrix(read.table(file_paths$Phi))[1:500,1:500]
  # kin_names <- data.table::fread(file_paths$Phi_names)
  # dimnames(kin)[[1]] <- kin_names$V1[1:200]
  # dimnames(kin)[[2]] <- kin_names$V1[1:200]
  #
  # # browser()
  #
  # dat <- snpStats::read.plink(file_paths$bedfile)
  # x <- as(dat$genotypes, "numeric")[1:200,1:100]
  # np <- dim(x)
  # n <- np[[1]]
  # p <- np[[2]]

  # str(dat)
  # xst <- snpStats::snp.post.multiply(dat$genotypes, diag(p))
  # dimnames(xst)[[2]] <- colnames(x)
  # xst[1:10,1:10]

  # causal <- colnames(x)[1:10]
  # not_causal <- setdiff(colnames(x), causal)
  # beta <- rep(0, length = p)
  # beta[which(colnames(x) %in% causal)] <- rnorm(n = length(causal))
  # mu <- as.numeric(x %*% beta)

  P <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = eta * sigma2 * kin)
  E <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))
  # y <- mu + sigma * matrix(rnorm(nsim * n), n, nsim)
  # y <- b0 + mu + t(P) + t(E)
  y <- MASS::mvrnorm(1, mu = mu, Sigma = eta * sigma2 * kin + (1 - eta) * sigma2 * diag(n))
  # y <- 0 + mu + P + E

  return(list(y = y, x = X, causal = causal, beta = beta, kin = kin,
              x_lasso = x_lasso,
              file_paths = file_paths))

  # new_model(name = "penfam", label = sprintf("type = %s, eta = %s", type, eta),
  #           params = list(mu = mu, n = n, x = x, beta = beta, type = type, not_causal = not_causal,
  #                         kin = kin, b0 = b0, sigma2 = sigma2, eta = eta, causal = causal),
  #           simulate = function(mu, sigma2, eta, kin, n, nsim) {
  #             P <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = eta * sigma2 * kin)
  #             E <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))
  #             # y <- mu + sigma * matrix(rnorm(nsim * n), n, nsim)
  #             y <- b0 + mu + t(P) + t(E)
  #             return(split(y, col(y))) # make each col its own list element
  #           })
}




# this function is similar to the one used to simulate the data included in the ggmix package
# see data-raw/bnpsd-data.R. I copied over this function to the R/utils.R file, because i modified it
# a little bit so that it doenst plot anyhting. and also in case I remove it here. that function
# will also be included in the package to simulate data
make_ADmixed_model_not_simulator <- function(n, p_test, p_kinship, k, s, Fst, b0, beta_mean,
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
    n1 <- 200; n2 <- 200; n3 <- 200; n4 <- 200; n5 <- 200
    # here’s the labels (for simplicity, list all individuals of S1 first, then S2, then S3)
    labs <- c( rep.int("S1", n1), rep.int("S2", n2), rep.int("S3", n3),
               rep.int("S4", n4), rep.int("S5", n5))
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

    snps_design <- c(setdiff(cnames, snps_kinships), causal)
    # length(snps_design)
    # setdiff(cnames, snps_kinships) %>% length()
    not_causal <- setdiff(snps_design, causal)

    Xkinship <- Xall[,snps_kinships]
    Xtest <- Xall[,snps_design]

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
    causal <- sample(snps_design, ncausal, replace = FALSE)
    not_causal <- setdiff(snps_design, causal)

    Xkinship <- Xall[,snps_kinships]
    Xtest <- Xall[,snps_design]

    # now estimate kinship using popkin
    # PhiHat <- popkin::popkin(X, subpops, lociOnCols = TRUE)
    PhiHat <- popkin::popkin(Xkinship, lociOnCols = TRUE)
  }

  kin <- 2 * PhiHat
  eiK <- eigen(kin)
  if (any(eiK$values < 1e-5)) { eiK$values[ eiK$values < 1e-5 ] <- 1e-5 }
  PC <- sweep(eiK$vectors, 2, sqrt(eiK$values), "*")
  plot(eiK$values)
  plot(PC[,1],PC[,2], pch = 19, col = rep(RColorBrewer::brewer.pal(5,"Paired"), each = 200))

  np <- dim(Xtest)
  n <- np[[1]]
  p <- np[[2]]

  x_lasso <- cbind(Xtest,PC[,1:10])
  x_lasso[1:5,1:5]
  # kin <- snpStats::xxt(dat$genotypes)/p

  beta <- rep(0, length = p)
  # beta[which(colnames(Xtest) %in% causal)] <- runif(n = length(causal), beta_mean - 0.1, beta_mean + 0.1)
  beta[which(colnames(Xtest) %in% causal)] <- rnorm(n = length(causal))
  mu <- as.numeric(Xtest %*% beta)

  # browser()
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

  return(list(y = y, x = Xtest, causal = causal, beta = beta, kin = kin,
              not_causal = not_causal,
              x_lasso = x_lasso))

  # new_model(name = "penfam", label = sprintf("type = %s, eta = %s", type, eta),
  #           params = list(mu = mu, n = n, x = x, beta = beta, type = type, not_causal = not_causal,
  #                         kin = kin, b0 = b0, sigma2 = sigma2, eta = eta, causal = causal),
  #           simulate = function(mu, sigma2, eta, kin, n, nsim) {
  #             P <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = eta * sigma2 * kin)
  #             E <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))
  #             # y <- mu + sigma * matrix(rnorm(nsim * n), n, nsim)
  #             y <- b0 + mu + t(P) + t(E)
  #             return(split(y, col(y))) # make each col its own list element
  #           })
}






make_ADmixed_model_not_simulator_with_validation <- function(n_train, p_test, p_kinship, k, s, Fst, b0, beta_mean,
                                                             n_validation,
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

  # n_total <- n_train + n_validation
  # train_ind <- sample(1:n_total, n_train)

  FF <- 1:k # subpopulation FST vector, up to a scalar
  # s <- 0.5 # desired bias coefficient
  # Fst <- 0.1 # desired FST for the admixed individuals
  geography <- match.arg(geography)
  if (geography == "1d") {
    obj <- bnpsd::q1d(n = n_train, k = k, s = s, F = FF, Fst = Fst)
    Q <- obj$Q
    FF <- obj$F

    obj_val <- bnpsd::q1d(n = n_validation, k = k, s = s, F = FF, Fst = Fst)
    Q_val <- obj_val$Q
    FF_val <- obj_val$F

  } else if (geography == "ind") {

    n_train_k <- n_train / k
    n1 <- n_train_k; n2 <- n_train_k; n3 <- n_train_k; n4 <- n_train_k; n5 <- n_train_k
    # here’s the labels (for simplicity, list all individuals of S1 first, then S2, then S3)
    labs <- c( rep.int("S1", n1), rep.int("S2", n2), rep.int("S3", n3),
               rep.int("S4", n4), rep.int("S5", n5))
    # data dimensions infered from labs:
    length(labs) # number of individuals "n"
    # desired admixture matrix ("is" stands for "Independent Subpopulations")
    Q <- bnpsd::qis(labs)
    FF <- 1:k # subpopulation FST vector, unnormalized so far
    FF <- FF/popkin::fst(FF)*Fst # normalized to have the desired Fst


    n_val_k <- n_validation / k
    n1_val <- n_val_k; n2_val <- n_val_k; n3_val <- n_val_k; n4_val <- n_val_k; n5_val <- n_val_k
    # here’s the labels (for simplicity, list all individuals of S1 first, then S2, then S3)
    labs_val <- c( rep.int("S1", n1_val), rep.int("S2", n2_val), rep.int("S3", n3_val),
                   rep.int("S4", n4_val), rep.int("S5", n5_val))

    Q_val <- bnpsd::qis(labs_val)
    FF_val <- 1:k # subpopulation FST vector, unnormalized so far
    FF_val <- FF_val/popkin::fst(FF_val)*Fst # normalized to have the desired Fst

  } else if (geography == "circ") {
    obj <- bnpsd::q1dc(n = n_train, k = k, s = s, F = FF, Fst = Fst)
    Q <- obj$Q
    FF <- obj$F

    obj_val <- bnpsd::q1dc(n = n_validation, k = k, s = s, F = FF, Fst = Fst)
    Q_val <- obj_val$Q
    FF_val <- obj_val$F # this will be the same as FF
  }


  ncausal <- p_test * percent_causal
  # browser()
  if (percent_overlap == "100") {
# browser()
    total_snps_to_simulate <- p_test + p_kinship - ncausal
    # this contains all SNPs (X_{Testing}:X_{kinship})

    pAnc <- bnpsd::rpanc(total_snps_to_simulate) # random vector of ancestral allele frequencies (length= number of loci)
    B <- bnpsd::rpint(pAnc, FF) # matrix of intermediate subpop allele freqs
    P <- bnpsd::rpiaf(B = B,Q = Q)
    out <- bnpsd::rgeno(P)

    Xall <- t(out) # genotypes are columns, rows are subjects
    cnames <- paste0("X", 1:total_snps_to_simulate)
    colnames(Xall) <- cnames
    rownames(Xall) <- paste0("id", 1:n_train)

    B_val <- bnpsd::rpint(pAnc, FF_val) # matrix of intermediate subpop allele freqs (FF is same as FF_val)
    # head(B_val);head(B)
    # sigma <- 1 # dispersion parameter of intermediate subpops
    # Q <- q1d(n, k, sigma) # non-trivial admixture proportions
    # head(obj$Q);head(obj_val$Q) # also the same, but would change depending on if n_train is different from n_validation
    # tail(obj$Q);tail(obj_val$Q)
    P_val <- bnpsd::rpiaf(B = B_val, Q = Q_val)
    # P_val[1:5,1:5] ; P[1:5,1:5]
    out_val <- bnpsd::rgeno(P_val)

    Xall_val <- t(out_val) # genotypes are columns, rows are subjects
    colnames(Xall_val) <- cnames
    rownames(Xall_val) <- paste0("id", 1:n_validation)

    # tra <- gaston::as.bed.matrix(Xall)
    # tra@p %>% plot
    # val <- gaston::as.bed.matrix(Xall_val)
    # dev.off()
    # these plots show the train and val MAFs are assez correler
    # plot(tra@p, val@p)
    # cor(tra@p, val@p)

    # subpops <- ceiling( (1:n)/n*k )
    # table(subpops) # got k=10 subpops with 100 individuals each

    # Snps used for kinship
    snps_kinships <- sample(cnames, p_kinship, replace = FALSE)
    length(snps_kinships)

    # all causal snps are in kinship matrix
    causal <- sample(snps_kinships, ncausal, replace = FALSE)
    length(causal)
    snps_design <- c(setdiff(cnames, snps_kinships), causal)
    # length(snps_design)
    # setdiff(cnames, snps_kinships) %>% length()
    not_causal <- setdiff(snps_design, causal)

    Xkinship <- Xall[,snps_kinships]
    Xtest <- Xall[,snps_design]

    Xkinship_val <- Xall_val[,snps_kinships]
    Xtest_val <- Xall_val[,snps_design]

    # now estimate kinship using popkin
    # PhiHat <- popkin::popkin(X, subpops, lociOnCols = TRUE)
    PhiHat <- popkin::popkin(Xkinship, lociOnCols = TRUE)
    PhiHat_val <- popkin::popkin(Xkinship_val, lociOnCols = TRUE)
    # popkin::plotPopkin(list(PhiHat,PhiHat_val))

  } else if (percent_overlap == "0") {

    total_snps_to_simulate <- p_test + p_kinship
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

    snps_design <- setdiff(cnames, snps_kinships)
    # length(snps_design)
    # setdiff(cnames, snps_kinships) %>% length()
    causal <- sample(snps_design, ncausal, replace = FALSE)
    not_causal <- setdiff(snps_design, causal)

    Xkinship <- Xall[,snps_kinships]
    Xtest <- Xall[,snps_design]

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

  kin_val <- 2 * PhiHat_val
  eiK_val <- eigen(kin_val)
  if (any(eiK_val$values < 1e-5)) { eiK_val$values[ eiK_val$values < 1e-5 ] <- 1e-5 }
  PC_val <- sweep(eiK_val$vectors, 2, sqrt(eiK_val$values), "*")
  # plot(eiK_val$values)

  par(mfrow = c(1,2))
  plot(PC[,1],PC[,2], pch = 19, col = rep(RColorBrewer::brewer.pal(5,"Paired"), each = 200), main = "training set")
  plot(PC_val[,1],PC_val[,2], pch = 19, col = rep(RColorBrewer::brewer.pal(5,"Paired"), each = 200), main = "test set")

  np <- dim(Xtest)
  n <- np[[1]]
  p <- np[[2]]

  x_lasso <- cbind(Xtest,PC[,1:10])
  x_lasso[1:5,1:5]

  # kin <- snpStats::xxt(dat$genotypes)/p

  beta <- rep(0, length = p)
  # beta_val <- rep(0, length = p)
  # beta[which(colnames(Xtest) %in% causal)] <- runif(n = length(causal), beta_mean - 0.1, beta_mean + 0.1)
  beta[which(colnames(Xtest) %in% causal)] <- rnorm(n = length(causal))
  # beta_val[which(colnames(Xtest) %in% causal)] <- rnorm(n = length(causal))

  mu <- as.numeric(Xtest %*% beta)
  mu_val <- as.numeric(Xtest_val %*% beta)

  P <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = eta * sigma2 * kin)
  E <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))
  y <- b0 + mu + P + E

  P <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = eta * sigma2 * kin_val)
  E <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))
  y_val <- b0 + mu_val + P + E

  return(list(y = y, x = Xtest, causal = causal, beta = beta, kin = kin,
              y_val = y_val, x_val = Xtest_val,
              not_causal = not_causal,
              x_lasso = x_lasso))

  # new_model(name = "penfam", label = sprintf("type = %s, eta = %s", type, eta),
  #           params = list(mu = mu, n = n, x = x, beta = beta, type = type, not_causal = not_causal,
  #                         kin = kin, b0 = b0, sigma2 = sigma2, eta = eta, causal = causal),
  #           simulate = function(mu, sigma2, eta, kin, n, nsim) {
  #             P <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = eta * sigma2 * kin)
  #             E <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))
  #             # y <- mu + sigma * matrix(rnorm(nsim * n), n, nsim)
  #             y <- b0 + mu + t(P) + t(E)
  #             return(split(y, col(y))) # make each col its own list element
  #           })
}



make_INDmixed_model_not_simulator <- function(n, p, ncausal, k, s, Fst, b0, beta_mean, eta, sigma2) {

  # k:	Number of intermediate subpopulations
  # s: The desired bias coefficient, which specifies σ indirectly. Required if sigma is missing
  # F: The length-k vector of inbreeding coefficients (or FST's) of the intermediate subpopulations,
  # up to a scaling factor (which cancels out in calculations). Required if sigma is missing
  # Fst: The desired final FST of the admixed individuals. Required if sigma is missing
  # browser()
  # define population structure
  n1 <- 200; n2 <- 200; n3 <- 200; n4 <- 200; n5 <- 200
  # here’s the labels (for simplicity, list all individuals of S1 first, then S2, then S3)
  labs <- c( rep.int("S1", n1), rep.int("S2", n2), rep.int("S3", n3), rep.int("S4", n4), rep.int("S5", n5))
  # data dimensions infered from labs:
  length(labs) # number of individuals "n"
  # desired admixture matrix ("is" stands for "Independent Subpopulations")
  Q <- bnpsd::qis(labs)

  FF <- 1:k # subpopulation FST vector, unnormalized so far
  FF <- FF/popkin::fst(FF)*Fst # normalized to have the desired Fst
  # s <- 0.5 # desired bias coefficient
  # Fst <- 0.1 # desired FST for the admixed individuals
  # obj <- bnpsd::q1d(n = n, k = k, s = s, F = FF, Fst = Fst) # admixture proportions from 1D geography
  # Q <- obj$Q
  # FF <- obj$F
  out <- bnpsd::rbnpsd(Q, FF, p)
  X <- t(out$X) # genotypes are columns, rows are subjects
  dim(X)
  colnames(X) <- paste0("X",1:ncol(X))
  rownames(X) <- paste0("id",1:nrow(X))
  dim(X)
  X[1:5,1:5]
  subpops <- ceiling( (1:n)/n*k )
  table(subpops) # got k=10 subpops with 100 individuals each
  # now estimate kinship using popkin
  # PhiHat <- popkin::popkin(X, subpops = c(rep(1, n1), rep(2, n2), rep(3, n3), rep(4, n4), rep(5, n5)), lociOnCols = TRUE)
  PhiHat <- popkin::popkin(X, lociOnCols = TRUE)
  # PhiHat[1:5,1:5]
  kin <- 2 * PhiHat
  # kin <- PhiHat
  # kin <- gaston::GRM(gaston::as.bed.matrix(X), autosome.only = FALSE)
  # inbrDiag(PhiHat)
  # dev.off()
  # plotPopkin(list(PhiHat, PhiHat2))

  isPD <- all(eigen(kin)$values > 0)
  how_many_neg_eigenvalues <- sum(eigen(kin)$values <= 0)
  if (!isPD) {
    kinPD <- as(nearPD(kin)$mat,"matrix")
    kin <- kinPD
  }


  kin[1:5,1:5]
  dim(PhiHat)
  eiK <- eigen(kin)
  # all(rownames(as.matrix(x))==rownames(kin))
  # deal with a small negative eigen value
  plot(eiK$values)
  sum(eiK$values < 0)
  if (any(eiK$values < 0)) { eiK$values[ eiK$values < 0 ] <- 0 }
  PC <- sweep(eiK$vectors, 2, sqrt(eiK$values), "*")
  # dev.off()
  plot(eiK$values)
  plot(PC[,1],PC[,2], pch = 19, col = rep(RColorBrewer::brewer.pal(5,"Paired"), each = n1))
  # X <- t(out$X)
  # dim(X)
  #Phi;Phi_names;bedfile;causal_list
  # browser()

  # need to re-order
  # all(rownames(X)==rownames(kin))
  # all(rownames(X) %in% rownames(kin))
  # X <- X[match(rownames(kin), rownames(X)),]
  # all(rownames(X)==rownames(kin))


  np <- dim(X)
  n <- np[[1]]
  p <- np[[2]]



  x_lasso <- cbind(X,PC[,1:10])
  x_lasso[1:5,1:5]
  # kin <- snpStats::xxt(dat$genotypes)/p
  causal <- sample(colnames(X), ncausal, replace = FALSE)
  not_causal <- setdiff(colnames(X), causal)
  beta <- rep(0, length = p)
  beta[which(colnames(X) %in% causal)] <- runif(n = length(causal), beta_mean - 0.2, beta_mean + 0.2)
  # beta[which(colnames(x) %in% causal)] <- rnorm(n = length(causal))
  mu <- as.numeric(X %*% beta)


  P <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = eta * sigma2 * kin)
  E <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))

  # P <- MASS::mvrnorm(1, rep(0,n), 1.2625 * kin)

  # y <- mu + sigma * matrix(rnorm(nsim * n), n, nsim)
  # y <- b0 + mu + t(P) + t(E)
  # y <- MASS::mvrnorm(1, mu = mu, Sigma = eta * sigma2 * kin + (1 - eta) * sigma2 * diag(n))
  # y <- mu + P + rnorm(n, 0, 1)
  y <- 0 + mu + P + E

  return(list(y = y, x = X, causal = causal, beta = beta, kin = kin,
              not_causal = not_causal,
              x_lasso = x_lasso))

  # new_model(name = "penfam", label = sprintf("type = %s, eta = %s", type, eta),
  #           params = list(mu = mu, n = n, x = x, beta = beta, type = type, not_causal = not_causal,
  #                         kin = kin, b0 = b0, sigma2 = sigma2, eta = eta, causal = causal),
  #           simulate = function(mu, sigma2, eta, kin, n, nsim) {
  #             P <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = eta * sigma2 * kin)
  #             E <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))
  #             # y <- mu + sigma * matrix(rnorm(nsim * n), n, nsim)
  #             y <- b0 + mu + t(P) + t(E)
  #             return(split(y, col(y))) # make each col its own list element
  #           })
}



# not used gaston testing code --------------------------------------------

# gaston_mat
# head(gaston_mat@ped)
#
# ld.x <- gaston::LD(gaston_mat, c(1,ncol(gaston_mat)))
# standardize(gaston_mat)
#
# plot(gaston_mat@p, gaston_mat@sigma, xlim=c(0,1))
# t <- seq(0,1,length=101);
# lines(t, sqrt(2*t*(1-t)), col="red")
#
# y <- gaston::LD.thin(gaston_mat, threshold = 0.4, max.dist = 500e3)
# y
#
# gaston_mat <- set.stats(gaston_mat)
# head(gaston_mat@ped)
# head(gaston_mat@snps)
#
#
# gaston_mat@snps[1:5,1:100]
#
# gaston_mat <- gaston::read.bed.matrix(file_paths$bedfile_for_kinship)
# grm_gaston <- gaston::GRM(gaston_mat)
#
# standardize(gaston_mat) <- "p"
# # head(grm_gaston)
# # grm_gaston[1:5,1:5]
# X <- as.matrix(gaston_mat)
# head(gaston_mat@bed)
# any(is.na(gaston_mat@p))
# X[1:5,1:5]
# X[is.na(X)] <- 0
# X[1:5,1:5]
# colMeans(X) %>% round(100) %>% plot
# apply(X, 2, sd) %>% plot
# sqrt(2*gaston_mat@p*(1-gaston_mat@p)) %>% plot
#
# col(X) %>% round(100) %>% plot
# scale(X)[1:5,1:5]
# scale(X, center = 2*gaston_mat@p,scale = sqrt(2*gaston_mat@p*(1-gaston_mat@p)))[1:5,1:5]
#
# x <- as.bed.matrix(TTN.gen, TTN.fam, TTN.bim)
# X <- as.matrix(x)
# X[1:5,1:4]
# standardize(x) <- "p"
# as.matrix(x)[1:5, 1:4]
# colMeans(as.matrix(x)) %>% round(10) %>% plot
# apply(as.matrix(x), 2, sd) %>% plot
# pheatmap::pheatmap(grm_gaston, show_rownames = F, show_colnames = F, color = rev(viridis::viridis(100)))
#
#
# eiK <- eigen(grm_gaston)
# any(eiK$values < 0)
# eiK$values[ eiK$values < 0 ] <- 0
# plot(eiK$vectors[,1], eiK$vectors[,2])
#
# PC <- sweep(eiK$vectors, 2, sqrt(eiK$values), "*")
# plot(PC[,1],PC[,2])
# plot(PC[,2],PC[,3])
#
# all(eigen(grm_gaston)$values > 0)
# eigen(grm_gaston)$values[which(eigen(grm_gaston)$values < 0)]
#
# svd()


