## @knitr models


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




## @knitr models-not-used

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





make_ADmixed_model_not_simulator <- function(n, p, ncausal, k, s, Fst, b0, beta_mean, eta, sigma2) {

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
  obj <- bnpsd::q1d(n = n, k = k, s = s, F = FF, Fst = Fst) # admixture proportions from 1D geography
  Q <- obj$Q
  FF <- obj$F
  out <- bnpsd::rbnpsd(Q, FF, p)
  X <- t(out$X) # genotypes are columns, rows are subjects
  dim(X)
  colnames(X) <- paste0("X",1:p)
  rownames(X) <- paste0("id",1:n)
  dim(X)
  X[1:5,1:5]
  subpops <- ceiling( (1:n)/n*k )
  table(subpops) # got k=10 subpops with 100 individuals each
  # now estimate kinship using popkin
  PhiHat <- popkin::popkin(X, subpops, lociOnCols = TRUE)
  PhiHat[1:5,1:5]
  kin <- 2 *PhiHat

  kin[1:5,1:5]
  dim(PhiHat)
  eiK <- eigen(kin)
  # all(rownames(as.matrix(x))==rownames(kin))
  # deal with a small negative eigen value
  if (any(eiK$values < 0)) { eiK$values[ eiK$values < 0 ] <- 0 }
  PC <- sweep(eiK$vectors, 2, sqrt(eiK$values), "*")
  # dev.off()
  plot(eiK$values)
  plot(PC[,1],PC[,2], pch = 19, col = rep(RColorBrewer::brewer.pal(10,"Paired"), each = 10))
  # X <- t(out$X)
  # dim(X)
  #Phi;Phi_names;bedfile;causal_list
  # browser()

  np <- dim(X)
  n <- np[[1]]
  p <- np[[2]]


  x_lasso <- cbind(X,PC[,1:10])
  x_lasso[1:5,1:5]
  # kin <- snpStats::xxt(dat$genotypes)/p
  causal <- sample(colnames(X), ncausal, replace = FALSE)
  not_causal <- setdiff(colnames(X), causal)
  beta <- rep(0, length = p)
  beta[which(colnames(X) %in% causal)] <- runif(n = length(causal), beta_mean - 0.1, beta_mean + 0.1)
  # beta[which(colnames(x) %in% causal)] <- rnorm(n = length(causal))
  mu <- as.numeric(X %*% beta)


  P <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = eta * sigma2 * kin)
  E <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))
  # y <- mu + sigma * matrix(rnorm(nsim * n), n, nsim)
  # y <- b0 + mu + t(P) + t(E)
  y <- MASS::mvrnorm(1, mu = mu, Sigma = eta * sigma2 * kin + (1 - eta) * sigma2 * diag(n))
  # y <- 0 + mu + P + E

  return(list(y = y, x = X, causal = causal, beta = beta, kin = kin,
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
  n1 <- 100; n2 <- 100; n3 <- 100
  # here’s the labels (for simplicity, list all individuals of S1 first, then S2, then S3)
  labs <- c( rep.int("S1", n1), rep.int("S2", n2), rep.int("S3", n3) )
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
  PhiHat <- popkin::popkin(X, subpops = c(rep(1, n1), rep(2, n2), rep(3, n3)), lociOnCols = TRUE)
  # PhiHat2 <- popkin::popkin(X, lociOnCols = TRUE)
  # PhiHat[1:5,1:5]
  kin <- 2 *PhiHat
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
  plot(PC[,1],PC[,2], pch = 19, col = rep(RColorBrewer::brewer.pal(3,"Paired"), each = n1))
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


  # P <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = eta * sigma2 * kin)
  # E <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))

  P <- MASS::mvrnorm(1, rep(0,n), 1.2625 * kin)

  # y <- mu + sigma * matrix(rnorm(nsim * n), n, nsim)
  # y <- b0 + mu + t(P) + t(E)
  # y <- MASS::mvrnorm(1, mu = mu, Sigma = eta * sigma2 * kin + (1 - eta) * sigma2 * diag(n))
  y <- mu + P + rnorm(n, 0, 1)
  # y <- 0 + mu + P + E

  return(list(y = y, x = X, causal = causal, beta = beta, kin = kin,
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
