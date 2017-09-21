## @knitr models


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
                         list(Phi = "/home/sahir/git_repositories/penfam/datahydra/kinship/Kinship_other_3600_causal_400.rel",
                              Phi_names = "/home/sahir/git_repositories/penfam/datahydra/kinship/Kinship_other_3600_causal_400.rel.id",
                              bedfile = "/home/sahir/git_repositories/penfam/datahydra/bed/X_test_4k",
                              causal_list = "/home/sahir/git_repositories/penfam/datahydra/snplist/causal_SNP_from_X_test_400")
                       },
                       causal_2k = {
                         list(Phi = "/home/sahir/git_repositories/penfam/datahydra/kinship/Kinship_other_2000_causal_2000.rel",
                              Phi_names = "/home/sahir/git_repositories/penfam/datahydra/kinship/Kinship_other_2000_causal_2000.rel.id",
                              bedfile = "/home/sahir/git_repositories/penfam/datahydra/bed/X_test_4k",
                              causal_list = "/home/sahir/git_repositories/penfam/datahydra/snplist/causal_SNP_from_X_test_2000")
                       },
                       causal_4k = {
                         list(Phi = "/home/sahir/git_repositories/penfam/datahydra/kinship/Kinship_causal_4000.rel",
                              Phi_names = "/home/sahir/git_repositories/penfam/datahydra/kinship/Kinship_causal_4000.rel.id",
                              bedfile = "/home/sahir/git_repositories/penfam/datahydra/bed/X_test_4k",
                              causal_list = "/home/sahir/git_repositories/penfam/datahydra/snplist/causal_SNP_from_X_test_4000")
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





make_mixed_model_not_simulator <- function(b0, eta, sigma2, type) {
  file_paths <- switch(type,
                       causal_400 = {
                         list(Phi = "/home/sahir/git_repositories/penfam/datahydra/kinship/Kinship_other_3600_causal_400.rel",
                              Phi_names = "/home/sahir/git_repositories/penfam/datahydra/kinship/Kinship_other_3600_causal_400.rel.id",
                              bedfile = "/home/sahir/git_repositories/penfam/datahydra/bed/X_test_4k",
                              causal_list = "/home/sahir/git_repositories/penfam/datahydra/snplist/causal_SNP_from_X_test_400")
                       },
                       causal_2k = {
                         list(Phi = "/home/sahir/git_repositories/penfam/datahydra/kinship/Kinship_other_2000_causal_2000.rel",
                              Phi_names = "/home/sahir/git_repositories/penfam/datahydra/kinship/Kinship_other_2000_causal_2000.rel.id",
                              bedfile = "/home/sahir/git_repositories/penfam/datahydra/bed/X_test_4k",
                              causal_list = "/home/sahir/git_repositories/penfam/datahydra/snplist/causal_SNP_from_X_test_2000")
                       },
                       causal_4k = {
                         list(Phi = "/home/sahir/git_repositories/penfam/datahydra/kinship/Kinship_causal_4000.rel",
                              Phi_names = "/home/sahir/git_repositories/penfam/datahydra/kinship/Kinship_causal_4000.rel.id",
                              bedfile = "/home/sahir/git_repositories/penfam/datahydra/bed/X_test_4k",
                              causal_list = "/home/sahir/git_repositories/penfam/datahydra/snplist/causal_SNP_from_X_test_4000")
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
  y <- b0 + mu + P + E

  return(list(y = y, x = x, causal = causal, beta = beta, kin = kin, isPD = isPD,
              neg_eigen = how_many_neg_eigenvalues))

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
