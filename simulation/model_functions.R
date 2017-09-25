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





make_mixed_model_not_simulator <- function(b0, eta, sigma2, type, related = TRUE) {
  file_paths <- if (!related) {
    switch(type,
           causal_400 = {
             list(Phi = "/home/sahir/git_repositories/penfam/datahydra/kinship/Kinship_other_3600_causal_400.rel",
                  Phi_names = "/home/sahir/git_repositories/penfam/datahydra/kinship/Kinship_other_3600_causal_400.rel.id",
                  bedfile = "/home/sahir/git_repositories/penfam/datahydra/bed/X_test_4k",
                  bedfile_for_kinship = "/home/sahir/git_repositories/penfam/datahydra/bed/X_other_3600_causal_400",
                  causal_list = "/home/sahir/git_repositories/penfam/datahydra/snplist/causal_SNP_from_X_test_400")
           },
           causal_2k = {
             list(Phi = "/home/sahir/git_repositories/penfam/datahydra/kinship/Kinship_other_2000_causal_2000.rel",
                  Phi_names = "/home/sahir/git_repositories/penfam/datahydra/kinship/Kinship_other_2000_causal_2000.rel.id",
                  bedfile = "/home/sahir/git_repositories/penfam/datahydra/bed/X_test_4k",
                  bedfile_for_kinship = "/home/sahir/git_repositories/penfam/datahydra/bed/X_other_2000_causal_2000",
                  causal_list = "/home/sahir/git_repositories/penfam/datahydra/snplist/causal_SNP_from_X_test_2000")
           },
           causal_4k = {
             list(Phi = "/home/sahir/git_repositories/penfam/datahydra/kinship/Kinship_causal_4000.rel",
                  Phi_names = "/home/sahir/git_repositories/penfam/datahydra/kinship/Kinship_causal_4000.rel.id",
                  bedfile = "/home/sahir/git_repositories/penfam/datahydra/bed/X_test_4k",
                  bedfile_for_kinship = "/home/sahir/git_repositories/penfam/datahydra/bed/X_test_4k",
                  causal_list = "/home/sahir/git_repositories/penfam/datahydra/snplist/causal_SNP_from_X_test_4000")
           })
  } else {
    switch(type,
           causal_400 = {
             list(Phi = "/home/sahir/git_repositories/penfam/datahydra_related/kinship/Kinship_other_3600_causal_400.rel",
                  Phi_names = "/home/sahir/git_repositories/penfam/datahydra_related/kinship/Kinship_other_3600_causal_400.rel.id",
                  bedfile = "/home/sahir/git_repositories/penfam/datahydra_related/bed/X_test_4k",
                  bedfile_for_kinship = "/home/sahir/git_repositories/penfam/datahydra_related/bed/X_other_3600_causal_400",
                  causal_list = "/home/sahir/git_repositories/penfam/datahydra_related/snplist/causal_SNP_from_X_test_400")
           },
           causal_2k = {
             list(Phi = "/home/sahir/git_repositories/penfam/datahydra_related/kinship/Kinship_other_2000_causal_2000.rel",
                  Phi_names = "/home/sahir/git_repositories/penfam/datahydra_related/kinship/Kinship_other_2000_causal_2000.rel.id",
                  bedfile = "/home/sahir/git_repositories/penfam/datahydra_related/bed/X_test_4k",
                  bedfile_for_kinship = "/home/sahir/git_repositories/penfam/datahydra_related/bed/X_other_2000_causal_2000",
                  causal_list = "/home/sahir/git_repositories/penfam/datahydra_related/snplist/causal_SNP_from_X_test_2000")
           },
           causal_4k = {
             list(Phi = "/home/sahir/git_repositories/penfam/datahydra_related/kinship/Kinship_causal_4000.rel",
                  Phi_names = "/home/sahir/git_repositories/penfam/datahydra_related/kinship/Kinship_causal_4000.rel.id",
                  bedfile = "/home/sahir/git_repositories/penfam/datahydra_related/bed/X_test_4k",
                  bedfile_for_kinship = "/home/sahir/git_repositories/penfam/datahydra_related/bed/X_test_4k",
                  causal_list = "/home/sahir/git_repositories/penfam/datahydra_related/snplist/causal_SNP_from_X_test_4000")
           })
    }

  # file_paths <- NULL
  # file_paths$bedfile_for_kinship <- "/home/sahir/git_repositories/penfam/datahydra_related/bed/X_other_3600_causal_400"
  # file_paths$Phi_names <- "/home/sahir/git_repositories/penfam/datahydra_related/kinship/Kinship_other_3600_causal_400.rel.id"
  # file_paths$bedfile <- "/home/sahir/git_repositories/penfam/datahydra_related/bed/X_test_4k"
  # ==================================

# Kinship Matrix ----------------------------------------------------------

  gaston_mat <- gaston::read.bed.matrix(file_paths$bedfile_for_kinship)
  gaston::standardize(gaston_mat) <- "p"
  kin <- gaston::GRM(gaston_mat)

  isPD <- all(eigen(kin)$values > 0)
  how_many_neg_eigenvalues <- sum(eigen(kin)$values <= 0)
  if (!isPD) {
    kinPD <- as(Matrix::nearPD(kin)$mat,"matrix")
    dimnames(kinPD)[[1]] <- kin_names$V1
    dimnames(kinPD)[[2]] <- kin_names$V1
    kin <- kinPD
  }

  # this is the W matrix (i.e. the snps used to construct the kinship matrix)
  w <- as.matrix(gaston_mat)

  # In standardized matrices, the NA values are replaced by zeroes, which amount to impute the missing
  # genotypes by the mean genotype
  w[is.na(w)] <- 0


# X matrix ----------------------------------------------------------------

  dat <- gaston::read.bed.matrix(file_paths$bedfile)
  gaston::standardize(dat) <- "p"
  x <- as.matrix(dat)
  x[is.na(x)] <- 0

  np <- dim(x)
  n <- np[[1]]
  p <- np[[2]]

  causal <- data.table::fread(file_paths$causal_list, header = FALSE)$V1
  not_causal <- setdiff(colnames(x), causal)
  beta <- rep(0, length = p)
  beta[which(colnames(x) %in% causal)] <- runif(n = length(causal), 0.1, 0.3)
  mu <- as.numeric(x %*% beta)

  P <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = eta * sigma2 * kin)
  E <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))

  y <- b0 + mu + P + E

  return(list(y = y, x = x, causal = causal, beta = beta, kin = kin, isPD = isPD, w = w,
              neg_eigen = how_many_neg_eigenvalues))
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
