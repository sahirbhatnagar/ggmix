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



make_mixed_model_SSC <- function(b0, eta, sigma2, percent_causal, percent_overlap) {
  
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
  
  # kin <- as.matrix(read.table(file_paths$Phi))
  # kin_names <- data.table::fread(file_paths$Phi_names)
  # dimnames(kin)[[1]] <- kin_names$V1
  # dimnames(kin)[[2]] <- kin_names$V1
  # isPD <- all(eigen(kin)$values > 0)
  # how_many_neg_eigenvalues <- sum(eigen(kin)$values <= 0)
  # if (!isPD) {
  #   kinPD <- as(nearPD(kin)$mat,"matrix")
  #   dimnames(kinPD)[[1]] <- kin_names$V1
  #   dimnames(kinPD)[[2]] <- kin_names$V1
  #   kin <- kinPD
  # }
  # dat <- snpStats::read.plink(file_paths$bedfile)
  # x <- as(dat$genotypes, "numeric")
  # np <- dim(x)
  # n <- np[[1]]
  # p <- np[[2]]
  # pcaCars <- princomp(kin, cor = TRUE)
  # # # view objects stored in pcaCars
  # # names(pcaCars)
  # # # proportion of variance explained
  # # summary(pcaCars)
  # # # scree plot
  # # plot(pcaCars, type = "l")
  # # first 10 PCs
  # x_lasso <- cbind(x,pcaCars$scores[,1:10])
  # 
  # # kin <- snpStats::xxt(dat$genotypes)/p
  # 
  # causal <- data.table::fread(file_paths$causal_list, header = FALSE)$V1
  # not_causal <- setdiff(colnames(x), causal)
  # beta <- rep(0, length = p)
  # beta[which(colnames(x) %in% causal)] <- runif(n = length(causal), 0.1, 0.3)
  # mu <- as.numeric(x %*% beta)
  
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
  # kin <- snpStats::xxt(dat$genotypes)/p
  
  causal <- data.table::fread(file_paths$causal_list, header = FALSE)$V1
  not_causal <- setdiff(colnames(X), causal)
  beta <- rep(0, length = p)
  beta[which(colnames(X) %in% causal)] <- runif(n = length(causal), 0.9, 1.1)
  # beta[which(colnames(x) %in% causal)] <- rnorm(n = length(causal))
  mu <- as.numeric(X %*% beta)
  
  new_model(name = "ggmixSSCv2", label = sprintf("percent_causal = %s, percent_overlap = %s, eta = %s, sigma = %s", 
                                               percent_causal, percent_overlap, eta, sigma2),
            params = list(mu = mu, n = n, x = X, x_lasso = x_lasso, 
                          beta = beta, percent_causal = percent_causal, 
                          percent_overlap = percent_overlap,
                          not_causal = not_causal,
                          kin = kin, b0 = b0, sigma2 = sigma2, eta = eta, causal = causal),
            simulate = function(mu, sigma2, eta, kin, n, nsim) {
              P <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = eta * sigma2 * kin)
              E <- MASS::mvrnorm(nsim, mu = rep(0, n), Sigma = (1 - eta) * sigma2 * diag(n))
              # y <- mu + sigma * matrix(rnorm(nsim * n), n, nsim)
              y <- b0 + mu + t(P) + t(E)
              return(split(y, col(y))) # make each col its own list element
            })
}

make_mixed_model_not_simulator <- function(b0, eta, sigma2,  percent_causal, percent_overlap) {
  
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
  y <- 0 + mu + P + E
  
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
