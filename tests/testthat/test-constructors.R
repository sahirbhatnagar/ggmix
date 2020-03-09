context(strwrap("check full and low rank constructors have valid output using
        admixed dataset"))

data("admixed")

svdXkinship <- svd(admixed$Xkinship)
number_nonzero_eigenvalues <- 20

fullrank_kinship <- try(new_fullrank_kinship(x = admixed$xtrain, y = admixed$ytrain,
                                             kinship = admixed$kin_train),
                        silent = TRUE)

fullrank_K <- try(new_fullrank_K(x = admixed$xtrain, y = admixed$ytrain,
                                 K = admixed$Xkinship),
                  silent = TRUE)

fullrank_UD <- try(new_fullrank_UD(x = admixed$xtrain, y = admixed$ytrain,
                                   U = svdXkinship$u[,which(svdXkinship$d > 0), drop = F],
                                   D = svdXkinship$d[which(svdXkinship$d > 0)]),
                   silent = TRUE)


test_that("no error in full rank constructor functions", {

  expect_false(inherits(fullrank_kinship, "try-error"))
  expect_false(inherits(fullrank_K, "try-error"))
  expect_false(inherits(fullrank_UD, "try-error"))
  expect_is(fullrank_kinship, class = "fullrank")
  expect_is(fullrank_K, class = "fullrank")
  expect_is(fullrank_UD, class = "fullrank")

  # expect_equivalent(class(coef(fit_sim)), "dgCMatrix")
  # expect_equal(dim(coef(fit_sim))[[1]], dim(fit_sim$design)[[2]] + 1)
  # expect_equal(dim(coef(fit_sim))[[2]], sum(fit_sim$converged))

})


# kin <- gaston::GRM(gaston::as.bed.matrix(admixed$Xkinship), autosome.only = FALSE)
# need to have pacman::p_load(RSpectra)

lowrank_kinship <- try(new_lowrank_kinship(x = admixed$xtrain, y = admixed$ytrain,
                                       kinship = admixed$kin_train,
                                       n_nonzero_eigenvalues = number_nonzero_eigenvalues,
                                       n_zero_eigenvalues = nrow(admixed$kin_train) - number_nonzero_eigenvalues),
                       silent = TRUE)

lowrank_K <- try(new_lowrank_K(x = admixed$xtrain, y = admixed$ytrain,
                               K = admixed$Xkinship,
                               n_nonzero_eigenvalues = number_nonzero_eigenvalues,
                               n_zero_eigenvalues = min(nrow(admixed$Xkinship),
                                                        ncol(admixed$Xkinship)) -
                                                          number_nonzero_eigenvalues),
                 silent = TRUE)

lowrank_UD <- try(new_lowrank_UD(x = admixed$xtrain, y = admixed$ytrain,
                                 U = svdXkinship$u, D = svdXkinship$d,
                                 n_nonzero_eigenvalues = length(svdXkinship$d),
                                 n_zero_eigenvalues = 50),
                  silent = TRUE)



test_that("no error in low rank constructor functions", {

  expect_false(inherits(lowrank_kinship, "try-error"))
  expect_false(inherits(lowrank_K, "try-error"))
  expect_false(inherits(lowrank_UD, "try-error"))
  expect_is(lowrank_kinship, class = "lowrank")
  expect_is(lowrank_K, class = "lowrank")
  expect_is(lowrank_UD, class = "lowrank")

  # expect_equivalent(class(coef(fit_sim)), "dgCMatrix")
  # expect_equal(dim(coef(fit_sim))[[1]], dim(fit_sim$design)[[2]] + 1)
  # expect_equal(dim(coef(fit_sim))[[2]], sum(fit_sim$converged))

})


