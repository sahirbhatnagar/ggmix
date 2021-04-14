# test U,D, kinship, K, inputs (see ggmix function)

context(strwrap("check different inputs for the relationship matrix U, D
                K and kinship when passed to ggmix"))

data("admixed")

svdXkinship <- svd(admixed$Xkinship)
number_nonzero_eigenvalues <- 20



test_that("Expect error when user has not provided number of (non)zero
          eigenvalues to ggmix", {
  skip_on_cran()
  # use these when lowrank has been implemented
  # expect_error(ggmix(x = admixed$x, y = admixed$y, kinship = admixed$kin,
  #                    estimation = "low"),
  #              regexp = strwrap("when kinship is specified and estimation=\"low\",
  #                    n_nonzero_eigenvalues must be specified."))
  #
  # expect_error(ggmix(x = admixed$x, y = admixed$y, K = admixed$Xkinship,
  #              estimation = "low"),
  #              regexp = strwrap("when K is specified and estimation=\"low\",
  #                    n_nonzero_eigenvalues must be specified."))
  #
  # expect_error(ggmix(x = admixed$x, y = admixed$y, U = svdXkinship$u,
  #                    D = svdXkinship$d,
  #                    estimation = "low"),
  #              regexp = strwrap("n_zero_eigenvalues must be specified when U and D have
  #                  been provided and estimation=\"low\""))
  #
  # expect_warning(ggmix(x = admixed$x, y = admixed$y, K = admixed$Xkinship,
  #                    estimation = "low", n_nonzero_eigenvalues = number_nonzero_eigenvalues),
  #              regexp = sprintf("n_zero_eigenvalues missing. setting to %g",
  #                               min(nrow(admixed$Xkinship), ncol(admixed$Xkinship)) -
  #                                 number_nonzero_eigenvalues))
  #
  # expect_warning(ggmix(x = admixed$x, y = admixed$y, kinship = admixed$kin,
  #                      estimation = "low", n_nonzero_eigenvalues = number_nonzero_eigenvalues),
  #                regexp = sprintf("n_zero_eigenvalues missing. setting to %g",
  #                                 nrow(admixed$kin) -
  #                                   number_nonzero_eigenvalues))

  expect_warning(ggmix(x = admixed$xtrain, y = admixed$ytrain, kinship = admixed$kin_train,
                       estimation = "full",
                       U = svdXkinship$u[,which(svdXkinship$d>0),drop=FALSE], D = svdXkinship$d[which(svdXkinship$d>0)],
                       K = admixed$Xkinship),
                 regexp = strwrap("kinship, U, D and K arguments have all been specified."))

})

context(strwrap("check simulation function with different geography"))

test_that("test simulation function", {
  skip_on_cran()
  set.seed(2345)
  sim1d <- gen_structured_model(n = 100,
                                p_design = 50,
                                p_kinship = 5e2,
                                geography = "1d",
                                percent_causal = 0.10,
                                percent_overlap = "100",
                                k = 5, s = 0.5, Fst = 0.1, nPC = 5,
                                b0 = 0, eta = 0.1, sigma2 = 1)
  expect_length(sim1d, 21)
  expect_is(sim1d, "list")

  simcirc <- gen_structured_model(n = 100,
                                  p_design = 50,
                                  p_kinship = 5e2,
                                  geography = "circ",
                                  percent_causal = 0.10,
                                  percent_overlap = "100",
                                  k = 5, s = 0.5, Fst = 0.1, nPC = 5,
                                  b0 = 0, eta = 0.1, sigma2 = 1)
  expect_length(simcirc, 21)
  expect_is(simcirc, "list")

  simind <- gen_structured_model(n = 100,
                                 p_design = 50,
                                 p_kinship = 5e2,
                                 geography = "ind",
                                 percent_causal = 0.10,
                                 percent_overlap = "100",
                                 k = 5, s = 0.5, Fst = 0.1, nPC = 10,
                                 b0 = 0, eta = 0.1, sigma2 = 1)
  expect_length(simind, 21)
  expect_is(simind, "list")
})



