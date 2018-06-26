# test U,D, kinship, K, inputs (see ggmix function)

context(strwrap("check different inputs for the relationship matrix U, D
                K and kinship when passed to ggmix"))

data("admixed")

svdXkinship <- svd(admixed$Xkinship)
number_nonzero_eigenvalues <- 20



test_that("Expect error when user has not provided number of (non)zero
          eigenvalues to ggmix", {

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

  expect_warning(ggmix(x = admixed$x, y = admixed$y, kinship = admixed$kin,
                       estimation = "full",
                       U = svdXkinship$u, D = svdXkinship$d,
                       K = admixed$Xkinship),
                 regexp = strwrap("kinship, U, D and K arguments have all been specified."))

})
