# test ranef

context(strwrap("check random effects"))

data("admixed")

test_that("Check random effects for multiple s values", {

  fit <- ggmix(x = admixed$x, y = admixed$y, kinship = admixed$kin,
               estimation = "full")
  gicfit <- gic(fit)

  expect_s3_class(fit, "ggmix_fit")
  expect_s3_class(gicfit, "ggmix_gic")
  expect_error(ranef(fit), regexp = "This function should be used with")
  expect_length(ranef(gicfit), length(admixed$y))
  expect_equal(dim(ranef(gicfit, s = c(0.1,0.2))), c(length(admixed$y), 2))

})
