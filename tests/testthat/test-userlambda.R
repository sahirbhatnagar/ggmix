# test use defined lambda which is based on predict methods test

context(strwrap("check user defined lambda"))

data("admixed")

test_that("Check predict and coef methods with multiple s values", {
  skip_on_cran()
  fit <- ggmix(x = admixed$xtrain, y = admixed$ytrain, kinship = admixed$kin_train,
               estimation = "full", lambda = c(0.235, 0.55, 0.19, 0.25, 0.00235, 0.00003))
  gicfit <- gic(fit)
  
  expect_error(predict(gicfit), regexp = "You need to supply a value for")
  expect_length(drop(predict(gicfit, newx = admixed$xtest)), length(admixed$ytest))
  expect_equal(dim(predict(gicfit, newx = admixed$xtest, s = c(0.1,0.2))),
               c(length(admixed$ytest), 2))
  expect_equal(predict(gicfit, type = "coef"), coef(gicfit))
  expect_equal(nrow(predict(gicfit, type = "all")), ncol(admixed$xtrain) + 3)
  expect_is(predict(gicfit, type = "non"), class = "matrix")
  expect_is(predict(fit, type = "non", s = c(0.1,0.2)), class = "list")
  expect_s3_class(fit, "lassofullrank")
  expect_s3_class(gicfit, "ggmix_gic")
})
