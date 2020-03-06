# test predict methods

context(strwrap("check predict and coef methods"))

data("admixed")

test_that("Check predict and coef methods with multiple s values", {

  fit <- ggmix(x = admixed$xtrain, y = admixed$ytrain, kinship = admixed$kin_train,
               estimation = "full")
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
