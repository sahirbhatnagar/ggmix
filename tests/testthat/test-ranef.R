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


# Requires popkin package


test_that("Check predicted random effects on test set", {

  ind <- sample(1:nrow(admixed$x), size = floor(nrow(admixed$x)*0.80), replace = FALSE)
  xtrain <- admixed$x[ind,,drop=FALSE]
  xtest <- admixed$x[-ind,,drop=FALSE]

  ytrain <- admixed$y[ind]
  ytest <- admixed$y[-ind]

  Xall <- rbind(xtest, xtrain)
  cov_train <- 2 * popkin::popkin(xtrain, lociOnCols = TRUE)

  cov_all <- 2 * popkin::popkin(Xall, lociOnCols = TRUE)

  cov_test_train <- cov_all[1:nrow(xtest), (nrow(xtest)+1):ncol(cov_all)]

  fit_ggmix <- ggmix(x = xtrain, y = ytrain, kinship = cov_train, verbose = 1)
  bicGGMIX <- gic(fit_ggmix, an = log(length(ytrain)))
  yhat_test <- predict(bicGGMIX, s="lambda.min", newx = xtest, type = "individual", covariance = cov_test_train)

  expect_length(yhat_test, n = nrow(xtest))


})
