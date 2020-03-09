# test ranef

context(strwrap("check random effects"))

data("admixed")

test_that("Check random effects for multiple s values", {

  fit <- ggmix(x = admixed$xtrain, y = admixed$ytrain, kinship = admixed$kin_train,
               estimation = "full")
  gicfit <- gic(fit)

  expect_s3_class(fit, "ggmix_fit")
  expect_s3_class(gicfit, "ggmix_gic")
  expect_error(ranef(fit), regexp = "This function should be used with")
  expect_length(ranef(gicfit), length(admixed$ytrain))
  expect_equal(dim(ranef(gicfit, s = c(0.1,0.2))), c(length(admixed$ytrain), 2))

})


# Requires popkin package


test_that("Check predicted random effects on test set", {

  draw <- gen_structured_model(n = 500,
                                p_design = 200,
                                p_kinship = 1e4,
                                geography = "1d",
                                percent_causal = 0.05,
                                percent_overlap = "100",
                                k = 5, s = 0.5, Fst = 0.1, nPC = 10,
                                b0 = 0, eta = 0.1, sigma2 = 1)

  fit <- ggmix(x = draw[["xtrain"]],
               y = draw[["ytrain"]],
               kinship = draw[["kin_train"]])
  predmat <- predict(fit,
                     newx = draw[["xtune"]],
                     type = "individual",
                     covariance = draw[["kin_tune_train"]],
                     s = fit$lambda)
  cvmat <- apply((draw[["ytune"]] - predmat)^2, 2, mean)
  lambda_min_ggmix <- fit$result[which.min(cvmat), "Lambda"]

  # inidividual level prediction
  yhat <- predict(fit,
                  s = lambda_min_ggmix,
                  newx = draw[["xtest"]],
                  type = "individual",
                  covariance = draw[["kin_test_train"]])

  expect_length(yhat, n = nrow(draw[["xtest"]]))


})
