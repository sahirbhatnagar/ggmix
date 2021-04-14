# test KKT

context(strwrap("KKT checks"))
set.seed(1234)
draw <- gen_structured_model(n = 500,
                             p_design = 100,
                             p_kinship = 1e4,
                             geography = "1d",
                             percent_causal = 0.05,
                             percent_overlap = "100",
                             k = 5, s = 0.5, Fst = 0.1, nPC = 10,
                             b0 = 0, eta = 0.1, sigma2 = 1)

fit <- ggmix(x = draw[["xtrain"]],
             y = draw[["ytrain"]],
             kinship = draw[["kin_train"]],
             estimation = "full", epsilon = 1e-7, verbose = 0)

test_that("Check predict and coef methods with multiple s values", {
  # fit <- ggmix(x = admixed$xtrain, y = admixed$ytrain, kinship = admixed$kin_train,
  #              estimation = "full")
  skip_on_cran()
  gicfit <- gic(fit)

  coef_res <- coef(gicfit, type = "all")
  et <- coef_res["eta",]
  sigs <- coef_res["sigma2",]
  bet <- coef_res[-which(rownames(coef_res) %in% c("eta","sigma2")),]
  kkt <- kkt_check(eta = et, sigma2 = sigs, beta = bet,
                   eigenvalues = fit$ggmix_object$D, x = fit$ggmix_object$x,
                   y = fit$ggmix_object$y, nt = length(fit$ggmix_object$y),
                   lambda = gicfit$lambda.min, tol.kkt = 1e-2)
  expect_true(all(abs(kkt)[-1] < 0.02))

})



fit <- ggmix(x = draw[["xtrain"]],
             y = draw[["ytrain"]],
             kinship = draw[["kin_train"]],
             estimation = "full", epsilon = 1e-10, 
             verbose = 0, lambda = c(0.22150, 0.025, 0.01, 0.33))

test_that("Check predict and coef methods with multiple s values and user defined lambda", {
  # fit <- ggmix(x = admixed$xtrain, y = admixed$ytrain, kinship = admixed$kin_train,
  #              estimation = "full")
  skip_on_cran()
  gicfit <- gic(fit)
  
  coef_res <- coef(gicfit, type = "all")
  et <- coef_res["eta",]
  sigs <- coef_res["sigma2",]
  bet <- coef_res[-which(rownames(coef_res) %in% c("eta","sigma2")),]
  kkt <- kkt_check(eta = et, sigma2 = sigs, beta = bet,
                   eigenvalues = fit$ggmix_object$D, x = fit$ggmix_object$x,
                   y = fit$ggmix_object$y, nt = length(fit$ggmix_object$y),
                   lambda = gicfit$lambda.min, tol.kkt = 1e-2)
  expect_true(all(abs(kkt)[-1] < 0.02))
  
})
