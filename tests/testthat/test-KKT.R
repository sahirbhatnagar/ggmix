# test KKT

context(strwrap("KKT checks"))

data("admixed")

test_that("Check predict and coef methods with multiple s values", {
  fit <- ggmix(x = admixed$x, y = admixed$y, kinship = admixed$kin,
               estimation = "full")
  gicfit <- gic(fit)

  coef_res <- coef(gicfit, type = "all")
  et <- coef_res["eta",]
  sigs <- coef_res["sigma2",]
  bet <- coef_res[-which(rownames(coef_res) %in% c("eta","sigma2")),]
  kkt <- kkt_check(eta = et, sigma2 = sigs, beta = bet,
                   eigenvalues = fit$ggmix_object$D, x = fit$ggmix_object$x,
                   y = fit$ggmix_object$y, nt = length(fit$ggmix_object$y),
                   lambda = gicfit$lambda.min, tol.kkt = 1e-3)
  expect_true(all(abs(kkt) < 0.02))

})
