library(ggmix)
devtools::load_all()
library(glmnet)
library(gaston)
data("admixed")

## ---- ggmix ----
fit_ggmix <- ggmix(x = admixed$x, y = admixed$y, kinship = admixed$kin, verbose = 1)
bicGGMIX <- gic(fit_ggmix, an = log(length(admixed$y)))
yhat_ggmix <- predict(bicGGMIX, newx = admixed$x)
RMSE_ggmix <- l2norm(yhat_ggmix - admixed$y)



## ---- two-step ----
x1 <- cbind(rep(1, nrow(admixed$x)))
fit_lme <- gaston::lmm.aireml(Y = admixed$y, X = x1, K = admixed$kin)
gaston_resid <- admixed$y - (fit_lme$BLUP_omega + fit_lme$BLUP_beta)
fitglmnet <- glmnet::cv.glmnet(x = admixed$x, y = gaston_resid,
                               standardize = T, alpha = 1, intercept = T)
yhat_twostep <- predict(fitglmnet, newx = admixed$x, s = "lambda.min")
RMSE_twostep <- l2norm(yhat_twostep - admixed$y)

## ---- lasso ----
fit_glmnet <- cv.glmnet(x = admixed$x_lasso, y = admixed$y,
                        alpha = 1, standardize = T,
                        penalty.factor = c(rep(1, ncol(admixed$x)), rep(0,10)))
# extract only betas for SNPs
betas_lasso <- coef(fit_glmnet, s = "lambda.min")[1:(ncol(admixed$x)+1), , drop = F]
yhat_lasso <- cbind(1, admixed$x) %*% betas_lasso
RMSE_lasso <- l2norm(yhat_lasso - admixed$y)
