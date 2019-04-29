library(ggmix)
library(glmnet)
library(gaston)
pacman::p_load(caret)
data("admixed")

# devtools::install_github('StoreyLab/popkin', build_opts=c())
library(popkin)

## ---- setup-data ----

set.seed(123)
ind <- caret::createDataPartition(admixed$y, p = 0.7, list = FALSE)[,1]

xtrain <- admixed$x[ind,,drop=FALSE]
xtest <- admixed$x[-ind,,drop=FALSE]

ytrain <- admixed$y[ind]
ytest <- admixed$y[-ind]

Xall <- rbind(xtest, xtrain)
Phi <- 2 * popkin::popkin(admixed$x, lociOnCols = TRUE)
cov_train <- Phi[ind,ind]
cov_test_train <- Phi[-ind,ind]



popkin::plotPopkin(cov_train)
popkin::plotPopkin(cov_test_train)

set.seed(123)
ind <- caret::createDataPartition(y, p = 0.7, list = FALSE)[,1]




# ggmix -------------------------------------------------------------------


fit_ggmix <- ggmix(x = xtrain, y = ytrain, kinship = cov_train, verbose = 1)
bicGGMIX <- gic(fit_ggmix, an = log(length(ytrain)))
yhat_test_ggmix_individual <- predict(bicGGMIX, s="lambda.min", newx = xtest, type = "individual", covariance = cov_test_train)
yhat_test_ggmix_population <- predict(bicGGMIX, s="lambda.min", newx = xtest, type = "response")
(RMSE_ggmix_test_individual <- l2norm(yhat_test_ggmix_individual - ytest))
(RMSE_ggmix_test_population <- l2norm(yhat_test_ggmix_population - ytest))



## ---- two-step ----

x1 <- cbind(rep(1, nrow(xtrain)))
fit_lme <- gaston::lmm.aireml(Y = ytrain, X = x1, K = cov_train)
gaston_resid <- ytrain - (fit_lme$BLUP_omega + fit_lme$BLUP_beta)
fit_twostep <- glmnet::cv.glmnet(x = xtrain, y = gaston_resid,
                                 standardize = T, alpha = 1, intercept = T)
yhat_twostep <- predict(fit_twostep, newx = xtest, s = "lambda.min")
(RMSE_twostep <- l2norm(yhat_twostep - ytest))


## ---- lasso ----

PC <- prcomp(xtrain)
xtrain_lasso <- cbind(xtrain, PC$x[,1:10])

fit_glmnet <- cv.glmnet(x = xtrain_lasso,
                        y = ytrain,
                        alpha = 1,
                        standardize = T,
                        penalty.factor = c(rep(1, ncol(xtrain)), rep(0,10)))

xtest_pc <- predict(PC, newdata = xtest)
xtest_lasso <- cbind(xtest, xtest_pc[,1:10])

yhat_lasso <- predict(fit_glmnet, s="lambda.min", newx = xtest_lasso)
(RMSE_lasso <- l2norm(yhat_lasso - ytest))


rbind(RMSE = c(ggmix = RMSE_ggmix_test_individual, twostep = RMSE_twostep, lasso = RMSE_lasso),
      active = c(ggmix = length(nonzeroCoef(coef(bicGGMIX, s = "lambda.min"))),
                 twostep = length(nonzeroCoef(coef(fitglmnet, s = "lambda.min"))),
                 lasso = length(nonzeroCoef(coef(fit_glmnet, s = "lambda.min")))))


