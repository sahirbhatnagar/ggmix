# library(ggmix)
devtools::load_all()
devtools::document()
library(glmnet)
library(gaston)
pacman::p_load(caret)
data("admixed")

# devtools::install_github('StoreyLab/popkin', build_opts=c())
library(popkin)

## ---- ggmix ----

set.seed(1234)
ind <- caret::createDataPartition(admixed$y, p = 0.8, list = FALSE)[,1]
xtrain <- admixed$x[ind,,drop=FALSE]
xtest <- admixed$x[-ind,,drop=FALSE]

ytrain <- admixed$y[ind]
ytest <- admixed$y[-ind]

Xall <- rbind(xtest, xtrain)
cov_train <- 2 * popkin::popkin(xtrain, lociOnCols = TRUE)
dim(cov_train)

cov_all <- 2 * popkin::popkin(Xall, lociOnCols = TRUE)
dim(cov_all)

cov_test_train <- cov_all[1:nrow(xtest), (nrow(xtest)+1):ncol(cov_all)]

dim(cov_test_train)


fit_ggmix <- ggmix(x = xtrain, y = ytrain, kinship = cov_train, verbose = 1)
bicGGMIX <- gic(fit_ggmix, an = log(length(ytrain)))
plot(bicGGMIX)
coef(bicGGMIX, s = "lambda.min")
yhat_test <- predict(bicGGMIX, s="lambda.min", newx = xtest, type = "individual", covariance = cov_test_train)
yhat_test_population <- predict(bicGGMIX, s="lambda.min", newx = xtest, type = "response")
(RMSE_ggmix_test <- l2norm(yhat_test - ytest))
(RMSE_ggmix_test_pop <- l2norm(yhat_test_population - ytest))

cor(yhat_test, ytest)^2
plot(yhat_test, ytest)
abline(a=0,b=1)
plot(yhat_test_population, yhat_test)
plot(yhat_test_population, ytest)
sd(yhat_test)
sd(yhat_test_population)
abline(a=0,b=1)
cor(yhat_test_population, ytest)^2
# yhat_train <- predict(bicGGMIX, newx = xtrain) + ranef(bicGGMIX)
# RMSE_ggmix_train <- l2norm(yhat_train - ytrain)



## ---- two-step ----
x1 <- cbind(rep(1, nrow(xtrain)))
fit_lme <- gaston::lmm.aireml(Y = ytrain, X = x1, K = cov_train)
gaston_resid <- ytrain - (fit_lme$BLUP_omega + fit_lme$BLUP_beta)
fitglmnet <- glmnet::cv.glmnet(x = xtrain, y = gaston_resid,
                               standardize = T, alpha = 1, intercept = T)
yhat_twostep <- predict(fitglmnet, newx = xtest, s = "lambda.min")
(RMSE_twostep <- l2norm(yhat_twostep - ytest))

cor(yhat_twostep, ytest)^2

## ---- lasso ----
# eiK <- eigen(cov_train)

PC <- prcomp(xtrain)
PC$rotation %>% dim
PC$x %>% dim


# if (any(eiK$values < 1e-5)) {
#   eiK$values[ eiK$values < 1e-5 ] <- 1e-5
# }
#
# PC <- sweep(eiK$vectors, 2, sqrt(eiK$values), "*")

xtrain_lasso <- cbind(xtrain, PC$x[,1:10])
# xtrain_lasso2 <- cbind(xtrain, prc$x[,1:10])

# eiK$vectors %>% dim

fit_glmnet <- cv.glmnet(x = xtrain_lasso,
                        y = ytrain,
                        alpha = 1,
                        standardize = T,
                        penalty.factor = c(rep(1, ncol(xtrain)), rep(0,10)))
plot(fit_glmnet)

xtest %>% dim
xtest_pc <- predict(PC, newdata = xtest)
xtest_pc %>% dim
xtest_lasso <- cbind(xtest, xtest_pc[,1:10])


# extract only betas for SNPs
# betas_lasso <- coef(fit_glmnet, s = "lambda.min")[1:(ncol(xtrain)+1), , drop = F]
# yhat_lasso <- cbind(1, xtest) %*% betas_lasso

yhat_lasso <- predict(fit_glmnet, s="lambda.min", newx = xtest_lasso)
(RMSE_lasso <- l2norm(yhat_lasso - ytest))
cor(as.vector(yhat_lasso), ytest)^2


rbind(RMSE = c(ggmix = RMSE_ggmix_test, twostep = RMSE_twostep, lasso = RMSE_lasso),
active = c(ggmix = length(nonzeroCoef(coef(bicGGMIX, s = "lambda.min"))),
  twostep = length(nonzeroCoef(coef(fitglmnet, s = "lambda.min"))),
  lasso = length(nonzeroCoef(coef(fit_glmnet, s = "lambda.min")))))


