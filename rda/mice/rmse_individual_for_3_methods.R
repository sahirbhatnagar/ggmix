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



## ---- setupe-train/test/validate-split ----


spec = c(train = .6, test = .3, validate = .1)

g = sample(cut(
  seq(nrow(admixed$x)),
  nrow(admixed$x)*cumsum(c(0,spec)),
  labels = names(spec)
))
g %>% table

train_ind <- which(g == "train")
validate_ind <- which(g == "validate")
test_ind <- which(g == "test")
# res = split(admixed$x, g)

xtrain <- admixed$x[train_ind,,drop=FALSE]
xtest <- admixed$x[test_ind,,drop=FALSE]
xvalidate <- admixed$x[validate_ind,,drop=FALSE]

ytrain <- admixed$y[train_ind]
ytest <- admixed$y[test_ind]
yvalidate <- admixed$y[validate_ind]

help(admixed)
cov_train <- admixed$kin[train_ind,train_ind]
cov_test_train <- admixed$kin[test_ind,train_ind]
cov_validate_train <- admixed$kin[validate_ind,train_ind]

cov_train %>% dim
cov_test_train %>% dim
cov_validate_train %>% dim

# ggmix -------------------------------------------------------------------


fit_ggmix <- ggmix(x = xtrain, y = ytrain, kinship = cov_train, verbose = 1)
predmat <- predict(fit_ggmix, newx = xtest, type = "individual", covariance = cov_test_train, s = fit_ggmix$lambda)
cvmat <- apply((ytest - predmat)^2, 2, mean)

dev.off()
bicGGMIX <- gic(fit_ggmix, an = log(length(ytrain)))
par(mfrow=c(1,2))
plot(bicGGMIX)
plot(log(fit_ggmix$lambda), cvmat, pch = 19)
lambda.min.name = which.min(cvmat)
lambda.min = fit_ggmix$result[lambda.min.name, "Lambda"]
bicGGMIX$lambda.min

sum(admixed$causal %in% rownames(predict(fit_ggmix, s = lambda.min, type = "coef")[nonzeroCoef(predict(fit_ggmix, s = lambda.min, type = "coef")),,drop=F])) / length(admixed$causal)
rownames(predict(fit_ggmix, s = lambda.min, type = "coef")[nonzeroCoef(predict(fit_ggmix, s = lambda.min, type = "coef")),,drop=F]) %>% length()
sum(admixed$causal %in% rownames(predict(bicGGMIX, s = "lambda.min", type = "coef")[nonzeroCoef(predict(bicGGMIX, s = "lambda.min", type = "coef")),,drop=F])) / length(admixed$causal)
rownames(predict(bicGGMIX, s = "lambda.min", type = "coef")[nonzeroCoef(predict(bicGGMIX, s = "lambda.min", type = "coef")),,drop=F]) %>% length()


valdpred <- predict(fit_ggmix, newx = xvalidate, type = "individual", covariance = cov_validate_train, s = lambda.min)
l2norm(yvalidate-valdpred)
valdpred <- predict(fit_ggmix, newx = xvalidate, type = "individual", covariance = cov_validate_train, s = bicGGMIX$lambda.min)
l2norm(yvalidate-valdpred)

plot(yvalidate, valdpred, pch = 19)
abline(a=0, b=1)
cor(yvalidate, valdpred)^2








coef()
admixed$causal

bicGGMIX$lambda.min

identical(predict(fit_ggmix, newx = xtest, type = "individual", covariance = cov_test_train, s = bicGGMIX$lambda),
          predict(bicGGMIX, newx = xtest, type = "individual", covariance = cov_test_train, s = bicGGMIX$lambda))

ggmix:::predict.ggmix_fit()
help(ggmix:::predict.ggmix_gic)
yhat_test_ggmix_individual <- predict(bicGGMIX, s="lambda.min", newx = xtest, type = "individual", covariance = cov_test_train)


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

fit_glmnet <- glmnet(x = xtrain_lasso,
                     y = ytrain,
                     alpha = 1,
                     standardize = T,
                     penalty.factor = c(rep(1, ncol(xtrain)), rep(0,10)))

xtest_pc <- predict(PC, newdata = xtest)
xtest_lasso <- cbind(xtest, xtest_pc[,1:10])

yhat_lasso <- predict(fit_glmnet, newx = xtest_lasso)
cvmat <- apply((ytest - yhat_lasso)^2, 2, mean)
dev.off()
plot(log(fit_glmnet$lambda), cvmat, pch=19)

xvalidate_pc <- predict(PC, newdata = xvalidate)
xvalidate_lasso <- cbind(xvalidate, xvalidate_pc[,1:10])


valpred <- predict(fit_glmnet, newx = xvalidate_lasso, s = fit_glmnet$lambda[which.min(cvmat)])
(RMSE_lasso <- l2norm(valpred - yvalidate))

nz_coef_lasso <- coef(fit_glmnet, s = fit_glmnet$lambda[which.min(cvmat)])[nonzeroCoef(coef(fit_glmnet, s = fit_glmnet$lambda[which.min(cvmat)])),,drop=F]

setdiff(rownames(nz_coef_lasso), c("(Intercept)",paste0("PC",1:10))) %>% length()

sum(setdiff(rownames(nz_coef_lasso), c("(Intercept)",paste0("PC",1:10))) %in% admixed$causal)

rbind(RMSE = c(ggmix = RMSE_ggmix_test_individual, twostep = RMSE_twostep, lasso = RMSE_lasso),
      active = c(ggmix = length(nonzeroCoef(coef(bicGGMIX, s = "lambda.min"))),
                 twostep = length(nonzeroCoef(coef(fitglmnet, s = "lambda.min"))),
                 lasso = length(nonzeroCoef(coef(fit_glmnet, s = "lambda.min")))))


