# options(warnPartialMatchArgs = FALSE, warnPartialMatchDollar = TRUE, warnPartialMatchAttr = TRUE)
library(magrittr)
library(bit64)
library(data.table)
library(glmnet)
library(MASS)
# library(lmmlasso)

# Simulate Data -----------------------------------------------------------


rm(list = ls())
source("~/git_repositories/penfam/R/functions.R")
source("~/git_repositories/penfam/R/sim-data.R")



# lambda <- lambda_sequence(x = X, y = Y, phi = Phi, lambda_min_ratio = 1e-4)

res <- penfam(x = X, y = Y, phi = Phi, lambda_min_ratio = 0.001)
res$lambda_min
coef(res, s=res$lambda_min)
predict(res, type = "nonzero", s = res$lambda_min)
coef(res)

plot(res, xvar = "norm")
plot(res, xvar = "lambda")
plot(res, xvar = "dev")

plot(res, type = "BIC", sign.lambda = 1)



plot(coef(res, s = res$lambda_min), pch = 19, col = "red")
points(seq_along(c(b0,b, eta, sigma2)), c(b0,b, eta, sigma2), pch = 19, col = "blue")
legend("bottomleft",
       legend = c("Estimated", "Truth"),
       col = c("red","blue"),
       pch = c(19, 19),
       bg = "gray90")





# lambda sequence ---------------------------------------------------------

rm(list = ls())
library(microbenchmark)
source("~/git_repositories/penfam/R/functions.R")
source("~/git_repositories/penfam/R/sim-data.R")









# grplasso ----------------------------------------------------------------

pacman::p_load(grplasso)
View(grplasso:::grplasso.default)


# glmnet ------------------------------------------------------------------

fit <- glmnet::cv.glmnet(x = x, y = y)
plot(fit)
coef(fit)[(coef(fit)==0) , ]

beta_hats <- as.matrix(coef(fit, s = "lambda.1se"))
s0 <- beta_hats[beta_hats!=0,, drop=F]

# true positive rate
sum(paste0("V",1:20) %in% rownames(s0))/20

# true negative rate
sum(paste0("V",21:1000) %ni% rownames(s0))/980



# lmmlasso ----------------------------------------------------------------

data(classroomStudy)
head(classroomStudy)
fit1 <- lmmlasso(x=classroomStudy$X,y=classroomStudy$y,z=classroomStudy$Z,
                 grp=classroomStudy$grp,lambda=15,pdMat="pdIdent")
summary(fit1)
plot(fit1)


Xmat <- model.matrix(as.formula(paste("~", paste(colnames(x), collapse = "+"))), data = as.data.frame(x))
Xmat[1:5,1:5]
x[1:5,1:5]
fit1 <- lmmlasso(x = Xmat, y = y, z = as.matrix(rep(1, 600)), grp = 1:600, lambda = 29, pdMat = "pdIdent")
summary(fit1)
plot(fit1)

beta_hats <- as.matrix(fit1$coefficients)
s0 <- beta_hats[beta_hats!=0,, drop=F]

# true positive rate
sum(paste0("V",1:20) %in% rownames(s0))/20

# true negative rate
sum(paste0("V",21:1000) %ni% rownames(s0))/980




# quadform ----------------------------------------------------------------

pacman::p_load(emulator)
pacman::p_load(microbenchmark)

jj <- matrix(rnorm(100*1000),ncol = 1000)
M <- crossprod(jj,jj)
M.lower <- t(chol(M))
x <- matrix(rnorm(8),4,2)

jj.1 <- t(x) %*% M %*% x
jj.2 <- quad.form(M,x)
jj.3 <- quad.form(M.lower,x,chol=TRUE)
print(jj.1)
print(jj.2)
print(jj.3)




emulator::quad.form
cprod()
