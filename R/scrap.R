#'
#' # Trying low rank ---------------------------------------------------------
#'
#' rm(list = ls())
#' pacman::p_load(gaston)
#' pacman::p_load(glmnet)
#' pacman::p_load(magrittr)
#' pacman::p_load(snpStats)
#' source("~/git_repositories/penfam/R/fitting.R")
#' source("~/git_repositories/penfam/R/functions.R")
#' source("~/git_repositories/penfam/R/methods.R")
#' source("~/git_repositories/penfam/R/plot.R")
#' source("~/git_repositories/penfam/simulation/model_functions.R")
#' dat <- make_mixed_model_not_simulator(b0 = 1, eta = 0.3, sigma2 = 2, type = "causal_400", related = TRUE)
#' w_svd <- svd(dat$w)
#' U <- w_svd$u
#'
#' # we dived by p-1 because thats how the matrix was standardized
#' Lambda <- w_svd$d^2 / (ncol(dat$x)-1)
#' any(Lambda<1e-5)
#' Lambda[Lambda<1e-5] <- 1e-5
#'
#' res <- lowrank(x = dat$x[,1:500], y = dat$y,  d = Lambda, u = U)
#'
#' res <- gic.penfam(x = dat$x[,1:500], y = dat$y,  d = Lambda, u = U)
#'
#' res
#'
#' plot(res)
#' coef(res)
#' plot(res$penfam.fit, type="coef")
#'
#' t(res$coef)
#'
#' bic <- function(eta, sigma2, beta, eigenvalues, x, y, nt, c, df_lambda) {
#'
#'   -2 * log_lik(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt) + c * df_lambda
#'
#' }
#'
#'
#'
#'
#' bic(eta = res$eta, sigma2 = res$sigma2)
#'
#' lasso <- cv.glmnet(x = dat$x, y = dat$y, alpha = 1)
#' fitted(lasso)
#' plot(lasso)
#'
#' coef(lasso, s="lambda.min")[nonzeroCoef(coef(lasso, s="lambda.min")),,drop=F]
#' grep("rs", rownames(coef(lasso)[nonzeroCoef(coef(lasso)),,drop=F]), value = T) %in% dat$causal
#' grep("rs", rownames(coef(res)[nonzeroCoef(coef(res)),,drop=F]), value = T) %in% dat$causal
#'
#' coef(lasso)[nonzeroCoef(coef(lasso)),,drop=F]
#' coef(res)[nonzeroCoef(coef(res)),,drop=F]
#'
#' dat$beta[dat$beta!=0]
#'
#' plot(coef(lass0)[-1], dat$beta)
#' abline(a=0, b=1)
#' all(colnames(dat$x)==rownames(coef(lass0)[-1]))
#'
#' coef(res, s = res$lambda_min) %>% head
#' plot(coef(resenet, s = resenet$lambda_min)[c(-1,-4002,-4003),,drop=T], dat$beta)
#' abline(a=0, b=1)
#' all(colnames(dat$x)==rownames(coef(lass0)[-1]))
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' # this needs work, as lambda.min is different
#' resenet <- lowrank(x = dat$x, y = dat$y, w = dat$w, alpha = 1)
#'
#' dim(dat$x)
#'
#' par(mfrow=c(3,2))
#' plot(res)
#' # c("coef","BIC", "QQranef","QQresid", "predicted", "Tukey-Anscombe")
#' plot(res, type = "BIC")
#' plot(res, type = "QQranef")
#' plot(res, type = "QQresid")
#' plot(res, type = "predicted")
#' plot(res, type = "Tukey")
#'
#' plot(resenet, type = "BIC")
#' predict(res, type = "nonzero", s = resenet$lambda_min)
#'
#' grep("rs", rownames(predict(res, type = "nonzero", s = res$lambda_min)), value = T) %in% dat$causal
#' grep("rs", rownames(predict(resenet, type = "nonzero", s = resenet$lambda_min)), value = T) %in% dat$causal
#'
#' plot(resenet, type = "QQranef")
#' plot(resenet, type = "QQresid")
#' plot(resenet, type = "predicted")
#' plot(resenet, type = "Tukey")
#'
#'
#' plot(res, type = "")
#' dev.off()
#' predict(res, type = "nonzero", s = res$lambda_min)
#' dat$causal
#' res$eta
#' res$sigma2
#'
#' lass0 <- cv.glmnet(x = dat$x, y = dat$y, alpha = 0.5)
#' plot(lass0)
#'
#' coef(lass0)[nonzeroCoef(coef(lass0)),,drop=F]
#' grep("rs", rownames(coef(lass0)[nonzeroCoef(coef(lass0)),,drop=F]), value = T) %in% dat$causal
#'
#' dat$beta[dat$beta!=0]
#'
#' plot(coef(lass0)[-1], dat$beta)
#' abline(a=0, b=1)
#' all(colnames(dat$x)==rownames(coef(lass0)[-1]))
#'
#' coef(res, s = res$lambda_min) %>% head
#' plot(coef(resenet, s = resenet$lambda_min)[c(-1,-4002,-4003),,drop=T], dat$beta)
#' abline(a=0, b=1)
#' all(colnames(dat$x)==rownames(coef(lass0)[-1]))
#'
#'
#' # options(warnPartialMatchArgs = FALSE, warnPartialMatchDollar = TRUE, warnPartialMatchAttr = TRUE)
#' library(magrittr)
#' library(glmnet)
#' library(MASS)
#'
#' rm(list=ls())
#' source("~/git_repositories/penfam/R/fitting.R")
#' source("~/git_repositories/penfam/R/functions.R")
#' source("~/git_repositories/penfam/R/methods.R")
#' source("~/git_repositories/penfam/R/plot.R")
#' source("~/git_repositories/penfam/R/sim-data.R")
#'
#' x <- X
#' y <- Y
#' phi <- Phi
#' lambda_min_ratio <- 0.001
#' nlambda <- 100
#' #convergence criterion
#' epsilon <- 1e-7
#' maxit <- 1e6
#' an = log(log(600)) * log(600)
#' lambda <- 0.10
#' #======================================
#'
#' np <- dim(x)
#' n <- np[[1]]
#' p <- np[[2]]
#'
#' # add column of 1s to x for intercept
#' x <- cbind(rep(1, n), x)
#' # x[1:5, 1:5]
#'
#' phi_eigen <- eigen(phi)
#' # this is a N_T x N_T matrix
#' U <- phi_eigen$vectors
#'
#' # dim(U)
#' # vector of length N_T
#' Lambda <- phi_eigen$values
#'
#' utx <- crossprod(U, x)
#' uty <- crossprod(U, y)
#'
#' # get sequence of tuning parameters
#' lamb <- lambda_sequence(x = utx, y = uty, phi = phi, lambda_min_ratio = lambda_min_ratio)
#'
#' res <- penfam(x = X, y = Y, phi = Phi, lambda_min_ratio = 0.001,
#'               nlambda = 100,
#'               an = log(log(n)) * log(n))
#'
#'
#' #' Check of KKT
#' #' @param x should be U^T X, where U is the matrix of eigenvectors and X
#' #'   contains the first column of ones for the intercept. x should be a mtrix of
#' #'   dimension n x (p+1)
#' #' @param beta should include intercept as well. A vector of length p+1.
#' #' @param tol.beta Tolerance for determining if a coefficient is zero
#' #' @param tol.kkt Tolerance for determining if an entry of the subgradient is
#' #'   zero
#'
#' kkt_check <- function(eta, sigma2, beta, eigenvalues, x, y, nt,
#'                       lambda, tol.beta = 1e-5, tol.kkt = 0.1){
#'
#'   # x = utx
#'   # y = uty
#'   # tol.beta = 1e-5 ; tol.kkt = 0.1
#'   # eigenvalues = Lambda
#'   # nt = dim(x)[[1]]
#'   # p = dim(x)[[2]]-1
#'   # eta = lamb$eta
#'   # sigma2 = lamb$sigma2
#'   # beta = c(lamb$beta0, rep(0, p))
#'   # =======================
#'
#'   di <- 1 + eta * (eigenvalues - 1)
#'   wi <- (1 / sigma2) * diag(1 / di)
#'
#'   # KKT for beta0
#'   kkt_beta0 <- grr_beta0(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt)
#'
#'   # KKT for eta
#'   kkt_eta <- grr_eta(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt)
#'
#'   # KKT for sigma2
#'   kkt_sigma2 <- grr_sigma2(eta = eta, sigma2 = sigma2, beta = beta, eigenvalues = eigenvalues, x = x, y = y, nt = nt)
#'
#'   # KKT for beta
#'   g0 <- crossprod(x[,-1, drop = F], wi) %*% (y - x %*% beta)
#'
#'   g <- g0 - lambda * sign(beta)
#'
#'   gg <- g0/lambda
#'
#'   # which of the betas are non-zero, subject to the tolerance level for beta
#'   oo <- abs(bhatt) > tol.beta
#'
#'   # cat(
#'   #   c(
#'   #     max(abs(g[oo])) > tol.kkt,
#'   #     min(gg[!oo]) < -1 - tol.kkt,
#'   #     max(gg[!oo]) > 1 + tol.kkt
#'   #   ),
#'   #   fill = T)
#'
#'   kkt_beta <- c(max(abs(g[oo])) > tol.kkt,
#'                 min(gg[!oo]) < -1 - tol.kkt,
#'                 max(gg[!oo]) > 1 + tol.kkt)
#'
#'   return(list(kkt_beta0 = kkt_beta0, kkt_eta = kkt_eta,
#'               kkt_sigma2 = kkt_sigma2, kkt_beta = kkt_beta))
#'
#'
#' }
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' # Simulate Data -----------------------------------------------------------
#'
#'
#' rm(list = ls())
#' source("~/git_repositories/penfam/R/functions.R")
#' source("~/git_repositories/penfam/R/sim-data.R")
#'
#'
#'
#' # lambda <- lambda_sequence(x = X, y = Y, phi = Phi, lambda_min_ratio = 1e-4)
#'
#' res <- penfam(x = X, y = Y, phi = Phi, lambda_min_ratio = 0.001)
#' res$lambda_min
#' coef(res, s=res$lambda_min)
#' predict(res, type = "nonzero", s = res$lambda_min)
#' coef(res)
#'
#' plot(res, xvar = "norm")
#' plot(res, xvar = "lambda")
#' plot(res, xvar = "dev")
#'
#' plot(res, type = "BIC", sign.lambda = 1)
#'
#'
#'
#' plot(coef(res, s = res$lambda_min), pch = 19, col = "red")
#' points(seq_along(c(b0,b, eta, sigma2)), c(b0,b, eta, sigma2), pch = 19, col = "blue")
#' legend("bottomleft",
#'        legend = c("Estimated", "Truth"),
#'        col = c("red","blue"),
#'        pch = c(19, 19),
#'        bg = "gray90")
#'
#'
#'
#'
#'
#' # lambda sequence ---------------------------------------------------------
#'
#' rm(list = ls())
#' library(microbenchmark)
#' source("~/git_repositories/penfam/R/functions.R")
#' source("~/git_repositories/penfam/R/sim-data.R")
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' # grplasso ----------------------------------------------------------------
#'
#' pacman::p_load(grplasso)
#' View(grplasso:::grplasso.default)
#'
#'
#' # glmnet ------------------------------------------------------------------
#'
#' fit <- glmnet::cv.glmnet(x = x, y = y)
#' plot(fit)
#' coef(fit)[(coef(fit)==0) , ]
#'
#' beta_hats <- as.matrix(coef(fit, s = "lambda.1se"))
#' s0 <- beta_hats[beta_hats!=0,, drop=F]
#'
#' # true positive rate
#' sum(paste0("V",1:20) %in% rownames(s0))/20
#'
#' # true negative rate
#' sum(paste0("V",21:1000) %ni% rownames(s0))/980
#'
#'
#'
#' # lmmlasso ----------------------------------------------------------------
#'
#' data(classroomStudy)
#' head(classroomStudy)
#' fit1 <- lmmlasso(x=classroomStudy$X,y=classroomStudy$y,z=classroomStudy$Z,
#'                  grp=classroomStudy$grp,lambda=15,pdMat="pdIdent")
#' summary(fit1)
#' plot(fit1)
#'
#'
#' Xmat <- model.matrix(as.formula(paste("~", paste(colnames(x), collapse = "+"))), data = as.data.frame(x))
#' Xmat[1:5,1:5]
#' x[1:5,1:5]
#' fit1 <- lmmlasso(x = Xmat, y = y, z = as.matrix(rep(1, 600)), grp = 1:600, lambda = 29, pdMat = "pdIdent")
#' summary(fit1)
#' plot(fit1)
#'
#' beta_hats <- as.matrix(fit1$coefficients)
#' s0 <- beta_hats[beta_hats!=0,, drop=F]
#'
#' # true positive rate
#' sum(paste0("V",1:20) %in% rownames(s0))/20
#'
#' # true negative rate
#' sum(paste0("V",21:1000) %ni% rownames(s0))/980
#'
#'
#'
#'
#' # quadform ----------------------------------------------------------------
#'
#' pacman::p_load(emulator)
#' pacman::p_load(microbenchmark)
#'
#' jj <- matrix(rnorm(100*1000),ncol = 1000)
#' M <- crossprod(jj,jj)
#' M.lower <- t(chol(M))
#' x <- matrix(rnorm(8),4,2)
#'
#' jj.1 <- t(x) %*% M %*% x
#' jj.2 <- quad.form(M,x)
#' jj.3 <- quad.form(M.lower,x,chol=TRUE)
#' print(jj.1)
#' print(jj.2)
#' print(jj.3)
#'
#'
#'
#'
#' emulator::quad.form
#' cprod()
#'
#'
#' # checking lower rank approx ----------------------------------------------
#'
#' library(magrittr)
#' options(digits = 4, scipen = 999)
#' X <- matrix(rnorm(n=100*40), ncol = 40)
#' X <- matrix(rnorm(n=100*40), ncol = 100)
#' dim(X)
#' Y <- rnorm(100)
#' Y <- rnorm(40)
#' xtx <- tcrossprod(X)
#' eigX <- eigen(xtx)
#'
#' U <- eigX$vectors
#' eigX$values %>% length()
#'
#' crossprod(U)[1:5,1:5]
#' tcrossprod(U)[1:5,1:5]
#'
#' svdX <- svd(X, nv = 100)
#' dim(X)
#' svdX$u %>% dim
#' svdX$d %>% length()
#' svdX$v %>% dim
#' U <- svdX$u
#' dim(U)
#' U1 <- U[,1:40]
#' U2 <- U[,41:100]
#'
#' (t(U2) %*% U2)[1:5,1:5]
#' (t(U2) %*% U2) %>% dim
#' (U2 %*% t(U2))[1:5,1:5]
#' crossprod(U2)[1:5,1:5]
#'
#' round(tcrossprod(U2)[1:5,1:5] + tcrossprod(U1)[1:5,1:5], 2)
#'
#' svdX$v %>% dim
#'
#'
#' crossprod(U1,Y) %>% dim
#'
#'
#' U1 %*% crossprod(U1,Y) %>% dim
#'
#' dim(U1)
#'
#'
#' Matrix::rankMatrix(X)
#'
#'
#'
#' eta = 0.6
#' I_nk = diag(5)
#' I_nk - eta * I_nk
#'
#' all.equal((1 - eta) * I_nk, I_nk - eta * I_nk)
#'
#'
#' all.equal(solve((1 - eta) * I_nk), 1/(1 - eta) * I_nk)
#'
