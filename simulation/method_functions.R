## @knitr methods
l2norm <- function(x) sqrt(sum(x^2))
library(glmnet)

lasso <- new_method("lasso", "Lasso",
                    method = function(model, draw) {
                      fitglmnet <- cv.glmnet(x = model$x, y = draw, alpha = 1, standardize = F)
                      list(beta = coef(fitglmnet, s = "lambda.min")[-1,,drop=F],
                           yhat = predict(fitglmnet, newx = model$x, s = "lambda.min"),
                           nonzero = coef(fitglmnet)[nonzeroCoef(coef(fitglmnet)),,drop=F],
                           nonzero_names = setdiff(rownames(coef(fitglmnet)[nonzeroCoef(coef(fitglmnet)),,drop=F]),c("(Intercept)")),
                           y = draw)
                    })


lassoPC <- new_method("lassoPC", "Lasso with 10 PC",
                    method = function(model, draw) {
                      fitglmnet <- cv.glmnet(x = model$x_lasso, y = draw, alpha = 1, standardize = T)
                      list(beta = coef(fitglmnet, s = "lambda.min")[-1,,drop=F],
                           yhat = predict(fitglmnet, newx = model$x_lasso, s = "lambda.min"),
                           nonzero = coef(fitglmnet)[nonzeroCoef(coef(fitglmnet)),,drop=F],
                           nonzero_names = setdiff(rownames(coef(fitglmnet)[nonzeroCoef(coef(fitglmnet)),,drop=F]),c("(Intercept)")),
                           y = draw)
                    })

lassoPCpf <- new_method("lassoPCpf", "Lasso with 10 PC Penalty Factor",
                      method = function(model, draw) {
                        fitglmnet <- cv.glmnet(x = model$x_lasso, y = draw, alpha = 1, standardize = T,
                                               penalty.factor = c(rep(1, ncol(model$x)), rep(0,10)))
                        
                        model_error <- l2norm(model$mu - 
                                              model$x %*% coef(fitglmnet, s = "lambda.min")[2:(ncol(model$x)+1),,drop=F])

                        list(beta = coef(fitglmnet, s = "lambda.min")[-1,,drop=F],
                             model_error = model_error,
                             eta = NA,
                             sigma2 = NA,
                             yhat = predict(fitglmnet, newx = model$x_lasso, s = "lambda.min"),
                             nonzero = coef(fitglmnet, s = "lambda.min")[nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop=F],
                             nonzero_names = setdiff(rownames(coef(fitglmnet, s = "lambda.min")[nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop=F]),c("(Intercept)")),
                             y = draw)
                      })

# lassoPCeigenpf <- new_method("lassoPCeigenpf", "Lasso with 10 PC Eigen on Kinship Penalty Factor",
#                         method = function(model, draw) {
#                           kin <- model$kin
#                           evv <- eigen(kin, symmetric=TRUE)
#                           pcs <- evv$vectors[,1:10]
#                           x_lasso <- cbind(model$x, pcs)
#                           fitglmnet <- cv.glmnet(x = x_lasso, y = draw, alpha = 1, standardize = T,
#                                                  penalty.factor = c(rep(1, 4000), rep(0,10)))
#                           list(beta = coef(fitglmnet, s = "lambda.min")[-1,,drop=F],
#                                yhat = predict(fitglmnet, newx = model$x_lasso, s = "lambda.min"),
#                                nonzero = coef(fitglmnet)[nonzeroCoef(coef(fitglmnet)),,drop=F],
#                                nonzero_names = setdiff(rownames(coef(fitglmnet)[nonzeroCoef(coef(fitglmnet)),,drop=F]),c("(Intercept)")),
#                                y = draw)
#                         })

source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/ggmix/simulation/packages.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/ggmix/simulation/functions.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/ggmix/R/fitting.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/ggmix/R/functions.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/ggmix/R/methods.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/ggmix/R/plot.R")

PENFAM <- new_method("penfam", "Penfam",
                     method = function(model, draw) {
                       fit <- penfam(x = model$x,
                                     y = draw,
                                     phi = model$kin,
                                     thresh_glmnet = 1e-10,
                                     epsilon = 1e-5,
                                     fdev = 1e-7,
                                     alpha = 1,
                                     tol.kkt = 1e-3,
                                     nlambda = 100,
                                     # an = log(log(model$n)) * log(model$n),
                                     # an = log(log(1000)),
                                     an = log(length(draw)),
                                     # lambda_min_ratio  = ifelse(model$n < model$p, 0.01, 0.001),
                                     lambda_min_ratio  = 0.05,
                                     eta_init = 0.5,
                                     maxit = 100)
                       
                       model_error <- l2norm(model$mu - 
                                               model$x %*% coef(fit, s = fit$lambda_min)[2:(ncol(model$x)+1),,drop=F])
                       
                       list(beta = fit$beta[,fit$lambda_min,drop=F], #this doesnt have intercept and is a 1-col matrix
                            model_error = model_error,
                            nonzero = predict(fit, type = "nonzero", s = fit$lambda_min),
                            nonzero_names = setdiff(rownames(predict(fit, type = "nonzero", s = fit$lambda_min)), c("(Intercept)","eta","sigma2")),
                            yhat = fit$predicted[,fit$lambda_min],
                            eta = fit$eta[,fit$lambda_min],
                            sigma2 = fit$sigma2[,fit$lambda_min],
                            y = draw
                       )
                     })

library(coxme)
pacman::p_load(gaston)

TWOSTEP <- new_method("twostep", "Two Step",
                      method = function(model, draw) {
                        
                        # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                        # fit_lme <- coxme::lmekin(Y ~ 1 + (1|id), data = pheno_dat, varlist = model$kin)
                        # newy <- residuals(fit_lme)
                        # fitglmnet <- glmnet::cv.glmnet(x = model$x, y = newy, standardize = F, alpha = 1)
                        
                        # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                        x1 <- cbind(rep(1, nrow(model$x)))
                        fit_lme <- gaston::lmm.aireml(draw, x1, K = model$kin)
                        gaston_resid <- draw - (fit_lme$BLUP_omega + fit_lme$BLUP_beta) 
                        fitglmnet <- glmnet::cv.glmnet(x = model$x, y = gaston_resid, 
                                                       standardize = T, alpha = 1, intercept = T)
                        
                        model_error <- l2norm(model$mu - 
                                                model$x %*% coef(fitglmnet, s = "lambda.min")[2:(ncol(model$x)+1),,drop=F])
                        
                        list(beta = coef(fitglmnet, s = "lambda.min")[-1,,drop=F],
                             yhat = predict(fitglmnet, newx = model$x, s = "lambda.min"),
                             nonzero = coef(fitglmnet, s = "lambda.min")[nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop=F],
                             nonzero_names = setdiff(rownames(coef(fitglmnet, s = "lambda.min")[nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop=F]),c("(Intercept)")),
                             model_error = model_error,
                             eta = NA,
                             sigma2 = NA,
                             y = draw)
                      })









#' Make folds for cross validation
#'
#' Divides the indices \code{1:n} into \code{nfolds} random folds of about the same size.
#'
#' @param n sample size
#' @param nfolds number of folds
make_folds <- function(n, nfolds) {
  nn <- round(n / nfolds)
  sizes <- rep(nn, nfolds)
  sizes[nfolds] <- sizes[nfolds] + n - nn * nfolds
  b <- c(0, cumsum(sizes))
  ii <- sample(n)
  folds <- list()
  for (i in seq(nfolds))
    folds[[i]] <- ii[seq(b[i] + 1, b[i + 1])]
  folds
}

cv <- new_method_extension("cv", "cross validated",
                           method_extension = function(model, draw, out,
                                                       base_method) {
                             nfolds <- 5
                             alpha <- base_method@settings$alpha
                             err <- matrix(NA, ncol(out$beta), nfolds)
                             ii <- make_folds(model$n, nfolds)
                             for (i in seq_along(ii)) {
                               train <- model
                               train@params$x <- model@params$x[-ii[[i]], ]
                               train@params$n <- model@params$x[-ii[[i]], ]
                               train_draw <- draw[-ii[[i]]]

                               test <- model
                               test@params$x <- model@params$x[ii[[i]], ]
                               test@params$n <- model@params$x[ii[[i]], ]
                               test_draw <- draw[ii[[i]]]
                               fit <- base_method@method(model = train,
                                                         draw = train_draw,
                                                         alpha = alpha,
                                                         lambda = out$lambda)
                               yhat <- test$x %*% fit$beta
                               ll <- seq(ncol(yhat))
                               err[ll, i] <- colMeans((yhat - test_draw)^2)
                             }
                             m <- rowMeans(err)
                             se <- apply(err, 1, sd) / sqrt(nfolds)
                             imin <- which.min(m)
                             ioneserule <- max(which(m <= m[imin] + se[imin]))
                             list(err = err, m = m, se = se, imin = imin,
                                  ioneserule = ioneserule,
                                  beta = out$beta[, imin],
                                  yhat = model$x %*% out$beta[, imin],
                                  alpha = alpha)
                           })
