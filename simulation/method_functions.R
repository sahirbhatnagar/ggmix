## @knitr methods

source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/ggmix/simulation/packages.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/ggmix/simulation/functions.R")

# source("~/git_repositories/ggmix/simulation/packages.R")
# source("~/git_repositories/ggmix/simulation/functions.R")

# lasso <- new_method("lasso", "Lasso",
#                     method = function(model, draw) {
#                       fitglmnet <- cv.glmnet(x = draw[["Xtest"]], y = draw, alpha = 1, standardize = F)
#                       list(beta = coef(fitglmnet, s = "lambda.min")[-1,,drop=F],
#                            yhat = predict(fitglmnet, newx = draw[["Xtest"]], s = "lambda.min"),
#                            nonzero = coef(fitglmnet)[nonzeroCoef(coef(fitglmnet)),,drop=F],
#                            nonzero_names = setdiff(rownames(coef(fitglmnet)[nonzeroCoef(coef(fitglmnet)),,drop=F]),c("(Intercept)")),
#                            y = draw)
#                     })

lasso <- new_method("lasso", "lasso",
                      method = function(model, draw) {
                        fitglmnet <- glmnet::cv.glmnet(x = draw[["x_lasso"]], y = draw[["y"]],
                                               alpha = 1, standardize = T,
                                    penalty.factor = c(rep(1, ncol(draw[["Xtest"]])),
                                                       rep(0,10)))

                        nz_names <- setdiff(rownames(coef(fitglmnet, s = "lambda.min")[glmnet::nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop = F]),c("(Intercept)"))

                        model_error <- l2norm(draw[["mu"]] -
                          draw[["Xtest"]] %*% coef(fitglmnet, s = "lambda.min")[2:(ncol(draw[["Xtest"]]) + 1),,drop = F])
                        # defined in Bertsimas et al. 2016
                        #Best Subset Selection via a Modern Optimization Lens
                        prediction_error <- model_error^2 / l2norm(draw[["mu"]])^2

                        yhat <- predict(fitglmnet, newx = draw[["x_lasso"]], s = "lambda.min")
                        # yhat <- cbind(1, draw[["Xtest"]]) %*% coef(fitglmnet, s = "lambda.min")[c("(Intercept)",colnames(draw[["Xtest"]])),,drop = F]
                        error_var <- l2norm(yhat - draw[["y"]])^2 / (length(draw[["y"]]) - length(nz_names))

                        list(beta = coef(fitglmnet, s = "lambda.min")[-1,,drop = F],
                             model_error = model_error,
                             prediction_error = prediction_error,
                             eta = NA,
                             sigma2 = NA,
                             yhat = yhat,
                             nonzero = coef(fitglmnet, s = "lambda.min")[glmnet::nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop = F],
                             nonzero_names = nz_names,
                             y = draw[["y"]],
                             error_variance = error_var,
                             causal = draw[["causal"]],
                             not_causal = draw[["not_causal"]],
                             p = ncol(draw[["Xtest"]])
                        )
                      })

ggmixed <- new_method("ggmix", "ggmix",
                     method = function(model, draw) {
                       fit <- ggmix(x = draw[["Xtest"]],
                                    y = draw[["y"]],
                                    kinship = draw[["kin"]],
                                    verbose = 1, dfmax = 100)
                       hdbic <- gic(fit)

                       model_error <- l2norm(draw[["mu"]] -
                                               draw[["Xtest"]] %*% coef(hdbic)[2:(ncol(draw[["Xtest"]]) + 1),,drop = F])

                       # mse_value <- crossprod(predict(hdbic, newx = draw[["Xtest"]]) + ranef(hdbic) - draw[["y"]]) / length(draw[["y"]])

                       prediction_error <- model_error^2 / l2norm(draw[["mu"]])^2

                       list(beta = coef(hdbic)[2:(ncol(draw[["Xtest"]]) + 1),,drop = F], #this doesnt have intercept and is a 1-col matrix
                            model_error = model_error,
                            prediction_error = prediction_error,
                            nonzero = coef(hdbic, type = "nonzero"),
                            nonzero_names = setdiff(rownames(coef(hdbic, type = "nonzero")), c("(Intercept)","eta","sigma2")),
                            yhat = predict(hdbic, newx = draw[["Xtest"]]) + ranef(hdbic),
                            eta = coef(hdbic, type = "nonzero")["eta",],
                            sigma2 = coef(hdbic, type = "nonzero")["sigma2",],
                            error_variance = (1 - coef(hdbic, type = "nonzero")["eta",]) * coef(hdbic, type = "nonzero")["sigma2",],
                            y = draw[["y"]],
                            causal = draw[["causal"]],
                            not_causal = draw[["not_causal"]],
                            p = ncol(draw[["Xtest"]])
                       )
                     })


# this one uses residuals to compare to
twostep <- new_method("twostep", "two step",
                      method = function(model, draw) {

                        # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                        # fit_lme <- coxme::lmekin(Y ~ 1 + (1|id), data = pheno_dat, varlist = model$kin)
                        # newy <- residuals(fit_lme)
                        # fitglmnet <- glmnet::cv.glmnet(x = draw[["Xtest"]], y = newy, standardize = F, alpha = 1)

                        # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                        x1 <- cbind(rep(1, nrow(draw[["Xtest"]])))
                        fit_lme <- gaston::lmm.aireml(draw[["y"]], x1, K = draw[["kin"]])
                        gaston_resid <- draw[["y"]] - (fit_lme$BLUP_omega + fit_lme$BLUP_beta)
                        fitglmnet <- glmnet::cv.glmnet(x = draw[["Xtest"]], y = gaston_resid,
                                                       standardize = T, alpha = 1, intercept = T)

                        nz_names <- setdiff(rownames(coef(fitglmnet, s = "lambda.min")[glmnet::nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop = F]),c("(Intercept)"))

                        model_error <- l2norm(draw[["mu"]] -
                          draw[["Xtest"]] %*% coef(fitglmnet, s = "lambda.min")[2:(ncol(draw[["Xtest"]]) + 1),,drop = F])

                        prediction_error <- model_error^2 / l2norm(draw[["mu"]])^2

                        # this should be compared with gaston_resid
                        yhat <- predict(fitglmnet, newx = draw[["Xtest"]], s = "lambda.min")
                        error_var <- l2norm(yhat - gaston_resid)^2 / (length(draw[["y"]]) - length(nz_names))


                        list(beta = coef(fitglmnet, s = "lambda.min")[-1,,drop = F],
                             yhat = yhat,
                             nonzero = coef(fitglmnet, s = "lambda.min")[glmnet::nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop=F],
                             nonzero_names = nz_names,
                             model_error = model_error,
                             prediction_error = prediction_error,
                             eta = NA,
                             sigma2 = NA,
                             y = gaston_resid,
                             error_variance = error_var,
                             causal = draw[["causal"]],
                             not_causal = draw[["not_causal"]],
                             p = ncol(draw[["Xtest"]])
                        )
                      })

# this one uses the original y to compare to
twostepY <- new_method("twostepY", "two step Y",
                      method = function(model, draw) {
                        
                        # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                        # fit_lme <- coxme::lmekin(Y ~ 1 + (1|id), data = pheno_dat, varlist = model$kin)
                        # newy <- residuals(fit_lme)
                        # fitglmnet <- glmnet::cv.glmnet(x = draw[["Xtest"]], y = newy, standardize = F, alpha = 1)
                        
                        # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                        x1 <- cbind(rep(1, nrow(draw[["Xtest"]])))
                        fit_lme <- gaston::lmm.aireml(draw[["y"]], x1, K = draw[["kin"]])
                        gaston_resid <- draw[["y"]] - (fit_lme$BLUP_omega + fit_lme$BLUP_beta)
                        fitglmnet <- glmnet::cv.glmnet(x = draw[["Xtest"]], y = gaston_resid,
                                                       standardize = T, alpha = 1, intercept = T)
                        
                        nz_names <- setdiff(rownames(coef(fitglmnet, s = "lambda.min")[glmnet::nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop = F]),c("(Intercept)"))
                        
                        model_error <- l2norm(draw[["mu"]] -
                                                draw[["Xtest"]] %*% coef(fitglmnet, s = "lambda.min")[2:(ncol(draw[["Xtest"]]) + 1),,drop = F])
                        
                        prediction_error <- model_error^2 / l2norm(draw[["mu"]])^2
                        
                        yhat <- predict(fitglmnet, newx = draw[["Xtest"]], s = "lambda.min")
                        error_var <- l2norm(yhat - draw[["y"]])^2 / (length(draw[["y"]]) - length(nz_names))
                        
                        
                        list(beta = coef(fitglmnet, s = "lambda.min")[-1,,drop = F],
                             yhat = yhat,
                             nonzero = coef(fitglmnet, s = "lambda.min")[glmnet::nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop=F],
                             nonzero_names = nz_names,
                             model_error = model_error,
                             prediction_error = prediction_error,
                             eta = NA,
                             sigma2 = NA,
                             y = draw[["y"]],
                             error_variance = error_var,
                             causal = draw[["causal"]],
                             not_causal = draw[["not_causal"]],
                             p = ncol(draw[["Xtest"]])
                        )
                      })




# this one uses residuals to compare to and stores the variance components (VC)
twostepVC <- new_method("twostep", "two step",
                      method = function(model, draw) {
                        
                        # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                        # fit_lme <- coxme::lmekin(Y ~ 1 + (1|id), data = pheno_dat, varlist = model$kin)
                        # newy <- residuals(fit_lme)
                        # fitglmnet <- glmnet::cv.glmnet(x = draw[["Xtest"]], y = newy, standardize = F, alpha = 1)
                        
                        # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                        x1 <- cbind(rep(1, nrow(draw[["Xtest"]])))
                        fit_lme <- gaston::lmm.aireml(draw[["y"]], x1, K = draw[["kin"]])
                        gaston_resid <- draw[["y"]] - (fit_lme$BLUP_omega + fit_lme$BLUP_beta)
                        fitglmnet <- glmnet::cv.glmnet(x = draw[["Xtest"]], y = gaston_resid,
                                                       standardize = T, alpha = 1, intercept = T)
                        
                        nz_names <- setdiff(rownames(coef(fitglmnet, s = "lambda.min")[glmnet::nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop = F]),c("(Intercept)"))
                        
                        model_error <- l2norm(draw[["mu"]] -
                                                draw[["Xtest"]] %*% coef(fitglmnet, s = "lambda.min")[2:(ncol(draw[["Xtest"]]) + 1),,drop = F])
                        
                        prediction_error <- model_error^2 / l2norm(draw[["mu"]])^2
                        
                        # this should be compared with gaston_resid
                        yhat <- predict(fitglmnet, newx = draw[["Xtest"]], s = "lambda.min")
                        error_var <- l2norm(yhat - gaston_resid)^2 / (length(draw[["y"]]) - length(nz_names))
                        
                        
                        list(beta = coef(fitglmnet, s = "lambda.min")[-1,,drop = F],
                             yhat = yhat,
                             nonzero = coef(fitglmnet, s = "lambda.min")[glmnet::nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop=F],
                             nonzero_names = nz_names,
                             model_error = model_error,
                             prediction_error = prediction_error,
                             eta = fit_lme$tau, # this is the VC for kinship
                             sigma2 = fit_lme$sigma2, # this is VC for error
                             y = gaston_resid,
                             error_variance = error_var,
                             causal = draw[["causal"]],
                             not_causal = draw[["not_causal"]],
                             p = ncol(draw[["Xtest"]])
                        )
                      })

# this one uses the original y to compare to and stores the variance components (VC)
twostepYVC <- new_method("twostepY", "two step Y",
                       method = function(model, draw) {
                         
                         # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                         # fit_lme <- coxme::lmekin(Y ~ 1 + (1|id), data = pheno_dat, varlist = model$kin)
                         # newy <- residuals(fit_lme)
                         # fitglmnet <- glmnet::cv.glmnet(x = draw[["Xtest"]], y = newy, standardize = F, alpha = 1)
                         
                         # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                         x1 <- cbind(rep(1, nrow(draw[["Xtest"]])))
                         fit_lme <- gaston::lmm.aireml(draw[["y"]], x1, K = draw[["kin"]])
                         gaston_resid <- draw[["y"]] - (fit_lme$BLUP_omega + fit_lme$BLUP_beta)
                         fitglmnet <- glmnet::cv.glmnet(x = draw[["Xtest"]], y = gaston_resid,
                                                        standardize = T, alpha = 1, intercept = T)
                         
                         nz_names <- setdiff(rownames(coef(fitglmnet, s = "lambda.min")[glmnet::nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop = F]),c("(Intercept)"))
                         
                         model_error <- l2norm(draw[["mu"]] -
                                                 draw[["Xtest"]] %*% coef(fitglmnet, s = "lambda.min")[2:(ncol(draw[["Xtest"]]) + 1),,drop = F])
                         
                         prediction_error <- model_error^2 / l2norm(draw[["mu"]])^2
                         
                         yhat <- predict(fitglmnet, newx = draw[["Xtest"]], s = "lambda.min")
                         error_var <- l2norm(yhat - draw[["y"]])^2 / (length(draw[["y"]]) - length(nz_names))
                         
                         
                         list(beta = coef(fitglmnet, s = "lambda.min")[-1,,drop = F],
                              yhat = yhat,
                              nonzero = coef(fitglmnet, s = "lambda.min")[glmnet::nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop=F],
                              nonzero_names = nz_names,
                              model_error = model_error,
                              prediction_error = prediction_error,
                              eta = fit_lme$tau, # this is the VC for kinship
                              sigma2 = fit_lme$sigma2, # this is VC for error
                              y = draw[["y"]],
                              error_variance = error_var,
                              causal = draw[["causal"]],
                              not_causal = draw[["not_causal"]],
                              p = ncol(draw[["Xtest"]])
                         )
                       })