## @knitr methods

# source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/ggmix/simulation/packages.R")
# source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/ggmix/simulation/functions.R")

source("/home/sahir/git_repositories/ggmix/simulation/packages.R")
source("/home/sahir/git_repositories/ggmix/simulation/functions.R")

# lasso <- new_method("lasso", "Lasso",
#                     method = function(model, draw) {
#                       fitglmnet <- cv.glmnet(x = draw[["xtrain"]], y = draw, alpha = 1, standardize = F)
#                       list(beta = coef(fitglmnet, s = "lambda.min")[-1,,drop=F],
#                            yhat = predict(fitglmnet, newx = draw[["xtrain"]], s = "lambda.min"),
#                            nonzero = coef(fitglmnet)[nonzeroCoef(coef(fitglmnet)),,drop=F],
#                            nonzero_names = setdiff(rownames(coef(fitglmnet)[nonzeroCoef(coef(fitglmnet)),,drop=F]),c("(Intercept)")),
#                            y = draw)
#                     })

lasso <- new_method("lasso", "lasso",
                      method = function(model, draw) {
                        fitglmnet <- glmnet::cv.glmnet(x = draw[["xtrain_lasso"]],
                                                       y = draw[["ytrain"]],
                                                       alpha = 1,
                                                       standardize = T,
                                                       penalty.factor = c(rep(1, ncol(draw[["xtrain"]])), rep(0,10)))

                        nz_names <- setdiff(rownames(coef(fitglmnet, s = "lambda.min")[glmnet::nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop = F]),c("(Intercept)",paste0("PC",1:10)))

                        model_error <- l2norm(draw[["mu_train"]] -
                          draw[["xtrain"]] %*% coef(fitglmnet, s = "lambda.min")[2:(ncol(draw[["xtrain"]]) + 1),,drop = F])

                        # defined in Bertsimas et al. 2016
                        # Best Subset Selection via a Modern Optimization Lens
                        prediction_error <- model_error^2 / l2norm(draw[["mu_train"]])^2

                        # inidividual level predictions
                        yhat <- predict(fitglmnet, newx = draw[["xtest_lasso"]], s = "lambda.min")

                        # yhat <- cbind(1, draw[["xtest"]]) %*% coef(fitglmnet, s = "lambda.min")[c("(Intercept)",colnames(draw[["xtest"]])),,drop = F]
                        # yhat_train <- cbind(1, draw[["xtrain"]]) %*% coef(fitglmnet, s = "lambda.min")[c("(Intercept)",colnames(draw[["xtrain"]])),,drop = F]
                        yhat_train <- predict(fitglmnet, newx = draw[["xtrain_lasso"]], s = "lambda.min")

                        error_var <- l2norm(yhat_train - draw[["ytrain"]])^2 / (length(draw[["ytrain"]]) - length(nz_names))

                        list(beta = coef(fitglmnet, s = "lambda.min")[-1,,drop = F],
                             model_error = model_error,
                             prediction_error = prediction_error,
                             eta = NA,
                             sigma2 = NA,
                             yhat = yhat, # on the test set using principal components
                             nonzero = coef(fitglmnet, s = "lambda.min")[glmnet::nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop = F],
                             nonzero_names = nz_names,
                             ytrain = draw[["ytrain"]],
                             ytest = draw[["ytest"]],
                             error_variance = error_var,
                             causal = draw[["causal"]],
                             not_causal = draw[["not_causal"]],
                             p = ncol(draw[["xtrain"]])
                        )
                      })

ggmixed <- new_method("ggmix", "ggmix",
                     method = function(model, draw) {
                       fit <- ggmix(x = draw[["xtrain"]],
                                    y = draw[["ytrain"]],
                                    kinship = draw[["kin_train"]],
                                    verbose = 1, dfmax = 100)
                       # hdbic <- gic(fit, an = log(length(draw[["ytrain"]])))
                       hdbic <- gic(fit)

                       model_error <- l2norm(draw[["mu_train"]] -
                                               draw[["xtrain"]] %*% coef(hdbic)[2:(ncol(draw[["xtrain"]]) + 1),,drop = F])

                       # inidividual level prediction
                       yhat <- predict(hdbic, s="lambda.min", newx = draw[["xtest"]],
                                       type = "individual", covariance = draw[["kin_test_train"]])

                       # mse_value <- crossprod(predict(hdbic, newx = draw[["xtrain"]]) + ranef(hdbic) - draw[["ytrain"]]) / length(draw[["ytrain"]])

                       prediction_error <- model_error^2 / l2norm(draw[["mu_train"]])^2

                       list(beta = coef(hdbic)[2:(ncol(draw[["xtrain"]]) + 1),,drop = F], #this doesnt have intercept and is a 1-col matrix
                            model_error = model_error,
                            prediction_error = prediction_error,
                            nonzero = coef(hdbic, type = "nonzero"),
                            nonzero_names = setdiff(rownames(coef(hdbic, type = "nonzero")), c("(Intercept)","eta","sigma2")),
                            # yhat = predict(hdbic, newx = draw[["xtrain"]]) + ranef(hdbic),
                            yhat = yhat,
                            ytrain = draw[["ytrain"]],
                            ytest = draw[["ytest"]],
                            eta = coef(hdbic, type = "nonzero")["eta",],
                            sigma2 = coef(hdbic, type = "nonzero")["sigma2",],
                            error_variance = (1 - coef(hdbic, type = "nonzero")["eta",]) * coef(hdbic, type = "nonzero")["sigma2",],
                            y = draw[["ytrain"]],
                            causal = draw[["causal"]],
                            not_causal = draw[["not_causal"]],
                            p = ncol(draw[["xtrain"]])
                       )
                     })


# this one uses residuals to compare to
twostep <- new_method("twostep", "two step",
                      method = function(model, draw) {

                        # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                        # fit_lme <- coxme::lmekin(Y ~ 1 + (1|id), data = pheno_dat, varlist = model$kin)
                        # newy <- residuals(fit_lme)
                        # fitglmnet <- glmnet::cv.glmnet(x = draw[["xtrain"]], y = newy, standardize = F, alpha = 1)

                        # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                        x1 <- cbind(rep(1, nrow(draw[["xtrain"]])))
                        fit_lme <- gaston::lmm.aireml(draw[["ytrain"]], x1, K = draw[["kin"]])
                        gaston_resid <- draw[["ytrain"]] - (fit_lme$BLUP_omega + fit_lme$BLUP_beta)
                        fitglmnet <- glmnet::cv.glmnet(x = draw[["xtrain"]], y = gaston_resid,
                                                       standardize = T, alpha = 1, intercept = T)

                        nz_names <- setdiff(rownames(coef(fitglmnet, s = "lambda.min")[glmnet::nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop = F]),c("(Intercept)"))

                        model_error <- l2norm(draw[["mutrain"]] -
                          draw[["xtrain"]] %*% coef(fitglmnet, s = "lambda.min")[2:(ncol(draw[["xtrain"]]) + 1),,drop = F])

                        prediction_error <- model_error^2 / l2norm(draw[["mutrain"]])^2

                        # this should be compared with gaston_resid
                        yhat <- predict(fitglmnet, newx = draw[["xtrain"]], s = "lambda.min")
                        error_var <- l2norm(yhat - gaston_resid)^2 / (length(draw[["ytrain"]]) - length(nz_names))


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
                             p = ncol(draw[["xtrain"]])
                        )
                      })

# this one uses the original y to compare to
twostepY <- new_method("twostepY", "two step Y",
                      method = function(model, draw) {

                        # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                        # fit_lme <- coxme::lmekin(Y ~ 1 + (1|id), data = pheno_dat, varlist = model$kin)
                        # newy <- residuals(fit_lme)
                        # fitglmnet <- glmnet::cv.glmnet(x = draw[["xtrain"]], y = newy, standardize = F, alpha = 1)

                        # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                        x1 <- cbind(rep(1, nrow(draw[["xtrain"]])))
                        fit_lme <- gaston::lmm.aireml(draw[["ytrain"]], x1, K = draw[["kin_train"]])
                        gaston_resid <- draw[["ytrain"]] - (fit_lme$BLUP_omega + fit_lme$BLUP_beta)
                        fitglmnet <- glmnet::cv.glmnet(x = draw[["xtrain"]], y = gaston_resid,
                                                       standardize = T, alpha = 1, intercept = T)

                        nz_names <- setdiff(rownames(coef(fitglmnet, s = "lambda.min")[glmnet::nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop = F]),c("(Intercept)"))

                        model_error <- l2norm(draw[["mu_train"]] -
                                                draw[["xtrain"]] %*% coef(fitglmnet, s = "lambda.min")[2:(ncol(draw[["xtrain"]]) + 1),,drop = F])

                        prediction_error <- model_error^2 / l2norm(draw[["mu_train"]])^2

                        yhat <- predict(fitglmnet, newx = draw[["xtest"]], s = "lambda.min")
                        yhat_train <- predict(fitglmnet, newx = draw[["xtrain"]], s = "lambda.min")
                        error_var <- l2norm(yhat_train - draw[["ytrain"]])^2 / (length(draw[["ytrain"]]) - length(nz_names))


                        list(beta = coef(fitglmnet, s = "lambda.min")[-1,,drop = F],
                             yhat = yhat,
                             nonzero = coef(fitglmnet, s = "lambda.min")[glmnet::nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop=F],
                             nonzero_names = nz_names,
                             model_error = model_error,
                             prediction_error = prediction_error,
                             eta = NA,
                             sigma2 = NA,
                             y = draw[["ytrain"]],
                             ytest = draw[["ytest"]],
                             error_variance = error_var,
                             causal = draw[["causal"]],
                             not_causal = draw[["not_causal"]],
                             p = ncol(draw[["xtrain"]])
                        )
                      })




# this one uses residuals to compare to and stores the variance components (VC)
twostepVC <- new_method("twostep", "two step",
                      method = function(model, draw) {

                        # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                        # fit_lme <- coxme::lmekin(Y ~ 1 + (1|id), data = pheno_dat, varlist = model$kin)
                        # newy <- residuals(fit_lme)
                        # fitglmnet <- glmnet::cv.glmnet(x = draw[["xtrain"]], y = newy, standardize = F, alpha = 1)

                        # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                        x1 <- cbind(rep(1, nrow(draw[["xtrain"]])))
                        fit_lme <- gaston::lmm.aireml(draw[["ytrain"]], x1, K = draw[["kin"]])
                        gaston_resid <- draw[["ytrain"]] - (fit_lme$BLUP_omega + fit_lme$BLUP_beta)
                        fitglmnet <- glmnet::cv.glmnet(x = draw[["xtrain"]], y = gaston_resid,
                                                       standardize = T, alpha = 1, intercept = T)

                        nz_names <- setdiff(rownames(coef(fitglmnet, s = "lambda.min")[glmnet::nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop = F]),c("(Intercept)"))

                        model_error <- l2norm(draw[["mutrain"]] -
                                                draw[["xtrain"]] %*% coef(fitglmnet, s = "lambda.min")[2:(ncol(draw[["xtrain"]]) + 1),,drop = F])

                        prediction_error <- model_error^2 / l2norm(draw[["mutrain"]])^2

                        # this should be compared with gaston_resid
                        yhat <- predict(fitglmnet, newx = draw[["xtrain"]], s = "lambda.min")
                        error_var <- l2norm(yhat - gaston_resid)^2 / (length(draw[["ytrain"]]) - length(nz_names))


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
                             p = ncol(draw[["xtrain"]])
                        )
                      })

# this one uses the original y to compare to and stores the variance components (VC)
twostepYVC <- new_method("twostepY", "two step Y",
                       method = function(model, draw) {

                         # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                         # fit_lme <- coxme::lmekin(Y ~ 1 + (1|id), data = pheno_dat, varlist = model$kin)
                         # newy <- residuals(fit_lme)
                         # fitglmnet <- glmnet::cv.glmnet(x = draw[["xtrain"]], y = newy, standardize = F, alpha = 1)

                         # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                         x1 <- cbind(rep(1, nrow(draw[["xtrain"]])))
                         fit_lme <- gaston::lmm.aireml(draw[["ytrain"]], x1, K = draw[["kin"]])
                         gaston_resid <- draw[["ytrain"]] - (fit_lme$BLUP_omega + fit_lme$BLUP_beta)
                         fitglmnet <- glmnet::cv.glmnet(x = draw[["xtrain"]], y = gaston_resid,
                                                        standardize = T, alpha = 1, intercept = T)

                         nz_names <- setdiff(rownames(coef(fitglmnet, s = "lambda.min")[glmnet::nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop = F]),c("(Intercept)"))

                         model_error <- l2norm(draw[["mutrain"]] -
                                                 draw[["xtrain"]] %*% coef(fitglmnet, s = "lambda.min")[2:(ncol(draw[["xtrain"]]) + 1),,drop = F])

                         prediction_error <- model_error^2 / l2norm(draw[["mutrain"]])^2

                         yhat <- predict(fitglmnet, newx = draw[["xtrain"]], s = "lambda.min")
                         error_var <- l2norm(yhat - draw[["ytrain"]])^2 / (length(draw[["ytrain"]]) - length(nz_names))


                         list(beta = coef(fitglmnet, s = "lambda.min")[-1,,drop = F],
                              yhat = yhat,
                              nonzero = coef(fitglmnet, s = "lambda.min")[glmnet::nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop=F],
                              nonzero_names = nz_names,
                              model_error = model_error,
                              prediction_error = prediction_error,
                              eta = fit_lme$tau, # this is the VC for kinship
                              sigma2 = fit_lme$sigma2, # this is VC for error
                              y = draw[["ytrain"]],
                              error_variance = error_var,
                              causal = draw[["causal"]],
                              not_causal = draw[["not_causal"]],
                              p = ncol(draw[["xtrain"]])
                         )
                       })
