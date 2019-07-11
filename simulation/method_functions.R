## @knitr methods

# source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/ggmix/simulation/packages.R")
# source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/ggmix/simulation/functions.R")

source("/home/sahir/git_repositories/ggmix/simulation/packages.R")
source("/home/sahir/git_repositories/ggmix/simulation/functions.R")


# lasso -------------------------------------------------------------------

lasso <- new_method("lasso", "lasso",
                    method = function(model, draw) {
                      # browser()
                      fitglmnet <- glmnet::glmnet(x = draw[["xtrain_lasso"]],
                                                  y = draw[["ytrain"]],
                                                  alpha = 1,
                                                  nlambda = 50,
                                                  standardize = FALSE,
                                                  penalty.factor = c(rep(1, ncol(draw[["xtrain"]])), rep(0,10)))
                      
                      # gicp <- -2*fitglmnet$nulldev*fitglmnet$dev.ratio + log(log(length(draw[["ytrain"]])))*log(ncol(draw[["xtrain"]])) * (fitglmnet$df)
                      
                      # lambda_min_lasso <- fitglmnet$lambda[which.min(gicp)]
                      # plot(log(fitglmnet$lambda), gicp)
                      # fitglmnet$nulldev - deviance(fitglmnet) / this is the same as fitglmnet$nulldev*fitglmnet$dev.ratio
                      
                      yhat_lasso <- predict(fitglmnet, newx = draw[["xtest_lasso"]])
                      cvmat <- apply((draw[["ytest"]] - yhat_lasso)^2, 2, mean)
                      lambda_min_lasso <- fitglmnet$lambda[which.min(cvmat)]
                      
                      nz_names <- setdiff(rownames(coef(fitglmnet, s = lambda_min_lasso)[glmnet::nonzeroCoef(coef(fitglmnet, s = lambda_min_lasso)),,drop = F]),c("(Intercept)",paste0("PC",1:10)))
                      
                      model_error <- l2norm(draw[["mu_train"]] -
                                              draw[["xtrain"]] %*% coef(fitglmnet, s = lambda_min_lasso)[2:(ncol(draw[["xtrain"]]) + 1),,drop = F])
                      
                      # defined in Bertsimas et al. 2016
                      # Best Subset Selection via a Modern Optimization Lens
                      prediction_error <- model_error^2 / l2norm(draw[["mu_train"]])^2
                      
                      # inidividual level predictions
                      yhat <- predict(fitglmnet, newx = draw[["xvalidate_lasso"]], s = lambda_min_lasso)
                      
                      # yhat <- cbind(1, draw[["xtest"]]) %*% coef(fitglmnet, s = "lambda.min")[c("(Intercept)",colnames(draw[["xtest"]])),,drop = F]
                      # yhat_train <- cbind(1, draw[["xtrain"]]) %*% coef(fitglmnet, s = "lambda.min")[c("(Intercept)",colnames(draw[["xtrain"]])),,drop = F]
                      yhat_train <- predict(fitglmnet, newx = draw[["xtrain_lasso"]], s = lambda_min_lasso)
                      
                      error_var <- l2norm(yhat_train - draw[["ytrain"]])^2 / (length(draw[["ytrain"]]) - (length(nz_names)+1+10)) # add back intercept and PCs
                      
                      list(beta = coef(fitglmnet, s = lambda_min_lasso)[-1,,drop = F],
                           model_error = model_error,
                           prediction_error = prediction_error,
                           eta = NA,
                           sigma2 = NA,
                           yhat = yhat, # on the test set using principal components
                           nonzero = coef(fitglmnet, s = lambda_min_lasso)[glmnet::nonzeroCoef(coef(fitglmnet, s = lambda_min_lasso)),,drop = F],
                           nonzero_names = nz_names,
                           ytrain = draw[["ytrain"]],
                           ytest = draw[["ytest"]],
                           yvalidate = draw[["yvalidate"]],
                           error_variance = error_var,
                           causal = draw[["causal"]],
                           not_causal = draw[["not_causal"]],
                           p = ncol(draw[["xtrain"]])
                      )
                    })



# lassoCV -----------------------------------------------------------------


# this uses train/validate split only
lassoCV <- new_method(name = "lassoCV", label = "Lasso with CV on training",
                      method = function(model, draw) {
                        
                        fitglmnet <- glmnet::cv.glmnet(x = draw[["xtrain_lasso"]],
                                                       y = draw[["ytrain"]],
                                                       alpha = 1,
                                                       nlambda = 50,
                                                       standardize = FALSE,
                                                       penalty.factor = c(rep(1, ncol(draw[["xtrain"]])), rep(0,10)))
                        
                        lambda_min_lasso <- fitglmnet$lambda.min
                        
                        nz_names <- setdiff(rownames(coef(fitglmnet, s = lambda_min_lasso)[glmnet::nonzeroCoef(coef(fitglmnet, s = lambda_min_lasso)),,drop = F]),c("(Intercept)",paste0("PC",1:10)))
                        
                        # this is design matrix of selected variables along with column of 1's for intercept
                        x_nonzero_train <- cbind(1,draw[["xtrain"]][, nz_names, drop = FALSE])
                        x_nonzero_validate <- cbind(1,draw[["xvalidate"]][, nz_names, drop = FALSE])
                        
                        xtx <- crossprod(x_nonzero_train)
                        # add small ridge in case solution has more than n nonzeros:
                        diag(xtx) <- diag(xtx) + ifelse(ncol(x_nonzero_train)>nrow(x_nonzero_train),1e-4,0)
                        bb <- solve(xtx, crossprod(x_nonzero_train, draw[["ytrain"]]))
                        
                        # this is on the validation set
                        yhat <- x_nonzero_validate %*% bb
                        
                        yhat_train <- predict(fitglmnet, newx = draw[["xtrain_lasso"]], s = lambda_min_lasso)
                        error_var <- l2norm(yhat_train - draw[["ytrain"]])^2 / (length(draw[["ytrain"]]) - (length(nz_names)+1+10)) # add back intercept and PCs
                        
                        # this is the lasso fitted betas
                        beta_orig <- beta_refit <- coef(fitglmnet, s = "lambda.min")[c(colnames(draw[["xtrain"]])),,drop = F]
                        
                        # this is the refitted betas
                        beta_refit[match(rownames(bb)[-1], rownames(beta_orig)),,drop=FALSE] <- bb[-1]
                        
                        # length(draw$beta)
                        # dim(beta_orig)
                        # l2norm(beta_orig-draw$beta)
                        # l2norm(beta_refit-draw$beta)
                        # defined in Bertsimas et al. 2016
                        # Best Subset Selection via a Modern Optimization Lens
                        model_error <- l2norm(draw[["mu_train"]] - x_nonzero_train %*% bb) # this uses the refitted betas
                        prediction_error <- model_error^2 / l2norm(draw[["mu_train"]])^2
                        
                        list(beta_refit = beta_refit,
                             beta_truth = draw[["beta"]],
                             model_error = model_error,
                             prediction_error = prediction_error,
                             eta = NA,
                             sigma2 = NA,
                             yhat = yhat, # on the test set using principal components
                             nonzero = coef(fitglmnet, s = lambda_min_lasso)[glmnet::nonzeroCoef(coef(fitglmnet, s = lambda_min_lasso)),,drop = F],
                             nonzero_names = nz_names,
                             ytrain = draw[["ytrain"]],
                             # ytest = draw[["ytest"]],
                             yvalidate = draw[["yvalidate"]],
                             error_variance = error_var,
                             causal = draw[["causal"]],
                             not_causal = draw[["not_causal"]],
                             p = ncol(draw[["xtrain"]])
                        )
                      })



# lassonoPC ---------------------------------------------------------------

# this only uses PC in the training, but ignores for prediction
lassoNOPC <- new_method("lassoNOPC", "lassoNOPC",
                        method = function(model, draw) {
                          
                          # draw <- draws@draws$r1.1
                          
                          fitglmnet <- glmnet::glmnet(x = draw[["xtrain_lasso"]],
                                                      y = draw[["ytrain"]],
                                                      alpha = 1,
                                                      nlambda = 50,
                                                      standardize = FALSE,
                                                      penalty.factor = c(rep(1, ncol(draw[["xtrain"]])), rep(0,10)))
                          
                          # yhat_lasso <- predict(fitglmnet, newx = draw[["xtest_lasso"]])
                          yhat_lasso <- cbind(1, draw[["xtest"]]) %*% coef(fitglmnet)[c("(Intercept)",colnames(draw[["xtest"]])),,drop = F]
                          cvmat <- apply((draw[["ytest"]] - yhat_lasso)^2, 2, mean)
                          lambda_min_lasso <- fitglmnet$lambda[which.min(cvmat)]
                          
                          nz_names <- setdiff(rownames(coef(fitglmnet, s = lambda_min_lasso)[glmnet::nonzeroCoef(coef(fitglmnet, s = lambda_min_lasso)),,drop = F]),c("(Intercept)",paste0("PC",1:10)))
                          
                          model_error <- l2norm(draw[["mu_train"]] -
                                                  draw[["xtrain"]] %*% coef(fitglmnet, s = lambda_min_lasso)[2:(ncol(draw[["xtrain"]]) + 1),,drop = F])
                          
                          # defined in Bertsimas et al. 2016
                          # Best Subset Selection via a Modern Optimization Lens
                          prediction_error <- model_error^2 / l2norm(draw[["mu_train"]])^2
                          
                          # inidividual level predictions
                          # yhat <- predict(fitglmnet, newx = draw[["xvalidate_lasso"]], s = lambda_min_lasso)
                          yhat <- cbind(1, draw[["xvalidate"]]) %*% coef(fitglmnet, s = lambda_min_lasso)[c("(Intercept)",colnames(draw[["xvalidate"]])),,drop = F]
                          
                          # yhat <- cbind(1, draw[["xtest"]]) %*% coef(fitglmnet, s = "lambda.min")[c("(Intercept)",colnames(draw[["xtest"]])),,drop = F]
                          # yhat_train <- cbind(1, draw[["xtrain"]]) %*% coef(fitglmnet, s = "lambda.min")[c("(Intercept)",colnames(draw[["xtrain"]])),,drop = F]
                          yhat_train <- predict(fitglmnet, newx = draw[["xtrain_lasso"]], s = lambda_min_lasso)
                          
                          error_var <- l2norm(yhat_train - draw[["ytrain"]])^2 / (length(draw[["ytrain"]]) - (length(nz_names)+1+10)) # add back intercept and PCs
                          
                          list(beta = coef(fitglmnet, s = lambda_min_lasso)[-1,,drop = F],
                               model_error = model_error,
                               prediction_error = prediction_error,
                               eta = NA,
                               sigma2 = NA,
                               yhat = yhat, # on the validation set without principal components
                               nonzero = coef(fitglmnet, s = lambda_min_lasso)[glmnet::nonzeroCoef(coef(fitglmnet, s = lambda_min_lasso)),,drop = F],
                               nonzero_names = nz_names,
                               ytrain = draw[["ytrain"]],
                               ytest = draw[["ytest"]],
                               yvalidate = draw[["yvalidate"]],
                               error_variance = error_var,
                               causal = draw[["causal"]],
                               not_causal = draw[["not_causal"]],
                               p = ncol(draw[["xtrain"]])
                          )
                        })



# twostepYVC-with Cross Validation --------------------------------------------------------------

# this one uses the original y to compare to and stores the variance components (VC) (USE THIS ONE)
twostepYVCCV <- new_method("twostepYVCCV", "two step YVC with CV",
                           method = function(model, draw) {
                             
                             # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                             # fit_lme <- coxme::lmekin(Y ~ 1 + (1|id), data = pheno_dat, varlist = model$kin)
                             # newy <- residuals(fit_lme)
                             # fitglmnet <- glmnet::cv.glmnet(x = draw[["xtrain"]], y = newy, standardize = F, alpha = 1)
                             
                             # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                             
                             # browser()
                             x1 <- cbind(rep(1, nrow(draw[["xtrain"]])))
                             fit_lme <- gaston::lmm.aireml(draw[["ytrain"]], x1, K = draw[["kin_train"]])
                             gaston_resid <- draw[["ytrain"]] - (fit_lme$BLUP_omega + fit_lme$BLUP_beta)
                             
                             fitglmnet <- glmnet::cv.glmnet(x = draw[["xtrain"]],
                                                            y = gaston_resid,
                                                            nlambda = 50,
                                                            standardize = FALSE,
                                                            alpha = 1,
                                                            intercept = T)
                             
                             lambda_min_lasso <- fitglmnet$lambda.min
                             
                             nz_names <- setdiff(rownames(coef(fitglmnet, s = lambda_min_lasso)[glmnet::nonzeroCoef(coef(fitglmnet, s = lambda_min_lasso)),,drop = F]),c("(Intercept)",paste0("PC",1:10)))
                             
                             # this is design matrix of selected variables along with column of 1's for intercept
                             x_nonzero_train <- cbind(1, draw[["xtrain"]][, nz_names, drop = FALSE])
                             x_nonzero_validate <- cbind(1, draw[["xvalidate"]][, nz_names, drop = FALSE])
                             xtx <- crossprod(x_nonzero_train)
                             
                             # add small ridge in case solution has more than n nonzeros:
                             diag(xtx) <- diag(xtx) + ifelse(ncol(x_nonzero_train)>nrow(x_nonzero_train),1e-4,0)
                             bb <- solve(xtx, crossprod(x_nonzero_train, draw[["ytrain"]]))
                             # draw$beta[which(draw$beta!=0)]
                             yhat <- x_nonzero_validate %*% bb
                             # l2norm(yhat-draw$yvalidate)
                             
                             # inidividual level predictions
                             # yhat <- predict(fitglmnet, newx = draw[["xvalidate_lasso"]], s = lambda_min_lasso)
                             
                             # yhat <- cbind(1, draw[["xtest"]]) %*% coef(fitglmnet, s = "lambda.min")[c("(Intercept)",colnames(draw[["xtest"]])),,drop = F]
                             # yhat_train <- cbind(1, draw[["xtrain"]]) %*% coef(fitglmnet, s = "lambda.min")[c("(Intercept)",colnames(draw[["xtrain"]])),,drop = F]
                             yhat_train <- predict(fitglmnet, newx = draw[["xtrain"]], s = lambda_min_lasso)
                             
                             error_var <- l2norm(yhat_train - draw[["ytrain"]])^2 / (length(draw[["ytrain"]]) - (length(nz_names)+1+10)) # add back intercept and PCs
                             
                             # defined in Bertsimas et al. 2016
                             # Best Subset Selection via a Modern Optimization Lens
                             model_error <- l2norm(draw[["mu_train"]] - x_nonzero_train %*% bb) # this uses the refitted betas
                             prediction_error <- model_error^2 / l2norm(draw[["mu_train"]])^2
                             
                             # this is the lasso fitted betas
                             beta_orig <- beta_refit <- coef(fitglmnet, s = "lambda.min")[c(colnames(draw[["xtrain"]])),,drop = F]
                             
                             # this is the refitted betas
                             beta_refit[match(rownames(bb)[-1], rownames(beta_orig)),,drop=FALSE] <- bb[-1]
                             # l2norm(beta_orig - draw$beta)
                             # l2norm(beta_refit - draw$beta)
                             
                             list(beta_refit = beta_refit,
                                  beta_truth = draw[["beta"]],
                                  yhat = yhat,
                                  nonzero = coef(fitglmnet, s = lambda_min_lasso)[glmnet::nonzeroCoef(coef(fitglmnet, s = lambda_min_lasso)),,drop=F],
                                  nonzero_names = nz_names,
                                  model_error = model_error,
                                  prediction_error = prediction_error,
                                  eta = fit_lme$tau, # this is the VC for kinship
                                  sigma2 = fit_lme$sigma2, # this is VC for error
                                  ytrain = draw[["ytrain"]],
                                  yvalidate = draw[["yvalidate"]],
                                  error_variance = error_var,
                                  causal = draw[["causal"]],
                                  not_causal = draw[["not_causal"]],
                                  p = ncol(draw[["xtrain"]])
                             )
                           })


# twostepYVC --------------------------------------------------------------

# this one uses the original y to compare to and stores the variance components (VC)
# this also uses train/test/validate split
twostepYVC <- new_method("twostepYVC", "two step YVC",
                         method = function(model, draw) {
                           
                           # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                           # fit_lme <- coxme::lmekin(Y ~ 1 + (1|id), data = pheno_dat, varlist = model$kin)
                           # newy <- residuals(fit_lme)
                           # fitglmnet <- glmnet::cv.glmnet(x = draw[["xtrain"]], y = newy, standardize = F, alpha = 1)
                           
                           # pheno_dat <- data.frame(Y = draw, id = rownames(model$kin))
                           x1 <- cbind(rep(1, nrow(draw[["xtrain"]])))
                           fit_lme <- gaston::lmm.aireml(draw[["ytrain"]], x1, K = draw[["kin_train"]])
                           gaston_resid <- draw[["ytrain"]] - (fit_lme$BLUP_omega + fit_lme$BLUP_beta)
                           
                           fitglmnet <- glmnet::glmnet(x = draw[["xtrain"]], y = gaston_resid,
                                                       nlambda = 50,
                                                       standardize = FALSE, alpha = 1, intercept = T)
                           
                           yhat_lasso <- predict(fitglmnet, newx = draw[["xtest"]])
                           cvmat <- apply((draw[["ytest"]] - yhat_lasso)^2, 2, mean)
                           lambda_min_lasso <- fitglmnet$lambda[which.min(cvmat)]
                           
                           nz_names <- setdiff(rownames(coef(fitglmnet, s = lambda_min_lasso)[glmnet::nonzeroCoef(coef(fitglmnet, s = lambda_min_lasso)),,drop = F]),c("(Intercept)"))
                           
                           model_error <- l2norm(draw[["mu_train"]] -
                                                   draw[["xtrain"]] %*% coef(fitglmnet, s = lambda_min_lasso)[2:(ncol(draw[["xtrain"]]) + 1),,drop = F])
                           
                           prediction_error <- model_error^2 / l2norm(draw[["mu_train"]])^2
                           
                           yhat <- predict(fitglmnet, newx = draw[["xvalidate"]], s = lambda_min_lasso)
                           yhat_train <- predict(fitglmnet, newx = draw[["xtrain"]], s = lambda_min_lasso)
                           error_var <- l2norm(yhat_train - draw[["ytrain"]])^2 / (length(draw[["ytrain"]]) - (length(nz_names)+1)) # add back intercept
                           
                           
                           list(beta = coef(fitglmnet, s = lambda_min_lasso)[-1,,drop = F],
                                yhat = yhat,
                                nonzero = coef(fitglmnet, s = lambda_min_lasso)[glmnet::nonzeroCoef(coef(fitglmnet, s = lambda_min_lasso)),,drop=F],
                                nonzero_names = nz_names,
                                model_error = model_error,
                                prediction_error = prediction_error,
                                eta = fit_lme$tau, # this is the VC for kinship
                                sigma2 = fit_lme$sigma2, # this is VC for error
                                ytrain = draw[["ytrain"]],
                                ytest = draw[["ytest"]],
                                yvalidate = draw[["yvalidate"]],
                                error_variance = error_var,
                                causal = draw[["causal"]],
                                not_causal = draw[["not_causal"]],
                                p = ncol(draw[["xtrain"]])
                           )
                         })


# ggmixed -----------------------------------------------------------------


# this use train/test/validate splits
ggmixed <- new_method("ggmix", "ggmix",
                      method = function(model, draw) {
                        # browser()
                        fit <- ggmix(x = draw[["xtrain"]],
                                     y = draw[["ytrain"]],
                                     kinship = draw[["kin_train"]],
                                     verbose = 0,
                                     nlambda = 50,
                                     # dfmax = 500,
                                     eta_init = runif(1, min = 0.05, max = 0.95))
                        # hdbic <- gic(fit, an = log(length(draw[["ytrain"]])))
                        # hdbic <- gic(fit)
                        # plot(hdbic)
                        # sum(setdiff(rownames(coef(hdbic, type = "nonzero")), c("(Intercept)","eta","sigma2")) %in%
                        #       draw[["causal"]])
                        # hdbic$lambda.min
                        
                        predmat <- predict(fit,
                                           newx = draw[["xtest"]],
                                           type = "individual",
                                           covariance = draw[["kin_test_train"]],
                                           s = fit$lambda)
                        cvmat <- apply((draw[["ytest"]] - predmat)^2, 2, mean)
                        lambda_min_ggmix <- fit$result[which.min(cvmat), "Lambda"]
                        
                        model_error <- l2norm(draw[["mu_train"]] -
                                                draw[["xtrain"]] %*% predict(fit, s = lambda_min_ggmix, type = "coef")[2:(ncol(draw[["xtrain"]]) + 1),,drop = F])
                        
                        # inidividual level prediction
                        yhat <- predict(fit,
                                        s = lambda_min_ggmix,
                                        newx = draw[["xvalidate"]],
                                        type = "individual",
                                        covariance = draw[["kin_validate_train"]])
                        
                        # yhat_hdbic <- predict(fit,
                        #                       s = hdbic$lambda.min,
                        #                       newx = draw[["xvalidate"]],
                        #                       type = "individual",
                        #                       covariance = draw[["kin_validate_train"]])
                        #
                        # l2norm(yhat-draw[["yvalidate"]])
                        # l2norm(yhat_hdbic-draw[["yvalidate"]])
                        # mse_value <- crossprod(predict(hdbic, newx = draw[["xtrain"]]) + ranef(hdbic) - draw[["ytrain"]]) / length(draw[["ytrain"]])
                        
                        prediction_error <- model_error^2 / l2norm(draw[["mu_train"]])^2
                        
                        list(beta = predict(fit, s = lambda_min_ggmix, type = "coef")[2:(ncol(draw[["xtrain"]]) + 1),,drop = F], #this doesnt have intercept and is a 1-col matrix
                             model_error = model_error,
                             prediction_error = prediction_error,
                             nonzero = coef(fit, type = "nonzero", s = lambda_min_ggmix),
                             nonzero_names = setdiff(rownames(coef(fit, type = "nonzero", s = lambda_min_ggmix)), c("(Intercept)","eta","sigma2")),
                             # yhat = predict(hdbic, newx = draw[["xtrain"]]) + ranef(hdbic),
                             yhat = yhat,
                             ytrain = draw[["ytrain"]],
                             ytest = draw[["ytest"]],
                             yvalidate = draw[["yvalidate"]],
                             eta = coef(fit, type = "nonzero", s = lambda_min_ggmix)["eta",],
                             sigma2 = coef(fit, type = "nonzero", s = lambda_min_ggmix)["sigma2",],
                             error_variance = (1 - coef(fit, type = "nonzero", s = lambda_min_ggmix)["eta",]) * coef(fit, type = "nonzero", s = lambda_min_ggmix)["sigma2",],
                             y = draw[["ytrain"]],
                             causal = draw[["causal"]],
                             not_causal = draw[["not_causal"]],
                             p = ncol(draw[["xtrain"]])
                        )
                      })



# ggmixHDBIC --------------------------------------------------------------

# this use HDBIC to select best model from training set only, and then predicts on validation set.
ggmixedHDBIC <- new_method("ggmixHDBIC", "ggmixHDBIC",
                           method = function(model, draw) {
                             
                             fit <- ggmix(x = draw[["xtrain"]],
                                          y = draw[["ytrain"]],
                                          kinship = draw[["kin_train"]],
                                          verbose = 0,
                                          nlambda = 50,
                                          eta_init = runif(1, min = 0.05, max = 0.95))
                             
                             hdbic <- gic(fit)
                             
                             nonzero_names <- setdiff(rownames(coef(fit, type = "nonzero", s = hdbic$lambda.min)), c("(Intercept)","eta","sigma2"))
                             
                             # this is design matrix of selected variables along with column of 1's for intercept
                             x_nonzero_train <- cbind(1,draw[["xtrain"]][, nonzero_names, drop = FALSE])
                             x_nonzero_validate <- cbind(1,draw[["xvalidate"]][, nonzero_names, drop = FALSE])
                             xtx <- crossprod(x_nonzero_train)
                             
                             # add small ridge in case solution has more than n nonzeros:
                             diag(xtx) <- diag(xtx) + ifelse(ncol(x_nonzero_train)>nrow(x_nonzero_train),1e-4,0)
                             bb <- solve(xtx, crossprod(x_nonzero_train, draw[["ytrain"]]))
                             
                             # this gives same result
                             # cbind(lm.fit(x = x_nonzero_train, y = draw[["ytrain"]])$coef)
                             # this is on validation set
                             yhat <- x_nonzero_validate %*% bb
                             
                             model_error <- l2norm(draw[["mu_train"]] - x_nonzero_train %*% bb)
                             prediction_error <- model_error^2 / l2norm(draw[["mu_train"]])^2
                             
                             # this is the ggmix fitted betas
                             beta_orig <- beta_refit <- coef(hdbic)[-1,,drop=F]
                             
                             # this is the refitted betas
                             beta_refit[match(rownames(bb)[-1], rownames(beta_orig)),,drop=FALSE] <- bb[-1]
                             
                             # some checks to makes sure everything lines up
                             # all(rownames(beta_refit)==colnames(draw$xtrain))
                             # l2norm(draw$beta - beta_refit)
                             # plot(draw$beta,beta_refit)
                             # abline(a=0, b=1)
                             # beta_refit[match(rownames(bb)[-1], rownames(beta_orig)),,drop=FALSE]
                             # beta_orig[match(rownames(bb)[-1], rownames(beta_orig)),,drop=FALSE]
                             # colnames(draw$xtrain)[match(draw$causal, colnames(draw$xtrain))]
                             # cbind(colnames(draw$xtrain)[match(draw$causal, colnames(draw$xtrain))],draw$beta[match(draw$causal, colnames(draw$xtrain))])
                             # draw$beta
                             # all(rownames(coef(hdbic))[-1] == colnames(draw$xtrain))
                             
                             list(beta_refit = beta_refit, # this doesnt have intercept and is a 1-col matrix of refitted betas using OLS
                                  beta_truth = draw[["beta"]],
                                  model_error = model_error,
                                  prediction_error = prediction_error,
                                  nonzero = coef(fit, type = "nonzero", s = hdbic$lambda.min),
                                  nonzero_names = nonzero_names,
                                  yhat = yhat,
                                  ytrain = draw[["ytrain"]],
                                  yvalidate = draw[["yvalidate"]],
                                  eta = coef(fit, type = "nonzero", s = hdbic$lambda.min)["eta",],
                                  sigma2 = coef(fit, type = "nonzero", s = hdbic$lambda.min)["sigma2",],
                                  error_variance = (1 - coef(fit, type = "nonzero", s = hdbic$lambda.min)["eta",]) * coef(fit, type = "nonzero", s = hdbic$lambda.min)["sigma2",],
                                  y = draw[["ytrain"]],
                                  causal = draw[["causal"]],
                                  not_causal = draw[["not_causal"]],
                                  p = ncol(draw[["xtrain"]])
                             )
                           })




# ggmixedHDBICLMM ---------------------------------------------------------


ggmixedHDBICLMM <- new_method("ggmixHDBICLMM", "ggmixHDBIC with LMM for prediction",
                              method = function(model, draw) {
                                # browser()
                                fit <- ggmix(x = draw[["xtrain"]],
                                             y = draw[["ytrain"]],
                                             kinship = draw[["kin_train"]],
                                             verbose = 2,
                                             nlambda = 50,
                                             # dfmax = 500,
                                             eta_init = runif(1, min = 0.05, max = 0.95))
                                # hdbic <- gic(fit, an = log(length(draw[["ytrain"]])))
                                hdbic <- gic(fit)
                                # plot(hdbic)
                                # sum(setdiff(rownames(coef(hdbic, type = "nonzero")), c("(Intercept)","eta","sigma2")) %in%
                                #       draw[["causal"]])
                                # hdbic$lambda.min
                                
                                # predmat <- predict(fit,
                                #                    newx = draw[["xtest"]],
                                #                    type = "individual",
                                #                    covariance = draw[["kin_test_train"]],
                                #                    s = fit$lambda)
                                # cvmat <- apply((draw[["ytest"]] - predmat)^2, 2, mean)
                                # lambda_min_ggmix <- fit$result[which.min(cvmat), "Lambda"]
                                
                                # inidividual level prediction
                                # yhat <- predict(fit,
                                #                 s = lambda_min_ggmix,
                                #                 newx = draw[["xvalidate"]],
                                #                 type = "individual",
                                #                 covariance = draw[["kin_validate_train"]])
                                #
                                # yhat <- predict(fit,
                                #                 s = hdbic$lambda.min,
                                #                 newx = draw[["xvalidate"]],
                                #                 type = "individual",
                                #                 covariance = draw[["kin_validate_train"]])
                                
                                # l2norm(yhat-draw[["yvalidate"]])
                                # l2norm(yhat_hdbic-draw[["yvalidate"]])
                                # mse_value <- crossprod(predict(hdbic, newx = draw[["xtrain"]]) + ranef(hdbic) - draw[["ytrain"]]) / length(draw[["ytrain"]])
                                
                                nonzero_names <- setdiff(rownames(coef(fit, type = "nonzero", s = hdbic$lambda.min)), c("(Intercept)","eta","sigma2"))
                                
                                # this is design matrix of selected variables along with column of 1's for intercept
                                x_nonzero_train <- cbind(1,draw[["xtrain"]][, nonzero_names, drop = FALSE])
                                x_nonzero_validate <- cbind(1,draw[["xvalidate"]][, nonzero_names, drop = FALSE])
                                
                                # x1 <- cbind(rep(1, nrow(draw[["xtrain"]])))
                                fit_lme <- gaston::lmm.aireml(Y = draw[["ytrain"]], X = x_nonzero_train, K = draw[["kin_train"]])
                                
                                # eta
                                eta_fit <- fit_lme$tau / (fit_lme$tau + fit_lme$sigma2)
                                
                                #sigma2
                                # fit_lme$sigma2
                                # yhat <- fit_lme$tau * draw[["kin_validate_train"]] %*% fit_lme$Py
                                
                                # fit_lme$BLUP_beta %>% length()
                                # ncol(x_nonzero_train)
                                
                                yhat <- x_nonzero_validate %*% fit_lme$BLUP_beta +
                                  ggmix:::bi_future_lassofullrank(
                                    eta = eta_fit,
                                    # eta = coef(fit, type = "nonzero", s = hdbic$lambda.min)["eta",],
                                    beta = as.matrix(fit_lme$BLUP_beta),
                                    eigenvalues = fit[["ggmix_object"]][["D"]],
                                    eigenvectors = fit[["ggmix_object"]][["U"]],
                                    x = cbind(1, fit[["ggmix_object"]][["x"]][,nonzero_names, drop = FALSE]), # these are the transformed x
                                    y = fit[["ggmix_object"]][["y"]], # these are the transformed y
                                    covariance = draw[["kin_validate_train"]])
                                
                                model_error <- l2norm(draw[["mu_train"]] -
                                                        x_nonzero_train %*% as.matrix(fit_lme$BLUP_beta))
                                prediction_error <- model_error^2 / l2norm(draw[["mu_train"]])^2
                                
                                # xtx <- crossprod(x_nonzero_train)
                                # # add small ridge in case solution has more
                                # # than n nonzeros:
                                # diag(xtx) <- diag(xtx) + 1e-4
                                # bb <- solve(xtx,
                                #             crossprod(x_nonzero_train, draw[["ytrain"]]))
                                # # draw$beta[which(draw$beta!=0)]
                                # yhat <- x_nonzero_validate %*% bb
                                # # l2norm(yhat-draw$yvalidate)
                                
                                list(beta = predict(fit, s = hdbic$lambda.min, type = "coef")[2:(ncol(draw[["xtrain"]]) + 1),,drop = F], #this doesnt have intercept and is a 1-col matrix
                                     model_error = model_error,
                                     prediction_error = prediction_error,
                                     nonzero = coef(fit, type = "nonzero", s = hdbic$lambda.min),
                                     nonzero_names = nonzero_names,
                                     # yhat = predict(hdbic, newx = draw[["xtrain"]]) + ranef(hdbic),
                                     yhat = yhat,
                                     ytrain = draw[["ytrain"]],
                                     ytest = draw[["ytest"]],
                                     yvalidate = draw[["yvalidate"]],
                                     eta = eta_fit,
                                     sigma2 = fit_lme$sigma2,
                                     error_variance = (1 - eta_fit) * fit_lme$sigma2,
                                     y = draw[["ytrain"]],
                                     causal = draw[["causal"]],
                                     not_causal = draw[["not_causal"]],
                                     p = ncol(draw[["xtrain"]])
                                )
                              })


refit <- new_method_extension(name = "refit", label = "refitted",
                              method_extension = function(model, draw, out,
                                                          base_method) {
                                beta <- apply(out$beta, 2, function(b) {
                                  ii <- b != 0
                                  if (sum(ii) == 0)
                                    return(b)
                                  xtx <- crossprod(model$x[, ii])
                                  # add small ridge in case solution has more
                                  # than n nonzeros:
                                  diag(xtx) <- diag(xtx) + 1e-4
                                  bb <- solve(xtx,
                                              crossprod(model$x[, ii], draw))
                                  b[ii] <- bb
                                  return(b)
                                })
                                return(list(beta = beta,
                                            yhat = model$x %*% beta,
                                            lambda = out$lambda,
                                            df = rep(NA, ncol(beta))))
                              })



###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%###########################################

# this one uses residuals to compare to (DO NOT USE)
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

# this one uses the original y to compare to (DO nOT USE)
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




# this one uses residuals to compare to and stores the variance components (VC) (DO NOT USE)
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

