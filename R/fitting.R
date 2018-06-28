#' Fit a linear mixed model with lasso or elasticnet regularization
#'
#' @description this is an older version of the code. currently not being used.
#' use \link{\code{lowrank}} instead
#' @param an is the constant used in the BIC calculation. The default choice is
#'   from Wang et al. (2009)
#' @references Wang, H., Li, B. and Leng, C. (2009) Shrinkage tuning parameter
#'   selection with a diverging number of parameters.J. R. Statist. Soc. B, 71,
#'   671–683.
# ggmix_old <- function(x, y, phi, lambda = NULL,
#                       lambda_min_ratio = ifelse(n < p, 0.01, 0.001),
#                       nlambda = 100,
#                       eta_init = 0.5,
#                       maxit = 100,
#                       fdev = 1e-4,
#                       alpha = 1, # elastic net mixing param. 1 is lasso, 0 is ridge
#                       thresh_glmnet = 1e-10, # this is for glmnet
#                       epsilon = 1e-5, # this is for ggmix
#                       an = log(log(n)) * log(n),
#                       tol.kkt = 1e-9) {
#
#   # rm(list=ls())
#   # source("~/git_repositories/ggmix/R/fitting.R")
#   # source("~/git_repositories/ggmix/R/functions.R")
#   # source("~/git_repositories/ggmix/R/methods.R")
#   # source("~/git_repositories/ggmix/R/plot.R")
#   # source("~/git_repositories/ggmix/R/sim-data.R")
#   # x <- X
#   # y <- Y
#   # phi <- Phi
#   # lambda_min_ratio <- ifelse(n < p, 0.01, 0.001)
#   # nlambda <- 100
#   # #convergence criterion
#   # thresh_glmnet <- 1e-10
#   # epsilon <- 1e-5
#   # tol.kkt <- 1e-5
#   # maxit <- 100
#   # eta_init <- 0.4
#   # an = log(log(600)) * log(600)
#   # lambda <- 0.10
#   # exact = F
#   # fdev <- 1e-4
#   # ======================================
#
#   this.call <- match.call()
#   if (!is.matrix(x)) {
#     stop("x has to be a matrix")
#   }
#   if (any(is.na(x))) {
#     stop("Missing values in x not allowed!")
#   }
#
#   np <- dim(x)
#   if (is.null(np) | (np[2] <= 1)) {
#     stop("x should be a matrix with 2 or more columns")
#   }
#
#   n <- np[[1]]
#   p <- np[[2]]
#
#   # add column of 1s to x for intercept
#   x <- cbind(beta0 = 1, x)
#   x[1:5, 1:5]
#
#   phi_eigen <- eigen(phi)
#   # this is a N_T x N_T matrix
#   U <- phi_eigen$vectors
#
#   # dim(U)
#   # vector of length N_T
#   Lambda <- phi_eigen$values
#
#   # Phi Inverse (used for prediction of random effects)
#   D_inv <- diag(1 / Lambda)
#   Phi_inv <- U %*% D_inv %*% t(U)
#
#   utx <- crossprod(U, x)
#   uty <- crossprod(U, y)
#
#   # get sequence of tuning parameters
#   lamb <- lambda_sequence(
#     x = utx, y = uty, eigenvalues = Lambda, nlambda = nlambda,
#     lambda_min_ratio = lambda_min_ratio,
#     eta_init = eta_init, epsilon = epsilon,
#     tol.kkt = tol.kkt
#   )
#
#   lambda_max <- lamb$sequence[[1]]
#
#   tuning_params_mat <- matrix(lamb$sequence, nrow = 1, ncol = nlambda, byrow = T)
#   dimnames(tuning_params_mat)[[1]] <- list("lambda")
#   dimnames(tuning_params_mat)[[2]] <- paste0("s", seq_len(nlambda))
#   lambda_names <- dimnames(tuning_params_mat)[[2]]
#
#   # coefficient_mat <- matrix(nrow = p + 3,
#   #                           ncol = nlambda,
#   #                           dimnames = list(c(paste0("beta",0:p), "eta","sigma2"),
#   #                                           lambda_names))
#
#   coefficient_mat <- matrix(
#     nrow = p + 3,
#     ncol = nlambda,
#     dimnames = list(
#       c(colnames(x), "eta", "sigma2"),
#       lambda_names
#     )
#   )
#
#   randomeff_mat <- matrix(
#     nrow = n,
#     ncol = nlambda,
#     dimnames = list(
#       c(paste0("Subject", 1:n)),
#       lambda_names
#     )
#   )
#
#   fitted_mat <- matrix(
#     nrow = n,
#     ncol = nlambda,
#     dimnames = list(
#       c(paste0("Subject", 1:n)),
#       lambda_names
#     )
#   )
#
#   predicted_mat <- matrix(
#     nrow = n,
#     ncol = nlambda,
#     dimnames = list(
#       c(paste0("Subject", 1:n)),
#       lambda_names
#     )
#   )
#
#   resid_mat <- matrix(
#     nrow = n,
#     ncol = nlambda,
#     dimnames = list(
#       c(paste0("Subject", 1:n)),
#       lambda_names
#     )
#   )
#
#   out_print <- matrix(NA,
#     nrow = nlambda, ncol = 10,
#     dimnames = list(
#       lambda_names,
#       c(
#         "Df",
#         # "%Dev",
#         "Deviance",
#         "Lambda",
#         "BIC",
#         "kkt_beta0",
#         "kkt_eta",
#         "kkt_sigma2",
#         "kkt_beta_nonzero",
#         "kkt_beta_subgr", "converged"
#       )
#     )
#   )
#
#   pb <- progress::progress_bar$new(
#     format = "  fitting over all tuning parameters [:bar] :percent eta: :eta",
#     total = nlambda, clear = FALSE, width = 90
#   )
#   pb$tick(0)
#
#   beta_init <- matrix(0, nrow = p + 1, ncol = 1)
#
#   # initial value for eta from lambda_sequence results
#   eta_init <- lamb$eta
#
#   # closed form solution value for sigma^2
#   sigma2_init <- (1 / n) * sum(((uty - utx %*% beta_init)^2) / (1 + eta_init * (Lambda - 1)))
#
#
#   for (LAMBDA in lambda_names) {
#
#     # LAMBDA <- "s1"
#     # ===========================
#
#     lambda_index <- which(LAMBDA == lambda_names)
#     lambda <- tuning_params_mat["lambda", LAMBDA][[1]]
#
#     # iteration counter
#     k <- 0
#
#     # to enter while loop
#     converged <- FALSE
#
#     while (!converged && k < maxit) {
#       Theta_init <- c(drop(beta_init), eta_init, sigma2_init)
#
#       # observation weights
#       di <- 1 + eta_init * (Lambda - 1)
#       wi <- (1 / sigma2_init) * (1 / di)
#
#       # fit beta
#       beta_next_fit <- glmnet::glmnet(
#         x = utx,
#         y = uty,
#         family = "gaussian",
#         weights = wi,
#         alpha = alpha,
#         penalty.factor = c(0, rep(1, p)),
#         standardize = FALSE,
#         intercept = FALSE,
#         lambda = c(.Machine$double.xmax, lambda),
#         thresh = thresh_glmnet
#       )
#
#       beta_next <- beta_next_fit$beta[, 2, drop = FALSE]
#
#       # fit eta
#       eta_next <- optim(
#         par = eta_init,
#         fn = fr_eta,
#         gr = grr_eta,
#         method = "L-BFGS-B",
#         control = list(fnscale = 1),
#         lower = 0.01,
#         upper = 0.99,
#         sigma2 = sigma2_init,
#         beta = beta_next,
#         eigenvalues = Lambda,
#         x = utx,
#         y = uty,
#         nt = n
#       )$par
#
#       # fit sigma (closed form)
#       sigma2_next <- (1 / n) * sum(((uty - utx %*% beta_next)^2) / (1 + eta_next * (Lambda - 1)))
#
#       k <- k + 1
#       # print(k)
#       Theta_next <- c(drop(beta_next), eta_next, sigma2_next)
#
#       converged <- crossprod(Theta_next - Theta_init) < epsilon
#       # converged <- max(abs(Theta_next - Theta_init) / (1 + abs(Theta_next))) <
#
#       beta_init <- beta_next
#       eta_init <- eta_next
#       sigma2_init <- sigma2_next
#
#       # message(sprintf("lambda = %s, iter = %f, l2 norm squared of Theta: %f \n log-lik: %f", LAMBDA, k, crossprod(Theta_next - Theta_init),
#       #                 log_lik(eta = eta_next, sigma2 = sigma2_next, beta = beta_next, eigenvalues = Lambda,x = utx, y = uty, nt = n)))
#     }
#
#     if (!converged) message(sprintf("algorithm did not converge for %s", LAMBDA))
#     # converged observation weights
#     di <- 1 + eta_next * (Lambda - 1)
#     wi <- (1 / sigma2_next) * (1 / di)
#
#     # nulldev <- drop(crossprod(uty - utx[,1, drop = F] %*% beta_next[1]))
#     # deviance <- drop(crossprod(uty - utx %*% beta_next))
#     # devRatio <- drop(1 - deviance/nulldev)
#
#     deviance <- log_lik(
#       eta = eta_next, sigma2 = sigma2_next, beta = beta_next,
#       eigenvalues = Lambda, x = utx, y = uty, nt = n
#     )
#
#
#     # the minus 1 is because our intercept is actually the first coefficient
#     # that shows up in the glmnet solution.
#     df <- length(glmnet::nonzeroCoef(beta_next)) - 1
#     print(df)
#     bic_lambda <- bic(
#       eta = eta_next, sigma2 = sigma2_next, beta = beta_next,
#       eigenvalues = Lambda, x = utx, y = uty, nt = n,
#       c = an, df_lambda = df
#     )
#
#     kkt_lambda <- kkt_check(
#       eta = eta_next, sigma2 = sigma2_next, beta = beta_next,
#       eigenvalues = Lambda, x = utx, y = uty, nt = n,
#       lambda = lambda, tol.kkt = tol.kkt
#     )
#
#
#     out_print[LAMBDA, ] <- c(
#       if (df == 0) 0 else df,
#       # devRatio,
#       deviance,
#       lambda,
#       bic_lambda,
#       kkt_lambda, converged
#     )
#
#     coefficient_mat[, LAMBDA] <- Theta_next
#
#
#     # prediction of random effects
#     # bi <- drop(eta_next * Phi %*% (y - x %*% beta_next)) / di
#     D_tilde_inv <- diag(1 / di)
#     V_inv <- U %*% D_tilde_inv %*% t(U)
#
#     bi <- drop(solve((1 / eta_next) * Phi_inv + V_inv) %*% V_inv %*% (y - x %*% beta_next))
#
#     # predicted values
#     yi_hat <- drop(x %*% beta_next) + bi
#
#     # fitted values
#     xbhat <- yi_hat - bi
#
#     # residuals
#     ri <- drop(y) - yi_hat
#
#     # bi <- drop(eta_next * Phi %*% (uty - utx %*% beta_next)) / di
#     # qqnorm(bi)
#     # abline(a = 0, b = 1, col = "red")
#     # plot(density(bi))
#
#     randomeff_mat[, LAMBDA] <- bi
#     fitted_mat[, LAMBDA] <- xbhat
#     predicted_mat[, LAMBDA] <- yi_hat
#     resid_mat[, LAMBDA] <- ri
#
#     deviance_change <- abs((out_print[lambda_index, "Deviance"] - out_print[lambda_index - 1, "Deviance"]) / out_print[lambda_index, "Deviance"])
#     message(sprintf("Deviance change = %.6f", deviance_change))
#
#     # this check: length(deviance_change) > 0 is for the first lambda since deviance_change returns numeric(0)
#     if (length(deviance_change) > 0) {
#       if (deviance_change < fdev | df > 50) break
#     }
#
#     pb$tick()
#   }
#
#
#   lambda_min <- out_print[which.min(out_print[, "BIC"]), "Lambda"]
#   id_min <- names(which(out_print[, "Lambda"] == lambda_min))
#
#   # if there is early stopping due to fdev, remove NAs
#   out_print <- out_print[complete.cases(out_print), ]
#
#   # get names of lambdas for which a solution was obtained
#   lambdas_fit <- rownames(out_print)
#
#   out <- list(
#     result = out_print,
#     x = x,
#     y = y,
#     coef = coefficient_mat[, lambdas_fit, drop = F],
#     b0 = coefficient_mat["beta0", lambdas_fit],
#     beta = as(coefficient_mat[colnames(x)[-1], lambdas_fit, drop = FALSE], "dgCMatrix"),
#     eta = coefficient_mat["eta", lambdas_fit, drop = FALSE],
#     sigma2 = coefficient_mat["sigma2", lambdas_fit, drop = FALSE],
#     nlambda = length(lambdas_fit),
#     randomeff = randomeff_mat[, lambdas_fit, drop = FALSE],
#     fitted = fitted_mat[, lambdas_fit, drop = FALSE],
#     predicted = predicted_mat[, lambdas_fit, drop = FALSE],
#     residuals = resid_mat[, lambdas_fit, drop = FALSE],
#     cov_names = colnames(x),
#     lambda_min = id_min,
#     lambda_min_value = lambda_min
#   )
#   # beta = beta_next,
#   # eta = eta_next,
#   # sigma2 = sigma2_next,
#   # df = nonzeroCoef(coef(beta_next_fit)),
#   # active = coef(beta_next_fit)[nonzeroCoef(coef(beta_next_fit)),, drop = F])
#   out$call <- this.call
#   class(out) <- "ggmix"
#   return(out)
# }


#' Fit a linear mixed model with lasso or elasticnet regularization
#'
#' @description By the time the user gets to this point in the function they
#'   should have already computed the singular values and the eigenvectors. this
#'   is to allow flexibility with this software, so that if for example its
#'   easier to those things elsewhere they can, else they will be computed in a
#'   wrapper function
#'
#' @seealso \code{\link[glmnet]{glmnet}}
#'
# w_svd <- svd(w)

# this is a N_T x N_T matrix
# U <- w_svd$u

# we dived by p-1 because thats how the matrix was standardized
# Lambda <- w_svd$d^2 / (p-1)
# Lambda[Lambda<1e-5] <- 1e-5

# This part shows the equivalence between the eigenvectors and eigenvalues
# from the kinship matrix vs. the SNP matrix used to construct the kinship
# U_w=U
# plot(U_w[,1],U_w[,2])
# phi_eigen <- eigen(phi)
# U_kinship <- phi_eigen$vectors
# plot(U_kinship[,1], U_w[,1])
# abline(a=0, b=-1)
# dim(U)
# vector of length N_T
# Lambda <- phi_eigen$values
# plot(Lambda, w_svd$d^2 / (p-1))
# abline(a=0,b=1,col="red")
# all.equal(Lambda, w_svd$d^2 / (p-1))

# lowrank <- function(x, y, d, u,
#                     # w,
#                     # an = log(log(n)) * log(p),
#                     lambda = NULL,
#                     lambda_min_ratio = ifelse(n < p, 0.01, 0.0001),
#                     nlambda = 100,
#                     eta_init = 0.5,
#                     maxit = 100,
#                     fdev = 1e-5,
#                     standardize = FALSE,
#                     alpha = 1, # elastic net mixing param. 1 is lasso, 0 is ridge
#                     thresh_glmnet = 1e-8, # this is for glmnet
#                     epsilon = 1e-4 # this is for ggmix
#                     # tol.kkt = 1e-9
# ) {
#
#   # rm(list=ls())
#   # pacman::p_load(gaston)
#   # pacman::p_load(glmnet)
#   # pacman::p_load(magrittr)
#   # pacman::p_load(snpStats)
#   # source("~/git_repositories/ggmix/R/fitting.R")
#   # source("~/git_repositories/ggmix/R/functions.R")
#   # source("~/git_repositories/ggmix/R/methods.R")
#   # source("~/git_repositories/ggmix/R/plot.R")
#   # # source("~/git_repositories/ggmix/R/sim-data.R")
#   # source("~/git_repositories/ggmix/simulation/model_functions.R")
#   # dat <- make_mixed_model_not_simulator(b0 = 1, eta = 0.3, sigma2 = 2, type = "causal_400", related = TRUE)
#   #
#   # x <- dat$x
#   # y <- dat$y
#   # w <- dat$w
#   # phi <- dat$kin
#   # # lambda_min_ratio <- ifelse(n < p, 0.01, 0.0001)
#   # lambda_min_ratio <- 0.01
#   # nlambda <- 100
#   # #convergence criterion
#   # thresh_glmnet <- 1e-7
#   # epsilon <- 1e-7
#   # # tol.kkt <- 1e-5
#   # maxit <- 100
#   # eta_init <- 0.5
#   # an = log(log(1069)) * log(4000)
#   # alpha = 1
#   # fdev <- 1e-4
#   # ======================================
#
#   this.call <- match.call()
#   if (!is.matrix(x)) {
#     stop("x has to be a matrix")
#   }
#   if (any(is.na(x))) {
#     stop("Missing values in x not allowed!")
#   }
#
#   np <- dim(x)
#   if (is.null(np) | (np[2] <= 1)) {
#     stop("x should be a matrix with 2 or more columns")
#   }
#
#   n <- np[[1]]
#   p <- np[[2]]
#
#   # browser()
#   # add column of 1s to x for intercept
#   vnames <- colnames(x)
#   if (is.null(vnames)) {
#     vnames <- paste("V", seq(p), sep = "")
#     colnames(x) <- vnames
#   }
#
#   x <- cbind(beta0 = 1, x)
#
#   utx <- crossprod(u, x)
#   uty <- crossprod(u, y)
#
#   # get sequence of tuning parameters
#   lamb <- lambda_sequence(
#     x = utx,
#     y = uty,
#     eigenvalues = Lambda,
#     nlambda = nlambda,
#     lambda_min_ratio = lambda_min_ratio,
#     eta_init = eta_init,
#     epsilon = epsilon
#   )
#   # browser()
#   lambda_max <- lamb$sequence[[1]]
#
#   lamb$sequence[[1]] <- .Machine$double.xmax
#
#   # tuning_params_mat <- matrix(c(lambda_max,lamb$sequence[-1]), nrow = 1, ncol = nlambda, byrow = T)
#   tuning_params_mat <- matrix(lamb$sequence, nrow = 1, ncol = nlambda, byrow = T)
#   dimnames(tuning_params_mat)[[1]] <- list("lambda")
#   dimnames(tuning_params_mat)[[2]] <- paste0("s", seq_len(nlambda))
#   lambda_names <- dimnames(tuning_params_mat)[[2]]
#
#   coefficient_mat <- matrix(
#     nrow = p + 3,
#     ncol = nlambda,
#     dimnames = list(
#       c(
#         colnames(x),
#         "eta", "sigma2"
#       ),
#       lambda_names
#     )
#   )
#
#   out_print <- matrix(NA,
#     nrow = nlambda, ncol = 8,
#     dimnames = list(
#       lambda_names,
#       c(
#         "Df",
#         "%Dev",
#         "Deviance",
#         "Lambda",
#         "saturated_loglik",
#         "model_loglik",
#         "intercept_loglik",
#         "converged"
#       )
#     )
#   )
#
#   # randomeff_mat <- matrix(nrow = n,
#   #                         ncol = nlambda,
#   #                         dimnames = list(c(paste0("Subject",1:n)),
#   #                                         lambda_names))
#
#   # fitted_mat <- matrix(nrow = n,
#   #                      ncol = nlambda,
#   #                      dimnames = list(c(paste0("Subject",1:n)),
#   #                                      lambda_names))
#
#   # predicted_mat <- matrix(nrow = n,
#   #                         ncol = nlambda,
#   #                         dimnames = list(c(paste0("Subject",1:n)),
#   #                                         lambda_names))
#
#   # resid_mat <- matrix(nrow = n,
#   #                     ncol = nlambda,
#   #                     dimnames = list(c(paste0("Subject",1:n)),
#   #                                     lambda_names))
#
#
#
#   pb <- progress::progress_bar$new(
#     format = "  fitting over all tuning parameters [:bar] :percent eta: :eta",
#     total = nlambda, clear = FALSE, width = 90
#   )
#   pb$tick(0)
#   # pb <- progress::progress_bar$new(total = nlambda)
#   # pb$tick(0)
#
#   # this includes the intercept
#   beta_init <- matrix(0, nrow = p + 1, ncol = 1)
#
#   # initial value for eta from lambda_sequence results
#   # this gives issues if lamb sequence returns eta at the boundary and
#   # fails kkt check
#   # eta_init <- lamb$eta
#
#   # closed form solution value for sigma^2
#   # sigma2_init <- (1 / n) * sum(((uty - utx %*% beta_init) ^ 2) / (1 + eta_init * (Lambda - 1)))
#
#   Lambda <- d
#
#   sigma2_init <- sigma2(
#     n = n,
#     x = utx,
#     y = uty,
#     eta = eta_init,
#     beta = beta_init,
#     eigenvalues = Lambda
#   )
#
#   # pb = txtProgressBar(min = 0, max = nlambda, initial = 0)
#
#   for (LAMBDA in lambda_names) {
#
#     # LAMBDA <- "s1"
#     # ===========================
#
#     lambda_index <- which(LAMBDA == lambda_names)
#     lambda <- tuning_params_mat["lambda", LAMBDA][[1]]
#
#     # iteration counter
#     k <- 0
#
#     # to enter while loop
#     converged <- FALSE
#
#     while (!converged && k < maxit) {
#       Theta_init <- c(as.vector(beta_init), eta_init, sigma2_init)
#
#       # observation weights
#       di <- 1 + eta_init * (Lambda - 1)
#       wi <- (1 / sigma2_init) * (1 / di)
#
#       # fit beta
#       beta_next_fit <- glmnet::glmnet(
#         x = utx,
#         y = uty,
#         family = "gaussian",
#         weights = wi,
#         alpha = alpha,
#         penalty.factor = c(0, rep(1, p)),
#         standardize = FALSE,
#         intercept = FALSE,
#         lambda = c(.Machine$double.xmax, lambda),
#         thresh = thresh_glmnet
#       )
#
#       # print(coef(beta_next_fit)[nonzeroCoef(coef(beta_next_fit)),])
#       beta_next <- beta_next_fit$beta[, 2, drop = FALSE]
#
#       # fit eta
#       eta_next <- optim(
#         par = eta_init,
#         fn = fr_eta,
#         gr = grr_eta,
#         method = "L-BFGS-B",
#         control = list(fnscale = 1),
#         lower = 0.1,
#         upper = 0.90,
#         sigma2 = sigma2_init,
#         beta = beta_next,
#         eigenvalues = Lambda,
#         x = utx,
#         y = uty,
#         nt = n
#       )$par
#
#       # fit sigma (closed form)
#       # sigma2_next <- (1 / n) * sum(((uty - utx %*% beta_next) ^ 2) / (1 + eta_next * (Lambda - 1)))
#       sigma2_next <- sigma2(n = n, x = utx, y = uty, beta = beta_next, eta = eta_next, eigenvalues = Lambda)
#
#       k <- k + 1
#
#       Theta_next <- c(as.vector(beta_next), eta_next, sigma2_next)
#
#       converged <- (crossprod(Theta_next - Theta_init) < epsilon)[1, 1]
#
#       beta_init <- beta_next
#       eta_init <- eta_next
#       sigma2_init <- sigma2_next
#     }
#
#     if (!converged) message(sprintf("algorithm did not converge for %s", LAMBDA))
#
#     # converged observation weights
#     # di <- 1 + eta_next * (Lambda - 1)
#     # wi <- (1 / sigma2_next) * (1 / di)
#
#     # a parameter for each observation
#     saturated_loglik <- log_lik(
#       eta = eta_next,
#       sigma2 = sigma2_next,
#       beta = 1,
#       eigenvalues = Lambda,
#       x = uty,
#       y = uty,
#       nt = n
#     )
#     # print(saturated_loglik)
#     # intercept only model
#     intercept_loglik <- log_lik(
#       eta = eta_next,
#       sigma2 = sigma2_next,
#       beta = beta_next[1, , drop = FALSE],
#       eigenvalues = Lambda,
#       x = utx[, 1, drop = FALSE],
#       y = uty,
#       nt = n
#     )
#     # print(intercept_loglik)
#     # model log lik
#     model_loglik <- log_lik(
#       eta = eta_next,
#       sigma2 = sigma2_next,
#       beta = beta_next,
#       eigenvalues = Lambda,
#       x = utx,
#       y = uty,
#       nt = n
#     )
#     # print(model_loglik)
#
#     deviance <- 2 * (saturated_loglik - model_loglik)
#     nulldev <- 2 * (saturated_loglik - intercept_loglik)
#     devratio <- 1 - deviance / nulldev
#
#     # the minus 1 is because our intercept is actually the first coefficient
#     # that shows up in the glmnet solution.
#     df <- length(glmnet::nonzeroCoef(beta_next)) - 1 + 2
#
#     # bic_lambda <- bic(eta = eta_next, sigma2 = sigma2_next, beta = beta_next,
#     #                   eigenvalues = Lambda, x = utx, y = uty, nt = n,
#     #                   c = an, df_lambda = df)
#
#     # kkt_lambda <- kkt_check(eta = eta_next, sigma2 = sigma2_next, beta = beta_next,
#     #                         eigenvalues = Lambda, x = utx, y = uty, nt = n,
#     #                         lambda = lambda, tol.kkt = tol.kkt)
#
#     out_print[LAMBDA, ] <- c(
#       if (df == 0) 0 else df,
#       devratio,
#       deviance,
#       lambda,
#       saturated_loglik,
#       model_loglik,
#       intercept_loglik,
#       # bic_lambda,
#       # kkt_lambda,
#       converged
#     )
#
#     coefficient_mat[, LAMBDA] <- Theta_next
#
#
#     # prediction of random effects
#     # bi <- drop(eta_next * Phi %*% (y - x %*% beta_next)) / di
#
#     # Phi Inverse (used for prediction of random effects)
#     # D_inv <- diag(1 / Lambda)
#     # Phi_inv <- U %*% D_inv %*% t(U)
#
#     # D_tilde_inv <- diag(1 / di)
#     # V_inv <- U %*% D_tilde_inv %*% t(U)
#
#     # bi <- as.vector(solve((1 / eta_next) * Phi_inv + V_inv) %*% U %*% D_inv %*% (uty - utx %*% beta_next))
#     # bi <- as.vector(U %*% diag(1 / (1/di + 1/(eta_next*Lambda))) %*% t(U) %*% U %*% D_tilde_inv %*% (uty - utx %*% beta_next))
#
#     # predicted values (this contains the intercept)
#     # yi_hat <- as.vector(x %*% beta_next) + bi
#
#     # fitted values
#     # xbhat <- yi_hat - bi
#
#     # residuals
#     # ri <- drop(y) - yi_hat
#
#     # bi <- drop(eta_next * Phi %*% (uty - utx %*% beta_next)) / di
#     # qqnorm(bi)
#     # abline(a = 0, b = 1, col = "red")
#     # plot(density(bi))
#
#     # randomeff_mat[,LAMBDA] <- bi
#     # fitted_mat[,LAMBDA] <- xbhat
#     # predicted_mat[,LAMBDA] <- yi_hat
#     # resid_mat[,LAMBDA] <- ri
#
#     deviance_change <- abs((out_print[lambda_index, "%Dev"] - out_print[lambda_index - 1, "%Dev"]) / out_print[lambda_index, "%Dev"])
#     # message(sprintf("Deviance change = %.6f", deviance_change))
#
#     # this check: length(deviance_change) > 0 is for the first lambda since deviance_change returns numeric(0)
#     if (length(deviance_change) > 0) {
#       if (deviance_change < fdev) break
#     }
#
#     pb$tick()
#
#     # setTxtProgressBar(pb,lambda_index)
#   }
#
#
#   # lambda_min <- out_print[which.min(out_print[,"BIC"]),"Lambda"]
#   # id_min <- names(which(out_print[,"Lambda"] == lambda_min))
#
#   # if there is early stopping due to fdev, remove NAs
#   out_print <- out_print[complete.cases(out_print), ]
#
#   # get names of lambdas for which a solution was obtained
#   lambdas_fit <- rownames(out_print)
#   out_print[1, "Lambda"] <- lambda_max
#   out <- list(
#     result = out_print,
#     x = x,
#     y = y,
#     utx = utx,
#     uty = uty,
#     u = u,
#     lambda = out_print[, "Lambda"],
#     eigenvalues = Lambda,
#     coef = coefficient_mat[, lambdas_fit, drop = F],
#     b0 = coefficient_mat["beta0", lambdas_fit],
#     beta = as(coefficient_mat[colnames(x)[-1], lambdas_fit, drop = FALSE], "dgCMatrix"),
#     df = out_print[lambdas_fit, "Df"],
#     eta = coefficient_mat["eta", lambdas_fit, drop = FALSE],
#     sigma2 = coefficient_mat["sigma2", lambdas_fit, drop = FALSE],
#     nlambda = length(lambdas_fit),
#     # randomeff = randomeff_mat[, lambdas_fit, drop = FALSE],
#     # fitted = fitted_mat[, lambdas_fit, drop = FALSE],
#     # predicted = predicted_mat[, lambdas_fit, drop = FALSE],
#     # residuals = resid_mat[, lambdas_fit, drop = FALSE],
#     cov_names = colnames(x) # ,
#     # lambda_min = id_min,
#     # lambda_min_value = lambda_min
#   )
#   # beta = beta_next,
#   # eta = eta_next,
#   # sigma2 = sigma2_next,
#   # df = nonzeroCoef(coef(beta_next_fit)),
#   # active = coef(beta_next_fit)[nonzeroCoef(coef(beta_next_fit)),, drop = F])
#   out$call <- this.call
#   class(out) <- "ggmix"
#   return(out)
# }
