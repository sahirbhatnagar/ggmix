#' @param an is the constant used in the BIC calculation. The default choice is
#'   from Wang et al. (2009)
#' @param exact If exact=FALSE, then the results from glmnet uses linear
#'   interpolation to make predictions for values of lambda that do not coincide
#'   with those used in the fitting algorithm. While this is often a good
#'   approximation, it can sometimes be a bit coarse. With exact=TRUE (default),
#'   the model is refit before predictions are made. See
#'   \code{\link[glmnet]{coef.glmnet}} for more details.
#' @references Wang, H., Li, B. and Leng, C. (2009) Shrinkage tuning parameter
#'   selection with a diverging number of parameters.J. R. Statist. Soc. B, 71,
#'   671–683.
penfam <- function(x, y, phi, lambda = NULL,
                   lambda_min_ratio  = ifelse(n < p, 0.01, 0.001),
                   nlambda = 100,
                   eta_init = 0.5,
                   maxit = 100,
                   thresh_glmnet = 1e-10, # this is for glmnet
                   epsilon = 1e-5, # this is for penfam
                   an = log(log(n)) * log(n),
                   tol.kkt = 1e-9,
                   exact = T) {

  # rm(list=ls())
  # source("~/git_repositories/penfam/R/fitting.R")
  # source("~/git_repositories/penfam/R/functions.R")
  # source("~/git_repositories/penfam/R/methods.R")
  # source("~/git_repositories/penfam/R/plot.R")
  # source("~/git_repositories/penfam/R/sim-data.R")
  # x <- X
  # y <- Y
  # phi <- Phi
  # lambda_min_ratio <- ifelse(n < p, 0.01, 0.001)
  # nlambda <- 100
  # #convergence criterion
  # thresh_glmnet <- 1e-10
  # epsilon <- 1e-5
  # tol.kkt <- 1e-5
  # maxit <- 100
  # eta_init <- 0.4
  # an = log(log(600)) * log(600)
  # lambda <- 0.10
  # exact = T
  #======================================

  this.call <- match.call()
  if (!is.matrix(x))
    stop("x has to be a matrix")
  if (any(is.na(x)))
    stop("Missing values in x not allowed!")

  np <- dim(x)
  if (is.null(np) | (np[2] <= 1))
    stop("x should be a matrix with 2 or more columns")

  n <- np[[1]]
  p <- np[[2]]

  # add column of 1s to x for intercept
  x <- cbind(1, x)
  x[1:5, 1:5]

  phi_eigen <- eigen(phi)
  # this is a N_T x N_T matrix
  U <- phi_eigen$vectors

  # dim(U)
  # vector of length N_T
  Lambda <- phi_eigen$values

  # Phi Inverse (used for prediction of random effects)
  D_inv <- diag(1 / Lambda)
  Phi_inv <- U %*% D_inv %*% t(U)

  utx <- crossprod(U, x)
  uty <- crossprod(U, y)

  # get sequence of tuning parameters
  (lamb <- lambda_sequence(x = utx, y = uty, eigenvalues = Lambda, nlambda = nlambda,
                           lambda_min_ratio = lambda_min_ratio,
                           eta_init = eta_init, epsilon = epsilon,
                           tol.kkt = tol.kkt))

  tuning_params_mat <- matrix(lamb$sequence, nrow = 1, ncol = nlambda, byrow = T)
  dimnames(tuning_params_mat)[[1]] <- list("lambda")
  dimnames(tuning_params_mat)[[2]] <- paste0("s",seq_len(nlambda))
  lambda_names <- dimnames(tuning_params_mat)[[2]]

  coefficient_mat <- matrix(nrow = p + 3,
                            ncol = nlambda,
                            dimnames = list(c(paste0("beta",0:p), "eta","sigma2"),
                                            lambda_names))

  randomeff_mat <- matrix(nrow = n,
                            ncol = nlambda,
                            dimnames = list(c(paste0("Subject",1:n)),
                                            lambda_names))

  fitted_mat <- matrix(nrow = n,
                       ncol = nlambda,
                       dimnames = list(c(paste0("Subject",1:n)),
                                       lambda_names))

  resid_mat <- matrix(nrow = n,
                       ncol = nlambda,
                       dimnames = list(c(paste0("Subject",1:n)),
                                       lambda_names))

  out_print <- matrix(NA, nrow = nlambda, ncol = 11,
                      dimnames = list(lambda_names,
                                      c("Df","%Dev","Lambda","BIC",
                                        "kkt_beta0",
                                        "kkt_eta",
                                        "kkt_sigma2",
                                        "kkt_beta_nonzero",
                                        "kkt_beta_subgr", "sum_wi", "converged")))

  pb <- progress::progress_bar$new(
    format = "  fitting over all tuning parameters [:bar] :percent eta: :eta",
    total = 100, clear = FALSE, width= 90)
  pb$tick(0)


  for (LAMBDA in lambda_names) {

    # LAMBDA <- "s1"
    # ===========================

    lambda_index <- which(LAMBDA == lambda_names)
    lambda <- tuning_params_mat["lambda", LAMBDA][[1]]


    if (lambda_index == 1) {
      # initial values for beta

      # this shows that fitting entire sequence of lambda is faster
      # tu <- microbenchmark::microbenchmark(
      # beta_init_fit <- glmnet(x = utx,
      #                         y = uty,
      #                         thresh = thresh_glmnet,
      #                         family = "gaussian",
      #                         penalty.factor = c(0, rep(1, p)),
      #                         standardize = FALSE,
      #                         intercept = FALSE,
      #                         lambda = lambda),
      # beta_init_fit <- glmnet(x = utx,
      #                         y = uty,
      #                         thresh = thresh_glmnet,
      #                         family = "gaussian",
      #                         penalty.factor = c(0, rep(1, p)),
      #                         standardize = FALSE,
      #                         intercept = FALSE,
      #                         lambda = NULL)
      # , times = 10)

      beta_init_fit <- glmnet(x = utx,
                              y = uty,
                              thresh = thresh_glmnet,
                              family = "gaussian",
                              penalty.factor = c(0, rep(1, p)),
                              standardize = FALSE,
                              intercept = FALSE,
                              lambda = NULL)

      # coef(beta_init_fit) %>% head
      # plot(beta_init_fit)
      # coef(beta_init_fit)[nonzeroCoef(coef(beta_init_fit)), , drop = F]

      # remove intercept since V1 is intercept
      beta_init <- coef(beta_init_fit, s = lambda, exact = exact,
                        x = utx,
                        y = uty,
                        thresh = thresh_glmnet,
                        family = "gaussian",
                        penalty.factor = c(0, rep(1, p)),
                        standardize = FALSE,
                        intercept = FALSE)[-1, , drop = F]
      # head(beta_init)
      # coef.approx = coef(beta_init_fit, s = lambda, exact = F)
      # coef.exact = coef(beta_init_fit, s = lambda, exact = T)
      # plot(coef.approx, coef.exact, xlab="not exact", ylab = "exact")
      # abline(a=0, b=1, col="red")
      #
      # nonzeroCoef(coef.approx) %>% length()
      # nonzeroCoef(coef.exact) %>% length()
      #
      # sum(nonzeroCoef(coef.approx) %in% nonzeroCoef(coef.exact)) / length(nonzeroCoef(coef.approx))
      # sum(nonzeroCoef(coef.exact) %in% nonzeroCoef(coef.approx)) / length(nonzeroCoef(coef.exact))

      # initial value for eta from lambda_sequence results
      eta_init <- lamb$eta

      # closed form solution value for sigma^2
      sigma2_init <- (1 / n) * sum(((uty - utx %*% beta_init) ^ 2) / (1 + eta_init * (Lambda - 1)))
    } else {

      # warm start
      beta_init <- beta_next
      eta_init <- eta_next
      sigma2_init <- sigma2_next
    }

    #iteration counter
    k <- 0

    # to enter while loop
    converged <- FALSE

    while (!converged && k < maxit) {

      Theta_init <- c(drop(beta_init), eta_init, sigma2_init)

      # observation weights
      di <- 1 + eta_init * (Lambda - 1)
      wi <- (1 / sigma2_init) * (1 / di)
      # wi <- (1 / sigma2_init) * (1 / (1 + eta_init * (Lambda - 1)))
      # length(wi)
      # plot(wi)
      # are all weights positive?
      # all(wi > 0)

      # fit beta
      beta_next_fit <- glmnet(x = utx,
                              y = uty,
                              family = "gaussian",
                              weights = wi,
                              penalty.factor = c(0, rep(1, p)),
                              standardize = FALSE,
                              intercept = FALSE,
                              lambda = NULL,
                              thresh = thresh_glmnet)
      # coef(beta_next_fit)[nonzeroCoef(coef(beta_next_fit)),, drop = F]

      beta_next <- coef(beta_next_fit, s = lambda, exact = exact,
                        x = utx,
                        y = uty,
                        family = "gaussian",
                        weights = wi,
                        penalty.factor = c(0, rep(1, p)),
                        standardize = FALSE,
                        intercept = FALSE,
                        thresh = thresh_glmnet)[-1, , drop = F]
      # plot(beta_next)

      # fit eta
      eta_next <- optim(par = eta_init,
                        fn = fr_eta,
                        gr = grr_eta,
                        method = "L-BFGS-B",
                        control = list(fnscale = 1),
                        lower = 0.01,
                        upper = 0.99,
                        sigma2 = sigma2_init,
                        beta = beta_next,
                        eigenvalues = Lambda,
                        x = utx,
                        y = uty,
                        nt = n)$par

      # fit sigma (closed form)
      sigma2_next <- (1 / n) * sum(((uty - utx %*% beta_next) ^ 2) / (1 + eta_next * (Lambda - 1)))

      k <- k + 1

      Theta_next <- c(drop(beta_next), eta_next, sigma2_next)

      converged <- crossprod(Theta_next - Theta_init) < epsilon
      # converged <- max(abs(Theta_next - Theta_init) / (1 + abs(Theta_next))) <

      beta_init <- beta_next
      eta_init <- eta_next
      sigma2_init <- sigma2_next

      # message(sprintf("lambda = %s, iter = %f, l2 norm squared of Theta: %f \n log-lik: %f", LAMBDA, k, crossprod(Theta_next - Theta_init),
      #                 log_lik(eta = eta_next, sigma2 = sigma2_next, beta = beta_next, eigenvalues = Lambda,x = utx, y = uty, nt = n)))

    }

    if (!converged) message(sprintf("algorithm did not converge for %s", LAMBDA))
    # converged observation weights
    di <- 1 + eta_next * (Lambda - 1)
    wi <- (1 / sigma2_next) * (1 / di)

    nulldev <- drop(crossprod(uty - utx[,1, drop = F] %*% beta_next[1]))
    deviance <- drop(crossprod(uty - utx %*% beta_next))
    devRatio <- drop(1 - deviance/nulldev)

    # the minus 1 is because our intercept is actually the first coefficient
    # that shows up in the glmnet solution.
    df <- length(glmnet::nonzeroCoef(beta_next)) - 1

    bic_lambda <- bic(eta = eta_next, sigma2 = sigma2_next, beta = beta_next,
                      eigenvalues = Lambda, x = utx, y = uty, nt = n,
                      c = an, df_lambda = df)

    kkt_lambda <- kkt_check(eta = eta_next, sigma2 = sigma2_next, beta = beta_next,
                            eigenvalues = Lambda, x = utx, y = uty, nt = n,
                            lambda = lambda, tol.kkt = tol.kkt)


    out_print[LAMBDA,] <- c(if (df == 0) 0 else df,
                            devRatio,
                            lambda,
                            bic_lambda,
                            kkt_lambda, sum(wi), converged)

    coefficient_mat[,LAMBDA] <- Theta_next


    # prediction of random effects
    # bi <- drop(eta_next * Phi %*% (y - x %*% beta_next)) / di
    D_tilde_inv <- diag(1 / di)
    V_inv <- U %*% D_tilde_inv %*% t(U)

    bi <- drop(solve((1 / eta_next) * Phi_inv + V_inv) %*% V_inv %*% (y - x %*% beta_next))

    # fitted values
    yi_hat <- drop(x %*% beta_next) + bi

    # residuals
    ri <- drop(y) - yi_hat

    # bi <- drop(eta_next * Phi %*% (uty - utx %*% beta_next)) / di
    # qqnorm(bi)
    # abline(a = 0, b = 1, col = "red")
    # plot(density(bi))

    randomeff_mat[,LAMBDA] <- bi
    fitted_mat[,LAMBDA] <- yi_hat
    resid_mat[,LAMBDA] <- ri

    pb$tick()
  }


  lambda_min <- out_print[which.min(out_print[,"BIC"]),"Lambda"]
  id_min <- names(which(out_print[,"Lambda"] == lambda_min))

  out <- list(result = out_print,
              x = x,
              y = y,
              coef = coefficient_mat,
              b0 = coefficient_mat["beta0",],
              beta = as(coefficient_mat[paste0("beta",1:p),, drop = FALSE],"dgCMatrix"),
              eta = coefficient_mat["eta",,drop = FALSE],
              sigma2 = coefficient_mat["sigma2",,drop = FALSE],
              nlambda = nlambda,
              randomeff = randomeff_mat,
              fitted = fitted_mat,
              residuals = resid_mat,
              cov_names = paste0("beta",1:p),
              lambda_min = id_min,
              lambda_min_value = lambda_min)
  # beta = beta_next,
  # eta = eta_next,
  # sigma2 = sigma2_next,
  # df = nonzeroCoef(coef(beta_next_fit)),
  # active = coef(beta_next_fit)[nonzeroCoef(coef(beta_next_fit)),, drop = F])
  out$call <- this.call
  class(out) <- "penfam"
  return(out)

}
