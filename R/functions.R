#' Calculate Sequence of Tuning Parameters For
#'
#' @description Function to calculate the sequence of tuning parameters based on
#'   the design matrix \code{x} and the response variable {y}. This is used in
#'   the \code{\link{shim_once}} function to calculate the tuning parameters
#'   applied to the main effects
#' @inheritParams lowrank
#' @param weights Separate penalty factors can be applied to each coefficient.
#'   This is a number that multiplies lambda to allow differential shrinkage,
#'   and can be used to apply adaptive LASSO. Can be 0 for some variables, which
#'   implies no shrinkage, and that variable is always included in the model.
#'   Default is 1 for all variables (and implicitly infinity for variables
#'   listed in exclude). Note: the penalty factors are internally rescaled to
#'   sum to nvars, and the lambda sequence will reflect this change.
#' @param lambda.factor The factor for getting the minimal lambda in lambda
#'   sequence, where \code{min(lambda) = lambda.factor * max(lambda).
#'   max(lambda)} is the smallest value of lambda for which all coefficients are
#'   zero. The default depends on the relationship between \code{N} (the number
#'   of rows in the matrix of predictors) and \code{p} (the number of
#'   predictors). If \code{N > p}, the default is \code{1e-6}, close to zero. If
#'   \code{N < p}, the default is \code{0.01}. A very small value of
#'   lambda.factor will lead to a saturated fit.
#' @param nlambda the number of lambda values - default is 100.
#' @param scale_x should the columns of x be scaled - default is FALSE
#' @param center_y should y be mean centered - default is FALSE
#' @return numeric vector of length \code{q}
#' @details The maximum lambda is calculated using the following inequality:
#'   \deqn{(N*w_j)^-1 | \sum x_ij y_i | \le \lambda_max}
#'
#'   The minimum lambda is given by lambda.factor*lambda_max. The sequence of
#'   nlambda values are decreasing from lambda_max to lambda_min on the log
#'   scale.
#'
#'   The penalty factors are internally rescaled to sum to the number of
#'   predictor variables in glmnet. Therefore, to get the correct sequence of
#'   lambdas when there are weights, this function first rescales the weights
#'   and then calclated the sequence of lambdas.
#'
#'   This formula is taken from section 2.5 of the \code{glmnet} paper in the
#'   Journal of Statistical Software (see references for details)
#'
# lambda_sequence <- function(x, y, eigenvalues, weights = NULL,
#                             lambda_min_ratio,
#                             epsilon = 1e-14,
#                             tol.kkt = 1e-9,
#                             eta_init = 0.5,
#                             nlambda = 100, scale_x = F, center_y = F) {
#
#   np <- dim(x)
#   n <- np[[1]]
#
#   # assuming the first column is the intercept, so we subtract 1
#   p <- np[[2]] - 1
#
#   # weights
#   di <- 1 + eta_init * (eigenvalues - 1)
#   # di_inverse <- diag(1 / di)
#
#   # initial value for beta0
#   beta0_init <- (sum(x[, 1] * y / di)) / (sum(x[,1] ^ 2 / di))
#
#   # this includes all other betas which are 0 by definition
#   beta_init <- as.matrix(c(beta0_init, rep(0,p)))
#   # sum is faster
#   # microbenchmark::microbenchmark(
#   #   mat = drop((t(utx0) %*% di_inverse %*% uty) / (t(utx0) %*% di_inverse %*% utx0)),
#   #     sum = (sum(utx0 * uty / di)) / (sum(utx0 ^ 2 / di)), times = 1000
#   # )
#
#   # closed form for sigma^2
#   sigma2_init <- (1 / n) * sum((y - x %*% beta_init) ^ 2 / di)
#
#   # sum version is faster
#   # mb <- microbenchmark(
#   #   mat = (1 / n) * t(uty - beta0_init * utx0) %*% di_inverse %*% (uty - beta0_init * utx0),
#   #   sum = (1 / n) * sum ((uty - beta0_init * utx0)^2 / (1 + eta_init * (eigenvalues - 1))),
#   #   times = 1000)
#   # ggplot2::autoplot(mb)
#
#   #iteration counter
#   k <- 0
#
#   # to enter while loop
#   converged <- FALSE
#
#   while (!converged) {
#
#     Theta_init <- c(beta_init, eta_init, sigma2_init)
#
#     # fit eta
#     eta_next <- optim(par = eta_init,
#                       fn = fr_eta,
#                       gr = grr_eta,
#                       method = "L-BFGS-B",
#                       control = list(fnscale = 1),
#                       lower = .01,
#                       upper = .99,
#                       sigma2 = sigma2_init,
#                       beta = beta_init,
#                       eigenvalues = eigenvalues,
#                       x = x,
#                       y = y,
#                       nt = n)$par
#
#     # weights
#     di <- 1 + eta_next * (eigenvalues - 1)
#
#     # di_inverse <- diag(1 / di)
#
#     # next value for beta0
#     beta0_next <- (sum(x[,1] * y / di)) / (sum(x[,1] ^ 2 / di))
#
#     beta_next <- as.matrix(c(beta0_next, rep(0, p)))
#
#     # closed form for sigma^2
#     sigma2_next <- (1 / n) * sum((y - x %*% beta_next) ^ 2 / di)
#
#     k <- k + 1
#
#     Theta_next <- c(beta_next, eta_next, sigma2_next)
#
#     converged <- crossprod(Theta_next - Theta_init) < epsilon
#     # converged <- max(abs(Theta_next - Theta_init) / (1 + abs(Theta_next))) < epsilon
#
#     message(sprintf("l2 norm squared of Theta_k+1 - Theta_k: %f \n log-lik: %f",
#                     crossprod(Theta_next - Theta_init),
#                     log_lik(eta = eta_next, sigma2 = sigma2_next, beta = beta_next,
#                             eigenvalues = eigenvalues,
#                             x = x, y = y, nt = n)))
#
#     beta0_init <- beta0_next
#     beta_init <- beta_next
#     eta_init <- eta_next
#     sigma2_init <- sigma2_next
#
#   }
#
#   # eta_next
#   # sigma2_next
#   di <- 1 + eta_next * (eigenvalues - 1)
#   wi <- (1 / sigma2_next) * (1 / di)
#   if (any(wi < 0)) stop("weights are negative")
#
#   # scale the weights to sum to nvars
#   # wi_scaled <- as.vector(wi) / sum(as.vector(wi)) * n
#
#   # wi_scaled <- as.vector(wi) * n
#
#   # lambda.max <- max(abs(colSums((wi * x[,-1]) * drop(y - x %*% beta_next))))
#
#   # this gives the same answer (see paper for details)
#   # we divide by sum(wi) here and not in glmnet because the sequence is determined
#   # on the log scale
#   # lambda.max <- max(abs(colSums(((1 / sum(wi_scaled)) * (wi_scaled * x[,-1]) * drop(y - x %*% beta_next)))))
#
#   lambda.max <- max(abs(colSums(((1 / sum(wi)) * (wi * x[,-1]) * drop(y - x %*% beta_next)))))
#
#   # lambda.max <- lambda.max * sum(wi)
#   # (x[,-1, drop = F]) %>% dim
#   # a <- colSums(x[,-1, drop = F]^2 * wi)
#   # b <- colSums(sweep(x[,-1, drop = F]^2, MARGIN = 1, wi, '*'))
#   # all(a == b)
#
#
#   kkt <- kkt_check(eta = eta_next, sigma2 = sigma2_next, beta = beta_next,
#                    eigenvalues = eigenvalues, x = x, y = y, nt = n,
#                    lambda = lambda.max, tol.kkt = tol.kkt)
#   # message(kkt)
#   out <- list(sequence = rev(exp(seq(log(lambda_min_ratio * lambda.max), log(lambda.max), length.out = nlambda))),
#               eta = eta_next, sigma2 = sigma2_next, beta0 = beta0_next, kkt = kkt)
#
# }

# bic <- function(eta, sigma2, beta, eigenvalues, x, y, nt, c, df_lambda) {
#
#   -2 * log_lik(eta = eta, sigma2 = sigma2, beta = beta,
#                eigenvalues = eigenvalues, x = x, y = y, nt = nt) + c * df_lambda
#
# }


#' An alternative to \code{summaryRprof()}
#'
#' \code{proftools} parses a profiling file and prints an easy-to-understand
#' table showing the most time-intensive function calls.
#'
#' Line numbers are included if \code{Rprof()} was run with
#' \code{line.numbering=TRUE}. If it was run with \code{memory.profiling=TRUE},
#' this function will probably break.
#'
#' Below the table are printed any files identified if line numbering is true,
#' the total time recorded by \code{Rprof()}, and the "parent call".  The
#' parent call consists of the parent call stack of all the call stacks in the\
#' table. Note that this is the parent call stack of only the printed lines,
#' not of all stacks recorded by \code{Rprof()}. This makes the table easier to read and fit into the console.
#'
#' @param file A profiling file generated by \code{Rprof()}
#' @param lines The number of lines (call stacks) you want returned. Lines are
#' printed from most time-intensive to least.
# proftable <- function(file, lines = 10) {
#   profdata <- readLines(file)
#   interval <- as.numeric(strsplit(profdata[1L], "=")[[1L]][2L]) / 1e+06
#   filelines <- grep("#File", profdata)
#   files <- profdata[filelines]
#   profdata <- profdata[-c(1, filelines)]
#   total.time <- interval * length(profdata)
#   ncalls <- length(profdata)
#   profdata <- gsub("\\\"| $", "", profdata)
#   calls <- lapply(profdata, function(x) rev(unlist(strsplit(x, " "))))
#   stacktable <- as.data.frame(table(sapply(calls, function(x) paste(x, collapse = " > "))) / ncalls * 100, stringsAsFactors = FALSE)
#   stacktable <- stacktable[order(stacktable$Freq[], decreasing = TRUE), 2:1]
#   colnames(stacktable) <- c("PctTime", "Call")
#   stacktable <- head(stacktable, lines)
#   shortcalls = strsplit(stacktable$Call, " > ")
#   shortcalls.len <- range(sapply(shortcalls, length))
#   parent.call <- unlist(lapply(seq(shortcalls.len[1]), function(i) Reduce(intersect, lapply(shortcalls,"[[", i))))
#   shortcalls <- lapply(shortcalls, function(x) setdiff(x, parent.call))
#   stacktable$Call = sapply(shortcalls, function(x) paste(x, collapse = " > "))
#   if (length(parent.call) > 0) {
#     parent.call <- paste(paste(parent.call, collapse = " > "), "> ...")
#   } else {
#     parent.call <- "None"
#   }
#   frac <- sum(stacktable$PctTime)
#   attr(stacktable, "total.time") <- total.time
#   attr(stacktable, "parent.call") <- parent.call
#   attr(stacktable, "files") <- files
#   attr(stacktable, "total.pct.time") <- frac
#   print(stacktable, row.names=FALSE, right=FALSE, digits=3)
#   if(length(files) > 0) {
#     cat("\n")
#     cat(paste(files, collapse="\n"))
#     cat("\n")
#   }
#   cat(paste("\nParent Call:", parent.call))
#   cat(paste("\n\nTotal Time:", total.time, "seconds\n"))
#   cat(paste0("Percent of run time represented: ", format(frac, digits=3)), "%")
#
#   invisible(stacktable)
# }
