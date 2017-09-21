# This is the main simulator file
rm(list = ls())
setwd("/home/sahir/git_repositories/penfam/simulation/")
pacman::p_load(simulator) # this file was created under simulator version 0.2.0
source("/home/sahir/git_repositories/penfam/simulation/model_functions.R")
source("/home/sahir/git_repositories/penfam/simulation/method_functions.R")
source("/home/sahir/git_repositories/penfam/simulation/eval_functions.R")

## @knitr init

name_of_simulation <- "Penfam simulation for IGES"

## @knitr main

# nsim needs to be at least 2
message("Generating data...")
sim <- new_simulation("iges-2017-penfam-v2", "IGES 2017") %>%
  generate_model(make_mixed_model, b0 = 2, sigma2 = 1,
                 # type = as.list(c(1,2,3)),
                 type = c("causal_400"),
                 eta = 0.6) %>%
  simulate_from_model(nsim = 2, index = 1:5)

message("Done generating data...")

message("Running simulations...")
sim <- run_method(sim, list(lasso, PENFAM, TWOSTEP),
                  parallel = list(socket_names = 5, libraries = c("glmnet","magrittr","MASS","progress","Matrix","coxme")))

message("Done simulations...")
sim <- sim %>% evaluate(list(rmse, tpr, fpr, r2))

save_simulation(sim)

sim <- simulator::load_simulation("iges-2017-penfam-v2")
# # sim %>% evaluate(list(sqrerr))
sim %>% evaluate(list(rmse))

# # sim %>% evaluate(list(muy))
# # warnings()
# #
# plot_eval(sim, metric_name = "eta") + panel_border()
# plot_eval(sim, metric_name = "sigma2") + panel_border()
plot_eval(sim, metric_name = "rmse") + panel_border()
# plot_eval(sim, metric_name = "r2") + panel_border()
#
plot_eval(sim, metric_name = "time") + panel_border()
plot_eval(sim, metric_name = "fpr") + panel_border() #+ ylim(c(0,1))
plot_eval(sim, metric_name = "tpr") + panel_border() #+ ylim(c(0,1))
plot_eval(sim, metric_name = "r2") + panel_border() #+ ylim(c(0,1))




sim %>% evaluate(list(tpr,fpr))
df <- as.data.frame(evals(sim))

ggplot(data = df, aes(x=fpr, y=tpr, color = Method)) + geom_point(size=3) +
  theme(legend.position = "bottom") + xlab("False positive rate") + ylab("True positive rate") + xlim(c(0.002,0.062)) +
  ylim(c(0,0.25))
ggsave(filename = "/home/sahir/Dropbox/PhD/Year4/IGES/poster/tpr_fpr.png")

sim %>% evaluate(list(rmse))
plot_eval(sim, metric_name = "rmse") + panel_border()

df <- melt(as.data.frame(evals(sim)@evals))
df


sim %>%
  tabulate_eval(metric_name = "rmse",
                format_args = list(nsmall = 3, digits = 0),
                output_type = "latex")


#
#
# head(df)

#
#
#
# sim %>% evaluate(list(rmse))
# sim %>% plot_eval(metric_name = "rmse") + panel_border()
#
#
# sim %>% evaluate(list(tpr,fpr))
# sim %>% plot_evals("fpr","tpr")
#
# sim %>% evaluate(list(sqrerr))
# sim %>% plot_eval(metric_name = "sqrerr") + panel_border()
# evals(sim)@evals
# #
# #
# # plot_eval(sim, metric_name = "muy") + panel_border()
# # plot_eval(sim, metric_name = "rmse",
# #              type = "raw", main = "p = 1")
# #
# # unname()
# # df <- as.data.frame(evals(sim))
# # head(df)
# # melt(as.data.frame(evals(sim)@evals))
# #
# # # sim_res <- sim %>%
# # #   evaluate(list(sqrerr, nnz, best_sqrerr))
# #
# # ## @knitr plots
# #
# # plot_eval_by(sim, "hisloss", varying = "prob")
# #
# # ## @knitr tables
# #
# # tabulate_eval(sim, "herloss", output_type = "markdown",
# #               format_args = list(digits = 1))
#
# # model@params$kin[1:10,1:10]
# # model@params$x %>% dim
# # draws@draws$r1.1
# # draws@draws$r1.2
# #
# dat <- make_mixed_model_not_simulator(2,.6,1,"causal_400")
# sum(dat$y>0)
# plot(dat$y)
# range(dat$y)
# plot(density(dat$y))
# # sum(draws@draws$r1.1>0)
# # sum(draws@draws$r1.2>0)
# # plot(dat$y, draws@draws$r1.1)
# #
# system.time(
# fit <- penfam(x = dat$x,
#               y = dat$y,
#               phi = dat$kin,
#               thresh_glmnet = 1e-10,
#               epsilon = 1e-5,
#               fdev = 1e-7,
#               alpha = 1,
#               tol.kkt = 1e-3,
#               nlambda = 100,
#               # an = log(log(1000)) * log(1000),
#               an = log(log(1000)),
#               # an = log(1000),
#               # lambda_min_ratio  = ifelse(model$n < model$p, 0.01, 0.001),
#               lambda_min_ratio  = 0.05,
#               eta_init = 0.5,
#               maxit = 100)
# )
# #
# plot(fit, type = "BIC")
# sum(rownames(coef(fit, s = fit$lambda_min)[nonzeroCoef(coef(fit, s = fit$lambda_min)),,drop=F]) %in% dat$causal)
# rownames(coef(fit, s = fit$lambda_min)[nonzeroCoef(coef(fit, s = fit$lambda_min)),,drop=F]) %>% length()
# sqrt(crossprod(fit$predicted[,fit$lambda_min] - dat$y))
# #
# #
# #
# #
# fitglmnet <- cv.glmnet(x = dat$x, y = dat$y, alpha = 1, standardize = F)
# plot(fitglmnet)
# sum(rownames(coef(fitglmnet)[nonzeroCoef(coef(fitglmnet)),,drop=F] ) %in% dat$causal)
# coef(fitglmnet)[nonzeroCoef(coef(fitglmnet)),,drop=F]
# coef(fitglmnet)[nonzeroCoef(coef(fitglmnet)),,drop=F] %>% dim
# sqrt(crossprod(predict(fitglmnet, newx = dat$x) - dat$y))
