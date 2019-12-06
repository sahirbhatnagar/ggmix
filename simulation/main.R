######################################
# R Source code file for simulation of ggmix paper
# Notes:
# This is the main simulator file
# I ran this on tmux by copy pasting into an rsession on hydra.
# parallel doesnt work in rstudio
# use this code to save the results. turns out you cant run simulator in parallel "by hand"
# you need to load the simulator object prior to adding new simulations to the object.
#
# Since the simulations take a long time to run, I save the results to a .rds file
# and then create figures based on that
# Author: Sahir Bhatnagar
# Created: 2018
# Updated: Dec 4, 2019
# To address referee comments. Need to add simulation where n is approx equal to p
#####################################



# This is the main simulator file
# rm(list = ls())

# setwd("/home/sahir/git_repositories/ggmix/simulation/")
pacman::p_load(simulator) # this file was created under simulator version 0.2.0
source("/home/sahir/git_repositories/ggmix/simulation/model_functions.R")
source("/home/sahir/git_repositories/ggmix/simulation/method_functions.R")
source("/home/sahir/git_repositories/ggmix/simulation/eval_functions.R")

# source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/ggmix/simulation/model_functions.R")
# source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/ggmix/simulation/method_functions.R")
# source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/ggmix/simulation/eval_functions.R")

## @knitr init

# name_of_simulation <- "thesis-ggmix-july1" this had only one eta value
# name_of_simulation <- "thesis-ggmix-july3" # this has more than 1 eta value 0.1,0.2,0.3,0.4
# name_of_simulation <- "thesis-ggmix-july12" # this has percent causal 0,0.01, and eta=0.1, 0.5
# name_of_simulation <- "ggmix-mar5" # this has percent causal 0,0.01, and eta=0.1, 0.5
# name_of_simulation <- "ggmix-apr29" # this has train/test split 50/50 split
# name_of_simulation <- "ggmix-may7" # this has train/test/validation split 60/20/20 split
# name_of_simulation <- "ggmix-jul8" # this has train/test/validation split 60/20/20 split but with dfmax not specified in ggmix
# name_of_simulation <- "ggmix-jul10" # this has train/validation split 80/20 split and HDBIC for ggmix ,with lm on active snps for prediction, and TPR at fixed FPR of 5%

## @knitr main


# I used this on Decemerb 4th 2019, because I wanted to re-use the simulated data,
# but apply the three methods with the TPR calculated at a fixed TPR of 5%
# For the n=p simulation, I need to create a new simulatior object. See below.
sim <- simulator::load_simulation(name_of_simulation, dir = "simulation/")
# sim <- sim %>%
#   run_method(list(twostepYVCCV, lassoCV, ggmixedHDBIC),
#              # run_method(list(lasso, ggmixed, twostepYVC, lassoNOPC),
#              parallel = list(socket_names = 40,
#                              libraries = c("glmnet","magrittr","MASS","Matrix","coxme","gaston","ggmix","popkin","bnpsd"))) %>%
#   evaluate(list(modelerror, prederror,tpr, fpr, nactive, eta, sigma2, tprFPR5, nactiveFPR5,
#                 correct_sparsity,mse, errorvariance, estimationerror))
# save_simulation(sim)
# as.data.frame(evals(sim))
# ls()
# # nsim needs to be at least 2
#
# # sim <- new_simulation(name_of_simulation, "Thesis 2018", dir = "simulation/") %>%
# #   generate_model(make_mixed_model_SSC, b0 = 0, sigma2 = list(2, 4),
# #                  eta = list(0.1, 0.5),
# #                  percent_causal = 1,
# #                  beta_mean = list(0.2, 0.6, 1),
# #                  percent_overlap = list("0","100"),
# #                  vary_along = c("sigma2","eta","beta_mean","percent_overlap")) %>%
# #   simulate_from_model(nsim = 4, index = 1:50) %>%
# #   run_method(list(lasso, ggmix, twostep),
# #              parallel = list(socket_names = 35,
# #                              libraries = c("glmnet","magrittr","MASS","progress","Matrix","coxme","gaston")))
# #
# # sim <- sim %>% evaluate(list(modelerror, tpr, fpr, nactive, eta, sigma2))
# # save_simulation(sim)
# # sim <- new_simulation(name_of_simulation, "jul-10", dir = "simulation/") %>%
#
# sim <- new_simulation(name_of_simulation, "dec-4", dir = "simulation/") %>%
#   generate_model(make_ADmixed_model_train_validate,
#                  b0 = 1,
#                  sigma2 = 1,
#                  beta_mean = 1,
#                  k = 10,
#                  s = 0.5,
#                  Fst = 0.1,
#                  geography = "1d",
#                  n = 1000, # 80/20 split
#                  p_design = 5000,
#                  p_kinship = 10000,
#                  percent_causal = list(0, 0.01),
#                  percent_overlap = list("0","100"),
#                  eta = list(0.1, 0.3),
#                  vary_along = c("percent_overlap","percent_causal","eta")
#
#                  # n = 1000, # 80/20 split
#                  # p_design = 5000,
#                  # p_kinship = 10000,
#                  # percent_causal = 0.01,
#                  # percent_overlap = "0",
#                  # eta = 0.1#,
#                  # vary_along = c("percent_overlap","percent_causal","eta")
#
#                  # n = 500, # 60/20/20 split
#                  # p_design = 500,
#                  # p_kinship = 1000
#                  # eta = 0.50,
#                  # geography = "1d",
#                  # percent_causal = 0,
#                  # percent_overlap = "0",
#                  # n = 800, # 50% train, 50% test
#                  # p_design = 500,
#                  # p_kinship = 1000
#
#   ) %>%
#   # simulate_from_model(nsim = 5, index = 1:40) %>%
#   simulate_from_model(nsim = 5, index = 1:40) %>%
#   run_method(list(twostepYVCCV, lassoCV, ggmixedHDBIC),
#              # run_method(list(lasso, ggmixed, twostepYVC, lassoNOPC),
#              parallel = list(socket_names = 40,
#                              libraries = c("glmnet","magrittr","MASS","Matrix","coxme","gaston","ggmix","popkin","bnpsd"))) %>%
#   evaluate(list(modelerror, prederror,tpr, fpr, nactive, eta, sigma2, tprFPR5, nactiveFPR5,
#                 correct_sparsity,mse, errorvariance, estimationerror))
# save_simulation(sim)
# as.data.frame(evals(sim))
# ls()



# n=p simulation ----------------------------------------------------------
# December 5th, 2019
# this has train/validation split 80/20 split and HDBIC for ggmix ,with lm on active snps for prediction, and n = p simulation for referee comments
name_of_simulation <- "ggmix-dec5"
sim <- new_simulation(name_of_simulation, "dec-5", dir = "simulation/") %>%
  generate_model(make_ADmixed_model_train_validate,
                 b0 = 1,
                 sigma2 = 1,
                 beta_mean = 1,
                 k = 10,
                 s = 0.5,
                 Fst = 0.1,
                 geography = "1d",
                 n = 1000, # 80/20 split
                 p_design = 5000,
                 p_kinship = 1000,
                 percent_causal = list(0, 0.01),
                 percent_overlap = list("0","100"),
                 eta = list(0.1, 0.3),
                 vary_along = c("percent_overlap","percent_causal","eta")

                 # n = 1000, # 80/20 split
                 # p_design = 5000,
                 # p_kinship = 10000,
                 # percent_causal = 0.01,
                 # percent_overlap = "0",
                 # eta = 0.1#,
                 # vary_along = c("percent_overlap","percent_causal","eta")

                 # n = 500, # 60/20/20 split
                 # p_design = 500,
                 # p_kinship = 1000
                 # eta = 0.50,
                 # geography = "1d",
                 # percent_causal = 0,
                 # percent_overlap = "0",
                 # n = 800, # 50% train, 50% test
                 # p_design = 500,
                 # p_kinship = 1000

  ) %>%
  # simulate_from_model(nsim = 5, index = 1:40) %>%
  simulate_from_model(nsim = 5, index = 1:40) %>%
  run_method(list(twostepYVCCV, lassoCV, ggmixedHDBIC),
             # run_method(list(lasso, ggmixed, twostepYVC, lassoNOPC),
             parallel = list(socket_names = 40,
                             libraries = c("glmnet","magrittr","MASS","Matrix","coxme","gaston","ggmix","popkin","bnpsd"))) %>%
  evaluate(list(modelerror, prederror,tpr, fpr, nactive, eta, sigma2, tprFPR5, nactiveFPR5,
                correct_sparsity,mse, errorvariance, estimationerror))
save_simulation(sim)
as.data.frame(evals(sim))
ls()


sim <- simulator::load_simulation(name_of_simulation, dir = "simulation/")
df <- as.data.frame(evals(sim))
saveRDS(df, file = "simulation/simulation_results/dec_5_2019_results.rds") # this has train/validate split only. and doing re-fit for MSE and HDBIC for ggmix and TPR at FPR 5 % and n=1000=kinship=1000

# make sure most recent version of ggmix is installed for the parallel code to work!!!!#!#@$!@$!@$
# remotes::install_github('sahirbhatnagar/ggmix', ref = "validate")

# ggmixed@method(draw = draws(sim)@draws$r1.2)
#
# tt <- ggmixed@method(draw = draws(sim)@draws$r1.1)
# draws(sim)@draws$r1.2$xtrain
# ls()
# tt$causal
# tt$nonzero
# sim <- sim %>% run_method(list(lasso, ggmixed))
# tp <- simulator::load_draws(dir = "simulation/", model_name = "ggmix_05_07_2019")
# sim <- sim %>% run_method(list(lasso, ggmixed))#,#, twostep, twostepY),
#              parallel = list(socket_names = 8,
#                              libraries = c("glmnet","magrittr","MASS","Matrix","coxme","gaston","ggmix","popkin","bnpsd"))

# save_simulation(sim)
# sim <- sim %>% evaluate(list(modelerror, prederror,tpr, fpr, nactive, eta, sigma2,
#                              correct_sparsity,mse, errorvariance))
# save_simulation(sim)
# as.data.frame(evals(sim))
# ls()
#
# sim <- sim %>% run_method(list(twostepYVC),
#                   parallel = list(socket_names = 40,
#                                   libraries = c("glmnet","magrittr","MASS","Matrix","coxme","gaston","ggmix","popkin","bnpsd"))) %>%
#   evaluate(list(modelerror, prederror,tpr, fpr, nactive, eta, sigma2,
#                 correct_sparsity,mse, errorvariance))
#
# sim <- sim %>% run_method(list(lasso, lassoNOPC, ggmixed),
#                           parallel = list(socket_names = 40,
#                                           libraries = c("glmnet","magrittr","MASS","Matrix","coxme","gaston","ggmix","popkin","bnpsd"))) %>%
#   evaluate(list(modelerror, prederror,tpr, fpr, nactive, eta, sigma2,
#                 correct_sparsity,mse, errorvariance))
#
# save_simulation(sim)
# ls()
#
# sim <- load_simulation(name = name_of_simulation, dir = "/home/sahir/git_repositories/ggmix/simulation/")
# sim %>% subset_simulation(methods = c("lasso","ggmix")) %>% plot_eval("mse")
# sim %>% plot_eval("tprFPR5")
# sim %>% plot_eval("nactiveFPR5")
# sim2 <- sim %>% subset_simulation(percent_overlap == "100" & percent_causal == 0.01 & eta == 0.10)
# name_of_simulation <- "ggmix-jul8"
# label = "jul-8"
# sim2 <- simulator::relabel(sim2, label = label)
# sim2 <- simulator::rename(sim2, name = name_of_simulation)
# save_simulation(sim2)
# sim2

# name_of_simulation <- "ggmix-jul8"
# sim2 <- load_simulation(name = name_of_simulation, dir = "/home/sahir/git_repositories/ggmix/simulation/")
#
# sim2 <- sim2 %>% run_method(list(lasso, lassoNOPC, ggmixed),
#                           parallel = list(socket_names = 40,
#                                           libraries = c("glmnet","magrittr","MASS","Matrix","coxme","gaston","ggmix","popkin","bnpsd"))) %>%
#   evaluate(list(modelerror, prederror,tpr, fpr, nactive, eta, sigma2,
#                 correct_sparsity,mse, errorvariance))
#
# save_simulation(sim2)
# ls()

rm(sim)
sim <- load_simulation(name = name_of_simulation, dir = "/home/sahir/git_repositories/ggmix/simulation/")
# sim <- sim %>%
#   run_method(list(ggmixedHDBIC),
#              parallel = list(socket_names = 40,
#                              libraries = c("glmnet","magrittr","MASS","Matrix","coxme","gaston","ggmix","popkin","bnpsd"))) %>%
#   evaluate(list(modelerror, prederror,tpr, fpr, nactive, eta, sigma2,
#                 correct_sparsity,mse, errorvariance))
# save_simulation(sim)
# as.data.frame(evals(sim))
# ls()

# sim %>% subset_simulation(methods = c("lasso","ggmix","lassoNOPC","ggmixHDBIC","lassoCV","twostepYVCCV")) %>% plot_eval("mse")
p1 <- sim %>% plot_eval("mse")

sim %>% subset_simulation(methods = c("lasso","ggmix","lassoNOPC","ggmixHDBIC")) %>% plot_eval("fpr")
sim %>% subset_simulation(methods = c("lasso","ggmix","lassoNOPC","ggmixHDBIC")) %>% plot_eval("fpr")
sim %>% subset_simulation(methods = c("lasso","ggmix","lassoNOPC","ggmixHDBIC","lassoCV")) %>% plot_eval("correct_sparsity")
sim %>% plot_eval("tpr")
sim %>% plot_eval("fpr")
sim %>% subset_simulation(methods = c("ggmix","ggmixHDBIC","lassoCV")) %>% plot_eval("errorvar")
sim %>% plot_eval("eta")
sim %>% plot_eval("time")
sim %>% plot_eval("nactive")
sim %>% tabulate_eval("nactive")
sim %>% plot_eval("me")
sim %>% subset_simulation(methods = c("lasso","ggmix","lassoNOPC","ggmixHDBIC")) %>% plot_eval("prederror")
sim %>% subset_simulation(methods = c("ggmix","ggmixHDBIC")) %>% plot_eval("eta")
sim %>% subset_simulation(methods = c("ggmix","ggmixHDBIC")) %>% plot_eval("errorvar")
sim %>% subset_simulation(methods = c("ggmix","ggmixHDBIC")) %>% plot_eval("sigma2")
sim %>% plot_eval("estimationerror")

plot_eval(sim, "tpr")
plot_eval(sim, "fpr")
plot_eval(sim, "correct_sparsity")

# save results ------------------------------------------------------------
sim <- load_simulation(name = name_of_simulation, dir = "/home/sahir/git_repositories/ggmix/simulation/")
df <- as.data.frame(evals(sim))
# sim %>% subset_simulation(methods = c("ggmix","lasso"))
# saveRDS(df, file = "simulation/simulation_results/may_02_2019_results.rds")
# saveRDS(df, file = "simulation/simulation_results/may_05_2019_results.rds") # this has lasso1se
# saveRDS(df, file = "simulation/simulation_results/may_06_2019_results.rds") # this has lasso1se + proper variance components for twostep, but im not using lasso1se
saveRDS(df, file = "simulation/simulation_results/jul_10_2019_results_v2.rds") # this has train/validate split only. and doing re-fit for MSE and HDBIC for ggmix and TPR at FPR 5 %
df %>% filter(Method=="twostepYVCCV")

simulator::tabulate_eval(sim, "mse")

# sim <- sim %>%
#   run_method(list(lasso, ggmixed, twostepY),
#              parallel = list(socket_names = 20,
#                              libraries = c("glmnet","magrittr","MASS","Matrix",
#                                            "coxme","gaston","popkin","bnpsd", "ggmix")))
# save_simulation(sim)
# sim <- sim %>%
#   evaluate(list(modelerror, prederror,tpr, fpr, nactive, eta, sigma2,
#                 correct_sparsity,mse, errorvariance))
# save_simulation(sim)
# as.data.frame(evals(sim))
# ls()

# ns <- seq(500,4000, by = 500)
# res <- vector("numeric", length = length(ns))

# for(jj in seq_along(ns)) {
# dat <- make_ADmixed_model_not_sim(b0 = 0, sigma2 = 1,
#                                   eta = 0.1,
#                                   # n = ns[jj],
#                                   n = 1000,
#                                   p_test = 1000,
#                                   beta_mean = 0.5,
#                                   # p_test = 500,
#                                   p_kinship = 10000,
#                                   geography = "ind",
#                                   # geography = "circ",
#                                   percent_causal = 0.01,
#                                   percent_overlap = "0",
#                                   # percent_overlap = "100",
#                                   k = 5, s = 0.5, Fst = 0.1)
#
# pheno_dat <- data.frame(Y = dat$y, id = paste0("ID",1:length(dat$y)))
# x1 <- cbind(rep(1, nrow(dat$x)))
# fit <- gaston::lmm.aireml(dat$y, x1, K = dat$kin)
# kin <- gaston::as.bed.matrix(dat$Xkinship)
# gaston::standardize(kin) <- "p"
# kin[1:5,1:5]
# kins <- gaston::GRM(kin)
# dim(kins)
# kins <- crossprod(dat$Xkinship)/ncol(dat$Xkinship)
# dim(kins)
# kins[1:5,1:5]
# popkin::plotPopkin(list(dat$kin, kins))
# fit <- gaston::lmm.aireml(dat$y, x1, K = kins)
#
# res[jj] <- fit$tau




# simdata -----------------------------------------------------------------
devtools::load_all()

library(ggmix)
library(glmnet)
library(gaston)
data("admixed")

## ---- ggmix ----
fit_ggmix <- ggmix(x = admixed$x, y = admixed$y, kinship = admixed$kin, verbose = 1)
bicGGMIX <- gic(fit_ggmix, an = log(length(admixed$y)))
yhat_ggmix <- predict(bicGGMIX, newx = admixed$x)
RMSE_ggmix <- l2norm(yhat_ggmix - admixed$y)

## ---- two-step ----
x1 <- cbind(rep(1, nrow(admixed$x)))
fit_lme <- gaston::lmm.aireml(Y = admixed$y, X = x1, K = admixed$kin)
gaston_resid <- admixed$y - (fit_lme$BLUP_omega + fit_lme$BLUP_beta)
fitglmnet <- glmnet::cv.glmnet(x = admixed$x, y = gaston_resid,
                               standardize = T, alpha = 1, intercept = T)
yhat_twostep <- predict(fitglmnet, newx = admixed$x, s = "lambda.min")
RMSE_twostep <- l2norm(yhat_twostep - admixed$y)

## ---- lasso ----
fit_glmnet <- cv.glmnet(x = admixed$x_lasso, y = admixed$y,
                        alpha = 1, standardize = T,
                        penalty.factor = c(rep(1, ncol(admixed$x)), rep(0,10)))
# extract only betas for SNPs
betas_lasso <- coef(fit_glmnet, s = "lambda.min")[1:(ncol(admixed$x)+1), , drop = F]
yhat_lasso <- cbind(1, admixed$x) %*% betas_lasso
RMSE_lasso <- l2norm(yhat_lasso - admixed$y)

# karim data --------------------------------------------------------------






pacman::p_load(gaston)
data("karim")
Phi <- 2 * karim$kin1
P = MASS::mvrnorm(1, rep(0,600), karim$s.g * Phi)
y <- 1 + karim$G %*% karim$b + P + rnorm(600,0,karim$s.e)

pheno_dat <- data.frame(Y = y, id = paste0("ID",1:length(y)))
x1 <- cbind(rep(1, nrow(karim$G)))
fit <- gaston::lmm.aireml(y, x1, K = Phi)

gaston_resid <- y - (fit$BLUP_omega + fit$BLUP_beta)
hist(gaston_resid)
fitglmnet <- glmnet::cv.glmnet(x = karim$G,
                               y = gaston_resid,
                               standardize = T, alpha = 1, intercept = T)
plot(fitglmnet)
pacman::p_load(magrittr)
head(coef(fitglmnet))
yhat <- predict(fitglmnet, newx = karim$G, s = "lambda.min")
plot(yhat, y)
karim$s.g
karim$s.e
fit_2$tau
fit_2$sigma2
karim$h.tot

# need an ID variable
dat <- data.frame(Y = y, x=1, id = 1:600)

# provide the kinship matrix
gfit1 <- coxme::lmekin(Y ~ x + (1|id), data=dat, varlist=Phi)
gfit1


karim$h.tot



fit <- ggmix(x = karim$G,
             y = y,
             kinship = Phi,
             verbose = 2)
hdbic <- gic(fit, an = log(600))
help(gic)
plot(hdbic)
dev.off()

# error variance
coef(hdbic, type = "nonzero")["sigma2",] * (1-coef(hdbic, type = "nonzero")["eta",])

#kinship variance
coef(hdbic, type = "nonzero")["sigma2",] * (coef(hdbic, type = "nonzero")["eta",])



# }
# plot(ns, res)

# plot(dat$y - (fit$BLUP_omega + fit$BLUP_beta),
#      newy)
# all.equal(dat$y - (fit$BLUP_omega + fit$BLUP_beta),
#      newy)
# abline(a=0,b=1)
gaston_resid <- dat$y - (fit$BLUP_omega + fit$BLUP_beta)
hist(gaston_resid)
fitglmnet <- glmnet::cv.glmnet(x = dat$x, y = gaston_resid, standardize = T, alpha = 1, intercept = T)
plot(fitglmnet)



res$causal
res$not_causal %>% length()
res$Xtest %>% colnames() %>% length
all(colnames(res$Xtest) == res$not_causal)
# res$x_lasso
hist(res$y)
res$kin %>% dim



# analyze results ---------------------------------------------------------

source("simulation/packages.R")
df <- readRDS("simulation/simulation_results/june_29_2018_results.rds")
df <- as.data.frame(evals(sim))
df <- df %>% separate(Model, into = c("simnames","b0","eta","Fst","geography","k","n","pkinship","pdesign","percentcausal","percentoverlap","s","sigma2"),
                      sep = "/")
DT <- as.data.table(df)

trop <- RSkittleBrewer::RSkittleBrewer("trop")

cbbPalette <- c("#8720B6","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
trop <- RSkittleBrewer::RSkittleBrewer("trop")

gg_sy <- theme(legend.position = "bottom", axis.text = element_text(size = 20),
               axis.title = element_text(size = 20), legend.text = element_text(size = 20),
               legend.title = element_text(size = 20),plot.title = element_text(size = 20) )

appender <- function(string) TeX(paste(string))

DT[percentoverlap=="percent_overlap_0", p_overlap := 0]
DT[percentoverlap=="percent_overlap_100", p_overlap := 100]
DT[, table(percentoverlap, p_overlap)]

DT[geography=="geography_ind", geo := "block"]
DT[geography=="geography_circ", geo := "circular"]
DT[geography=="geography_1d", geo := "1D"]
DT[, table(geography)]

DT[, scenario:= as.numeric(as.character(stringr::str_extract_all(parameterIndex, "\\d", simplify = T)))]
DT$scenario %>% table
DT[, scen:=ifelse(scenario==1,"Strong Hierarchy",ifelse(scenario==2, "Weak Hierarchy", ifelse(scenario==3,"Interactions Only",ifelse(scenario==4, "Strong Hierarchy (Linear)", ifelse(scenario==5, "Main Effects Only", "Linear v2")))))]
DT$scen %>% table
DT[, scen:=factor(scen, levels = c("Strong Hierarchy", "Weak Hierarchy","Interactions Only","Strong Hierarchy (Linear)","Linear v2","Main Effects Only"))]
DT$scen %>% table



sim %>%
  plot_eval(metric_name = "correct_sparsity")
sim %>%
  plot_eval(metric_name = "prederror")# + ylim(c(0, 10))
sim %>%
  plot_eval(metric_name = "me")# + ylim(c(0, 10))
# simulator::draws(sim)
# we added an extra simulation scenario here.
# sim <- simulator::load_simulation(name_of_simulation, dir = "simulation/")
#
# mref <-  generate_model(make_model = make_mixed_model_SSC, b0 = 0, sigma2 = 4,
#                         eta = 0.10,
#                         percent_causal = 1,
#                         percent_overlap = list("0","100"),
#                         vary_along = "percent_overlap")
#
# dref <- simulate_from_model(mref, nsim = 4, index = 1:50)
# sim <- simulator::add(sim, mref)
# sim <- simulator::add(sim, dref)
# oref <- simulator::run_method(dref, list(lassoPCpf, PENFAM, TWOSTEP),
#                               parallel = list(socket_names = 35,
#                                               libraries = c("glmnet","magrittr","MASS","progress","Matrix","coxme","gaston")))
# sim <- simulator::add(sim, oref)
# simulator::save_simulation(sim)
# sim
#
#

## @knitr load-results

sim <- simulator::load_simulation(name_of_simulation, dir = "simulation/")


## @knitr plots

sim %>%
  subset_simulation(methods = c("ggmix","lasso")) %>%
  plot_eval(metric_name = "mse", scales = "free") + panel_border() #+
# ggplot2::geom_hline(yintercept = 0.5)

sim %>%
  subset_simulation(methods = c("ggmix")) %>%
  plot_eval(metric_name = "sigma2") + panel_border()

# appender <- function(string) TeX(paste(string))
# df <- as.data.frame(evals(sim))
# df <- df %>% separate(Model, into = c("simnames","beta0","eta","percent_causal","percent_overlap","sigma2"),
#                       sep = "/")
#
# DT <- as.data.table(df, stringsAsFactors = FALSE)
# DT[Method=="lassoPCpf", Method := "lasso"]
# DT[Method=="penfam", Method := "ggmix"]
# DT[Method=="twostep", Method := "2 step"]
# DT[, table(Method, useNA = "al")]
# DT[, Method := droplevels(Method)]
# DT[, table(Method, useNA = "al")]
# DT[, scen := case_when(percent_overlap == "percent_overlap_0" ~ "No overlap",
#                        percent_overlap == "percent_overlap_100" ~ "100% overlap")]
# DT[,table(scen, useNA = "al")]
# pacman::p_load_gh("hrbrmstr/hrbrthemes")
# pacman::p_load(ggrepel)
# # pacman::p_load(Cairo)
# pacman::p_load(extrafont)
# extrafont::loadfonts()
#
#
# df_tpr_nactive <- DT[, c("Method","scen","tpr","nactive")] %>%
#   group_by(Method, scen) %>%
#   summarise(mean.tpr = mean(tpr, na.rm = TRUE), sd.tpr = sd(tpr, na.rm = TRUE),
#             mean.nactive = mean(nactive, na.rm = TRUE), sd.nactive = sd(nactive, na.rm = TRUE))
#
# p1_tpr_nactive <- ggplot(data = df_tpr_nactive, aes(x = mean.nactive, y = mean.tpr, color = Method, label = Method)) +
#   geom_point(size = 2.1) +
#   geom_text_repel(
#     data = subset(df_tpr_nactive, mean.nactive < 40),
#     nudge_x      = 20,
#     direction    = "y",
#     hjust        = 0,
#     segment.size = 0.2
#   ) +
#   geom_text_repel(
#     data = subset(df_tpr_nactive, mean.nactive >= 40),
#     nudge_x      = 5,
#     direction    = "y",
#     hjust        = 0,
#     segment.size = 0.2
#   ) +
#   geom_errorbar(aes(ymin = mean.tpr - sd.tpr, ymax = mean.tpr + sd.tpr), size = 1.1) +
#   geom_errorbarh(aes(xmin = mean.nactive - sd.nactive, xmax = mean.nactive + sd.nactive), size = 1.1) +
#   facet_rep_wrap(~scen, scales = "free", ncol = 2,
#                  repeat.tick.labels = 'left',
#                  labeller = as_labeller(appender,
#                                         default = label_parsed)) +
#   # scale_color_brewer(palette = "Dark2")+
#   scale_color_manual(values=RColorBrewer::brewer.pal(12, "Paired")[-11], guide=guide_legend(ncol=3)) +
#   labs(x="Number of active variables", y="True Positive Rate",
#        title="True Positive Rate vs. Number of Active Variable (Mean +/- 1 SD)",
#        subtitle="Based on 200 simulations",
#        caption="") +
#   theme_ipsum_rc(axis_title_just = "bt") +
#   theme(legend.position = "right",
#         legend.text=element_text(size=14),
#         strip.text = element_text(size=14))
#
# reposition_legend(p1_mse_nactive, 'center', panel='panel-2-3')



## @knitr tpr

tabulate_eval(sim, "tpr", output_type = "markdown",
              caption = "Mean True Positive Rate (Standardar Error) over 200 simulations",
              format_args = list(nsmall = 3, digits = 0))

## @knitr fpr

tabulate_eval(sim, "fpr", output_type = "markdown",
              format_args = list(digits = 4, nsmall=3))


## @knitr nactive

tabulate_eval(sim, "nactive", output_type = "markdown",
              format_args = list(digits = 3, nsmall=3))

## @knitr model-error

tabulate_eval(sim, "me", output_type = "markdown",
              format_args = list(digits = 3, nsmall=3))

## @knitr pred-error

tabulate_eval(sim, "me", output_type = "markdown",
              format_args = list(digits = 3, nsmall=3))

## @knitr not-used


# # sim %>% evaluate(list(muy))
# # warnings()
# #
# plot_eval(sim, metric_name = "eta") + panel_border()
# plot_eval(sim, metric_name = "sigma2") + panel_border()
# plot_eval(sim, metric_name = "rmse") + panel_border()
# # plot_eval(sim, metric_name = "r2") + panel_border()
# #
# plot_eval(sim, metric_name = "time") + panel_border()
# plot_eval(sim, metric_name = "fpr") + panel_border() #+ ylim(c(0,1))
# plot_eval(sim, metric_name = "tpr") + panel_border() #+ ylim(c(0,1))
# plot_eval(sim, metric_name = "r2") + panel_border() #+ ylim(c(0,1))
#
#
#
#
# sim %>% evaluate(list(tpr,fpr))
# df <- as.data.frame(evals(sim))
# dfn <- df
# dfn$Method <- as.character(dfn$Method)
# dfn[dfn$Method=="penfam","Method"] <- "ggmix"
#
# pdf("~/Dropbox/jobs/hec/talk/tpr_fpr.pdf")#, width = 11, height = 8)
# ggplot(data = dfn[dfn$Method!="twostep",], aes(x=fpr, y=tpr, color = Method)) + geom_point(size=4) +
#   theme(legend.position = "bottom", axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"),
#         legend.text = element_text(size=16), legend.title = element_text(size=16)) +
#   xlab("False positive rate") + ylab("True positive rate") + xlim(c(0.002,0.062)) +
#   ylim(c(0,0.25))
# dev.off()
#
# pdf("~/Dropbox/jobs/hec/talk/tpr_fpr_all_causal400.pdf")#, width = 11, height = 8)
# p1 <- ggplot(data = dfn, aes(x=fpr, y=tpr, color = Method)) + geom_point(size=4) +
#   theme(legend.position = "bottom", axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"),
#         legend.text = element_text(size=16), legend.title = element_text(size=16)) +
#   xlab("False positive rate") + ylab("True positive rate") + xlim(c(0.002,0.062)) +
#   ylim(c(0,0.25))
# dev.off()
#
#
# pacman::p_load(gridExtra)
# df <- data.frame(ggmix = "26.7 (1.06)", lasso = "32.8 (0.87)", twostep = "4784.1 (0.19)")
# rownames(df) <- "mean RMSE (sd)"
# t1 <- grid.table(df)
# t1 <- tableGrob(df)
# dev.off()
#
# pdf("~/Dropbox/jobs/hec/talk/tpr_fpr_all_causal400.pdf")#, width = 11, height = 8)
# plot_grid(p1,t1, ncol = 1, rel_heights =  c(1,0.2))
# dev.off()
#
#
# png("~/Dropbox/jobs/hec/talk/tpr_fpr_all_causal400.png", width = 11, height = 8, units = "in", res = 150)
# plot_grid(p1,t1, ncol = 1, rel_heights =  c(1,0.2))
# dev.off()
#
#
#
# theme(axis.text=element_text(size=12),
#       axis.title=element_text(size=14,face="bold"))
#
# ggsave(filename = "/home/sahir/Dropbox/PhD/Year4/IGES/poster/tpr_fpr.png")
# ggsave(filename = "/home/sahir/Dropbox/jobs/hec/talk/tpr_fpr.png")
#
#
# sim %>% evaluate(list(rmse))
# plot_eval(sim, metric_name = "rmse") + panel_border()
#
# df <- melt(as.data.frame(evals(sim)@evals))
# df
#
#
# sim %>%
#   tabulate_eval(metric_name = "rmse",
#                 format_args = list(nsmall = 3, digits = 0),
#                 output_type = "latex")
#
#
# #
# #
# # head(df)
#
# #
# #
# #
# # sim %>% evaluate(list(rmse))
# # sim %>% plot_eval(metric_name = "rmse") + panel_border()
# #
# #
# # sim %>% evaluate(list(tpr,fpr))
# # sim %>% plot_evals("fpr","tpr")
# #
# # sim %>% evaluate(list(sqrerr))
# # sim %>% plot_eval(metric_name = "sqrerr") + panel_border()
# # evals(sim)@evals
# # #
# # #
# # # plot_eval(sim, metric_name = "muy") + panel_border()
# # # plot_eval(sim, metric_name = "rmse",
# # #              type = "raw", main = "p = 1")
# # #
# # # unname()
# # # df <- as.data.frame(evals(sim))
# # # head(df)
# # # melt(as.data.frame(evals(sim)@evals))
# # #
# # # # sim_res <- sim %>%
# # # #   evaluate(list(sqrerr, nnz, best_sqrerr))
# # #
# # # ## @knitr plots
# # #
# # # plot_eval_by(sim, "hisloss", varying = "prob")
# # #
# # # ## @knitr tables
# # #
# # # tabulate_eval(sim, "herloss", output_type = "markdown",
# # #               format_args = list(digits = 1))
# #
# # # model@params$kin[1:10,1:10]
# # # model@params$x %>% dim
# # # draws@draws$r1.1
# # # draws@draws$r1.2
# # #
# dat <- make_mixed_model_not_simulator(b0 = 1, eta = 0.5, sigma2 = 4, percent_causal = 1, percent_overlap = "100")
dat <- make_ADmixed_model_not_simulator(n = 200, p = 5000, ncausal = 20, k = 3, s = 0.5, Fst = 0.1, b0 = 0, beta_mean = 1,
                                        eta = 0.5, sigma2 = 4)
dat <- make_INDmixed_model_not_simulator(n = 200, p = 5000, ncausal = 10, k = 3, s = 0.5, Fst = 0.1, b0 = 0,
                                         beta_mean = 1, eta = 0.5, sigma2 = 4)
# dat$kin[1:25,1:25]
# pheno_dat <- data.frame(Y = dat$y, id = rownames(dat$kin))
pheno_dat <- data.frame(Y = dat$y, id = paste0("ID",1:length(dat$y)))
# all(names(dat$y)==rownames(dat$kin))
# hist(dat$y)
# head(pheno_dat)
# fitting with no intercept, i.e., 0 + (1|id) returns error
# fit_lme <- coxme::lmekin(Y ~ 1 + (1|id), data = pheno_dat, varlist = dat$kin)
# newy <- residuals(fit_lme)
#
# all(rownames(dat$x)==rownames(dat$kin))
x1 <- cbind(rep(1, nrow(dat$x)))
fit <- gaston::lmm.aireml(dat$y, x1, K = dat$kin)
# plot(dat$y - (fit$BLUP_omega + fit$BLUP_beta),
#      newy)
# all.equal(dat$y - (fit$BLUP_omega + fit$BLUP_beta),
#      newy)
# abline(a=0,b=1)
gaston_resid <- dat$y - (fit$BLUP_omega + fit$BLUP_beta)
hist(gaston_resid)
fitglmnet <- glmnet::cv.glmnet(x = dat$x, y = gaston_resid, standardize = T, alpha = 1, intercept = T)
plot(fitglmnet)
# yhat <- predict(fitglmnet, newx = dat$x_lasso, s = "lambda.min")
# as.numeric(sqrt(crossprod(gaston_resid - yhat)))
(nonz2step <- setdiff(rownames(coef(fitglmnet, s = "lambda.min")[nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop=F]),c("(Intercept)")))
(tpr2step <- length(intersect(nonz2step, dat$causal))/length(dat$causal))
dim(dat$x)
dim(dat$x_lasso)
fitglmnet2 <- glmnet::cv.glmnet(x = dat$x_lasso, y = dat$y, standardize = T, alpha = 1, intercept = T,
                                penalty.factor = c(rep(1, 5000), rep(0, 10)))
plot(fitglmnet2)
# yhat2 = predict(fitglmnet2, newx = dat$x_lasso, s = "lambda.min")
# as.numeric(sqrt(crossprod(dat$y - yhat2)))
(nonzlasso <- setdiff(rownames(coef(fitglmnet2, s = "lambda.min")[nonzeroCoef(coef(fitglmnet2, s = "lambda.min")),,drop=F]),c("(Intercept)","")))
(tprlasso <- length(intersect(nonzlasso, dat$causal))/length(dat$causal))

l2norm <- function(x) sqrt(sum(x^2))

l2norm(dat$x %*% dat$beta - dat$x %*% coef(fitglmnet2, s = "lambda.min")[2:(ncol(dat$x)+1),,drop=F])


coef(fitglmnet2, s = "lambda.min")[1,,drop=F]
### two step
length(nonz2step) ; tpr2step

### lasso
length(nonzlasso) ; tprlasso




all(names(dat$y)==rownames(dat$x))
all(rownames(dat$x)==rownames(dat$kin))
fit <- penfam(x = dat$x,
              y = dat$y,
              phi = dat$kin,
              # thresh_glmnet = 1e-10,
              # epsilon = 1e-5,
              # fdev = 1e-7,
              # alpha = 1,
              # tol.kkt = 1e-3,
              # nlambda = 100,
              # an = log(log(length(dat$y))) * log(length(dat$y)),
              # an = log(log(length(dat$y))),
              # an = log(length(dat$y)),
              # lambda_min_ratio  = ifelse(model$n < model$p, 0.01, 0.001),
              # lambda_min_ratio  = 0.01,
              eta_init = 0.5,
              maxit = 100)


# try re-arranigning termsof kinship to match x matrix, try hgp real data
plot(fit, "BIC")
coef(fit)[match(dat$causal, rownames(coef(fit))),]
nonzero = predict(fit, type = "nonzero", s = fit$lambda_min)
nonzero_names = setdiff(rownames(predict(fit, type = "nonzero", s = fit$lambda_min)), c("(Intercept)","eta","sigma2"))
length(intersect(nonzero_names, dat$causal))/length(dat$causal)
length(nonzero_names)

l2norm(dat$x %*% dat$beta - dat$x %*% coef(fit, s = fit$lambda_min)[2:(ncol(dat$x)+1),,drop=F])
#
#
#
# # dat <- make_mixed_model_not_simulator(2,.6, 1, percent_causal = 2.5, percent_overlap = "0")
# dat$kin %>% colnames
# dat$kin %>% rownames
# dat$file_paths
# kins <- snpStats::read.plink(dat$file_paths$X_Phi)
# kins
#
# x <- gaston::read.bed.matrix(dat$file_paths$X_Phi)
# slotNames(x)
# x@snps$chr %>% table
# x@sigma
#
# plot(x@p, x@sigma, xlim=c(0,1))
# t <- seq(0,1,length=101);
# lines(t, sqrt(2*t*(1-t)), col="red")
# plot(2*x@p, x@mu)
# abline(a=0,b=1, col = "red")
#
# as.matrix(x)[1:5,1:5]
# standardize(x) <- "p"
# x@standardize_mu_sigma
# as.matrix(x)[1:5,1:5]
# X <- as.matrix(x)
# # any(is.na(X[,2,drop=F]))
#
# kin <- gaston::GRM(x, autosome.only = FALSE)
# kin[1:5,1:5]
# all(complete.cases(kin))
# eiK <- eigen(kin)
#
# # deal with a small negative eigen value
# eiK$values[ eiK$values < 0 ] <- 0
# any(eiK$values < 0)
# PC <- sweep(eiK$vectors, 2, sqrt(eiK$values), "*")
# dim(PC)
# plot(PC[,1], PC[,2])
# PC[,1:10]
#
# rownames(kin)
# head(pheno_dat)
# hist(dat$y)
# dat$kin[1:5,1:5]
#
# # fit_lme <- coxme::lmekin(Y ~ 1 + (1|id), data = pheno_dat, varlist = xxmat)
#
# pacman::p_load(snpStats)
# xxmat <- snpStats::xxt(kins$genotypes, correct.for.missing = TRUE)
# evv <- eigen(xxmat, symmetric=TRUE)
# pcs <- evv$vectors[,1:5]
# evals <- evv$values[1:5]
# plot(evv$values)
# str(xxmat)
# xxmat[1:5,1:5]/4000
# kinPD <- as(nearPD(xxmat)$mat,"matrix")
# dimnames(kinPD)[[1]] <- dat$file_paths$Phi_names
# dimnames(kinPD)[[2]] <- kin_names$V1
# kin <- kinPD
#
# dat$kin[1:5,1:5]
# all(eigen(dat$kin)$values > 0)
# all(evv$values > 0)
# # pheno_dat <- data.frame(Y = dat$y, id = rownames(dat$kin))
#
# pheno_dat <- data.frame(Y = dat$y, id = rownames(xxmat))
# head(pheno_dat)
#
# # fitting with no intercept, i.e., 0 + (1|id) returns error
# fit_lme <- coxme::lmekin(Y ~ 1 + (1|id), data = pheno_dat, varlist = dat$kin)
# fit_lme <- coxme::lmekin(Y ~ 1 + (1|id), data = pheno_dat, varlist = xxmat)
# pheatmap::pheatmap(dat$kin)
# dat$kin[1:5,1:5]
# fit_lme %>% names
# newy <- residuals(fit_lme)
# plot( newy, dat$y)
# cor(newy,dat$y)
# abline(a=0,b=1)
# fit_lme %>% str
# coef(fit_lme)
# newy <- residuals(fit_lme)
# hist(newy)
# fitglmnet <- glmnet::cv.glmnet(x = dat$x_lasso, y = newy, standardize = F, alpha = 1, intercept = F,
#                                penalty.factor = c(rep(1, 4000), rep(0, 10)))
# plot(fitglmnet)
# yhat = predict(fitglmnet, newx = dat$x_lasso, s = "lambda.min")
# as.numeric(sqrt(crossprod(newy - yhat)))
# nonz <- setdiff(rownames(coef(fitglmnet, s = "lambda.min")[nonzeroCoef(coef(fitglmnet, s = "lambda.min")),,drop=F]),c("(Intercept)"))
# length(intersect(nonz, dat$causal))/length(dat$causal)
# coef(fitglmnet, s = "lambda.min")[1:3,]
#
# as.numeric(sqrt(crossprod(dat$y - (yhat+coef(fit_lme)$fixed+coef(fit_lme)$random$id))))
# as.numeric(sqrt(crossprod(dat$y - (yhat+coef(fit_lme)$fixed))))
#
# fitglmnet2 <- glmnet::cv.glmnet(x = dat$x_lasso, y = dat$y, standardize = T, alpha = 1, intercept = T,
#                                 penalty.factor = c(rep(1, 4000), rep(0, 10)))
# plot(fitglmnet2)
# yhat2 = predict(fitglmnet2, newx = dat$x_lasso, s = "lambda.min")
# as.numeric(sqrt(crossprod(dat$y - yhat2)))
# nonz2 <- setdiff(rownames(coef(fitglmnet2, s = "lambda.min")[nonzeroCoef(coef(fitglmnet2, s = "lambda.min")),,drop=F]),c("(Intercept)"))
# length(intersect(nonz2, dat$causal))/length(dat$causal)
#
#
# fitglmnetwithint <- glmnet::cv.glmnet(x = dat$x, y = newy, standardize = T, alpha = 1)
# plot(fitglmnetwithint)
#
# fit <- penfam(x = dat$x,
#               y = dat$y,
#               phi = dat$kin,
#               thresh_glmnet = 1e-10,
#               epsilon = 1e-5,
#               fdev = 1e-4,
#               alpha = 1,
#               tol.kkt = 1e-3,
#               nlambda = 100,
#               # an = log(log(model$n)) * log(model$n),
#               an = log(log(length(dat$y))),
#               # lambda_min_ratio  = ifelse(model$n < model$p, 0.01, 0.001),
#               lambda_min_ratio  = 0.05,
#               eta_init = 0.5,
#               maxit = 100)
#
# nonzero = predict(fit, type = "nonzero", s = fit$lambda_min)
# nonzero_names = setdiff(rownames(predict(fit, type = "nonzero", s = fit$lambda_min)), c("(Intercept)","eta","sigma2"))
# plot(fit)
#
# require(snpStats)
# data(for.exercise)
# controls <- rownames(subject.support)[subject.support$cc==0]
# use <- seq(1, ncol(snps.10), 10)
# ctl.10 <- snps.10[controls,use]
#
#
# ###################################################
# ### code chunk number 2: xxt-matrix
# ###################################################
# xxmat <- xxt(ctl.10, correct.for.missing=FALSE)
# xxmat[1:5,1:5]
#
# ###################################################
# ### code chunk number 3: eigen
# ###################################################
# evv <- eigen(xxmat, symmetric=TRUE)
# pcs <- evv$vectors[,1:5]
# evals <- evv$values[1:5]
# evals
#
#
#
# # dat$causal
# # dat$x %>% dim
# # dat$y %>% hist
# # dat$kin %>% dim
# # sum(dat$y>0)
# # plot(dat$y)
# # range(dat$y)
# # plot(density(dat$y))
# # # sum(draws@draws$r1.1>0)
# # # sum(draws@draws$r1.2>0)
# # # plot(dat$y, draws@draws$r1.1)
# # #
# # system.time(
# # fit <- penfam(x = dat$x,
# #               y = dat$y,
# #               phi = dat$kin,
# #               thresh_glmnet = 1e-10,
# #               epsilon = 1e-5,
# #               fdev = 1e-7,
# #               alpha = 1,
# #               tol.kkt = 1e-3,
# #               nlambda = 100,
# #               # an = log(log(1000)) * log(1000),
# #               an = log(log(1000)),
# #               # an = log(1000),
# #               # lambda_min_ratio  = ifelse(model$n < model$p, 0.01, 0.001),
# #               lambda_min_ratio  = 0.05,
# #               eta_init = 0.5,
# #               maxit = 100)
# # )
# # #
# # plot(fit, type = "BIC")
# # sum(rownames(coef(fit, s = fit$lambda_min)[nonzeroCoef(coef(fit, s = fit$lambda_min)),,drop=F]) %in% dat$causal)
# # rownames(coef(fit, s = fit$lambda_min)[nonzeroCoef(coef(fit, s = fit$lambda_min)),,drop=F]) %>% length()
# # sqrt(crossprod(fit$predicted[,fit$lambda_min] - dat$y))
# # #
# # #
# # #
# # #
# # fitglmnet <- cv.glmnet(x = dat$x, y = dat$y, alpha = 1, standardize = F)
# # plot(fitglmnet)
# # sum(rownames(coef(fitglmnet)[nonzeroCoef(coef(fitglmnet)),,drop=F] ) %in% dat$causal)
# # coef(fitglmnet)[nonzeroCoef(coef(fitglmnet)),,drop=F]
# # coef(fitglmnet)[nonzeroCoef(coef(fitglmnet)),,drop=F] %>% dim
# # sqrt(crossprod(predict(fitglmnet, newx = dat$x) - dat$y))
