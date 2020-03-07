# source("manuscript/bin/setup.R")
## ---- simulation-results ----

# df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/june_29_2018_results.rds")
# df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/july_1_2018_results.rds")
# df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/july_1_2018_results_with_twostepY.rds")
# df <- df %>% separate(Model,
#                       into = c("simnames","b0","eta_p","Fst","geography","k","n",
#                                "pkinship","ptest","percentcausal",
#                                "percentoverlap","s","sigma2_p"),
#                       sep = "/")

# df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/july_12_2018_results_with_null_model_VC.rds")
# df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/july_12_2018_results_with_null_model_VC_lasso_has_proper_MSE.rds")
# df <- readRDS("C:/Users/sahir/Documents/git_repositories/ggmix/simulation/simulation_results/july_12_2018_results_with_null_model_VC.rds")
# df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/may_02_2019_results.rds")
# df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/may_05_2019_results.rds") # this has lasso1se
# df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/may_06_2019_results.rds")
# df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/may_07_2019_results.rds")
# df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/jul_10_2019_results.rds") # this was used in first submission to plos genetics
df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/jul_10_2019_results_v2.rds") # this has TPR at FPR 5%, otherwise same jul_10_2019_results.rds

df <- df %>% separate(Model,
                      into = c("simnames","b0","beta_mean","eta_p","Fst","geography","k","n",
                               "pkinship","ptest","percentcausal",
                               "percentoverlap","s","sigma2_p"),
                      sep = "/")

DT <- data.table::as.data.table(df, stringsAsFactors = FALSE)

# DT[, table(percentoverlap, p_overlap)]
DT[, p_overlap := case_when(percentoverlap == "percent_overlap_0" ~ "No causal SNPs in Kinship",
                            percentoverlap == "percent_overlap_100" ~ "All causal SNPs in Kinship")]
DT[, p_causal := case_when(percentcausal == "percent_causal_0" ~ "Null model",
                            percentcausal == "percent_causal_0.01" ~ "1% of SNPs are causal")]

DT[, p_causal := factor(p_causal, levels = c("Null model","1% of SNPs are causal"))]

# DT[,table(geography)]
# DT[geography == "geography_ind", structure := "block"]
# DT[geography == "geography_circ", structure := "circular"]
## --PATCH- # to keep code with as little change as possible, im renaming 1D to block
## because everything below is for "block". Even though, the may7th results are for 1D structure.
DT[geography == "geography_1d", structure := "block"]
# DT[, table(geography, structure)]
# DT[, structure := factor(structure, levels = c("block","1D","circular"))]
DT[, structure := factor(structure, levels = c("block"))]
DT[, eta_p := case_when(eta_p == "eta_0.1" ~ "10% Heritability",
                        eta_p == "eta_0.3" ~ "30% Heritability")]
# DT[, table(eta_p)]
# use twostepY, which compares to the original Y
# DT[, table(Method)]
# DT <- DT[Method %ni% c("twostep","twostepY","lasso1se")]
# DT[Method == "twostepY", Method := "twostep"]
DT[Method == "twostepYVCCV", Method := "twostep"]
DT[Method == "lassoCV", Method := "lasso"]
DT[Method == "ggmixHDBIC", Method := "ggmix"]
# DT <- DT[Method != "lasso"]
# DT[Method == "lasso1se", Method := "lasso"]
# DT[,table(Method)]
# DT[, Method := factor(Method, levels = c("twostep","twostepY","lasso","ggmix"))]
DT[, Method := factor(Method, levels = c("twostep","lasso","ggmix"))]
# DT[, table(percentcausal,p_causal)]
# DT[, Method := factor(Method, levels = c("twostep","twostepY","lasso","ggmix"))]
# DT[, table(Method)]

DT[Method == "twostep", errorvar := sigma2]
DT[Method == "twostep", tau := eta]
DT[Method == "twostep", eta := tau/(tau + sigma2)]

DT[, me2 := (1/1000) * me^2]
# DT[Method == "twostep"]
appender <- function(string) latex2exp::TeX(paste(string))

theme_box <- theme_ipsum_rc(axis_title_just = "bt",axis_title_size = 16, axis = T, ticks = F, axis_col = "black") +
  theme(legend.position = "bottom",title = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        strip.text = element_text(size = 18))


## ---- table-of-results ----

DT$RMSE <- sqrt(DT$mse)
# tt <- DT %>% tidyr::pivot_longer(cols = c("me","prederror","tpr","fpr","nactive",
#                                     "eta","sigma2","tprFPR5","nactiveFPR5",
#                                     "correct_sparsity","mse","RMSE","errorvar","estimationerror",
#                                     "time","tau","sigma2","me2"),
#                                  names_to = "metric")

tt <- DT %>%
  group_by(Method, eta_p, p_overlap, p_causal) %>%
  summarise(RMSE = qwraps2::mean_sd(RMSE, denote_sd = "paren"),
            TPR = qwraps2::mean_sd(tprFPR5, denote_sd = "paren"),
            Heritability = qwraps2::mean_sd(eta, denote_sd = "paren"),
            Errorvar = qwraps2::mean_sd(errorvar, denote_sd = "paren"),
            Estimationerror = qwraps2::mean_sd(estimationerror, denote_sd = "paren"),
            Nactive = qwraps2::median_iqr(nactive, digits = 0))

tt[tt$Method=="lasso","Heritability"] <- "--"
# tt[tt$p_causal=="Null model","TPR"] <- "--"


pt <- tt %>% pivot_longer(cols = c("RMSE","TPR","Heritability","Errorvar","Estimationerror","Nactive"),
                    names_to = "metric") %>%
  unite(col = "type", p_causal,p_overlap,eta_p) %>%
  mutate(type = factor(type,
                       levels = c("Null model_No causal SNPs in Kinship_10% Heritability",
                                  "Null model_No causal SNPs in Kinship_30% Heritability",
                                  "Null model_All causal SNPs in Kinship_10% Heritability",
                                  "Null model_All causal SNPs in Kinship_30% Heritability",
                                  "1% of SNPs are causal_No causal SNPs in Kinship_10% Heritability",
                                  "1% of SNPs are causal_No causal SNPs in Kinship_30% Heritability",
                                  "1% of SNPs are causal_All causal SNPs in Kinship_10% Heritability",
                                  "1% of SNPs are causal_All causal SNPs in Kinship_30% Heritability")),
         metric = factor(metric,
                         labels = c("TPR at FPR=5%","Model Size","RMSE","Estimation Error","Error Variance","Heritability"),
                         levels = c("TPR","Nactive","RMSE","Estimationerror","Errorvar","Heritability"))
         ) %>%
  pivot_wider(id_cols = c("metric","Method"), names_from = "type", values_from = "value") %>%
  arrange(metric, Method) %>%
  select("metric","Method","Null model_No causal SNPs in Kinship_10% Heritability",
         "Null model_No causal SNPs in Kinship_30% Heritability",
         "Null model_All causal SNPs in Kinship_10% Heritability",
         "Null model_All causal SNPs in Kinship_30% Heritability",
         "1% of SNPs are causal_No causal SNPs in Kinship_10% Heritability",
         "1% of SNPs are causal_No causal SNPs in Kinship_30% Heritability",
         "1% of SNPs are causal_All causal SNPs in Kinship_10% Heritability",
         "1% of SNPs are causal_All causal SNPs in Kinship_30% Heritability")


## ---- print-sim-table ----

kable(pt, "latex", booktabs = T, align = c("l","l","c","c","c","c","c","c","c","c"),
      caption = c("Mean (standard deviation) from 200 simulations stratified by the number of causal SNPs (null, 1\\%), the overlap between causal SNPs and kinship matrix (no overlap, all causal SNPs in kinship), and true heritability (10\\%, 30\\%).
                  For all simulations, sample size is $n=1000$, the number of covariates is $p=5000$, and the number of SNPs used to estimate the kinship matrix is $k=10000$.
                  TPR at FPR=5\\% is the true positive rate at a fixed false positive rate of 5\\%.
                  Model Size ($|\\widehat{S}_{\\hat{\\lambda}}|$) is the number of selected variables in the training set using the high-dimensional BIC for \\texttt{ggmix} and 10-fold cross validation for \\texttt{lasso} and \\texttt{twostep}.
                  RMSE is the root mean squared error on the test set.
                  Estimation error is the squared distance between the estimated and true effect sizes.
                  Error variance ($\\sigma^2$) for \\texttt{twostep} is estimated from an intercept only LMM with a single random effect and is modeled explicitly in \\ggmix. For the \\texttt{lasso} we use $\\protect\\frac{1}{n - |\\widehat{S}_{\\hat{\\lambda}}|} \\protect||\\bY - \\bX \\widehat{\\bbeta}_{\\hat{\\lambda}}||_2^2$~\\citep{reid2016study} as an estimator for $\\sigma^2$.
                  Heritability ($\\eta$) for \\texttt{twostep} is estimated as $\\sigma_g^2 / (\\sigma_g^2 + \\sigma_e^2)$ from an intercept only LMM with a single random effect where $\\sigma_g^2$ and $\\sigma_e^2$ are the variance components for the random effect and error term, respectively. $\\eta$ is explictly modeled in \\ggmix. There is no positive way to calculate $\\eta$ for the \\texttt{lasso} since we are using a PC adjustment."),
      col.names = c("Metric","Method",rep(c("10%","30%"),4))) %>%
  kable_styling(latex_options = "striped", full_width = TRUE,
                position = "center",font_size = 7, stripe_index = c(1:3,7:9,13:15)) %>%
  add_header_above(c(" "," ","No overlap" = 2, "All causal SNPs\nin kinship" = 2, "No overlap" = 2, "All causal SNPs\nin kinship" = 2)) %>%
  add_header_above(c(" "," ","Null model" = 4, "1% Causal SNPs" = 4)) %>%
  column_spec(1, bold=T) %>%
  collapse_rows(columns = 1, latex_hline = "major", valign = "middle") %>%
  footnote(general = "Median (Inter-quartile range) is given for Model Size.",
           general_title = "Note:") #%>%
  # footnote(general = "Model Size is the number of selected variables using the high-dimensional BIC for \texttt{ggmix} and 10-fold cross validation for \texttt{lasso} and \texttt{twostep}.")
  #          # number = c("", "Footnote 2; "),
  #          # alphabet = c("Footnote A; ", "Footnote B; "),
  #          # symbol = c("Footnote Symbol 1; ", "Footnote Symbol 2"),
  #          general_title = "Note: ", number_title = "Type I: ",
  #          alphabet_title = "Type II: ", symbol_title = "Type III: ",
  #          footnote_as_chunk = T, title_format = c("italic", "underline")
  # )


## ---- table-of-results-for-n-equal-to-k ----

df2 <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/dec_5_2019_results.rds") # this has n=1000=k=1000

df2 <- df2 %>% separate(Model,
                      into = c("simnames","b0","beta_mean","eta_p","Fst","geography","k","n",
                               "pkinship","ptest","percentcausal",
                               "percentoverlap","s","sigma2_p"),
                      sep = "/")

DT2 <- data.table::as.data.table(df2, stringsAsFactors = FALSE)

# DT[, table(percentoverlap, p_overlap)]
DT2[, p_overlap := case_when(percentoverlap == "percent_overlap_0" ~ "No causal SNPs in Kinship",
                            percentoverlap == "percent_overlap_100" ~ "All causal SNPs in Kinship")]
DT2[, p_causal := case_when(percentcausal == "percent_causal_0" ~ "Null model",
                           percentcausal == "percent_causal_0.01" ~ "1% of SNPs are causal")]

DT2[, p_causal := factor(p_causal, levels = c("Null model","1% of SNPs are causal"))]

# DT2[,table(geography)]
# DT[geography == "geography_ind", structure := "block"]
# DT[geography == "geography_circ", structure := "circular"]
## --PATCH- # to keep code with as little change as possible, im renaming 1D to block
## because everything below is for "block". Even though, the may7th results are for 1D structure.
DT2[geography == "geography_1d", structure := "block"]
# DT[, table(geography, structure)]
# DT[, structure := factor(structure, levels = c("block","1D","circular"))]
DT2[, structure := factor(structure, levels = c("block"))]
DT2[, eta_p := case_when(eta_p == "eta_0.1" ~ "10% Heritability",
                        eta_p == "eta_0.3" ~ "30% Heritability")]
# DT[, table(eta_p)]
# use twostepY, which compares to the original Y
# DT[, table(Method)]
# DT <- DT[Method %ni% c("twostep","twostepY","lasso1se")]
# DT[Method == "twostepY", Method := "twostep"]
DT2[Method == "twostepYVCCV", Method := "twostep"]
DT2[Method == "lassoCV", Method := "lasso"]
DT2[Method == "ggmixHDBIC", Method := "ggmix"]
# DT <- DT[Method != "lasso"]
# DT[Method == "lasso1se", Method := "lasso"]
# DT[,table(Method)]
# DT[, Method := factor(Method, levels = c("twostep","twostepY","lasso","ggmix"))]
DT2[, Method := factor(Method, levels = c("twostep","lasso","ggmix"))]
# DT[, table(percentcausal,p_causal)]
# DT[, Method := factor(Method, levels = c("twostep","twostepY","lasso","ggmix"))]
# DT[, table(Method)]

DT2[Method == "twostep", errorvar := sigma2]
DT2[Method == "twostep", tau := eta]
DT2[Method == "twostep", eta := tau/(tau + sigma2)]

DT2[, me2 := (1/1000) * me^2]


DT2$RMSE <- sqrt(DT2$mse)
# tt <- DT %>% tidyr::pivot_longer(cols = c("me","prederror","tpr","fpr","nactive",
#                                     "eta","sigma2","tprFPR5","nactiveFPR5",
#                                     "correct_sparsity","mse","RMSE","errorvar","estimationerror",
#                                     "time","tau","sigma2","me2"),
#                                  names_to = "metric")

tt2 <- DT2 %>%
  group_by(Method, eta_p, p_overlap, p_causal) %>%
  summarise(RMSE = qwraps2::mean_sd(RMSE, denote_sd = "paren"),
            TPR = qwraps2::mean_sd(tprFPR5, denote_sd = "paren"),
            Heritability = qwraps2::mean_sd(eta, denote_sd = "paren"),
            Errorvar = qwraps2::mean_sd(errorvar, denote_sd = "paren"),
            Estimationerror = qwraps2::mean_sd(estimationerror, denote_sd = "paren"),
            Nactive = qwraps2::median_iqr(nactive, digits = 0))

tt2[tt2$Method=="lasso","Heritability"] <- "--"



pt2 <- tt2 %>% pivot_longer(cols = c("RMSE","TPR","Heritability","Errorvar","Estimationerror","Nactive"),
                          names_to = "metric") %>%
  unite(col = "type", p_causal,p_overlap,eta_p) %>%
  mutate(type = factor(type,
                       levels = c("Null model_No causal SNPs in Kinship_10% Heritability",
                                  "Null model_No causal SNPs in Kinship_30% Heritability",
                                  "Null model_All causal SNPs in Kinship_10% Heritability",
                                  "Null model_All causal SNPs in Kinship_30% Heritability",
                                  "1% of SNPs are causal_No causal SNPs in Kinship_10% Heritability",
                                  "1% of SNPs are causal_No causal SNPs in Kinship_30% Heritability",
                                  "1% of SNPs are causal_All causal SNPs in Kinship_10% Heritability",
                                  "1% of SNPs are causal_All causal SNPs in Kinship_30% Heritability")),
         metric = factor(metric,
                         labels = c("TPR at FPR=5%","Model Size","RMSE","Estimation Error","Error Variance","Heritability"),
                         levels = c("TPR","Nactive","RMSE","Estimationerror","Errorvar","Heritability"))
  ) %>%
  pivot_wider(id_cols = c("metric","Method"), names_from = "type", values_from = "value") %>%
  arrange(metric, Method) %>%
  select("metric","Method","Null model_No causal SNPs in Kinship_10% Heritability",
         "Null model_No causal SNPs in Kinship_30% Heritability",
         "Null model_All causal SNPs in Kinship_10% Heritability",
         "Null model_All causal SNPs in Kinship_30% Heritability",
         "1% of SNPs are causal_No causal SNPs in Kinship_10% Heritability",
         "1% of SNPs are causal_No causal SNPs in Kinship_30% Heritability",
         "1% of SNPs are causal_All causal SNPs in Kinship_10% Heritability",
         "1% of SNPs are causal_All causal SNPs in Kinship_30% Heritability")


## ---- print-sim-table-for-n-equal-to-k ----

kable(pt2, "latex", booktabs = T, align = c("l","l","c","c","c","c","c","c","c","c"),
      caption = c("Mean (standard deviation) from 200 simulations stratified by the number of causal SNPs (null, 1\\%), the overlap between causal SNPs and kinship matrix (no overlap, all causal SNPs in kinship), and true heritability (10\\%, 30\\%).
                  For all simulations, sample size is $n=1000$, the number of covariates is $p=5000$, and the number of SNPs used to estimate the kinship matrix is $k=1000$.
                  TPR at FPR=5\\% is the true positive rate at a fixed false positive rate of 5\\%.
                  Model Size ($|\\widehat{S}_{\\hat{\\lambda}}|$) is the number of selected variables in the training set using the high-dimensional BIC for \\texttt{ggmix} and 10-fold cross validation for \\texttt{lasso} and \\texttt{twostep}.
                  RMSE is the root mean squared error on the test set.
                  Estimation error is the squared distance between the estimated and true effect sizes.
                  Error variance ($\\sigma^2$) for \\texttt{twostep} is estimated from an intercept only LMM with a single random effect and is modeled explicitly in \\ggmix. For the \\texttt{lasso} we use $\\protect\\frac{1}{n - |\\widehat{S}_{\\hat{\\lambda}}|} \\protect||\\bY - \\bX \\widehat{\\bbeta}_{\\hat{\\lambda}}||_2^2$~\\citep{reid2016study} as an estimator for $\\sigma^2$.
                  Heritability ($\\eta$) for \\texttt{twostep} is estimated as $\\sigma_g^2 / (\\sigma_g^2 + \\sigma_e^2)$ from an intercept only LMM with a single random effect where $\\sigma_g^2$ and $\\sigma_e^2$ are the variance components for the random effect and error term, respectively. $\\eta$ is explictly modeled in \\ggmix. There is no positive way to calculate $\\eta$ for the \\texttt{lasso} since we are using a PC adjustment."),
      col.names = c("Metric","Method",rep(c("10%","30%"),4))) %>%
  kable_styling(latex_options = "striped", full_width = TRUE,
                position = "center",font_size = 7, stripe_index = c(1:3,7:9,13:15)) %>%
  add_header_above(c(" "," ","No overlap" = 2, "All causal SNPs\nin kinship" = 2, "No overlap" = 2, "All causal SNPs\nin kinship" = 2)) %>%
  add_header_above(c(" "," ","Null model" = 4, "1% Causal SNPs" = 4)) %>%
  column_spec(1, bold=T) %>%
  collapse_rows(columns = 1, latex_hline = "major", valign = "middle") %>%
  footnote(general = "Median (Inter-quartile range) is given for Model Size.",
           general_title = "Note:") #%>%
# general_title = "Note:") %>%
# footnote(general = "Model Size is the number of selected variables using the high-dimensional BIC for \texttt{ggmix} and 10-fold cross validation for \texttt{lasso} and \texttt{twostep}.")
#          # number = c("", "Footnote 2; "),
#          # alphabet = c("Footnote A; ", "Footnote B; "),
#          # symbol = c("Footnote Symbol 1; ", "Footnote Symbol 2"),
#          general_title = "Note: ", number_title = "Type I: ",
#          alphabet_title = "Type II: ", symbol_title = "Type III: ",
#          footnote_as_chunk = T, title_format = c("italic", "underline")
# )




## ---- plot-kinship-sim ----

dat <- lapply(#list("ind","1d","circ"),
              list("1d"),
              function(i)
                ggmix::gen_structured_model(n = 1200,
                                            p_design = 5000,
                                            p_kinship = 10000,
                                            geography = i,
                                            percent_causal = 0.01,
                                            percent_overlap = "100",
                                            k = 10, s = 0.5, Fst = 0.1,
                                            b0 = 0,
                                            eta = 0.1, sigma2 = 1,
                                            train_tune_test = c(0.99,.005,0.005))
)
# str(dat)
# dev.off()

# par(omi = c(0.3,0.3,0.3,0.3))
# popkin::plotPopkin(x = list(dat[[1]]$kin, dat[[2]]$kin, dat[[3]]$kin),
#                    titles = c("block", "1D", "circular"),
#                    marPad = 0.05)
# popkin::plot_popkin(kinship = list(dat[[1]]$kin_train),
#                     titles = c("Empirical Kinship Matrix with 1D Structure"),
#                     marPad = 0.05)
popkin::plot_popkin(kinship = list(dat[[1]]$kinship)#,
                    # col = RColorBrewer::display.brewer.all(),
                    # titles = c("Empirical Kinship Matrix with 1D Structure"),
                    # mar_pad = 0.05
                    )
# popkin::plot_popkin(kinship = list(dat[[1]]$coancestry),
#                     # titles = c("Empirical Kinship Matrix with 1D Structure"),
#                     marPad = 0.05)

## ---- plot-pc-sim ----

xlabs <- "1st principal component"
ylabs <- "2nd principal component"


# par(mfrow = c(1,3))
plot(dat[[1]]$PC[,1],dat[[1]]$PC[,2],
     pch = 19, col = rep(RColorBrewer::brewer.pal(10,"Paired"), each = table(dat[[1]]$subpops)[1]),
     xlab = xlabs, ylab = ylabs, cex = 0.5)#,
     # main = "1D Structure")
# plot(dat[[2]]$x_lasso[,5001],dat[[2]]$x_lasso[,5002],
#      pch = 19, col = rep(RColorBrewer::brewer.pal(5,"Paired"), each = 200),
#      xlab = xlabs, ylab = ylabs,
#      main = "1D Geography")
# plot(dat[[3]]$x_lasso[,5001],dat[[3]]$x_lasso[,5002],
#      pch = 19, col = rep(RColorBrewer::brewer.pal(5,"Paired"), each = 200),
#      xlab = xlabs, ylab = ylabs,
#      main = "Circular Geography")

## ---- plot-kinship-pc-combined-sim ----

# not used in manuscript
dev.off()
m <- rbind(c(1, 2,3), c(4,5,6))
layout(m)
layout.show(m)
layout(m)
par(mar = c(3, 3, 0, 0))
for (i in 1:3) plot(1, 1, type = "n")

def.par <- par(no.readonly = TRUE) # save default, for resetting...

## divide the device into two rows and two columns
## allocate figure 1 all of row 1
## allocate figure 2 the intersection of column 2 and row 2
dev.off()
par()$mar
layout(matrix(c(1,2,3,4,5,6,7,7), 2, 4, byrow = TRUE))
par(omi = c(0.3,0.3,0.3,0.3), pty="s")
popkin::plotPopkin(x = list(dat[[1]]$kin, dat[[2]]$kin, dat[[3]]$kin),
                   titles = c("block", "1D", "circular"),
                   legMar = c(5,0.2,4,2),
                   marPad = 0.05, addLayout = FALSE)
## show the regions that have been allocated to each plot
layout.show(7)



## divide device into two rows and two columns
## allocate figure 1 and figure 2 as above
## respect relations between widths and heights
nf <- layout(matrix(c(1,1,0,2), 2, 2, byrow = TRUE), respect = TRUE)
layout.show(nf)

## ---- plot-correct-sparsity-sim-null-model ----

p1_cs <- ggplot(DT[p_causal == "Null model"][p_overlap == "No causal SNPs in Kinship"][structure == "block"],
                aes(Method, correct_sparsity, fill = Method)) +
  ggplot2::geom_boxplot() +
  facet_rep_grid(p_overlap ~ eta_p, scales = "fixed",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,3,4,2)]) +
  labs(x = "", y = "",
       title = "Correct Sparsity Results for the Null Model",
       subtitle = "Based on 200 simulations",
       caption = "") +
  theme_box

p1_cs


## ---- plot-estimation-error-sim-null-model ----

p1_esterror <- ggplot(DT[p_causal == "Null model"][p_overlap == "No causal SNPs in Kinship"][structure == "block"],
                aes(Method, estimationerror, fill = Method)) +
  ggplot2::geom_boxplot() +
  facet_rep_grid(p_overlap ~ eta_p, scales = "fixed",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,3,4,2)]) +
  labs(x = "", y = "",
       title = "Estimation Error Results for the Null Model",
       subtitle = "Based on 200 simulations",
       caption = "") +
  theme_box

p1_esterror

# reposition_legend(p1_mse, 'center', panel='panel-2-3')

## ---- plot-6in1-10percentHerit-1pcausal-allinkinship ----

pm_cs <- ggplot(DT[p_causal != "Null model"][structure == "block"][eta_p == "10% Heritability"][p_overlap == "All causal SNPs in Kinship"],
                aes(Method, correct_sparsity, fill = Method)) +
  ggplot2::geom_boxplot() +
  # facet_rep_grid(p_overlap ~ eta_p, scales = "fixed",
  #                repeat.tick.labels = 'left',
  #                labeller = as_labeller(appender,
  #                                       default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,3,4,2)]) +
  labs(x = "", y = "Correct sparsity",
       # title = "Correct Sparsity results for the Model with 1% Causal SNPs",
       subtitle = "(A)"
       # caption = ""
       ) +
  theme_box + theme(legend.position = "none")

pm_esterror <- ggplot(DT[p_causal != "Null model"][structure == "block"][eta_p == "10% Heritability"][p_overlap == "All causal SNPs in Kinship"],
                aes(Method, estimationerror, fill = Method)) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  # facet_rep_grid(p_overlap ~ eta_p, scales = "fixed",
  #                repeat.tick.labels = 'left',
  #                labeller = as_labeller(appender,
  #                                       default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,3,4,2)]) +
  labs(x = "",
       y = latex2exp::TeX("Estimation error"),
       # y = latex2exp::TeX("Estimation error $(||\\hat{\\beta} - \\beta_{truth}||_2^2)$"),
       # title = "Correct Sparsity results for the Model with 1% Causal SNPs",
       subtitle = "(B)"
       # caption = ""
  ) +
  theme_box + theme(legend.position = "none")+ scale_y_continuous(limits = quantile(DT$estimationerror, c(0.1, 0.9)))

# PATCH DOING RMSE instead of MSE, because that how RDA was done by Tianyuan
df_me_nactive <- DT[structure == "block"][p_causal != "Null model", c("Method", "eta_p", "p_overlap", "nactive", "mse")][eta_p == "10% Heritability"][p_overlap == "All causal SNPs in Kinship"] %>%
  group_by(Method, eta_p, p_overlap) %>%
  summarise(mean.me = mean(sqrt(mse), na.rm = TRUE), sd.me = sd(sqrt(mse), na.rm = TRUE),
            mean.nactive = mean(nactive, na.rm = TRUE), sd.nactive = sd(nactive, na.rm = TRUE))

pm_mse_nactive <- ggplot(data = df_me_nactive,
                        aes(x = mean.nactive, y = mean.me, color = Method, label = Method)) +
  geom_point(size = 2.1) +
  geom_text_repel(
    data = subset(df_me_nactive, Method == "ggmix"),
    nudge_x      = 15,
    nudge_y = 0.15,
    # size = 8,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
    geom_text_repel(
      data = subset(df_me_nactive, Method == "lasso"),
      nudge_x      = -90,
      nudge_y = -0.10,
      # size = 8,
      direction    = "y",
      hjust        = 0,
      segment.size = 0.2
    ) +
  geom_text_repel(
    data = subset(df_me_nactive, mean.nactive >= 287),
    nudge_x      = 35,
    nudge_y = 0.10,
    # size = 8,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_errorbar(aes(ymin = mean.me - sd.me, ymax = mean.me + sd.me, width=5), size = 1.1) +
  geom_errorbarh(aes(xmin = mean.nactive - sd.nactive, xmax = mean.nactive + sd.nactive, height = 0.02), size = 1.1) +
  # facet_rep_grid(p_overlap ~ eta_p, scales = "free",
                 # repeat.tick.labels = 'left',
                 # labeller = as_labeller(appender,
                                        # default = label_parsed)) +
  scale_color_manual(values = cbbPalette[c(7,3,4,2)], guide = guide_legend(ncol=3)) +
  labs(x = "Number of active variables", y = "Root mean squared prediction error",
       # title = "Mean Squared Error vs. Number of Active Variable (Mean +/- 1 SD) for Model with 1% Causal SNPs",
       subtitle = "(C)",
       caption = "mean +/- 1 standard deviation"
       ) +
  theme_box+ theme(legend.position = "none") +
  scale_x_continuous(limits = c(0,400), breaks = seq(0,400,100))



# DT[, table(Method)]
df_me_nactive <- DT[structure == "block"][Method %in% c("lasso","ggmix")][p_causal != "Null model", c("Method", "eta_p", "p_overlap", "nactive", "mse")][eta_p == "10% Heritability"][p_overlap == "All causal SNPs in Kinship"] %>%
  group_by(Method, eta_p, p_overlap) %>%
  summarise(mean.me = mean(sqrt(mse), na.rm = TRUE), sd.me = sd(sqrt(mse), na.rm = TRUE),
            mean.nactive = mean(nactive, na.rm = TRUE), sd.nactive = sd(nactive, na.rm = TRUE))

pm_mse_nactive_zoom <- ggplot(data = df_me_nactive,
                         aes(x = mean.nactive, y = mean.me, color = Method, label = Method)) +
  geom_point(size = 2.1) +
  geom_text_repel(
    data = subset(df_me_nactive, mean.nactive < 150),
    nudge_x      = 11,
    nudge_y = 0.04,
    # size = 8,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_text_repel(
    data = subset(df_me_nactive, mean.nactive >= 150),
    nudge_x      = 11,
    nudge_y = .06,
    # size = 8,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_errorbar(aes(ymin = mean.me - sd.me, ymax = mean.me + sd.me, width = 5), size = 1.1) +
  geom_errorbarh(aes(xmin = mean.nactive - sd.nactive, xmax = mean.nactive + sd.nactive, height = 0.007), size = 1.1) +
  # facet_rep_grid(p_overlap ~ eta_p, scales = "free",
  # repeat.tick.labels = 'left',
  # labeller = as_labeller(appender,
  # default = label_parsed)) +
  scale_color_manual(values = cbbPalette[c(3,4,2)], guide = guide_legend(ncol=3)) +
  labs(x = "Number of active variables", y = "Root mean squared prediction error",
       # title = "Mean Squared Error vs. Number of Active Variable (Mean +/- 1 SD) for Model with 1% Causal SNPs",
       subtitle = "(C)",
       caption = "mean +/- 1 standard deviation"
  ) +
  theme_box+ theme(legend.position = "none") +
  scale_x_continuous(limits = c(0,350))


# dummy2 <- data.frame(eta_p = c("10% Heritability", "30% Heritability"), Z = c(0.1, 0.5))
dummy2 <- data.frame(eta_p = c("10% Heritability"), Z = c(0.1))

pm_eta <- ggplot(DT[structure == "block"][p_causal != "Null model"][Method %in% c("twostep","ggmix")][eta_p == "10% Heritability"][p_overlap == "All causal SNPs in Kinship"],
                 aes(Method, eta, fill = Method)) +
  ggplot2::geom_boxplot() +
  # facet_rep_grid(p_overlap ~ eta_p, scales = "free",
                 # repeat.tick.labels = 'left',
                 # labeller = as_labeller(appender,
                                        # default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,4)]) +
  labs(x = "", y = TeX("Estimated heritability $(\\hat{\\eta})$"),
       # title = TeX("Estimated Heritability for the Model with 1% Causal SNPs"),
       subtitle = "(E)",
       caption = "horizontal dashed line is the true value") +
  theme_ipsum_rc(axis_title_just = "bt",axis_title_size = 16, axis = T, ticks = F, axis_col = "black") +
  theme(legend.position = "none",title = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 18)) +
  geom_hline(data = dummy2, aes(yintercept = Z), linetype = 2, col = "#2f4f4f")

# dummy2 <- data.frame(eta_p = c("10% Heritability", "30% Heritability"), Z = c((1 - 0.1), (1 - 0.3)))
dummy2 <- data.frame(eta_p = c("10% Heritability"), Z = c((1 - 0.1)))

pm_errorvar <- ggplot(DT[structure == "block"][p_causal != "Null model"][eta_p == "10% Heritability"][p_overlap == "All causal SNPs in Kinship"],
                      aes(Method, errorvar, fill = Method)) +
  ggplot2::geom_boxplot() +
  # facet_rep_grid(p_overlap ~ eta_p, scales = "free",
  #                repeat.tick.labels = 'left',
  #                labeller = as_labeller(appender,
  #                                       default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,3,4)]) +
  labs(x = "", y = TeX("Estimated error variance $(\\hat{\\sigma^2})$"),
       # title = TeX("Estimated Error Variance for the Model with 1% Causal SNPs"),
       subtitle = "(F)",
       caption = "horizontal dashed line is the true value") +
  theme_ipsum_rc(axis_title_just = "bt",axis_title_size = 16, axis = T, ticks = F, axis_col = "black") +
  theme(legend.position = "none",title = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        strip.text = element_text(size = 18)) +
  geom_hline(data = dummy2, aes(yintercept = Z), linetype = 2, col = "#2f4f4f")


df_tpr_fpr <- DT[structure == "block"][p_causal != "Null model"][, c("Method", "eta_p", "p_overlap", "tpr", "fpr")][eta_p == "10% Heritability"][p_overlap == "All causal SNPs in Kinship"] %>%
  group_by(Method, eta_p, p_overlap) %>%
  summarise(mean.tpr = mean(tpr, na.rm = TRUE), sd.tpr = sd(tpr, na.rm = TRUE),
            mean.fpr = mean(fpr, na.rm = TRUE), sd.fpr = sd(fpr, na.rm = TRUE))

pm_tpr_fpr <- ggplot(data = df_tpr_fpr, aes(x = mean.fpr, y = mean.tpr, color = Method, label = Method)) +
  geom_point(size = 2.1) +
  # geom_text_repel() +
  geom_text_repel(
    data = subset(df_tpr_fpr, mean.fpr < 0.02),
    nudge_x      = 0.002,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_text_repel(
    data = subset(df_tpr_fpr, Method == "lasso"),
    nudge_x      = -0.012,
    nudge_y = 0.049,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_text_repel(
    data = subset(df_tpr_fpr, Method == "twostep"),
    nudge_x      = 0.007,
    nudge_y = 0.03,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_errorbar(aes(ymin = mean.tpr - sd.tpr, ymax = mean.tpr + sd.tpr, width=0.0005), size = 1.1) +
  geom_errorbarh(aes(xmin = mean.fpr - sd.fpr, xmax = mean.fpr + sd.fpr, height = 0.009), size = 1.1) +
  # facet_rep_grid(p_overlap ~ eta_p, scales = "free",
  #                repeat.tick.labels = 'left',
  #                labeller = as_labeller(appender,
  #                                       default = label_parsed)) +
  scale_color_manual(values = cbbPalette[c(7,3,4)], guide=guide_legend(ncol=3)) +
  labs(x="False positive rate", y="True positive rate",
       # title="True Positive Rate vs. False Positive Rate (Mean +/- 1 SD) for the Model with 1% Causal SNPs",
       subtitle="(D)",
       caption="mean +/- 1 standard deviation") +
  theme_box + scale_y_continuous(limits = c(0.6,1), breaks = seq(0.6,1, 0.1)) +
  scale_x_continuous(limits = c(0,.06), breaks = seq(0,0.06, 0.02)) + theme(legend.position = "none")


pm_cs +
  pm_esterror+
  pm_mse_nactive+
  # pm_mse_nactive_zoom +
  pm_tpr_fpr +
  pm_eta +
  pm_errorvar
# dev.off()
#
#
# cowplot::plot_grid(pm_mse_nactive,
#                    pm_mse_nactive_zoom,
#                    pm_cs,
#                    pm_eta,
#                    pm_errorvar,
#                    pm_tpr_fpr, nrow = 3)
# library(cowplot)
# plot_to_gtable(pm_eta)
# plot_to_gtable(pm_errorvar)
# plot_to_gtable(pm_cs)
# plot_to_gtable(pm_mse_nactive)
# plot_to_gtable(pm_mse_nactive_zoom)
# plot_to_gtable(pm_tpr_fpr)

## ---- plot-correct-sparsity-sim-1p-causal ----

p1_cs <- ggplot(DT[p_causal != "Null model"][structure == "block"],
                aes(Method, correct_sparsity, fill = Method)) +
  ggplot2::geom_boxplot() +
  facet_rep_grid(p_overlap ~ eta_p, scales = "fixed",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,3,4,2)]) +
  labs(x = "", y = "",
       title = "Correct Sparsity results for the Model with 1% Causal SNPs",
       subtitle = "Based on 200 simulations",
       caption = "") +
  theme_box


p1_cs


## ---- plot-estimation-error-sim-1p-causal ----

p1_esterror_1p <- ggplot(DT[p_causal != "Null model"][structure == "block"],
                aes(Method, estimationerror, fill = Method)) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  facet_rep_grid(p_overlap ~ eta_p, scales = "fixed",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,3,4,2)]) +
  labs(x = "", y = "",
       title = "Estimation Error results for the Model with 1% Causal SNPs",
       subtitle = "Based on 200 simulations",
       caption = "") +
  theme_box+ scale_y_continuous(limits = quantile(DT$estimationerror, c(0.1, 0.9)))


p1_esterror_1p


## ---- plot-me-nactive-sim-null ----

df_me_nactive <- DT[p_overlap == "No causal SNPs in Kinship"][structure == "block"][p_causal == "Null model", c("Method", "eta_p", "p_overlap", "nactive", "mse")] %>%
  group_by(Method, eta_p, p_overlap) %>%
  summarise(mean.me = mean(mse, na.rm = TRUE), sd.me = sd(mse, na.rm = TRUE),
         mean.nactive = mean(nactive, na.rm = TRUE), sd.nactive = sd(nactive, na.rm = TRUE))

p1_me_nactive <- ggplot(data = df_me_nactive,
                        aes(x = mean.nactive, y = mean.me, color = Method, label = Method)) +
  geom_point(size = 2.1) +
  geom_text_repel(
    data = subset(df_me_nactive, mean.nactive < 100),
    nudge_x      = 5,
    nudge_y = 9,
    # size = 8,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_text_repel(
    data = subset(df_me_nactive, mean.nactive >= 100),
    nudge_x      = 5,
    nudge_y = 9,
    # size = 8,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_errorbar(aes(ymin = mean.me - sd.me, ymax = mean.me + sd.me), size = 1.1) +
  geom_errorbarh(aes(xmin = mean.nactive - sd.nactive, xmax = mean.nactive + sd.nactive), size = 1.1) +
  facet_rep_grid(p_overlap ~ eta_p, scales = "free",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_color_manual(values = cbbPalette[c(7,3,4,2)], guide = guide_legend(ncol=3)) +
  labs(x = "Number of active variables", y = "Model Error",
       title = "Model Error vs. Number of Active Variable (Mean +/- 1 SD) for Null Model",
       subtitle = "Based on 200 simulations",
       caption = "") +
  theme_box

p1_me_nactive


## ---- plot-me-nactive-sim-1p-causal ----

df_me_nactive <- DT[structure == "block"][p_causal != "Null model", c("Method", "eta_p", "p_overlap", "nactive", "me")] %>%
  group_by(Method, eta_p, p_overlap) %>%
  summarise(mean.me = mean(me, na.rm = TRUE), sd.me = sd(me, na.rm = TRUE),
            mean.nactive = mean(nactive, na.rm = TRUE), sd.nactive = sd(nactive, na.rm = TRUE))

p1_me_nactive <- ggplot(data = df_me_nactive,
                        aes(x = mean.nactive, y = mean.me, color = Method, label = Method)) +
  geom_point(size = 2.1) +
  geom_text_repel(
    data = subset(df_me_nactive, mean.nactive < 100),
    nudge_x      = 5,
    nudge_y = 9,
    # size = 8,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_text_repel(
    data = subset(df_me_nactive, mean.nactive >= 100),
    nudge_x      = 5,
    nudge_y = 9,
    # size = 8,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_errorbar(aes(ymin = mean.me - sd.me, ymax = mean.me + sd.me), size = 1.1) +
  geom_errorbarh(aes(xmin = mean.nactive - sd.nactive, xmax = mean.nactive + sd.nactive), size = 1.1) +
  facet_rep_grid(p_overlap ~ eta_p, scales = "free",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_color_manual(values = cbbPalette[c(7,3,4,2)], guide = guide_legend(ncol=3)) +
  labs(x = "Number of active variables", y = "Model Error",
       title = "Model Error vs. Number of Active Variable (Mean +/- 1 SD) for Model with 1% Causal SNPs",
       subtitle = "Based on 200 simulations",
       caption = "") +
  theme_box

p1_me_nactive

## ---- plot-mse-nactive-sim-null ----

df_mse_nactive <- DT[p_overlap == "No causal SNPs in Kinship"][structure == "block"][p_causal == "Null model", c("Method", "eta_p", "p_overlap", "nactive", "mse")] %>%
  group_by(Method, eta_p, p_overlap) %>%
  summarise(mean.me = mean(mse, na.rm = TRUE), sd.me = sd(mse, na.rm = TRUE),
            mean.nactive = mean(nactive, na.rm = TRUE), sd.nactive = sd(nactive, na.rm = TRUE))

p1_mse_nactive <- ggplot(data = df_mse_nactive,
                        aes(x = mean.nactive, y = mean.me, color = Method, label = Method)) +
  geom_point(size = 2.1) +
  geom_text_repel(
    data = subset(df_mse_nactive, mean.nactive < 100),
    nudge_x      = 5,
    nudge_y = 9,
    # size = 8,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_text_repel(
    data = subset(df_mse_nactive, mean.nactive >= 100),
    nudge_x      = 5,
    nudge_y = 9,
    # size = 8,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_errorbar(aes(ymin = mean.me - sd.me, ymax = mean.me + sd.me), size = 1.1) +
  geom_errorbarh(aes(xmin = mean.nactive - sd.nactive, xmax = mean.nactive + sd.nactive), size = 1.1) +
  facet_rep_grid(p_overlap ~ eta_p, scales = "free",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_color_manual(values = cbbPalette[c(7,3,4,2)], guide = guide_legend(ncol=3)) +
  labs(x = "Number of active variables", y = "Mean Squared Error",
       title = "Mean Squared Error vs. Number of Active Variable (Mean +/- 1 SD) for Null Model",
       subtitle = "Based on 200 simulations",
       caption = "") +
  theme_box

p1_mse_nactive


## ---- plot-mse-nactive-sim-1p-causal ----

df_me_nactive <- DT[structure == "block"][p_causal != "Null model", c("Method", "eta_p", "p_overlap", "nactive", "mse")] %>%
  group_by(Method, eta_p, p_overlap) %>%
  summarise(mean.me = mean(mse, na.rm = TRUE), sd.me = sd(mse, na.rm = TRUE),
            mean.nactive = mean(nactive, na.rm = TRUE), sd.nactive = sd(nactive, na.rm = TRUE))

p1_me_nactive <- ggplot(data = df_me_nactive,
                        aes(x = mean.nactive, y = mean.me, color = Method, label = Method)) +
  geom_point(size = 2.1) +
  geom_text_repel(
    data = subset(df_me_nactive, mean.nactive < 100),
    nudge_x      = 5,
    nudge_y = 9,
    # size = 8,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_text_repel(
    data = subset(df_me_nactive, mean.nactive >= 100),
    nudge_x      = 5,
    nudge_y = 9,
    # size = 8,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_errorbar(aes(ymin = mean.me - sd.me, ymax = mean.me + sd.me), size = 1.1) +
  geom_errorbarh(aes(xmin = mean.nactive - sd.nactive, xmax = mean.nactive + sd.nactive), size = 1.1) +
  facet_rep_grid(p_overlap ~ eta_p, scales = "free",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_color_manual(values = cbbPalette[c(7,3,4,2)], guide = guide_legend(ncol=3)) +
  labs(x = "Number of active variables", y = "Mean Squared Error",
       title = "Mean Squared Error vs. Number of Active Variable (Mean +/- 1 SD) for Model with 1% Causal SNPs",
       subtitle = "Based on 200 simulations",
       caption = "") +
  theme_box

p1_me_nactive


## ---- plot-mse-nactive-sim-1p-causal-zoom-in ----

df_me_nactive <- DT[structure == "block"][Method %in% c("lasso","ggmix")][p_causal != "Null model", c("Method", "eta_p", "p_overlap", "nactive", "mse")] %>%
  group_by(Method, eta_p, p_overlap) %>%
  summarise(mean.me = mean(mse, na.rm = TRUE), sd.me = sd(mse, na.rm = TRUE),
            mean.nactive = mean(nactive, na.rm = TRUE), sd.nactive = sd(nactive, na.rm = TRUE),
            median.me = median(mse, na.rm = TRUE),
            median.nactive = median(nactive, na.rm = TRUE))

p1_me_nactive <- ggplot(data = df_me_nactive,
                        aes(x = mean.nactive, y = mean.me, color = Method, label = Method)) +
  geom_point(size = 2.1) +
  geom_text_repel(
    data = subset(df_me_nactive, mean.nactive < 100),
    nudge_x      = 2,
    nudge_y = 2,
    # size = 8,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_text_repel(
    data = subset(df_me_nactive, mean.nactive >= 100),
    nudge_x      = 2,
    nudge_y = 2,
    # size = 8,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_errorbar(aes(ymin = mean.me - sd.me, ymax = mean.me + sd.me), size = 1.1) +
  geom_errorbarh(aes(xmin = mean.nactive - sd.nactive, xmax = mean.nactive + sd.nactive), size = 1.1) +
  facet_rep_grid(p_overlap ~ eta_p, scales = "free",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_color_manual(values = cbbPalette[c(3,4,2)], guide = guide_legend(ncol=3)) +
  labs(x = "Number of active variables", y = "Mean Squared Error",
       title = "Mean Squared Error vs. Number of Active Variable (Mean +/- 1 SD) for Model with 1% Causal SNPs",
       subtitle = "Based on 200 simulations",
       caption = "") +
  # ylim(c(0,10)) +
  theme_box

p1_me_nactive






## ---- plot-prederror-nactive-sim ----

# not-used in manuscript
df_prederror_nactive <- DT[, c("Method", "structure", "p_overlap", "nactive", "prederror")] %>%
  group_by(Method, structure, p_overlap) %>%
  summarise(mean.prederror = mean(prederror, na.rm = TRUE), sd.prederror = sd(prederror, na.rm = TRUE),
            mean.nactive = mean(nactive, na.rm = TRUE), sd.nactive = sd(nactive, na.rm = TRUE))

p1_prederror_nactive <- ggplot(data = df_prederror_nactive, aes(x = mean.nactive, y = mean.prederror, color = Method, label = Method)) +
  geom_point(size = 2.1) +
  geom_text_repel(
    data = subset(df_prederror_nactive, mean.nactive < 100),
    nudge_x      = 5,
    nudge_y = 9,
    # size = 8,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_text_repel(
    data = subset(df_prederror_nactive, mean.nactive >= 100),
    nudge_x      = 5,
    nudge_y = 9,
    # size = 8,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_errorbar(aes(ymin = mean.prederror - sd.prederror, ymax = mean.prederror + sd.prederror), size = 1.1) +
  geom_errorbarh(aes(xmin = mean.nactive - sd.nactive, xmax = mean.nactive + sd.nactive), size = 1.1) +
  facet_rep_grid(p_overlap ~ structure, scales = "free",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  # scale_color_brewer(palette = "Dark2")+
  scale_color_manual(values = cbbPalette[c(7,3,4)], guide = guide_legend(ncol=3)) +
  labs(x = "Number of active variables", y = "Prediction Error",
       title = "Model Error vs. Number of Active Variable (Mean +/- 1 SD)",
       subtitle = "Based on 200 simulations",
       caption = "") +
  theme_box

p1_prederror_nactive


## ---- plot-tpr-sim ----
# not used in manuscript
p1_tpr <- ggplot(DT, aes(method, tpr, fill = method)) +
  ggplot2::geom_boxplot() +
  facet_rep_wrap(~scen, scales = "free", ncol = 2,
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(12, "Paired")[-11], guide=guide_legend(ncol=3)) +
  labs(x="", y="True Positive Rate",
       title="True Positive Rate",
       subtitle="Based on 200 simulations",
       caption="") +
  theme_ipsum_rc() + theme(legend.position = "right", axis.text.x = element_text(angle = 25, hjust = 1),
                           legend.text=element_text(size=14),
                           strip.text = element_text(size=14))

reposition_legend(p1_tpr, 'center', panel='panel-2-3')


## ---- plot-fpr-sim ----
# not used in manuscript
# DT[, table(time)]
DT[fpr==1, fpr:=NA] # this is skewing the plots, so we remove them and set to NA
p1_fpr <- ggplot(DT, aes(method, fpr, fill = method)) +
  ggplot2::geom_boxplot() +
  facet_rep_wrap(~scen, scales = "free", ncol = 2,
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(12, "Paired")[-11], guide=guide_legend(ncol=3)) +
  labs(x="", y="False Positive Rate",
       title="False Positive Rate",
       subtitle="Based on 200 simulations",
       caption="") +
   theme_ipsum_rc() + theme(legend.position = "right", axis.text.x = element_text(angle = 25, hjust = 1),
                            legend.text=element_text(size=14),
                            strip.text = element_text(size=14))


reposition_legend(p1_fpr, 'center', panel='panel-2-3')

## ---- plot-nactive-sim-null-model ----

p1_nactive <- ggplot(DT[structure == "block"][p_causal == "Null model"][p_overlap == "No causal SNPs in Kinship"], aes(Method, nactive, fill = Method)) +
  ggplot2::geom_boxplot() +
  facet_rep_grid(p_overlap ~ eta_p, scales = "fixed",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,3,4)]) +
  labs(x = "", y = "",
       title = "Number of Active Variables for Null Model",
       subtitle = "Based on 200 simulations",
       caption = "a variable is active if its estimated coefficient is non-zero") +
  theme_box

p1_nactive


## ---- plot-nactive-sim-1p-causal ----

p1_nactive <- ggplot(DT[structure == "block"][p_causal != "Null model"], aes(Method, nactive, fill = Method)) +
  ggplot2::geom_boxplot() +
  facet_rep_grid(p_overlap ~ eta_p, scales = "fixed",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,3,4)]) +
  labs(x = "", y = "",
       title = "Number of Active Variables for Model with 1% Causal SNPs",
       subtitle = "Based on 200 simulations",
       caption = "a variable is active if its estimated coefficient is non-zero") +
  theme_box

p1_nactive

## ---- plot-eta-sim-null-model ----

dummy2 <- data.frame(eta_p = c("10% Heritability", "30% Heritability"), Z = c(0.1, 0.3))
p1_eta <- ggplot(DT[structure == "block"][p_causal == "Null model"][Method %in% c("twostep","ggmix")][p_overlap == "No causal SNPs in Kinship"],
                 aes(Method, eta, fill = Method)) +
  ggplot2::geom_boxplot() +
  facet_rep_grid(p_overlap ~ eta_p, scales = "fixed",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,4)]) +
  labs(x = "", y = TeX("$\\hat{\\eta}$"),
       title = TeX("Estimated Heritability for the Null Model"),
       subtitle = "Based on 200 simulations",
       caption = "horizontal dashed line is the true value") +
  theme_ipsum_rc(axis_title_just = "bt",axis_title_size = 20, axis = T, ticks = F, axis_col = "black") +
  theme(legend.position = "none",title = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        strip.text = element_text(size = 18)) +
  geom_hline(data = dummy2, aes(yintercept = Z), linetype = 2, col = "#2f4f4f")

p1_eta

## ---- plot-eta-sim-1p-causal ----

dummy2 <- data.frame(eta_p = c("10% Heritability", "30% Heritability"), Z = c(0.1, 0.3))

p1_eta <- ggplot(DT[structure == "block"][p_causal != "Null model"][Method %in% c("twostep","ggmix")],
                 aes(Method, eta, fill = Method)) +
  ggplot2::geom_boxplot() +
  facet_rep_grid(p_overlap ~ eta_p, scales = "free",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,4)]) +
  labs(x = "", y = TeX("$\\hat{\\eta}$"),
       title = TeX("Estimated Heritability for the Model with 1% Causal SNPs"),
       subtitle = "Based on 200 simulations",
       caption = "horizontal dashed line is the true value") +
  theme_ipsum_rc(axis_title_just = "bt",axis_title_size = 20, axis = T, ticks = F, axis_col = "black") +
  theme(legend.position = "none",title = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        strip.text = element_text(size = 18)) +
geom_hline(data = dummy2, aes(yintercept = Z), linetype = 2, col = "#2f4f4f")

p1_eta



## ---- plot-errorvar-sim-null-model ----

dummy2 <- data.frame(eta_p = c("10% Heritability", "30% Heritability"), Z = c((1 - 0.1), (1 - 0.3)))

p1_errorvar <- ggplot(DT[structure == "block"][p_causal == "Null model"][p_overlap == "No causal SNPs in Kinship"],
                    aes(Method, errorvar, fill = Method)) +
  ggplot2::geom_boxplot() +
  facet_rep_grid(p_overlap ~ eta_p, scales = "fixed",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,3,4)]) +
  labs(x = "", y = "Error Variance",
       title = TeX("Estimated Error Variance for the Null Model"),
       subtitle = "Based on 200 simulations",
       caption = "horizontal dashed line is the true value") +
  theme_ipsum_rc(axis_title_just = "bt",axis_title_size = 20, axis = T, ticks = F, axis_col = "black") +
  theme(legend.position = "none",title = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        strip.text = element_text(size = 18)) +
  geom_hline(data = dummy2, aes(yintercept = Z), linetype = 2, col = "#2f4f4f")

p1_errorvar


## ---- plot-errorvar-sim-1p-causal ----

dummy2 <- data.frame(eta_p = c("10% Heritability", "30% Heritability"), Z = c((1 - 0.1), (1 - 0.3)))
p1_errorvar <- ggplot(DT[structure == "block"][p_causal != "Null model"],
                    aes(Method, errorvar, fill = Method)) +
  ggplot2::geom_boxplot() +
  facet_rep_grid(p_overlap ~ eta_p, scales = "free",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,3,4)]) +
  labs(x = "", y = "Error Variance",
       title = TeX("Estimated Error Variance for the Model with 1% Causal SNPs"),
       subtitle = "Based on 200 simulations",
       caption = "horizontal dashed line is the true value") +
  theme_ipsum_rc(axis_title_just = "bt",axis_title_size = 20, axis = T, ticks = F, axis_col = "black") +
  theme(legend.position = "none",title = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        strip.text = element_text(size = 18)) +
  geom_hline(data = dummy2, aes(yintercept = Z), linetype = 2, col = "#2f4f4f")

p1_errorvar


## ---- plot-tpr-fpr-sim-null-model ----

df_tpr_fpr <- DT[structure == "block"][p_causal == "Null model"][p_overlap == "No causal SNPs in Kinship"][, c("Method", "eta_p", "p_overlap", "tpr", "fpr")] %>%
  group_by(Method, eta_p, p_overlap) %>%
  summarise(mean.tpr = mean(tpr, na.rm = TRUE), sd.tpr = sd(tpr, na.rm = TRUE),
            mean.fpr = mean(fpr, na.rm = TRUE), sd.fpr = sd(fpr, na.rm = TRUE))

p1_tpr_fpr <- ggplot(data = df_tpr_fpr, aes(x = mean.fpr, y = mean.tpr, color = Method, label = Method)) +
  geom_point(size = 2.1) +
  # geom_text_repel() +
  geom_text_repel(
    data = subset(df_tpr_fpr, mean.fpr < 0.02),
    nudge_x      = 0.002,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_text_repel(
    data = subset(df_tpr_fpr, mean.fpr > 0.02),
    nudge_x      = 0.005,
    nudge_y = 0.015,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_errorbar(aes(ymin = mean.tpr - sd.tpr, ymax = mean.tpr + sd.tpr), size = 1.1) +
  geom_errorbarh(aes(xmin = mean.fpr - sd.fpr, xmax = mean.fpr + sd.fpr), size = 1.1) +
  facet_rep_grid(p_overlap ~ eta_p, scales = "free",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_color_manual(values = cbbPalette[c(7,3,4)], guide=guide_legend(ncol=3)) +
  labs(x="False positive rate", y="True positive rate",
       title="True Positive Rate vs. False Positive Rate (Mean +/- 1 SD) for the Null Model",
       subtitle="Based on 200 simulations",
       caption="") +
theme_box

p1_tpr_fpr



## ---- plot-tpr-fpr-sim-1p-causal ----

df_tpr_fpr <- DT[structure == "block"][p_causal != "Null model"][, c("Method", "eta_p", "p_overlap", "tpr", "fpr")] %>%
  group_by(Method, eta_p, p_overlap) %>%
  summarise(mean.tpr = mean(tpr, na.rm = TRUE), sd.tpr = sd(tpr, na.rm = TRUE),
            mean.fpr = mean(fpr, na.rm = TRUE), sd.fpr = sd(fpr, na.rm = TRUE))

p1_tpr_fpr <- ggplot(data = df_tpr_fpr, aes(x = mean.fpr, y = mean.tpr, color = Method, label = Method)) +
  geom_point(size = 2.1) +
  # geom_text_repel() +
  geom_text_repel(
    data = subset(df_tpr_fpr, mean.fpr < 0.02),
    nudge_x      = 0.002,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_text_repel(
    data = subset(df_tpr_fpr, mean.fpr > 0.02),
    nudge_x      = 0.005,
    nudge_y = 0.015,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_errorbar(aes(ymin = mean.tpr - sd.tpr, ymax = mean.tpr + sd.tpr), size = 1.1) +
  geom_errorbarh(aes(xmin = mean.fpr - sd.fpr, xmax = mean.fpr + sd.fpr), size = 1.1) +
  facet_rep_grid(p_overlap ~ eta_p, scales = "free",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_color_manual(values = cbbPalette[c(7,3,4)], guide=guide_legend(ncol=3)) +
  labs(x="False positive rate", y="True positive rate",
       title="True Positive Rate vs. False Positive Rate (Mean +/- 1 SD) for the Model with 1% Causal SNPs",
       subtitle="Based on 200 simulations",
       caption="") +
  theme_box

p1_tpr_fpr



## ---- plot-mse-sim ----

p1_mse <- ggplot(DT[Method %in% c("lasso","ggmix")], aes(Method, mse, fill = Method)) +
  ggplot2::geom_boxplot() +
  facet_rep_grid(p_overlap ~ structure, scales = "free",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,3,4)]) +
  labs(x = "", y = "",
       title = "Mean Squared Error",
       subtitle = "Based on 200 simulations",
       caption = "") +
  theme_box

p1_mse



## ---- plot-runtime-sim-null-model ----

# dummy2 <- data.frame(eta_p = c("10% Heritability", "30% Heritability"), Z = c((1 - 0.1), (1 - 0.3)))
p1_errorvar <- ggplot(DT[structure == "block"][p_causal == "Null model"][p_overlap == "No causal SNPs in Kinship"],
                      aes(Method, log(time), fill = Method)) +
  ggplot2::geom_boxplot() +
  facet_rep_grid(p_overlap ~ eta_p, scales = "free",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,3,4)]) +
  labs(x = "", y = "log run time (seconds)",
       title = TeX("Log Run Time (seconds) for the Null Model"),
       subtitle = "Based on 200 simulations") +
  theme_ipsum_rc(axis_title_just = "bt",axis_title_size = 20, axis = T, ticks = F, axis_col = "black") +
  theme(legend.position = "none",title = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        strip.text = element_text(size = 18)) #+
# geom_hline(data = dummy2, aes(yintercept = Z), linetype = 2, col = "#2f4f4f")

p1_errorvar




## ---- plot-runtime-sim-1p-causal ----

# dummy2 <- data.frame(eta_p = c("10% Heritability", "30% Heritability"), Z = c((1 - 0.1), (1 - 0.3)))
p1_errorvar <- ggplot(DT[structure == "block"][p_causal != "Null model"],
                      aes(Method, log(time), fill = Method)) +
  ggplot2::geom_boxplot() +
  facet_rep_grid(p_overlap ~ eta_p, scales = "free",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,3,4)]) +
  labs(x = "", y = "log run time (seconds)",
       title = TeX("Log Run Time (seconds) for the Model with 1% Causal SNPs"),
       subtitle = "Based on 200 simulations") +
  theme_ipsum_rc(axis_title_just = "bt",axis_title_size = 20, axis = T, ticks = F, axis_col = "black") +
  theme(legend.position = "none",title = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        strip.text = element_text(size = 18)) #+
  # geom_hline(data = dummy2, aes(yintercept = Z), linetype = 2, col = "#2f4f4f")

p1_errorvar


## ---- plot-fpr-tpr-boxplot-sim ----
# not-used
# DT_tpr_fpr <- melt(DT[, .(method, tpr, fpr, scen)], id.vars = c("method","scen"))
# p1_tpr <- ggplot(DT_tpr_fpr, aes(x = method, y = value, fill = variable)) +
#   ggplot2::geom_boxplot() +
#   facet_rep_wrap(~scen, scales = "free", ncol = 2,
#                  repeat.tick.labels = 'left',
#                  labeller = as_labeller(appender,
#                                         default = label_parsed)) +
#   scale_fill_manual(values=RColorBrewer::brewer.pal(11, "Paired"), guide=guide_legend(ncol=3)) +
#   labs(x="", y="Test Set MSE",
#        title="Simulation Study Results: Test Set MSE",
#        subtitle="Based on 200 simulations",
#        caption="") +
#   theme_ipsum_rc() + theme(legend.position = "right", axis.text.x = element_text(angle = 25, hjust = 1))
# reposition_legend(p1_tpr, 'center', panel='panel-2-3')


