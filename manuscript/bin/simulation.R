
## ---- simulation-results ----

df <- readRDS(here::here("manuscript/data/simulation_results.rds"))

df <- df %>% separate(Model,
  into = c(
    "simnames", "b0", "beta_mean", "eta_p", "Fst", "geography", "k", "n",
    "pkinship", "ptest", "percentcausal",
    "percentoverlap", "s", "sigma2_p"
  ),
  sep = "/"
)

DT <- data.table::as.data.table(df, stringsAsFactors = FALSE)

DT[, p_overlap := case_when(
  percentoverlap == "percent_overlap_0" ~ "No causal SNPs in Kinship",
  percentoverlap == "percent_overlap_100" ~ "All causal SNPs in Kinship"
)]
DT[, p_causal := case_when(
  percentcausal == "percent_causal_0" ~ "Null model",
  percentcausal == "percent_causal_0.01" ~ "1% of SNPs are causal"
)]

DT[, p_causal := factor(p_causal, levels = c("Null model", "1% of SNPs are causal"))]


## --PATCH- # to keep code with as little change as possible, im renaming 1D to block
## because everything below is for "block". Even though, the may7th results are for 1D structure.
DT[geography == "geography_1d", structure := "block"]
DT[, structure := factor(structure, levels = c("block"))]
DT[, eta_p := case_when(
  eta_p == "eta_0.1" ~ "10% Heritability",
  eta_p == "eta_0.3" ~ "30% Heritability"
)]
DT[Method == "twostepYVCCV", Method := "twostep"]
DT[Method == "lassoCV", Method := "lasso"]
DT[Method == "ggmixHDBIC", Method := "ggmix"]
DT[, Method := factor(Method, levels = c("twostep", "lasso", "ggmix"))]
DT[Method == "twostep", errorvar := sigma2]
DT[Method == "twostep", tau := eta]
DT[Method == "twostep", eta := tau / (tau + sigma2)]

DT[, me2 := (1 / 1000) * me^2]
appender <- function(string) latex2exp::TeX(paste(string))

# theme_box <- theme_ipsum_rc(axis_title_just = "bt",axis_title_size = 16, axis = T, ticks = F, axis_col = "black") +
#   theme(legend.position = "bottom",title = element_text(size = 20),
#         axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16),
#         axis.text.y = element_text(size = 16),
#         legend.text = element_text(size = 16), legend.title = element_text(size = 16),
#         strip.text = element_text(size = 18))


## ---- table-of-results ----

DT$RMSE <- sqrt(DT$mse)

tt <- DT %>%
  group_by(Method, eta_p, p_overlap, p_causal) %>%
  summarise(
    RMSE = qwraps2::mean_sd(RMSE, denote_sd = "paren"),
    TPR = qwraps2::mean_sd(tprFPR5, denote_sd = "paren"),
    Heritability = qwraps2::mean_sd(eta, denote_sd = "paren"),
    Errorvar = qwraps2::mean_sd(errorvar, denote_sd = "paren"),
    Estimationerror = qwraps2::mean_sd(estimationerror, denote_sd = "paren"),
    Nactive = qwraps2::median_iqr(nactive, digits = 0)
  )

tt[tt$Method == "lasso", "Heritability"] <- "--"

pt <- tt %>%
  pivot_longer(
    cols = c("RMSE", "TPR", "Heritability", "Errorvar", "Estimationerror", "Nactive"),
    names_to = "metric"
  ) %>%
  unite(col = "type", p_causal, p_overlap, eta_p) %>%
  mutate(
    type = factor(type,
      levels = c(
        "Null model_No causal SNPs in Kinship_10% Heritability",
        "Null model_No causal SNPs in Kinship_30% Heritability",
        "Null model_All causal SNPs in Kinship_10% Heritability",
        "Null model_All causal SNPs in Kinship_30% Heritability",
        "1% of SNPs are causal_No causal SNPs in Kinship_10% Heritability",
        "1% of SNPs are causal_No causal SNPs in Kinship_30% Heritability",
        "1% of SNPs are causal_All causal SNPs in Kinship_10% Heritability",
        "1% of SNPs are causal_All causal SNPs in Kinship_30% Heritability"
      )
    ),
    metric = factor(metric,
      labels = c("TPR at FPR=5%", "Model Size", "RMSE", "Estimation Error", "Error Variance", "Heritability"),
      levels = c("TPR", "Nactive", "RMSE", "Estimationerror", "Errorvar", "Heritability")
    )
  ) %>%
  pivot_wider(id_cols = c("metric", "Method"), names_from = "type", values_from = "value") %>%
  arrange(metric, Method) %>%
  select(
    "metric", "Method", "Null model_No causal SNPs in Kinship_10% Heritability",
    "Null model_No causal SNPs in Kinship_30% Heritability",
    "Null model_All causal SNPs in Kinship_10% Heritability",
    "Null model_All causal SNPs in Kinship_30% Heritability",
    "1% of SNPs are causal_No causal SNPs in Kinship_10% Heritability",
    "1% of SNPs are causal_No causal SNPs in Kinship_30% Heritability",
    "1% of SNPs are causal_All causal SNPs in Kinship_10% Heritability",
    "1% of SNPs are causal_All causal SNPs in Kinship_30% Heritability"
  )


## ---- print-sim-table ----

kable(pt,
  booktabs = T, align = c("l", "l", "c", "c", "c", "c", "c", "c", "c", "c"),
  caption = c("Mean (standard deviation) from 200 simulations stratified by the number of causal SNPs (null, 1\\%), the overlap between causal SNPs and kinship matrix (no overlap, all causal SNPs in kinship), and true heritability (10\\%, 30\\%).
                  For all simulations, sample size is $n=1000$, the number of covariates is $p=5000$, and the number of SNPs used to estimate the kinship matrix is $k=10000$.
                  TPR at FPR=5\\% is the true positive rate at a fixed false positive rate of 5\\%.
                  Model Size ($|\\widehat{S}_{\\hat{\\lambda}}|$) is the number of selected variables in the training set using the high-dimensional BIC for \\texttt{ggmix} and 10-fold cross validation for \\texttt{lasso} and \\texttt{twostep}.
                  RMSE is the root mean squared error on the test set.
                  Estimation error is the squared distance between the estimated and true effect sizes.
                  Error variance ($\\sigma^2$) for \\texttt{twostep} is estimated from an intercept only LMM with a single random effect and is modeled explicitly in \\ggmix. For the \\texttt{lasso} we use $\\frac{1}{n - |\\widehat{S}_{\\hat{\\lambda}}|} ||\\bY - \\bX \\widehat{\\bbeta}_{\\hat{\\lambda}}||_2^2$~\\citep{reid2016study} as an estimator for $\\sigma^2$.
                  Heritability ($\\eta$) for \\texttt{twostep} is estimated as $\\sigma_g^2 / (\\sigma_g^2 + \\sigma_e^2)$ from an intercept only LMM with a single random effect where $\\sigma_g^2$ and $\\sigma_e^2$ are the variance components for the random effect and error term, respectively. $\\eta$ is explictly modeled in \\ggmix. There is no positive way to calculate $\\eta$ for the \\texttt{lasso} since we are using a PC adjustment."),
  col.names = c("Metric", "Method", rep(c("10%", "30%"), 4))
) %>%
  kable_styling(
    latex_options = "striped", full_width = TRUE,
    position = "center", font_size = 7, stripe_index = c(1:3, 7:9, 13:15)
  ) %>%
  add_header_above(c(" ", " ", "No overlap" = 2, "All causal SNPs\nin kinship" = 2, "No overlap" = 2, "All causal SNPs\nin kinship" = 2)) %>%
  add_header_above(c(" ", " ", "Null model" = 4, "1% Causal SNPs" = 4)) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 1, latex_hline = "major", valign = "middle") %>%
  footnote(
    general = "Median (Inter-quartile range) is given for Model Size.",
    general_title = "Note:"
  )


## ---- plot-kinship-sim ----

dat <- lapply(
  list("1d"),
  function(i) {
    ggmix::gen_structured_model(
      n = 1200,
      p_design = 5000,
      p_kinship = 10000,
      geography = i,
      percent_causal = 0.01,
      percent_overlap = "100",
      k = 10, s = 0.5, Fst = 0.1,
      b0 = 0,
      eta = 0.1, sigma2 = 1,
      train_tune_test = c(0.99, .005, 0.005)
    )
  }
)

popkin::plot_popkin(kinship = list(dat[[1]]$kinship))


## ---- plot-pc-sim ----

xlabs <- "1st principal component"
ylabs <- "2nd principal component"

plot(dat[[1]]$PC[, 1], dat[[1]]$PC[, 2],
  pch = 19, col = rep(RColorBrewer::brewer.pal(10, "Paired"), each = table(dat[[1]]$subpops)[1]),
  xlab = xlabs, ylab = ylabs, cex = 0.5
)
