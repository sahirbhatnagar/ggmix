## ---- simulation-results ----

# df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/june_29_2018_results.rds")
# df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/july_1_2018_results.rds")
# df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/july_1_2018_results_with_twostepY.rds")
# df <- df %>% separate(Model,
#                       into = c("simnames","b0","eta_p","Fst","geography","k","n",
#                                "pkinship","ptest","percentcausal",
#                                "percentoverlap","s","sigma2_p"),
#                       sep = "/")

df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/july_12_2018_results_with_null_model_VC.rds")
# df <- readRDS("C:/Users/sahir/Documents/git_repositories/ggmix/simulation/simulation_results/july_12_2018_results_with_null_model_VC.rds")
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
DT[geography == "geography_ind", structure := "block"]
DT[geography == "geography_circ", structure := "circular"]
DT[geography == "geography_1d", structure := "1D"]
# DT[, table(geography, structure)]
DT[, structure := factor(structure, levels = c("block","1D","circular"))]
DT[, eta_p := case_when(eta_p == "eta_0.1" ~ "10% Heritability",
                        eta_p == "eta_0.5" ~ "50% Heritability")]
# DT[, table(eta_p)]
# use twostepY, which compares to the original Y
# DT[, table(Method)]
DT <- DT[Method != "twostep"]
DT[Method == "twostepY", Method := "twostep"]
# DT[, Method := factor(Method, levels = c("twostep","twostepY","lasso","ggmix"))]
DT[, Method := factor(Method, levels = c("twostep","lasso","ggmix"))]
# DT[, table(percentcausal,p_causal)]
# DT[, Method := factor(Method, levels = c("twostep","twostepY","lasso","ggmix"))]
# DT[, table(Method)]

DT[Method == "twostep", errorvar := sigma2]
# DT[Method == "twostep"]
appender <- function(string) latex2exp::TeX(paste(string))

theme_box <- theme_ipsum_rc(axis_title_just = "bt",axis_title_size = 20, axis = T, ticks = F, axis_col = "black") +
  theme(legend.position = "bottom",title = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        strip.text = element_text(size = 18))


## ---- plot-kinship-sim ----

dat <- lapply(list("ind","1d","circ"),
              function(i)
                ggmix::gen_structured_model(n = 1000,
                                            p_test = 5000,
                                            p_kinship = 10000,
                                            geography = i,
                                            percent_causal = 0.01,
                                            percent_overlap = "100",
                                            k = 5, s = 0.5, Fst = 0.1,
                                            b0 = 0,
                                            eta = 0.1, sigma2 = 1)
)
# str(dat)
# dev.off()

par(omi = c(0.3,0.3,0.3,0.3))
# popkin::plotPopkin(x = list(dat[[1]]$kin, dat[[2]]$kin, dat[[3]]$kin),
#                    titles = c("block", "1D", "circular"),
#                    marPad = 0.05)
popkin::plotPopkin(x = list(dat[[1]]$kin),
                   titles = c("Empirical Kinship Matrix with Block Structure"),
                   marPad = 0.05)


## ---- plot-pc-sim ----

xlabs <- "1st principal component"
ylabs <- "2nd principal component"

# par(mfrow = c(1,3))
plot(dat[[1]]$x_lasso[,5001],dat[[1]]$x_lasso[,5002],
     pch = 19, col = rep(RColorBrewer::brewer.pal(5,"Paired"), each = 200),
     xlab = xlabs, ylab = ylabs,
     main = "Block Structure")
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
       caption = TeX("$\\eta$ = 10%")) +
  theme_box

p1_cs
# reposition_legend(p1_mse, 'center', panel='panel-2-3')

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
       caption = TeX("$\\eta$ = 10%")) +
  theme_box


p1_cs


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
            mean.nactive = mean(nactive, na.rm = TRUE), sd.nactive = sd(nactive, na.rm = TRUE))

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

dummy2 <- data.frame(eta_p = c("10% Heritability", "50% Heritability"), Z = c(0.1, 0.5))
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

dummy2 <- data.frame(eta_p = c("10% Heritability", "50% Heritability"), Z = c(0.1, 0.5))

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

dummy2 <- data.frame(eta_p = c("10% Heritability", "50% Heritability"), Z = c((1 - 0.1), (1 - 0.5)))

p1_errorvar <- ggplot(DT[Method %in% c("twostep","ggmix")][structure == "block"][p_causal == "Null model"][p_overlap == "No causal SNPs in Kinship"],
                    aes(Method, errorvar, fill = Method)) +
  ggplot2::geom_boxplot() +
  facet_rep_grid(p_overlap ~ eta_p, scales = "fixed",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,4)]) +
  labs(x = "", y = TeX("$\\hat{\\sigma}^2$"),
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

dummy2 <- data.frame(eta_p = c("10% Heritability", "50% Heritability"), Z = c((1 - 0.1), (1 - 0.5)))
p1_errorvar <- ggplot(DT[Method %in% c("twostep","ggmix")][structure == "block"][p_causal != "Null model"],
                    aes(Method, errorvar, fill = Method)) +
  ggplot2::geom_boxplot() +
  facet_rep_grid(p_overlap ~ eta_p, scales = "free",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,4)]) +
  labs(x = "", y = TeX("$\\hat{\\sigma}^2$"),
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

p1_mse <- ggplot(DT, aes(Method, mse, fill = Method)) +
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


