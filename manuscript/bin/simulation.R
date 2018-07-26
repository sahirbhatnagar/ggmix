## ---- simulation-results ----

# df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/june_29_2018_results.rds")
# df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/july_1_2018_results.rds")
# df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/july_1_2018_results_with_twostepY.rds")
# df <- df %>% separate(Model,
#                       into = c("simnames","b0","eta_p","Fst","geography","k","n",
#                                "pkinship","ptest","percentcausal",
#                                "percentoverlap","s","sigma2_p"),
#                       sep = "/")

df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/july_12_2018_results_with_null_model.rds")
df <- df %>% separate(Model,
                      into = c("simnames","b0","beta_mean","eta_p","Fst","geography","k","n",
                               "pkinship","ptest","percentcausal",
                               "percentoverlap","s","sigma2_p"),
                      sep = "/")

DT <- data.table::as.data.table(df, stringsAsFactors = FALSE)

# DT[, table(percentoverlap, p_overlap)]
DT[, p_overlap := case_when(percentoverlap == "percent_overlap_0" ~ "No causal SNPs in Kinship",
                            percentoverlap == "percent_overlap_100" ~ "All causal SNPs in Kinship")]

DT[geography == "geography_ind", structure := "block"]
DT[geography == "geography_circ", structure := "circular"]
DT[geography == "geography_1d", structure := "1D"]
# DT[, table(geography, structure)]
DT[, structure := factor(structure, levels = c("block","1D","circular"))]
DT[, Method := factor(Method, levels = c("twostep","twostepY","lasso","ggmix"))]

# DT[, Method := factor(Method, levels = c("twostep","twostepY","lasso","ggmix"))]
# DT[, table(Method)]
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
popkin::plotPopkin(x = list(dat[[1]]$kin, dat[[2]]$kin, dat[[3]]$kin),
                   titles = c("block", "1D", "circular"),
                   marPad = 0.05)


## ---- plot-pc-sim ----

xlabs <- "1st principal component"
ylabs <- "2nd principal component"

par(mfrow = c(1,3))
plot(dat[[1]]$x_lasso[,5001],dat[[1]]$x_lasso[,5002],
     pch = 19, col = rep(RColorBrewer::brewer.pal(5,"Paired"), each = 200),
     xlab = xlabs, ylab = ylabs,
     main = "Block Structure")
plot(dat[[2]]$x_lasso[,5001],dat[[2]]$x_lasso[,5002],
     pch = 19, col = rep(RColorBrewer::brewer.pal(5,"Paired"), each = 200),
     xlab = xlabs, ylab = ylabs,
     main = "1D Geography")
plot(dat[[3]]$x_lasso[,5001],dat[[3]]$x_lasso[,5002],
     pch = 19, col = rep(RColorBrewer::brewer.pal(5,"Paired"), each = 200),
     xlab = xlabs, ylab = ylabs,
     main = "Circular Geography")

## ---- plot-kinship-pc-combined-sim ----

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

## ---- plot-correct-sparsity-sim ----

p1_cs <- ggplot(DT[percentcausal=="percent_causal_0"], aes(Method, correct_sparsity, fill = Method)) +
    ggplot2::geom_boxplot() +
    facet_rep_grid(p_overlap ~ structure, scales = "fixed",
                   repeat.tick.labels = 'left',
                   labeller = as_labeller(appender,
                                          default = label_parsed)) +
    scale_fill_manual(values = cbbPalette[c(7,3,4,2)]) +
    labs(x = "", y = "",
         title = "Correct Sparsity",
         subtitle = "Based on 200 simulations",
         caption = "") +
  theme_box



p1_cs
# reposition_legend(p1_mse, 'center', panel='panel-2-3')

## ---- plot-me-nactive-sim ----
DT[percentcausal=="percent_causal_0.01"]
df_me_nactive <- DT[percentcausal=="percent_causal_0.01", c("Method", "structure", "p_overlap", "nactive", "me")] %>%
  group_by(Method, structure, p_overlap) %>%
  summarise(mean.me = mean(me, na.rm = TRUE), sd.me = sd(me, na.rm = TRUE),
         mean.nactive = mean(nactive, na.rm = TRUE), sd.nactive = sd(nactive, na.rm = TRUE))

p1_me_nactive <- ggplot(data = df_me_nactive, aes(x = mean.nactive, y = mean.me, color = Method, label = Method)) +
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
  facet_rep_grid(p_overlap ~ structure, scales = "free",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_color_manual(values = cbbPalette[c(7,3,4,2)], guide = guide_legend(ncol=3)) +
  labs(x = "Number of active variables", y = "Model Error",
       title = "Model Error vs. Number of Active Variable (Mean +/- 1 SD)",
       subtitle = "Based on 200 simulations",
       caption = "") +
  theme_box

p1_me_nactive


## ---- plot-mse-nactive-sim ----

df_mse_nactive <- DT[percentcausal=="percent_causal_0", c("Method", "structure", "p_overlap", "nactive", "mse")] %>%
  group_by(Method, structure, p_overlap) %>%
  summarise(mean.mse = mean(mse, na.rm = TRUE), sd.mse = sd(mse, na.rm = TRUE),
            mean.nactive = mean(nactive, na.rm = TRUE), sd.nactive = sd(nactive, na.rm = TRUE))

p1_mse_nactive <- ggplot(data = as.data.table(df_mse_nactive)[Method!="twostepY"],
                         aes(x = mean.nactive, y = mean.mse, color = Method, label = Method)) +
  geom_point(size = 2.1) +
  # geom_text_repel(
  #   data = subset(df_mse_nactive, mean.nactive < 100),
  #   nudge_x      = 5,
  #   nudge_y = 9,
  #   # size = 8,
  #   direction    = "y",
  #   hjust        = 0,
  #   segment.size = 0.2
  # ) +
  # geom_text_repel(
  #   data = subset(df_mse_nactive, mean.nactive >= 100),
  #   nudge_x      = 5,
  #   nudge_y = 9,
  #   # size = 8,
  #   direction    = "y",
  #   hjust        = 0,
  #   segment.size = 0.2
  # ) +
  geom_errorbar(aes(ymin = mean.mse - sd.mse, ymax = mean.mse + sd.mse), size = 1.1) +
  geom_errorbarh(aes(xmin = mean.nactive - sd.nactive, xmax = mean.nactive + sd.nactive), size = 1.1) +
  facet_rep_grid(p_overlap ~ structure, scales = "free",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_color_manual(values = cbbPalette[c(7,3,4,2)], guide = guide_legend(ncol=3)) +
  labs(x = "Number of active variables", y = "Mean Squared Error",
       title = "Mean Squared Error vs. Number of Active Variable (Mean +/- 1 SD)",
       subtitle = "Based on 200 simulations",
       caption = "") +
  theme_box

p1_mse_nactive


## ---- plot-prederror-nactive-sim ----

# not-used
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
# not used
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
#not used
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

## ---- plot-nactive-sim ----

p1_nactive <- ggplot(DT, aes(Method, nactive, fill = Method)) +
  ggplot2::geom_boxplot() +
  facet_rep_grid(p_overlap ~ structure, scales = "fixed",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(7,3,4)]) +
  labs(x = "", y = "",
       title = "Number of Active Variables",
       subtitle = "Based on 200 simulations",
       caption = "") +
  theme_box

p1_nactive


## ---- plot-eta-sim ----

DT$eta_p
p1_eta <- ggplot(DT[Method=="ggmix"][percentcausal=="percent_causal_0"][eta_p=="eta_0.1"],
                 aes(Method, eta, fill = Method)) +
  ggplot2::geom_boxplot() +
  facet_rep_grid(p_overlap ~ structure, scales = "fixed",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(4)]) +
  labs(x = "", y = TeX("$\\hat{\\eta}$"),
       title = TeX("Estimated $\\eta$ Parameter in the ggmix Model"),
       subtitle = "Based on 200 simulations",
       caption = "") +
  theme_ipsum_rc(axis_title_just = "bt",axis_title_size = 20, axis = T, ticks = F, axis_col = "black") +
  theme(legend.position = "none",title = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        strip.text = element_text(size = 18))

p1_eta

## ---- plot-sigma2-sim ----

p1_sigma2 <- ggplot(DT[Method=="ggmix"][percentcausal=="percent_causal_0"][eta_p=="eta_0.5"],
                    aes(Method, sigma2, fill = Method)) +
  ggplot2::geom_boxplot() +
  facet_rep_grid(p_overlap ~ structure, scales = "fixed",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values = cbbPalette[c(4)]) +
  labs(x = "", y = TeX("$\\hat{\\sigma}^2$"),
       title = TeX("Estimated $\\sigma^2$ Parameter in the ggmix Model"),
       subtitle = "Based on 200 simulations",
       caption = "") +
  theme_ipsum_rc(axis_title_just = "bt",axis_title_size = 20, axis = T, ticks = F, axis_col = "black") +
  theme(legend.position = "none",title = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        strip.text = element_text(size = 18))

p1_sigma2

## ---- plot-tpr-fpr-sim ----

df_tpr_fpr <- DT[, c("Method", "structure", "p_overlap", "tpr", "fpr")] %>%
  group_by(Method, structure, p_overlap) %>%
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
  facet_rep_grid(p_overlap ~ structure, scales = "free",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_color_manual(values = cbbPalette[c(7,3,4)], guide=guide_legend(ncol=3)) +
  labs(x="False positive rate", y="True positive rate",
       title="True Positive Rate vs. False Positive Rate (Mean +/- 1 SD)",
       subtitle="Based on 200 simulations",
       caption="") +
theme_box

p1_tpr_fpr

# reposition_legend(p1_tpr_fpr, 'center', panel='panel-2-3')

## ---- plot-errorvar-sim ----

df_errorvar <- DT[, c("Method", "structure", "p_overlap", "errorvar")] %>%
  group_by(p_overlap, structure, Method) %>%
  summarise(mean = mean(errorvar, na.rm = TRUE),
            sd = sd(errorvar, na.rm = TRUE)) #%>%
  # mutate(value = sprintf("%0.2f (%0.2f)", mean,sd))

p1_errorvar <- ggplot(data = df_errorvar, aes(Method, mean, color = Method, label = Method)) +
  geom_point(size = 2.1) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), size = 1.1, width = 0.51) +
  geom_hline(yintercept = (1-0.1)*1, linetype = 2, col = "#2f4f4f") +
  # geom_errorbarh(aes(xmin = mean.fpr - sd.fpr, xmax = mean.fpr + sd.fpr), size = 1.1) +
  facet_rep_grid(p_overlap ~ structure, scales = "free",
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_color_manual(values = cbbPalette[c(7,3,4)], guide=guide_legend(ncol=3)) +
  labs(x="", y="Error Variance",
       title="Error Variance (Mean +/- 1 SD)",
       subtitle="Based on 200 simulations. Dotted line represents the true error variance of 0.90.",
       caption="") +
  theme_box

p1_errorvar


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


