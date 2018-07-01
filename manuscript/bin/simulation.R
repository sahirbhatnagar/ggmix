## ---- simulation-results ----

df <- readRDS("/home/sahir/git_repositories/ggmix/simulation/simulation_results/june_29_2018_results.rds")
df <- df %>% separate(Model,
                      into = c("simnames","b0","eta","Fst","geography","k","n",
                               "pkinship","ptest","percentcausal",
                               "percentoverlap","s","sigma2"),
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
DT[, Method := factor(Method, levels = c("twostep","lasso","ggmix"))]
# DT[, table(Method)]
appender <- function(string) latex2exp::TeX(paste(string))


## ---- plot-correct-sparsity-sim ----

p1_cs <- ggplot(DT, aes(Method, correct_sparsity, fill = Method)) +
    ggplot2::geom_boxplot() +
    # gg_sy +
    # facet_rep_wrap(~scen, scales = "free", ncol = 2,
    #                repeat.tick.labels = 'left',
    #                labeller = as_labeller(appender,
    #                                       default = label_parsed)) +
    facet_rep_grid(p_overlap ~ structure, scales = "fixed",
                   repeat.tick.labels = 'left',
                   labeller = as_labeller(appender,
                                          default = label_parsed)) +
    scale_fill_manual(values = cbbPalette[c(7,3,4)]) +
    labs(x = "", y = "",
         title = "Correct Sparsity",
         subtitle = "Based on 200 simulations",
         caption = "") +
    # panel_border()+
    # background_grid()+
    theme_ipsum_rc() +
  theme(legend.position = "bottom",title = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        strip.text = element_text(size = 18))

# , legend.text=element_text(size=18)

p1_cs
# reposition_legend(p1_mse, 'center', panel='panel-2-3')

## ---- plot-me-nactive-sim ----

df_me_nactive <- DT[, c("Method", "structure", "p_overlap", "nactive", "me")] %>%
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
  # scale_color_brewer(palette = "Dark2")+
  scale_color_manual(values = cbbPalette[c(7,3,4)], guide = guide_legend(ncol=3)) +
  labs(x = "Number of active variables", y = "Model Error",
       title = "Model Error vs. Number of Active Variable (Mean +/- 1 SD)",
       subtitle = "Based on 200 simulations",
       caption = "") +
  theme_ipsum_rc(axis_title_just = "bt",axis_title_size = 16) +
    theme(legend.position = "bottom",title = element_text(size = 20),
          axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
          legend.text = element_text(size = 16), legend.title = element_text(size = 16),
          strip.text = element_text(size = 18))

p1_me_nactive


## ---- plot-prederror-nactive-sim ----

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
  theme_ipsum_rc(axis_title_just = "bt",axis_title_size = 16) +
  theme(legend.position = "bottom",title = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        strip.text = element_text(size = 18))

p1_prederror_nactive


## ---- plot-tpr-sim ----

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
  theme_ipsum_rc() +
  theme(legend.position = "bottom",title = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        strip.text = element_text(size = 18))

p1_nactive




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
  theme_ipsum_rc(axis_title_just = "bt",axis_title_size = 16) +
  theme(legend.position = "bottom",title = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        strip.text = element_text(size = 18))

p1_tpr_fpr

# reposition_legend(p1_tpr_fpr, 'center', panel='panel-2-3')


## ---- plot-fpr-tpr-boxplot-sim ----

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


