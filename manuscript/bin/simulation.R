## ---- simulation-results ----

df <- readRDS("/home/sahir/git_repositories/sail/my_sims/simulation_results/may_15_2018_results.rds")
# df <- readRDS("C:/Users/sahir/Documents/git_repositories/sail/my_sims/simulation_results/apr_25_2018_results.rds")
df <- df %>% separate(Model, into = c("simnames","betaE","corr","lambda.type","n","p","parameterIndex","SNR_2"),
                      sep = "/")

DT <- as.data.table(df, stringsAsFactors = FALSE)
DT <- DT[parameterIndex != "parameterIndex_4"]
# DT[parameterIndex=="parameterIndex_1", table(Method)] %>% names %>% dput
# DT[, table(Method)]

DT[Method=="Adaptivesail", Method := "adaptive sail"]
# DT[Method=="Adaptivesailweak", Method := "Asail weak"]
DT[Method=="Adaptivelasso", Method := "adaptive lasso"]
DT[Method=="sailweak", Method := "sail weak"]
# DT[Method=="sail", Method := "sail strong"]
# DT[Method=="Asail", Method := "Asail strong"]
DT[Method=="linearsail", Method := "linear sail"]
DT[, Method := droplevels(Method)]
# DT[, table(Method)]
# DT[, table(Method)] %>% names %>% dput
appender <- function(string) TeX(paste(string))

# DT[, method := factor(Method, levels = c("lasso","Alasso","lassoBT", "GLinternet", "HierBasis", "SPAM", "gamsel",
#                                        "sail strong", "Asail strong", "sail weak", "Asail weak"))]

DT[, method := factor(Method, levels = c("lasso","adaptive lasso","lassoBT", "GLinternet", "HierBasis", "SPAM", "gamsel",
                                         "sail", "adaptive sail",  "sail weak", "linear sail"))]

# DT[, table(method)]
# DT[, table(parameterIndex)]
DT[, scenario:= as.numeric(as.character(stringr::str_extract_all(parameterIndex, "\\d", simplify = T)))]
# DT[, table(scenario)]
DT[, scenario := replace(scenario, which(scenario==6), 4)]
DT[, scen := case_when(scenario==1 ~ "1a) Strong Hierarchy",
                       scenario==2 ~ "1b) Weak Hierarchy",
                       scenario==3 ~ "1c) Interactions Only",
                       scenario==4 ~ "2) Linear Effects",
                       scenario==5 ~ "3) Main Effects Only")]
DT[, scen := factor(scen, levels = c("1a) Strong Hierarchy", "1b) Weak Hierarchy","1c) Interactions Only","2) Linear Effects", "3) Main Effects Only"))]
# DT$scen %>% table
#Truth obeys strong hierarchy (parameterIndex = 1)
#Truth obeys weak hierarchy (parameterIndex = 2)
#Truth only has interactions (parameterIndex = 3)
#Truth is linear (parameterIndex = 4)
#Truth only has main effects (parameterIndex = 5)


## ---- plot-mse-sim ----

p1_mse <- ggplot(DT, aes(method, mse, fill = method)) +
    ggplot2::geom_boxplot() +
    # gg_sy +
    # facet_rep_wrap(~scen, scales = "free", ncol = 2,
    #                repeat.tick.labels = 'left',
    #                labeller = as_labeller(appender,
    #                                       default = label_parsed)) +
    facet_rep_wrap(~scen, scales = "free", ncol = 2,
                   repeat.tick.labels = 'left',
                   labeller = as_labeller(appender,
                                          default = label_parsed)) +
    scale_fill_manual(values=RColorBrewer::brewer.pal(12, "Paired")[-11], guide=guide_legend(ncol=3)) +
    # ggplot2::labs(y = "Test Set MSE", title = "") + xlab("") +
    labs(x="", y="Test Set MSE",
         title="Test Set MSE",
         subtitle="Based on 200 simulations",
         caption="") +
    # panel_border()+
    # background_grid()+
    theme_ipsum_rc() + theme(legend.position = "right", axis.text.x = element_text(angle = 25, hjust = 1),
                             legend.text=element_text(size=14),
                             strip.text = element_text(size=14))

# , legend.text=element_text(size=18)


reposition_legend(p1_mse, 'center', panel='panel-2-3')

## ---- plot-mse-nactive-sim ----


df_mse_nactive <- DT[, c("method","scen","mse","nactive")] %>%
  group_by(method, scen) %>%
  summarise(mean.mse = mean(mse, na.rm = TRUE), sd.mse = sd(mse, na.rm = TRUE),
         mean.nactive = mean(nactive, na.rm = TRUE), sd.nactive = sd(nactive, na.rm = TRUE)) %>%
  mutate(scen = case_when(scen == "1a) Strong Hierarchy" ~ "1a) Strong Hierarchy (|S_0| = 7)",
                          scen == "1b) Weak Hierarchy" ~ "1b) Weak Hierarchy (|S_0| = 5)",
                          scen == "1c) Interactions Only" ~ "1c) Interactions Only (|S_0| = 2)",
                          scen == "2) Linear Effects" ~ "2) Linear Effects (|S_0| = 7)",
                          scen == "3) Main Effects Only" ~ "3) Main Effects Only (|S_0| = 5)")) %>%
  mutate(scen = factor(scen, levels = c("1a) Strong Hierarchy (|S_0| = 7)",
                                        "1b) Weak Hierarchy (|S_0| = 5)",
                                        "1c) Interactions Only (|S_0| = 2)",
                                        "2) Linear Effects (|S_0| = 7)",
                                        "3) Main Effects Only (|S_0| = 5)")))

p1_mse_nactive <- ggplot(data = df_mse_nactive, aes(x = mean.nactive, y = mean.mse, color = method, label = method)) +
  geom_point(size = 2.1) +
  geom_text_repel(
    data = subset(df_mse_nactive, mean.nactive < 100),
    nudge_x      = 40,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_text_repel(
    data = subset(df_mse_nactive, mean.nactive >= 100),
    nudge_x      = 5,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_errorbar(aes(ymin = mean.mse - sd.mse, ymax = mean.mse + sd.mse), size = 1.1) +
  geom_errorbarh(aes(xmin = mean.nactive - sd.nactive, xmax = mean.nactive + sd.nactive), size = 1.1) +
  facet_rep_wrap(~scen, scales = "free", ncol = 2,
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  # scale_color_brewer(palette = "Dark2")+
  scale_color_manual(values=RColorBrewer::brewer.pal(12, "Paired")[-11], guide=guide_legend(ncol=3)) +
  labs(x="Number of active variables", y="Test Set MSE",
       title="Test Set MSE vs. Number of Active Variable (Mean +/- 1 SD)",
       subtitle="Based on 200 simulations",
       caption="") +
  theme_ipsum_rc(axis_title_just = "bt") +
  theme(legend.position = "right",
        legend.text=element_text(size=14),
        strip.text = element_text(size=14))

reposition_legend(p1_mse_nactive, 'center', panel='panel-2-3')

# pacman::p_load(psych)
# affect.mat2 <- describeBy(DT[, c("mse","nactive")], group = list(DT$method, DT$scen), mat = TRUE)
# dev.off()
# par(family="serif")
# error.crosses(affect.mat2[c(56:66),],
#               affect.mat2[c(1:11),],
#               labels=unique(affect.mat2$group1),
#               xlab="Number of Active Variables",
#               main = "ADNI Data: Means (+/- 1 SD) from 200 Train/Validate/Test Splits",
#               sd = TRUE,
#               cex.lab = 1.4,
#               cex.axis = 1.4,
#               cex.main = 1.5,
#               # xlim = c(0, 34),
#               ylab="Test Set MSE",
#               colors = RColorBrewer::brewer.pal(11, "Paired"),
#               pch=16,cex=2)

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

p1_nactive <- ggplot(DT, aes(method, nactive, fill = method)) +
  ggplot2::geom_boxplot() +
  facet_rep_wrap(~scen, scales = "free", ncol = 2,
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(12, "Paired")[-11], guide=guide_legend(ncol=3)) +
  labs(x="", y="Number of active variables",
       title="Number of active variables",
       subtitle="Based on 200 simulations",
       caption="") +
  theme_ipsum_rc() + theme(legend.position = "right", axis.text.x = element_text(angle = 25, hjust = 1),
                           legend.text=element_text(size=14),
                           strip.text = element_text(size=14))


reposition_legend(p1_nactive, 'center', panel='panel-2-3')




## ---- plot-tpr-fpr-sim ----

DT[fpr==1, fpr:=NA] # this is skewing the plots, so we remove them and set to NA
df_tpr_fpr <- DT[, c("method","scen","tpr","fpr")] %>%
  group_by(method, scen) %>%
  summarise(mean.tpr = mean(tpr, na.rm = TRUE), sd.tpr = sd(tpr, na.rm = TRUE),
            mean.fpr = mean(fpr, na.rm = TRUE), sd.fpr = sd(fpr, na.rm = TRUE)) %>%
  mutate(scen = case_when(scen == "1a) Strong Hierarchy" ~ "1a) Strong Hierarchy (|S_0| = 7)",
                          scen == "1b) Weak Hierarchy" ~ "1b) Weak Hierarchy (|S_0| = 5)",
                          scen == "1c) Interactions Only" ~ "1c) Interactions Only (|S_0| = 2)",
                          scen == "2) Linear Effects" ~ "2) Linear Effects (|S_0| = 7)",
                          scen == "3) Main Effects Only" ~ "3) Main Effects Only (|S_0| = 5)")) %>%
  mutate(scen = factor(scen, levels = c("1a) Strong Hierarchy (|S_0| = 7)",
                                        "1b) Weak Hierarchy (|S_0| = 5)",
                                        "1c) Interactions Only (|S_0| = 2)",
                                        "2) Linear Effects (|S_0| = 7)",
                                        "3) Main Effects Only (|S_0| = 5)")))

p1_tpr_fpr <- ggplot(data = df_tpr_fpr, aes(x = mean.fpr, y = mean.tpr, color = method, label = method)) +
  geom_point(size = 2.1) +
  # geom_text_repel() +
  geom_text_repel(
    data = subset(df_tpr_fpr, mean.fpr < 0.05),
    nudge_x      = 0.02,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_text_repel(
    data = subset(df_tpr_fpr, mean.fpr > 0.05),
    nudge_x      = 0.00,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2
  ) +
  geom_errorbar(aes(ymin = mean.tpr - sd.tpr, ymax = mean.tpr + sd.tpr), size = 1.1) +
  geom_errorbarh(aes(xmin = mean.fpr - sd.fpr, xmax = mean.fpr + sd.fpr), size = 1.1) +
  facet_rep_wrap(~scen, scales = "free", ncol = 2,
                 repeat.tick.labels = 'left',
                 labeller = as_labeller(appender,
                                        default = label_parsed)) +
  scale_color_manual(values=RColorBrewer::brewer.pal(12, "Paired")[-11], guide=guide_legend(ncol=3)) +
  labs(x="False positive rate", y="True positive rate",
       title="True Positive Rate vs. False Positive Rate (Mean +/- 1 SD)",
       subtitle="Based on 200 simulations",
       caption="") +
  theme_ipsum_rc(axis_title_just = "bt") +
  theme(legend.position = "right",
        legend.text=element_text(size=14),
        strip.text = element_text(size=14))

reposition_legend(p1_tpr_fpr, 'center', panel='panel-2-3')


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


