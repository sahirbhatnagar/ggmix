## ---- GAW20-prediction-RMSE-activeVariable ----

root <- "/home/sahir/git_repositories/ggmix/RealData/"
load(paste0(root, "GAW20-RMSE-Nactive.RData"))

ggplot(RMSEACTIVE, aes(x=meanACTIVE,y=meanRMSE,colour=Method,label=Method)) + geom_point(shape=18,size=5) +
  #ylim(0.29,0.38) +
  #xlim(0,40) +
  geom_errorbar(aes(ymin=meanRMSE - sdRMSE,ymax=meanRMSE + sdRMSE,width=0.5), size = 1.1) +
  geom_errorbarh(aes(xmin=meanACTIVE,xmax=upper90ACTIVE,height=0), size = 1.1) +
  geom_errorbarh(aes(xmin=meanACTIVE,xmax=upper95ACTIVE,height=0), size = 1.1, lty = 3) +
  scale_color_manual(values = cbbPalette[c(7,3,4,2)], guide = guide_legend(ncol=3)) +
  labs(x = "Number of active variables", y = "Root mean squared prediction error",
       title = "Prediction Root Mean Squared Error vs. Number of Active Variables",
       subtitle = "Based on five-fold cross validation of 200 GAW20 simulations",
       caption = "") + theme_bw() + geom_text(aes(label=Method),hjust=-0.5,vjust=-0.5) +
  theme(legend.position = "bottom",title = element_text(size = 20),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        strip.text = element_text(size = 18),
        axis.ticks=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black',size=0.3)) +
  scale_x_continuous(breaks = seq(1,51,5))




## ---- GAW20-R2 ----

root <- here::here("manuscript/data/")
load(paste0(root, "Chr1AroundCausal.ld.RData"))
Gaw <- as.matrix(Gaw)

ha <- rowAnnotation(foo = anno_mark(at = 309,labels = "rs9661059"))
har2 <- HeatmapAnnotation(`R-square` = anno_points(as.vector(Gaw[309,])))
ComplexHeatmap::Heatmap(Gaw,cluster_columns = F,cluster_rows = F,show_row_names = F,show_column_names = F,right_annotation = ha,top_annotation = har2,col = col_fun,name = "R squared")
