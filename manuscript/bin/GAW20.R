## ---- GAW20-prediction-RMSE-activeVariable ----

root <- here::here("manuscript/data/")
load(paste0(root, "GAW20-RMSE-Nactive.RData"))

knitr::kable(RMSEACTIVE)


## ---- GAW20-R2 ----

root <- here::here("manuscript/data/")
load(paste0(root, "Chr1AroundCausal.ld.RData"))
Gaw <- as.matrix(Gaw)
col_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))


ha <- rowAnnotation(foo = anno_mark(at = 309, labels = "rs9661059"))
har2 <- HeatmapAnnotation(`R-square` = anno_points(as.vector(Gaw[309, ])))
ComplexHeatmap::Heatmap(Gaw,
  cluster_columns = F, cluster_rows = F,
  show_row_names = F, show_column_names = F,
  right_annotation = ha, top_annotation = har2,
  col = col_fun, name = "R squared"
)
