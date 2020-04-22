## ---- packages ----

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(ggmix)
pacman::p_load(ComplexHeatmap) # for heatmaps
pacman::p_load(circlize) # for heatmaps
pacman::p_load(ggplot2) # for UKB imputed plot
pacman::p_load(latex2exp)
# pacman::p_load(truncnorm)
# pacman::p_load(lemon)
pacman::p_load(qwraps2) # used for simulation results table
pacman::p_load(kableExtra) # used for simulation results table
pacman::p_load(magrittr)
pacman::p_load(tidyr)
pacman::p_load_gh("thomasp85/patchwork")
pacman::p_load(simulator)
pacman::p_load(dplyr)
pacman::p_load(data.table)
pacman::p_load(here)
# pacman::p_load(roxygen2)
pacman::p_load(popkin)
pacman::p_load(bnpsd)
# pacman::p_load_gh("hrbrmstr/hrbrthemes")
# pacman::p_load(ggrepel)
# pacman::p_load(Cairo)
# pacman::p_load(extrafont)
# extrafont::loadfonts()

## ---- globals ----

cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
gg_sy <- ggplot2::theme(legend.position = "bottom", axis.text = element_text(size = 20),
                        axis.title = element_text(size = 20), legend.text = element_text(size = 20),
                        legend.title = element_text(size = 20))


