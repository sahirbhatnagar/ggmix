## ---- packages ----

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

# pacman::p_load_gh("sahirbhatnagar/sail", dependencies = FALSE)
devtools::load_all("/home/sahir/git_repositories/sail/")
# devtools::load_all("C:/Users/sahir/Documents/git_repositories/sail")
pacman::p_load(ggplot2)
pacman::p_load(doMC)
registerDoMC(cores = 8)
pacman::p_load(latex2exp)
pacman::p_load(truncnorm)
pacman::p_load(lemon)
pacman::p_load(magrittr)
pacman::p_load(simulator)
pacman::p_load(dplyr)
pacman::p_load(tidyr)
pacman::p_load(data.table)
pacman::p_load(psych)
# pacman::p_load(cowplot)
pacman::p_load_gh("hrbrmstr/hrbrthemes")
pacman::p_load(ggrepel)
pacman::p_load(Cairo)
pacman::p_load(extrafont)
extrafont::loadfonts()
# pacman::p_load(multipanelfigure)

## ---- globals ----

cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
trop <- RSkittleBrewer::RSkittleBrewer("trop")
gg_sy <- theme(legend.position = "bottom", axis.text = element_text(size = 20),
               axis.title = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 20))
