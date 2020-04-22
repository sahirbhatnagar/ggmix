## ---- UKB-Figure ----
root <- here::here("manuscript/data/")
load(paste0(root, "UKB-plot-data.RData"))
p1 <- ggplot(valstats, aes(x=index,y=RMSE,color=method)) +
  geom_line(size = 1.5) +
  theme_bw() +
  scale_color_manual(values = cbbPalette[c(7,3,4)]) +
  theme(legend.position = "bottom", axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15)) +
  labs(color="",x=expression(lambda~index),y="RMSE in model selection set") +
  ylim(0.8,1.3)

p2 <- ggplot(teststats, aes(x=log2(coef),y=RMSE,color=method)) + ylim(0.9,1) + xlim(9,12) +
  geom_point(size = 5) + xlab(expression(log[2](Number~of~active~variables))) +
  ylab("RMSE in test set") +
  labs(color="") +
  theme_bw()  +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(angle = 0,hjust=0.5,vjust = 1)) +
  scale_color_manual(values = cbbPalette[c(7,3,4,2)])

p1+p2+plot_annotation(tag_levels = 'A', tag_prefix = "(", tag_suffix = ")")


## ---- UKB-chromosome-distribution ----

ukb_df <- data.frame(Chromosome = factor(1:22),
                     freq = c(493,989,754,645,522,1515,520,388,338,190,141,658,204,300,543,138,292,231,128,962,3,46),
                     genotyped = c(65,95,62,78,55,284,48,42,47,20,24,71,24,14,77,27,56,20,20,99,1,4))


ukb_long <- do.call(rbind, lapply(1:nrow(ukb_df), function(i) {
  data.frame(Chromosome = rep(ukb_df$Chromosome[i], each = ukb_df$freq[i]),
             Genotyped = c(rep("Genotyped", ukb_df$genotyped[i]), rep("Imputed", ukb_df$freq[i] - ukb_df$genotyped[i]))
  )
}
))

ggplot(ukb_long, aes(x=Chromosome, fill = Genotyped)) +
  geom_bar() +
  ylab("Number of SNPs used in model") +
  labs(color="") +
  theme_bw()  +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(angle = 0,hjust=0.5,vjust=0)) +
  scale_fill_manual(values = cbbPalette[c(7,3)])
