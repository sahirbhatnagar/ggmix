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

## ---- Mice-comparison-fixTPR ----
root <- "/home/sahir/git_repositories/ggmix/RealData/"
load(paste0(root, "Mice-200Bootstrap.RData"))
load(paste0(root, "mice.RData"))


ggmixpie <- c(ggmixfail, 200 - ggmixfail)
ggmixlab <- c("Failure","Success")

lassopie <- c(lassofail, 200 - lassofail)
lassolab <- c("Failure","Success")

twosteppie <- c(twostepfail, 200 - twostepfail)
twosteplab <- c("Failure","Success")

#Twostep

layout(matrix(c(1,1:11,1,1,12:21,22,22:32,22,22,33:42,43,43:53,43,43,54:63), 6, 12, byrow = TRUE))
par(mar=c(2,2,2,2))

plotGenome <- genotype[,1:3]
plotGenome$count <- NA

pie(twosteppie,labels = paste0(prop.table(twosteppie)*100,"%"),cex=2,main = "(a) twostep",cex.main = 3,col = c("grey",cbbPalette[7]), radius = 0.8)

for (j in 1:nrow(plotGenome)) {
  plotGenome$count[j] <- length(which(twostepCoef==plotGenome$marker[j]))
}

for (chr in 1:20) {
  subdat <- plotGenome[plotGenome$chr==chr,]
  if (chr!=20) {
    if (!chr %in% c(1,11)) {
      plot(subdat$cM,subdat$count,type="l",xlab=NA,ylab=NA,ylim=c(0,200-twostepfail),main=paste0("Chr",chr),cex.main =2,yaxt='n',cex.axis=1.5)
      abline(h=(200-twostepfail)/2,lty=2,col="red")
    }
    else {
      plot(subdat$cM,subdat$count,type="l",xlab=NA,ylab=NA,ylim=c(0,200-twostepfail),main=paste0("Chr",chr),cex.main =2,cex.axis=1.5)
      abline(h=(200-twostepfail)/2,lty=2,col="red")
    }
  }
  else {
    plot(subdat$cM,subdat$count,type="l",xlab=NA,ylab=NA,ylim=c(0,200-twostepfail),main="ChrX",cex.main=2,yaxt='n',cex.axis=1.5)
    abline(h=(200-twostepfail)/2,lty=2,col="red")
  }
}


#LASSO

plotGenome <- genotype[,1:3]
plotGenome$count <- NA

pie(lassopie,labels = paste0(prop.table(lassopie)*100,"%"),cex=2,main = "(b) lasso",cex.main=3,col = c("grey",cbbPalette[3]), radius = 0.8)

for (j in 1:nrow(plotGenome)) {
  plotGenome$count[j] <- length(which(lassoCoef==plotGenome$marker[j]))
}

for (chr in 1:20) {
  subdat <- plotGenome[plotGenome$chr==chr,]
  if (chr!=20) {
    if (!chr %in% c(1,11)) {
      plot(subdat$cM,subdat$count,type="l",xlab=NA,ylab=NA,ylim=c(0,200-lassofail),main=paste0("Chr",chr),cex.main =2,yaxt='n',cex.axis=1.5)
      abline(h=(200-lassofail)/2,lty=2,col="red")
    }
    else {
      plot(subdat$cM,subdat$count,type="l",xlab=NA,ylab=NA,ylim=c(0,200-lassofail),main=paste0("Chr",chr),cex.main =2,cex.axis=1.5)
      abline(h=(200-lassofail)/2,lty=2,col="red")
    }
  }
  else {
    plot(subdat$cM,subdat$count,type="l",xlab=NA,ylab=NA,ylim=c(0,200-lassofail),main="ChrX",cex.main=2,yaxt='n',cex.axis=1.5)
    abline(h=(200-lassofail)/2,lty=2,col="red")
  }
}

#GGMIX

plotGenome <- genotype[,1:3]
plotGenome$count <- NA

pie(ggmixpie,labels = paste0(prop.table(ggmixpie)*100,"%"),cex=2,main = "(c) ggmix",cex.main=3,col = c("grey",cbbPalette[4]), radius = 0.8)

for (j in 1:nrow(plotGenome)) {
  plotGenome$count[j] <- length(which(ggmixCoef==plotGenome$marker[j]))
}

for (chr in 1:20) {
  subdat <- plotGenome[plotGenome$chr==chr,]
  if (chr!=20) {
    if (!chr %in% c(1,11)) {
      plot(subdat$cM,subdat$count,type="l",xlab=NA,ylab=NA,ylim=c(0,200-ggmixfail),main=paste0("Chr",chr),cex.main =2,yaxt='n',cex.axis=1.5)
      abline(h=(200-ggmixfail)/2,lty=2,col="red")
    }
    else {
      plot(subdat$cM,subdat$count,type="l",xlab=NA,ylab=NA,ylim=c(0,200-ggmixfail),main=paste0("Chr",chr),cex.main =2,cex.axis=1.5)
      abline(h=(200-ggmixfail)/2,lty=2,col="red")
    }
  }
  else {
    plot(subdat$cM,subdat$count,type="l",xlab=NA,ylab=NA,ylim=c(0,200-ggmixfail),main="ChrX",cex.main=2,yaxt='n',cex.axis=1.5)
    abline(h=(200-ggmixfail)/2,lty=2,col="red")
  }
}


## ---- UKB-Figure ----
root <- "/home/sahir/git_repositories/ggmix/RealData/"
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
# ukb_df$freq %>% sum
# ukb_df$genotyped %>% sum

ukb_long <- do.call(rbind, lapply(1:nrow(ukb_df), function(i) {
  data.frame(Chromosome = rep(ukb_df$Chromosome[i], each = ukb_df$freq[i]),
             Genotyped = c(rep("Genotyped", ukb_df$genotyped[i]), rep("Imputed", ukb_df$freq[i] - ukb_df$genotyped[i]))
  )
}
))

# ukb_long$Chromosome <- as.numeric(ukb_long$Chromosome)
# str(ukb_long)
# levels(ukb_long$Genotyped)
# ukb_long$Genotyped %>% table
# ukb_long$Chromosome %>% table
# table(ukb_long$Chromosome, ukb_long$Genotyped) %>% rowSums()

# str(ukb_long)

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
