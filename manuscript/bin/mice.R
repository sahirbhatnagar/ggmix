## ---- Mice-comparison-fixTPR ----
root <- here::here("manuscript/data/")
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



## ---- Mice-R2 ----

root <- here::here("manuscript/data/")
load(paste0(root,"markerCorrelationBCG.RData"))
col_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))

chrname <- substr(rownames(corMice),1,4)
chrcount <- c(63,43,33,49,31,34,29,32,32,27,41,28,27,23,29,21,26,17,18,16)
chrlabel <- c()
for (i in 1:20) {
  chrlabel <- c(chrlabel,rep(i,chrcount[i]))
}
chrlabel[chrlabel==20] <- "X"

ha <- rowAnnotation(foo = anno_mark(at = c(which(rownames(corMice)=="D1Mit435"),which(rownames(corMice)=="D11Mit119")),labels=c("D1Mit435","D11Mit119")))
hachr <- columnAnnotation(Chr = chrlabel)

ComplexHeatmap::Heatmap(corMice,cluster_columns = F,cluster_rows = F,show_row_names = F,show_column_names = F,right_annotation = ha,top_annotation = hachr,col = col_fun,name = "R squared")
