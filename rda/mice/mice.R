setwd("/mnt/GREENWOOD_SCRATCH/tianyuan.lu/BCG/")

library(ggmix)
library(gaston)
library(glmnet)

load("mice.RData")
load("~/Dropbox/mcgill/students/tianyuan/results/mice/mice.RData")
STRAND <- rs42

markerinfo <- genotype[,c(1:3)]

pureGeno <- genotype[,-c(1:3)]
pureGeno <- t(pureGeno)

response <- STRAND$load
sampCount <- as.matrix(table(STRAND$rcs))

kinGeno <- as.data.frame(matrix(NA,ncol = nrow(genotype)))
colnames(kinGeno) <- genotype$marker
for (i in 1:nrow(sampCount)) {
  nrep <- sampCount[i,1]
  nstrain <- as.numeric(rownames(sampCount)[i])
  for (j in 1:nrep) {
    kinGeno <- rbind.data.frame(kinGeno,pureGeno[nstrain,])
  }
}
kinGeno <- kinGeno[-1,]
rownames(kinGeno) <- rownames(STRAND)

predictor <- as.matrix(cbind.data.frame(kinGeno,STRAND$sex))
kinGeno <- as.matrix(kinGeno)

DT <- as.bed.matrix(kinGeno)
standardize(DT) <- "p"
kin <- GRM(DT, autosome.only = FALSE)

ggmixFit <- ggmix(x=predictor,y=response,K=kinGeno,verbose = 1)
#ggmixFit <- ggmix(x=predictor,y=response,kinship=kin,verbose = 1)
bicGGMIX <- gic(ggmixFit,an = log(length(admixed$y)))
coef(bicGGMIX,type = "nonzero")

x1 <- cbind(rep(1,nrow(predictor)))
fit_lme <- lmm.aireml(response, x1, kin)
gaston_resid <- response - (fit_lme$BLUP_omega + fit_lme$BLUP_beta)

fitglmnet <- glmnet(x = predictor, y = gaston_resid, standardize = T, alpha = 1, intercept = T)

tLL <- fitglmnet$nulldev - deviance(fitglmnet)
k <- fitglmnet$df
n <- fitglmnet$nobs
BICtwostep <- log(n)*k - tLL
BICtwostep

coefMat <- as.matrix(fitglmnet$beta)
leastBIC <- coefMat[,which.min(BICtwostep)]
nz_names_twostep <- leastBIC[leastBIC!=0]
nz_names_twostep

pcKin <- prcomp(kin, scale. = T, center = T)
topthreePCs <- pcKin$rotation[,c(1:3)]
rownames(topthreePCs) <- rownames(kin)

pdwithPC <- as.matrix(cbind.data.frame(topthreePCs,predictor))
fitlasso <- glmnet(x = pdwithPC, y = response,
                   standardize = T, intercept = T,
                   alpha=1, penalty.factor = c(rep(0,3), rep(1,ncol(pdwithPC)-3)))

tLL <- fitlasso$nulldev - deviance(fitlasso)
k <- fitlasso$df
n <- fitlasso$nobs
BIClasso <- log(n)*k - tLL
BIClasso

coefMat <- as.matrix(fitlasso$beta)
leastBIC <- coefMat[,which.min(BIClasso)]
nz_names_lasso <- leastBIC[leastBIC!=0]
nz_names_lasso
