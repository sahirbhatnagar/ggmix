pacman::p_load_gh('StoreyLab/popkin')
pacman::p_load_gh('StoreyLab/bnpsd') # load this after popkin
pacman::p_load('lfa')
pacman::p_load('RColorBrewer')

# this shows that UKB doesnt have much relatedness ------------------------

X <- BEDMatrix::BEDMatrix("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/penfam/simulation/data/bed/X_other_8k")
Phi <- popkin(X)
# set outer margin for axis labels (left and right are non-zero)
par(oma=c(0,1.5,0,3))
# set inner margin for subpopulation labels (bottom and left are non-zero), add padding
par(mar=c(1,1,0,0)+0.2)
# now plot!
plotPopkin(Phi)

par(oma=c(0,1.5,0,3))
par(mar=c(1,1,0,0)+0.2)
plotPopkin(inbrDiag(Phi))

Phi[1:5,1:5]

# using bnspd package to simulate data ------------------------------------


# dimensions of data/model
n <- 100 # number of individuals (NOTE this is 10x less than in publication!)
k <- 10 # number of intermediate subpops
# define population structure
FF <- 1:k # subpopulation FST vector, up to a scalar
s <- 0.5 # desired bias coefficient
Fst <- 0.1 # desired FST for the admixed individuals
obj <- q1d(n, k, s=s, F=FF, Fst=Fst) # admixture proportions from 1D geography
Q <- obj$Q
FF <- obj$F
# get pop structure parameters of the admixed individuals
Theta <- coanc(Q,FF) # the coancestry matrix
diag(Theta)

pacman::p_functions("bnpsd")
# verify that we got the desired Fst!
Fst2 <- fst(Q,FF)
Fst2
# this should also equal Fst
inbr <- diag(Theta)
popkin::fst(inbr)
s2 <- mean(Theta)/Fst
s2
# visualize the per-subpopulation inbreeding coefficients (FSTs)
par(mar=c(4,4,0,0)+0.2) # shrink default margins
colIS <- brewer.pal(k, "Paired") # indep. subpop. colors
barplot(FF, col=colIS, names.arg=1:k, ylim=c(0,1),
        xlab='Subpopulation', ylab='Inbreeding coeff.')
dev.off()
par(mar=c(1,4,0,0)+0.2) # shrink default margins
barplot(t(Q), col=colIS, border=NA, space=0, ylab='Admixture prop.')
mtext('Individuals', 1)

dev.off()
# Visualize the coancestry matrix using "popkin"!
# set outer margin for axis labels (left and right are non-zero)
par(oma=c(0,1.5,0,3))
# zero inner margin (plus padding) because we have no individual or subpopulation labels
par(mar=c(0,0,0,0)+0.2)
# now plot!
popkin::plotPopkin(Theta)

diag(Theta)

m <- 10000 # number of loci in simulation (NOTE this is 30x less than in publication!)
# draw all random Allele Freqs (AFs) and genotypes
# reuse the previous F,Q
out <- rbnpsd(Q, FF, m)
X <- out$X # genotypes
P <- out$P # IAFs (individual-specific AFs)
B <- out$B # intermediate AFs
pAnc <- out$Pa # ancestral AFs
par(mar=c(4,4,0,0)+0.2) # shrink default margins for these figures
# inspect distribution of ancestral AFs (~ Uniform(0.01,0.5))
hist(pAnc, xlab='Ancestral AF', main='', xlim=c(0,1))

# distribution of intermediate population AFs
# (all subpopulations combined)
# (will be more dispersed toward 0 and 1 than ancestral AFs)
hist(B, xlab="Intermediate Subpopulation AF", main="", xlim=c(0,1))

# distribution of IAFs (admixed individuals)
# (admixture reduces differentiation, so these resemble ancestral AFs a bit more)
hist(P, xlab="Individual-specific AF", main="", xlim=c(0,1))

dev.off()
# genotype distribution of admixed individuals
barplot(table(X), xlab="Genotypes", ylab="Frequency", col="white")


# for best estimates, group individuals into subpopulations using the geography
# this averages more individuals in estimating the minimum kinship
subpops <- ceiling( (1:n)/n*k )
table(subpops) # got k=10 subpops with 100 individuals each

# now estimate kinship using popkin
PhiHat <- popkin(X, subpops)
PhiHat2 <- popkin(X)
diag(PhiHat)
diag(PhiHat)*2
head(PhiHat)
# replace diagonal with inbreeding coeffs. to match coancestry matrix
ThetaHat <- inbrDiag(PhiHat)
diag(ThetaHat)

# Visualize the coancestry matrix using "popkin"!
# set outer margin for axis labels (left and right are non-zero)
par(oma=c(0,1.5,0,3))
# increase inner top margin for panel titles
par(mar=c(0,0,2.5,0)+0.2)
# now plot!
x <- list(Theta, ThetaHat)
x <- list(PhiHat, 2*PhiHat2)
titles <- c("Truth", "Estimate")
plotPopkin(x, titles)

pacman::p_funs("BEDMatrix")
tt <- gaston::as.bed.matrix(t(X))
kin <- gaston::GRM(tt, autosome.only = FALSE)
kin[1:5,1:5]
2*PhiHat[1:5,1:5]
2*PhiHat2[1:5,1:5]

dev.off()


par(oma=c(0,1.5,0,3))
# increase inner top margin for panel titles
par(mar=c(0,0,2.5,0)+0.2)
plotPopkin(list(kin, 2*PhiHat, 2*PhiHat2), c("gaston kin","2*popkin subpops","2*popkin"))


all(complete.cases(kin))
eiK <- eigen(kin)
# all(rownames(as.matrix(x))==rownames(kin))
# deal with a small negative eigen value
any(eiK$values < 0)
sum(eiK$values < 0)
eiK$values[ eiK$values < 0 ] <- 0
PC <- sweep(eiK$vectors, 2, sqrt(eiK$values), "*")
plot(eiK$values)
plot(PC[,2],PC[,5])
