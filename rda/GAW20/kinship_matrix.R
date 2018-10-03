

#load("~/GAW manuscript/Cov_confounder_714_samples.RData")  # 714 subjects have pre-treatment methylation data
load("~/git_repositories/ggmix/rda/GAW20/Cov_confounder_714_samples.RData")  # 714 subjects have pre-treatment methylation data

#----------------------------------------------------------------------------------#
# Covariate of primary interest: difference of HDL and Trig between visit 4 and 2
#----------------------------------------------------------------------------------#

mycov = data.frame(hdl.diff =Cov.confder$hdl.post - Cov.confder$hdl.pre, trr.diff= log(Cov.confder$trr.post) - log(Cov.confder$trr.pre))
rownames(mycov) = Cov.confder$GAWSUBJ

#--------------------------------------------------------------------------------#
# Confounders
#--------------------------------------------------------------------------------#

myconf = Cov.confder[,c(15, 2:10)]
myconf = myconf[,-c(7,9)] # delete the fast time 1 and 3; because methylation and lipid measures (we used) were taken at visit 2 and 4.
rownames(myconf) = Cov.confder$GAWSUBJ


#--------------------------------------------------------------------------------#
# Kinship matrix
#--------------------------------------------------------------------------------#
library(kinship2)

PED  = read.csv (file = "~/share/share/PROJECTS/GAW/GAW20/Real Data Package/real_phenos_covars/PED.csv.gz")


PED[is.na(PED$sex),"sex"] = rep(3, sum(is.na(PED$sex)))
rel = PED[!is.na(PED$zygosity),c("GAWSUBJ", "GAWTWIN", "zygosity")]

colnames(rel) =c("id1", "id2", "code")
gaw.ped=pedigree( id = PED$GAWSUBJ, dadid=PED$GAWDAD, momid=PED$GAWMOM,
                  sex=  PED$sex, affected = rep(1, nrow(PED)), status=rep(1, nrow(PED)), relation  = rel)
gaw.kin=kinship(gaw.ped)

# kinship for the 714 subjects
sample.used <- Cov.confder$GAWSUBJ
id.use <- which(colnames(gaw.kin) %in% sample.used)

kin.714 <- gaw.kin[id.use, id.use]


