######################################
# R Source code file for creating simulated dataset to be included in the ggmix package
# Simulated data with admixture population using bnpsd package
# Author: Sahir Bhatnagar
# Created: June 20, 2018
# Updated:
#####################################


pacman::p_load_gh('StoreyLab/popkin')
pacman::p_load_gh('StoreyLab/bnpsd')
pacman::p_load(MASS)
devtools::load_all()

set.seed(2345)
admixed <- gen_structured_model(n = 500,
                                p_test = 200,
                                p_kinship = 1e4,
                                geography = "1d",
                                percent_causal = 0.05,
                                percent_overlap = "100",
                                k = 5, s = 0.5, Fst = 0.1,
                                b0 = 0, beta_mean = 1,
                                eta = 0.1, sigma2 = 1)

usethis::use_data(admixed, overwrite = TRUE)

# phi_eigen <- eigen(dat$kin)
# U_kinship <- phi_eigen$vectors
# Lambda <- phi_eigen$values
# which(Lambda < 1e-5)
# if (any(Lambda < 1e-5)) Lambda[which(Lambda < 1e-5)] <- 1e-05

# popkin::plotPopkin(dat$kin)
# res <- gic.penfam(x = dat$x, y = dat$y,  d = Lambda, u = U_kinship, an = log(length(dat$y)))
# dev.off()
# plot(res)
# res$penfam.fit$result
# (nonzero = res$penfam.fit$coef[,res$lambda.min.name,drop = F][nonzeroCoef(res$penfam.fit$coef[,res$lambda.min.name,drop = F]),,drop = F])
# (nonzero_names = setdiff(rownames(nonzero), c("beta0","eta","sigma2")))
# length(intersect(nonzero_names, dat$causal))/length(dat$causal)
# length(nonzero_names)
# res$penfam.fit$sigma2[,res$lambda.min.name]
# res$penfam.fit$eta[,res$lambda.min.name]


