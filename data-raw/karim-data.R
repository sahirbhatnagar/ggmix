######################################
# R Source code file for creating simulated dataset to be included in the ggmix package
# Simulated data from Karim
# Author: Sahir Bhatnagar
# Created: June 13, 2018
# Updated:
#####################################

#' Ci-joint les données dont je t’ai parlé cette après-midi. C’est une liste
#' “data” avec
#'
#' data$b c’est les 1000 betas avec 10 sans non-nuls; kin1: matrice kinship 600
#' x 600;
#' s.e la variance de l’erreur
#' s.g la polygenic variance
#' G: matrice 600 x
#' 1000 SNPs, avec environ 800 communs et 200 rares;
#'
#' If you simulate data using this scenario:
#'
#' Phi <- 2 * kin1
#' P = mvrnorm(1, rep(0,600), s.g*Phi)
#' y <- intercept + G%*%b+ P + rnorm(600,0,s.e)
#'
#' then the QTL heritability of y will be 8% (i.e. the 10 SNPs will explain 8%
#' of the trait’s total heritability), and the trait total heritability is set
#' to be 60%.

load(file = "data-raw/Scenario0.RData")
karim <- data
usethis::use_data(karim, overwrite = TRUE)
