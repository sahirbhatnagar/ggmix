context(strwrap("check full and low rank constructors have valid output using
        admixed dataset"))

data("admixed")

svdXkinship <- svd(admixed$Xkinship)
number_nonzero_eigenvalues <- 20

fullrank_kinship <- new_fullrank_kinship(kinship = admixed$kin)
fullrank_K <- new_fullrank_K(K = admixed$Xkinship)
fullrank_UD <- new_fullrank_UD(U = svdXkinship$u, D = svdXkinship$d)

# kin <- gaston::GRM(gaston::as.bed.matrix(admixed$Xkinship), autosome.only = FALSE)

lowrank_kinship <- new_lowrank_kinship(kinship = admixed$kin,
                                       n_nonzero_eigenvalues = 20,
                                       n_zero_eigenvalues = nrow(admixed$kin) - 20)
lowrank_K <- new_lowrank_K(K = admixed$Xkinship)
lowrank_UD <- new_lowrank_UD(U = svdXkinship$u, D = svdXkinship$d)

