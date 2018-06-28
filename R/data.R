#' Karim's Simulated Data
#'
#' A simulated dataset with a kinship matrix
#'
#' @format A list with 6 elements: \describe{ \item{b}{vector of length 1000
#'   representing the true regression coefficients. 10 non-zero coefficients,
#'   the rest are 0.}\item{kin1}{the true kinship matrix}\item{s.g}{polygenic
#'   variance, set to be 1.26}\item{s.e}{error variance, set to be
#'   1}\item{h.tot}{the total trait heritability. Set to be 60%.}\item{G}{matrix
#'   of genotypes of dimension 600 x 1000 SNPs, with approximately 800 commun
#'   and 200 rare SNPs}}
#' @details If you simulate data using the scenario provided in the example,
#'   then the QTL heritability of y will be 8% (i.e. the 10 SNPs will explain 8%
#'   of the trait’s total heritability), and the trait total heritability is set
#'   to be 60%
#' @examples
#' data(karim)
#' # Simulate a response using the genotype matrix and the kinship matrix
#' Phi <- 2 * karim$kin1
#' intercept <- 1
#' P <- MASS::mvrnorm(1, rep(0,600), karim$s.g * Phi)
#' y <- intercept + karim$G %*% karim$b + P + rnorm(600,0,karim$s.e)
"karim"


#' Simulated Dataset with 1D Geography
#'
#' A simualted dataset to show the utility of this package
#'
#' @details The code used to simulate the data is available at
#'   \url{https://github.com/sahirbhatnagar/ggmix/blob/master/data-raw/bnpsd-data.R}.
#'    See \code{\link{gen_structured_model}} for more details on the output and
#'   how the function used to simulate the data.
#' @format A list with the following elements \describe{\item{y}{simulated
#'   response vector} \item{x}{simulated design matrix, 500 x 200 matrix}
#'   \item{causal}{character vector of the names of the causal SNPs}
#'   \item{beta}{the vector of true regression coefficients. 5% of the SNPs in x
#'   are causal} \item{kin}{2 times the estimated kinship} \item{Xkinship}{the
#'   matrix of SNPs used to estimate the kinship matrix, 500 x 10,000 matirx}
#'   \item{not_causal}{character vector of the non-causal SNPs}
#'   \item{causal_positive}{character vector of the causal SNPs with positive
#'   regression coefficient} \item{causal_negative}{character vector of the
#'   causal SNPs with negative regression coefficient}\item{x_lasso}{the design
#'   matrix which also includes \code{nPC} principal components} }
#' @references Ochoa, Alejandro, and John D. Storey. 2016a. "FST And Kinship for
#'   Arbitrary Population Structures I: Generalized Definitions." bioRxiv
#'   doi:10.1101/083915.
#' @references Ochoa, Alejandro, and John D. Storey. 2016b. "FST And Kinship for
#'   Arbitrary Population Structures II: Method of Moments Estimators." bioRxiv
#'   doi:10.1101/083923.
#' @examples
#' data(admixed)
#' str(admixed)
"admixed"
