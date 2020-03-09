#' Karim's Simulated Data
#'
#' A simulated dataset with a kinship matrix
#'
#' @format A list with 6 elements: \describe{ \item{b}{vector of length 1000
#'   representing the true regression coefficients. 10 non-zero coefficients,
#'   the rest are 0.}\item{kin1}{the true kinship matrix}\item{s.g}{polygenic
#'   variance, set to be 1.26}\item{s.e}{error variance, set to be
#'   1}\item{h.tot}{the total trait heritability. Set to be 60%.}\item{G}{matrix
#'   of genotypes of dimension 600 x 1000 SNPs, with approximately 800 common
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
#' A simulated dataset to show the utility of this package
#'
#' @details The code used to simulate the data is available at
#'   \url{https://github.com/sahirbhatnagar/ggmix/blob/master/data-raw/bnpsd-data.R}.
#'    See \code{\link{gen_structured_model}} for more details on the output and
#'   how the function used to simulate the data.
#' @return A list with the following elements \describe{\item{ytrain}{simulated
#'   response vector for training set} \item{ytune}{simulated response vector
#'   for tuning parameter selection set} \item{ytest}{simulated response vector
#'   for test set} \item{xtrain}{simulated design matrix for training
#'   set}\item{xtune}{simulated design matrix for tuning parameter selection
#'   set}\item{xtest}{simulated design matrix for testing set}
#'   \item{xtrain_lasso}{simulated design matrix for training set for lasso
#'   model. This is the same as xtrain, but also includes the nPC principal
#'   components} \item{xtune_lasso}{simulated design matrix for tuning parameter
#'   selection set for lasso model. This is the same as xtune, but also includes
#'   the nPC principal components}\item{xtest}{simulated design matrix for
#'   testing set for lasso model. This is the same as xtest, but also includes
#'   the nPC principal components} \item{causal}{character vector of the names
#'   of the causal SNPs} \item{beta}{the vector of true regression coefficients}
#'   \item{kin_train}{2 times the estimated kinship for the training set
#'   individuals} \item{kin_tune_train}{The covariance matrix between the tuning
#'   set and the training set individuals} \item{kin_test_train}{The covariance
#'   matrix between the test set and training set individuals}
#'   \item{Xkinship}{the matrix of SNPs used to estimate the kinship matrix}
#'   \item{not_causal}{character vector of the non-causal SNPs} \item{PC}{the
#'   principal components for population structure adjustment} }
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
