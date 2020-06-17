#' An example of collating results from \code{poolne_estim}, in R
#'
#' \code{genomalicious} allows collation of output files from Gautier et al.'s (2013) program
#' \code{poolne_estim}. \cr
#' \cr
#' Read counts from replicate sequencing events are processed by \code{genomalicious::poolne_estim_inputs},
#' which provides the inputs for \code{poolne_estim}. The results of the \code{poolne_estim}
#' analysis can then be imported into R with the function \code{genomalicious::poolne_estim_outputs}. \cr
#' \cr
#' The product of this workflow is exemplified in this dataset. \cr
#'
#' @usage data(data_PoolPi)
#'
#' @format A data table with 32 rows and 10 columns.
#'
#' @details
#' These columns would be created by \code{genomalicious::poolne_estim_outputs}:
#' \itemize{
#'   \item \code{MRK} = The marker ID, assigned by \code{poolne_estim}.
#'   \item \code{PI} = The posterior mean estimate of the population Ref allele frequency,
#'     as estimated by \code{poolne_estim}.
#'   \item \code{SD} = The standard deviation of the the population Ref allele frequency,
#'     as estimated by \code{poolne_estim}.
#'   \item \code{POOL} = The population pool ID.
#'   \item \code{CHROM} = The chromosome ID.
#'   \item \code{LOCUS} = The locus ID.
#' }
#' The following columns have been added as metadata:
#' #' \itemize{
#'   \item \code{REF} = The reference allele.
#'   \item \code{ALT} = The alternate allele.
#'   \item \code{INDS} = The number of diploid individuals pooled.
#' }
#'
#' @references Gautier et al. (2013) Estimation of population allele frequencies from
#' next-generation sequencing data: pool-versus individual-based genotyping.
#' Molecular Ecology, 22(14), 3766–3779.
#'
#' @docType data
#' @keywords datasets
#' @name data_PoolPi
NULL
