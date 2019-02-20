#' An example of collated data after analysis by \code{poolne_estim}
#'
#' \code{pgposer} allows collation of output files from Gautier et al.'s (2013) programme
#' \code{poolne_estim}. \cr
#' \cr
#' Read counts from replicate sequencing events are processed by \code{PoolNeEstim_gen_inputs},
#' which provides the inputs for \code{poolne_estim}. The results of the \code{poolne_estim}
#' analysis can then be pulled into \code{R} with the function \code{poolne_estim_outputs}. \cr
#' \cr
#' The product of this workflow is exemplified in this dataset. \cr
#' \cr
#' Data of this structure forms one of the key inputs for a number of analysis functions
#' provided by \code{pgposer}, e.g.: \cr
#' \itemize{
#'   \item \code{sim_freq}
#'   \item \code{pca_pools}
#'   \item \code{dapc_pools}
#'   \item \code{Fst_pools_pairs}
#'   \item \code{Fst_pools_demes}
#' }
#'
#' @usage data(pgposerPi)
#'
#' @format A data table with 6 columns and 15 rows.
#'
#' @details
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
#'
#' @references Gautier et al. (2013) Estimation of population allele frequencies from
#' next-generation sequencing data: pool-versus individual-based genotyping.
#' Molecular Ecology, 22(14), 3766–3779.
#'
#' @docType data
#' @keywords datasets
#' @name pgposerPi
NULL
