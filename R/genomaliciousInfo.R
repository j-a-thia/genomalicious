#' An example of pool-seq sample info
#'
#' A dataset containing the population pool ID and the number of individuals
#' used in a pool-seq experiment. \cr
#' \cr
#' Data of this structure would be used as the \code{dat} parameter in the
#' function \code{poolne_estim_input}, which prepares pooled replicate sequencing
#' data for analysis by Gautier et al.'s (2013) \code{poolne_estim}.
#'
#' @usage data(pgposerInfo)
#'
#' @format A data table with 2 columns and 2 rows.
#'
#' @details
#' \itemize{
#'   \item \code{POOL} = The population pool ID.
#'   \item \code{INDS} = The number of diploid individuals that contributed toward
#'     the population pool.
#' }
#'
#' @references Gautier et al. (2013) Estimation of population allele frequencies from
#' next-generation sequencing data: pool-versus individual-based genotyping.
#' Molecular Ecology, 22(14), 3766–3779.
#'
#' @docType data
#' @keywords datasets
#' @name pgposerInfo
NULL
