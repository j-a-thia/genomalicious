#' An example of pool-seq read info
#'
#' A dataset containing read count information from a pool-seq experiment with
#' replicate library preparations. \cr
#' \cr
#' Data of this structure would be used as the \code{dat} parameter in the
#' function \code{poolne_estim_input}, which prepares pooled replicate sequencing
#' data for analysis by Gautier et al.'s (2013) \code{poolne_estim}.
#'
#' @usage data(pgposerReads)
#'
#' @format A data table with 6 columns and 18 rows.
#'
#' @details
#' \itemize{
#'   \item \code{POOL} = The population pool ID.
#'   \item \code{SAMPLE} = The sample replicate ID. Has values of 1-3, which
#'     represent one of three replicate sequencing samples of the same pool of individuals.
#'   \item \code{CHROM} = The chromosome ID.
#'   \item \code{LOCUS} = The locus ID.
#'   \item \code{RO} = The number of Ref read counts at the locus.
#'   \item \code{DP} = The depth (total reads observed) at the locus.
#' }
#'
#' @references Gautier et al. (2013) Estimation of population allele frequencies from
#' next-generation sequencing data: pool-versus individual-based genotyping.
#' Molecular Ecology, 22(14), 3766–3779.
#'
#' @docType data
#' @keywords datasets
#' @name pgposerReads
NULL
