#' An example pooled allele frequency data set
#'
#' A long-format data table containing individual genotypes for 4 populations, each with 30
#' individuals, and 200 biallelic SNP loci. Simulated with FastSimCoal2 (Excoffier et al., 2013).
#'
#' @usage data(data_PoolFreqs)
#'
#' @format A data table with 800 rows and 11 columns.
#'
#' @details
#' These columns would be created by \code{genomalicious::poolne_estim_outputs}:
#' \itemize{
#'   \item \code{CHROM} = The contig ID.
#'   \item \code{POS} = The SNP position.
#'   \item \code{LOCUS} = The unique locus ID.
#'   \item \code{ALT} = The alternate allele.
#'   \item \code{REF} = The reference allele.
#'   \item \code{POOL} = The pool ID.
#'   \item \code{FREQ} = The frequency of the alternate allele.
#'   \item \code{DP} = The total read depth.
#'   \item \code{AO} = The alternate read counts.
#'   \item \code{RO} = The reference read counts.
#'   \item \code{INDS} = The number of pooled diploid individuals.
#' }
#'
#' @references Gautier et al. (2013) Estimation of population allele frequencies from
#' next-generation sequencing data: pool-versus individual-based genotyping.
#' Molecular Ecology, 22(14), 3766–3779.
#'
#' @docType data
#' @keywords datasets
#' @name data_PoolFreqs
NULL
