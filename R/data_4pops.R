#' An example of multi-population genotype data
#'
#' A dataset containing individual genotypes for 4 populations, each with 30
#' individuals, and 1205 biallelic SNP loci. Simulated with FastSimCoal2 (Excoffier et al., 2013).
#'
#' @usage data(data_4pops)
#'
#' @format A data table with 144600 rows and 5 columns.
#'
#' @details
#' \itemize{
#'   \item \code{POP} = The population ID.
#'   \item \code{SAMPLE} = The sample ID: Ind[Pop]_[Number].
#'   \item \code{CHROM} = The contig ID.
#'   \item \code{LOCUS} = The unique locus ID.
#'   \item \code{GT} = The genotype ID, allele separated by '/'. 0 = Ref allele, 1 = Alt allele.
#' }
#'
#' @references Excoffier et al. (2013) Robust demographic inference from genomic SNP data.
#' PLoS Genetics: 9(10), e1003905.
#'
#' @docType data
#' @keywords datasets
#' @name data_4pops
NULL
