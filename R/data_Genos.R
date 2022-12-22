#' An example genotype data set
#'
#' A long-format data table containing individual genotypes for 4 populations, each with 30
#' individuals, and 200 biallelic SNP loci. Simulated with FastSimCoal2 (Excoffier et al., 2013).
#'
#' @usage data(data_Genos)
#'
#' @format A data table with 24000 rows and 11 columns.
#'
#' @details
#' \itemize{
#'   \item \code{CHROM} = The contig ID.
#'   \item \code{POS} = The SNP position.
#'   \item \code{LOCUS} = The unique locus ID.
#'   \item \code{POP} = The population ID.
#'   \item \code{SAMPLE} = The sample ID: Ind[Pop]_[Number].
#'   \item \code{GT} = The genotype ID, coded as counts of the alternate allele.
#'   \item \code{DP} = The total read depth.
#'   \item \code{AO} = The alternate read counts.
#'   \item \code{RO} = The reference read counts.
#'   \item \code{ALT} = The alternate allele.
#'   \item \code{REF} = The reference allele.
#' }
#'
#' @references Excoffier et al. (2013) Robust demographic inference from genomic SNP data.
#' PLoS Genetics: 9(10), e1003905.
#'
#' @docType data
#' @keywords datasets
#' @name data_Genos
NULL
