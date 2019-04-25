#' Calculate diversity statistics with \code{IDIP} from a data table of genotypes
#'
#' @param snpDat Data table: a long-format data table with SNP genotypes
#' coded as per VCF specifications, e.g. ('0/0', '0/1', '1/1').
#' Three columns are required:
#' \enumerate{
#'    \item (1) The sampled individual ID (see param \code{sampCol}).
#'    \item (2) The locus ID (see param \code{locusCol}).
#'    \item (3) The genotype (see param \code{genoCol}).
#' }
#'
#' @param sampCol Character: The column name with the sampled individual information.
#'
#' @param locusCol Character: The column name with the locus information.
#'
#' @param genoCol Character: The column name with the genotype information.
#'
#' @examples
#' data(genomalicious4pops)
#' snpDat <- genomalicious4pops
#' @export
idip_DT <- function(snpDat, strucMat
                    , sampCol=NA, locusCol=NA, genoCol=NA){

}
