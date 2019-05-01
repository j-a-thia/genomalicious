#' Calculate diversity statistics with \code{IDIP} from a data table of genotypes
#'
#' @param snpDat Data table: a long-format data table with SNP genotypes
#' coded as per VCF specifications, e.g. ('0/0', '0/1', '1/1').
#' Three columns are required:
#' \enumerate{
#'    \item (1) The population ID (see param \code{popCol}).
#'    \item (2) The locus ID (see param \code{locusCol}).
#'    \item (3) The genotype (see param \code{genoCol}).
#' }
#'
#' @param strucMat Matrix: Hierarchy of "aggregates", minimum of 2 rows and
#' as many columns as the lowest units of the hierarchy. 1st row is the top most
#' level (e.g. ecosystem or metapopulation), with lower levels in proceding rows.
#'
#' @param popCol Character: The column name with the population information, i.e.
#' the lowest levels of the population hierarchy (last row of \code{strucMat}).
#'
#' @param locusCol Character: The column name with the locus information.
#'
#' @param genoCol Character: The column name with the genotype information.
#'
#' @return A matrix with columns representing loci, rows are the statistics
#' from the \code{IDIP} analysis.
#'
#' @examples
#' data(genomalicious_4pops)
#' snpDat <- genomalicious_4pops
#' strucMat <- matrix(c(rep('metapop', 4)
#'                , paste('Group', c(1,1,2,2))
#'                , paste('Pop', 1:4, sep='')), ncol=4, byrow=TRUE)
#' @export
idip_DTgenos <- function(snpDat, strucMat, popCol='POP'
                    , locusCol='LOCUS', genoCol='GT'){

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(libs in c('HierDpart', 'data.table')){ require(libs, character.only=TRUE) }

  # Stop if strucMat doesn't meet min rows
  if(nrow(struMat)<2){
    stop('Argument strucMat must be a matrix with at least 2 rows.')}

  # Population check
  uniqPops <- unique(snpDat$POP)
  if(sum(uniqPops %in% strucMat[nrow(strucMat),]) < length(uniqPops)){
    stop('All unique populations in snpDat must be in the last row (lowest level)
    in the strucMat matrix of population aggregate hierarchy.')
  }

  # Rename columns
  colnames(snpDat)[
    which(colnames(snpDat)%in%c(popCol, locusCol, genoCol))] <- c(
      'POP', 'LOCUS', 'GT')

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  # Unique loci
  uniqLoci <- unique(snpDat$LOCUS)

  # Run IDIP for each locus
  lociDivpar <- lapply(uniqLoci, function(locus){
              # For each locus, get Ref and Alt allele counts per population
              als <- lapply(uniqPops, function(pop){
                     allele_counts(snpDat[LOCUS==locus & POP==pop]$GT)
                     })
              als <- do.call('rbind', als)
              rownames(als) <- uniqPops
              # Run IDIP
              idip <- IDIP(abun=t(as.matrix(als)), struc=strucMat)
              return(idip)
  })
  lociDivpar <- do.call('cbind', lociDivpar)
  colnames(lociDivpar) <- uniqLoci

  return(lociDivpar)
}
