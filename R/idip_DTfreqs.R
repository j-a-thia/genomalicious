#' Calculate diversity statistics with \code{IDIP} from a data table of frequencies
#'
#' @param freqDat Data table: a long-format data table with SNP genotypes
#' allele frequencies. Three columns are required:
#' \enumerate{
#'    \item (1) The population ID (see param \code{popCol}).
#'    \item (2) The locus ID (see param \code{locusCol}).
#'    \item (3) The allele frequencies (see param \code{freqCol}).
#' }
#'
#' @param strucMat Matrix: Hierarchy of "aggregates", minimum of 2 rows and
#' as many columns as the lowest units of the hierarchy. 1st row is the top most
#' level (e.g. ecosystem or metapopulation), with lower levels in proceding rows.
#'
#' @param popCol Character: The column name with the population information, i.e.
#' the lowest levels of the population hierarchy (last row of \code{strucMat}).
#' Default = \code{'POP'}.
#'
#' @param locusCol Character: The column name with the locus information.
#' Default = \code{'LOCUS'}.
#'
#' @param freqCol Character: The column name with the allele frequency information.
#' Default = \code{'FREQ'}.
#'
#' @param alleleCol Character: If number of alleles > 2, then specify a column that
#' contains the allele information, otherwise, Default = \code{NA}. Ideally,
#' all alleles possible at a locus should recorded for each population.
#'
#'
#'
#' @return A matrix with columns representing loci, rows are the statistics
#' from the \code{IDIP} analysis.
#'
#' @examples
#' data(genomalicious_Freqs)
#' freqDat <- genomalicious_PoolPi
#' strucMat <- matrix(c(rep('metapop', 4)
#'                , paste('Group', c(1,1,2,2))
#'                , paste('Pop', 1:4, sep='')), ncol=4, byrow=TRUE)
#' idip_DTfreqs(freqDat=genomaliciousPi, strucMat=metaPop
#'              , popCol='POOL', locusCol='LOCUS', freqCol='PI')
#'
#' @export
idip_DTfreqs <- function(freqDat, strucMat, popCol='POP'
                         , locusCol='LOCUS', freqCol='FREQ'
                         , alleleCol=NA){

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(libs in c('HierDpart', 'data.table', 'tidyr')){ require(libs, character.only=TRUE) }

  # Stop if strucMat doesn't meet min rows
  if(nrow(struMat)<2){
    stop('Argument strucMat must be a matrix with at least 2 rows.')}

  # Population check
  uniqPops <- unique(freqDat$POP)
  if(sum(uniqPops %in% strucMat[nrow(strucMat),]) < length(uniqPops)){
    stop('All unique populations in freqDat must be in the last row (lowest level)
         in the strucMat matrix of population aggregate hierarchy.')
  }

  # Rename columns
  colnames(freqDat)[
    match(c(popCol, locusCol, freqCol), colnames(freqDat))] <- c(
      'POP', 'LOCUS', 'FREQ')

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  # Unique loci
  uniqLoci <- unique(freqDat$LOCUS)

  # Run IDIP for each locus
  lociDivpar <- lapply(uniqLoci, function(locus){
                sub <- freqDat[LOCUS==locus]
                sub[, ALT.FREQ:=1-FREQ]
                sub <- spread(sub[, c('POP', 'LOCUS', 'FREQ')]
                              , key='POP', value='FREQ')[, !'LOCUS']
                sub <- rbind(sub, 1-sub[1,])

    # Run IDIP
    idip <- IDIP(abun=as.matrix(sub), struc=strucMat)
    return(idip)
  })
  lociDivpar <- do.call('cbind', lociDivpar)
  colnames(lociDivpar) <- uniqLoci

  return(lociDivpar)
  }
