#' Generate a matrix of allele frequencies from genotypes
#'
#' Parses a data table of genotypes and returns a matrix of allele frequencies.
#'
#' @param dat Data table: It is expected that there are only two alleles, and therefore, only three possible genotypes:
#' 0/0, 0/1 (or 1/0), and 1/1, where the Ref allele is '0'. This data.table needs the following columns:
#' \enumerate{
#'    \item The population ID (see param \code{popCol}).
#'    \item The sample ID (see param \code{sampCol})
#'    \item The locus ID (see param \code{locusCol}).
#'    \item The genotypes (see param \code{genoCol}).
#' }
#' @param popCol Character: The column name with the sampled individual information. Default = \code{'POP'}.
#'
#' @param sampCol Character: The column name with the sampled individual information. Default = \code{'SAMPLE'}.
#'
#' @param locusCol Character: The column name with the locus information. Default = \code{'LOCUS'}.
#'
#' @param genoCol Character: The column name with the genotype information. Default = \code{'GT'}.
#'
#' @param returnMat Logical: Should an allele frequency matrix be returned?
#' Default = \code{TRUE}, but if \code{FALSE}, will return a data table.
#'
#' @return Returns a matrix of allele frequencies for the Alt allele, or
#' a data table with the columns \code{$POP}, \code{$SAMPLE}, \code{$LOCUS}
#' and \code{$ALT.FREQ}.
#'
#' @examples
#' # Import genotype data
#' data(data_4pops)
#' data_4pops
#'
#' # Convert to frequency matrix
#' genos2freqs(data_4pops)
#'
#' @export
genos2freqs <- function(dat
                        , popCol='POP'
                        , sampCol='SAMPLE'
                        , locusCol='LOCUS'
                        , genoCol='GT'
                        , returnMat=TRUE){
  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  # Check the class of dat
  if(!'data.table' %in% class(dat)){ stop("Argument dat isn't a data table.")}

  # Check that all columns are specified correctly
  if(sum(c(popCol, sampCol, locusCol, genoCol) %in% colnames(dat)) != 4){
    stop('Arguments popCol, sampCol, locusCol, and genoCol must be
         columns in dat. See ?genos2freqs.')
  }

  # Reassign column names
  colReass <- match(c(popCol, sampCol, locusCol, genoCol), colnames(dat))
  colnames(dat)[colReass] <- c('POP', 'SAMPLE', 'LOCUS', 'GT')

  # Convert to integer counts
  genoClass <- class(dat$GT)
  if(genoClass=='character'){
    dat[, GT:=genoscore_converter(GT)]
  } else if(genoClass=='numeric'){
    dat[, GT:=as.integer(GT)]
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Get the Alt allele frequency
  freqDt <- dat[, sum(GT)/(length(GT)*2), by=c('POP', 'LOCUS')]
  setnames(freqDt, 'V1', 'ALT.FREQ')

  # If a matrix of allele frequencies desired?
  if(returnMat==TRUE){
    freqMat <- as.matrix(spread(freqDt, key=LOCUS, value=ALT.FREQ)
                  , rownames='POP')
    return(freqMat)
  } else{ return(freqDt) }

  # ............ END
}
