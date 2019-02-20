#' Calculate reference allele frequencies form read counts
#'
#' Parses a data table of read counts and calculates Ref allele frequencies \cr
#' as the observed reference reads divided by the total number of reads. This is \cr
#' primarily meant for pool-seq data.
#'
#' @param dat Data table: The read counts. Requires the following columns: \cr
#' \enumerate{
#'    \item \code{$POOL} = the population pool ID.
#'    \item \code{$LOCUS} = the locus ID.
#'    \item \code{$DP} = the total read depth.
#'    \item \code{$RO} = the observed reference read counts.
#' }
#'
#' @return A matrix of allele freuqencies for each population. Loci in columns, populations
#' in rows, and Ref allele frequencies in the cells.
#'
#' @examples
#' data(pgposerReads)
#'
#' # Pull out a single replicate from the pgposer demo dataset.
#' rep1Reads <- pgposerReads[grep(pattern='Rep1', x=pgposerReads$SAMPLE),]
#' rep1Reads
#'
#' # Add a population pool column:
#' rep1Reads$POOL <- unlist(lapply(strsplit(rep1Reads$SAMPLE, '_'), function(X){ X[1] }))
#'
#' # Generate a matrix of Ref allele frequencies
#' freqMat <- reads2freqs(rep1Reads)
#' freqMat
#'
#' @export
reads2freqs <- function(dat){
  # BEGIN ............

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(libs in c('data.table', 'tidyr')){ require(libs, character.only=TRUE)}

  # Check the class of dat.
  if(!'data.table' %in% class(dat)){ stop("Argument dat isn't a data table.")}

  # Check that the correct columns are in dat.
  if(length(which((c('POOL', 'LOCUS', 'DP', 'RO') %in% colnames(dat))==FALSE)) > 0){
    stop("Argument dat needs the columns $POOL, $LOCUS, $DP, and $RO.")
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Calculate Ref frequencies as RO/DP
  dat[, RF:=RO/DP]

  # Convert data from long to wide format
  freqDT <- tidyr::spread(dat[, c('POOL', 'LOCUS', 'RF')], key=LOCUS, value=RF)

  # Turn freqDT from data.table to a matrix, dropping the POOL column
  freqMat <- as.matrix(freqDT[, !'POOL', with=FALSE])

  # Give rownames the population POOL IDs
  rownames(freqMat) <- freqDT$POOL

  # Return the matrix
  return(freqMat)

  # ............ END
}
