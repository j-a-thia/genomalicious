#' Generate a matrix of allele frequencies from genotypes
#'
#' Parses a data table of genotypes and returns a matrix of allele frequencies.
#'
#' @param dat Data table: It is expected that there are only two alleles, and therefore, only three possible genotypes:
#' 0/0, 0/1 (or 1/0), and 1/1, where the Ref allele is '0'. This data.table needs the following columns:
#' \code{POP}, the population ID; \code{SAMPLE}, the individual ID; \code{LOCUS}, the locus ID;
#' and \code{GT}, the genotype.
#'
#' @return Returns a matrix of allele frequencies for the Ref allele (coded as '0' in the genotype)
#'
#' @examples
#' # Import genotype data
#' data(genomaliciousGenos)
#' genomaliciousGenos
#'
#' # Convert to frequency matrix
#' genos2freqs(genomaliciousGenos)
#'
#' @export
genos2freqs <- function(dat){
  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  # Check the class of dat
  if(!'data.table' %in% class(dat)){ stop("Argument dat isn't a data table.")}

  # Check that demes has the right columns.
  if(length(which((c('POP', 'SAMPLE', 'LOCUS', 'GT') %in% colnames(dat))==FALSE)) > 0){
    stop("Argument dat needs the columns $POP, $SAMPLE, $LOCUS, and $GT.")
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Split the data based on LOCUS, then iterate through each LOCUS.
  popFreqs <- lapply(split(dat, dat$LOCUS), function(L){

    # An empty matrix that will store the population frequencies
    # for the Ref allele at the Lth LOCUS.
    locusFreqs <- matrix(NA, 0, 1)
    colnames(locusFreqs) <- L$LOCUS[1]

    # Iterate through each Pth population
    for(P in unique(L$POP)){
      # Subset the data by POP, obtain all alleles
      pop.als <- unlist(strsplit(L[POP==P]$GT, '/'))
      # What frequency are the Ref allele?
      pop.ref.freq <- length(which(pop.als=='0')) / length(pop.als)
      # Adjust the object's structure
      pop.ref.freq <- as.matrix(pop.ref.freq)
      colnames(pop.ref.freq) <- L$LOCUS[1]
      rownames(pop.ref.freq) <- P
      # Row bind the population data to the LOCUS matrix.
      locusFreqs <- rbind(locusFreqs, pop.ref.freq)
    }
    return(locusFreqs)
  })

  # Column bind the loci.
  return(do.call('cbind', popFreqs))

  # ............ END
}
