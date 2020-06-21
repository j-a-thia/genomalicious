#' Filter loci by minor allele frequency (MAF)
#'
#' Parses a matrix of allele frequencies to determine which loci conform to the
#' desired minor allele frequency.
#'
#' @param dat Matrix or data.table: Default expectation is that user is supplying a matrix of allele
#' frequencies; i.e., loci in columns, populations in rows, and allele frequencies in cells (i.e.
#' \code{type=='freqs'}). \n\n
#' Alternatively, #' a data table of genotypes (i.e. \code{type=='genos'}). Loci are
#' expected to be biallelic. Recorded either as '/' separated alleles (0/0, 0/1 1/1),
#' or as counts of the Alt allele (0, 1, 2). For a data.table, expecting the columns:
#' \enumerate{
#'    \item The sample ID (see param \code{sampCol})
#'    \item The locus ID (see param \code{locusCol}).
#'    \item The genotypes (see param \code{genoCol}).
#' }
#'
#' @param maf Numeric: The minor allele frequency. E.g. 0.05 will filter for 5%, which will remove
#' a locus if its frequency is < 0.05 or > 0.95.
#'
#' @param type Character: Default = 'freqs', expected that \code{dat} is a matrix of allele frequencies.
#' Alternatively, if \code{dat} is a data table of of genotypes, set \code{type} to 'genos'.
#'
#' @param popCol Character: The column name with the sampled population information. Default = \code{'POP'}.
#' Only needed when \code{dat} is a long-format data table of genotypes.
#'
#' @param sampCol Character: The column name with the sampled individual information. Default = \code{'SAMPLE'}.
#' Only needed when \code{dat} is a long-format data table of genotypes.
#'
#' @param locusCol Character: The column name with the locus information. Default = \code{'LOCUS'}.
#' Only needed when \code{dat} is a long-format data table of genotypes.
#'
#' @param genoCol Character: The column name with the genotype information. Default = \code{'GT'}.
#' Only needed when \code{dat} is a long-format data table of genotypes.
#'
#' @return Returns an integer vector of column numbers in \code{dat} that conform
#' to the MAF value specified. These values can then be used to filter the allele frequency matrix.
#'
#' @examples
#' # MATRIX OF ALLELE FREQUENCIES
#' data(data_Freqs)
#'
#' # Filter for MAF=0.20
#' filter_maf(data_Freqs, maf=0.20, type='freqs')
#'
#' # LONG TABLE OF GENOTYPES
#' data(data_4pops)
#'
#' # Filter for MAF=0.05
#' filter_maf(data_4pops, maf=0.05, type='genos')
#'
#'
#' @export
filter_maf <- function(dat, maf=0.05, type='freqs', sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT'){

  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  libs <- library(data.table)
  for(L in libs){ require(L, character.only=TRUE)}

  # Make sure the class of dat matches the type of data specified.
  if('matrix'%in%class(dat) & type!='freqs'){
    stop("Argument dat is a matrix. Check this is a matrix of allele frequencies
         and set argument type to 'freqs'. See ?filter_maf.")
  }
  if ("data.table" %in% class(dat) & type!="genos") {
    stop("Argument dat is a data table. Check this is a data.table of genotypes
         and set argument type to 'genos'. See ?filter_maf.")
  }

  # Check columns for data.table
  if('data.table' %in% class(dat)){
    if(sum(sampCol, locusCol, genoCol) %in% colnames(dat)!=3){
      stop("Not all specified columns (sampCol, locusCol, genoCol) are in data.table dat. See ?filter_maf.")
    }
  }

  # Check that the MAF is between 0 and 1.
  if(maf < 0 | maf > 1){
    stop("Argument maf needs to be a numeric between 0 and 1. See ?filter_maf.")
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Specify the min and max MAF
  minF <- maf
  maxF <- 1 - maf

  # If the input is a matrix of allele frequencies (columns = loci, rows = pops)
  if(type=='freqs'){
    test <- apply(dat, 2, function(f){
      if(min(f) >= minF & max(f) <= maxF){ return('Yes')
      } else {
        return('No')
      }
    })
    return(names(which(test=='Yes')))
  }

  # If the input if a data.table of individuals and genotypes.
  if(type=='genos'){
    # Reassign column names
    colReass <- match(c(sampCol, locusCol, genoCol), colnames(dat))
    colnames(dat)[colReass] <- c('SAMPLE', 'LOCUS', 'GT')

    # Frequencies
    freqs <- genos2freqs(dat)

    # Test for MAF
    test <- apply(freqs, 2, function(f){
      if(min(f) >= minF & max(f) <= maxF){ return('Yes')
      } else {
        return('No')
      }
    })
    return(names(which(test=='Yes')))
  }


}
