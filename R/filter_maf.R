#' Filter loci by minor allele frequency (MAF)
#'
#' Parses a matrix of allele frequencies to determine which loci conform to the
#' desired minor allele frequency.
#'
#' @param dat Matrix or data table: Default expectation is that user is supplying a matrix of Ref allele
#' frequencies; i.e., loci in columns, populations in rows, and allele frequencies in cells. Alternatively,
#' a data table of genotypes can be supplied and allele frequencies will be calculated. It is expected that
#' there are only two alleles, and therefore, only three possible genotypes: 0/0, 0/1 (or 1/0), and 1/1, where
#' the Ref allele is '0'. This data table needs the following columns: \code{$POP}, the population ID; \code{$IND},
#' the individual ID; \code{$LOCUS}, the locus ID; and \code{$GENO}, the genotype.
#'
#' @param maf Numeric: The minor allele frequency. E.g. 0.05 will filter for 5%, which will remove
#' a locus if its frequency is < 0.05 or > 0.95.
#'
#' @param type Character: Default = 'freqs', expected that \code{dat} is a matrix of allele frequencies.
#' Alternatively, if \code{dat} is a data table of of genotypes, set \code{type} to 'genos'.
#'
#' @return Returns an integer vector of column numbers in \code{dat} that conform
#' to the MAF value specified. These values can then be used to filter the allele frequency matrix.
#'
#' @examples
#'
#'
#' @export
filter_maf <- function(dat, maf=0.05, type='freqs'){

  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  libs <- library(data.table)
  for(L in libs){ require(L, character.only=TRUE)}

  # Make sure the class of dat matches the type of data specified.
  if('matrix'%in%class(dat) & type!='freqs'){
    stop("Argument dat is a matrix. Check this is a matrix of allele frequencies
         and set argument type to 'freqs'.")
  }
  if ("data.table" %in% class(dat) & type!="genos") {
    stop("Argument dat is a data table. Check this is a data table of genotypes\n
         and set argument type to 'genos'.")
  }

  # Check that the MAF is between 0 and 1.
  if(maf < 0 | maf > 1){
    stop("Argument maf needs to be a numeric between 0 and 1.")
  }

  # Check that all the correct columns are in dat.
  if(class(dat)=='data.table'){
    if(length(which((c('POP', 'IND', 'LOCUS', 'GENO') %in% colnames(dat))==FALSE)) > 0){
      stop("Argument dat needs the columns $POP, $IND, $LOCUS, and $GENO.")
    }
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

  # If the input if a data.table of indiviudals and genotypes.
  # Needs a column $POP (pop info), $IND (individual ID), $LOCUS
  if(type=='genos'){
    F <- genos2freqs(dat)

    test <- apply(F, 2, function(f){
      if(min(f) >= minF & max(f) <= maxF){ return('Yes')
      } else {
        return('No')
      }
    })
    return(names(which(test=='Yes')))
  }


}
