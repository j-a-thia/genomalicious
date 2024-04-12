#' Retrieve amino acid names
#'
#' Retrieve amino acid names based on an input string, following IUPAC names.
#'
#' @param X.in Character: a vector of amino acids to retrieve names. Can be
#' either single letter, three letter, or full. See Details.
#'
#' @param X.type Character: a single value, one of "single", "three", or "full".
#' The type of the input for \code{X.in}.
#'
#' @param Y.type Character: a single value, one of "single", "three", or "full".
#' The type of the output you want to return for \code{X.in}.
#'
#' @param showTable Logical: Show the amino acid name table instead of performing
#' the name retrieval? Default = FALSE.
#'
#' @details Input the string of amino acid in either single letter (e.g., 'A'),
#' three letter (e.g., 'Ala') or full (e.g., 'Alanin') with \code{X.in}. Specify
#' what the input format is with \code{X.type}. Specify the desired format of
#' name you want back with \code{Y.type}.
#'
#' @returns Returns a character vector of amino acid names in desired output format.
#'
#' @examples
#' library(genoalicious)
#'
#' # Vector of amino acids in single letter form
#' aa_vec <- c('A','W','Y','N','K','Q')
#'
#' # Show the table
#' amino_acid_names(showTable=TRUE)
#'
#' # Match to three letters
#' aa3_vec <- amino_acid_names(X.in=aa_vec, X.type='single', Y.type='three')
#'
#' aa3_vec
#'
#' # Match three letters to full
#' aafull_vec <- amino_acid_names(X.in=aa3_vec, X.type='three', Y.type='full')
#'
#' aafull_vec
#'
#' @export

amino_acid_names <- function(X.in, X.type, Y.type, showTable=FALSE){

  require(data.table)

  aaTab <- data.table(
    single=c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'),
    three=c('Ala','Cys','Asp','Glu','Phe','Gly','His','Ile','Lys','Leu','Met','Asn','Pro','Gln','Arg','Ser','Thr','Val','Trp','Tyr'),
    full=c('Alanine','Cysteine','Aspartic acid','Glutamic acid','Phenylalanine','Glycine','Histidine','Isoleucine','Lysine','Leucine','Methionine','Asparagine','Proline','Glutamine','Arginine','Serine','Threonine','Valine','Tryptophan','Tyrosine')
  )

  type.vec <- colnames(aaTab)

  if(showTable==TRUE){
    return(aaTab)
  } else if(showTable==FALSE){
    # Check X.type is specified correctly.
    if(length(X.type)!=1){
      stop('Argument `X.type` must be a single character value. See ?amino_acid_names.')
    }
    if(!X.type %in% type.vec){
      stop('Argument `X.type` must be one of "single", "three", or "full". See ?amino_acid_names.')
    }
    # Check Y.type is specified correctly.
    if(length(Y.type)!=1){
      stop('Argument `Y.type` must be a single character value. See ?amino_acid_names.')
    }
    if(!Y.type %in% type.vec){
      stop('Argument `Y.type` must be one of "single", "three", or "full". See ?amino_acid_names.')
    }
    # Match X.in to X.type in the amino acid table and return as Y.type.
    return(aaTab[[Y.type]][match(X.in, aaTab[[X.type]])])
  }
}

