#' Convert between genotype scores (separated alleles vs counts)
#'
#' Assumes biallelic genotypes and can interchange between separated
#' alleles ('0/0', '0/1', '1/1') or Alt allele counts (0, 1, 2).
#'
#' @param dat Character/Integer: A vector of genotypes. If the class
#' is \code{'character'}, then assumes separated alleles and converts into
#' allele counts. The opposite is true, if class is \code{'integer'}, will
#' convert into separated alleles.
#'
#' @details If \code{dat} is a vector of characters, then missing values should
#' take the form of './.'. Otherwise, if \code{dat} is a vector of integers,
#' missing values should take the form of \code{NA}.
#'
#' @examples
#' library(genomalicious)
#' genoscore_converter(c('0/0', '0/1', '1/1'))
#' genoscore_converter(c(0, 1, 2))
#'
#' @export
genoscore_converter <- function(dat){
  if(!class(dat) %in% c('character', 'integer', 'numeric')){
    stop('Class of argument dat is wrong: see ?genoscore_converter')
  }

  if(class(dat)=='numeric'){ dat <- as.integer(dat)}

  if(class(dat)=='character' & sum(dat %in% c('./.', '0/0', '0/1', '1/1'))!=length(dat)){
    stop('Character genotypes are not coded properly: see ?genoscore_converter')
  }

  if(class(dat)=='integer' & sum(dat %in% c(NA, 0, 1, 2))!=length(dat)){
    stop('Integer genotypes are not coded properly: see ?genoscore_converter')
  }

  if(class(dat)=='character'){
    gt <- lapply(dat, function(x){
          if(x=='./.'){ return(NA)
          } else{
            x <- unlist(strsplit(x, '/'))
            return(sum(as.integer(x)))
          }
    })
    return(unlist(gt))
  }

  if(class(dat)=='integer'){
    gt <- lapply(dat, function(x){
              if(is.na(x)){ return('./.')
              } else if(x==0){ return('0/0')
              } else if(x==1){ return('0/1')
              } else if(x==2){ return('1/1')
              }
            })
    return(unlist(gt))
  }
}
