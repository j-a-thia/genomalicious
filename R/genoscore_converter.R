#' Convert between genotype scores (separated alleles vs counts)
#'
#' Assumes biallelic genotypes and can interchange between separated
#' alleles ('0/0', '0/1', '1/1') and Ref allele counts (0, 1, 2).
#'
#' @param dat Character/Integer: A vector of genotypes. If the class
#' is \code{'character'}, then assumes separated alleles and converts into
#' allele counts. The opposite is true if class is \code{'integer'}, will
#' convert into separated alleles.
#'
#' @export
genoscore_converter <- function(dat){
  if(!class(dat) %in% c('character', 'integer')){
    stop('Class of argument dat is wrong.')
  }

  if(class(dat)=='character'){
    gt <- lapply(strsplit(dat, split='/', fixed=TRUE), function(x){sum(as.integer(x))})
    return(unlist(gt))
  }

  if(class(dat)=='integer'){
    gt <- lapply(dat, function(x){
              if(x==2){ return('0/0')
              } else if(x==1){ return('0/1')
              } else if(x==0){ return('1/1')
                    }
            })
    return(unlist(gt))
  }
}
