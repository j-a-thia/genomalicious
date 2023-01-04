#' Subfunction for calculating F-statistics
#'
#' Calculates F-statistics from genotype or allele frequency data in a
#' long-format data table as per Weir & Cockerham (1984). Assumes that all
#' values are for a single biallelic SNP locus.
#'
#' Note, this function is not exported.
#'
#' @param dat Data table: In long-format, requires columns \code{$POP},
#' \code{$SAMPLE}, \code{$LOCUS} and \code{$GT} for genotypes;
#' requires \code{$POP}, \code{$LOCUS}, \code{$FREQ}, and \code{$INDS} for
#' allele frequencies
#'
#' @param type Character: either \code{"genos"} or \code{"freqs"}.
#'
#' @param fstat Character: a vector containing the F-statistics to calculate
#' when the data are genotypes. Must include one of \code{"FST"},
#' \code{"FIS"}, or \code{"FIT"}.
#'
#' @returns Returns a list with \code{$genome}, the genome-wide F-statistic, and
#' \code{$locus}, the per locus F-statistic.
#'
#' @examples
#' data_Genos %>%
#'    fstat_subfun(., type='genos', fstat=c('FST','FIS','FIT'))
#'
#' left_join(data_PoolFreqs, data_PoolInfo) %>%
#'    setnames(., 'POOL', 'POP') %>%
#'    fstat_subfun(., type='freqs')

fstat_subfun <- function(dat, type, fstat=NULL){
  require(data.table); require(tidyverse)

  # Empty list
  result <- list(genome=NULL, locus=NULL)

  # If genotype data
  if(type=='genos'){
    # Summarise
    D <- list(
      dat[, .(P=sum(GT)/(length(GT)*2)), by=c('LOCUS','POP')],
      dat[, .(N=length(unique(SAMPLE))), by=c('LOCUS','POP')],
      dat[, .(H=sum(GT==1)/length(GT)), by=c('LOCUS','POP')]
    ) %>%
      reduce(full_join, by=c('LOCUS','POP'))

    num.pops <- D$POP %>% unique %>% length

    # Variance components
    D.varcomp <- D[, fstat_varcomps(pi=P, ni=N, r=num.pops, hi=H, hetStand=TRUE), by=LOCUS]

    # Iterate through requested F-statistics
    for(f in fstat){
      if(f=='FST'){
        result$genome[[f]] <- D.varcomp[, sum(A)/sum(A+B+C)]
        result$locus[[f]] <- D.varcomp[, .(FST=A/(A+B+C)), by=LOCUS]
      } else if(f=='FIS'){
        result$genome[[f]] <- D.varcomp[, 1-(sum(C)/sum(B+C))]
        result$locus[[f]] <- D.varcomp[, .(FIS=1-(C/sum(B+C))), by=LOCUS]
      } else if(f=='FIT'){
        result$genome[[f]] <- D.varcomp[, 1-(sum(C)/sum(A+B+C))]
        result$locus[[f]] <- D.varcomp[, .(FIT=1-(C/sum(A+B+C))), by=LOCUS]
      }
    }

    # Output
    result <- list(
      genome=as.data.table(do.call('cbind', result$genome)),
      locus=result$locus %>% reduce(full_join, by='LOCUS')
    )
  }

  # If allele frequency data
  if(type=='freqs'){
    # Number of populations
    num.pops <- dat$POP %>% unique %>% length

    # Variance components
    D.varcomp <- dat[, fstat_varcomps(pi=FREQ, ni=INDS, r=num.pops, hetStand=FALSE), by=LOCUS]

    # Calculations
    result$genome <- D.varcomp[, .(FST=sum(MSP-MSG)/sum(MSP+((Nc-1)*MSG)))]
    result$locus <- D.varcomp[, .(FST=(MSP-MSG)/(MSP+((Nc-1)*MSG))), by=LOCUS]
  }

  # Output
  return(result)
}
