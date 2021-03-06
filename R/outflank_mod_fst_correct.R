#' @title Modificatin of OutFLANK's \code{WC_FST_FiniteSample_Diploids_2Alleles()}
#'
#' @description Skip the need for diploid genotypes, i.e. just use allele frequencies,
#' and produce the same output as \code{WC_FST_FiniteSample_Diploids_2Alleles()}.
#' Calculate FST with a sample size correction.
#'
#' @param alleleFreqs Numeric: A vector of the Reference allele frequency, with
#' the length == number of populations.
#'
#' @param sampSize Integer: A vector of sample sizes, with the length == number
#' of populations and the order corresponding to the those in \code{alleleFreqs}.
#'
#' @param Ho Integer: An optional vector of the heterozygote frequency,
#' with the length == number of populations and the order corresponding
#' to the those in \code{alleleFreqs}. Default = \code{NULL}, see Details.
#'
#' @details By using just allele frequencies, there is no information about
#' the observed frequency of heterozygotes, Ho. The function will make the
#' naive assumption that Ho = He, if Ho is not specified.
#'
#' @return Returns a list of values related to FST:
#'  \itemize{
#'  \item   He:  the expected heterozygosity of the locus
#'  \item   FST:  Fst (with sample size correction)
#'  \item 	T1: The numerator of the Fst calculation (a from Weir and Cockerham 1984)
#'  \item   T2: The denominator of the Fst calculation (a+b+c from Weir and Cockerham 1984)
#'  }

outflank_mod_fst_correct <-function(alleleFreqs, sampSize, Ho){

  # Determine the variance in allele frequencies
  sample_sizes = sampSize
  n_ave = mean(sample_sizes)
  n_pops = length(alleleFreqs) #r
  r = n_pops
  n_c = (n_pops*n_ave - sum(sample_sizes^2)/(n_pops*n_ave))/(n_pops-1)
  p_freqs = alleleFreqs
  p_ave = sum(sample_sizes*p_freqs)/(n_ave*n_pops)
  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)

  if(s2==0){return(1); break}

  # Heterozygosity
  He <- 1-sum(p_ave^2, (1-p_ave)^2)

  if(is.null(Ho)==FALSE){ h_freqs <- Ho
  } else{ h_freqs = He }

  h_ave = sum(sample_sizes*h_freqs)/(n_ave*n_pops)

  # This is where sample size correction occurs
  a <- n_ave/n_c*(s2 - 1/(n_ave-1)*(p_ave*(1-p_ave)-((r-1)/r)*s2-(1/4)*h_ave))

  b <- n_ave/(n_ave-1)*(p_ave*(1-p_ave) - (r-1)/r*s2 - (2*n_ave - 1)/(4*n_ave)*h_ave)

  c <- 1/2*h_ave

  FST <- a/(a+b+c)
  return(list(He=He,FST=FST, T1=a, T2=(a+b+c)))
}
