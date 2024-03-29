#' Perform a PCoA (principal coordinates analysis) on population allele frequencies
#'
#' Takes a long-format data table of allele frequencies, calculates pairwise FST,
#' and conducts a PCoA of these pairwise FST values using R's \code{ape::pcoa()}
#' function (be sure to cite Paradis & Schliep 2019).
#'
#' @param dat Data table: A long-format data table containing allele
#' frequency estimates for biallelic SNP loci. Requires the following columns:
#' \enumerate{
#'   \item The population ID (see param \code{popCol}).
#'   \item The locus ID (see param \code{locusCol}).
#'   \item The allele frequency column column (see param \code{genoCol}).
#'   \item The number of pooled diploids used to generate the allele frequency
#'   estimates (see param \code{indsCol})
#' }
#'
#' @param popCol Character: The column name with the population information.
#' Default is \code{'POP'}.
#'
#' @param locusCol Character: The column name with the locus information.
#' Default is \code{'LOCUS'}.
#'
#' @param freqCol Character: The column name with the allele freuqency information.
#' Default is \code{'FREQ'}.
#'
#' @param indsCol Character: The column name with the number of individuals
#' contributing to the allele freuqency estimate. Default is \code{indsCol}.
#'
#' @returns Returns a list-like object of class \code{eigen}, as per \code{ape::pcoa},
#' with an additional long-format data table of genome-wide pairiwse FST
#' estimates in an index \code{$fst}.
#'
#' @references Paradis & Schliep (2019) Bioinformatics. DOI: 10.1093/bioinformatics/bty633
#'
#' @examples
#' library(genomalicious)
#' data(data_PoolFreqs)
#' data(data_PoolInfo)
#'
#' # Note columns in data_PoolFreqs
#' colnames(data_PoolFreqs)
#'
#' # We need to add in the number of diploid individuals, $INDS
#' newFreqData <- left_join(data_PoolFreqs, data_PoolInfo)
#' head(newFreqData)
#'
#' # Fit the PCoA
#' PCOA <-pcoa_freqs(newFreqData)
#'
#' # Object class
#' class(PCOA)
#'
#' # Axes
#' PCOA$vectors
#'
#' # Eigenvalues
#' PCOA$values
#'
#' # Pairwise FST in long-format data table
#' PCOA$fst
#'
#' @export
pcoa_freqs <- function(
    dat, popCol='POP', locusCol='LOCUS', freqCol='FREQ', indsCol='INDS'){

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  require(ape); require(data.table); require(tidyverse)

  if(sum(c(popCol,locusCol,freqCol,indsCol) %in% colnames(dat))!=4){
    stop(
      'Argument `dat` requires columns `popCol`, `locusCol`, `freqCol`, and
        `indsCol` to estimate F-statistics from allele frequencies. See ?fstat_calc.')
  }

  dat <- as.data.table(dat) %>%
    copy %>%
    setnames(., old=c(popCol,locusCol,freqCol,indsCol), c('POP','LOCUS','FREQ','INDS'))

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  pair.fst <- dat %>%
    fstat_calc(., type='freqs', pairwise=TRUE, global=FALSE)

  pair.mat <- pair.fst$genome %>%
    pairwiseMat2DT(., flip=TRUE, X1='POP1', X2='POP2', Y='FST')

  PCOA <- pcoa(pair.mat)

  PCOA$fst <- pair.fst[['genome']]

  return(PCOA)
}

