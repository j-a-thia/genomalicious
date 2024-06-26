#' Allele frequency correlations across genomic windows
#'
#' This function takes a data table of population allele frequencies
#' across different chromosomes/contigs and estimates correlations between
#' loci within specified windows. Loci within these windows can be
#' randomly subsmapled to reduce computational time and memory cost.
#'
#' @param dat Data table: Contains populations allele frequencies.
#'
#' @param chromCol Character: The column with chromosome information.
#' Default is 'CHROM'.
#'
#' @param posCol Character: The column with position information.
#' Default is 'POS'.
#'
#' @param locusCol Character: The column with locus ID.
#' Default is 'LOCUS'.
#'
#' @param popCol Character: The column with population ID.
#' Default is 'POP'.
#'
#' @param freqCol Character: The column with frequencies of a focal allele.
#' Default is 'FREQ'.
#'
#' @param windowSize Integer: The window size to use. Default is 10000.
#'
#' @param numLoci Integer: The number of loci to subsample per window. You may
#' need to adjust this value depending on the density of markers in your dataset.
#' Default is 0, which triggers all loci to be used. Specify >0 if you want
#' to control the subsampling.
#'
#' @returns Returns a long-format data table of pairwise correlations between
#' loci within windows, within chromosomes. Has the columns:
#' \enumterate{
#'  \item \code{$CHROM}, the chromosome ID.
#'  \item \code{$WINDOW}, the window ID.
#'  \item \code{$LOCUS.1}, the locus 1 ID.
#'  \item \code{$CHROM}, the locus 2 ID.
#'  \item \code{$DIST}, the distance between loci in base pairs.
#'  \item \code{$CORR}, the Pearson correlation between locus allele frequencies.
#' }
#'
#' @export
#'
#' @examples
#' library(genomalicious)
#'
#' # RADseq pool-seq dataset
#' data(data_PoolFreqs)
#'
#' # Note, only subset of contigs have multiple alleles, so only
#' # a few pairwise loci are returned.
#' allele_freq_corr_windows(data_PoolFreqs)

allele_freq_corr_windows <- function(
    dat, chromCol='CHROM', posCol='POS', locusCol='LOCUS',
    popCol='POP', freqCol='FREQ', windowSize=10000, numLoci=0
){
  # -------------------------------------------+
  # CHECKS
  # -------------------------------------------+
  require(data.table); require(tidyverse)

  # Are columns present?
  check.cols <- c(chromCol, posCol, locusCol, popCol, freqCol)
  if(sum(check.cols %in% colnames(dat))!=5){
    stop(
      'Check that all columns specified in arguments `chromCol`, `posCol`, `locusCol`, `popCol` and `freqCol` are present in `dat`. See ?allele_freq_corr_windows.'
    )
  }

  # Reassign columns
  dat <- dat %>%
    copy %>%
    setnames(., check.cols, c('CHROM','POS','LOCUS','POP','FREQ'))

  # -------------------------------------------+
  # MAIN CODE
  # -------------------------------------------+
  uniq.chrom <- dat$CHROM %>% unique

  dist.cor.list <- list()

  # Iterate over chromosomes
  for(chrom in uniq.chrom){
    # Subset chromosome
    dat.chrom <- dat[CHROM==chrom]

    # Max sequence position
    seq.max <- max(dat.chrom$POS)

    # Table of windows
    win.tab <- data.table(START=seq(1, seq.max, windowSize)) %>%
      .[, END:=START+windowSize] %>%
      .[, WINDOW:=1:.N]

    # Table of unique loci
    loc.tab <- dat.chrom[, c('CHROM','POS','LOCUS')] %>%
      unique()

    # Iterate over windows
    for(win in 1:max(win.tab$WINDOW)){
      # Get the unique loci in window
      uniq.loc.win <- loc.tab %>%
        .[POS>=win.tab$START[win] & POS<=win.tab$END[win]] %>%
        .[['LOCUS']]

      # If there are loci in the window, get correlations
      if(length(uniq.loc.win)>0){
        # Loci to subset?
        if(numLoci==0){
          work.loc.win <- uniq.loc.win
        } else if(numLoci < length(uniq.loc.win)){
          work.loc.win <- sample(uniq.loc.win, size=numLoci, replace=FALSE)
        } else if(numLoci > length(uniq.loc.win)){
          work.loc.win <- uniq.loc.win
        }

        # Loci distances
        loc.dist.win <- dat.chrom[LOCUS %in% work.loc.win] %>%
          .[,c('LOCUS','POS')] %>%
          unique %>%
          as.data.frame() %>%
          column_to_rownames(., 'LOCUS') %>%
          dist(., method='euclidean') %>%
          as.matrix() %>%
          genomalicious::pairwiseMat2DT(., X1='LOCUS.1', X2='LOCUS.2', Y='DIST') %>%
          data.table(CHROM=chrom, WINDOW=win, .)

        # Loci correlations
        loc.cor.win <- dat.chrom[LOCUS %in% work.loc.win] %>%
          DT2Mat_freqs(., popCol='POP', locusCol='LOCUS', freqCol='FREQ') %>%
          cor(., method='pearson') %>%
          genomalicious::pairwiseMat2DT(., X1='LOCUS.1', X2='LOCUS.2', Y='CORR') %>%
          data.table(CHROM=chrom, WINDOW=win, .)

        # Save
        dist.cor.list[[chrom]][[win]] <- left_join(
          loc.dist.win, loc.cor.win,
          by=c('CHROM','WINDOW','LOCUS.1','LOCUS.2')
        ) %>%
          as.data.table() %>%
          .[LOCUS.1!=LOCUS.2]
      }
    }
  }

  lapply(dist.cor.list, function(X) do.call('rbind', X)) %>%
    do.call('rbind', .) %>%
    return()
}
