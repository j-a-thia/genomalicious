#' Generate an allele frequency data table
#'
#' Takes a data.table of genotypes or allele counts and calculates the allele
#' frequency for each allele. Can be used for multiallelic datasets.
#'
#' @param dat Data.table: Long-format data table of variants, e.g., as read in
#' with \code{genomalicious::vcf2DT}.
#'
#' @param type Character: Two modes, one of "genos" for individual genotype data,
#' or "counts" of allele in populations.
#'
#' @param sampCol Character: The column with sample ID information. Default is "SAMPLE".
#' Only needed if \code{type=="genos"}.
#'
#' @param popCol Character: The column with population ID information. Default is "POP".
#'
#' @param locusCol Character: The column with locus ID information. Default is "LOCUS".
#'
#' @param genoCol Character: The column with genotype information. Default is "GT".
#' Only needed if \code{type=="genos"}. Genotypes must be in character format where
#' alleles are separated by the delimiter, "/". For example, "0/1" is one Ref and
#' one Alt allele 1; "2/2" is two Alt allele 2.
#'
#' @param countCol Character: The column with allele count information for all alleles.
#' For example, in pool-seq of populations, the number of read counts for each allele.
#' Default is "COUNTS". Only needed if \code{type=="counts"}. Counts should be separated
#' by commas, with the Ref allele first. E.g., "20,60,4" would indicate 20, 60, and 4
#' counts of the Ref allele, Alt allele 1, and Alt allele 2, respectively.
#'
#' @param indsCol Character: The column with the number of sampled individuals per
#' population. Default is "INDS".
#'
#' @details This function assumes no missing values. For \code{type=="genos"}, all
#' sampled individuals must have a genotype value for each locus.
#' For \code{type=="counts"}, all sampled populations must have count data for each
#' locus. You could impute for individuals, or drop loci with missing data for
#' for individual or population datasets.
#'
#' Note, when \code{type=="counts"}, the allele frequencies are based on the
#' proportion of counts per allele relative to the total number of observed counts
#' at a locus. However, this function will align the total sample number of
#' sequenced individuals against the counts.
#'
#' @returns Returns a long format data table with the following columns:
#' \enumerate{
#'  \item \code{$POP}, the population ID column.
#'  \item \code{$LOCUS}, the locus ID column.
#'  \item \code{$ALLELE}, the allele ID column (0 is Ref, and each subsequent
#'  Alt allele is 1 -> n alleles).
#'  \item \code{$COUNTS}, the number of observations of the allele: the number of
#'  individuals for genotype data, or the number of counts (e.g., reads) for
#'  population count data.
#'  \item \code{$INDS}, the number of individuals sampled per population.
#'  \item \code{$FREQ}, the estimated allele frequency.
#'  \item \code{$HET}, the proportion of heterozygotes, calculated directly from
#'  genotype data, or estimated as the expected heteroygosity for population
#'  allele frequencies. Assumes diploid organisms.
#' }
#'
#' @examples
#' library(genomalicious)
#'
#' # Import biallelic SNPs as genotypes or population counts
#' data(data_Genos)
#' data(data_PoolFreqs)
#'
#' # On genotypes, convert the $GT values to characters.
#' dat_gt <- data_Genos %>%
#'   copy %>%
#'   .[, GT:=as.character(GT)] %>%
#'   .[GT==0, GT:='0/0'] %>%
#'   .[GT==1, GT:='0/1'] %>%
#'   .[GT==2, GT:='1/1']
#'
#' print(dat_gt)
#'
#' allele_freqs_DT(dat=dat_gt, type='genos')
#'
#' # On counts, need to make a $COUNTS column, and add in 30 individuals
#' # per locus per population in a new $INDS column.
#' dat_counts <- data_PoolFreqs %>%
#'   copy %>%
#'   .[, COUNTS:=paste(RO,AO,sep=',')] %>%
#'   .[, INDS:=30]
#'
#' print(dat_counts)
#'
#' allele_freqs_DT(dat=dat_counts, type='counts')
#' @export
allele_freqs_DT <- function(
    dat, type, sampCol='SAMPLE', popCol='POP', locusCol='LOCUS',
    genoCol='GT', countCol='COUNTS', indsCol='INDS'
){

  # ----------------------------------------+
  # ASSERTIONS + REFORMAT
  # ----------------------------------------+
  require(data.table); require(tidyverse)

  # Check type
  if(!type %in% c('genos','counts')){
    stop('Argument `type` must be one of "genos" or "counts". See ?allele_freq_DT.')
  }

  # Check columns and reassign
  if(type=='genos'){
    check.cols <- c(sampCol,popCol,locusCol,genoCol) %in% colnames(dat)
    if(length(check.cols)!=4){
      stop('Argument `type`=="genos": At least one of `sampCol`, `popCol`, `locusSol`, or `genoCol` is not in `dat`. See ?allele_freq_DT.')
    } else {
      dat <- dat %>% copy %>% as.data.table %>%
        setnames(
          .,
          old=c(sampCol,popCol,locusCol,genoCol),
          new=c('SAMPLE','POP','LOCUS','GT')
        )
    }
  } else if(type=='counts'){
    check.cols <- c(popCol,locusCol,countCol,indsCol) %in% colnames(dat)
    if(length(check.cols)!=4){
      stop('Argument `type`=="freqs": At least one of `popCol`, `locusCol`, `countCol`, or `indsCol` is not in `dat`. See ?allele_freq_DT.')
    } else{
      dat <- dat %>% copy %>% as.data.table %>%
        setnames(
          .,
          old=c(popCol,locusCol,countCol,indsCol),
          new=c('POP','LOCUS','COUNTS','INDS')
          )
    }
  }

  # ----------------------------------------+
  # MAIN CODE
  # ----------------------------------------+

  # Using genotypes or frequencies?
  if(type=='genos'){
    # Population information
    pop.samps <- unique(dat[,c('POP','SAMPLE')])
    pop.size <- pop.samps[, .(INDS=length(SAMPLE)), by=POP]

    # Split the genotypes
    dat.spl <- dat[, tstrsplit(GT,'/')] %>%
      cbind(dat[, c('SAMPLE','LOCUS')], .)

    # Heterozygote proportions per population and locus
    het.tab <- dat.spl[, .(HET=sum(V1==V2)), by=c('SAMPLE','LOCUS')] %>%
      left_join(., pop.samps) %>%
      .[, .(HET=sum(HET)/length(SAMPLE)), by=c('POP','LOCUS')]

    # Allele frequencies per population and locus
    freq.tab <- data.table::melt(
      dat.spl,
      id.vars=c('SAMPLE','LOCUS'),
      variable.name='VAR',
      value.name='ALLELE'
    ) %>%
      left_join(., pop.samps) %>%
      as.data.table %>%
      .[, .(COUNTS=length(VAR)), by=c('POP','LOCUS','ALLELE')]  %>%
      left_join(., pop.size) %>%
      as.data.table() %>%
      .[, FREQ:=COUNTS/(INDS*2)] %>%
      setorder(., POP, LOCUS)

    # Output
    left_join(freq.tab, het.tab) %>%
      as.data.table %>%
      return()
  } else if(type=='counts'){
    # Population information
    pop.sizes <- dat[, c('POP','INDS')] %>% unique

    # Split the allele counts, then get allele frequencies as
    # proportion of read counts.
    freq.tab <- dat[, tstrsplit(COUNTS, ',')] %>%
      cbind(dat[, c('POP','LOCUS')], .) %>%
      data.table::melt(
        .,
        id.vars=c('POP','LOCUS'),
        variable.name='ALLELE',
        value.name='COUNTS'
      ) %>%
      setorder(., POP, LOCUS) %>%
      .[!is.na(COUNTS)] %>%
      .[, COUNTS:=as.integer(COUNTS)] %>%
      .[, ALLELE:=as.integer(sub('V','',ALLELE))-1] %>%
      .[, FREQ:=COUNTS/sum(COUNTS), by=c('POP','LOCUS')] %>%
      left_join(x=., y=pop.sizes) %>%
      as.data.table()

    # Heterozygosity
    het.tab <- freq.tab %>%
      .[, .(HET=1-sum(FREQ^2)), by=c('POP','LOCUS')] %>%
      left_join(., pop.sizes) %>%
      as.data.table %>%
      .[, .(HET=(INDS/(INDS-1))*HET), by=c('LOCUS','POP')]

    # Output
    left_join(freq.tab, het.tab) %>%
    return()
  }
}
