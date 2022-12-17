#' Replace missing genotypes
#'
#' For each locus, missing genotypes are replaced with the most common
#' genotype. Can be done across all sampled individuals or by population.
#' Loci must be biallelic.
#'
#' NOTE: it is recommended that missing genotypes are imputed using
#' inferences of linkage and genotype likelihood. However, if you need
#' a quick-and-dirty approach, this function might be useful for
#' preliminary analyses, or if missing data is very low.
#'
#' @param dat Data table: A long data table, e.g. like that imported from
#' \code{vcf2DT}. Genotypes can be coded as '/' separated characters
#' (e.g. '0/0', '0/1', '1/1'), or integers as Alt allele counts (e.g. 0, 1, 2).
#' Must contain the following columns,
#' \enumerate{
#'   \item The sampled individuals (see param \code{sampCol}).
#'   \item The locus ID (see param \code{locusCol}).
#'   \item The genotype column (see param \code{genoCol}).
#' }
#'
#' @param sampCol Character: The column name with the sampled individual information.
#' Default is \code{'SAMPLE'}.
#'
#' @param locusCol Character: The column name with the locus information.
#' Default is \code{'LOCUS'}.
#'
#' @param genoCol Character: The column name with the genotype information.
#' Default is \code{'GT'}.
#'
#' @param popCol Character: An optional argument. The column name with the
#' population information. Default is \code{NULL}. If specified, genotype
#' replacement at each locus is done per population, not across all
#' sampled individuals.
#'
#' @details If genotypes are coded as characters, \code{NA} or \code{'./.'}
#' should be used to code missing genotypes. Otherwise if genotypes
#' are coded as integers, \code{NA} should code missing genotypes.
#' Whether the most common genotype is estimated across individuals or
#' for each population depends on parameterisation of \code{popCol}.
#'
#' @examples
#' library(genomalicious)
#'
#' data(data_4pops)
#'
#' D <- data_4pops %>% copy
#'
#' # Sites with missing data
#' D[sample(1:nrow(D), round(0.1*nrow(D)), FALSE), GT:=NA] %>%
#'  setnames(., 'GT', 'GT.MISS')
#'
#' # Replace across individuals
#' D.rep.inds <- replace_miss_genos(
#'    dat=D, sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT.MISS'
#' ) %>%
#'    setnames(., 'GT', 'GT.INDS')
#'
#' # Replace within populations
#' D.rep.pops <- replace_miss_genos(
#'    dat=D, sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT.MISS', popCol='POP'
#' ) %>%
#'    setnames(., 'GT', 'GT.POPS')
#'
#' # Tabulate comparisons between methods
#' compReplace <- left_join(
#'    data_4pops[, c('LOCUS','SAMPLE','POP','GT')],
#'    D[, c('LOCUS','SAMPLE','POP','GT.MISS')]
#' ) %>%
#' .[is.na(GT.MISS), !'GT.MISS'] %>%
#'    left_join(., D.rep.inds[,c('LOCUS','SAMPLE','POP','GT.INDS')]) %>%
#'    left_join(., D.rep.pops[,c('LOCUS','SAMPLE','POP','GT.POPS')])
#'
#' # Number of correct matches is slightly higher when using the most
#' # common genotype within populations
#' compReplace[GT==GT.INDS] %>% nrow
#' compReplace[GT==GT.POPS] %>% nrow
#'
#'
#' @export
replace_miss_genos <- function(
    dat, sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT', popCol=NULL
  ){

  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  require(data.table)

  # Get the class of the genotypes
  gtClass <- class(dat[[genoCol]])

  # Check that genotypes are characters or counts
  if(!gtClass %in% c('character', 'numeric', 'integer')){
    stop("Check that genotypes are coded as '/' separated characters or as
         counts of the Alt allele. See: ?replace_miss_genos")
  }

  # Turn missing data into NAs if genotypes are characters
  if(gtClass=='character'){
    dat[[genoCol]][dat[[genoCol]]=='./.'] <- NA
  }

  # Convert numeric genotypes into integers
  if(gtClass=='numeric'){
    dat[[genoCol]] <- as.integer(dat[[genoCol]])
  }

  # Rename columns
  if(is.null(popCol)){
    if(sum(colnames(dat) %in% c(sampCol, locusCol, genoCol))!=3){
      stop('Arguments `sampCol`, `locusCol`, and `genoCol` must all be
           column names in argument `dat`. See ?replace_miss_genos')
    }

    colnames(dat)[match(c(sampCol, locusCol, genoCol), colnames(dat))] <- c('SAMPLE','LOCUS','GT')
  } else if(!is.null(popCol)){
    if(sum(colnames(dat) %in% c(sampCol, locusCol, genoCol, popCol))!=4){
      stop('Arguments `sampCol`, `locusCol`, `genoCol`, and `popCol` must all be
           column names in argument `dat`. See ?replace_miss_genos')
    }

    colnames(dat)[match(c(sampCol, locusCol, genoCol, popCol), colnames(dat))] <- c('SAMPLE','LOCUS','GT','POP')
  }

  # --------------------------------------------+
  # Internal functions
  # --------------------------------------------+
  FUN_common_geno <- function(x){
    xclass <- class(x)
    tab <- sort(table(x))
    xcommon <- names(tab)[length(tab)]
    return(as.integer(xcommon))
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  if(is.null(popCol)){
    genoCommTab <- dat[, .(GT=FUN_common_geno(GT)), by=LOCUS]
  } else if(!is.null(popCol)){
    genoCommTab <- dat[, .(GT=FUN_common_geno(GT)), by=c('LOCUS','POP')]
  }

  datGenos <- dat[!is.na(GT)]
  datMiss <- dat[is.na(GT), !'GT']

  rbind(datGenos, left_join(datMiss, genoCommTab)) %>%
    setorder(., LOCUS, SAMPLE) %>%
    return()
  # .......... END
}
