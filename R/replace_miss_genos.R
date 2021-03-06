#' Replace missing genotypes
#'
#' For each locus, missing genotypes are replaced with either
#' the median or mode genotype. Can be done across all sampled individuals
#' or by population. Loci must be biallelic. NOTE: it is recommended
#' that missing genotypes are imputed, but if you need a
#' quick-and-dirty approach, this function might be useful for
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
#' @param methodReplace Character: The way in which missing genotypes
#' are replaced. Either by the code{'median'} (middle) genotype, or by
#' the \code{'mode'} (most common) genotype. Default = \code{'mode'}.
#'
#' @details If genotypes are coded as characters, \code{NA} or \code{'./.'}
#' should be used to code missing genotypes. Otherwise if genotypes
#' are coded as integers, \code{NA} should code missing genotypes.
#'
#' @examples
#' # Create a link to raw external datasets in genomalicious
#' genomaliciousExtData <- paste0(find.package('genomalicious'), '/inst/extdata')
#'
#' # Use this to create a path to the genomalicious ind-seq VCF
#' vcfPath <- paste0(genomaliciousExtData, '/data_indseq.vcf')
#'
#' # Import the VCF
#' indseq <- vcf2DT(vcfPath)
#' indseq
#'
#' # Add a population column
#' indseq[, POP:=sub(
#'         pattern='Ind'
#'         , replacement=''
#'         , x=strsplit(SAMPLE, '.', fixed=TRUE)[[1]][1])
#'     , by=c('LOCUS', 'SAMPLE') ]
#'
#' # Sites with missing data
#' genomiss <- which(indseq$GT=='./.')
#' genomiss
#' indseq[genomiss,]
#'
#' # Replace genotypes using the mode across all individuals
#' indseq_mode_all <- replace_miss_genos(
#'     dat=indseq
#'     , sampCol='SAMPLE'
#'     , locusCol='LOCUS'
#'     , popCol=NULL
#'     , methodReplace='mode')
#'
#' # Replace genotypes using the mode within each population
#' indseq_mode_pops <- replace_miss_genos(
#'     dat=indseq
#'     , sampCol='SAMPLE'
#'     , locusCol='LOCUS'
#'     , popCol='POP'
#'     , methodReplace='mode')
#'
#' # Replace genotypes using the median within each population
#' indseq_med_pops <- replace_miss_genos(
#'     dat=indseq
#'     , sampCol='SAMPLE'
#'     , locusCol='LOCUS'
#'     , popCol='POP'
#'     , methodReplace='median')
#'
#' # Compare between methods
#' cbind(indseq_mode_all[genomiss, GT]
#'     , indseq_mode_pops[genomiss, GT]
#'     , indseq_med_pops[genomiss, GT])
#'
#' @export
replace_miss_genos <- function(dat
                               , sampCol='SAMPLE'
                               , locusCol='LOCUS'
                               , genoCol='GT'
                               , popCol=NULL
                               , methodReplace='mode'){

  # BEGIN ...........

  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  require(data.table)

  # Check the replacement method is specified correctly
  if((methodReplace %in% c('median', 'mode'))==FALSE){
    stop("Argument `methodReplace` misspecififed. See: ?replace_miss_genos")
  }

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

  # --------------------------------------------+
  # Internal functions
  # --------------------------------------------+
  FUN_median_geno <- function(x){
    return(x[round(length(x)/2)])
  }

  FUN_mode_geno <- function(x){
    tab <- sort(table(x))
    return(names(tab)[length(tab)])
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Split the data by loci, or loci * population combinations
  if(is.null(popCol)==TRUE){
    dat <- split(dat, dat[[locusCol]])
  } else if(is.null(popCol)==FALSE){
    dat <- split(dat, dat[, c(locusCol, popCol), with=FALSE])
  }

  # Iterate through each X facet of data
  datRepl <- lapply(dat, function(X){
    # Get the genotypes, omit missing values
    x_gt <- sort(na.omit(X[[genoCol]]))

    # Determine which value to replace missing genotypes
    if(methodReplace=='median'){ x_repl <- FUN_median_geno(x_gt)
    } else if(methodReplace=='mode'){ x_repl <- FUN_mode_geno(x_gt) }

    # Replace missing genotypes
    X[[genoCol]][is.na(X[[genoCol]])] <- x_repl

    # Return updated value
    return(X)
  })

  return(do.call('rbind', datRepl))
  # .......... END
}
