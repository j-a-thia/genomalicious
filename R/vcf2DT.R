#' VCF file to data table
#'
#' Reads a VCF file and converts to a long format data table. Note, that whilst
#' the \code{data.table} object class is very memory efficient, very large genomic
#' datasets might take longer to read in, and/or be difficult to hold in
#' memory. Take your operating system and the size of your input dataset into
#' consideration when using this function.
#'
#' @param vcfFile Character: The path to the input VCF file.
#'
#' @param dropCols Character: Vector of column names from the VCF that you
#' want to drop from the output data table. Use this for any column that occurs
#' before the 'FORMAT' column in the original VCF file. Default = \code{NULL}.
#'
#' @param keepComments Logical: Should the VCF comments be kept?
#' Default = \code{FALSE}. See Details for parameterisation.
#'
#' @param keepInfo Logical: Should the VCF info for each locus be kept?
#' Default = \code{FALSE}.
#'
#' @details Firstly, it should be noted that while data tables are a really
#' excellent way of handling genotype and sequence read information in R,
#' they are not necessarily the most efficient way to do so for very large
#' genomic datasets. Take your operating system and/or dataset in mind before
#' using this function. Most RADseq datasets should be manageable, but
#' whole-genome data can be challenging if you do not have a lot of available
#' memory. You can always try loading in subsets (e.g., by chromosome or contigs)
#' of your dataset to see how feasible it is to load with this function.
#'
#' @return A \code{data.table} object is returned with all the columns contained in
#' the original VCF file with some additions:
#' \itemize{
#'     \item A column called \code{$LOCUS} is generated. This is the concatenation of the
#'              \code{$CHROM} and \code{$POS} column to form a locus ID.
#'     \item A column called \code{$SAMPLE} is generated. This contains the sample IDs that
#'              are the columns that follow the \code{$FORMAT} column in the original VCF.
#'     \item The items in the original \code{$FORMAT} column of the VCF are given their own columns.
#' } \cr
#' Note, for VCF files produced by Stacks, the $CHROM is given the same value
#' as the $ID column. \cr\cr
#' When \code{keepInfo==TRUE} and/or \code{keepComments==TRUE}, these are returned
#' as attributes. E.g., if the returned object is \code{vcfDT}, then you can
#' access Info and Comments (respectively) with: \code{attr(vcfDT, 'vcf_info')}
#' and \code{attr(vcfDT, 'vcf_comments')}.
#'
#' @examples
#' # Create a link to raw external datasets in genomalicious
#' genomaliciousExtData <- paste0(find.package('genomalicious'), '/extdata')
#'
#' # This command here shows you the VCF file that comes with genomalicious
#' list.files(path=genomaliciousExtData, pattern='indseq.vcf')
#'
#' # Use this to create a path to that file
#' vcfPath <- paste0(genomaliciousExtData, '/data_indseq.vcf')
#'
#' # You can read the file in as lines to see what it
#' # looks like:
#' readLines(vcfPath) %>%  head
#' readLines(vcfPath) %>%  tail
#'
#' # Now read it in as a data table
#' readVcf1 <- vcf2DT(vcfFile=vcfPath)
#' readVcf1 %>% print()
#'
#' # Read in VCF, but drop some columns,
#' # and keep comments and info.
#' readVcf2 <- vcf2DT(vcfPath
#'    , dropCols=c('QUAL')
#'    , keepComments=TRUE
#'    , keepInfo=TRUE)
#'
#' readVcf2 %>% print
#'
#' attr(readVcf2, 'vcf_comments')
#' attr(readVcf2, 'vcf_info')
#'
#' @export
vcf2DT <- function(vcfFile, dropCols=NULL, keepComments=FALSE, keepInfo=FALSE){

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # #### Libraries and assertions            ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  require(data.table); require(tidyverse)

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # #### Code: VCF to data table             ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # What is the position of the column heads?
  n <- 1
  headPos <- NULL
  while(length(headPos)==0){
    headPos <- grep('#CHROM', readLines(vcfFile, n=n), value=FALSE)
    if(length(headPos)==0){
      n <- n + 1000
    }
  }

  # Read file from header
  cat('(1/4) Reading in VCF as a data table', sep='\n')
  vcfDT <- fread(vcfFile, skip=headPos-1, sep='\t', header=TRUE)

  # Adjust header
  colnames(vcfDT) <- gsub(pattern='#', replace='', x=colnames(vcfDT))

  # Which columns are the sample? The ones after the FORMAT column.
  sampCols <- colnames(vcfDT)[(which(colnames(vcfDT)=='FORMAT')+1):ncol(vcfDT)]

  # Generate a $LOCUS column, place at the start of the data table
  cat('(2/4) Generating locus IDs', sep='\n')
  vcfDT[, LOCUS:=paste0(CHROM, '_', POS)]

  # Get the locus info as a vector and drop from data table
  if(keepInfo==TRUE){
    vcfInfo <- vcfDT$INFO
    names(vcfInfo) <- vcfDT$LOCUS
  }
  vcfDT <- vcfDT[, !'INFO']

  # Drop unwanted columns here to save memory
  if(is.null(dropCols)==FALSE){
    vcfDT <- vcfDT[, !dropCols, with=FALSE]
  }

  # Now convert the data from wide to long
  cat('(3/4) Converting from wide to long format', sep='\n')
  vcfDT <- data.table::melt(
    data=vcfDT,
    id.vars=colnames(vcfDT)[!colnames(vcfDT)%in%sampCols],
    measure.vars=sampCols,
    variable.name='SAMPLE',
    value.name='DATA') %>%
    as.data.table()

  # Make sure SAMPLE is a character
  vcfDT[, SAMPLE:=as.character(SAMPLE)]

  # Separate out the FORMAT data components into their own columns
  cat('(4/4) Parsing data for each sample', sep='\n')

  # ... Get the format names
  formatNames <- unlist(strsplit(vcfDT$FORMAT[1], split=':'))

  # ... Drop the format column to save space now
  vcfDT <- vcfDT[, !'FORMAT']

  # ... If the $DATA column is '.', add in NA
  vcfDT[DATA=='.', DATA:=NA]

  # ... Separate $DATA by $FORMAT names
  vcfDT <- vcfDT[, tstrsplit(DATA, ':', names=formatNames)] %>%
    cbind(vcfDT[, !'DATA'], .) %>%
    as.data.table()

  # .... Replace '.' values in the data columns with NA
  for(f in formatNames){
    vcfDT[[f]][vcfDT[[f]]=='.'] <- NA
  }

  # Make sure DP and RO integers
  if('DP' %in% colnames(vcfDT)){ vcfDT[, DP:=as.integer(DP)] }
  if('RO' %in% colnames(vcfDT)){ vcfDT[, RO:=as.integer(RO)] }

  # Attach header as an attribute, if specified.
  if(keepComments==TRUE){
    attr(vcfDT, 'vcf_comments') <- readLines(vcfFile, n=headPos-1)
  }

  # Attach info as an attribute if, if specified.
  if(keepInfo==TRUE){
    attr(vcfDT, 'vcf_info') <- vcfInfo
  }

  # Finish
  cat('All done! <3', '\n')

  # Return the data.table, drop any columns if specified.
  return(vcfDT)
}
