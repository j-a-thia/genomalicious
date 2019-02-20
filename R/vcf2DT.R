#' VCF file to data table
#'
#' Reads a VCF file and converts to a data table.
#'
#' @param vcfFile Character: The path to the input VCF file.
#'
#' @return A data.table with all the columns contained in the original VCF file with
#' some additions:
#' \itemize{
#'     \item A column called \code{LOUCS} is generated. This is the concatenation of the
#'              \code{CHROM} and \code{POS} column to form a locus ID.
#'     \item A column called \code{SAMPLE} is generated. This contains the sample IDs that
#'              are the columns that follow the \code{FORMAT} column in the original VCF.
#'     \item The items in the original \code{FORMAT} column of the VCF are given their own columns.
#' } \cr
#' Note, for VCF files produced by Stacks, the $CHROM is given the same value as the $ID column.
#'
#' @examples
#' # Create a link to raw external datasets in pgposer
#' pgposerExtData <- paste0(find.package('pgposer'), '/extdata')
#'
#' # This command here shows you the VCF file that comes with pgposer
#' list.files(pgposerExtData, pattern='vcf')
#'
#' # Use this to create a path to that file
#' vcfPath <- paste0(pgposerExtData, '/pgposer.vcf.txt')
#'
#' # You can read the file in as lines to see what it
#' # looks like:
#' readLines(vcfPath)
#'
#' # Now read it in as a data.table
#' vcf2DT(vcfPath)
#'
#' @export
vcf2DT <- function(vcfFile){

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # #### Libraries and assertions            ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  require(data.table)

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # #### Code                                ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # What is the position of the header?
  headPos <- grep('#CHROM', readLines(vcfFile), value=FALSE)

  # Read file from header
  vcfDT <- fread(vcfFile, skip=headPos-1, sep='\t', header=TRUE)

  # Adjust header
  colnames(vcfDT) <- gsub(pattern='#', replace='', x=colnames(vcfDT))

  # Which columns are the sample? The ones after the FORMAT column.
  sampCols <- (which(colnames(vcfDT)=='FORMAT')+1):ncol(vcfDT)

  # Now convert the data from wide to long
  vcfDT <- melt(data=vcfDT, id.vars=1:(sampCols[1]-1), measure.vars=sampCols, variable.name='SAMPLE', value.name='DATA')

  # Make sure SAMPLE is a character
  vcfDT$SAMPLE <- as.character(vcfDT$SAMPLE)

  # Generate a LOCUS column
  vcfDT[, LOCUS:=paste(CHROM, POS, sep='_')]

  # Extract the first value in the FORMAT column as a vector
  format <- unlist(strsplit(vcfDT$FORMAT[1], ':'))

  # Iterate through each index in 'format', then access the DATA column from vcfDT
  # and pull out that specific value.
  for(f in 1:length(format)){
    vcfDT[, format[f]] <- unlist(lapply(strsplit(vcfDT$DATA, ':'), function(d){ d[f] }))
  }

  # Drop FORMAT and DATA
  vcfDT <- vcfDT[, !c('FORMAT','DATA'), with=F]

  # Make sure DP, RO, and AO are integers
  for(i in c('DP', 'RO', 'AO')){
    if(i %in% colnames(vcfDT)){ vcfDT[[i]] <- as.integer(vcfDT[[i]]) }
  }

  # Return the data.table
  return(vcfDT)
}

