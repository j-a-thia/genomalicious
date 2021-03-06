#' VCF file to data table
#'
#' Reads a VCF file and converts to a long format data table.
#' Also provides functionality to go in reverse,
#' i.e. writing a VCF from a long format data table.
#'
#' @param vcfFile Character: The path to the input VCF file.
#'
#' @param dropCols Character: Vector of column names from the VCF that you
#' want to drop from the data table. Default = \code{NULL}.
#' Only relevant when argument \code{flip==FALSE}.
#'
#' @param keepComments Logical: Should the VCF comments be kept?
#' Default = \code{FALSE}. See Details for parameterisation.
#'
#' @param keepInfo: Logical: Should the VCF info for each locus be kept?
#' Default = \code{FALSE}. See Details for parameterisation.
#'
#' @param flip Logical: Go in the reverse direction, i.e. data table to VCF?
#' Default = \code{FALSE}, i.e. VCF to data table. See Details for parameterisation.
#'
#' @param dat Data table: Long format data table containing read and genotype information.
#' Ideally, a data table produced with \code{vcf2DT(..., flip=FALSE)}. This
#' argument is only relevant when \code{flip==TRUE} and should contain all the
#' columns that are specified in argument \code{vcfValues}, and optionally,
#' the attributes \code{'vcf_comments'} and \code{'vcf_info'}. See Details.
#'
#' @param vcfValues List: Information about relevant columns for constructing
#' a VCF file from a long format data table. Only relevant when \code{flip==TRUE}.
#' All these values must be specific column names in argument \code{dat}.
#' Structure with the following indices,
#' \enumerate{
#'     \item \code{$locusCol}: Character, single value specifying the column
#'         of locus names. Default = \code{'LOCUS'}.
#'     \item \code{$variantCols}: Character, 1+ values specifying information
#'         on variant sites. Default = \code{c('CHROM', 'POS', 'REF', 'ALT')}.
#'     \item \code{$formatCols}: Character, 1+ values specifying information
#'         on the format of sample data. Default = \code{c('DP', 'RO', 'AO')}.
#'    \item \code{$sampCol}: Character, a single value specifying the column
#'         of sample names. Default = \code{'SAMPLE'}.
#' }
#'
#' @details Firstly, it should be noted that while data tables are a really
#' excellent way of handling genotype and sequence read information in R,
#' they are not necessarily the most efficient way to do so. Importing VCFs
#' as data table (or the reverse, exporting data tables as VCFs), can take
#' a considerable amount of time if the number of loci and samples are large.
#' However, a bit of patience is worth it! \cr
#'
#' The direction of this function, i.e. importing VCFs or exporting VCFs,
#' is controlled by the argumnet \code{flip}. The remaineder of this section
#' will describe the nuances of parameterising import/export. \cr
#'
#' If \code{flip==FALSE}, this converts a VCF file into a long format data table.
#' Parameterisation can be simple, just specify the VCF file path with
#' \code{vcfFile}. More fine-scale control can be exerted by removing specific
#' columns using the argument \code{dropCols}. Note, however, that the FORMAT
#' and INFO columns from the VCF are always removed. The VCF INFO column contains
#' specifics of the called variants, which you may want to keep by specifying
#' the argument \code{keepInfo==TRUE}. Doing so saves INFO as an attribute, such
#' that if \code{dt} is the returned data table, \code{attr(dt, 'vcf_info')}
#' will return a vector of the original VCF INFO values for each locus.
#' If you want to keep the leading comments of the VCF (lines starting with '##'),
#' then specifying \code{keepComments==TRUE} saves these lines as an attribute,
#' accessible with \code{attr(dt, 'vcf_comments')}. Keeping the INFO column and
#' comment lines of the original VCF is important if you want to access these later,
#' particularly if you want to convert the data table back into a VCF.
#'
#' If \code{flip==TRUE}, this converts a long format data table into a VCF file.
#' To do so requires specification of the data table using argument \code{dat}
#' and the name of the target saved VCF file using \code{vcfFile}. Genotypes and
#' read counts of variants are stored in wide format in VCF files
#' (see Thia & Riginos 2019 BioRxiv for illustration). Each ith row is a locus
#' and specific details of variant and the jth samples are in columns. Specification
#' of where this data is stored in \code{dat} is managed by the argument
#' \code{vcfValues}, which is a list containing the indices \code{$loci},
#' $\code{variants}, \code{$format}, and \code{$samples}: see the Examples.
#' Note: (1) if you change \code{vcfValues} away from default values in any way,
#' all indices of \code{vcfValues} must be manually specified; and (2) at minimum
#' you should always keep the CHROM and POS columns, \code{vcfValues$variants=c('CHROM', 'POS')},
#' because these are import if you want to re-import the VCF with \code{vcf2DT()}.
#' The VCF file can be constructed with leading comment lines and the variant
#' INFO column, using the arguments \code{keepComments=TRUE} and \code{keepInfo=TRUE},
#' respectively. However, in order for this to occur, the data table, \code{dt},
#' specified in \code{dat}, must have the attributes \code{attr(dt, 'vcf_comments')}
#' and \code{attr(dt, 'vcf_info')}, respectively. These attributes can be created
#' when importing a VCF file with \code{vcf2DT()}, see above.
#'
#' Note, writing VCFs from data tables is very slow. Going from a long format to
#' wide format for many samples and loci is memory intensive. The current algorithm
#' deals with each unique locus individually and in order, writing on the fly
#' to avoid loading #' extremely large wide format datasets into memory.
#' Would gladly accept any suggestions on how to improve the algorithm!
#'
#' @return When \code{flip==FALSE}, i.e. converting a VCF to data table, a
#' \code{data.table} object is returned with all the columns contained in
#' the original VCF file with some additions:
#' \itemize{
#'     \item A column called \code{LOUCS} is generated. This is the concatenation of the
#'              \code{CHROM} and \code{POS} column to form a locus ID.
#'     \item A column called \code{SAMPLE} is generated. This contains the sample IDs that
#'              are the columns that follow the \code{FORMAT} column in the original VCF.
#'     \item The items in the original \code{FORMAT} column of the VCF are given their own columns.
#' } \cr
#' Note, for VCF files produced by Stacks, the $CHROM is given the same value
#' as the $ID column. \cr\cr
#' When \code{flip==TRUE}, i.e. converting a long format data table to VCF, a
#' text file in VCF format is written to the current working directory, or
#' the path specified in argument \code{vcfFile}.
#'
#' @references This & Riginos (2019) genomalicious: serving up a smorgasbord of
#' R functions for population genomic analyses. BioRxiv.
#'
#' @examples
#' # Create a link to raw external datasets in genomalicious
#' genomaliciousExtData <- paste0(find.package('genomalicious'), '/extdata')
#'
#' # This command here shows you the VCF file that comes with genomalicious
#' list.files(genomaliciousExtData, pattern='_poolseq.vcf')
#'
#' # Use this to create a path to that file
#' vcfPath <- paste0(genomaliciousExtData, '/genomalicious_poolseq.vcf')
#'
#' # You can read the file in as lines to see what it
#' # looks like:
#' readLines(vcfPath)
#'
#' # Now read it in as a data table
#' readVcf1 <- vcf2DT(vcfPath)
#' readVcf1
#'
#' # Read in VCF, but drop some columns,
#' # and keep comments and info.
#' readVcf2 <- vcf2DT(vcfPath
#'    , dropCols=c('FILTER', 'ID')
#'    , keepComments=TRUE
#'    , keepInfo=TRUE)
#'
#' readVcf2
#'
#' attr(readVcf2, 'vcf_comments')
#' attr(readVcf2, 'vcf_info')
#'
#' # Convert a data table back to a VCF.
#' # Code shows how you can exlude columns from the VCF.
#' # Here, will drop the $QUAL and $DP columns.
#' col_locus <- 'LOCUS'
#' col_var <- c('CHROM', 'POS', 'REF', 'ALT')
#' col_form <- c('RO', 'AO')
#' col_samp <- 'SAMPLE'
#' vcfValList <- list(loci=col_locus
#'     , variants=col_var
#'     , format=col_form
#'     , samples=col_samp)
#'
#' vcf2DT(vcfFile='dt2vcf_example.vcf'
#'     , keepComments=TRUE
#'     , keepInfo=TRUE
#'     , flip=TRUE
#'     , dat=readVcf2
#'     , vcfValues=vcfValList)
#'
#' readLines('dt2vcf_example.vcf')
#'
#' vcf2DT('dt2vcf_example.vcf', flip=FALSE)
#'
#' @export
vcf2DT <- function(vcfFile
                   , dropCols=NULL
                   , keepComments=FALSE
                   , keepInfo=FALSE
                   , flip=FALSE
                   , dat
                   , vcfValues=list(loci='LOCUS'
                         , variants=c('CHROM', 'POS', 'REF', 'ALT')
                         , format=c('GT', 'DP', 'RO', 'AO')
                         , samples='SAMPLE'))
{

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # #### Libraries and assertions            ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  require(data.table)

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # #### Code: VCF to data table             ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(flip==FALSE){
    # What is the position of the header?
    headPos <- grep('#CHROM', readLines(vcfFile), value=FALSE)

    # Read file from header
    cat('(1/4) Reading in VCF as a data table', sep='\n')
    vcfDT <- fread(vcfFile, skip=headPos-1, sep='\t', header=TRUE)

    # Adjust header
    colnames(vcfDT) <- gsub(pattern='#', replace='', x=colnames(vcfDT))

    # Generate a $LOCUS column, place at the start of the data table
    cat('(2/4) Generating locus IDs', sep='\n')
    vcfDT <- cbind(LOCUS=vcfDT[, paste(CHROM, POS, sep='_')], vcfDT)

    # Get the locus info as a vector and drop from data table
    vcfInfo <- vcfDT$INFO
    names(vcfInfo) <- vcfDT$LOCUS
    vcfDT <- vcfDT[, !'INFO']

    # Which columns are the sample? The ones after the FORMAT column.
    sampCols <- (which(colnames(vcfDT)=='FORMAT')+1):ncol(vcfDT)

    # Now convert the data from wide to long
    cat('(3/4) Converting from wide to long format', sep='\n')
    vcfDT <- melt(data=vcfDT, id.vars=1:(sampCols[1]-1), measure.vars=sampCols, variable.name='SAMPLE', value.name='DATA')

    # Make sure SAMPLE is a character
    vcfDT$SAMPLE <- as.character(vcfDT$SAMPLE)

    # Add the $FORMAT data into the data table
    cat('(4/4) Collecting and organising FORMAT data', sep='\n')

    # Split the $FORMAT column, identify NAs, rotate, rename columns, and bind.
    formatDat <- t(vcfDT[, strsplit(DATA, ':')])
    formatDat[which(formatDat=='.')] <- NA
    colnames(formatDat) <- unlist(strsplit(vcfDT$FORMAT[1], ':'))
    vcfDT <- cbind(vcfDT, as.data.table(formatDat))

    # Drop FORMAT and DATA
    vcfDT <- vcfDT[, !c('FORMAT','DATA'), with=FALSE]

    # Make sure DP, RO, and AO are integers
    for(i in c('DP', 'RO', 'AO')){
      if(i %in% colnames(vcfDT)){ vcfDT[[i]] <- as.integer(vcfDT[[i]]) }
    }

    # Attach header as an attribute, if specified.
    if(keepComments==TRUE){
      attr(vcfDT, 'vcf_comments') <- readLines(vcfFile, n=headPos-1)
    }

    # Attach info as an attribute if, if specified.
    if(keepInfo==TRUE){
      attr(vcfDT, 'vcf_info') <- vcfInfo
    }

    # Return the data.table, drop any columns if specified.
    if(is.null(dropCols)){
      return(vcfDT)
    } else{
      return(vcfDT[, !dropCols, with=FALSE])
    }

    cat('All done! <3', '\n')
  }

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # #### Code: Data table to VCF             ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(flip==TRUE){
    # Check that all columns in vcfValues are in dat.
    if(sum(unlist(vcfValues) %in% colnames(dat))!=length(unlist(vcfValues))){
      stop('Not all columns specified in `vcfValues` match `colnames(dat)`: See ?vcf2DT.')
    }

    # Check that CHROM and POS are specified in vcfValues$variants,
    # and put CHROM and POS first if they are not already.
    if(sum(c('CHROM', 'POS') %in% vcfValues$variants)!=2){
      stop("Argument `vcfValues$variants` must at minimum contain the values
           'CHROM' and 'POS': See ?vcf2DT.")
    } else{
      chrompos <- match(c('CHROM', 'POS'), vcfValues$variants)
      vcfValues$variants <- c(vcfValues$variants[chrompos]
                              , vcfValues$variants[-chrompos])
    }

    # Internal function
    FUN_format_paste <- function(x){
      x[is.na(x)] <- '.'
      return(paste(trimws(x), collapse=':'))
    }

    # The unquie loci and samples
    uniqLoci <- unique(dat[, c('CHROM', 'POS', 'LOCUS')])
    setorder(uniqLoci, CHROM, POS)

    uniqSamps <- sort(unique(dat$SAMPLE))

    # Locus counter
    countLoci <- round(quantile(1:nrow(uniqLoci), seq(0.05, 1, by=0.05)))
    names(countLoci) <- sub('%', '', names(countLoci))

    # Create a data table of variants (varDat).
    cat('(1/4) Collecting the variant data.', '\n')
    varDat <- unique(dat[,c(vcfValues$variants, vcfValues$loci), with=FALSE])

    # Should INFO and comments be kept?
    cat('(2/4) Collecting info and comments if specified.', '\n')
    # Extract VCF INFO column
    if(keepInfo==TRUE){
      outInfo <- attr(dat, 'vcf_info')
      if(is.null(outInfo)==TRUE){
        cat("   NOTE: `keepInfo==TRUE` but `attr(dat, 'vcf_info')` is NULL.", sep='\n')
      }
      varDat$INFO <- outInfo[varDat$LOCUS]
    } else{ outInfo <- NULL }

    # Extract VCF comments
    if(keepComments==TRUE){
      # Comments stored in attributes
      vcfComms <- attr(dat, 'vcf_comments')
      # The comments to be out put
      outComms <- lapply(vcfValues$format, function(f){
        vcfComms[grep(pattern=paste0('##FORMAT=<ID=', f), x=vcfComms)]
      })
      outComms <- unlist(outComms)

      if(is.null(outComms)==TRUE){
        cat("NOTE: `keepComments==TRUE` but `attr(dat, 'vcf_comments')` is NULL.", sep='\n')
      }

      # If INFO is desired, add to comments
      if(keepInfo==TRUE){
        outComms <- c(vcfComms[grep(pattern='##INFO', x=vcfComms)], outComms)
      }
    } else{ outComms <- NULL }

    # Create a data table of sample data (sampDat).
    # Will iterate through loci and write lines on the fly.
    cat('(3/4) Collecting the sample data.', '\n')
    sampDat <- dat[, c(vcfValues$loci, vcfValues$samples, vcfValues$format), with=FALSE]

    # Create the VCF
    cat('(4/4) Writing the VCF.', '\n')

    # Write info, comments, and column header.
    vcfHead <- paste0('#', paste(c(vcfValues$variants, 'FORMAT', uniqSamps), collapse='\t'))
    writeLines(c(outComms, vcfHead), vcfFile)

    # Write the VCF lines for each ith locus.
    cat('   % complete: ')
    for(i in 1:length(uniqLoci)){
      locus <- uniqLoci$LOCUS[i]

      # Make count if at checkpoint
      if(i %in% countLoci){ cat(names(countLoci[countLoci==i]), ' ') }

      # Subset data based on locus: get desired format values
      # and the sample names.
      locDT <- sampDat[LOCUS==locus, c(vcfValues$format, vcfValues$samples), with=FALSE]

      # Rotate into wide format such that samples are in columns,
      loc_samp_vals <- as.data.table(t(locDT[, !vcfValues$samples, with=FALSE]))

      # Iterate over all columns with .SDcols and apply
      # FUN_formate_paste() to each, which combines all format
      # values into a single string.
      loc_samp_vals <- loc_samp_vals[, lapply(.SD, FUN_format_paste)
                                     , .SDcols=colnames(loc_samp_vals)]

      # Label columns based on samples
      colnames(loc_samp_vals) <- locDT$SAMPLE

      # Reorganise to keep sample name consistency
      loc_samp_vals <- loc_samp_vals[, uniqSamps, with=FALSE]

      # Write the new line to file.
      vcfLine <- cbind(varDat[LOCUS==locus, !'LOCUS']
                       , FORMAT=paste(vcfValues$format, collapse=':')
                       , loc_samp_vals)

      fwrite(x=vcfLine, file=vcfFile, append=TRUE, sep='\t')
    }
    cat('\n')

    cat('All done! <3', '\n')
  }
}


