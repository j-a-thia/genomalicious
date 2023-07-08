#' Converts a BED file of annotations into a GenBank feature table.
#'
#' The GenBank feature table is a 5 column table that is used to submit
#' annotations to GenBank (www.ncbi.nlm.nih.gov/WebSub/html/help/feature-table.html).
#' This function takes a BED file of annotations and converts it into the feature
#' table format for GenBank submission.
#'
#' At this stage, the function's use is largely limited to simple annotations of
#' complete seqeunces and is restricted to genes, mRNA, tRNA, rRNA, and misc_features.
#' Coding of partial sequences and more complex annotations may be added in the future.
#'
#' @param bedFiles Character: A vector of BED file names to process.
#'
#' @param featureRef Data.table: This contains the list of features that are to
#' be mapped against the annotations in BED format. Requires 3 columns:
#' \enumerate{
#'    \item \code{$NAME}, the name of the feature, which matches the 4th column
#'    of the BED file.
#'    \item \code{$FEATURE}, the type of feature, e.g., 'gene', 'tRNA', or 'misc_feature'.
#'    \item \code{$DETAIL}, the details to accompany the feature, which could be
#'    the gene or product name, or notes to attach to misc_features.
#' }
#'
#' @param outFile Character: The output file name for the feature table.
#'
#' @returns Writes the feature table to file.
#'
#' @examples
#' library(genomalicious)
#'
#' # Create a link to raw external datasets in genomalicious
#' genomaliciousExtData <- paste0(find.package('genomalicious'), '/extdata')
#'
#' # We will create a feature table of gene and coding sequences for the protein-coding
#' # genes from the Bathygobius cocosensis mitogenome.
#' bcocoFtRefs <- fread(paste0(genomaliciousExtData, '/data_Bcocosensis_mito_features.csv'))
#' bcocoBedPath <- paste0(genomaliciousExtData, '/data_Bcocosensis_mito_annots.bed')
#'
#' # Take a look
#' bcocoFtRefs
#' bcocoBedPath
#' readLines(bcocoBedPath)
#'
#' # Run
#' bed2featuretab(bcocoBedPath, bcocoFtRefs, 'Bcoco_feature_table.txt')
#'
#' @export

bed2featuretab <- function(bedFiles, featureRef, outFile){
  require(data.table)
  require(tidyverse)

  result <- character()

  for(i in 1:length(bedFiles)){
    # Get the ith BED
    bedfile.i <- bedFiles[[i]]

    # Read in BED
    data.i <- fread(bedfile.i, header=FALSE) %>%
      setnames(., new=c('SAMPLE','START','END','NAME','QUAL','STRAND'))

    # BED files are base 0, and not inclusive of the last base.
    # So you need to add +1 to the start position only.
    data.i[,START:=START+1]

    # Join in the details for the features and make the individual
    # tables for each feature
    feat.i <- left_join(data.i, featureRef) %>%
      as.data.table %>%
      # Give each row an index
      .[, INDEX:=1:.N] %>%
      # Split on index
      split(., by='INDEX') %>%
      # Iterate over each index
      lapply(., function(X){
        # Get strand and assign start and end
        if(X$STRAND=='+'){
          st <- X$START; en <- X$END
        } else if(X$STRAND=='-'){
          st <- X$END; en <- X$START
        }
        # Get the feature
        ft <- X$FEATURE
        # Get the detail
        deet <- X$DETAIL
        # 1st column: start
        col1 <- c(st, '') %>% as.character
        # 2nd column: end
        col2 <- c(en, '') %>% as.character
        # 3rd column: gene, or empty for note
        col3 <- c(ft, '') %>% as.character
        # 4th column
        if(ft %in% c('misc_feature')){
          col4 <- c('', 'note')
        } else if(ft=='gene'){
          col4 <- c('', 'gene')
        } else if(ft %in% c('CDS','mRNA','tRNA','rRNA')){
          col4 <- c('', 'product')
        }
        # 5th column
        col5 <- c('', deet) %>% as.character
        # Combine
        cbind(col1, col2, col3, col4, col5)
      }) %>%
      # Merge into single table
      do.call('rbind', .) %>%
      # But now iterate over rows and con
      apply(., 1, function(Y){
        paste(Y, collapse='\t')
      })

    # Add in header for ith sample
    result <- c(result,c(paste0('>Feature ', data.i$SAMPLE[1]), feat.i), '')
  }

  writeLines(result, outFile)
}
