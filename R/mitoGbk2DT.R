#' Convert a mitogenome GENBANK file into a data table
#'
#' This function reads in a GENBANK file and outputs an R data.table object.
#' It has been written for the import of mitogenome data into R, but could be
#' used for importing other GENBANK files but is \strong{not rigorously tested}.
#'
#' @param genbankFile Character: the path to the GENBANK file.
#'
#' @param type_vec Character: a vector of features to extract. The default value
#' has been set to extract the major elements of mitochondrial genomes, but
#' it could be edited to extract other features. Default is
#' \code{c('gene', 'CDS', 'tRNA', 'rRNA', 'D-loop', 'misc_feature')}.
#'
#' @return Returns a data.table with the following columns: \cr
#' \enumerate{
#'   \item \code{$NAME} = The name of the genetic feature.
#'   \item \code{$TYPE} = The type of genetic feature.
#'   \item \code{$STRAND} = Is given a 1 for the reported strand, and a -1
#'   if it is on the complementary strand.
#'   \item \code{$START} = The starting base position.
#'   \item \code{$END} = The end base position.
#' }
#'
#' @examples
#' library(genomalicious)
#'
#' # Create a link to raw external datasets in genomalicious
#' genomaliciousExtData <- paste0(find.package('genomalicious'), '/extdata')
#'
#' # Read in a GENBANK file
#' gbk.read <- mitoGbk2DT(paste(genomaliciousExtData, 'data_Bcocosensis_mitogenome.gbk', sep='/'))
#' gbk.read
#'
#' @export
mitoGbk2DT <- function(
  genbankFile,
  type_vec=c('gene', 'CDS', 'tRNA', 'rRNA', 'D-loop', 'misc_feature')){
  require(data.table)
  require(tidyverse)
  require(stringr)

  # Read in lines
  gbk <- readLines(genbankFile)

  # Start and end of features
  ft.st <- grep('FEATURES', gbk)
  ft.end <- grep('ORIGIN', gbk)

  # Subset the features out of GENBANK file
  ft.data <- gbk[ft.st:(ft.end-1)]

  # The positions of the elements
  ele.pos <- lapply(paste0('  ', type_vec, '  '), grep, x=ft.data) %>%
    unlist() %>% sort()

  # List of elements
  elementList <- list()

  # Iterature through each ith element number
  for(i in 1:length(ele.pos)){
    # Lines for the ith element
    if(i != length(ele.pos)){
      x <- ele.pos[i]
      y <- ele.pos[i+1]-1
    } else{
      x <- ele.pos[i]
      y <- length(ft.data)
    }

    # Subset feature data for ith element
    ft.sub <- ft.data[x:y]

    # The type
    sub.type <- type_vec[
      lapply(type_vec, grep, x=ft.sub[1]) %>%
        lapply(., length) %>%
        unlist() %>%
        as.logical()]

    # The nucleotide start and stop bases
    sub.st.end <- strsplit(ft.sub[1], sub.type)[[1]] %>%
      trimws() %>%
      .[2]

    # This removes a join(), or order(), when the element spans both
    # ends of the sequence.
    if(grepl('join', sub.st.end)){
      sub.st.end <- str_split(sub.st.end,',')[[1]] %>%
        gsub('join', '', .) %>%
        gsub('\\(', '', .) %>%
        gsub('\\)', '', .)
    } else if(grepl('order', sub.st.end)){
      sub.st.end <- str_split(sub.st.end,',')[[1]] %>%
        gsub('order', '', .) %>%
        gsub('\\(', '', .) %>%
        gsub('\\)', '', .)
    }

    # The strand, -1 is 'complement' is present
    sub.strand <- if_else(
      lapply(sub.st.end, function(x){
        grepl('complement', x=x)}) %>%
        unlist() %>% sum() == 0,
      1, -1)

    # The strand is -1, remove 'complement'.
    if(sub.strand==-1){
      sub.st.end <- gsub('complement', '', sub.st.end) %>%
        gsub('\\(', '', .) %>%
        gsub('\\)', '', .)
    }

    # Adjust the start and stop bases
    sub.st.end <- strsplit(sub.st.end, '\\.\\.') %>%
      lapply(., gsub, pattern='>', replacement='') %>%
      lapply(., as.integer)

    # The gene
    if(sub.type %in% c('tRNA','rRNA','misc_feature')){
      sub.name <- ft.sub[grep('product=', ft.sub)] %>%
        strsplit(., '/product=') %>%
        unlist() %>%
        trimws() %>%
        paste(., collapse='') %>%
        gsub('\"', '', .)
    } else if(sub.type %in% c('gene', 'CDS')){
      sub.name <- ft.sub[grep('gene=', ft.sub)] %>%
        strsplit(., '/gene=') %>%
        unlist() %>%
        trimws() %>%
        paste(., collapse='') %>%
        gsub('\"', '', .)
    } else if(sub.type=='D-loop'){
      sub.name <- 'D-loop'
    }

    # Add to list
    elementList[[i]] <- data.table(
      NAME=sub.name, TYPE=sub.type, STRAND=sub.strand,
      START=sub.st.end %>% lapply(., function(x) x[1]) %>% unlist(),
      END=sub.st.end %>% lapply(., function(x) x[2]) %>% unlist()
    )
  }

  # Output
  elementList %>% do.call('rbind',. ) %>% return()
}
