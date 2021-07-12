#' Plot DNA alignments
#'
#' This function takes in an alignment of DNA sequences in a matrix format
#' and returns a plot of that alignment.
#'
#' @param alignMat Matrix: the DNA alignments. Rows are the taxa, columns are
#' the base positions, and cells contain the values. Values should be one of
#' 'A', 'T', 'G', 'C', '?', 'N', '-'.
#'
#' @param nuc_colours Character: named vector of colour values for nucleotides.
#' For example, \code{c('A'='blue', 'T'='yellow', 'G'='green', 'C'='red')}.
#' Default is NA, which will produce automated colours.
#'
#' @param other_colours Character: named vector colour values for '?' or 'N'
#' values. For example, \code{c('?'='black', 'N'='black')}.
#' Default is NA, which will produce automated colours.
#'
#' @param gap_colour Chracter: named vector colour for gaps. For examples,
#' \code{c('-'='white')}. #' Default is NA, which will produce
#' automated colours.
#'
#' @param border_colour Character: single value, the colour of borders around
#' base positions. Default is NA, no assigned border colour.
#'
#' @param pos_vec Integer: a vector of base positions to label plot.
#' Default is NA, which will produced automated values.
#'
#' @param show_bases Logical: whether the values of bases should be plotted.
#' Default is TRUE.
#'
#' @param show_legend Logical: whether a legend should be plotted.
#' Default is FALSE.
#'
#' @examples
#' # Create a link to raw external datasets in genomalicious
#' genomaliciousExtData <- paste0(find.package('genomalicious'), '/extdata')
#'
#' # Path to the demo FASTA file
#' fasta <- paste0(genomaliciousExtData, '/data_COI.fasta')
#'
#' # Multi sequence alignmnet of demo COI data.
#' co1 <- align_many_genes_dna(fasta, gene.names='COI')
#'
#' # Plot base positions from 200:250
#' align_plot_dna(as.matrix(co1$COI$align)[, 200:250])
#'
#' @export

align_plot_dna<- function(
  alignMat, nuc_colours=NA, other_colours=NA, gap_colour=NA,
  border_colour=NA, pos_vec=NA, show_bases=TRUE, show_legend=FALSE){

  require(tidyverse)
  require(data.table)

  # Is input a matrix?
  if(!'matrix' %in% class(alignMat)){
    stop('Argument alignMat is not of class matrix. See ?align_plot_dna.')
  }

  # If colours not specified
  if(is.na(nuc_colours)){
    nuc_colours <- c(
      'A'='royalblue2',
      'T'='gold',
      'G'='green3',
      'C'='firebrick2')
  }

  if(is.na(other_colours)){
    other_colours <- c('?'='black', 'N'='black')
  }

  if(is.na(gap_colour)){
    gap_colour <- c('-'='white')
  }

  # Combine colours
  colour_map <- c(nuc_colours, gap_colour, other_colours)

  # Convert matrix into data.table
  alignDT <- alignMat %>%
    t() %>%
    as.data.table() %>%
    .[,POS:=1:.N] %>%
    melt(., id.vars='POS', variable.name='TAXA', value.name='X')

  # Maximum base position
  pos.max <- alignDT$POS %>% max

  # If a vector of base positions is not specified
  if(is.na(pos_vec)){ pos_vec <- seq(1, pos.max, 5) }

  # Create alignment
  alignGG <- (ggplot(alignDT, aes(x=POS, y=TAXA,fill=X))
              + theme(
                panel.background = element_blank()
              )
              + geom_tile(
                colour=border_colour
              )
              + scale_fill_manual(values=colour_map)
              + scale_x_continuous(breaks=pos_vec)
              + labs(x='Base position', y='Taxa', fill=NULL)
  )

  # If base values desired
  if(show_bases==TRUE){
    alignGG <- alignGG + geom_text(mapping=aes(label=X))
  }

  # If the legend should be removed
  if(show_legend==FALSE){
    alignGG <- alignGG + guides(fill=FALSE)
  }

  # Output
  return(alignGG)

}

