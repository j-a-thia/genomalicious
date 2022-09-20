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
#' @param pos_vec Integer: a vector of base positions to plot, subsets alignment.
#' Default is NULL, which will plot all bases.
#'
#' @param show_bases Logical: whether the values of bases should be plotted.
#' Default is TRUE.
#'
#' @param base_size Numeric: the text size used for nucleotide bases. Default is 3.
#'
#' @param show_legend Logical: whether a legend should be plotted.
#' Default is FALSE.
#'
#' @examples
#' library(genomalicious)
#'
#' # Create a link to raw external datasets in genomalicious
#' genomaliciousExtData <- paste0(find.package('genomalicious'), '/extdata')
#'
#' # Path to the demo FASTA file
#' fasta <- paste0(genomaliciousExtData, '/data_COI.fasta')
#'
#' # Multi sequence alignmnet of demo COI data.
#' aln <- align_many_genes_dna(fasta, gene.names='COI')
#'
#' # Plot base positions from 200:250
#' align_plot_dna(as.matrix(aln$COI$align), pos_vec=200:250)
#'
#' # Custom colours, no text
#' align_plot_dna(
#'    as.matrix(aln$COI$align),
#'    nuc_colours=c(`A`='#ce0073',`T`='#e46adf',`G`='royalblue',`C`='#08c7e0'),
#'    pos_vec=200:250,
#'    show_bases=FALSE
#' )
#'
#' @export

align_plot_dna<- function(
  alignMat, nuc_colours=NULL, other_colours=NULL, gap_colour=NA,
  border_colour=NA, pos_vec=NULL, show_bases=TRUE, base_size=3, show_legend=FALSE){

  require(tidyverse)
  require(data.table)

  # Is input a matrix?
  if(!'matrix' %in% class(alignMat)){
    stop('Argument alignMat is not of class matrix. See ?align_plot_dna.')
  }

  # If colours not specified
  if(is.null(nuc_colours)){
    nuc_colours <- c(
      'A'='royalblue2',
      'T'='gold',
      'G'='green3',
      'C'='firebrick2')
  }

  if(is.null(other_colours)){
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

  # If pos_vec specified, subset
  if(!is.null(pos_vec)){
    alignDT <- alignDT[POS %in% pos_vec]
  }

  # Maximum base position
  pos.max <- alignDT$POS %>% max

  # Create alignment
  alignGG <- (ggplot(alignDT, aes(x=POS, y=TAXA,fill=X))
              + theme(
                panel.background = element_blank(),
                axis.ticks.length= unit(2,'mm')
              )
              + geom_tile(
                colour=border_colour
              )
              + scale_fill_manual(values=colour_map)
              + scale_x_continuous(
                breaks=~round(unique(pretty(.))),
                expand=c(0,0)
                )
              + scale_y_discrete(
                expand=c(0,0)
              )
              + labs(x='Base position', y='Taxa', fill=NULL)
  )

  # If base values desired
  if(show_bases==TRUE){
    alignGG <- alignGG + geom_text(mapping=aes(label=X), size=base_size)
  }

  # If the legend should be removed
  if(show_legend==FALSE){
    alignGG <- alignGG + theme(legend.position='none')
  }

  # Output
  return(alignGG)

}

