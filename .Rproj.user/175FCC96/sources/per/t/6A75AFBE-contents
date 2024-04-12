#' Plot amino acid alignments
#'
#' This function takes in an alignment of amino acid sequences in a matrix format
#' and returns a plot of that alignment.
#'
#' @param alignMat Matrix: the DNA alignments. Rows are the taxa, columns are
#' the base positions, and cells contain the values. Values should be one of
#' the 20 amino acids in single letter form (e.g., 'A', 'V'), '?', '-'.
#'
#' @param aa_colours Character: named vector of colour values for nucleotides.
#' For example, \code{c('K'='blue', 'D'='red'...)}. You need a colour assigned
#' to each amino acid. Default is NA, which will produce automated colours based
#' on whether an amino acid is small non-polar, hydrophobic, polar, negatively
#' charged, or positively charged.
#'
#' @param other_colours Character: named vector colour values for '?'.
#' For example, \code{c('?'='black')}. Default is NA, which will produce
#' automated colours.
#'
#' @param gap_colour Chracter: named vector colour for gaps. For examples,
#' \code{c('-'='white')}. Default is NA, which will produce
#' automated colours.
#'
#' @param border_colour Character: single value, the colour of borders around
#' base positions. Default is NA, no assigned border colour.
#'
#' @param pos_vec Integer: a vector of base positions to plot, subsets alignment.
#' Default is NULL, which will plot all bases.
#'
#' @param show_residues Logical: whether the values of residues should be plotted.
#' Default is TRUE.
#'
#' @param residue_size Numeric: the text size used for amino acid residues. Default is 3.
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
#' fasta <- paste0(genomaliciousExtData, '/data_COI_aa.fasta')
#'
#' # Multi sequence alignmnet of demo COI data.
#' aln <- align_many_genes_aa(fasta, gene.names='COI')
#'
#' # Plot base positions from 200:210
#' align_plot_aa(as.matrix(aln$COI$align), pos_vec=200:210)
#'
#' # Default colours, no text, and borders
#' align_plot_aa(
#'    as.matrix(aln$COI$align),
#'    aa_colours=c(
#'     # Small non-polar
#'     'G'='#F8C641','A'='#F8C641','S'='#F8C641','T'='#F8C641',
#'     # Hydrophobic
#'     'C'='grey80','V'='grey80','I'='grey80','L'='grey80','P'='grey80','F'='grey80','Y'='grey80','M'='grey80','W'='grey80',
#'     # Polar
#'     'N'='#08c7e0','Q'='#08c7e0','H'='#08c7e0',
#'     # Negatively charged
#'     'D'='#ce0073','E'='#ce0073',
#'     # Positively charged
#'     'K'='#59A3FF','R'='#59A3FF'
#'     ),
#'    border_colour='grey20',
#'    pos_vec=200:210,
#'    show_residues=FALSE
#' )
#'
#' @export

align_plot_aa<- function(
    alignMat, aa_colours=NULL, other_colours=NULL, gap_colour=NA,
    border_colour=NA, pos_vec=NULL, show_residues=TRUE, residue_size=3, show_legend=FALSE){

  require(tidyverse)
  require(data.table)

  # Is input a matrix?
  if(!'matrix' %in% class(alignMat)){
    stop('Argument alignMat is not of class matrix. See ?align_plot_aa.')
  }

  if(is.null(aa_colours)){
    aa_colours <- c(
      # Small non-polar
      'G'='#F8C641',
      'A'='#F8C641',
      'S'='#F8C641',
      'T'='#F8C641',
      # Hydrophobic
      'C'='#59D863',
      'V'='#59D863',
      'I'='#59D863',
      'L'='#59D863',
      'P'='#59D863',
      'F'='#59D863',
      'Y'='#59D863',
      'M'='#59D863',
      'W'='#59D863',
      # Polar
      'N'='#DB8BFF',
      'Q'='#DB8BFF',
      'H'='#DB8BFF',
      # Negatively charged
      'D'='#FF5858',
      'E'='#FF5858',
      # Positively charged
      'K'='#59A3FF',
      'R'='#59A3FF'
    )
  }

  if(is.null(other_colours)){
    other_colours <- c('?'='black')
  }

  if(is.na(gap_colour)){
    gap_colour <- c('-'='white')
  }

  other_colours <- c('?'='black', 'N'='black')

  # Combine colours
  colour_map <- c(aa_colours, gap_colour, other_colours)

  # Convert matrix into data.table
  alignDT <- alignMat %>%
    t() %>%
    as.data.table() %>%
    .[,POS:=1:.N] %>%
    melt(., id.vars='POS', variable.name='TAXA', value.name='X') %>%
    .[, TAXA:=factor(TAXA,levels=rownames(alignMat))]

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
              + labs(x='Alignment position', y='Taxa', fill=NULL)
  )

  # If base values desired
  if(show_residues==TRUE){
    alignGG <- alignGG + geom_text(mapping=aes(label=X), size=residue_size)
  }

  # If the legend should be removed
  if(show_legend==FALSE){
    alignGG <- alignGG + theme(legend.position='none')
  }

  # Output
  return(alignGG)

}

