#' Plot a gene map
#'
#' Uses base positions information to make a linear plot of genetic features
#' to produce a gene map. Separates into "gene features" (plotted as large blocks
#' on the main chromosome) and "extra features" (plotted as small bars offset from
#' the main chromosome).
#'
#' @param mapDT Data.table: genetic feature information. Requires the columns: \cr
#' \enumerate{
#'    \item $NAME = Character, the name of the genetic feature.
#'    \item $TYPE = Character, the type of genetic feature.
#'    \item $STRAND = Integer, the starnd, either 1 or -1.
#'    \item $START = Integer, the starting base position.
#'    \item $END = Integer, the ending base position.
#' }
#'
#' @param genome_len Integer: the genome length. Default = NULL. If unspecified,
#' will be assigned the final base pair of the last genetic feature in \code{mapDT}.
#'
#' @param gene_colour Character: a vector of colours to plot genes. Each item
#' is a colour, with the gene accessible through \code{names(gene_colour)}.
#' See Details for parameterisation.
#'
#' @param gene_type Character: a vector of values present in \code{mapDT$TYPE}
#' that will be plotted as large coloured bars. For plotting purposes,
#' "gene features". Default = \code{c('gene', 'rRNA')}.
#' See Details for parameterisation.
#'
#' @param extra_type Character: a vector of values present in \code{mapDT$TYPE}
#' that will be plotted as small grey bars. For plotting purposes, "extra features".
#' Default = \code{c('tRNA', 'D-loop')}. See Details for parameterisation.
#'
#' @param plot_xmax Numeric: a single value, the maximum x-axis limit.
#' Default is the genome length, as per \code{genome_len}.
#'
#' @param extra_ypos Numeric: a single value, the starting y-axis position for
#' extra features. See Details for parameterisation.
#'
#' @param plot_ymax Numeric: a single value, the maximum y-axis limit.
#' See Details for parameterisation.
#'
#' @param gene_txt_size Integer: a single value, the size for gene feature labels.
#' Default is 4.
#'
#' @param extra_txt_size Integer: a single value, the size for extra feature labels.
#' Default is 4.
#'
#' @param font Character: a single value, the font family to use.
#' Default is \code{'Arial'}.
#'
#' @param gene_border Character: a single value, the colour for borders around
#' gene features. Default is NA, no border.
#'
#' @param gene_size Numeric: a single value, the thickness of borders around
#' gene features, if a colour if specified in \code{gene_border}. Default is 1.
#'
#' @details There are two major features plotted, "gene features" and
#' "extra features". These names are just for convention: gene features are
#' plotted as large coloured bars in center of the plot on the main "chromosome",
#' whereas extra features are plotted as small grey bars above/below the
#' gene features, offset from the main chromosome. Anything could be plotted as
#' a gene or extra feature, and these are specified through \code{gene_type}
#' and \code{extra_type}.
#'
#' The name of the genetic feature being plotted is the value of \code{mapDT$NAME}.
#' This value is effectively evaluated as a mathematical expression to allow
#' italics for gene names and mixed formatting in gene names. The internal function
#' call is the evaluation of values by \code{geom_text(..., parse=TRUE)} and
#' \code{geom_text_repel(..., parse=TRUE)}.
#'
#' The value of \code{mapDT$STRAND} dictates the position of the coloured bars.
#' A value of 1 places "genes" on the top of the genomic strand, whereas a value
#' of -1 places "genes" below the genomic strand.
#'
#' The colour of the gene features is specified through \code{gene_colour} as
#' a named vector. If there are two genes, 'COX1' and 'COX2', specification of
#' their colours can be done like so: \code{c(COX1='pink', COX2='blue')}.
#' If colours are not specified, one colour is automatically assigned to each
#' unique "gene".
#'
#' The value of \code{extra_ypos} specifies that distance of the extra features
#' from the gene features. Set larger if things are looking squashed.
#' Additionally, \code{plot_ymax} sets the maximal plotting area, so set this
#' value larger if things are not fitting well.
#'
#' @return Returns a gg object.
#'
#' @examples
#' library(genomalicious)
#'
#' # Create a link to raw external datasets in genomalicious
#' genomaliciousExtData <- paste0(find.package('genomalicious'), '/extdata')
#'
#' # Read in a GENBANK file of the Bathygobius cocosensis mitogenome
#' gbk.read <- mitoGbk2DT(paste(genomaliciousExtData, 'data_Bcocosensis.gbk', sep='/'))
#' head(gbk.read)
#'
#' # Subset out the "CDS" types and plot genes, rRNA, tRNA, and D-loop.
#' # Rename rRNAs for nicer plotting. Because $NAME is evaluated by the
#' # expression() function, it is useful to put single quotations around characters
#' to have them read as characters internally by gene_map_plot().
#' gbk.read[TYPE!='CDS'] %>%
#' .[NAME=='12S ribosomal RNA', NAME:='12S'] %>%
#' .[NAME=='16S ribosomal RNA', NAME:='16S'] %>%
#' .[, NAME:=paste0("'", NAME, "'")] %>%
#' gene_map_plot(mapDT=., genome_len=16692, extra_txt_size=3)
#'
#' # Plot just the COX genes and the D-loop as "gene features" with
#' # custom colours and a border. Again, not the use of single quotes nested in
#' double quotes, which will match up to the edited gene $NAME column below.
#' gene.col.vec <- c(
#' "'COX1'"='royalblue',
#' "'COX2'"='firebrick3',
#' "'COX3'"='mediumpurple2',
#' "'CYTB'"='plum3',
#' "'D-loop'"='grey40')
#'
#' # Subset focal genes, add quotes to ensure characters are parsed as characters.
#' gbk.read[NAME %in% c('COX1','COX2','COX3','CYTB','D-loop')] %>%
#' .[, NAME:=paste0("'", NAME, "'")] %>%
#'   gene_map_plot(
#'     mapDT=., genome_len=16692,
#'     gene_type=c('gene', 'D-loop'), gene_colour=gene.col.vec,
#'     extra_type=NULL, gene_border='black')
#'
#' # It is possible to parse characters without the double quotes, but note how
#' # the '-' character in 'D-loop' has been parsed as a minus symbol.
#' gbk.read[NAME %in% c('COX1','COX2','COX3','CYTB','D-loop')] %>%
#'   gene_map_plot(
#'     mapDT=., genome_len=16692, gene_type=c('gene', 'D-loop'),
#'     extra_type=NULL, gene_border='black')
#'
#' @export

gene_map_plot <- function(
  mapDT, genome_len=NULL, gene_colour=NULL,
  gene_type=c('gene', 'rRNA'), extra_type=c('tRNA', 'D-loop'),
  plot_xmax=genome_len, extra_ypos=3, plot_ymax=5,
  gene_txt_size=4, extra_txt_size=4, font='Arial',
  gene_border=NA, gene_size=2
){

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   ENVIRONMENT AND CHECKS   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  require(data.table); require(tidyverse); require(ggrepel)

  # Check position variables
  if(length(extra_ypos)!=1){
    stop('Argument extra_ypos must be length 1. See ?gene_map_plot.')
  }
  if(length(plot_ymax)!=1){
    stop('Argument plot_ymax must be length 1. See ?gene_map_plot.')
  }

  # Internal function
  FUN_gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  # Subset just the required columns
  mapDT <- mapDT %>%
    as.data.table %>%
    .[,c('START', 'END', 'TYPE','NAME', 'STRAND')]

  # Get the plotted genes data
  genesDT <- mapDT[TYPE %in% gene_type] %>%
    .[, Y.MAX:=if_else(STRAND==1, 1, 0)] %>%
    .[, Y.MIN:=if_else(STRAND==1, 0, -1)] %>%
    .[, X.MID:=(START+END)/2]

  # Get the extra plotted items
  extraDT <- mapDT[TYPE %in% extra_type] %>%
    .[, Y.MAX:=if_else(STRAND==1, 1, 0)] %>%
    .[, Y.MIN:=if_else(STRAND==1, 0, -1)] %>%
    .[, X.MID:=(START+END)/2]

  # Set gene colours if not specified
  if(is.null(gene_colour)){
    uniq_genes <- genesDT$NAME %>% unique() %>% sort()
    gene_colour <- FUN_gg_color_hue(length(uniq_genes))
    names(gene_colour) <- uniq_genes
  } else{
    if(FALSE %in% c(genesDT$NAME %in% names(gene_colour))){
      stop('Argument gene_colour has been manually specified, but all the genes
      to be plot are not accessible in names(gene_colour). See ?gene_map_plot.')
    }
  }

  # Set genome length if not specified
  if(is.null(genome_len)){
    genome_len <- mapDT$END %>% max()
  }

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####    PLOT GENE FEATURES   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # Base plot
  ggMito <- (
    ggplot()
    # Theme
    + theme(
      axis.line.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.title.y=element_blank(),
      axis.ticks.length.x=unit(2, 'mm'),
      axis.title.x=element_blank(),
      legend.position='none',
      panel.grid.minor.y=element_blank(),
      panel.grid.major.y=element_blank(),
      text=element_text(family=font)
    )
    # Gene features
    + geom_rect(
      data=genesDT,
      mapping=aes(xmin=START, xmax=END, ymin=Y.MIN, ymax=Y.MAX, fill=NAME),
      colour=gene_border, size=gene_size
    )
    # The central genome line
    + geom_rect(
      data=data.table(START=0, END=genome_len, Y.MIN=-0.1, Y.MAX=0.1),
      mapping=aes(xmin=START, xmax=END, ymin=Y.MIN, ymax=Y.MAX),
    )
    # Colours
    + scale_colour_manual(values=gene_colour)
    + scale_fill_manual(values=gene_colour)
    # Axis limits
    + xlim(0, plot_xmax)
    + ylim(-plot_ymax,plot_ymax)
  )
  # Text for genes features on positive and negative strands
  if(nrow(genesDT[STRAND==1]) > 0){
    ggMito <- (
      ggMito
      + geom_text(
        data=genesDT[STRAND==1],
        mapping=aes(x=X.MID, y=1.2, label=NAME, colour=NAME),
        angle=45, hjust='left', size=gene_txt_size, family=font, parse=TRUE
      )
    )
  }
  if(nrow(genesDT[STRAND==-1]) > 0){
    ggMito <- (
      ggMito
      + geom_text(
        data=genesDT[STRAND==-1],
        mapping=aes(x=X.MID, y=-1.2, label=NAME, colour=NAME),
        angle=45, hjust='right', size=gene_txt_size, family=font, parse=TRUE
      )
    )
  }

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   PLOT EXTRA FEATURES   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(nrow(extraDT)>0){
    # Subset out extra features on positive vs negative strand
    extraPos <- extraDT[STRAND==1]
    extraNeg <- extraDT[STRAND==-1]

    # Vector of colours: alternating between successive features
    extra.colour <- c('grey20', 'grey50')
    extra.pos.colour <- extra.colour[1:nrow(extraPos)%%2 + 1]
    extra.neg.colour <- extra.colour[1:nrow(extraNeg)%%2 + 1]

    # Add in features and their labels
    ggMito <- (
      ggMito
      + geom_rect(
        data=extraPos,
        mapping=aes(xmin=START, xmax=END, ymin=extra_ypos, ymax=extra_ypos+0.3),
        fill=extra.pos.colour, colour=NA
      )
      + geom_rect(
        data=extraNeg,
        mapping=aes(xmin=START, xmax=END, ymin=-extra_ypos, ymax=-extra_ypos-0.3),
        fill=extra.neg.colour, colour=NA
      )
      + geom_text_repel(
        data=extraPos,
        mapping=aes(x=X.MID, y=extra_ypos+0.3, label=NAME),
        angle=45, min.segment.length = unit(0, 'lines'),
        ylim=c(extra_ypos+0.5,plot_ymax), nudge_y=extra_ypos+0.5,
        size=extra_txt_size, family=font, parse=TRUE, segment.size=0.25
      )
      + geom_text_repel(
        data=extraNeg,
        mapping=aes(x=X.MID, y=-extra_ypos-0.3, label=NAME),
        angle=45, min.segment.length = unit(0, 'lines'),
        ylim=c(-plot_ymax,-extra_ypos-0.5), nudge_y=-extra_ypos-0.5,
        size=extra_txt_size, family=font, parse=TRUE, segment.size=0.25
      )
    )
  }

  # >>>>>>>>>>>>>>>>>>
  ####   OUTPUT   ####
  # >>>>>>>>>>>>>>>>>>
  return(ggMito)
}
