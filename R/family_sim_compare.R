#' Compare simulated and observed inferred genetic relationships
#'
#' This function can help compare observed estimates of relatedness to expected
#' values of relatedness for different familial relationships given a set of
#' loci and population alelle frequencies. It returns a data table combining
#' the simulated and observed estiamtes and a plot of density curves overlaying
#' the distribution of observed relatedness estimates on top of simulated values.
#'
#' The function takes the output of \code{family_sim_genos} and additional calculations
#' of the genetic relationiship matrix (GRM) that are performed by the user on
#' the simulated and observed individuals.
#'
#' @param simGRM Matrix: The simulated GRM using the individuals generated
#' from \code{family_sim_genos}. See Details.
#'
#' @param obsGRM Matrix: The observed GRM using the observed individuals.
#' See Details.
#'
#' @param numSims Integer: The number of simulated individuals for each family
#' relationship. Default is 100.
#'
#' @param look Character: The look of the plot. Default = \code{'ggplot'}, the
#' typical gray background with gridlines produced by \code{ggplot2}. Alternatively,
#' when set to \code{'classic'}, produces a base R style plot.
#'
#' @param legendPos Character: Where should the legend be positioned? Default is
#' \code{'top'}, but could also be one of, \code{'right'}, \code{'bottom'},
#' \code{'left'}, or \code{'none'}.
#'
#' @param curveAlpha Numeric: A value between 0 and 1 to set the transparency of
#' the density curves. Default = 0.7.
#'
#' @param curveFill Character: A vector of colours to fill density curves,
#' but is an optional argument. Default = \code{NULL}. If specified, must be
#' a length of 5, with colours corresponding to the levels 'Unrelated', 'Cousins',
#' 'Half-siblings', 'Siblings', and 'Observed', in that order.
#'
#' @param curveOutline Character: A vector of colours to for density curve outlines,
#' but is an optional argument. Default = \code{NULL}. If specified, must be
#' a length of 5, with colours corresponding to the levels 'Unrelated', 'Cousins',
#' 'Half-siblings', 'Siblings', and 'Observed', in that order.
#'
#' @param facetFamily Logical: Should the plot be faceted by family relationship?
#' Default is FALSE.
#'
#' @details The GRMs for arguments \code{simGRM} and \code{obsGRM} need to be
#' created by the user with whatever program they want to use to calculate
#' pairwise relatedness among individuals. The same function call should be used
#' ob both datasets. The GRM should be a square matrix with the relatedness of
#' individuals to themselves on the diagonal, and their relatedness to other
#' individuals on the off-diagonal.
#'
#' @returns Returns two objects. The first is data.table with the following columns:
#' \enumerate{
#'    \item \code{$SIM}, the simulation number for pairs of simulated individuals,
#'    or 'NA' for the pairs of observed individuals.
#'    \item \code{$SAMPLE1}, the sample ID for the first individual.
#'    \item \code{$SAMPLE2}, the sample ID for the second individual.
#'    \item \code{$FAMILY}, the familial relationship for simulated individuals.
#'    \item \code{$RELATE}, the estimated relatedness.
#' }
#'
#' The second is a ggplot object which plots density curves for the estimated
#' relatedness values calculated for simulated pairs of unrelated individuals,
#' cousins, half-siblings, and siblings, with the observed relatedness values
#' from the user's dataset overlayed. Dashed lines are used to demarcate the
#' theoretical expected relatedness values for unrelated individuals (0),
#' cousins (0.125), half-siblings (0.25), and siblings (0.5).
#'
#' @examples
#' library(genomalicious)
#' data(data_Genos)
#'
#' # Subset Pop1 genotypes
#' genosPop1 <- data_Genos[POP=='Pop1', c('SAMPLE', 'LOCUS', 'GT')]
#'
#' # Get the allele frequencies for Pop1
#' freqsPop1 <- genosPop1[, .(FREQ=sum(GT)/(length(GT)*2)), by=LOCUS]
#'
#' # Simulate 100 families
#' simFamily <- family_sim_genos(
#'    freqData=freqsPop1,
#'    locusCol='LOCUS',
#'    freqCol='FREQ',
#'    numSims=100
#' )
#'
#' ### THE OBSERVED GENOMIC RELATIONSHIPS MATRIX
#' library(sommer)
#'
#' # Note, for sommer, we have to adjust genotyeps to range from -1 to 1.
#'
#' # A genotype matrix for the focal pairs
#' obsGenosMat <- genosPop1 %>%
#'   copy %>%
#'   .[, GT:=GT-1] %>%
#'   DT2Mat_genos()
#'
#' # Calculate the GRM
#' obsGRM <- sommer::A.mat(obsGenosMat)
#'
#' ### THE SIMULATED GENOMIC RELATIONSHIPS MATRIX
#' # Convert simulated families into a genotype matrix
#' simGenosMat <- simFamily$focal.pairs %>%
#'   copy %>%
#'   .[, GT:=GT-1] %>%
#'   DT2Mat_genos()
#'
#' # Calculate the GRM
#' simGRM <- sommer::A.mat(simGenosMat)
#'
#' ### COMPARE THE OBSERVED AND SIMULATED
#' relComp <- family_sim_compare(
#'    simGRM=simGRM,
#'    obsGRM=obsGRM,
#'    look='classic'
#' )
#'
#' # The data
#' relComp$data
#'
#' # Simulated dataset
#' relComp$data[!is.na(SIM)]
#'
#' # The observed dataset
#' relComp$data[is.na(SIM)]
#'
#' # Plot of relatedness values. Dashed lines denote relatedness
#' # values of 0, 0.0625, 0.125, 0.25, and 0.5, which are the theoretical
#' # expectations for unrelated individuals, half-cousins, cousins, half-siblings,
#' # and siblings, respectively.
#' # You will note a large variance are the expected values, which
#' # is not surprising for this very small SNP dataset (200 loci).
#' relComp$plot
#'
#' # Take a look at the "known" relationships in the observed dataset using the
#' # offspring we created.
#' # Note, whilst values are centered on their theoretical expectations,
#' # you will probably find the "observed" relatedness might be lower or even higher.
#' relComp$data[FAMILY=='Half-siblings']$RELATE %>% summary()
#' relComp$data[FAMILY=='Siblings']$RELATE %>% summary()
#' relComp$data[FAMILY=='Parent-offspring']$RELATE %>% summary()
#'
#' @export

family_sim_compare <- function(
    simGRM, obsGRM, plotColous=NULL, look='ggplot', legendPos='right',
    curveAlpha=0.7, curveFill=NULL, curveOutline=NULL, facetFamily=FALSE
){
  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  require(data.table); require(ggplot2); require(tidyverse)

  # Default colours
  default.colours <- c("#DE3232", "#FF9510", "#E8CE00", "#54C567", "#02CEF2", "#4F74E6","#DA5EE4")

  # Adjust values for curveOutline if NULL or if length != 7
  if(is.null(curveOutline)){
    curveOutline <- default.colours
  }

  if(length(curveOutline)<7){
    curveOutline <- rep('grey20', 7)
  }

  # Check that fill colours are specified correctly if NULL, and
  # if not NULL, make sure there are only 7 colours
  if(!is.null(curveFill)){
    if(length(curveFill)!=7){
      stop('Argument `curveFill` must be a character vector of 7 colours. See ?family_sim_plot.')
    }
  }

  if(is.null(curveFill)){
    curveFill <- default.colours
  }

  # Legend position
  if(!legendPos %in% c('top','right','left','bottom')){
    stop('Argument `legendPos` must be on of "top", "right", "left", and "bottom". See ?family_sim_plot.')
  }

  # Set the plot theme by look
  if(!look %in% c('ggplot','classic')){
    stop('Argument `look` must be one of "ggplot" or "classic". See ?family_sim_plot.')
  }

  if(look=='ggplot'){
    plotTheme <- theme_gray() + theme(legend.position=legendPos, axis.ticks.length = unit(0.2, 'cm'))
  } else if(look=='classic'){
    plotTheme <- theme_bw() + theme(
      panel.grid.major=element_blank()
      , panel.grid.minor=element_blank()
      , text=element_text(colour='black')
      , legend.position=legendPos
      , axis.ticks.length=unit(0.2, 'cm'))
  }

  # Make sure facetFamily is logical
  if(class(facetFamily)!='logical'){
    stop('Argument `facetFamily` must be a logical object. See ?family_sim_plot.')
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+

  # Simulated relatedness values in long-format
  simRel <- simGRM %>%
    pairwiseMat2DT(., X1='SAMPLE1', X2='SAMPLE2', Y='RELATE') %>%
    .[SAMPLE1!=SAMPLE2]

  simRel[, SIM1:=sub('_.*', '', SAMPLE1)]
  simRel[, SIM2:=sub('_.*', '', SAMPLE2)]

  simRel[SIM1!=SIM2, FAMILY:='Unrelated']

  simRel[
    SIM1 == SIM2 &
      grepl('_unrel_',SAMPLE1)==TRUE &
      grepl('_unrel_',SAMPLE2)==TRUE,
    FAMILY:='Unrelated'
  ]

  simRel[
    SIM1 == SIM2 &
      grepl('_po_',SAMPLE1)==TRUE &
      grepl('_po_',SAMPLE2)==TRUE,
    FAMILY:='Parent-offspring'
  ]

  simRel[
    SIM1 == SIM2 &
      grepl('_sib_',SAMPLE1)==TRUE &
      grepl('_sib_',SAMPLE2)==TRUE,
    FAMILY:='Siblings'
  ]

  simRel[
    SIM1 == SIM2 &
      grepl('_hsib_',SAMPLE1)==TRUE &
      grepl('_hsib_',SAMPLE2)==TRUE,
    FAMILY:='Half-siblings'
  ]

  simRel[
    SIM1 == SIM2 &
      grepl('_cuz_',SAMPLE1)==TRUE &
      grepl('_cuz_',SAMPLE2)==TRUE,
    FAMILY:='Cousins'
  ]

  simRel[
    SIM1 == SIM2 &
      grepl('_hcuz_',SAMPLE1)==TRUE &
      grepl('_hcuz_',SAMPLE2)==TRUE,
    FAMILY:='Half-cousins'
  ]

  simRel[SIM1 == SIM2 & is.na(FAMILY), FAMILY:='Unrelated']

  simRel <- simRel[, SIM:=paste0(SIM1, '|', SIM2)] %>%
    .[, !c('SIM1','SIM2')]

  # Observed relatedness values in long-format.
  obsRel <- obsGRM %>%
    pairwiseMat2DT(., X1='SAMPLE1', X2='SAMPLE2', Y='RELATE') %>%
    .[SAMPLE1!=SAMPLE2] %>%
    data.table(SIM=NA, . ,FAMILY='Observed')

  # Make the data
  rel_data <- rbind(obsRel, simRel) %>%
    as.data.table %>%
    .[, FAMILY:=factor(FAMILY, levels=c('Observed','Unrelated','Half-cousins','Cousins','Half-siblings','Siblings','Parent-offspring'))]

  # Make the plot
  min.val.x <- round(min(rel_data$RELATE)/0.1) * 0.1

  rel_gg <- ggplot(rel_data, aes(x=RELATE, fill=FAMILY, colour=FAMILY)) +
    plotTheme +
    geom_density(alpha=curveAlpha,position="identity") +
    geom_vline(xintercept=c(0,0.0625,0.125,0.25,0.5), linetype='longdash') +
    scale_colour_manual(values=curveOutline) +
    scale_fill_manual(values=curveFill) +
    scale_x_continuous(breaks=seq(min.val.x, 1, 0.1)) +
    labs(x='Relatedness', y='Density', fill='', colour='')

  # Facet the plot?
  if(facetFamily==TRUE){
  rel_gg <- rel_gg + facet_wrap(~FAMILY, ncol=1, nrow=7)
  }

  # Return a list
  list(data=rel_data, plot=rel_gg) %>% return()
}

