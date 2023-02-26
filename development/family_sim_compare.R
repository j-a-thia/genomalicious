library(AGHmatrix)

simGenoMat <- simFamily %>%
  dcast(., SAMPLE~LOCUS, value.var='GT') %>%
  as.data.frame() %>%
  column_to_rownames(., 'SAMPLE') %>%
  as.matrix()

simGRM <- Gmatrix(simGenoMat, method='Yang', ploidy=2)

obsGenoMat <- data_Genos %>%
  .[POP=='Pop1'] %>%
  DT2Mat_genos(., sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT')

obsGRM <- Gmatrix(obsGenoMat, method='Yang', ploidy=2)




family_sim_compare <- function(
    simFamily, simGRM, obsGRM, look='ggplot', legendPos='right',
    curveAlpha=0.7, curveFill=NULL, curveOutline=NULL
){
  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  require(data.table); require(ggplot2); require(tidyverse)

  # Adjust values for curveOutline if NULL or if length != 5
  if(is.null(curveOutline)){
    curveOutline <- rep('grey20', 5)
  }

  if(length(curveOutline)<5){
    curveOutline <- rep(curveOutline[1], 5)
  }

  # Check that fill colours are specified correctly if NULL, and
  # if not NULL, make sure there are only 5 colours
  if(!is.null(curveFill)){
    if(length(curveFill)!=5){
      stop('Argument `curveFill` must be a character vector of 5 colours. See ?family_sim_plot.')
    }
  }

  if(is.null(curveFill)){
    curveFill <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3")
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

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Number of sims
  numSims <- simFamily$SIM %>% max

  # Compile the relatedness data table from simulated and observed GRMs.
  simRel <- lapply(1:numSims, function(sim){
    # Line up the known pairs. Note, that samples with "G3.2" are half-siblings
    # with "G3.3" and cousins with "G3.4". But these relationships are not
    # included to keep things balanced.
    unrel <- c(paste0('S',sim,c('_UR.1','_UR.2')))
    sibs <- c(paste0('S',sim,c('_G3.1','_G3.2')))
    halfsibs <- c(paste0('S',sim,c('_G3.1','_G3.3')))
    cousins <- c(paste0('S',sim,c('_G3.1','_G3.4')))

    data.table(
      SIM=sim,
      SAMPLE1=c(unrel[1],sibs[1],halfsibs[1],cousins[1]),
      SAMPLE2=c(unrel[2],sibs[2],halfsibs[2],cousins[2]),
      FAMILY=c('Unrelated','Siblings','Half-siblings','Cousins'),
      RELATE=c(
        simGRM[unrel[1],unrel[2]],
        simGRM[sibs[1],sibs[2]],
        simGRM[halfsibs[1],halfsibs[2]],
        simGRM[cousins[1],cousins[2]]
      )
    )
  }) %>%
    do.call('rbind', .)

  obsRel <- combn(colnames(obsGRM),2) %>%
    apply(., 2, function(x){
      data.table(SIM=NA, SAMPLE1=x[1], SAMPLE2=x[2], FAMILY='Observed', RELATE=obsGRM[x[1],x[2]])
      }) %>%
    do.call('rbind', .)

  # Make the data
  rel_data <- rbind(obsRel, simRel) %>%
    as.data.table %>%
    .[, FAMILY:=factor(FAMILY, levels=c('Unrelated','Cousins','Half-siblings','Siblings','Observed'))]

  # Make the plot
  rel_gg <- ggplot(rel_data, aes(x=RELATE, fill=FAMILY, colour=FAMILY)) +
    plotTheme +
    geom_density(alpha=curveAlpha,position="identity") +
    geom_vline(xintercept=c(0,0.125,0.25,0.5), linetype='longdash') +
    scale_colour_manual(values=curveOutline) +
    scale_fill_manual(values=curveFill) +
    scale_x_continuous(breaks=seq(-0.1, 1, 0.1)) +
    labs(x='Relatedness', y='Density', fill='', colour='')

  # Return a list
  list(data=rel_data, plot=rel_gg) %>% return()
}
