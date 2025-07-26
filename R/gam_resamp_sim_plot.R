#' Resampling-simulation framework for genomic animal models: Summary plot
#'
#' @description This is one of a set of functions that support use of a resampling-simulation
#' framework for assessing the robustness of genomic animal models. These functions
#' form a pipeline:
#' \enumerate{
#' \item \code{gam_resamp_sim_snps}
#' \item \code{gam_resamp_sim_phenos}
#' \item \code{gam_resamp_sim_stats}
#' \item \code{gam_resamp_sim_plot}
#' }
#' The \code{gam_resamp_sim_plot} function is used to generate a summary plot
#' of the observed and simulated estimates. You will need to
#' compile the estimates from your analyses into data.tables that will be fed
#' into this function. Estimates could be anything (heritability, additive variance,
#' residual environmental variance, model comparison P-values, etc).
#'
#' @param obsEstimDT Data.table: The observed estimates across SNP sets.
#' Requires the columns:
#' \enumerate{
#' \item \code{$SET}: The SNP set ID.
#' \item \code{$ESTIM}: The estimated value.
#' }
#'
#' @param simEstimDT Data.table: The simulated estimates across SNP sets. Note,
#' for each SNP set, you will need to take the median across nested runs.
#' Requires the columns:
#' \enumerate{
#' \item \code{$AMAT}: The A-matrix scenario.
#' \item \code{$SET}: The SNP set ID.
#' \item \code{$HERIT}: The simulated heritability level.
#' \item \code{$ESTIM}: The estimated value.
#' }
#'
#' @param xLabel Character: The label for the x-axis. Default is 'Estimate'.
#'
#' @param yLabel Character: The label for the y-axis. Default is 'Number of SNP sets'.
#'
#' @param AmtxColours Character: A vector of three colours for plotting the A-matrix
#' scenarios. Values are in this order: 'Observed', 'Simulated (without QTLs)',
#' 'Simulated (with QTLs)'. Default is NULL, so if unspecified, automatic colours
#' will be used.
#'
#' @param boxYpos Numeric: A vector of three values for positioning boxes for boxplots
#' below underneath histograms. These values must be negative. Values are in this
#' order: 'Observed', 'Simulated (without QTLs)', 'Simulated (with QTLs)'.
#' Default is NULL, so if unspecified, automatic values will be used.
#'
#' @param boxWith Numeric: A single value, the width of boxes for boxplots.
#' Default is 1.5.
#'
#' @returns Returns a ggplot object. For each heritability level simulated,
#' the observed and simulated estimates will be illustrated with histograms,
#' with boxplots summarising the distribution below them.
#'
#'@export
#'
gam_resamp_sim_plot <- function(
    obsEstimDT, simEstimDT, xLabel='Estimate', yLabel='Number of SNP sets',
    AmtxColours=NULL, boxYpos=NULL, boxWidth=1.5){
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   LIBRARIES AND ASSERTIONS   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  require(genomalicious)

  # Specify colours
  if(is.null(AmtxColours)){
    AmtxColours <- c('Observed'='#AD85E2','Simulated (without QTLs)'='#B5D84E','Simulated (with QTLs)'='#019E1E')
  } else{
    AmtxColours <- AmtxColours[1:3]
    names(AmtxColours) <- c('Observed','Simulated (without QTLs)','Simulated (with QTLs)')
  }

  # Specify y positions
  if(is.null(boxYpos)){
    boxYpos <- c('Observed'=-3,'Simulated (without QTLs)'=-7,'Simulated (with QTLs)'=-11)
  } else{
    boxYpos <- boxYpos[1:3]
    names(boxYpos) <- c('Observed','Simulated (without QTLs)','Simulated (with QTLs)')
  }

  # Observed data columns
  check.obs <- sum(c('SET','ESTIM') %in% colnames(obsEstimDT))

  if(check.obs<2){
    stop(
      'Argument `obsEstimDT` needs columns "SET" and "ESTIM". See ?gam_resamp_sim_plot.'
    )
  }

  # Simulated data columns
  check.sim <- sum(c('AMAT','SET','HERIT','ESTIM') %in% colnames(simEstimDT))

  if(check.sim<4){
    stop(
      'Argument `simEstimDT` needs columns "AMAT", "SET", "HERIT" and "ESTIM". See ?gam_resamp_sim_plot.'
    )
  }

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   SUMMARISE STATISTICS FOR BOXPLOTS   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  summStats <- gam_resamp_sim_stats(obsEstimDT = obsEstimDT, simEstimDT = simEstimDT)

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   MAKE THE PLOT   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # Internal theme
  my_gg_theme <- theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(size=10, face='bold'),
      legend.text = element_text(size=10),
      axis.text.y = element_text(angle=90, hjust=0.5),
      axis.ticks.length = unit(1.5,'mm'),
      legend.position = 'top',
      plot.title = element_text(hjust=0.5, size=10, face='bold')
    )

  # Prepare for plotting
  Amtx.levels <- c('Observed','Simulated (without QTLs)','Simulated (with QTLs)')

  obsEstimDT$AMAT <- 'Observed'
  obsEstimDT[, AMAT:=factor(AMAT, levels=Amtx.levels)]

  simEstimDT[, AMAT:=factor(AMAT, levels=Amtx.levels)]
  simEstimDT[, HERIT.PLOT:=paste0('bolditalic(h)^bold(2)~bold("= ', HERIT,' (sim)")')]

  obsStats <- summStats$obs %>%
    copy %>%
    .[, Y.POS:=boxYpos['Observed']] %>%
    setnames(., old=c('2.5%','25%','50%','75%','97.5%'), new=c('MIN','LOW','MID','HIGH','MAX'))
  simStats <- summStats$sim %>%
    copy %>%
    .[AMAT=='Simulated (without QTLs)', Y.POS:=boxYpos['Simulated (without QTLs)']] %>%
    .[AMAT=='Simulated (with QTLs)', Y.POS:=boxYpos['Simulated (with QTLs)']] %>%
    setnames(., old=c('2.5%','25%','50%','75%','97.5%'), new=c('MIN','LOW','MID','HIGH','MAX')) %>%
    .[, HERIT.PLOT:=paste0('bolditalic(h)^bold(2)~bold("= ', HERIT,' (sim)")')]

  # Plot
  ggplot() +
    my_gg_theme +
    # Histograms
    geom_histogram(
      data=obsEstimDT,
      mapping=aes(x=ESTIM, fill=AMAT, colour=AMAT),
      alpha = 0.6, position='identity'
    ) +
    geom_histogram(
      data=simEstimDT,
      mapping=aes(x=ESTIM, fill=AMAT, colour=AMAT),
      alpha = 0.6, position='identity'
    ) +
    # Whiskers: The 2.5% and 97.5% percentiles
    geom_segment(
      data=obsStats,
      mapping=aes(x=MIN, xend=MAX, y=Y.POS, yend=Y.POS, colour=AMAT)
    ) +
    geom_segment(
      data=simStats,
      mapping=aes(x=MIN, xend=MAX, y=Y.POS, yend=Y.POS, colour=AMAT)
    ) +
    # Box: 25% and 75%
    geom_rect(
      data=obsStats,
      mapping=aes(
        xmin=LOW, xmax=HIGH,
        ymin=Y.POS-boxWidth, ymax=Y.POS+boxWidth,
        fill=AMAT
      )
    ) +
    geom_rect(
      data=simStats,
      mapping=aes(
        xmin=LOW, xmax=HIGH,
        ymin=Y.POS-boxWidth, ymax=Y.POS+boxWidth,
        fill=AMAT
      )
    ) +
    # Median
    geom_errorbarh(
      data=obsStats,
      mapping=aes(xmin = MID, xmax=MID, y=Y.POS),
      stat = "identity", height = boxWidth*2, colour='grey20', linewidth=1
    ) +
    geom_errorbarh(
      data=simStats,
      mapping=aes(xmin = MID, xmax=MID, y=Y.POS),
      stat = "identity", height = boxWidth*2, colour='grey20', linewidth=1
    ) +
    # Plot visuals
    facet_wrap(~HERIT.PLOT, labeller=label_parsed) +
    scale_fill_manual(values=AmtxColours) +
    scale_colour_manual(values=AmtxColours) +
    labs(
      x=xLabel, y=yLabel,
      fill='Dataset',
      colour='Dataset'
    )
}
