#' Resampling-simulation framework for genomic animal models: Summary statistics
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
#' The \code{gam_resamp_sim_stats} function is used to summarise observed and
#' simulated estimates and produce statistics for reporting. You will need to
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
#' @returns Returns a list with the following indexed slots:
#' \code{$obs}: The observed summary statistics, with columns:
#' \enumerate{
#' \item \code{$AMAT}: The A-matrix scenario ('Observed').
#' \item \code{$2.5%}: The 2.5% percentile across SNP sets.
#' \item \code{$25%}: The 25% percentile across SNP sets.
#' \item \code{$50%}: The 50% percentile across SNP sets.
#' \item \code{$75%}: The 75% percentile across SNP sets.
#' \item \code{$97.5%}: The 97.5% percentile across SNP sets.
#' }
#' \code{$sim}: The simulated summary statistics, with the columns:
#' \enumerate{
#' \item \code{$AMAT}: The A-matrix scenario ('with QTLs' and 'without QTLs').
#' \item \code{$HERIT}: The heritability levels simulated.
#' \item \code{$2.5%}: The 2.5% percentile across SNP sets.
#' \item \code{$25%}: The 25% percentile across SNP sets.
#' \item \code{$50%}: The 50% percentile across SNP sets.
#' \item \code{$75%}: The 75% percentile across SNP sets.
#' \item \code{$97.5%}: The 97.5% percentile across SNP sets.
#' }
#' \code{$comp}: The comparison between observed and simulated estimates, with the columns:
#' \enumerate{
#' \item \code{$AMAT}: The A-matrix scenario ('with QTLs' and 'without QTLs').
#' \item \code{$HERIT}: The heritability levels simulated.
#' \item \code{$Pct(S>Omed)}: Percent of simulated estimates > observed median.
#' \item \code{$Pct(S<Omed)}: Percent of simulated estimates < observed median.
#' \item \code{$Pct(Ovl[S,O])}: Percent of overlap between simulated and observed estimates.
#' \item \code{$75%}: The 75% percentile across SNP sets.
#' \item \code{$97.5%}: The 97.5% percentile across SNP sets.
#' }
#'
#'@export
#'
gam_resamp_sim_stats <- function(obsEstimDT, simEstimDT){
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   LIBRARIES AND ASSERTIONS   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  require(genomalicious); require(overlapping)

  # Observed data columns
  check.obs <- sum(c('SET','ESTIM') %in% colnames(obsEstimDT))

  if(check.obs<2){
    stop(
      'Argument `obsEstimDT` needs columns "SET" and "ESTIM". See ?gam_resamp_sim_stats.'
    )
  }

  # Simulated data columns
  check.sim <- sum(c('AMAT','SET','HERIT','ESTIM') %in% colnames(simEstimDT))

  if(check.sim<4){
    stop(
      'Argument `simEstimDT` needs columns "AMAT", "SET", "HERIT" and "ESTIM". See ?gam_resamp_sim_stats.'
    )
  }

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   LIBRARIES AND ASSERTIONS   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # Observed summary statistics
  obsStats <- obsEstimDT[, .(
    `2.5%`=quantile(ESTIM, 0.025),
    `25%`=quantile(ESTIM, 0.25),
    `50%`=quantile(ESTIM, 0.5),
    `75%`=quantile(ESTIM, 0.75),
    `97.5%`=quantile(ESTIM, 0.975)
  )] %>%
    data.table(AMAT='Observed', .)

  # Simulated summary statistics
  simStats <- simEstimDT[, .(
    `2.5%`=quantile(ESTIM, 0.025),
    `25%`=quantile(ESTIM, 0.25),
    `50%`=quantile(ESTIM, 0.5),
    `75%`=quantile(ESTIM, 0.75),
    `97.5%`=quantile(ESTIM, 0.975)
  ), by=c('AMAT','HERIT')]

  # Comparison summary statistics
  compStats <- left_join(
    simEstimDT %>%
      .[, .(`Pct(S>Omed)`=sum(ESTIM>(obsStats$`50%`))/length(ESTIM) * 100), by=c('AMAT','HERIT')],
    simEstimDT %>%
      .[, .(`Pct(S<Omed)`=sum(ESTIM<(obsStats$`50%`))/length(ESTIM) * 100), by=c('AMAT','HERIT')]
  ) %>%
    left_join(
      .,
      simEstimDT %>%
        .[, .(`Pct(Ovl[S,O])`=overlapping::overlap(list(obsEstimDT$ESTIM, ESTIM), type='2')$OV*100), by=c('AMAT','HERIT')] %>%
        .[, `Pct(Ovl[S,O])`:=round(`Pct(Ovl[S,O])`,2)]
    )

  # Output
  list(obs=obsStats, sim=simStats, comp=compStats)
}
