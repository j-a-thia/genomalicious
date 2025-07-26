#' Resampling-simulation framework for genomic animal models: SNP set
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
#' The \code{gam_resamp_sim_snps} function is used to generate a suite of
#' resampled SNP sets to explore how different compositions of genetic markers
#' influence the performance of genomics animal models.
#'
#' @param uniqLocusDT Data.table: Unique loci assigned to linkage blocks.
#' Requires the columns:
#' \enumerate{
#' \item \code{$LINK.BLOCK}: the linkage block ID.
#' \item \code{$LOCUS}, the locus ID (biallelic SNPs).
#'
#' @param seqType Character: The type of sequencing data used. One of
#' \code{'wgs'} for whole-genome sequencing (linkage blocks are genomic windows),
#' or \code{'rrs'} for reduced-representation sequencing (linkage blocks are
#' contigs, or genomic fragments).
#'
#' @param numSnpSets Integer: The number of unique resampled SNP sets to draw.
#' One random SNP will be drawn per linkage block.
#'
#' @param reserveBlocksRRS Numeric: The number of linkage blocks to reserve for the
#' for reduced-representation sequencing datasets (\code{seqType=='rrs'}). If you
#' intend to use the full pipeline, then you need to specify this value to be >0.
#' These reserved linkage blocks will be used to simulate unsampled genomic regions
#' in the function, \code{gam_resamp_sim_phenos}. You should \code{reserveBlocksRRS}
#' to be the maximum number of QTLs you would like to simulate, because the default
#' behaviour of \code{gam_resamp_sim_phenos} is to draw a maximum of one locus
#' per reserved linkage block.
#'
#' @returns Returns a list with \code{numSnpSets} number of slots. Each slot
#' contains the following indexed items:
#' \enumerate{
#' \item \code{$focal}: Data.table with the randomly sampled SNPs (per linkage
#'    block) that will be used in this SNP set, with columns \code{$SET},
#'    \code{$LINK.BLOCK}, and \code{$LOCUS}. This will be present for both
#'    \code{seqType == 'wgs'} and \code{seqType == 'rrs'}. For \code{seqType == 'rrs'},
#'    a subset of \code{reserveBlocksRRS} linkage blocks will be dropped randomly
#'    and reserved for simulation downstream.
#' \item \code{$reserved}: Data.table with the reserved linkage blocks for the
#'    SNP set when \code{seqType == 'rrs'}, with columns \code{$SET} and \code{$LINK.BLOCk}.
#' }
#'
#'@export

gam_resamp_sim_snps <- function(
    uniqLocusDT, seqType, numSnpSets, reserveBlocksRRS=NULL){

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   LIBRARIES AND ASSERTIONS   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  require(genomalicious)

  check.cols <- sum(c('LINK.BLOCK','LOCUS') %in% colnames(uniqLocusDT))

  # Columns
  if(check.cols<2){
    stop(
      'Arguments `uniqLocusDT` must have column "LINK.BLOCK" and "LOCUS". See ?gam_resamp_sim_snps.')
  }

  # Dropping linkage blocks for reduced-representation
  if(seqType=='rrs' & is.null(reserveBlocksRRS)){
    stop(
      'Argument `seqType`=="rrs", but the number of linkage blocks to drop per SNP set has not been specified by `reserveBlocksRRS`. See ?gam_resamp_sim_snps.'
    )
  }

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   MAKE SNP SETS   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # If a whole-genome sequencing dataset
  if(seqType=='wgs'){
    snpSetList <- lapply(1:numSnpSets, function(i){
      D.loci <- uniqLocusDT[, .(LOCUS=sample(LOCUS,1,FALSE)), by='LINK.BLOCK']
      list(focal=data.table(SET=i, D.loci), reserved=NULL)
    })
  # If a reduced representation sequencing dataset
  } else if(seqType=='rrs'){
    snpSetList <- lapply(1:numSnpSets, function(i){
      # Reserve linkage blocks for later QTL simulation
      link.block.reserve <- uniqLocusDT %>%
        .[, sample(unique(LINK.BLOCK), reserveBlocksRRS, FALSE)]

      D.loci <- uniqLocusDT %>%
        .[!LINK.BLOCK %in% link.block.reserve,] %>%
        .[,.(LOCUS=sample(LOCUS,1,FALSE)), by='LINK.BLOCK']

      list(focal=data.table(SET=i, D.loci), reserved=data.table(SET=i, LINK.BLOCK=link.block.reserve))
    })
  }

  return(snpSetList)
}
