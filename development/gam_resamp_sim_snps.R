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
#' 1. \code{$LINK.BLOCK}, the linkage block ID.
#' 2. \code{$LOCUS}, the locus ID (biallelic SNPs).
#'
#' @param seqType Character: The type of sequencing data used. One of
#' \code{'wgs'} for whole-genome sequencing (linkage blocks are genomic windows),
#' or \code{'rrs'} for reduced-representation sequencing (linkage blocks are
#' contigs, or genomic fragments).
#'
#' @param numSnpSets Integer: The number of unique resampled SNP sets to draw.
#' One random SNP will be drawn per linkage block.
#'
#' @param dropBlocksRRS Numeric: The number of linkage blocks to drops for the
#' for reduced-representation sequencing datasets (\code{seqType=='rrs'}). If you
#' intend to use the full pipeline, then you need to specify this value to be >0.
#' These dropped linkage blocks will be used to simulate unsampled genomic regions
#' in the function, \code{gam_resamp_sim_phenos}. You should \code{dropBlocksRRS}
#' to be the maximum number of QTLs you would like to simulate, because the default
#' behaviour of \code{gam_resamp_sim_phenos} is to draw a maximum of one locus
#' per reserved linkage block.

gam_resamp_sim_snps <- function(
    uniqLocusDT, seqType, numSnpSets, dropBlocksRRS=NULL){

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
  if(seqType=='rrs' & is.null(dropBlocksRRS)){
    stop(
      'Argument `seqType`=="rrs", but the number of linkage blocks to drop per SNP set has not been specified by `dropBlocksRRS`. See ?gam_resamp_sim_snps.'
    )
  }

  if(dropBlocksRRS==0){
    warning(
      'Argument `dropBlocksRRS` has been set to value of 0. This will cause issues for downstream simulations. See See ?gam_resamp_sim_snps.'
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
        .[, sample(unique(LINK.BLOCK), dropBlocksRRS, FALSE)]

      D.loci <- uniqLocusDT %>%
        .[!LINK.BLOCK %in% link.block.reserve,] %>%
        .[,.(LOCUS=sample(LOCUS,1,FALSE)), by='LINK.BLOCK']

      list(focal=data.table(SET=i, D.loci), reserved=data.table(SET=i, LINK.BLOCK=link.block.reserve))
    })
  }

  return(snpSetList)
}
