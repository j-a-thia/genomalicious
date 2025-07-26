#' Resampling-simulation framework for genomic animal models: Phenotypes
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
#' The \code{gam_resamp_sim_phenos} function is used to generate a suite of
#' phenotypes for different SNP sets to explore how different compositions of
#' genetic markers influence the performance of genomics animal models.
#'
#' @param snpSetList List: The output from the function \code{gam_resamp_sim_snps}.
#'
#' @param genosDT Data.table: The observed genotype data for each individual.
#' Requires the columns:
#' \enumerate{
#' \item \code{$LINK.BLOCK}: The linkage block ID.
#' \item \code{$LOCUS}: The locus ID.
#' \item \code{$SAMPLE}: The sample ID.
#' \item \code{$GT}: The genotype score (0, 1, or 2)
#' }
#'
#' @param runsPerSet Integer: The number of simulation runs to perform for each
#' SNP set. Recommend a minimum of 3-5. Each simulation will randomly draw a
#' `qtlNumber` of loci to act as QTLs. So you can think of this argument as
#' determining the number of resampled nested QTL sets to simulate within each
#' SNP set.
#'
#' @param seqType Character: The type of sequencing data used. One of
#' \code{'wgs'} for whole-genome sequencing (linkage blocks are genomic windows),
#' or \code{'rrs'} for reduced-representation sequencing (linkage blocks are
#' contigs, or genomic fragments).
#'
#' @param qtlNumber Integer: The number of loci to draw as simulated QTLs.
#'
#' @param h2Levels Numeric: A vector of heritability levels to simulate. Must be
#' between 0 and 1.
#'
#' @param effectSizes Numeric: A vector of QTL effect sizes to draw from in
#' simulations. All values must be positive. If you want to simulate a uniform
#' effect size, simply specify a single value of 1. However, you can use a distribution
#' of values if you want to simulate variable effect sizes across QTLs.
#'
#' @param numCores Integer: Number of cores to run simulations in parallel.
#'
#' @details This function is effectively calling \code{family_sim_qtl} under the
#' hood. If you want more information about that function, see \code{?family_sim_qtl}.
#'
#' Two scenarios are simulated with respect to whether QTLs will be present
#' in the calculated A-matrix (genomic relationship matrix): 'with QTLs' and 'without QTLs'.
#'
#' In the 'with QTLs' scenario, loci in the focal SNP set are randomly drawn to
#' use as simulated QTLs. This represents a 'best case' context where the focal
#' SNP set captures all the loci contributing to a traits additive genetic variance.
#'
#' In the contrasting scenario, 'without QTLs', loci that are not in the focal
#' SNP set are randomly drawn to simulate as QTLs. For \code{seqType=='wgs'},
#' these can be any loci found throughout the genome that is not in the focal SNP set.
#' For \code{seqType=='rrs'}, loci in the reserved contigs are are used, and one
#' locus is drawn per contig. Therefore, if \code{seqType=='rrs'}, you must make
#' sure that \code{qtlNumber} is the same as the number of reserved SNP sets.
#'
#' @returns Returns a list with one slot for each SNP set. Each slot contains
#' a data.table with the following columns:
#' \enumerate{
#' \item \code{$AMAT}: The A-matrix scenario.
#' \item \code{$SET}: The SNP set ID.
#' \item \code{$HERIT}: The simulated heritability level.
#' \item \code{$RUN}: The simulated run ID.
#' \item \code{$SAMPLE}: The sample ID.
#' \item \code{$G}: The genetic value.
#' \item \code{$E}: The residual environmental value.
#' \item \code{$P}: The phenotypic value (G+E).
#' }
#'
#'@export
#'
gam_resamp_sim_phenos <- function(
    snpSetList, genosDT, runsPerSet, seqType,
    qtlNumber, h2Levels, effectSizes, numCores=1){
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   LIBRARIES AND ASSERTIONS   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  require(genomalicious); require(doSNOW)

  # Check columns in genoDT
  check.cols <- sum(c('LINK.BLOCK','LOCUS','SAMPLE','GT') %in% colnames(genosDT))
  if(check.cols<4){
    stop(
      'Argument `genoDT` must have columns "LINK.BLOCK", "LOCUS", "SAMPLE", and "GT". See ?gam_resamp_sim_phenos.'
    )
  }

  # Number of SNP sets
  num.sets <- length(snpSetList)

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   SIMULATE PHENOTYPES   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # Data table of parameters for phenotype simulation
  paramDT <- CJ(SET=1:num.sets, HERIT=h2Levels, RUN=1:runsPerSet) %>%
    .[, PARAM:=1:.N]

  paramSplitList <- split(paramDT, by=c('SET','HERIT'))

  my.cluster <- makeCluster(numCores, type = "SOCK")
  registerDoSNOW(my.cluster)

  pb <- txtProgressBar(max=length(paramSplitList), style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)

  phenoSimsList <- foreach(Dparam=paramSplitList, .options.snow=opts) %dopar% {
    require(genomalicious)

    # Get set ID and heritability
    set.id <- Dparam$SET[1]
    herit <- Dparam$HERIT[1]

    # Focal loci, focal linkage blocks, reserved linkage blocks (only RRS)
    loci.focal <- snpSetList[[set.id]]$focal$LOCUS
    link.focal <- snpSetList[[set.id]]$focal$LINK.BLOCK
    link.reserved <- snpSetList[[set.id]]$reserved$LINK.BLOCK

    # Focal genotypes for the SNP set
    genoSubFocal <- genosDT[LOCUS %in% loci.focal]

    # RRS only: reserved genotypes for the SNP set
    if(seqType=='rrs'){
      genoSubReserved <- genosDT[LINK.BLOCK %in% link.reserved]
    }

    # WGS only: unique reserved loci
    if(seqType=='wgs'){
      uniqLociReserved <- unique(genosDT[!LOCUS %in% loci.focal]$LOCUS)
    }

    # Simulate phenotypes where loci are in the focal set
    phenoS1 <- lapply(Dparam$RUN, function(run.id){
      family_sim_qtl(
        famGenos = genoSubFocal, numLoci = qtlNumber,
        additiveVar = herit, environVar = 1-herit, effectSizes = effectSizes
      ) %>%
        .[['trait']] %>%
        data.table(AMAT='Simulated (with QTLs)', SET=set.id, HERIT=herit, RUN=run.id, .)
    })

    # Simulate phenotypes where loci are not in the focal set
    phenoS0 <- lapply(Dparam$RUN, function(run.id){
      if(seqType=='wgs'){
        # Draw loci
        loci.reserved <- sample(uniqLociReserved, qtlNumber, FALSE)
        # Simulate trait and phenotypes
        Dtrait <- family_sim_qtl(
          famGenos = genosDT[LOCUS %in% loci.reserved], qtlNames = loci.reserved,
          additiveVar = herit, environVar = 1-herit, effectSizes = effectSizes
        ) %>%
          .[['trait']] %>%
          data.table(AMAT='Simulated (without QTLs)', SET=set.id, HERIT=herit, RUN=run.id, .)
      } else if(seqType=='rrs'){
        # Draw loci
        loci.reserved <- genoSubReserved[, sample(LOCUS,1,FALSE), by=LINK.BLOCK] %>%
          .[['V1']] %>%
          sample(., qtlNumber, replace=FALSE)
        # Simulate trait and phenotypes
        Dtrait <- family_sim_qtl(
          famGenos = genoSubReserved[LOCUS %in% loci.reserved], qtlNames = loci.reserved,
          additiveVar = herit, environVar = 1-herit, effectSizes = effectSizes
        ) %>%
          .[['trait']] %>%
          data.table(AMAT='Simulated (without QTLs)', SET=set.id, HERIT=herit, RUN=run.id, .)
      }
    })

    # Return the phenotyoes
    rbind(rbindlist(phenoS1), rbindlist(phenoS0))
  }; stopCluster(my.cluster); gc()

  # Outout
  return(phenoSimsList)
}
