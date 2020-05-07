#' Convert \code{fastsimcoal2} SNPs into diploid genotypes
#'
#' Takes in a \code{.arp} file produced from a \code{fastsimcoal2} (Excoffier et al. 2013)
#' simulation of SNP data and produces diploid genotypes from the sampled haploid genomes.
#' NOTE: This function was developed for even sample sizes (i.e. total number
#' of haploid genomes per population is divisible by 2) and for simulations where
#' all loci are polymorphic. Users should proceed with caution if their simulated data
#' does not fit these specifications.
#'
#' @param arpFile Character: The name of the \code{.arp} file produced from
#' the \code{fastsimecoal2} simulations.
#'
#' @param padIDs Logical: Should ID numbers be padded with '0'? Useful for
#' sorting character values. Defalt is TRUE.
#'
#' @details Diploid genotypes are created by combining pairs of haploid genomes
#' from each population (i.e. "Sample" 1 -> n in \code{fastsimcoal2}).
#'
#' @return Returns a data table with the following columns:
#' \enumerate{
#'    \item \code{$POP} Population ID (simualted samples, "Sample" 1 -> n;
#'    takes form \code{pop_[ID]}).
#'    \item \code{$SAMPLE} Individual ID (the diploid individual from each
#'    population; takes form \code{ind_[pop]_[ID]}).
#'    \item \code{$LOCUS} The locus ID (takes form \code{locus_[chrom]_[ID]}).
#'    \item \code{$CHROM} The chromosome ID (takes form \code{chrom_[ID]})
#'    \item \code{$GT} The genotype, scored as '/' separated alleles, where 0
#'    is ancestral, and 1 is derived.
#' }
#'
#' @references
#' Excoffier et al. (2013) Robust demographic inference from genomic and SNP data. PLoS Genetics.
#' @export
#'
fastsimcoal2_snps2dips <- function(arpFile, padIDs=TRUE){
  library(data.table); library(stringr)

  # Read data in line by line.
  dataIn <- readLines(arpFile)

  # Number of chromosomes
  numChrom <- dataIn[grep(pattern='#Number of independent chromosomes', x=dataIn)]
  numChrom <- as.integer(trimws(strsplit(numChrom, ':')[[1]][2]))

  # Get the line number for chromosome polymorphisms
  polyLines <- grep(pattern='polymorphic positions on chromosome', x=dataIn)

  # Make locus names using the chromosome and SNP number
  chrom_loci <- unlist(
    lapply(polyLines, function(n.line){
      chrom <- trimws(strsplit(dataIn[n.line], 'chromosome')[[1]][2])

      loci <- trimws(strsplit(dataIn[n.line+1], ',')[[1]])
      loci <- gsub(pattern='#', replacement='', x=loci)

      if(padIDs==TRUE){
        chrom <- str_pad(chrom, nchar(numChrom), 'left', '0')
        loci <- str_pad(loci, nchar(max(as.integer(loci))), 'left', '0')
      }

      return(paste0('locus_', chrom, '_', loci))
    }))

  # Get the number of samples (populations)
  numSamps <- dataIn[grep(pattern='NbSamples=', x=dataIn)]
  numSamps <- strsplit(numSamps, '=')[[1]][2]

  # Get the lines that sample data begins
  sampLines <- grep(pattern='SampleName=', x=dataIn)

  # For each nth sample (population), create diploid genoypes
  # by combining pairs of haploid genomes
  diploids <- lapply(1:numSamps, function(n.samp){
    # For the nth sample, what line does it start on?
    st.samp <- sampLines[n.samp]

    # How big is the (haploid) sample size?
    samp.size <- strsplit(dataIn[st.samp+1], '=')[[1]][2]

    # Extract the haploid genomes for the nth sample
    hap.genomes <- dataIn[st.samp + 2 + 1:samp.size]
    hap.genomes <- strsplit(hap.genomes, '\t')

    # How many items per sampled genome (SNPs are the last item)
    n.items <- length(hap.genomes[[1]])

    # Get the SNP data for each genome
    hap.genomes <- unlist(
      lapply(hap.genomes, function(genome){ trimws(genome[n.items]) })
    )

    # Convert string into individual alleles for haploid genomes
    hap.genomes <- do.call('rbind'
                           , lapply(hap.genomes, function(genome){
                             alleles <- strsplit(genome, '')[[1]]
                             return(matrix(alleles, nrow=1, ncol=length(alleles)))
                           }))

    # Get the genotypes for all diploid individuals in this sample.
    # Simply combine two haploid genomes to make a genotype
    samp.gts <- do.call('rbind'
                        , lapply(seq(1, nrow(hap.genomes), by=2), function(gt.first){
                          gts <- apply(hap.genomes[c(gt.first, gt.first+1),], 2, function(x){
                            paste(sort(x), collapse='/') })
                          gts <- matrix(gts, nrow=1, ncol=length(gts))
                          return(gts)
                        }))

    # Rename the columns to be loci
    colnames(samp.gts) <- chrom_loci

    # Make a data table in long-format
    if(padIDs==TRUE){
      pop.name <- paste0('pop_', str_pad(n.samp, nchar(numSamps), 'left', '0'))
      samp.name <- paste0('ind_', n.samp, '_', str_pad(1:nrow(samp.gts), nchar(nrow(samp.gts)), 'left', '0'))
    } else{
      pop.name <- paste0('pop_', n.samp)
      samp.name <- paste0('ind_', n.samp, '_', 1:nrow(samp.gts))
    }

    samp.gts <- data.table(POP=pop.name, SAMPLE=samp.name, as.data.table(samp.gts))
    samp.gts <- melt(samp.gts, id.vars=c('POP', 'SAMPLE'), variable.name='LOCUS', value.name='GT')

    return(samp.gts)
  })

  # Combine population list items into single data.table
  diploids <- do.call('rbind', diploids)

  # Make sure loci are coded as characters.
  diploids[, LOCUS:=as.character(LOCUS)]

  # Create a chromosome column
  diploids[, CHROM:=paste0('chrom_', strsplit(LOCUS, split='_')[[1]][2])]

  return(diploids)
}
