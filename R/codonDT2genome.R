#' Align a codon data.table to genomic positions
#'
#' Takes a codon table produced from \code{dna2codonDT} and aligns genomic
#' positions against exons and codons, using a set of annotations. Note that 
#' this function is designed for a single gene. To apply this function to
#' multiple genes, you will need to use a loop.
#'
#' @param codonDT Data.table: The output from \code{dna2codonDT}. Note, this
#' must have been generated using the exon functionality. It could be the
#' compressed or uncompressed version (see \code{?dna2codonDT}).
#' The required columns are:
#' \enumerate{
#'    \item \code{$EXON}: Integer, the exon ID.
#'    \item \code{$CODON}: Integer, the codon ID.
#'    \item \code{$NUC.GENE}: Integer/Character, the nucleotide ID within the
#'    exon coding sequence. An integer class if an uncompressed codon table, but
#'    a character class if a compressed codon table (see \code{?dna2codonDT}).
#' }
#'
#' @param genomeAnnot Data.table: The exon coding sequence annotations, e.g., from
#' a GFF3 annotations file. Each row should be an exon. The required columns are:
#' \enumerate{
#'    \item \code{$CHROM}: Character, the chromosome ID.
#'    \item \code{$START}: Integer, the starting position (on forward strand).
#'    \item \code{$END}: Integer, the end position (on forward strand).
#'    \item \code{$STRAND}: Character, the strand position, one of \code{'+'} (forward),
#'    or \code{'-'} (reverse). Note, that if the gene is on the reverse strand,
#'    the "end" of the exon (relative to the forward strand) is actually the
#'    "start" of the exon (relative to the reverse strand).
#'    \item \code{$EXON}: Integer, the exon ID.
#' }
#'
#' @param isCompressed Logical: Is the codon data table provided in \code{codonDT}
#' compressed? See \code{dna2codonDT}.
#'
#' @returns Returns the codon table but with the additional column \code{$CHROM}
#' and \code{$POS}, the chromosome and the genomic positions, respectively.
#'
#' @export

codonDT2genome <- function(codonDT, genomeAnnot, isCompressed){
  require(data.table)
  require(tidyverse)

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   ASSERTIONS   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  check_columns_codonDT <- sum(c('CODON','NUC.GENE','EXON') %in% colnames(codonDT))
  if(check_columns_codonDT!=3){
    stop('Argument `codonDT` requires columns $CODON, $NUC.GENE, and $EXON. See ?codonDT2genome.')
  }

  check_columns_genomeAnnot <- sum(c('CHROM','START','END','STRAND') %in% colnames(genomeAnnot))
  if(check_columns_genomeAnnot!=4){
    stop('Argument `genomeAnnot` requires columns $CHROM, $START, $END, and $STRAND See ?codonDT2genome.')
  }

  if(class(isCompressed)!='logical'){
    stop('Argument `isCompressed` should be a logical class. See ?codonDT2genome.')
  }

  check_compress_nucs <- sum(grepl('\\|', codonDT$NUC.GENE)) > 0
  if(isCompressed==TRUE & check_compress_nucs==FALSE){
    stop('The column $NUC.GENE in argument `codonDT` does not appear to be compressed as expected, given argument `isCompressed==TRUE`. Please check that specifications are correct. See ?codonDT2genome.')
  }
  if(isCompressed==FALSE & check_compress_nucs==TRUE){
    stop('The column $NUC.GENE in argument `codonDT` appears to be compressed, which is unexpected, given argument `isCompressed==FALSE`. Please check that specifications are correct. See ?codonDT2genome.')
  }

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   EXECUTE   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # Number of exons and strand position of the gene
  uniq.exons <- sort(unique(genomeAnnot$EXON))
  strand <- genomeAnnot$STRAND[1]

  # Create a table of nucleotide positions in the genome, relative to exons.
  nucGenome <- lapply(uniq.exons, function(i){
    st <- genomeAnnot[EXON==i]$START
    en <- genomeAnnot[EXON==i]$END
    chrom <- genomeAnnot[EXON==i]$CHROM
    if(strand=='+'){seq_i <- seq(st,en,1)
    } else if(strand=='-'){seq_i <- rev(seq(st,en,1))}
    data.table(CHROM=chrom,POS=seq_i,EXON=i)
  }) %>%
    do.call('rbind', .) %>%
    as.data.table %>%
    .[, NUC.GENE:=1:.N]

  if(isCompressed==FALSE){
    result <- left_join(nucGenome,codonDT)
  }

  if(isCompressed==TRUE){
    # For this, rename the $NUC.GENE column
    setnames(nucGenome, 'NUC.GENE', 'NUC.GENE.LONG')

    # Make a table of nucleotides within the gene with respect to codons
    nucGene <- codonDT %>%
      .[, .(NUC.GENE.LONG=str_split(NUC.GENE,'\\|')[[1]]),by=c('EXON','CODON','NUC.GENE')]
    nucGene[, NUC.GENE.LONG:=as.integer(NUC.GENE.LONG)]

    # Remake the table of nucleotides in the genome with respect to exons + codons
    nucGenome <- left_join(
      nucGenome,
      nucGene[, c('NUC.GENE','NUC.GENE.LONG','CODON')]
    )

    # Collapse the chromosomes
    chromCollapse <- nucGenome %>%
      .[, .(CHROM=paste(sort(unique(CHROM)),collapse='|')), by=c('CODON','NUC.GENE')]

    # Collapse the genomic positions
    if(strand=='+'){
      posCollapse <- nucGenome %>%
        .[, .(POS=paste(sort(POS),collapse='|')), by=c('CODON','NUC.GENE')]
    } else if(strand=='-'){
      posCollapse <- nucGenome %>%
        .[, .(POS=paste(rev(sort(POS)),collapse='|')), by=c('CODON','NUC.GENE')]
    }

    # The final result
    result <- left_join(codonDT, chromCollapse) %>%
      left_join(., posCollapse)
  }

  # Output
  return(result)
}
