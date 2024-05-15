#' Calculate heterozygosity from genotypes or allele frequencies
#'
#' This function takes in a long format data table of genotypes or allele
#' frequencies (allelic counts) and calculates heterozygosity. Heterozygosity
#' is calculated as either the SNP-wise heterozygosity (at each SNP, excluding
#' monomorphic sites), or as the genomic heterozygosity (the SNP-wise
#' heterozygosity standardised for total sites assayed, including monomorphic
#' sites). For indivdiual genotypes, you can use multiallelic data.
#' For population allele frequencies, you must code as biallelic reference and
#' alternate alleles.
#'
#' @param snpData Data table: Genotypes for individuals or frquencies (read coiunts)
#' for populations in long-format. See Details for parameterisation.
#'
#' @param extraData Data table: Extra info for individual or populations.
#' See Details for parameterisation.
#'
#' @param type Character: One of \code{'genos'} or \code{'freqs'}, to calculate
#' heterozygosity on genotype or frequencies, respectively.
#'
#' @param method Charachter: One of \code{'genomic'} or \code{'snpwise'} to
#' calculate genomic or SNP-wise heterozygosity, respectively.
#'
#' @param chromCol Character: The column with chromosome ID.
#' Default is \code{'CHROM'}.
#'
#' @param posCol Character: The column with the positional info.
#' Default is \code{'POS'}.
#'
#' @param locusCol Character: The column with locus ID.
#' Default is \code{'LOCUS'}.
#'
#' @param genoCol Character: The column with the genotype info.
#' Default is \code{'GT'}. Genotypes should be scored as alleles separated by
#' '/', e.g., '0/0', '0/1', '1/1', etc.
#'
#' @param roCol Character: The column with reference allele counts.
#' Default is \code{'RO'}.
#'
#' @param aoCol Character: The column with alternate allele counts.
#' Default is \code{'AO'}.
#'
#' @param sampCol Character: The column name with the sampled individual ID.
#' Default is \code{'SAMPLE'}.
#'
#' @param popCol Character: The column name with population ID.
#' Default is \code{'POP'}.
#'
#' @param covCol Character: The column with the number of genomic sites covered
#' per chromosome. Default is \code{'COV.SITES'}.
#'
#' @param indsCol Character: The column name with the number of pooled individuals
#' per population per chromosome. Default is \code{'INDS'}.
#'
#' @details The genomic heterozygosity (also known as the autosomal
#' heterozygosity) has been demonstrated to be the more accurate and robust
#' measure of heterozygosity (Schmidt et al. 2021). SNP-wise heterozygosity can
#' suffer from sampling biases (filtering, missing data, sample size, etc).
#' Note, estimates of genomic heterozygosity are orders of magnitude less than
#' SNP-wise heterozygosity (so do not be alarmed if you see very small values!).
#'
#' Heterozygsity for population pools is calculated using the method from
#' Ferretti et al. (2013). You can calculate genomic heterozygosity for population
#' pools (standardising by total covered sites), or SNP-wise heterozygosity too
#' (standardising by the number of polymorphic sites observed).
#'
#' You must specify both \code{type} and \code{method}. This will dicate the
#' required column needed for \code{snpData} and \code{extraData}.
#'
#' If \code{type=='genos'} and \code{method=='genomic'}:
#' \enumerate{
#'    \item \code{snpData} requires columns specified in:
#'    \code{chromCol}, \code{posCol}, \code{locusCol} and \code{sampCol}
#'
#'    \item \code{extraData} requires columns specified in:
#'    \code{chromCol}, \code{sampCol}, and \code{covCol}.
#' }
#'
#' If \code{type=='genos'} and \code{method=='snpwise'}
#' \enumerate{
#'    \item \code{snpData} requires columns specified in:
#'    \code{chromCol}, \code{posCol}, \code{locusCol} and \code{sampCol}
#'
#'    \item \code{extraData} will NOT be used.
#' }
#'
#' If \code{type=='freqs'} and \code{method=='genomic'}:
#' \enumerate{
#'    \item \code{snpData} requires columns specified in:
#'    \code{chromCol}, \code{posCol}, \code{locusCol}, \code{popCol},
#'    \code{roCol}, and \code{aoCol}.
#'
#'    \item \code{extraData} requires columns specified in:
#'    \code{chromCol}, \code{sampCol}, \code{covCol}, and \code{indsCol}.
#' }
#'
#' If \code{type=='freqs'} and \code{method=='snpwise'}:
#' \enumerate{
#'    \item \code{snpData} requires columns specified in:
#'    \code{chromCol}, \code{posCol}, \code{locusCol}, \code{popCol},
#'    \code{roCol}, and \code{aoCol}.
#'
#'    \item \code{extraData} requires columns specified in:
#'    \code{chromCol}, \code{sampCol}, and \code{indsCol}.
#' }
#'
#' @return Returns a data table of heterozygosity estimates per sample or population.
#' Genomic heterozygosity is reported per chromosome, SNP-wise heterozygosity is
#' reported across all SNPs.
#'
#' @references
#' Ferretti et al. (2013) Molecular Ecology. DOI: 10.1111/mec.12522 \cr
#' Schmidt et al. (2021) Methods in Ecology and Evolution. DOI: 10.1111/2041-210X.13659
#'
#' @examples
#' library(genomalicious)
#'
#' data(data_Genos)
#' data(data_PoolFreqs)
#'
#' # Convert genos to characters
#' data_Genos[, GT:=genoscore_converter(GT)]
#'
#' # Make extra data for the samples and populatin pools
#' extraSampInfo <- CJ(
#'   SAMPLE=unique(data_Genos$SAMPLE),
#'   CHROM=unique(data_Genos$CHROM),
#'   COV.SITES=150
#'   )
#'
#' extraPoolInfo <- CJ(
#'   POP=unique(data_PoolFreqs$POP),
#'   CHROM=unique(data_PoolFreqs$CHROM),
#'   COV.SITES=150,
#'   INDS=30
#'   )
#'
#' # Genomic heterozygosity of individuals, per chromosome/contig
#' het_calc(data_Genos, extraSampInfo, type='genos', method='genomic')
#'
#' # SNP-wise heterozygosity of indivdiuals, SNP-wise
#' het_calc(data_Genos, extraSampInfo, type='genos', method='snpwise')
#'
#' # Genomic heterozygosity of population pools, per chromosome/contig
#' het_calc(data_PoolFreqs, extraPoolInfo, type='freqs', method='genomic')
#'
#' # Genomic heterozygosity of population pools, SNP-wise
#' het_calc(data_PoolFreqs, extraPoolInfo, type='freqs', method='snpwise')
#' @export

het_calc <- function(
    snpData, extraData, method, type, chromCol='CHROM', posCol='POS', locusCol='LOCUS',
    sampCol='SAMPLE', genoCol='GT', popCol='POP', roCol='RO', aoCol='AO',
    covCol='COV.SITES', indsCol='INDS'
){
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   CHECK AND ENVIRONMENT   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # Check type
  if(!type %in% c('genos','freqs')){
    stop('Argument `type` must be one of "genos" or "freqs". See ?het_calc.')
  }

  # Check method
  if(!method %in% c('genomic','snpwise')){
    stop('Argument `method` must be one of "genomic" or "snpwise". See ?het_calc.')
  }

  # Check the columns in snpData and extraData
  if(type=='genos'){
    # Genomic heterozygosity
    if(method=='genomic'){
      # SNP data
      check.col.snp <- c(chromCol, posCol, locusCol, sampCol, genoCol)
      if(length(check.col.snp)!=5){
        stop('Argument `type`=="genos" and `method`=="genomic". All columns `chromCol`, `posCol`, `locusCol` and `sampCol` must be in `snpData`. See ?het_calc.')
      }
      snpData <- copy(snpData) %>%
        setnames(., check.col.snp, c('CHROM','POS','LOCUS','SAMPLE','GT'))
      # Extra data
      check.col.extra <- c(chromCol, sampCol, covCol)
      if(length(check.col.extra)!=3){
        stop('Argument `type`=="genos" and `method`=="genomic". All columns `chromCol`, `sampCol`, and `covCol` must be in `extraData`. See ?het_calc.')
      }
      extraData <- copy(extraData) %>%
        setnames(., check.col.extra, c('CHROM','SAMPLE','COV.SITES'))
    }
    # SNP-wise heterozygosity
    if(method=='snpwise'){
      # SNP data
      check.col.snp <- c(chromCol, posCol, locusCol, sampCol, genoCol)
      if(length(check.col.snp)!=5){
        stop('Argument `type`=="genos" and `method`=="snpwise". All columns `chromCol`, `posCol`, `locusCol` and `sampCol` must be in `snpData`. See ?het_calc.')
      }
      snpData <- copy(snpData) %>%
        setnames(., check.col.snp, c('CHROM','POS','LOCUS','SAMPLE','GT'))
    }
  }

  if (type=='freqs'){
    # Genomic heterozygosity
    if(method=='genomic'){
      # SNP data
      check.col.snp <- c(chromCol, posCol, locusCol, popCol, roCol, aoCol)
      if(length(check.col.snp)!=6){
        stop('Argument `type`=="freqs" and `method`=="genomic". All columns `chromCol`, `posCol`, `locusCol`, `popCol`, `roCol` and `aoCol` must be in `snpData`. See ?het_calc.')
      }
      snpData <- copy(snpData) %>%
        setnames(., check.col.snp, c('CHROM','POS','LOCUS','POP','RO','AO'))
      # Extra data
      check.col.extra <- c(chromCol, popCol, covCol, indsCol)
      if(length(check.col.extra)!=4){
        stop('Argument `type`=="freqs" and `method`=="genomic". All columns `chromCol`, `sampCol`, `covCol`, and `indsCol` must be in `extraData`. See ?het_calc.')
      }
      extraData <- copy(extraData) %>%
        setnames(., check.col.extra, c('CHROM', 'POP', 'COV.SITES', 'INDS'))
    }
    # SNP-wise heterozygosity
    if(method=='snpwise'){
      # SNP data
      check.col.snp <- c(chromCol, posCol, locusCol, popCol, roCol, aoCol)
      if(length(check.col.snp)!=6){
        stop('Argument `type`=="freqs" and `method`=="snpwise". All columns `chromCol`, `posCol`, `locusCol`, `popCol`, `roCol` and `aoCol` must be in `snpData`. See ?het_calc.')
      }
      snpData <- copy(snpData) %>%
        setnames(., check.col.snp, c('CHROM', 'POS', 'LOCUS', 'POP', 'RO', 'AO'))
      # Extra data
      check.col.extra <- c(chromCol, popCol, indsCol)
      if(length(check.col.extra)!=3){
        stop('Argument `type`=="freqs" and `method`=="snpwise". All columns `chromCol`, `sampCol`, and `indsCol` must be in `extraData`. See ?het_calc.')
      }
      extraData <- copy(extraData) %>%
        setnames(., check.col.extra, c('CHROM', 'POP', 'INDS'))
    }
  }

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   GENOMIC HETEROZYGOSITY   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # For genomic heterozygosity, iterate over combinations of
  # chromosome and populations in the extra data table.
  if(method=='genomic'){
    # Using individual genotypes
    if(type=='genos'){
      snpData[, ALLELE.1:=sub('/.*','',GT), by=c('SAMPLE','LOCUS')]
      snpData[, ALLELE.2:=sub('.*/','',GT), by=c('SAMPLE','LOCUS')]
      snpData[, HET:=(ALLELE.1!=ALLELE.2), by=c('SAMPLE','LOCUS')]

      result <- snpData[, .(HET=sum(HET)), by=c('CHROM','SAMPLE')] %>%
        merge.data.table(., extraData) %>%
        .[, HET:=HET/COV.SITES] %>%
        .[, c('CHROM','SAMPLE','HET')]
    }
    # Using population allele frequencies
    if(type=='freqs'){
      snpData <- snpData %>%
        .[, RO:=as.integer(RO)] %>%
        .[, AO:=as.integer(AO)] %>%
        .[, DP:=AO+RO]

      result <- lapply(1:nrow(extraData), function(i){
        chrom <- extraData$CHROM[i]
        pop <- extraData$POP[i]
        cov.sites <- extraData$COV.SITES[i]
        inds <- extraData$INDS[i]

        allele.frac <- snpData[POP==pop] %>%
          .[, NUMER:=2*(AO*(DP-AO)), by=LOCUS] %>%
          .[, DENOM:=DP*(DP-1)] %>%
          .[, sum(NUMER/DENOM)]

        het <- (inds/(inds-1)) * (1/cov.sites) * allele.frac

        data.table(CHROM=chrom, POP=pop, HET=het)
      }) %>%
        do.call('rbind',.)
    }
  }

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ####   SNP-WISE HETEROZYGOSITY   ####
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  if(method=='snpwise'){
    # Using individual allele frequencies
    if(type=='genos'){
      snpData[, ALLELE.1:=sub('/.*','',GT), by=c('SAMPLE','LOCUS')]
      snpData[, ALLELE.2:=sub('.*/','',GT), by=c('SAMPLE','LOCUS')]
      snpData[, HET:=(ALLELE.1!=ALLELE.2), by=c('SAMPLE','LOCUS')]
      result <- snpData[, .(HET=sum(HET)/length(HET)), by=c('SAMPLE')]
    }
    # Using population allele frequencies
    if(type=='freqs'){
      snpData <- snpData %>%
        .[, RO:=as.integer(RO)] %>%
        .[, AO:=as.integer(AO)] %>%
        .[, DP:=AO+RO]

      result <- lapply(1:nrow(extraData), function(i){
        chrom <- extraData$CHROM[i]
        pop <- extraData$POP[i]
        inds <- extraData$INDS[i]
        num.sites <- snpData[POP==pop] %>% nrow

        allele.frac <- snpData[POP==pop] %>%
          .[, NUMER:=2*(AO*(DP-AO)), by=LOCUS] %>%
          .[, DENOM:=DP*(DP-1)] %>%
          .[, sum(NUMER/DENOM)]

        het <- (inds/(inds-1)) * (1/num.sites) * allele.frac

        data.table(CHROM=chrom, POP=pop, HET=het)
      }) %>%
        do.call('rbind',.)
    }
  }

  # Output
  return(result)
}
