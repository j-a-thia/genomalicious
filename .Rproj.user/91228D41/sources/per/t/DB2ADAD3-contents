# CODE FOR BUILDING PACKAGE #

# Good website:
#   http://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html

# Developer libraries
libs <- c('devtools', 'roxygen2', 'testthat', 'knitr', 'data.table', 'tidyverse')
for(L in libs){require(L, character.only=TRUE)}

# Make documents
roxygenise(clean=TRUE) # Sometimes this throws an error?
roxygenise()

# Load currently installed genomalicious
library(genomalicious)

# Make pool info data
data_PoolInfo <- data.table(POOL=paste0('Pop', 1:4), INDS=30)
save(data_PoolInfo, file='data/data_PoolInfo.RData')

# Make pool reads data
data_PoolReads <- fread('inst/extdata/data_PoolReads.csv')
save(data_PoolReads, file='data/data_PoolReads.RData')
refalt <- unique(data_PoolReads[, c('LOCUS', 'REF', 'ALT')])

# Make pool ref allele (pi) frequency estimate data
data_PoolFreqs <- fread('inst/extdata/data_PoolFreqs.csv')
data_PoolFreqs <- cbind(data_PoolFreqs, refalt[match(data_PoolFreqs$LOCUS, refalt$LOCUS), c('REF', 'ALT')]
                         , INDS=data_PoolInfo$INDS[match(data_PoolFreqs$POOL, data_PoolInfo$POOL)]
                         )
save(data_PoolFreqs, file='data/data_PoolFreqs.RData')

# Make a wide matrix of allele frequencies
data_FreqsMat <- as.data.frame(
                          spread(fread('inst/extdata/data_PoolFreqs.csv')[, c('POOL', 'LOCUS', 'FREQ')]
                            , key=LOCUS, value=FREQ))
rownames(data_FreqsMat) <- data_FreqsMat$POOL
save(data_FreqsMat, file='data/data_FreqsMat.RData')

# Make the 4 pop genotype dataset
fsc_genos <- fread('inst/extdata/fsc2_radseq_sim_1_1.gen', skip=1)
fsc_head <- colnames(fread('inst/extdata/fsc2_radseq_sim_1_1.gen', nrow=0))

fsc_tab <- fsc_genos[, 1:(ncol(fsc_genos)-1)] %>%
  setnames(., new=fsc_head) %>%
  setnames(
    .,
    old=c('Chrom','Pos','Anc_all','Der_all'),
    new=c('CHROM','POS','REF','ALT')
  ) %>%
  melt(
    .,
    id.vars=c('CHROM','POS','REF','ALT'),
    variable.name='SAMPLE',
    value.name='GT'
  ) %>%
  .[, SAMPLE:=as.character(SAMPLE)] %>%
  .[, SAMPLE:=sub('A_', 'Ind', SAMPLE)] %>%
  .[, POP:=sub('Ind', 'Pop', sub('\\_.*', '', SAMPLE))] %>%
  .[, CHROM:=paste0('Contig', CHROM)] %>%
  .[, LOCUS:=paste(CHROM, POS, sep='_')] %>%
  .[, c('CHROM','POS','LOCUS','POP','SAMPLE','GT')] %>%
  print

fsc_tab[, length(unique(LOCUS)), by=CHROM]$V1 %>%  table

keep.loci <- fsc_tab[LOCUS %in% filter_maf(fsc_tab, type='genos', maf=0.05)]$LOCUS %>%
  unique() %>% .[1:500]

fsc_tab[LOCUS %in% keep.loci, length(unique(LOCUS)), by=CHROM]$V1 %>%  table

fsc_tab[GT=='1/0', GT:='0/1']
data_4pops <- fsc_tab[LOCUS %in% keep.loci]

data_4pops %>%  pca_genos(., popCol='POP') %>%  pca_plot

save(data_4pops, file='data/data_4pops.RData')

# Make a VCF from 4 pop genotypes
sub4pops <- data_4pops[
  LOCUS %in% sample(x=unique(data_4pops$LOCUS), size=8, replace=FALSE)
  & SAMPLE %in% data_4pops[, unique(SAMPLE)[1:8], by=POP]$V1, ]

sub4pops[, DP:=rnbinom(nrow(sub4pops), mu=30, size=2)]

hist(sub4pops$DP)

sub4pops <-
  apply(sub4pops, 1, function(xx){
    dp <- as.integer(xx['DP'])
    if(xx['GT']=='0/0'){
      ro <- dp; ao <- 0
    } else if(xx['GT']=='0/1'){
      ro <- sum(rbinom(dp/2, size=1, prob=1))
      ao <- dp - ro
    } else{
      ro <- 0; ao <- dp
    }
    return(data.table(RO=ro, AO=ao))
  }) %>%
  do.call('rbind', .) %>%
  cbind(sub4pops, .) %>%
  as.data.table()

vcf4pops <- do.call('rbind', lapply(unique(sub4pops$LOCUS), function(locus){
  X <- sub4pops[LOCUS==locus]
  alleles <- sample(size=2, x=c('G', 'A', 'T', 'C'), replace=FALSE)

  col_names <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', X$SAMPLE)

  vcfRow <- matrix('.', nrow=1, ncol=9+nrow(X), dimnames=list(NULL, col_names))

  vcfRow[,'CHROM'] <- X$CHROM[1]
  vcfRow[1,'POS'] <- X$POS[1]
  vcfRow[1,'REF'] <- alleles[1]
  vcfRow[1,'ALT'] <- alleles[2]
  vcfRow[1,'INFO'] <- paste0('DP=', sum(X$DP))
  vcfRow[1,'QUAL'] <- 30
  vcfRow[1,'FORMAT'] <- 'GT:DP:RO:AO'

  for(samp in unique(X$SAMPLE)){
    xsamp <- unlist(X[SAMPLE==samp, c('GT', 'DP', 'RO', 'AO')])
    vcfRow[, samp] <- paste(xsamp, collapse=':')
  }

  vcfRow <- as.data.table(vcfRow)

  setnames(vcfRow, 'CHROM', '#CHROM')

  return(vcfRow)
}))

fwrite(vcf4pops, 'inst/extdata/data_indseq.vcf', sep='\t', quote=FALSE)


