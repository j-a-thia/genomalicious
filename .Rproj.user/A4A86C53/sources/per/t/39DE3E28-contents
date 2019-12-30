# CODE FOR BUILDING PACKAGE #

# Good website:
#   http://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html

# Developer libraries
libs <- c('devtools', 'roxygen2', 'testthat', 'knitr', 'data.table', 'tidyr')
for(L in libs){require(L, character.only=TRUE)}

# Make documents
roxygenise(clean=TRUE) # Sometimes this throws an error?
roxygenise()

# Load currently installed genomalicious
library(genomalicious)

# Make pool info data
genomalicious_PoolInfo <- data.table(POOL=paste0('Pop', 1:4), INDS=30)
save(genomalicious_PoolInfo, file='data/genomalicious_PoolInfo.RData')

# Make pool reads data
genomalicious_PoolReads <- fread('inst/extdata/genomalicious_PoolReads.csv')
save(genomalicious_PoolReads, file='data/genomalicious_PoolReads.RData')
refalt <- unique(genomalicious_PoolReads[, c('LOCUS', 'REF', 'ALT')])

# Make pool ref allele (pi) frequency estimate data
genomalicious_PoolPi <- fread('inst/extdata/genomalicious_PoolPi.csv')
genomalicious_PoolPi <- cbind(genomalicious_PoolPi, refalt[match(genomalicious_PoolPi$LOCUS, refalt$LOCUS), c('REF', 'ALT')]
                         , INDS=genomalicious_PoolInfo$INDS[match(genomalicious_PoolPi$POOL, genomalicious_PoolInfo$POOL)]
                         )
save(genomalicious_PoolPi, file='data/genomalicious_PoolPi.RData')

# Make a wide matrix of allele frequencies
genomalicious_Freqs <- as.data.frame(
                          spread(fread('inst/extdata/genomalicious_PoolPi.csv')[, c('POOL', 'LOCUS', 'PI')]
                            , key=LOCUS, value=PI))
rownames(genomalicious_Freqs) <- genomalicious_Freqs$POOL
genomalicious_Freqs <- as.matrix(genomalicious_Freqs[, -which(colnames(genomalicious_Freqs)=='POOL')])
save(genomalicious_Freqs, file='data/genomalicious_Freqs.RData')

# genomalicious_FreqsLong <- as.data.table(genomalicious_Freqs)
# genomalicious_FreqsLong$POP <- rownames(genomalicious_Freqs)
# genomalicious_FreqsLong <- melt(genomalicious_FreqsLong, id='POP', variable.name='LOCUS', value.name='FREQ')
# save(genomalicious_FreqsLong, file='data/genomalicious_FreqsLong.RData')

# Make the 4 pop genotype dataset
genomalicious_4pops <- readRDS('inst/extdata/genomalicious_FastSimCoal_30Diploids.RDS')
setnames(genomalicious_4pops, c('IND', 'GENO'), c('SAMPLE', 'GT'))
genomalicious_4pops$GT[genomalicious_4pops$GT=='1/0'] <- '0/1'
goodloci <- filter_maf(genomalicious_4pops, 0.01, 'genos')
genomalicious_4pops <- genomalicious_4pops[LOCUS %in% goodloci]
save(genomalicious_4pops, file='data/genomalicious_4pops.RData')

# Make a VCF from 4 pop genotypes
data("genomalicious_4pops")

sub4pops <- genomalicious_4pops[
  LOCUS %in% sample(x=unique(genomalicious_4pops$LOCUS), size=8, replace=FALSE)
  & SAMPLE %in% genomalicious_4pops[, unique(SAMPLE)[1:8], by=POP]$V1, ]

sub4pops[, DP:=rnbinom(nrow(sub4pops), mu=30, size=1)]

hist(sub4pops$DP)

sub4pops[, RO:=sum(rbinom(DP, size=1, prob=0.5)), by=c('SAMPLE', 'LOCUS')]
sub4pops[, AO:=DP-RO]

badDP4pops <- which(sub4pops$DP==0 | sub4pops$DP==1)
sub4pops$GT[badDP4pops] <- './.'

vcf4pops <- do.call('rbind', lapply(unique(sub4pops$LOCUS), function(locus){
  X <- sub4pops[LOCUS==locus]
  alleles <- sample(size=2, x=c('G', 'A', 'T', 'C'), replace=FALSE)

  col_names <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', X$SAMPLE)

  vcfRow <- matrix('.', nrow=1, ncol=9+nrow(X), dimnames=list(NULL, col_names))

  vcfRow[,'CHROM'] <- X$CHROM[1]
  vcfRow[1,'POS'] <- strsplit(locus, '_')[[1]][3]
  vcfRow[1,'REF'] <- alleles[1]
  vcfRow[1,'ALT'] <- alleles[2]
  vcfRow[1,'INFO'] <- paste0('DP=', sum(X$DP))
  vcfRow[1,'QUAL'] <- 30
  vcfRow[1,'FORMAT'] <- 'GT:DP:RO:AO'

  for(samp in unique(X$SAMPLE)){
    xsamp <- unlist(X[SAMPLE==samp, c('GT', 'DP', 'RO', 'AO')])
    vcfRow[, samp] <- paste(xsamp, collapse=':')
  }

  return(vcfRow)
}))

fwrite(vcf4pops, 'inst/extdata/genomalicious_indseq.vcf', sep='\t')

##This is a toy dataset for the R package genomalicious - it emulates a VCF file for individually genotyped data
##INFO=<ID=DP,Number=1,Type=Integer,Description='The total depth across samples'>
##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description='The total depth in a sample'>
##FORMAT=<ID=RO,Number=1,Type=Integer,Description='The reference allele counts in a sample'>
##FORMAT=<ID=AO,Number=1,Type=Integer,Description='The alternate allele counts in a sample'>


