# CODE FOR BUILDING PACKAGE #

# Good website:
#   http://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html

# Developer libraries
libs <- c('devtools', 'roxygen2', 'testthat', 'knitr', 'data.table', 'tidyverse')
for(L in libs){require(L, character.only=TRUE)}

# Make all documents
roxygenise(clean=TRUE)

# Make just those documents that have changed
roxygenise()

# Load currently installed genomalicious
library(genomalicious)

# # Make the 4 pop genotype dataset
# fsc_genos <- fread('inst/extdata/fsc2_sim_radseq_1_1.gen', skip=1)
# fsc_head <- colnames(fread('inst/extdata/fsc2_sim_radseq_1_1.gen', nrow=0))
#
# fsc_tab <- fsc_genos[, 1:(ncol(fsc_genos)-1)] %>%
#   setnames(., new=fsc_head) %>%
#   setnames(
#     .,
#     old=c('Chrom','Pos','Anc_all','Der_all'),
#     new=c('CHROM','POS','REF','ALT')
#   ) %>%
#   melt(
#     .,
#     id.vars=c('CHROM','POS','REF','ALT'),
#     variable.name='SAMPLE',
#     value.name='GT'
#   ) %>%
#   .[, SAMPLE:=sub('G_', 'Ind', SAMPLE)] %>%
#   .[, POP:=sub('Ind', 'Pop', sub('\\_.*', '', SAMPLE))] %>%
#   .[, CHROM:=paste0('Contig', CHROM)] %>%
#   .[, LOCUS:=paste(CHROM, POS, sep='_')] %>%
#   .[, c('CHROM','POS','LOCUS','POP','SAMPLE','GT')] %>%
#   print
#
# fsc_tab[, length(unique(LOCUS)), by=CHROM]$V1 %>%  table
#
# keep.loci <- fsc_tab[LOCUS %in% filter_maf(fsc_tab, type='genos', maf=0.05)]$LOCUS %>%
#   unique() %>% .[1:200]
#
# fsc_tab[LOCUS %in% keep.loci, length(unique(LOCUS)), by=CHROM]$V1 %>%  table
#
# data_Genos <- fsc_tab[LOCUS %in% keep.loci] %>%
#   .[, DP:=rnbinom(n=1, size=15, prob=0.3), by=c('CHROM','SAMPLE')] %>%
#   .[, AO:=rbinom(n=1, size=DP, prob=GT/2), by=c('CHROM','SAMPLE')] %>%
#   .[, RO:=DP-AO]
#
# locs_Genos <- data_Genos[, c('LOCUS')] %>% unique()
# locs_Genos <- cbind(
#   locs_Genos,
#   sapply(1:nrow(locs_Genos), function(i){
#     sample(c('T','A','C','G'), 2, replace=FALSE)
#   }) %>%
#     t() %>%
#     as.data.table() %>%
#     setnames(., new=c('ALT','REF'))
# )
#
# data_Genos <- left_join(data_Genos, locs_Genos)
#
# data_Genos %>%  pca_genos(., popCol='POP') %>%  pca_plot
#
# save(data_Genos, file='data/data_Genos.RData')
#
# # Make pool data
# data_PoolFreqs <- data_Genos %>%
#   .[, .(FREQ=sum(GT)/(length(GT)*2)), by=c('CHROM','POS','LOCUS','ALT','REF','POP')] %>%
#   .[, DP:=rnbinom(1, 20, prob=0.2), by=c('POP','CHROM','LOCUS')] %>%
#   .[, AO:=rbinom(n=1, size=DP, prob=FREQ), by=c('POP','CHROM','LOCUS')] %>%
#   .[, RO:=DP-AO, by=c('POP','CHROM','LOCUS')] %>%
#   .[, FREQ:=AO/DP, by=c('POP','CHROM','LOCUS')] %>%
#   .[, POOL:=as.integer(sub('Pop', '', POP))]
# save(data_PoolFreqs, file='data/data_PoolFreqs.RData')
#
# # Make pool info data
# data_PoolInfo <- data.table(POOL=1:4, INDS=30)
# save(data_PoolInfo, file='data/data_PoolInfo.RData')
#
# # Make a frequency matrix
# data_FreqsMat <- data_PoolFreqs %>%
#   dcast(., POOL~LOCUS, value.var='FREQ') %>%
#   as.data.frame() %>%
#   column_to_rownames(., 'POOL') %>%
#   as.matrix()
# save(data_FreqsMat,file='data/data_FreqsMat.RData')
#
# # Make VCFs
# vcf_head <- "##This is a toy dataset for the R package genomalicious - it emulates a VCF file
# ##INFO=<ID=DP,Number=1,Type=Integer,Description='The total depth across samples'>
# ##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>
# ##FORMAT=<ID=DP,Number=1,Type=Integer,Description='The total depth in a sample'>
# ##FORMAT=<ID=RO,Number=1,Type=Integer,Description='The reference allele counts in a sample'>
# ##FORMAT=<ID=AO,Number=1,Type=Integer,Description='The alternate allele counts in a sample'>"
#
# # Make a VCF from 4 pop genotypes
# vcf4pops <- do.call('rbind', lapply(unique(data_Genos$LOCUS), function(locus){
#   X <- data_Genos[LOCUS==locus]
#   X[, GT:=genoscore_converter(GT)]
#
#   vcfRow <- data.table(
#     CHROM=X$CHROM[1],
#     POS=X$POS[1],
#     REF=X$REF[1],
#     ALT=X$ALT[1],
#     INFO=paste0('DP=', sum(X$DP)),
#     QUAL=30,
#     FORMAT='GT:DP:RO:AO',
#     X[, paste(GT,DP,RO,AO,sep=":"), by=c('LOCUS','SAMPLE')] %>%
#       dcast(., LOCUS ~ SAMPLE, value.var='V1') %>%
#       .[, !'LOCUS']
#   )
#
#   setnames(vcfRow, 'CHROM', '#CHROM')
#
#   return(vcfRow)
# }))
#
# fwrite(vcf4pops, 'inst/extdata/data_indseq.vcf', sep='\t', quote=FALSE)
# readLines('inst/extdata/data_indseq.vcf') %>%
#   c(vcf_head, .) %>%
#   writeLines(., 'inst/extdata/data_indseq.vcf')
#
# vcfPools <- vcf4pops <- do.call('rbind', lapply(unique(data_PoolFreqs$LOCUS), function(locus){
#   X <- data_PoolFreqs[LOCUS==locus] %>%
#     copy %>%
#     .[, GT:='0/1']
#
#   vcfRow <- data.table(
#     CHROM=X$CHROM[1],
#     POS=X$POS[1],
#     REF=X$REF[1],
#     ALT=X$ALT[1],
#     INFO=paste0('DP=', sum(X$DP)),
#     QUAL=30,
#     FORMAT='GT:DP:RO:AO',
#     X[, paste(GT,DP,RO,AO,sep=":"), by=c('LOCUS','POOL')] %>%
#       dcast(., LOCUS ~ POOL, value.var='V1') %>%
#       .[, !'LOCUS']
#   )
#
#   setnames(vcfRow, 'CHROM', '#CHROM')
#
#   return(vcfRow)
# }))
#
# fwrite(vcfPools, 'inst/extdata/data_poolseq.vcf', sep='\t', quote=FALSE)
# readLines('inst/extdata/data_poolseq.vcf') %>%
#   c(vcf_head, .) %>%
#   writeLines(., 'inst/extdata/data_poolseq.vcf')

