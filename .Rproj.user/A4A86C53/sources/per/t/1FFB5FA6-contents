# Good website:
#   http://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html

libs <- c('devtools', 'roxygen2', 'testthat', 'knitr')
for(L in libs){require(L, character.only=T)}

roxygenise('./', clean=TRUE)

library(data.table)

library(genomalicious)

genomalicious_PoolInfo <- data.table(POOL=paste0('Pop', 1:4), INDS=30)
save(genomalicious_PoolInfo, file='data/genomalicious_PoolInfo.RData')

genomalicious_PoolReads <- fread('inst/extdata/genomalicious_PoolReads.csv')
save(genomalicious_PoolReads, file='data/genomalicious_PoolReads.RData')
refalt <- unique(genomalicious_PoolReads[, c('LOCUS', 'REF', 'ALT')])

genomalicious_PoolPi <- fread('inst/extdata/genomalicious_PoolPi.csv')
genomalicious_PoolPi <- cbind(genomalicious_PoolPi, refalt[match(genomalicious_PoolPi$LOCUS, refalt$LOCUS), c('REF', 'ALT')]
                         , INDS=genomalicious_PoolInfo$INDS[match(genomalicious_PoolPi$POOL, genomalicious_PoolInfo$POOL)]
                         )
save(genomalicious_PoolPi, file='data/genomalicious_PoolPi.RData')


# genomaliciousGenos <- fread('inst/extdata/genomaliciousGenos.csv')
# save(genomaliciousGenos, file='data/genomaliciousGenos.RData')

genomalicious_Freqs <- readRDS('inst/extdata/genomalicious_Freqs.RDS')
save(genomalicious_Freqs, file='data/genomalicious_Freqs.RData')

genomalicious_FreqsLong <- as.data.table(genomalicious_Freqs)
genomalicious_FreqsLong$POP <- rownames(genomalicious_Freqs)
genomalicious_FreqsLong <- melt(genomalicious_FreqsLong, id='POP', variable.name='LOCUS', value.name='FREQ')
save(genomalicious_FreqsLong, file='data/genomalicious_FreqsLong.RData')

genomalicious_4pops <- readRDS('inst/extdata/genomalicious_FastSimCoal_30Diploids.RDS')
setnames(genomalicious_4pops, c('IND', 'GENO'), c('SAMPLE', 'GT'))
genomalicious_4pops$GT[genomalicious_4pops$GT=='1/0'] <- '0/1'
goodloci <- filter_maf(genomalicious_4pops, 0.01, 'genos')
genomalicious_4pops <- genomalicious_4pops[LOCUS %in% goodloci]
save(genomalicious_4pops, file='data/genomalicious_4pops.RData')
