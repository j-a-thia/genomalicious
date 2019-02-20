# Good website:
#   http://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html

libs <- c('devtools', 'roxygen2', 'testthat', 'knitr')
for(L in libs){require(L, character.only=T)}

library(data.table)

genomaliciousInfo <- data.table(POOL=paste0('Pop', 1:4), INDS=30)
save(genomaliciousInfo, file='data/genomaliciousInfo.RData')

genomaliciousReads <- fread('inst/extdata/genomaliciousReads.csv')
save(genomaliciousReads, file='data/genomaliciousReads.RData')

genomaliciousPi <- fread('inst/extdata/genomaliciousPi.csv')
save(genomaliciousPi, file='data/genomaliciousPi.RData')

genomaliciousGenos <- fread('inst/extdata/genomaliciousGenos.csv')
save(genomaliciousGenos, file='data/genomaliciousGenos.Rdata')

fwrite(pgposerPi, 'inst/extdata/genomaliciousPi.csv')
fwrite(pgposerReads, 'inst/extdata/genomaliciousReads.csv')
fwrite(pgposerGenos, 'inst/extdata/genomaliciousGenos.csv')

for(oldfile in list.files('./inst/extdata/')){
  newfile <- gsub('pgposerDemo', 'genomalicious', oldfile)
  file.rename(from=paste0('./inst/extdata/', oldfile), to=paste0('./inst/extdata/', newfile))
}
