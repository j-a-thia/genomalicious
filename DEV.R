# Good website:
#   http://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html

libs <- c('devtools', 'roxygen2', 'testthat', 'knitr')
for(L in libs){require(L, character.only=T)}

library(data.table)
