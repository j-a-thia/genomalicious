% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genoscore_converter.R
\name{genoscore_converter}
\alias{genoscore_converter}
\title{Convert between genotype scores (separated alleles vs counts)}
\usage{
genoscore_converter(dat)
}
\arguments{
\item{dat}{Character/Integer: A vector of genotypes. If the class
is \code{'character'}, then assumes separated alleles and converts into
allele counts. The opposite is true, if class is \code{'integer'}, will
convert into separated alleles.}
}
\description{
Assumes biallelic genotypes and can interchange between separated
alleles ('0/0', '0/1', '1/1') or Alt allele counts (0, 1, 2).
}
\details{
If \code{dat} is a vector of characters, then missing values should
take the form of './.'. Otherwise, if \code{dat} is a vector of integers,
missing values should take the form of \code{NA}.
}
\examples{
library(genomalicious)
genoscore_converter(c('0/0', '0/1', '1/1'))
genoscore_converter(c(0, 1, 2))

}
