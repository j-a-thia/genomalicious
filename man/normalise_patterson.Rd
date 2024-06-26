% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalise_patterson.R
\name{normalise_patterson}
\alias{normalise_patterson}
\title{Perform the Patterson et al. (2006) normalisation to a genotype matrix}
\usage{
normalise_patterson(dat)
}
\arguments{
\item{dat}{Matrix: Counts of Alt allele per individual (0, 1, or 2)
, e.g. the C(i,j) table described in Patterson et al. (2006).}
}
\value{
Returns a matrix with the same dimensions as \code{dat}, but with
genotypes normalised as per Patterson et al. (2006).
}
\description{
Perform the Patterson et al. (2006) normalisation to a genotype matrix
}
\references{
Patterson et al. (2006) Population structure and eigenanalysis. PLOS Genetics.
}
