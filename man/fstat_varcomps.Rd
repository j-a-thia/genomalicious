% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fstat_calc.R
\name{fstat_varcomps}
\alias{fstat_varcomps}
\title{Calculate variance components for F-statistics}
\usage{
fstat_varcomps(pi, ni, r, hi)
}
\arguments{
\item{pi}{Numeric: A vector of allele frequencies for each population.}

\item{ni}{Numeric: A vector of sample sizes (number of diploid individuals)
for each population.}

\item{r}{Integer: A single value, the number of populations.}

\item{hi}{Numeric: A vector of observed heterozygosities for each population.}
}
\value{
Returns a data table with three columns \code{$A}, \code{$B},
\code{$C}, which correspond to the 'a', 'b' and 'c' variance components
described in Weir & Cockheram (1984).
}
\description{
Takes a list of values of allele frequencies, sample sizes, number of
populations, and heterozygosity, and returns variance components, as
per Weir & Cockerham (1984). Assumes that all values are for a single
biallelic SNP locus. \cr\cr Note, this function is not exported and is an
internal function for \code{fstat_calc}.
}
\references{
Weir & Cockerham (1984) Evolution. DOI: 10.1111/j.1558-5646.1984.tb05657.x
Weir et al. (2002) Annals of Human Genetics. DOI: 10.1146/annurev.genet.36
}
\keyword{internal}
