% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poolfstat_DT2Fst.R
\name{poolfstat_DT2Fst}
\alias{poolfstat_DT2Fst}
\title{Calculate FST with \code{poolfstat} from a data table of read counts}
\usage{
poolfstat_DT2Fst(dat, pool.info, method = "Anova")
}
\arguments{
\item{dat}{Data table: Contains read counts, e.g. like that been
produced by the function \code{vcf2DT}. Must contain all the following columns:
\enumerate{
   \item \code{$CHROM} The chromosome (contig) ID.
   \item \code{$POS} The variant position on the chromosome.
   \item \code{$REF} The reference allele.
   \item \code{$ALT} The alternate allele.
   \item \code{$POOL} The pool ID.
   \item \code{$AO} The number of reads supporting the alternate allele.
   \item \code{$RO} The number of reads supporting the reference allele.
}}

\item{pool.info}{Data table: Contains the sample sample sizes (number of diploids) for
for each unique pool listed in \code{dat$POOL}. Requires two columns:
\enumerate{
   \item \code{$POOL} The pools listed in \code{dat$POOL}.
   \item \code{$INDS} The number of diploid individuals for the pools.
}}

\item{method}{Character: Either 'Anova' (default) or 'Identity'. Passed to \code{method} argument
in \code{poolfstat::computeFST}.}
}
\value{
Returns a list with two indices: \code{$Fst} is the calculated FST among the
pools using a function call of \code{poolfstat::computeFST}, whereas \code{$pooldat} is the
\code{poolfstat} object used to generate said FST values.
}
\description{
Takes a data table of read counts and creates an object of class
\code{poolfstat}. The FST for the pools in the data table is calculated using
the function \code{poolfstat::computeFST}. Also requires pool size information.
}
\examples{
library(genomalicious)

# Load in the pool metadata and reads
data(data_PoolInfo)
data(data_PoolFreqs)

# Pool info
data_PoolInfo

# Pool reads in $DP, $AO, and $RO
data_PoolFreqs[, c('POOL','DP','AO','RO')]

# Calculate FST using poolfstat
Y <- poolfstat_DT2Fst(data_PoolFreqs, data_PoolInfo)

# Output is a list
class(Y)

# Outout from poolfstat::computeFST
Y$Fst

# The pooldata class object, generated from data table of pooled reads
class(Y$pooldat)
Y$pooldat

}
