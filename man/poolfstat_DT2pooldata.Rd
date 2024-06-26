% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poolfstat_DT2pooldata.R
\name{poolfstat_DT2pooldata}
\alias{poolfstat_DT2pooldata}
\title{Convert a data table of read counts into a pooldata object}
\usage{
poolfstat_DT2pooldata(dat, flip = FALSE, poolInfo)
}
\arguments{
\item{dat}{Data table: A long-format data table of read counts for pool-seq
samples of biallelic SNP loci to be converted into a pooldata object.
Must contain all the following columns:
\enumerate{
   \item \code{$CHROM} The chromosome (contig) ID.
   \item \code{$POS} The variant position on the chromosome.
   \item \code{$REF} The reference allele.
   \item \code{$ALT} The alternate allele.
   \item \code{$LOCUS} The locus ID.
   \item \code{$POOL} The pool ID.
   \item \code{$AO} The number of reads supporting the alternate allele.
   \item \code{$RO} The number of reads supporting the reference allele.
}
Alternatively, a pooldata object to be converted into a data table.
See param \code{flip} with respect to specifying the direction of the conversion.}

\item{flip}{Logical: Should the function be reversed? Default is \code{FALSE}.
When \code{flip==FALSE}, a data table is converted into a pooldata object.
When \code{flip==TRUE}, a pooldata object is converted into a data table.}

\item{poolInfo}{Data table: Contains the sample sample sizes (number of diploids) for
for each unique pool listed in \code{dat$POOL}. Requires two columns:
\enumerate{
   \item \code{$POOL} The pools listed in \code{dat$POOL}.
   \item \code{$INDS} The number of diploid individuals for the pools.
}
You only need to specify \code{poolInfo} when \code{flip==FALSE}, that is,
transforming a data table itno a pooldata object.}
}
\value{
When \code{flip==FALSE}, returns a pooldata object as per
the \code{poolfstat} package. Alternatively, when \code{flip==TRUE}, returns
a data table.
}
\description{
Used to takes a data table of pool-seq read counts and creates an object of class
"pooldata", as per the \code{poolfstat} package (Hivert et al. 2018).
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

# Make pooldata object
pooldataObj <-poolfstat_DT2pooldata(
   data_PoolFreqs,
   flip=FALSE,
   poolInfo=data_PoolInfo
)

class(pooldataObj)

pooldataObj

# And go back to data table
pooldataTab <- poolfstat_DT2pooldata(pooldataObj, flip=TRUE)

head(pooldataTab)

}
\references{
Hivert et al. (2018) Genetics. DOI: 10.1534/genetics.118.300900
}
