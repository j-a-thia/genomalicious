% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_space_loci.R
\name{filter_space_loci}
\alias{filter_space_loci}
\title{Filter loci based on their spacing}
\usage{
filter_space_loci(
  dat,
  chromCol = "CHROM",
  posCol = "POS",
  locusCol = "LOCUS",
  stepSize
)
}
\arguments{
\item{dat}{Data.table: Contains the information on the genomic context of each
locus, that is, their position and which chromosome/contig they reside on. You can
pass this function a genotype data.table (e.g., as produced from \code{vcf2DT()}),
because it will subset only the unique chromosome, position, and locus information.
Must contain the columns:
\enumerate{
   \item The chromosome/contig ID (see param \code{chromCol}).
   \item The positional information (see param \code{posCol}).
   \item The locus ID (see param \code{locusCol}).
}}

\item{chromCol}{Character: The column name with the chromosome information.
Default = \code{'CHROM'}.}

\item{posCol}{Character: The column name with the position information.
Default = \code{'POS'}.}

\item{locusCol}{Character: The column name with the locus name information.
Default = \code{'LOCUS'}.}

\item{stepSize}{Integer: the size of steps between loci.}
}
\value{
Returns a vector of loci IDs to be kept.
}
\description{
A function to space loci based on a particular step size.
For each chromosome (contig), starting from the first locus, each
subsequent locus is assessed with respect to whether it is further
than the step size from the previous locus. If a locus is not >= the
step size, it will be skipped, and the next locus evaluated.
}
