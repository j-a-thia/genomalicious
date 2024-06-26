% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adegenet_DT2genX.R
\name{adegenet_DT2genX}
\alias{adegenet_DT2genX}
\title{Convert a long-format data table of genotypes into a genind/genlight object}
\usage{
adegenet_DT2genX(
  dat,
  genX,
  sampCol = "SAMPLE",
  locusCol = "LOCUS",
  genoCol = "GT",
  popCol = NULL
)
}
\arguments{
\item{dat}{Data.table: Contains samples and loci in rows, with a
separate column endcoding genotypes, where each allele is separated
by a '/', or as a count of ref alleles (0, 1 or 2). Assumes biallelic loci.
Code missing data as './.' for characters, and \code{NA} for Ref allele counts.}

\item{genX}{Character: One of \code{'genind'} or \code{'genlight'}.}

\item{sampCol}{Character: Column with sample information. Default = \code{'SAMPLE'}.}

\item{locusCol}{Character: Column with locus information. Default = \code{'LOCUS'}.}

\item{genoCol}{Character: Column with genotype information. Default = \code{'GT'}.}

\item{popCol}{Character: Column with population information. Optional. Default is NULL.}
}
\value{
A \code{genind} object, with the slot \code{pop} slot filled if argumnet
\code{popCol} is specified.
}
\description{
For a data table with samples and loci as rows (long-format data table),
create a genind/genlight object, as per the \code{adegenet} package
(Jombart 2011 Bioinformatics). Data must be biallelic SNPs.
}
\examples{
library(genomalicious)

data(data_Genos)

data_Genos

# Genind without and with populations
adegenet_DT2genX(data_Genos, genX='genind')
adegenet_DT2genX(data_Genos, genX='genind', popCol='POP')

# Genlight without and with populations
adegenet_DT2genX(data_Genos, genX='genlight')
adegenet_DT2genX(data_Genos, genX='genlight', popCol='POP')

}
\references{
Jombart (2008) adegenet 1.3-1: new tools for the analysis of genome-wide SNP data. Bioinformatics.
}
