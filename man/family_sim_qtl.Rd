% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/family_sim_qtl.R
\name{family_sim_qtl}
\alias{family_sim_qtl}
\title{Simulate quantitative trait from simulated families}
\usage{
family_sim_qtl(
  famGenos,
  numLoci = NULL,
  qtlLoci = NULL,
  additiveVar,
  environVar
)
}
\arguments{
\item{famGenos}{Data.table: Long-format, must contain the columns
\code{$SAMPLE} of sample IDs, \code{$LOCUS} of locus IDs, and \code{$GT} of
genotypes.}

\item{numLoci}{Integer: The number of loci contributing to the trait. These loci
are drawn from random in the dataset, \code{famGenos}. Therefore, it is
essential the \code{numLoci} is <= the total unique loci in \code{famGenos}.
Alternate parameterisation is to specify the names of the loci to use with \code{qtlLoci}.}

\item{qtlLoci}{Character: A character vector of loci contributing to the trait.
These loci must be present in the dataset, \code{famGenos}.
Alternate parameterisation is to specify the number of  loci to use with \code{numLoci}.}

\item{additiveVar}{Numeric: The additive genetic variance.}

\item{environVar}{Numeric: The environmental variance.}
}
\value{
A list with the indexes \code{$trait} and \code{$loci} is returned, containing
information on the traits and the underpinning loci.

The index \code{$trait} contains a data.table with the following columns:
\enumerate{
   \item \code{$SAMPLE}, the sample IDs.
   \item \code{$G}, the additive genetic values, with variance \emph{Va}.
   \item \code{$E}, the environmental values, with variance \emph{Ve}.
   \item \code{$P}, the phenoypic values, \emph{G} + \emph{E}.
}

The index \code{$loci} contains a data.table with the following columns:
\enumerate{
   \item \code{$LOCUS}, the locus IDs.
   \item \code{$FREQ}, the frequency of the focal allele.
   \item \code{$A.VAL}, the additive genetic value for this locus.
}
}
\description{
This function is designed to pair with the \code{family_sim_genos} function
to simulate a quantitative trait from a set of simulated families.
However, it could also be used to simulate a quantitative trait from observed data.
The function is currently only designed to simulate a single trait controlled
by \emph{n} loci with phenotypic variance, \emph{Vp}, the sum of the
additive, \emph{Va}, and environmental, \emph{Ve}, variances.
}
\details{
Note, the current implementation of this function generates a specific
additive + environmental variance for a trait, given a set of allele frequencies
and the number/identify of loci of interest. The phenotypes and genotypic values
simulated are only relevant to the input sample only.

A total of \emph{n} == \code{numLoci} loci from \code{famGenos} are drawn
at random to underpin the quantitative trait. It is assumed that all loci contribute
equally to the trait, i.e., \emph{Va}/\emph{n}. The simulation starts by first
fitting genotype AA with a genetic value of +1, Aa with a value of 0,
and aa with a value of -1. Because the additive genetic variance
is also dependent on alelle frequencies: \cr\cr
\emph{Va} = 2\emph{pq}[\emph{a} + \emph{d}(\emph{q} - \emph{p})]^2 \cr\cr
The genetic values are then modified relative to the allele frequencies to ensure
the total additivie genetic variance in the sample sum to the specified \emph{Va}.
The phenotype per individual is the sum of their genotypic values plus randomly
drawn environmental deviation (\code{rnorm(.., mean=0, sd=sqrt(Ve))}).
}
\examples{
library(genomalicious)
data(data_Genos)

# Subset Pop1 genotypes
genosPop1 <- data_Genos[POP=='Pop1', c('SAMPLE', 'LOCUS', 'GT')]

# Get the allele frequencies for Pop1
freqsPop1 <- genosPop1[, .(FREQ=sum(GT)/(length(GT)*2)), by=LOCUS]

# Simulate 100 familial relationships of each class
simFamily <- family_sim_genos(
   freqData=freqsPop1,
   locusCol='LOCUS',
   freqCol='FREQ',
   numSims=100,
   returnParents=TRUE,
   returnPedigree=TRUE
)

# Take a look at the focal pairs
simFamily$focal.pairs

# Take a look at the parentals
simFamily$parents

# Take a look at the pedigree
simFamily$pedigree

### THE OBSERVED GENOMIC RELATIONSHIPS MATRIX
library(AGHmatrix)

# A genotype matrix for the focal pairs
obsGenosMat <- genosPop1 \%>\% DT2Mat_genos()

# Calculate the GRM
obsGRM <- Gmatrix(obsGenosMat, method='Yang', ploidy=2)

### THE SIMULATED GENOMIC RELATIONSHIPS MATRIX
# Convert simulated families into a genotype matrix
simGenosMat <- DT2Mat_genos(simFamily$focal.pairs)

# Calculate the GRM
simGRM <- Gmatrix(simGenosMat, method='Yang', ploidy=2)

### COMPARE THE OBSERVED AND SIMULATED
relComp <- family_sim_compare(
   simGRM=simGRM,
   obsGRM=obsGRM,
   look='classic'
)

# The data
relComp$data

# Simulated dataset
relComp$data[!is.na(SIM)]

# The observed dataset
relComp$data[is.na(SIM)]

# Plot of relatedness values. Dashed lines denote relatedness
# values of 0, 0.0625, 0.125, 0.25, and 0.5, which are the theoretical
# expectations for unrelated individuals, half-cousins, cousins, half-siblings,
# and siblings/parent-offspring, respectively.
# You will note a large variance in the expected values, which
# is not surprising for this very small SNP dataset (200 loci).
relComp$plot

### SIMULATE A QUANTITATIVE TRAIT

# Combine the focal pairs and parentals, and simualte a trait controlled by
# 100 loci with Va = 1, and Ve = 1.
simQTL <- family_sim_qtl(
  famGenos=rbind(simFamily$focal.pairs, simFamily$parents),
  numLoci=100, additiveVar=1, environVar=1
  )

# The trait values
simQTL$trait

# The locus values
simQTL$loci

}
\references{
Falconer and MacKay (1996) Introduction to Quantitative Genetics, 3rd Ed., chapter 8, page 129.
}
