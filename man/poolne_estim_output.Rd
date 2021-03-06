% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poolne_estim_output.R
\name{poolne_estim_output}
\alias{poolne_estim_output}
\title{Merge \code{poolne_estim} outputs}
\usage{
poolne_estim_output(stat=c('pi', 'ne'), datDir=NA, lociDir=NA)
}
\arguments{
\item{stat}{Character: One of 'pi' (population allele frequency estimates) or 'ne' (effective pool
size estimates). Default = NA.}

\item{datDir}{Character: A directory to search for \code{summary_pi.out} or \code{summary_ne_eps.out} files.
Default = NA. If unspecified, function assumes files are in the current working directory.}

\item{lociDir}{Character: A directory to search for \code{Loci.txt} files. Default = NA.
If unspecified, function assumes files are in the current working directory. You do not need
to specify this if \code{stat = 'pi'}.}
}
\value{
If \code{stat = 'pi'}, a data.table with the following columns is returned: \cr
\enumerate{
  \item \code{$MRK} = The marker number, as used in \code{poolne_estim}.
  \item \code{$PI} = The \code{poolne_estim} posterior mean of the population Ref allele
    frequency, pi.
  \item \code{$SD} = The \code{poolne_estim} estimate of the standard deviation in \code{PI}.
  \item \code{$POOL} = The population pool ID.
  \item \code{$RUN} = The run ID.
  \item \code{$CHROM} = The chromosome ID.
  \item \code{$LOCUS} = The locus ID.
} \cr
If \code{stat = 'ne'}, a data.table with the following columns is returned: \cr
\enumerate{
  \item \code{$POOL} = The population pool ID.
  \item \code{$RUN} = The run ID.
  \item \code{$SAMPLE} = The replicate sample ID, as used in \code{poolne_estim}.
  \item \code{$NEuntr} = The untruncated \code{poolne_estim} posterior mean effective pool-size, ne.
  \item \code{$SD.NEuntr} = The standard deviation of \code{NEraw}.
  \item \code{$NE} = The truncated \code{poolne_estim} posterior mean effective pool-size, ne.
  \item \code{$SD.NE} = The standard deviation of \code{NE}.
  \item \code{$E} = The \code{poolne_estim} posterior mean of the experimental error, epsilon.
  \item \code{$SD.E} = The standard deviation of \code{E}.
  \item \code{$PER} = Percent acceptance.
  \item \code{$ACCRATE} = The acceptance rate of the chain.
}
}
\description{
Searches a user specified directory for \code{summary_pi.out} or \code{summary_ne_eps.out} files,
which are the output files from \code{poolne_estim} containing estimates of population allele frequencies
and effective pool size (respectively), merging these outputs into a single data.table.
}
\details{
If merging \code{summary_pi.out} files, \code{Loci.txt} files are also required, which have
the names of the loci in the same order as the rows in the \code{summary_pi.out} file.
The \code{Loci.txt} files are produced by the function \code{poolne_estim_input}. \cr
\cr
Files are expected to follow the naming convention: \code{[runID]_[poolID]_summary_pi.out},
\code{[runID]_[poolID]_summary_ne_eps.out} and \code{[runID]_[poolID]_Loci.txt}.
The \code{[runID]} \code{[poolID]} portion is used to assign run ID and the population pool ID, respectively.
}
\examples{
# Create a link to raw external datasets in genomalicious
genomaliciousExtData <- paste0(find.package('genomalicious'), '/extdata')

# Use list.files() to show the summary_pi.out, Loci.txt and
# summary_ne_eps.out text files
list.files(genomaliciousExtData, pattern='summary_pi.out')
list.files(genomaliciousExtData, pattern='Loci.txt')
list.files(genomaliciousExtData, pattern='summary_ne_eps.out')

# Merge the outputs
pi.data <- poolne_estim_output(stat='pi', datDir=genomaliciousExtData, lociDir=genomaliciousExtData)
ne.data <- poolne_estim_output(stat='ne', datDir=genomaliciousExtData)

}
