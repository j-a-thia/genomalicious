#' Merge \code{poolne_estim} outputs
#'
#' Searches a user specified directory for \code{summary_pi.out} or \code{summary_ne_eps.out} files,
#' which are the output files from \code{poolne_estim} containing estimates of population allele frequencies
#' and effective pool size (respectively), merging these outputs into a single data.table.
#'
#' If merging \code{summary_pi.out} files, \code{Loci.txt} files are also required, which have
#' the names of the loci in the same order as the rows in the \code{summary_pi.out} file.
#' The \code{Loci.txt} files are produced by the function \code{poolne_estim_input}. \cr
#' \cr
#' Files are expected to follow the naming convention: \code{[runID]_[poolID]_summary_pi.out},
#' \code{[runID]_[poolID]_summary_ne_eps.out} and \code{[runID]_[poolID]_Loci.txt}.
#' The \code{[runID]} \code{[poolID]} portion is used to assign run ID and the population pool ID, respectively.
#'
#' @usage poolne_estim_output(stat=c('pi', 'ne'), datDir=NA, lociDir=NA)
#'
#' @param stat Character: One of 'pi' (population allele frequency estimates) or 'ne' (effective pool
#' size estimates). Default = NA.
#'
#' @param datDir Character: A directory to search for \code{summary_pi.out} or \code{summary_ne_eps.out} files.
#' Default = NA. If unspecified, function assumes files are in the current working directory.
#'
#' @param lociDir Character: A directory to search for \code{Loci.txt} files. Default = NA.
#' If unspecified, function assumes files are in the current working directory. You do not need
#' to specify this if \code{stat = 'pi'}.
#'
#' @return If \code{stat = 'pi'}, a data.table with the following columns is returned: \cr
#' \enumerate{
#'   \item \code{$MRK} = The marker number, as used in \code{poolne_estim}.
#'   \item \code{$PI} = The \code{poolne_estim} posterior mean of the population Ref allele
#'     frequency, pi.
#'   \item \code{$SD} = The \code{poolne_estim} estimate of the standard deviation in \code{PI}.
#'   \item \code{$POOL} = The population pool ID.
#'   \item \code{$RUN} = The run ID.
#'   \item \code{$CHROM} = The chromosome ID.
#'   \item \code{$LOCUS} = The locus ID.
#' } \cr
#' If \code{stat = 'ne'}, a data.table with the following columns is returned: \cr
#' \enumerate{
#'   \item \code{$POOL} = The population pool ID.
#'   \item \code{$RUN} = The run ID.
#'   \item \code{$SAMPLE} = The replicate sample ID, as used in \code{poolne_estim}.
#'   \item \code{$NEuntr} = The untruncated \code{poolne_estim} posterior mean effective pool-size, ne.
#'   \item \code{$SD.NEuntr} = The standard deviation of \code{NEraw}.
#'   \item \code{$NE} = The truncated \code{poolne_estim} posterior mean effective pool-size, ne.
#'   \item \code{$SD.NE} = The standard deviation of \code{NE}.
#'   \item \code{$E} = The \code{poolne_estim} posterior mean of the experimental error, epsilon.
#'   \item \code{$SD.E} = The standard deviation of \code{E}.
#'   \item \code{$PER} = Percent acceptance.
#'   \item \code{$ACCRATE} = The acceptance rate of the chain.
#' }
#'
#' @examples
#' # Create a link to raw external datasets in genomalicious
#' genomaliciousExtData <- paste0(find.package('genomalicious'), '/extdata')
#'
#' # Use list.files() to show the summary_pi.out, Loci.txt and
#' # summary_ne_eps.out text files
#' list.files(genomaliciousExtData, pattern='summary_pi.out')
#' list.files(genomaliciousExtData, pattern='Loci.txt')
#' list.files(genomaliciousExtData, pattern='summary_ne_eps.out')
#'
#' # Merge the outputs
#' pi.data <- poolne_estim_output(stat='pi', datDir=genomaliciousExtData, lociDir=genomaliciousExtData)
#' ne.data <- poolne_estim_output(stat='ne', datDir=genomaliciousExtData)
#'
#' @export
poolne_estim_output <- function(stat=NA, datDir=NA, lociDir=NA){
  # BEGIN ............
  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  require(data.table)

  # Check that stat is one of specified values.
  if((stat %in% c('pi', 'ne')==FALSE)){
    stop("Argument stat must take the value of 'pi' or 'ne'.")
  }

  if(is.na(datDir)==TRUE){ datDir <- getwd()}
  if(is.na(datDir)==TRUE){ lociDir <- getwd()}

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  if(stat=='pi'){
    pi.ls <- list.files(path=datDir, pattern='summary_pi.out')
    loci.ls <- list.files(path=lociDir, pattern='Loci.txt')

    # Make a datatable, merging the 'summary_pi.out' files.
    allPools.dt <- data.table(matrix(data=NA, nrow=0, ncol=5))
    colnames(allPools.dt) <- c('MRK', 'PI_M', 'PI_SD', 'POOL', 'RUN')

    # Iterate through the pi files
    for(pi in pi.ls){
      # Load pool
      pi.dt <- fread(paste0(datDir,'/',pi))
      # If there are >9999 markers, poolne_etim doesn't record these; need to manually replace.
      pi.dt$MRK <- 1:nrow(pi.dt)
      # Extract pool name and run ID
      pi.dt$POOL <- strsplit(pi,'_',fixed=TRUE)[[1]][2]
      pi.dt$RUN <- strsplit(pi,'_',fixed=TRUE)[[1]][1]
      # Bind to full dataset
      allPools.dt <- rbind(allPools.dt, pi.dt)
    }

    # DO a quick readjustment of colnames for allPools.dt
    colnames(allPools.dt) <- c('MRK', 'PI', 'SD', 'POOL', 'RUN')

    # Make a datatable, merging the 'Loci.txt' files.
    allLoci.dt <- data.table(matrix(data=NA, nrow=0, ncol=5))
    colnames(allLoci.dt) <- c('CHROM','LOCUS','MRK','POOL', 'RUN')

    # Iterate through the loci files
    for(loci in loci.ls){
      loci.dt <- fread(paste0(lociDir,'/',loci))
      loci.dt$MRK <- 1:nrow(loci.dt)
      loci.dt$POOL <- strsplit(loci,'_',fixed=TRUE)[[1]][2]
      loci.dt$RUN <- strsplit(pi,'_',fixed=TRUE)[[1]][1]
      allLoci.dt <- rbind(allLoci.dt, loci.dt)
    }

    # Match LOCUS/MRK in allLoci.dt to MRK in allPools.dt
    # First create a list of POOL names
    pools <- as.list(unique(allPools.dt$POOL))
    # Then iterate through each Xth POOL name
    # out.ls contains individual datatables for each Xth POOL
    out.ls <- lapply(pools, function(X){
      # Subset allLoci.dt and allPools.dt by Xth POOL
      loci.sub <- subset(allLoci.dt, POOL==X)
      allpools.sub <- subset(allPools.dt, POOL==X)
      # Match values of CHROM/LOCUS by MRK in allpools.sub with MRK in
      # in loci.sub. CHROM/LOCUS values are appended to allpools.sub, which is returned.
      allpools.sub <- cbind(allpools.sub
                            ,loci.sub[,c('CHROM','LOCUS')][match(allpools.sub$MRK,loci.sub$MRK)])
      return(allpools.sub)
    })

    # rbind the datatables in out.ls, and return the complete datatable.
    return(do.call('rbind', out.ls))
  }

  if(stat=='ne'){
    ne.ls <- list.files(path=datDir, pattern='summary_ne_eps.out')

    # Make a data.table, merging the 'summary_ne_eps.out' files.
    allPools.dt <- data.table(matrix(data=NA, nrow=0, ncol=11))
    colnames(allPools.dt) <- c('POOL', 'RUN', 'SAMPLE', 'NEuntr', 'SD.NEuntr', 'NE', 'SD.NE', 'E', 'SD.E'
                               , 'PER', 'ACCRATE')

    # Iterate through the pi files
    for(ne in ne.ls){
      # Load pool
      ne.dt <- fread(paste0(datDir, '/', ne))
      # Extract pool name and run ID
      ne.dt <- data.table(POOL=strsplit(ne, '_', fixed=TRUE)[[1]][2]
                          , RUN=strsplit(ne, '_', fixed=TRUE)[[1]][1]
                          , ne.dt
      )
      # Correct column names
      colnames(ne.dt) <- c('POOL', 'RUN', 'SAMPLE', 'NEuntr', 'SD.NEuntr', 'NE', 'SD.NE', 'E', 'SD.E'
                           , 'PER', 'ACCRATE')
      # Bind to full dataset
      allPools.dt <- rbind(allPools.dt, ne.dt)
    }

    # Return data
    return(allPools.dt)
  }
  # ............ END
}

