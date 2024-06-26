#' Generate \code{bayescan} input files
#'
#' Generates an input file for Bayescan (Foll & Gaggiotti, 2008) from pooled allele frequencies. \cr
#'
#' @param dat Data table: The biallelic SNP data. Requires all of the following columns: \cr
#'              \enumerate{
#'                \item The population pool ID. \cr
#'                \item The locus ID. \cr
#'                \item The population Ref allele frequency. \cr
#'              }
#'
#' @param pool.info Data table: The population pool metadata. Requires all of the following columns: \cr
#'              \enumerate{
#'                \item \code{$POOL} = The population pool ID. \cr
#'                \item \code{$INDS} = The number of diploid individuals in each pooled library.
#'              }
#'
#' @param file.bayescan Character: The name of the input file for \code{bayescan}/
#'
#' @param file.loci Character: A file that contains an information about the locus ID in
#' the \code{bayescan} input file (i.e., \code{file.bayescan}).
#'
#' @param poolCol Character: Population pool ID. Default = \code{'POOL'}
#'
#' @param locusCol Character: Locus ID. Default = \code{'LOCUS'}
#'
#' @param freqCol Character: The reference allele frequency. Default = \code{'FREQ'}.
#'
#' @details The allele counts in the Bayescan input file generated reflect the number of haploid
#' genomes pooled. E.g. if 20 individuals were pooled, i.e. 40 genomes, and the Ref allele was
#' estimated at a frequency of 0.7, the counts would be Ref=28 and Alt=12. NOTE: The estimated number
#' of indivdiuals carrying an allele is rounded to nearest integer (e.g. 1.5 = 2, and 1.4 = 1),
#' with the exception when the number of individauls is < 1 but > 0, in which it is always rounded to 1.
#'
#' @references Foll & Gaggiotti (2008) A genome scan method to identify selected loci appropriate
#' for both dominant and codominant markers: A Bayesian perspective. Genetics 180: 977-993.
#'
#' @examples
#' library(genomalicious)
#'
#' data(data_PoolFreqs)
#' data(data_PoolInfo)
#'
#' bayescan_inputs_pools(dat=data_PoolFreqs,
#' pool.info=data_PoolInfo,
#' file.bayescan='Bayescan_input.txt',
#' file.loci='Bayescan_loci.txt',
#' poolCol='POOL',
#' locusCol='LOCUS',
#' freqCol='FREQ')
#'
#' @export
bayescan_inputs_pools <- function(dat, pool.info, file.bayescan, file.loci
                                 , poolCol, locusCol, freqCol) {

  # BEGIN ...........
  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  for(lib in c('data.table')){ require(lib, character.only=TRUE) }

  # Check class of dat.
  if(!'data.table' %in% class(dat)){ stop("Argument dat isn't a data table.")}

  # Check class of pool.info.
  if(!'data.table' %in% class(pool.info)){ stop("Argument pool.info isn't a data table") }

  # Test for the necessary columns in dat.
  if(sum(c(poolCol, locusCol, freqCol) %in% colnames(dat))!= 3){
    stop("Argument dat needs the columns specified by the arguments: poolCol, locusCol, freqCol")
  }

  # Test for the necessary columns in pool.info.
  if(sum(c('POOL', 'INDS') %in% colnames(pool.info))!=2){
    stop("Argument pool.info needs the columns $POOL and $INDS.")
  }

  # --------------------------------------------+
  # Code
  # --------------------------------------------+
  # Reassign population, locus, and frequency columns
  colReass <- match(c(poolCol, locusCol, freqCol), colnames(dat))
  colnames(dat)[colReass] <- c('POOL', 'LOCUS', 'P')

  # Create some character objects to insert into input file
  num.loci <- paste0('[loci]=',length(unique(dat$LOCUS)))
  num.pops <- paste0('[populations]=',length(unique(dat$POOL)))

  # Reduce data columns then split the data on $POOL.
  dat.spl <- lapply(split(dat[,c('POOL','LOCUS','P')], dat$POOL), function(X){
    setorder(X, LOCUS)
    return(X)
  })

  # Iterate through each Xth population and make a Bayescan-friendly data table.
  # The observed REF and ALT alleles counts are derived from the estimated values in $P,
  # with respect to the sampled number of genomes (2*diploid individuals).
  BS.ls <- lapply(dat.spl, function(X){
    pool <- X$POOL[1]

    X$INDS <- pool.info[POOL==pool]$IND

    X$REF.COUNT <- apply(X[, c('P', 'INDS')], 1, function(Y){
                    p <- Y[['P']]
                    inds <- Y[['INDS']]
                    ref.count <- p * (inds * 2)
                    if(p != 0 & ref.count < 1){ ref.count <- 1
                    } else if(p != 1 & ref.count < 1){ ref.count <- 1
                    } else{ ref.count <- round(ref.count)
                    }
                    return(ref.count)
                  })

    X[, ALT.COUNT:=(INDS*2)-REF.COUNT]

    # BS.dt is a data table for each population in the Bayescan (BS) format.
    # There are 5 columns: 1=the locus number; 2=total sampled genomes (2xdiploid individuals);
    # 3=the number of alleles at the locus (which is 2 because SNPs are biallelic);
    # 4=Ref allele count; 5=Alt allele count.
    BS.dt <-data.table(MARKER=1:nrow(X), LOCUS=X$LOCUS, INDS=X$INDS, ALLELES=2
                       , REF=X$REF.COUNT, ALT=X$ALT.COUNT)
    return(BS.dt)
  })

  # Write (append) lines to txt file.
  # First open the file
  # Next append the number of loci and number of pops
  # Then iterate through each j-th population in BS.ls and write data to file
  file.create(file.bayescan)

  for(i in c(num.loci,'',num.pops,'')){ write(i,file=file.bayescan,append=TRUE) }

  for(j in 1:length(names(BS.ls))){
    write(paste0('[pop]=',j),file=file.bayescan,append=TRUE)
    write.table(BS.ls[[j]][,-'LOCUS'],file=file.bayescan,append=TRUE,col.names=FALSE,row.names=FALSE)
    write('',file=file.bayescan,append=TRUE)
  }

  # Create a file of marker/locus order
  fwrite(x=BS.ls[[1]][,c('MARKER','LOCUS')], file=file.loci, sep='\t')

  # ........... END
}
