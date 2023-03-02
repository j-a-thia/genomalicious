#' @title  Calculate F-statistics from genotypes or allele frequencies
#'
#' @description Function takes a genotypes or allele frequencies in a long-format data table
#' and calculates Weir & Cockerham's F-statistics (Weir & Cockerham, 1984).
#' Permutations can be used to test statistical significance of FST in
#' genotype data sets. See Details for more information.
#'
#' @param dat Data table: A long-format data table of biallelic genotypes,
#' coded as '/' separated alleles ('0/0', '0/1', '1/1'), or counts
#' of the Alt alleles (0, 1, 2, respectively). Alternatively, a long-format
#' data table of allele frequencies.
#' Columns required for \strong{both} genotypes and allele frequencies:
#' \enumerate{
#'    \item The population ID (see param \code{popCol}).
#'    \item The locus ID (see param \code{locusCol}).
#' }
#' Columns required only for \strong{genotypes}:
#' \enumerate{
#'    \item The sample ID (see param \code{sampCol}).
#'    \item The genotypes (see param \code{genoCol}).
#' }
#' Columns required only for \strong{allele frequencies}:
#' \enumerate{
#'    \item The allele frequencies (see param \code{freqCol}).
#'    \item The number of individuals used to obtain the allele frequency
#'    estimate (see param \code{indsCol}).
#' }
#'
#' @param type Character: One of \code{'genos'} or \code{'freqs'}, to calculate
#' F-statistics from genotype or allele frequency data, respectively.
#'
#' @param fstat Character: A vector of F-statistics to calculate. This is only
#' applicable for genotype data. Must include one of \code{"FST"},
#' \code{"FIS"}, or \code{"FIT"}.
#'
#' @param popCol Character: The column name with the population information.
#' Default is \code{'POP'}.
#'
#' @param sampCol Character: The column name with the sampled individual information.
#' Default is \code{'SAMPLE'}.
#'
#' @param locusCol Character: The column name with the locus information.
#' Default is \code{'LOCUS'}.
#'
#' @param genoCol Character: The column name with the genotype information.
#' Default is \code{'GT'}.
#'
#' @param freqCol Character: The column name with the allele freuqency information.
#' Default is \code{'FREQ'}.
#'
#' @param indsCol Character: The column name with the number of individuals
#' contributing to the allele freuqency estimate. Default is \code{indsCol}.
#'
#' @param global Logical: Calculate global F-statistics across all populations?
#' Default is \code{TRUE}. Note, \code{global==TRUE} and \code{pairwise==TRUE}
#' is not permitted.
#'
#' @param pairwise Logical: Calculate pairwise FST between populations?
#' Default is \code{FALSE}. Note, \code{global==TRUE} and \code{pairwise==TRUE}
#' is not permitted.
#'
#' @param permute Logical: Should permutations be performed to test the
#' statistical significance of F-statistics? Default is \code{FALSE}.
#' Can only be performed on genotype data, i.e., \code{type=='genos'}.
#'
#' @param numPerms Integer: The number of permutations to perform. Default is 100.
#'
#' @param numCores Integer: The number of cores to use when running permutations. Default is 1.
#'
#' @details With genotype data, the F-statistics FST, FIS, and FIT can be calculated.
#' Only FST can be calculated from allele frequency data.
#'
#' F-statistics from genotype data are calculated from the variance components
#' 'a', 'b', and 'c', which have been standardised for observed heterozygosity.
#' FST from allele frequencies is calculated using estimates of MSP
#' (mean squares for populations), MSG (mean squares for gametes) and
#' nc (the weighted mean sample size).
#'
#' Permutation tests involve random shuffling of individuals
#' among populations, recalculating F-statistics, and testing the hypothesis
#' that the permuted F-statistic > observed F-statistic. The p-value represents the
#' proportion of permutation that were TRUE to this expression. That is,
#' if no permuted values are greater than the observed, p=0. Likewise, if all
#' the permuted values are greater than the observed, p=1.
#'
#' @return A list is returned with three indexes.
#'
#' The first index is \code{$genome}, the genome-wide F-statistics. If global
#' estimates were requested, \code{global==TRUE}, then this is just a single row;
#' the estimates across all populations. If pairwise esimates were requested,
#' \code{pairwise==TRUE}, then there are \code{$POP1} and \code{$POP2}, which
#' represent two populations tested.
#'
#' The second index is \code{$locus}, the locus-specific F-statistics. This is
#' a data table with a \code{$LOCUS} column for global estimates at each locus.
#' when \code{global==TRUE}. If pairwise population estimates have been requested,
#' \code{pairwise==TRUE}, then there are \code{$POP1} and \code{$POP2}, which
#' represent the two populations tested.
#'
#' The third index is \code{$permute}, the permutation results. This index will
#' be \code{NULL} when frequencies are used, i.e., \code{type=='freqs'}, and
#' will only contain data if \code{type=='genos'} and \code{permute==TRUE}.
#' \code{$permute} is itself a list, with two subindexes:
#' \enumerate{
#'    \item \code{$fstat}: The permuted F-statistics. If \code{global==TRUE}, then
#'    this will simply be a single row of global estimates. If \code{pairwise==TRUE},
#'    then this will be a data table with columns \code{$POP1}, \code{$POP2}, and
#'    a column for each F-statistic.
#'
#'    \item \code{$pval}: The permuted p-values. This is a long-format data table.
#'    If \code{global==TRUE}, then there are two column: \code{$STAT}, which
#'    contains the F-statistic; and \code{$PVAL}, which contains the global
#'    permuted p-value. If \code{pairwise==TRUE}, then there will two additional
#'    columns, \code{$POP1} and \code{$POP2}.
#' }
#'
#' @references Weir & Cockerham (1984) Evolution. DOI: 10.1111/j.1558-5646.1984.tb05657.x
#' Weir et al. (2002) Annals of Human Genetics. DOI: 10.1146/annurev.genet.36
#'
#' @examples
#' library(genomalicious)
#'
#' data(data_Genos)
#' data(data_PoolFreqs)
#' data(data_PoolInfo)
#'
#' ### Genotypes, global estimates for FST, FIS and FIT, with 50 permutations
#' genoGlobal <- fstat_calc(
#'    data_Genos,
#'    type='genos',
#'    fstat=c('FST','FIS','FIT'),
#'    global=TRUE,
#'    permute=TRUE,
#'    numPerms=50)
#'
#' # Genome-wide estimate
#' genoGlobal$genome
#'
#' # Locus specific estimates
#' genoGlobal$locus %>% head
#'
#' # Permutations
#' genoGlobal$permute$fstat
#' genoGlobal$permute$pval
#'
#' ### Genotypes, pairwise estimates for FST, FIS and FIT, no permutations
#' genoPairs <- fstat_calc(
#'    data_Genos,
#'    type='genos',
#'    fstat=c('FST','FIS','FIT'),
#'    global=FALSE,
#'    pairwise=TRUE
#'    )
#'
#' # Genome-wide estimates
#' genoPairs$genome %>% head
#'
#' # Locus specific estimates
#' genoPairs$locus %>% head
#'
#' # Permutations (there are none!)
#' genoPairs$permute
#'
#' ### Genotypes, pairwise, FST only, with 10 permutations
#' genoPairsPerms <- fstat_calc(
#'    data_Genos,
#'    type='genos',
#'    fstat='FST',
#'    global=FALSE,
#'    pairwise=TRUE,
#'    permute=TRUE,
#'    numPerms=10
#'    )
#'
#' # Permuted FST
#' genoPairsPerms$permute$fstat %>% head
#'
#' # Permuted p-value
#' genoPairsPerms$permute$pval %>% head
#'
#' ### Frequencies, pairwise (can only do FST)
#' # Note columns in data_PoolFreqs
#' colnames(data_PoolFreqs)
#'
#' # We need to add in the number of diploid individuals, $INDS
#' newFreqData <- left_join(data_PoolFreqs, data_PoolInfo)
#' head(newFreqData)
#'
#' # There is no $POP column, but there is a $POOL column that contains
#' # the population pool ID. Need to manually specify this.
#' freqPairs <- fstat_calc(
#'    newFreqData,
#'    type='freqs',
#'    global=FALSE,
#'    pairwise=TRUE
#' )
#'
#' freqPairs
#'
#' @export

fstat_calc <- function(
    dat, type, fstat=NULL, popCol='POP', sampCol='SAMPLE', locusCol='LOCUS',
    genoCol='GT', freqCol='FREQ', indsCol='INDS', global=TRUE, pairwise=FALSE,
    permute=FALSE, numPerms=100, numCores=1
  ){

  # --------------------------------------------+
  # Assertions and environment
  # --------------------------------------------+
  require(data.table); require(tidyverse)

  # Make sure dat is a data.table
  dat <- as.data.table(dat)

  # Has type been specified correctly?
  if(!type %in% c('genos','freqs')){
    stop('Argument `type` must be one of "genos" or "freqs".')
  }

  # Check that global and pairwise are logicals
  if(class(global)!='logical'){
    stop('Argument `global` must be a logical value. See ?fstat_calc.')
  }

  if(class(pairwise)!='logical'){
    stop('Argument `pairwise` must be a logical value. See ?fstat_calc.')
  }

  # Cannot specify bpth global and pairwise
  if(global==TRUE & pairwise==TRUE){
    stop('Arguments `global` and `pairwise` cannot both be TRUE. See ?fstat_calc.')
  }

  # Must specify one of global or pairwise
  if(global==FALSE & pairwise==FALSE){
    stop('One of arguments `global` and `pairwise` must be TRUE. See ?fstat_calc.')
  }

  # Check that permute is logical
  if(class(permute)!='logical'){
    stop('Argument `permute` must be a logical value. See ?fstat_calc.')
  }

  # If permute is specified, then check numPerms is >0 and is an integer
  if(permute==TRUE){
    if(numPerms<1){
      stop('Argument `numPerms` must be an integer >0. See ?fstat_calc.')
    }

    numPerms <- as.integer(numPerms)
  }

  # Checks specific to each type
  if(type=='genos'){
    # Check all necessary columns are present and reassign values.
    if(sum(c(popCol,sampCol,locusCol,genoCol) %in% colnames(dat))!=4){
      stop(
        'Argument `dat` requires columns `popCol`, `sampCol`, `locusCol` and
        `genoCol` to estimate F-statistics from genotypes. See ?fstat_calc.')
    }

    # If genotypes are characters, convert to integers.
    if(class(dat$GT)=='character'){
      dat$GT <- genoscore_converter(dat$GT)
    }

    # If no F-statistics have been specified for genotypes, stop.
    if(is.null(fstat)){
      stop(
      'Argument `fstat` cannot be NULL when requesting F-statistics
      from genotypes. See ?fstat.')
    }

    # Make sure the specified value in fstat are correct
    if(sum(!fstat %in% c('FST','FIS','FIT'))!=0){
      stop(
        'Argument `fstat` must only contain "FST", "FIS", and/or "FIT" when
        requesting F-statistics from genotypes. See ?fstat.'
      )
    }
  }

  if(type=='freqs'){
    if(sum(c(popCol,locusCol,freqCol,indsCol) %in% colnames(dat))!=4){
      stop(
        'Argument `dat` requires columns `popCol`, `locusCol`, `freqCol`, and
        `indsCol` to estimate F-statistics from allele frequencies. See ?fstat_calc.')
    }

    if(is.null(fstat)!=TRUE){
      stop(
        'Argument `fstat` must be NULL when requesting FST from allele
        frequencies. See ?fstat_calc.'
      )
    }
  }

  # --------------------------------------------+
  # Code for genotypes
  # --------------------------------------------+

  if(type=='genos'){
    # Rename columns
    col.index <- match(c(popCol,sampCol,locusCol,genoCol), colnames(dat))
    colnames(dat)[col.index] <- c('POP','SAMPLE','LOCUS','GT')

    # START: FST among all populations
    if(global==TRUE){
      cat('F-statistic calculation on genotype data, global estimate', '\n')

      # Basic result
      result <- fstat_subfun(dat, type='genos', fstat=fstat)

      # For genotype permutations
      if(permute == TRUE){
        cat('Performing permutations', '\n')
        # Permuted global F-statistics
        result$permute$fstat <- fstat_permute(
          dat=dat, fstat=fstat, numPerms=numPerms, numCores=numCores
          )
        # Permuted p-values
        result$permute$pval <- lapply(fstat, function(fx){
          p <- sum(result$permute$fstat[[fx]] > result$permute$genome[[fx]]) / numPerms
          data.table(STAT=fx, PVAL=p)
        }) %>%
          do.call('rbind', .)
      } else{
        result$permute <- NULL
      }
      # END: FST among all populations
    }

    # START: FST between pairwise populations
    if(pairwise==TRUE){
      cat('F-statistic calculation on genotype data, pairwise estimates', '\n')
      pair.comp <- combn(unique(dat$POP), 2) %>%
        t() %>%
        as.data.table() %>%
        setnames(., new=c('POP1','POP2'))

      num.pairs <- nrow(pair.comp)

      result <- list(genome=list(), locus=list(), permute=list(fstat=NULL, pval=NULL))

      # START: Iterate over i pairs
      for(i in 1:nrow(pair.comp)){
        cat('Estimates for pair:', i, '/', num.pairs, '\n')
        pop1 <- pair.comp$POP1[i]
        pop2 <- pair.comp$POP2[i]

        # Calculate genome-wide and locus-wise F-statistics
        pair.calc <- dat[POP %in% c(pop1, pop2)] %>%
          fstat_subfun(., type='genos', fstat=fstat) %>%
          lapply(., function(D){
            data.table(POP1=pop1, POP2=pop2, D)
          })

        # Add genome-wide values into main result
        result$genome[[i]] <- pair.calc$genome
        result$locus[[i]] <- pair.calc$locus

        # If permuting
        if(permute==TRUE){
          # Perform the permutation
          pair.perm <- fstat_permute(
            dat[POP %in% c(pop1, pop2)],
            fstat=fstat, numPerms=numPerms, numCores=numCores
          ) %>%
            data.table(POP1=pop1, POP2=pop2, .)

          # Get the p-value
          pair.pval <- lapply(fstat, function(fx){
            p <- sum(pair.perm[[fx]] > pair.calc$genome[[fx]]) / numPerms
            data.table(STAT=fx, PVAL=p)
          }) %>%
            do.call('rbind', .) %>%
            data.table(POP1=pop1, POP2=pop2, .)

          # Add permutations into main results
          result$permute$fstat[[i]] <- pair.perm
          result$permute$pval[[i]] <- pair.pval
        } else{
          result$permute <- NULL
        }

        rm(i)
        # END: Iterate over i pairs
      }

      # Combine results across indexed pairs into single data table
      result$genome <- do.call('rbind', result$genome)

      result$locus <- do.call('rbind', result$locus)

      if(!is.null(result$permute$fstat)){
        result$permute$fstat <- do.call('rbind', result$permute$fstat)
      }

      if(!is.null(result$permute$pval)){
        result$permute$pval <- do.call('rbind', result$permute$pval) %>%
          setorder(., 'STAT')
      }

      # END: FST between population pairs
    }
  }

  # --------------------------------------------+
  # Code for allele frequencies
  # --------------------------------------------+
  if(type=='freqs'){
    # Rename columns
    col.index <- match(c(popCol,locusCol,freqCol,indsCol), colnames(dat))
    colnames(dat)[col.index] <- c('POP','LOCUS','FREQ','INDS')

    # START: FST among all populations
    if(global==TRUE){
      cat('F-statistic calculation on frequency data, global estimate', '\n')
      result <- fstat_subfun(dat, type='freqs')

      # No pemutations
      result$permute <- NULL

      # END: FST among all populations
    }

    # START: FST between pairwise populations
    if(pairwise==TRUE){
      cat('F-statistic calculation on frequency data, global estimate', '\n')

      pair.comp <- combn(unique(dat$POP), 2) %>%
        t() %>%
        as.data.table() %>%
        setnames(., new=c('POP1','POP2'))

      num.pairs <- nrow(pair.comp)

      result <- list(genome=list(), locus=list())

      # START: Iterate over i population pairs
      for(i in 1:nrow(pair.comp)){
        pop1 <- pair.comp$POP1[i]
        pop2 <- pair.comp$POP2[i]

        pair.calc <- dat[POP %in% c(pop1, pop2)] %>%
          fstat_subfun(., type='freqs') %>%
          lapply(., function(D){
            data.table(POP1=pop1, POP2=pop2, D)
          })

        result$genome[[i]] <- pair.calc$genome
        result$locus[[i]] <- pair.calc$locus
        # END: Iteration over i population pairs
      }

      # Combine results
      result$genome <- do.call('rbind', result$genome)
      result$locus <- do.call('rbind', result$locus)
      result$permute <- NULL
      # END: FST between pairwise populations
    }
  }

  # --------------------------------------------+
  # Output results
  # --------------------------------------------+
  return(result)
}

#' @title  Calculate variance components for F-statistics
#'
#' @description Takes a list of values of allele frequencies, sample sizes, number of
#' populations, and heterozygosity, and returns variance components, as
#' per Weir & Cockerham (1984). Assumes that all values are for a single
#' biallelic SNP locus. \cr\cr Note, this function is not exported and is an
#' internal function for \code{fstat_calc}.
#'
#' @keywords internal
#'
#' @param pi Numeric: A vector of allele frequencies for each population.
#'
#' @param ni Numeric: A vector of sample sizes (number of diploid individuals)
#' for each population.
#'
#' @param r Integer: A single value, the number of populations.
#'
#' @param hetStand Logical: Should the estimates be standardised for observed
#' heterozygosity?
#'
#' @param hi Numeric: A vector of observed heterozygosities for each population.
#'
#' @returns Returns a data table. If \code{hetStand==FALSE}, then the list has three
#' columns: \code{$MSP}, mean sqaures for populations; \code{$MSG}, mean squares
#' for gametes; \code{$Nc}, a sample size constant. If \code{hetStand==TRUE},
#' then the three columns are \code{$A}, \code{$B}, \code{$C}, which
#' correspond to the 'a', 'b' and 'c' variance components described in
#' Weir & Cockheram (1984).
#'
#' @references Weir & Cockerham (1984) Evolution. DOI: 10.1111/j.1558-5646.1984.tb05657.x
#' Weir et al. (2002) Annals of Human Genetics. DOI: 10.1146/annurev.genet.36
#'
#' @examples
#' fstat_varcomps(pi=c(0.5,0.25), ni=c(20,20), r=2, hetStand=FALSE)
#'
#' fstat_varcomps(pi=c(0.5,0.25), ni=c(20,20), r=2, hi=c(0.05,0.375), hetStand=TRUE)
fstat_varcomps <- function(pi, ni, r, hi=NULL, hetStand=FALSE){
  require(data.table); require(tidyverse)

  if(hetStand==FALSE){
    # Mean weighted allele frequency
    p.mean <- sum(ni * pi)/sum(ni)

    # Mean squares variance components
    msp <- (1/(r-1)) * sum(ni * (pi - p.mean)^2)
    msg <- (1/sum(ni-1)) * sum(ni * pi * (1-pi))

    # Sample size correction factor
    nc <- (1/(r-1)) * ( sum(ni) - (sum(ni^2)/sum(ni)) )

    # Output as data.table
    return(data.table(MSP=msp, MSG=msg, Nc=nc))
  }

  if(hetStand==TRUE){
    # The mean sample size
    n.mean <- sum(ni/r)

    # The sample size scaling parameter
    nc <- (r*n.mean - sum((ni^2)/(r*n.mean))) / (r-1)

    # The average sample allele frequency
    p.mean <- sum((ni*pi)/(r*n.mean))

    # The variance in allele frequencies
    s2 <- sum( (ni*(pi-p.mean)^2)/((r-1)*n.mean) )

    # The average heterozygosity
    h.mean <- sum( (ni*hi)/(r*n.mean) )

    # The a, b, and c components
    a <- (nc/n.mean) * (s2 - (1/(n.mean-1))*((p.mean*(1-p.mean)) - (s2*(r-1)/r) - (0.25*h.mean)))

    b <- (n.mean/(n.mean-1)) * ((p.mean*(1-p.mean)) - (s2*(r-1)/r) - (h.mean*((2*n.mean-1)/(4*n.mean))))

    c <- 0.5 * h.mean

    # Return as list
    return(data.table(A=a, B=b, C=c))
  }
}

#' @title  Subfunction for calculating F-statistics
#'
#' @description Calculates F-statistics from genotype or allele frequency data in a
#' long-format data table as per Weir & Cockerham (1984). Assumes that all
#' values are for a single biallelic SNP locus. \cr\cr Note, this function is
#' not exported and is an internal function for \code{fstat_calc}.
#'
#' @keywords internal
#'
#' @param dat Data table: In long-format, requires columns \code{$POP},
#' \code{$SAMPLE}, \code{$LOCUS} and \code{$GT} for genotypes;
#' requires \code{$POP}, \code{$LOCUS}, \code{$FREQ}, and \code{$INDS} for
#' allele frequencies
#'
#' @param type Character: either \code{"genos"} or \code{"freqs"}.
#'
#' @param fstat Character: a vector containing the F-statistics to calculate
#' when the data are genotypes. Must include one of \code{"FST"},
#' \code{"FIS"}, or \code{"FIT"}.
#'
#' @returns Returns a list with \code{$genome}, the genome-wide F-statistic, and
#' \code{$locus}, the per locus F-statistic.
#'
#' @examples
#' data_Genos %>%
#'    fstat_subfun(., type='genos', fstat=c('FST','FIS','FIT'))
#'
#' left_join(data_PoolFreqs, data_PoolInfo) %>%
#'    setnames(., 'POOL', 'POP') %>%
#'    fstat_subfun(., type='freqs')
fstat_subfun <- function(dat, type, fstat=NULL){
  require(data.table); require(tidyverse)

  # Empty list
  result <- list(genome=NULL, locus=NULL)

  # If genotype data
  if(type=='genos'){
    # Summarise
    D <- list(
      dat[, .(P=sum(GT)/(length(GT)*2)), by=c('LOCUS','POP')],
      dat[, .(N=length(unique(SAMPLE))), by=c('LOCUS','POP')],
      dat[, .(H=sum(GT==1)/length(GT)), by=c('LOCUS','POP')]
    ) %>%
      reduce(full_join, by=c('LOCUS','POP'))

    num.pops <- D$POP %>% unique %>% length

    # Variance components
    D.varcomp <- D[, fstat_varcomps(pi=P, ni=N, r=num.pops, hi=H, hetStand=TRUE), by=LOCUS]

    # Iterate through requested F-statistics
    for(f in fstat){
      if(f=='FST'){
        result$genome[[f]] <- D.varcomp[, sum(A)/sum(A+B+C)]
        result$locus[[f]] <- D.varcomp[, .(FST=A/(A+B+C)), by=LOCUS]
      } else if(f=='FIS'){
        result$genome[[f]] <- D.varcomp[, 1-(sum(C)/sum(B+C))]
        result$locus[[f]] <- D.varcomp[, .(FIS=1-(C/sum(B+C))), by=LOCUS]
      } else if(f=='FIT'){
        result$genome[[f]] <- D.varcomp[, 1-(sum(C)/sum(A+B+C))]
        result$locus[[f]] <- D.varcomp[, .(FIT=1-(C/sum(A+B+C))), by=LOCUS]
      }
    }

    # Output
    result <- list(
      genome=as.data.table(do.call('cbind', result$genome)),
      locus=result$locus %>% reduce(full_join, by='LOCUS')
    )
  }

  # If allele frequency data
  if(type=='freqs'){
    # Number of populations
    num.pops <- dat$POP %>% unique %>% length

    # Variance components
    D.varcomp <- dat[, fstat_varcomps(pi=FREQ, ni=INDS, r=num.pops, hetStand=FALSE), by=LOCUS]

    # Calculations
    result$genome <- D.varcomp[, .(FST=sum(MSP-MSG)/sum(MSP+((Nc-1)*MSG)))]
    result$locus <- D.varcomp[, .(FST=(MSP-MSG)/(MSP+((Nc-1)*MSG))), by=LOCUS]
  }

  # Output
  return(result)
}

#' @title Permute the genome-wide F-statistics for genotype data
#'
#' @description Uses permutation of individuals among populations to generate a null
#' distribution of F-statistics. Assumes biallelic SNP genotypes. \cr\cr
#' Note, this function is not exported and is an internal function for
#' \code{fstat_calc}.
#'
#' @keywords internal
#'
#' @param dat Data table: In long-format, requires columns \code{$POP},
#' \code{$SAMPLE}, \code{$LOCUS} and \code{$GT} for genotypes;
#' requires \code{$POP}, \code{$LOCUS}, \code{$FREQ}, and \code{$INDS} for
#' allele frequencies
#'
#' @param fstat Character: A vector containing the F-statistics to calculate
#' when the data are genotypes. Must include one of \code{"FST"},
#' \code{"FIS"}, or \code{"FIT"}.
#'
#' @param numPerms Integer: The number of permutations to perform.
#'
#' @param numCores Integer: The number of cores to use. Default = 1. When 1,
#' this will be a single core processes using lapply. If >1, then this will
#' utilise the doParallel package to run on a mini cluster.
#'
#' @returns Returns a data table of permuted F-statistics for each permutation.
#'
#' @examples
#' data(data_Genos)
#'
#' fstat_permute(data_Genos, 'FST', 100, 1)
#'
#' fstat_permute(data_Genos, c('FST','FIT','FIS'), 100, 3)
fstat_permute <- function(dat, fstat, numPerms, numCores=1){
  require('tidyverse'); require('doParallel'); require('data.table')

  pop.tab <- dat[, c('POP','SAMPLE')] %>% unique
  num.samps <- pop.tab %>% nrow

  # If single core
  if(numCores==1){
    result <- lapply(1:numPerms, function(i){
      pop.tab.i <- pop.tab %>% copy
      pop.tab.i$SAMPLE <- sample(pop.tab.i$SAMPLE, size=num.samps, replace=FALSE)

      dat[, c('SAMPLE','LOCUS','GT')] %>%
        left_join(., pop.tab.i, by='SAMPLE') %>%
        fstat_subfun(., type='genos', fstat=fstat) %>%
        .[['genome']] %>%
        data.table(PERM=i, .)
    }) %>%
      do.call('rbind',.)
  }

  # If multi core
  if(numCores>1){
    my.cluster <- makeCluster(numCores)
    registerDoParallel(my.cluster)

    result <- foreach(i=1:numPerms) %dopar% {
      require('genomalicious')

      pop.tab.i <- pop.tab %>% copy
      pop.tab.i$SAMPLE <- sample(pop.tab.i$SAMPLE, size=num.samps, replace=FALSE)

      dat[, c('SAMPLE','LOCUS','GT')] %>%
        left_join(., pop.tab.i, by='SAMPLE') %>%
        fstat_subfun(., type='genos', fstat=fstat) %>%
        .[['genome']] %>%
        data.table(PERM=i, .)
    } %>%
      do.call('rbind',.)

    stopCluster(my.cluster)
  }

  return(result)
}


