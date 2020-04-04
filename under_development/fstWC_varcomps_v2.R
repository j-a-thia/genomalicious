#### This should integrate fine into fstWC_freqs, but fstWC_genos will
#### need to be modified to incorporate.

data("genomalicious_4pops")

dat <- NULL
samp_size <- NULL
pop_id <- NULL
input_type <- NULL


# Setup For genotypes
mat4pops <- DT2Mat_genos(genomalicious_4pops, 'SAMPLE', 'LOCUS', 'GT')

dat <- mat4pops[, ]
dat <- apply(dat, 2, genoscore_converter)

samp_size <- matrix(rep(30, nrow(dat)), nrow=nrow(dat), ncol=ncol(dat))
pop_id <- c(rep('pop1',30), rep('pop2',30), rep('pop3',30), rep('pop4',30))
input_type <- 'genos' # or 'freqs'

# Setup for freqs
dat <- genomalicious_Freqs

samp_size <- matrix(rep(30, nrow(dat)), nrow=nrow(dat), ncol=ncol(dat))
colnames(samp_size) <- colnames(dat)
pop_id <- NULL
input_type <- 'freqs'

# dat Matrix: Allele frequencies for populations, or biallelic genotypes of
# of individuals scores as integer counts of the Alt allele (0, 1, 2).
# Populations or individuals are in rows, loci are in columns.
# Column names are loci IDs and must exactly match those in argument \code{samp_size}.

# input_type Character: One of two possible values: 'genos', calcualte
# variance components from genotype matrix, or 'freqs', calculate variance
# components from an allele frequency matrix.

# samp_size Matrix: The sample size for each locus in each population.
# Allows for different sample sizes at each locus, for example, if there
# is missing data. Populations in rows, loci in columns. Rows must be
# in the same order as rows in \code{dat}. Column names are loci IDs
# and must all occur in \code{dat}.

# pop_id  Charater: A vector of population IDs. Default is NULL and
# is only required if \code{input_type==TRUE}.
# Must be the same length as \code{nrow(dat)} and the order of values
# must match the order of rows in \code{dat}.

# DETAILS:
# It is important to labels rows as population ID OR specify the poopulation
# ID in \code{pop_vec}. If you don't, the function will assume every unique
# row ID is a population ID.

require(data.table)

# Check input type
if(!input_type %in% c('genos', 'freqs')){
  stop("Argument input type is incorrectly specified, see ?fstWC_varcomps")
}

# If there are no locus names, stop.
if(length(colnames(dat))!=ncol(dat)){
  stop("All colnames of argument dat (the loci) need to named, see ?fstWC_varcomps")
}

# If using genotypes, but no population IDs specified...
if(input_type=='genos' & is.null(pop_id)==TRUE){
  stop("Argument input_type=='genos', but pop_id is not specified, see ?fstWC_varcomps")
}

# If using frequencies, but no sample size matrix
if(input_type=='freqs' & is.null(samp_size)==TRUE){
  stop("Argument input_type=='freqs', but samp_size is not specified, see ?fstWC_varcomps")
}

# If using frequencies:
# (1) Check columns are equal for dat and samp_size,
# (2) Make sure column names are the same for dat and samp_size,
# (3) For good measure, make columns be in the same order.
if(input_type=='freqs'){
  if(ncol(dat)!=ncol(samp_size)){
    stop("The number columns (loci) must be equal in dat and samp_size, see ?fstWC_varcomps")
  }
  if(sum(colnames(samp_size) %in% colnames(dat))!=ncol(dat)){
    stop("All column names in dat must be in samp_size, see ?fstWC_varcomps")
  }

  samp_size <- samp_size[,colnames(dat)]
}

# Get the locus names
lociNames <- colnames(dat)

if(input_type=='freqs'){
  # The number of populations
  r <- nrow(dat)

  # Iterate over each locus
  lociVar <- lapply(lociNames, function(locus){
    # Allele frequency and sample size for each ith population
    pi <- dat[, locus]
    ni <- samp_size[, locus]

    # Mean weighted allele frequency
    p.mean <- sum(ni * pi)/sum(ni)

    # Mean squares variance components
    msp <- (1/(r-1)) * sum(ni * (pi - p.mean)^2)
    msg <- (1/sum(ni-1)) * sum(ni * pi * (1-pi))

    # Sample size correction factor
    nc <- (1/(r-1)) * ( sum(ni) - (sum(ni^2)/sum(ni)) )

    # Return the locus specific parameters
    return(data.table(LOCUS=locus, NUMER=msp-msg, DENOM=msp+(nc-1)*msg))
  })
}

if(input_type=='genos'){
  popNames <- sort(unique(pop_id))

  # The number of populations
  r <- length(unique(pop_id))

  # Iterate over each locus
  lociVar <- lapply(lociNames, function(locus){
    # Get the allele frequencies, observed heterozygosity,
    # and sample size for each ith population
    gtDT <- na.omit(data.table(POP=pop_id, GT=dat[, locus]))
    piDT <- gtDT[, .( FREQ=sum(GT)/(length(GT)*2) ), by=POP]
    niDT <- gtDT[, .( SAMPS=length(GT) ), by=POP]
    hiDT <- gtDT[, .( HO=sum(GT==1)/length(GT) ), by=POP]
    pi <- piDT$FREQ
    names(pi) <- piDT$POP
    ni <- niDT$SAMPS
    names(ni) <- niDT$POP
    hi <- hiDT$HO
    names(hi) <- hiDT$POP

    # Reorder, just to be sure
    ni <- ni[names(pi)]
    hi <- hi[names(pi)]

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

    # Return as numerator and denominator
    return(data.table(LOCUS=locus, NUMER=a, DENOM=a+b+c))
  })
}

# Merge list items together
lociVar <- do.call('rbind', lociVar)

# Return data table
return(lociVar)
