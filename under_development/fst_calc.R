FUN_varcomps_freqs <- function(pi, ni, r){
  # Mean weighted allele frequency
  p.mean <- sum(ni * pi)/sum(ni)

  # Mean squares variance components
  msp <- (1/(r-1)) * sum(ni * (pi - p.mean)^2)
  msg <- (1/sum(ni-1)) * sum(ni * pi * (1-pi))

  # Sample size correction factor
  nc <- (1/(r-1)) * ( sum(ni) - (sum(ni^2)/sum(ni)) )

  # Output as data.table
  data.table(NUMER=msp-msg, DENOM=msp+(nc-1)*msg)
}

FUN_varcomps_genos <- function(ni, pi, hi, r){
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
  return(data.table(NUMER=a, DENOM=a+b+c))
}

FUN_fst_genos <- function(D){
  left_join(
  D[, .(FREQ=sum(GT)/(length(GT)*2)), by=c('POP','LOCUS')],
  D[, .(HO=sum(GT==1)/length(GT)), by=c('POP','LOCUS')]
) %>%
  left_join(
    ., D[, .(N=length(unique(SAMPLE))), by=c('POP','LOCUS')]
  ) %>%
  .[, FUN_varcomps_genos(ni=N, pi=FREQ, hi=HO, r=length(unique(POP))), by=LOCUS] %>%
  .[, FST:=NUMER/DENOM]
}

# For genotypes
dat <- data_4pops
fstLocus <- FUN_fst_genos(dat)

# For genotype permutations
popPermTab <- dat[, c('POP','SAMPLE')] %>% unique

# For frequencies
dat <- data_PoolFreqs %>% copy %>% setnames(., old='POOL', new='POP')
fstLocus <- dat %>%
  .[, FUN_varcomps_freqs(pi=FREQ, ni=INDS, r=length(unique(POP))), by=LOCUS]

# Global FST
fstGlobal <- fstLocus[, sum(NUMER)/sum(DENOM)]



