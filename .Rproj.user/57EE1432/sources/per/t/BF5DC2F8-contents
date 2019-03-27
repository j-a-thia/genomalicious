FUN_format_dadi_poolreads <- function(dat){
  # Quickly format data into the Dadi format following Gutenkunst methood of diffusion
  # approximation for demographic inference:
  # http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000695
  # This method is appropriate for pooled RADseq data where each population is represented
  # by a single sample and reads represent the allele observations.
  #
  # INPUTS:
  #   dat   (data.table)    The input data. It is expected that the data is in long
  #                           format; i.e. populations are all in a single column in rows.
  #                           The data table must contain the following column names:
  #
  #                             LOCUS, a unique locus ID.
  #                             SAMPLE, a population ID.
  #                             RO, the reference read observations (counts).
  #                             AO, the alternate read observations (counts).
  #                             REF, the reference allele state.
  #                             ALT, the alternate allele state.
  #
  # OUTPUTS:
  #   A data table is returned in the Dadi format:
  #   http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000695
  #
  # BEGIN ............

  library(data.table); library(tidyr)

  # Split the data by LOCUS
  dat.spl <- split(dat, as.factor(dat$LOCUS))

  # Use spread() from tidyr to put SAMPLE in columns and RO/AO in cells
  allele1 <- spread(dat[,c('LOCUS','SAMPLE','RO')], SAMPLE, RO)
  allele2 <- spread(dat[,c('LOCUS','SAMPLE','AO')], SAMPLE, AO)

  # Extract the allelic states and remove duplicates rows
  states <- dat[,c('LOCUS','REF','ALT')]
  states <- states[which(!duplicated(states)),]

  # Add in some '-' into the allelic states, this is a Dadi requirement.
  for(i in 1:nrow(states)){
    states$REF[i] <- paste0('-',states$REF[i],'-')
    states$ALT[i] <- paste0('-',states$ALT[i],'-')
  }

  # Sort the order of individual datatables by LOCUS
  setorder(allele1, LOCUS); setorder(allele2, LOCUS); setorder(states, LOCUS)

  # Combine all individual data tables into single object
  dadi.dt <- data.table(states[,2:3], Allele1=gsub('-', '', states$REF)
                        , allele1[,2:ncol(allele1)], Allele2=gsub('-', '', states$ALT), allele2[,2:ncol(allele2)]
                        , Locus=states$LOCUS)

  # Return Dadi format
  return(dadi.dt)
  # ............ END
}
