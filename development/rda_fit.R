library(genomalicious)

data("data_4pops")

Xdat <- data_4pops %>%
  copy %>%
  .[, c('SAMPLE','POP')] %>%
  unique()

Ydat <- data_4pops

predFormula <- 'POP'

sampCol <- 'SAMPLE'

locusCol <- 'LOCUS'

genoCol <- 'GT'

type='genos'

# --------------------------------------------+
# Libraries and assertions
# --------------------------------------------+
for(lib in c('data.table', 'dplyr', 'vegan')){ require(lib, character.only = TRUE)}

# Check that scaling is specified
if(!scaling %in% c('covar', 'corr', 'patterson', 'none')){
  stop('Argument `scaling`` is invalid. See: ?pca_genos')
}

# Get the class of the genotypes
gtClass <- class(Ydat[[genoCol]])

# Check that genotypes are characters or counts
if(!gtClass %in% c('character', 'numeric', 'integer')){
  stop("Check that genotypes are coded as '/' separated characters or as
         counts of the Alt allele. See: ?pca_genos")
}

# Convert characters of separated alleles to counts
if(gtClass=='character'){
  Ydat[[genoCol]] <- genoscore_converter(Ydat[[genoCol]])
}

# Convert numeric allele counts to integers
if(gtClass=='numeric'){
  dat[[genoCol]] <- as.integer(dat[[genoCol]])
}

# --------------------------------------------+
# Code
# --------------------------------------------+
# Convert to a genotype matrix
if(type=='genos'){
  Ymat <- DT2Mat_genos(Ydat, sampCol=sampCol, genoCol=genoCol, locusCol=locusCol)
}

samp.vec <- rownames(Ymat)

Xdat <- Xdat[match(samp.vec, Xdat[[sampCol]]),]

rdaFormula <- paste0('rda(Ymat ~ ', predFormula, ', data=Xdat, scale=FALSE)')

RDAobj <- eval(parse(text=rdaFormula))

