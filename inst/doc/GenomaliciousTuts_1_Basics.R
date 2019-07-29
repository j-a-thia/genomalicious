## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE, warning=FALSE----------------------------------------
library(genomalicious)

## ----eval=TRUE-----------------------------------------------------------
data("genomalicious_Freqs")

# Take a look at the structure of this data, the dim() function reports the dimensons.
dim(genomalicious_Freqs)

# Print the data to screen.
genomalicious_Freqs

## ------------------------------------------------------------------------
# Create a link to raw external datasets in genomalicious
genomaliciousExtData <- paste0(find.package('genomalicious'), '/extdata')

# This command here shows you the VCF file that comes with genomalicious
list.files(genomaliciousExtData, pattern='_poolseq.vcf')

# Use this to create a path to that file
vcfPath <- paste0(genomaliciousExtData, '/genomalicious_poolseq.vcf')

# The value of vcfPath with depend on your system
vcfPath

## ------------------------------------------------------------------------
# Read in and print the first 10 lines of the demo VCF
head(readLines(vcfPath), 10)

## ------------------------------------------------------------------------
# Import VCF into R using the path name
poolSnps <- vcf2DT(vcfPath)

# First 8 rows of the imported SNP data
head(poolSnps, 8)

# The imported data is stored as a data.table object
class(poolSnps)

## ------------------------------------------------------------------------
# Column vector for sample and loci
head(poolSnps$SAMPLE)

head(poolSnps$LOCUS)

## ---- warning=FALSE------------------------------------------------------
# Subset rows by depth.
poolSnps[DP > 110,]

# You can apply functions to columns, for example,
# take the mean depth.
poolSnps[, mean(DP)]

# You can even specify subsets of the data that you want to
# apply the function to, for example, take the mean depth
# for each locus.
poolSnps[, mean(DP), by=LOCUS]

# You can also use this feature of column manipulation to
# apply a function to the data and create a new column using
# the ':=' notation. For example, add a new column
# containting the frequency of reads that contain the
# reference allele.
poolSnps[, R.FREQ:=RO/DP]
head(poolSnps, 4)

## ------------------------------------------------------------------------
# Converting to a matrix
matDat <- as.matrix(poolSnps)
class(matDat)
head(matDat, 4)

# Converting to a data frame
dfDat <- as.data.frame(poolSnps)
class(dfDat)
head(dfDat, 4)

## ------------------------------------------------------------------------
# Convert long format data table of allele frequencies to a matrix
freqMat <- DT2Mat_freqs(poolSnps, popCol='SAMPLE', locusCol='LOCUS', freqCol='R.FREQ')
freqMat

# Check the class
class(freqMat)

# Sample names are stored in the rows of the matrix
rownames(freqMat)

# Loci names are stored in the columns of the matrix
colnames(freqMat)

## ------------------------------------------------------------------------
freqDT <- DT2Mat_freqs(freqMat, popCol='SAMPLE', locusCol='LOCUS', freqCol='R.FREQ', flip=TRUE)
head(freqDT, 8)

## ------------------------------------------------------------------------
# Import the dataset
data(genomalicious_4pops)
indSnps <- genomalicious_4pops

# Have a look at its structure
head(indSnps)

# The number of unique samples, per population
indSnps[, length(unique(SAMPLE)), by=POP]

# The number of unique loci
length(unique(indSnps$LOCUS))

## ------------------------------------------------------------------------
# Practise run with simple vectors
genoscore_converter(c('0/0', '0/1', '1/1'))

genoscore_converter(c(0L, 1, 2))

# Now manipulate the data table of genotypes
indSnps[, GT:=genoscore_converter(GT)]

head(indSnps, 8)

## ------------------------------------------------------------------------
# Make a matrix of genotypes in separated format
genoMat.sep <- DT2Mat_genos(indSnps, sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT', genoScore='sep')

genoMat.sep[1:8, 1:4]

# Make a matrix of genotypes in counts format
genoMat.counts <- DT2Mat_genos(indSnps, sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT', genoScore='counts')

genoMat.counts[1:8, 1:4]

# Convert counts genotype matrix into long format data table
# of separated alleles
genoDT.sep <- DT2Mat_genos(genoMat.counts, sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT', genoScore='sep', flip=TRUE)

head(genoDT.sep, 8)

# Convert separated genotype matrix in long format data table
# of reference allele counts
genoDT.counts <- DT2Mat_genos(genoMat.counts, sampCol='SAMPLE', locusCol='LOCUS', genoCol='GT', genoScore='counts', flip=TRUE)

head(genoDT.counts, 8)

