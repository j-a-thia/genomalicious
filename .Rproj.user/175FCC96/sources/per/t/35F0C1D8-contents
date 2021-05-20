#' Generate input file for \code{MrBayes}
#'
#' Takes a list of data matrices and converts into a NEXUS file that can
#' be read by \code{MrBayes}. Can be used for multiple data partitions and
#' data types.
#'
#' @param dat.list List: Each index contains a matrix for a data partition.
#' Could be a matrix of morphological characters, DNA alignments for a gene, etc.
#' Row names are taxa, column names are the characters, with cells containing
#' the character states per taxa.
#'
#' @param type.vec Character: A vector of the data partition types (the matrices
#' in \code{dat.list}), as per \code{MrBayes} "datatype" requirements:
#' \code{"Standard"}, \code{"DNA"}, \code{"RNA"}, \code{"Protein"}, or
#' \code{"Restriction"}. Note that the length of \code{type.vec} must be the same
#' as the length (number of data partitions) in \code{dat.list}.
#'
#' @param set.vec Character: A vector of the data partition sets (the matrices
#' in \code{dat.list}), effectively a set of user defined names/IDs for the
#' data partitions. Used for generating the "charset" specifications in
#' \code{MrBayes}.
#'
#' @param out.file Character: A single value, the name of the NEXUS file to save.
#'
#' @details Primary functionality is to convert data generated in R into NEXUS.
#' Calculates data ranges and sets up partitions. Addition of \codes{MrBayes}
#' run commands are not supported.
#'
#' @return Writes a NEXUS file.
#'
#' @export
#'
mrbayes_input <- function(dat.list, type.vec, set.vec, out.file){
  require(tidyverse); require(data.table)

  # Check that there are the same number of rows in all data matrices
  if(lapply(dat.list, nrow) %>% unlist() %>% unique() %>% length() > 1){
    stop('The matrices in argument dat.list do not have the same number of rows.
       See ?mrbayes_input.')
  }

  # Check the rownames are congruent across datasets
  dat.rownames <- lapply(dat.list, rownames)

  if(length(dat.rownames)>1){
    dat.rowcheck <- apply(combn(1:length(dat.list),2), 2, function(xy){
      x <- dat.rownames[[xy[1]]]
      y <- dat.rownames[[xy[2]]]

      return(FALSE %in% (x==y))
    })
  } else{
    dat.rowcheck <- FALSE
  }

  if(TRUE %in% dat.rowcheck){
    stop('At least one of the matrices in argument dat.list has rownames that
       are not congruent with other matrices in dat.list. See ?mrbayes_input.')
  }

  # Make sure the data type is specified properly
  if(FALSE %in% (type.vec %in% c('Standard','DNA','RNA','Protein','Restriction'))){
    stop("Argument type.vec must contain only 'Standard', 'DNA', 'RNA', 'Protein'
       or 'Restriction'. See ?mrbayes_input.")
  }

  # Make sure the data type and character set vectors are specified properly.
  if(length(type.vec)!=length(dat.list)){
    stop('Argument type.vec must be same length as dat.list. See ?mrbayes_input.')
  }

  if(length(set.vec)!=length(dat.list)){
    stop('Argument set.vec must be same length as dat.list. See ?mrbayes_input.')
  }

  # Get number of taxa
  num.taxa <- length(dat.rownames[[1]])

  # Get the partition sizes
  part.size <- lapply(dat.list, ncol) %>% unlist()

  # Get the total number of characters
  num.chars <- part.size %>% sum()

  # Unique data types
  type.uniq <- unique(type.vec)

  # If a mixed data type, get the breaks of each data partition and position ranges
  if(length(type.uniq)>1){
    # Get the partition breakpoints
    part.n <- lapply(1:length(dat.list), function(i){
      data.table(PART=i, TYPE=type.vec[i], SET=set.vec[[i]], CHAR=1:ncol(dat.list[[i]]))
    }) %>%
      do.call('rbind', .) %>%
      .[, POS:=1:.N] %>%
      .[, RANGE:='0'] %>%
      .[CHAR==1] %>%
      setorder(., PART)

    # Get partition ranges as character vector
    type.ranges <- lapply(1:nrow(part.n), function(i){
      if(i!=nrow(part.n)){
        x <- paste0(part.n$TYPE[i],':',part.n$POS[i],'-',part.n$POS[i+1]-1)
      } else{
        x <- paste0(part.n$TYPE[i],':',part.n$POS[i],'-',num.chars)
      }
    }) %>%
      unlist()

    # Add ranges into partition table
    part.n$RANGE <- gsub("^[^:]*:", '', type.ranges)

    # The character vector for type to write out
    type.out <- type.ranges %>%
      paste(., collapse=',') %>%
      paste0('mixed(', ., ')')
  } else{
    type.out <- type.uniq
  }

  # Initiate vector to hold text for output
  outText <- c('#NEXUS', '\n', 'Begin data;')

  # Add dimensions
  outText <- c(
    outText,
    paste(
      '\tDimensions',
      paste0('ntax=', num.taxa),
      paste0('nchar=', num.chars),
      ';',
      sep=' ')
  )

  # Add format
  outText <- c(
    outText,
    paste(
      '\tFormat',
      paste0('datatype=', type.out),
      'interleave=yes gap=- missing=?;',
      sep=' ')
  )

  # Add in data matrices
  datmat.out <- lapply(dat.list, function(partX){
    lapply(rownames(partX), function(spe){
      paste0(spe, '\t', partX[spe,] %>% paste(., collapse=''))
    }) %>%
      unlist() %>%
      c(., '')
  }) %>%
    unlist()

  outText <- c(outText, '\tMatrix', '', datmat.out, '\t;', 'End;', '\n')

  # Add Mr Bayes partitions
  charset.out <- lapply(1:nrow(part.n), function(i){
    paste0('\tcharset ', part.n$SET[i], ' = ', part.n$RANGE[i], ';')
  }) %>%
    unlist()

  outText <- c(
    outText,
    'begin mrbayes;', '',
    charset.out,
    paste0(
      '\tpartition favored = ',
      nrow(part.n),
      ': ',
      paste(set.vec, collapse=', '),';'),
    '',
    'end;'
  )

  # Write out
  writeLines(outText, out.file)
}
