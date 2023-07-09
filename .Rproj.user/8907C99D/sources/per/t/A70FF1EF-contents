#' Correct uneven gaps in an assembled genome
#'
#' Uneven gap sizes can form during genome assembly as different
#' assembling, polishing, and gap closing tools pile on top of each
#' other. This makes submission of a genome assembly problematic,
#' because varying gap sizes are ambiguous, and ideally we want unknown
#' gaps to have a standardised size. \cr\cr
#' This function takes a genome assembly and corrects gap sizes given
#' a user-define threshold. It is assumed that any gap larger than the
#' threshold is an unknown gap.
#'
#' @param genomeSS DNAStringSet: the genome assembly as a 'DNAStringSet'
#' object from the \code{Biostrings} package.
#'
#' @param threshold Integer: the minimum string of 'N's to consider a gap.
#'
#' @param correct Integer: the standardised number of 'N's for gaps of unknown size.
#'
#' @param numCores Integer: the number of cores for parallel processing. Default is 2.
#'
#' @returns Returns the genome assembly back as a DNAStringSet object with all
#' gaps meeting the threshold size standardized to the corrected size.
#' Sequences with corrected gaps are labelled with '_correct_gaps'.
#'
#' @export
assembly_correct_gaps <- function(genomeSS, threshold, correct, numCores=2){
  require(Biostrings); require(data.table); require(doParallel); require(tidyverse)

  # Quick checks
  if(!'DNAStringSet' %in% class(genomeSS)){
    stop('Argument `genomeSS` must be a "DNAStringSet" class. See ?assembly_correct_gaps.')
  }

  if(class(threshold)!='integer' | threshold<1){
    stop('Argument `treshold` must be an "integer" class with a values >0. See ?assembly_correct_gaps.')
  }

  if(class(correct)!='integer' | correct<1){
    stop('Argument `correct` must be an "integer" class with a values >0. See ?assembly_correct_gaps.')
  }

  if(class(numCores)!='integer' | numCores<1){
    stop('Argument `numCores` must be an "integer" class with a values >0. See ?assembly_correct_gaps.')
  }

  # Set up cluster
  my_cluster <- makeCluster(numCores)

  registerDoParallel(my_cluster)

  # Iterate through each ith genomic sequence
  resultsList <- foreach(i=1:length(genomeSS)) %dopar% {
    require(Biostrings); require(data.table); require(tidyverse)
    # Convert the sequence into a long-format data table
    seq_i_tab <- genomeSS[i] %>%
      as.matrix() %>%
      t() %>%
      as.data.table() %>%
      setnames(., new='DNA')

    # Add in an index and test whether the bases are Ns
    seq_i_tab[, INDEX:=1:.N]
    seq_i_tab[, IS.N:=DNA=='N']

    # Subset and sort indexes with Ns. Gap IDs are currently unknown.
    n_i_tab <- seq_i_tab[IS.N==TRUE] %>%
      setorder(., INDEX) %>%
      .[, GAP:=0L]

    # Counter of unique gap IDs.
    uniq_gap_count <- 1L

    # Iterate through each jth index where there is an N.
    # Test whether the jth index is exactly j+1 the previous index.
    # If yes, then they are consecutive and part of the same gap.
    # If no, then they are different gaps, and the counter goes up
    for(j in 1:nrow(n_i_tab)){
      if(j == 1){
        n_i_tab$GAP[j] <- uniq_gap_count
      } else if(j > 1){
        index_j <- n_i_tab$INDEX[j]
        index_j_m1 <- n_i_tab$INDEX[j-1]
        index_consec <- index_j - 1 == (index_j_m1)

        if(index_consec==TRUE){
          n_i_tab$GAP[j] <- uniq_gap_count
        } else if(index_consec==FALSE){
          uniq_gap_count <- uniq_gap_count + 1
          n_i_tab$GAP[j] <- uniq_gap_count
        }
      }
    }

    # Summarise the size of unique gaps
    size_i_tab <- n_i_tab[, .(SIZE=length(INDEX)), by=GAP] %>%
      .[, KEEP:=SIZE>=threshold]

    keep_gaps_i <- size_i_tab[KEEP==TRUE]$GAP

    # Return original sequence, or sequence with standardised gaps
    if(nrow(size_i_tab)==0){
      result_i <- genomeSS[i]
    } else {
      # The N string to replace corrected gaps
      n_string <- paste(rep('N',correct), collapse='')

      # Starting site for each gap
      keep_i_tab <- n_i_tab[GAP %in% keep_gaps_i] %>%
        copy %>%
        .[, START:=INDEX==min(INDEX), by=GAP]

      # Add in the info for the gaps to keep
      seq_i_tab <- left_join(seq_i_tab, keep_i_tab)

      # Replace the string of the start site for each gap
      seq_i_tab[!is.na(GAP) & START==TRUE, DNA:=n_string]

      # Keep only sites that are 'A', 'T', 'G', 'C', or the first
      # site for each gap. Make sure the sites are in order then
      # paste them back together and convert into a string set.
      new_seq_i <- seq_i_tab[is.na(GAP) | (!is.na(GAP) & START==TRUE)] %>%
        setorder(., INDEX) %>%
        .[['DNA']] %>%
        paste(., collapse='') %>%
        DNAStringSet()

      # Add in the sequence name, but note that it is corrected.
      names(new_seq_i) <- paste0(names(genomeSS[i]), '_correct_gaps')

      # Add to results
      result_i <- new_seq_i
    }

    # Cleanup and return results
    gc()
    result_i
  }

  stopCluster(my_cluster)

  gc()

  do.call('c', resultsList)
}
