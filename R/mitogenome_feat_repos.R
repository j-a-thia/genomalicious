#' Reposition a table of mitogenome features
#'
#' Takes a data table of mitogenome features a reposisitions it relative to a
#' new start position. Could be used to reposition genomic features of any
#' circular chromosomes.
#'
#' @param featuresDT Data.table: Contains the genomic features. Requires only
#' two columns: \code{$START} and \code{$END}, which are the start and end base
#' positions for each a feature. Each row is a unique feature.
#'
#' @param feat_base_zero Logical: Do the base positions in \code{featuresDT$START}
#' and \code{featuresDT$END} have a base value of zero. In other words, does the
#' count start at 0 (the first base)? See Details.
#'
#' @param genome_len Integer: The genome length, assuming the circular chromosome
#' has been linearised. If unspecified, it takes on the value of \code{max(featuresDT$END)}.
#' It must be one the same scale as \code{feat_base_zero}. Default is NULL.
#'
#' @param new_start_pos Integer: The new starting position for the features.
#' Must be a value >= 2 and must be <= \code{genome_len}. It must be on the same scale
#' as \code{feat_base_zero}.
#'
#' @details Consider the sequence: ATATAGCCG. The first base 'A', would have have
#' a numbering value of 0 and 1 for base 0 and 1, respectively. Likewise, the
#' last base, 'G', would have a value of 8 or 9. \cr\cr
#'
#' Repositioning is done with respect to the initial base level. So to achieve
#' the split ATAT|AGCCG, \code{new_start_pos==4} for \code{feat_base_zero==TRUE},
#' and \code{new_start_pos==5} for \code{feat_base_zero==FALSE}.
#'
#' NOTE: In this iteration of the function, it is possible that a feature can
#' be split, such that it is separated into sections on the left and right end
#' of a linear sequence (by virtue of the circularity of the mitogenome).
#' However, once this split has been made, there is no current functionality to
#' reverse the output (i.e. to merge a feature that has been split into sections).
#'
#' @return Returns a data.table with same structure as \code{featuresDT}, but
#' features have new start and end positions as specified by the repositioning
#' value of \code{new_start_pos}. If the \code{feat_base_zero==TRUE}, then the
#' values are returned as base 0, otherwise they are returned as base 1.
#'
#' @examples
#' #' library(genomalicious)
#'
#' # Create a link to raw external datasets in genomalicious
#' genomaliciousExtData <- paste0(find.package('genomalicious'), '/extdata')
#'
#' # Read in a GENBANK file of the Bathygobius cocosensis mitogenome
#' gbk.read <- mitoGbk2DT(paste(genomaliciousExtData, 'data_Bcocosensis.gbk', sep='/'))
#' head(gbk.read)
#'
#' # Reposition to start at COX1
#' gbk.read[NAME=='COX1']
#'
#' reposCOX1 <- mitogenome_feat_repos(
#'   featuresDT=gbk.read,
#'   new_start_pos=5495,
#'   feat_base_zero=FALSE,
#'   genome_len=16692
#' )
#'
#' # Reposition in the middle of COX1
#' reposMidCOX1 <- mitogenome_feat_repos(
#'   featuresDT=gbk.read,
#'   new_start_pos=6272,
#'   feat_base_zero=FALSE,
#'   genome_len=16692
#' )
#'
#' @export

mitogenome_feat_repos <- function(featuresDT, new_start_pos, feat_base_zero, genome_len=NULL){
  require(data.table); require(tidyverse); require(stringr)

  if(is.null(featuresDT)==FALSE){
    if(sum(c('START', 'END') %in% colnames(featuresDT))!=2){
      stop('Argument featuresDT must have the columns $START and $END.
    See ?mitogenome_feat_repos.')
    }

    if(is.logical(feat_base_zero)==FALSE){
      stop('Argument feat_base_zero must be a logical value.
    See ?mitogenome_feat_repos.')
    }
  }

  # Make sure featuresDT is a data.table
  featuresDT <- as.data.table(featuresDT)

  # Check if features are base zero.
  if(feat_base_zero==TRUE){
    featuresDT <- featuresDT %>%
      .[, START:=START+1] %>%
      .[, END:=END+1]

    new_start_pos <- new_start_pos + 1
  }

  # Genome length checks
  if(genome_len < max(featuresDT$END)){
    stop('Argument genome_len cannot be less than max(featuresDT$END).
  See ?mitogenome_feat_repos')
  }

  if(is.null(genome_len)==TRUE){
    genome_len <- max(featuresDT$END)
  }

  # Start position checks
  if(new_start_pos < 2){
    stop('Argument new_start_pos must be an integer > 2. See ?mitogenome_feat_repos.')
  }

  if(new_start_pos > genome_len){
    stop('Argument new_start_pos must be < genome_len. See ?mitogenome_feat_repos.')
  }

  # Size of initial sequences to left (L) and right (R) of the new start position.
  size.init.L <- new_start_pos - 1
  size.init.R <- genome_len - size.init.L

  # Empty list to hold repositioned data
  featuresRepos <- list()

  # Iterate through each row of the original features data table
  for(f in 1:nrow(featuresDT)){
    # Get the original positions
    f.start <- featuresDT$START[f]
    f.end <- featuresDT$END[f]

    # Will the feature need to be split? I.e. is the original feature start <
    # new seq start, and is the original feature end > new seq start?
    f.is.split <- new_start_pos > f.start & new_start_pos < f.end

    # Code for handling different scenarios
    if(f.is.split==TRUE){
      # If the feature needs to be split

      # New left region
      f.start.new.L <- new_start_pos - size.init.L
      f.end.new.L <- f.end - size.init.L

      # New right region
      f.start.new.R <- f.start + size.init.R
      f.end.new.R <- (new_start_pos - 1) + size.init.R

      # Output
      f.out <- data.table(
        featuresDT[f, !c('START', 'END')],
        START=c(f.start.new.L, f.start.new.R),
        END=c(f.end.new.L, f.end.new.R),
        REPOS=c('Split, left region', 'Split, right region')
      )
    } else if(f.start > new_start_pos){
      # If the feature is shifted left, i.e. Original feature start > new seq start.
      f.out <- data.table(
        featuresDT[f, !c('START', 'END')],
        START=f.start - size.init.L,
        END=f.end - size.init.L,
        REPOS='Repositioned left'
      )
    } else if(f.end < new_start_pos){
      # If the feature is shifted right, i.e. Original feature end < new seq start.
      f.out <- data.table(
        featuresDT[f, !c('START', 'END')],
        START=f.start + size.init.R,
        END=f.end + size.init.R,
        REPOS='Repositioned right'
      )
    } else if(f.start == new_start_pos){
      # If the features starts at the new start position, i.e. original feature
      # start == new seq start feature.
      f.out <- data.table(
        featuresDT[f, !c('START', 'END')],
        START=f.start - size.init.L,
        END=f.end - size.init.L,
        REPOS='Repositioned left'
      )
    }

    # Add to list
    featuresRepos[[f]] <- f.out
  }

  # Bind as rows
  featuresRepos <- do.call('rbind', featuresRepos)

  # Adjust base positions to base zero if that was iniial format
  if(feat_base_zero==TRUE){
    featuresRepos <- featuresRepos %>%
      .[, START:=START-1] %>%
      .[, END:=END-1]
  }

  # Output
  return(featuresRepos)
}


