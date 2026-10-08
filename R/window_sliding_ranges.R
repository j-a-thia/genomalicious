#' Generate genomic sliding window ranges
#'
#' This function takes in a data.table of chromosome/contigs details and
#' returns overlapping window coordinates.
#'
#' @param dat Data.table: The chromosome/contig details. Requires the columns
#' \code{$CHROM} the sequence ID, the \code{WIDTH} the sequence length in bases.
#'
#' @param step_size Numeric: The step size for the sliding windows.
#'
#' @param flank_size Numeric: The flank size for the sliding windows.
#'
#' @details The function aims to generate overlapping sliding windows.
#' Argument \code{step_size} will dictate the spacing of each window's midpoint.
#' The width of the window is determined by argument \code{flank_size}, which
#' extends either side of the midpoint. An error will return if the
#' \code{flank_size} does not produce overlapping windows.
#'
#' @export
#'
#' @examples
#' library(genomalicious)
#'
#' chrom.tab <- data.table(CHROM=c('chrom1','chrom2'), WIDTH=c(100000,550000))
#'
#' slide.ranges <- window_sliding_ranges(chrom.tab, step_size=500, flank_size=400)
#'
#' print(slide.ranges)

window_sliding_ranges <- function(dat, step_size, flank_size) {
  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  library(data.table)

  stopifnot(
    is.data.table(dat),
    all(c("CHROM", "WIDTH") %in% names(dat)),
    length(step_size) == 1L, step_size > 0,
    length(flank_size) == 1L, flank_size > 0,
    all(dat$WIDTH > 0)
  )

  if (step_size >= 2 * flank_size) {
    stop("'step_size' must be less than 2 * 'flank_size' for windows to overlap.")
  }

  chromosomes <- unique(dat[, .(CHROM, WIDTH)])

  if (anyDuplicated(chromosomes$CHROM)) {
    stop("Each CHROM must have one unique WIDTH.")
  }

  # --------------------------------------------+
  # Generate the windows
  # --------------------------------------------+
  result <- rbindlist(lapply(seq_len(nrow(chromosomes)), function(i) {
    chrom <- chromosomes$CHROM[i]
    width <- chromosomes$WIDTH[i]

    midpoint <- if (width >= step_size) {
      seq(step_size, width, by = step_size)
    } else {
      width
    }

    data.table(
      CHROM = chrom,
      MIDPOINT = midpoint,
      START = pmax(1, midpoint - flank_size + 1),
      END = pmin(width, midpoint + flank_size)
    )
  }))

  result[, WINDOW.ID := .I]
  setcolorder(result, c("CHROM", "WINDOW.ID", "MIDPOINT", "START", "END"))

  return(result)
}
