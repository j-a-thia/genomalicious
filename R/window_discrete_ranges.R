#' Generate discrete genomic window ranges
#'
#' This function takes in a data.table of chromosome/contig details and
#' returns non-overlapping window coordinates.
#'
#' @param dat Data.table: Chromosome/contig details. Requires the columns
#' \code{CHROM}, the sequence ID, and \code{WIDTH}, the sequence length in bases.
#'
#' @param window_size Numeric: The intended total width of each window in bases.
#'
#' @details Windows are generated from the start of each chromosome and do not
#' overlap. If the chromosome width is not divisible by \code{window_size},
#' the final window is shorter.
#'
#' @export
#'
#' @examples
#' library(genomalicious)
#'
#' chrom.tab <- data.table(CHROM = c("chrom1", "chrom2"), WIDTH = c(100000, 550000))
#'
#' discrete.ranges <- window_discrete_ranges(chrom.tab, window_size = 50000)
#'
#' print(discrete.ranges)

window_discrete_ranges <- function(dat, window_size) {
  # --------------------------------------------+
  # Libraries and assertions
  # --------------------------------------------+
  library(data.table)

  stopifnot(
    is.data.table(dat),
    all(c("CHROM", "WIDTH") %in% names(dat)),
    length(window_size) == 1L,
    is.finite(window_size),
    window_size > 0,
    window_size %% 1 == 0,
    all(dat$WIDTH > 0)
  )

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

    start <- seq(1, width, by = window_size)
    end <- pmin(start + window_size - 1, width)

    data.table(
      CHROM = chrom,
      MIDPOINT = floor((start + end) / 2),
      START = start,
      END = end
    )
  }))

  result[, WINDOW.ID := .I]
  setcolorder(result, c("CHROM", "WINDOW.ID", "MIDPOINT", "START", "END"))

  return(result)
}
