#' Test the circularity of a genetic sequence
#'
#' A rough function to assess the circularity of a genetic sequence, for example,
#' an assembly of a bacterial chromosome or a eukaryotic organelle genome.
#'
#' @param query_seq Character: The query sequence.
#'
#' @param word_size Integer: The size of the words to search for. Default = 20.
#'
#' @param search_start Integer: The starting base position. Default = 1.
#'
#' @param step_size Integer: The sliding window step size. Default = 1.
#'
#' @param wiggle Integer: The "wiggle room" for the word size in case the
#' sliding window over-shoots the end of the sequence. Default = 2.
#'
#' @details A sliding window is used to assess for the presence of a replicated
#' character string in the query sequence. This replicated pattern is used to
#' infer the circular start and end points.
#'
#' The function \code{circle_cutter} can then be used to excise the single
#' non-replicated linear sequence.
#'
#' @return If the function finds a replicated character string, it will return
#' a data.table with the columns:
#' \enumerate{
#'    \item \code{$STEP}: Integer, the step number in the sliding window.
#'    \item \code{$WORD}: Character, the character string.
#'    \item \code{$START}: Integer, the start position.
#'    \item \code{$END}: Integer, the end position.
#'    \item \code{$SIZE}: Integer, the word size.
#' }
#'
#' @export
#'
#' @examples
#' x <- 'AATTGGCCACTATCTGCTAGCTAGCATAGCATCGATCAGCATGACGCGCAAAATTGGCC'
#'
#' # Find character motif that is repeated
#' motif_hits <- circularity_test(x, word_size = 8)
#' motif_seq <- substr(x, motif_hits[1,1], motif_hits[1,2])
#'
#' circle_cutter(query_seq = x, motif=motif_seq)
#'
circularity_test <- function(
    query_seq, word_size=20,
    search_start=1, step_size=1, wiggle=2){

  # Make sure the query is a character
  query_seq <- as.character(query_seq)

  # Query length
  query_len <- nchar(query_seq)

  # Start and end of the first word
  st <- search_start
  en <- st + word_size - 1

  # Track the steps
  step_track <- 0
  num_steps <- 1

  # START WHILE LOOP
  # Continue stepping through words
  while(step_track==0){
    hits <- NULL

    if(st < query_len & en <= query_len){
      # If both the start and end of the word are less then query length
      word.i <- substr(query_seq, st, en)

      hits <- str_locate_all(query_seq, word.i)[[1]]
    } else if(st > query_len){
      # Kill if start is greater than query length
      step_track <- 1
    } else if(en > query_len){
      # If the end if greater than query length, remove wiggle room and see
      # if it helps.
      en_wig <- en - wiggle
      if(st < query_len & en_wig < query_len){
        word.i <- substr(query_seq, st, en_wig)

        hits <- str_locate_all(query_seq, word.i)[[1]]
      } else{
        warning(
          paste0(
            '\n',
            'Got to end of query with no circularity detected. \n\n',
            'You might also consider changing the word size. \n\n',
            'The last word was too truncated to test: \n\n',
            'Start = ', st, ', End = ', en, ', Size = ', nchar(word.i), '\n\n',
            word.i
          )
        )
        hits <- matrix(NA, nrow=1, ncol=2)
        step_track <- 1
      }
    }

    # Were there multiple hits?
    if(nrow(hits)==1){
      st <- st + step_size
      en <- en + step_size
      num_steps <- num_steps + 1
    } else if(is.null(hits)){
      warning(
        paste0(
          '\n',
          'Got to end of query with no circularity detected. \n\n',
          'You might also consider changing the word size.'
        )
      )
      step_track <- 1
    } else{
      step_track <- 1
    }
  # END WHILE LOOP
  }

  return(hits)
}
