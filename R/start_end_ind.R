#' start and end points for each cluster
#'
#'
#' @param n: cluster size
#' @return A list of start and end points for each cluster
#' @export 
start_end_ind <- function(n){
  end <- cumsum(n)
  start <- end - n + 1
  return(list(start = start, end = end))
}