#' Computes values of a random walk.
#'
#' description
#'
#' @usage
#'
#' @param distance.matrix Numeric distance matrix created by \emph{DistanceMatrix}.
#' @param random.walk Data frame with a random walk created by \emph{GenerateRandomWalk}.
#' @return Data frame, the input random walk with a new column named \emph{distance}
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#' @export
ComputeDistances=function(distance.matrix, slotting){
  #get matrix values
  slotting$distances=distance.matrix[cbind(slotting$a, slotting$b)]
  return(slotting)
}
