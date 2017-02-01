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
ComputeDistances=function(distance.matrix, random.walk){
  #get matrix values
  random.walk$distances=distance.matrix[cbind(random.walk$a, random.walk$b)]
  return(random.walk)
}
