#' Manhattan distance between two vectors
#'
#' Computes the manhattan distance between two numeric vectors.
#'
#' @usage
#'
#' @param x Numeric vector.
#' @param y Numeric vector.
#' @return Numeric, a manhattan distance
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#' x=c(1,3,4,5)
#' y=c(1,3,4,5)
#' ManhattanDistance(x, y)
#' y=c(23, 45, 67, 24)
#' ManhattanDistance(x, y)
#' @export
ManhattanDistance=function(x, y){
  sum(abs(x - y))
}
