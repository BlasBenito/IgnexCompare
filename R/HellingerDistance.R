#' Hellinger distance between two vectors
#'
#' Computes the Hellinger distance between two numeric vectors.
#'
#' @usage
#'
#' @param x Numeric vector.
#' @param y Numeric vector.
#' @return Numeric, a Hellinger distance
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#' x=c(1,3,4,5)
#' y=c(1,3,4,5)
#' HellingerDistance(x, y)
#' y=c(23, 45, 67, 24)
#' HellingerDistance(x, y)
#' @export
HellingerDistance=function(x, y){
  sqrt(1/2 * sum(sqrt(x)-sqrt(y))^2)
}
