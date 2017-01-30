#' running description
#'
#' description
#'
#' @usage
#'
#' @param param Number of rows of the results table.
#' @return whatever
#' \itemize{
#'  \item{"parameter 1"}{Stuff}
#'  \item{"parameter 2"}{Stuff}
#' }
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#'
#' #generating random input data
#' results.table=GenerateResultsTable(10)
#' str(results.table)
#' @export
PlotDistanceMatrix=function(distance.matrix){

  require(fields)

  x.axis=as.numeric(rownames(distance.matrix))
  y.axis=as.numeric(colnames(manhattan.matrix))

  image.plot(x.axis, y.axis, distance.matrix, xlab="Sequence A", ylab="Sequence B", main="Manhattan dissimilarity.")

}
