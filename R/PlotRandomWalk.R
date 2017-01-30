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
PlotRandomWalk=function(distance.matrix, random.walk){

  require(gplots)

  PlotDistanceMatrix(distance.matrix)
  lines(random.walk$b, random.walk$a, lwd=2, col="black")

}
