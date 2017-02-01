#' Decreases matrix resolution.
#'
#' Decreases the resolution of a distance matrix. It outputs the 10th percentile or the minimum of the values of the original cells inside of every new coarser cell.
#'
#' @usage
#'
#' @param distance.matrix A distance matrix produced by \emph{DistanceMatrix}.
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#' @export
AggregateDistanceMatrix=function(distance.matrix){

  #checking for even dimensions and adding dummy rows and columns
  distance.matrix.cols=ncol(distance.matrix)
  distance.matrix.rows=nrow(distance.matrix)

  if (is.even(distance.matrix.cols) == FALSE){
    distance.matrix.even=cbind(distance.matrix, distance.matrix[,distance.matrix.cols])
  }

  if (is.even(distance.matrix.rows) == FALSE){
    distance.matrix.even=rbind(distance.matrix.even, distance.matrix.even[distance.matrix.rows,])
  }

  #dimensions
  distance.matrix.even.cols=ncol(distance.matrix.even)
  distance.matrix.even.rows=nrow(distance.matrix.even)

  #setting col and row names
  colnames(distance.matrix.even)=1:ncol(distance.matrix.even)
  rownames(distance.matrix.even)=1:nrow(distance.matrix.even)

  #finding least divisor
  divisors.columns=vector()
  divisors.rows=vector()

  for (i in 2:(min(c(distance.matrix.even.cols, distance.matrix.even.rows))/2)){

    #computing divisor
    divisor.columns=distance.matrix.even.cols/i
    divisor.rows=distance.matrix.even.rows/i

    #only storing integer divisors
    if (all.equal(divisor.columns, as.integer(divisor.columns), tolerance=0) == TRUE){divisors.columns=c(divisors.columns, divisor.columns)}

    if (all.equal(divisor.rows, as.integer(divisor.rows), tolerance=0) == TRUE){divisors.rows=c(divisors.rows, divisor.rows)}

  }

  #getting the divisor that is closest to the 10th of the row or col number
  window.cols=divisors.columns[which.min(abs(divisors.columns-10))]
  window.rows=divisors.rows[which.min(abs(divisors.rows-10))]

  #new empty matrix
  coarse.matrix.cols=distance.matrix.even.cols/window.cols
  coarse.matrix.rows=distance.matrix.even.rows/window.rows
  coarse.matrix=matrix(data=NA, nrow=coarse.matrix.rows, ncol=coarse.matrix.cols)
  coarse.matrix.cells=coarse.matrix.cols*coarse.matrix.rows

  #coordinates of the given windows
  min.col=seq(from=1, to=distance.matrix.even.cols, by=window.cols)
  max.col=seq(from=window.cols, to=distance.matrix.even.cols, by=window.cols)
  min.row=seq(from=1, to=distance.matrix.even.rows, by=window.rows)
  max.row=seq(from=window.rows, to=distance.matrix.even.rows, by=window.rows)

  #for column in the coarse matrix
  for (current.column in 1:coarse.matrix.cols){

    #for row in the coarse matrix
    for (current.row in 1:coarse.matrix.rows){

      #minimum value
      # coarse.matrix[current.row, current.column]=min(as.vector(distance.matrix.even[min.row[current.row]:max.row[current.row], min.col[current.column]:max.col[current.column]]))

      #10 percent quantile
      coarse.matrix[current.row, current.column]=quantile(as.vector(distance.matrix.even[min.row[current.row]:max.row[current.row], min.col[current.column]:max.col[current.column]]), probs=0.1)

    }

  }

  #col and row names
  colnames(coarse.matrix)=1:ncol(coarse.matrix)
  rownames(coarse.matrix)=1:nrow(coarse.matrix)

  return(coarse.matrix)

}
