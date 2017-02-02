#' Search minimum cost path using nearest neighbour.
#'
#' description
#'
#' @usage
#'
#' @param distance.matrix A distance matrix provided by \emph{DistanceMatrix}.
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#'
#' #generating random input data
#' results.table=GenerateResultsTable(10)
#' str(results.table)
#' @export
LeastCostNN=function(distance.matrix){

  current.col=1
  current.row=1
  last.col=ncol(distance.matrix)
  last.row=nrow(distance.matrix)
  path.rows=vector()
  path.cols=vector()
  path.rows[1]=current.row
  path.cols[1]=current.col

  #finding least cost path in the coarse matrix
  while (current.col!=last.col & current.row!=last.row){

    # print(paste("Current row = ", current.row, "; current column = ", current.col, sep=""))

    #MOVING WITHIN MATRIX BOUNDS
    if (current.row+1 <= last.row & current.col+1 <= last.col)
    {

      if (distance.matrix[current.row+1, current.col] < distance.matrix[current.row, current.col+1])
      {
        current.row=current.row + 1
      } else {
        current.col=current.col + 1
      }

      #writing path
      path.rows=c(path.rows, current.row)
      path.cols=c(path.cols, current.col)

    }

    #LAST STEP
    if (current.row==last.row & current.col < last.col){
      # print("Last row reached.")
      repeat{
        current.col=current.col+1
        path.rows=c(path.rows, current.row)
        path.cols=c(path.cols, current.col)
        if (current.col==last.col){break}
      }

    }

    if (current.row < last.row & current.col == last.col){
      # print("Last column reached.")
      repeat{
        current.row=current.row+1
        path.rows=c(path.rows, current.row)
        path.cols=c(path.cols, current.col)
        if (current.row==last.row){break}
      }

    }

  }

  #creating dataframe
  solution=data.frame(a=path.rows, b=path.cols)

  #computing cost
  solution=ComputeDistances(distance.matrix, solution)

  return(solution)

}
