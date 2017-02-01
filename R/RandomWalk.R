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
RandomWalk=function(distance.matrix){

  current.col=1
  current.row=1
  last.col=ncol(distance.matrix) #sequence B
  last.row=nrow(distance.matrix) #sequence A
  path.rows=vector()
  path.cols=vector()
  path.rows[1]=current.row
  path.cols[1]=current.col
  moves.left.col=last.col-1
  moves.left.row=last.row-1

  #finding least cost path in the coarse matrix
  repeat{

      #tossing a coin
      random=sample(c(1,0))

      #filter values if bounds are reached
      if (current.row == last.row & current.col < last.col){
        random=c(0,1)
      }

      if (current.col == last.col & current.row < last.row){
        random=c(1,0)
      }

      #moving
      current.row = current.row + random[1]
      current.col = current.col + random[2]

      #writing path
      path.rows = c(path.rows, current.row)
      path.cols = c(path.cols, current.col)

      if (current.row == last.row & current.col == last.col){break}

  }#end of repeat

  #creating dataframe
  solution=data.frame(a=path.rows, b=path.cols)

  return(solution)

}
