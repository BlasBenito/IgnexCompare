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
LeastCostNNRandom=function(distance.matrix, max.random.threshold=NULL){

  #starting values for the search
  current.col=1
  current.row=1
  last.col=ncol(distance.matrix)
  last.row=nrow(distance.matrix)
  path.rows=vector()
  path.cols=vector()
  path.rows[1]=current.row
  path.cols[1]=current.col
  # random.threshold=runif(1, 0, max.random.threshold)
  random.threshold=max.random.threshold

  #searching
  repeat{

    #random move or targeted move?
    if (runif(1) > random.threshold){

      #TARGETED

      #while we are not done...
      if (current.row < last.row & current.col < last.col){

        #move to neighboring cell with the lower distance value
        if (distance.matrix[current.row+1, current.col] < distance.matrix[current.row, current.col+1])

        #defining move to do
        {move=c(1,0)} else {move=c(0,1)}

      }

      #change move if bounds are reached
      if (current.row == last.row & current.col < last.col){move=c(0,1)}
      if (current.col == last.col & current.row < last.row){move=c(1,0)}

      #moving
      current.row = current.row + move[1]
      current.col = current.col + move[2]

      #writing path
      path.rows=c(path.rows, current.row)
      path.cols=c(path.cols, current.col)

    } else {

      #RANDOM MOVE

      #coin toss
      move=sample(c(1,0))

      #filter values if bounds are reached
      if (current.row == last.row & current.col < last.col){move=c(0,1)}
      if (current.col == last.col & current.row < last.row){move=c(1,0)}

      #moving
      current.row = current.row + move[1]
      current.col = current.col + move[2]

      #writing path
      path.rows = c(path.rows, current.row)
      path.cols = c(path.cols, current.col)

    }

    #stop if last row and col are reached
    if (current.row == last.row & current.col == last.col){break}

  }#end of repeat

  #creating dataframe
  solution=data.frame(a=path.rows, b=path.cols)

  #computing cost
  solution=ComputeDistances(distance.matrix, solution)

  return(solution)

}
