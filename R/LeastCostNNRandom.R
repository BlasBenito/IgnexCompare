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
LeastCostNNRandom=function(distance.matrix, random.threshold=NULL){

  if (random.threshold < 0 | random.threshold > 1){
    stop("random.threshold should be in the interval [0, 1].")
  }

  if(is.null(random.threshold)){random.threshold=0.2}

  current.col=1
  current.row=1
  last.col=ncol(distance.matrix)
  last.row=nrow(distance.matrix)
  path.rows=vector()
  path.cols=vector()
  path.rows[1]=current.row
  path.cols[1]=current.col

  repeat{

    #random number
    random=runif(1)

    #random move or targeted move?
    if (random > random.threshold){

      #TARGETED

      #selecting cell with lower value
      if (current.row < last.row & current.col < last.col){
        if (distance.matrix[current.row+1, current.col] < distance.matrix[current.row, current.col+1])
        {move=c(1,0)} else {move=c(0,1)}
      }

      #filter values if bounds are reached
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
