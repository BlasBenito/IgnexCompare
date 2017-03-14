#' Sequence Slotting
#'
#' It takes as inputs two pollen sequences produced by the function \emph{PrepareInputSequences}, and computes the best slotting.
#'
#' @usage
#'
#' @param sequences A list produced by the function \emph{PrepareInputSequences} or \emph{DistanceMatrix}
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#' #generating random input data
#' results.table=GenerateResultsTable(10)
#' str(results.table)
#'@export
LeastCost=function(cost){

  #dimensions
  cost.columns=ncol(cost)
  cost.rows=nrow(cost)

  #array to store travel costs
  cumulative.cost=matrix(nrow=cost.rows, ncol=cost.columns)
  rownames(cumulative.cost)=rownames(cost)
  colnames(cumulative.cost)=colnames(cost)

  #first value
  cumulative.cost[1,1]=cost[1,1]
  rownames(cumulative.cost)=rownames(cost)
  colnames(cumulative.cost)=colnames(cost)

  #initiating first column
  cumulative.cost[1, ] = cumsum(cost[1, ])

  #initiating the first row
  cumulative.cost[, 1] = cumsum(cost[, 1])

  #rest of the array
  for (column in 1:(cost.columns-1)){
    for (row in 1:(cost.rows-1)){

      #just for clarity
      next.row=row+1
      next.column=column+1

      #value of the next cell
      cumulative.cost[next.row, next.column] = min(cumulative.cost[row, next.column], cumulative.cost[next.row, column]) + cost[next.row, next.column]

    }
  }

  #distance
  solution=cumulative.cost[cost.rows, cost.columns]

  return(solution)
}
