#' Sequence Slotting
#'
#' It takes as inputs two pollen sequences produced by the function \emph{PrepareInputSequences}, and computes the best slotting.
#'
#' @usage
#'
#' @param cost A distance matrix produced by \emph{DistanceMatrix}.
#' @param diagonal Boolean. If FALSE, the cost is computed
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#' #generating random input data
#' results.table=GenerateResultsTable(10)
#' str(results.table)
#'@export
LeastCost=function(cost, diagonal){

  #setting diagonal if it's empty
  if(is.null(diagonal)){diagonal=FALSE}

  #dimensions
  cost.columns=ncol(cost)
  cost.rows=nrow(cost)

  #matrix to store cumulative cost
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

  #computing cumulative cost
  for (column in 1:(cost.columns-1)){
    for (row in 1:(cost.rows-1)){

      next.row=row+1
      next.column=column+1

      if(diagonal==TRUE){
        cumulative.cost[next.row, next.column] = min(cumulative.cost[row, next.column], cumulative.cost[next.row, column], cumulative.cost[row, column]) + cost[next.row, next.column]
      } else {
        cumulative.cost[next.row, next.column] = min(cumulative.cost[row, next.column], cumulative.cost[next.row, column]) + cost[next.row, next.column]
      }

    } #end of rows
  } #end of columns


  #BACKWARDS PASS TO COMPUTE BEST ALIGNMENT BETWEEN SEQUENCES
  # backwards.pass<-function(cost, cumulative.cost){

  #file to store the case pairings
  pairings<-data.frame(A=nrow(cumulative.cost), B=ncol(cumulative.cost), cost=cost[nrow(cumulative.cost), ncol(cumulative.cost)], cumulative.cost=cumulative.cost[nrow(cumulative.cost), ncol(cumulative.cost)])

  #defining coordinates of the focal cell
  focal.row<-pairings$A
  focal.column<-pairings$B

  #going through the matrix
  repeat{

    #defining values o focal row
    focal.cumulative.cost<-cumulative.cost[focal.row, focal.column]
    focal.cost<-cost[focal.row, focal.column]

    #SCANNING NEIGHBORS

    #dataframe with neighbors
    if(diagonal==TRUE){
    neighbors<-data.frame(A=c(focal.row-1, focal.row-1, focal.row), B=c(focal.column, focal.column-1, focal.column-1))
    } else {
      neighbors<-data.frame(A=c(focal.row-1, focal.row), B=c(focal.column, focal.column-1))
    }


    #removing neighbors with coordinates lower than 1 (out of bounds)
    neighbors[neighbors<1]<-NA
    neighbors<-na.omit(neighbors)
    if(nrow(neighbors)==0){break}

    #computing cost and cumulative cost values for the neighbors
    if(nrow(neighbors)>1){
      neighbors$cost<-diag(cost[neighbors$A, neighbors$B])
      neighbors$cumulative.cost<-diag(x=cumulative.cost[neighbors$A, neighbors$B])
    }else{
      neighbors$cost<-cost[neighbors$A, neighbors$B]
      neighbors$cumulative.cost<-cumulative.cost[neighbors$A, neighbors$B]
    }

    #getting the neighbor with a minimum cumulative.cost
    neighbors<-neighbors[which.min(neighbors$cumulative.cost), c("A", "B")]

    #temporal dataframe to rbind with pairings
    pairings.temp<-data.frame(A=neighbors$A, B=neighbors$B, cost=cost[neighbors$A, neighbors$B], cumulative.cost=cumulative.cost[neighbors$A, neighbors$B])

    #putting them together
    pairings<-rbind(pairings, pairings.temp)

    #decreasing values of focal.row and focal.column
    focal.row<-pairings[nrow(pairings), "A"]
    focal.column<-pairings[nrow(pairings), "B"]

  }#end of repeat

  output<-list()
  output$cumulative.distance<-cumulative.cost[cost.rows, cost.columns]
  output$pairings<-pairings

  return(output)
}



