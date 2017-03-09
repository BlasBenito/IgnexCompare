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
#' @export
SeqSlotClassic=function(sequences, psi.mode=NULL){

  if(is.null(psi.mode)){
    psi.mode="classic"
  }

  #initial checks
  if (is.null(sequences)){
    stop("WOOOSH! No input file provided, aborting.")
  }

  #checking if there is a distance matrix in the input object
  if ("distance.matrix" %not-in% names(sequences)){
    #message
    sequences=DistanceMatrix(sequences=sequences, method="manhattan")
  }

  #extracting objects from the input list
  cost=sequences$distance.matrix
  sum.distances.sequence.A=sequences$sum.distances.sequence.A
  sum.distances.sequence.B=sequences$sum.distances.sequence.B

  #input data
  cost=sequences$distance.matrix
  cost.columns=ncol(cost)
  cost.rows=nrow(cost)

  #array to store travel costs
  cumulative.cost=matrix(nrow=cost.rows, ncol=cost.columns)

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

  #COMPUTING PSI
  if (psi.mode=="classic"){
    solution.cost=(solution*2)+(cost[1,1]*2)
      sum.distances.sequences=sum.distances.sequence.A+sum.distances.sequence.B
      if (sum.distances.sequences != 0 & solution.cost !=0){
        psi = (solution.cost - sum.distances.sequences) / sum.distances.sequences
        if (psi < 0.0001){psi=0}
      } else {
        psi = NA
      }
  }

  if (psi.mode="modern"){
    psi=solution/((cost.rows+cost.columns)-1)
  }

  #printing result
  cat(paste("Psi value =", round(psi, 4), sep=" "), sep="\n")

  #WRITING RESULTS (for compatibility with the other seqslot functions)
  #####################################################################
  previous.names=names(sequences)

  sequences[[9]]=NA
  sequences[[10]]=NA
  sequences[[11]]=solution
  sequences[[12]]=NA
  sequences[[13]]=psi
  sequences[[14]]=NA

  names(sequences)=c(previous.names, "iterations", "lowest.costs", "best.slotting", "best.slotting.cost", "psi", "p.value")

  return(sequences)

}
