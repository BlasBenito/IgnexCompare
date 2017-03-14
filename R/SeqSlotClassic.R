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
SeqSlotClassic=function(sequences=NULL, compute.p.value=NULL){

  #initial checks
  if (is.null(sequences)){
    stop("WOOOSH! No input file provided, aborting.")
  }

  #checking if there is a distance matrix in the input object
  if (!is.matrix(sequences$distance.matrix)){
    sequences=DistanceMatrix(sequences=sequences, method="manhattan")
  }

  #default value for compute.p.value
  if (is.null(compute.p.value)){
    compute.p.value=FALSE
  }

  #extracting objects from the input list
  cost=sequences$distance.matrix
  sum.distances.sequence.A=sequences$sum.distances.sequence.A
  sum.distances.sequence.B=sequences$sum.distances.sequence.B

  #input data
  cost=sequences$distance.matrix

  #SLOTTING
  solution=LeastCost(cost)

  #COMPUTING PSI
  #psi.classic
    solution.cost=(solution*2)+(cost[1,1]*2)
      sum.distances.sequences=sum.distances.sequence.A+sum.distances.sequence.B
      if (sum.distances.sequences != 0 & solution.cost !=0){
        psi.classic = (solution.cost - sum.distances.sequences) / sum.distances.sequences
        if (psi.classic < 0.0001){psi.classic=0}
      } else {
        psi.classic = NA
      }

    #psi.modern
    psi.modern=solution/((nrow(cost)+ncol(cost))-1)

  #printing result
  cat(paste("Psi value =", round(psi.classic, 4), sep=" "), sep="\n")

  #COMPUTING P-VALUE
  if (compute.p.value==TRUE){

    #iterating
    for (i in 1:999){

    }

  }

  #WRITING RESULTS
  sequences$psi.classic=psi.classic
  sequences$psi.modern=psi.modern
  sequences$p.value=NA


  return(sequences)

}

