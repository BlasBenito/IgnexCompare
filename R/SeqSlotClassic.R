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
  cost=unlist(sequences$distance.matrix)
  sum.distances.sequence.A=unlist(sequences$sum.distances.sequence.A)
  sum.distances.sequence.B=unlist(sequences$sum.distances.sequence.B)

  #input data
  cost=sequences$distance.matrix

  #SLOTTING
  solution=LeastCost(cost)

  #COMPUTING PSI
  psi=ComputePsi(distance.matrix=cost, least.cost=solution, autosum.A=sum.distances.sequence.A, autosum.B=sum.distances.sequence.B)

  #printing result
  cat(paste("Psi value =", round(psi[1], 4), sep=" "), sep="\n")

  #WRITING RESULTS
  sequences$psi.classic=psi[1]
  sequences$psi.modern=psi[2]

  #COMPUTING P-VALUE
  ##################################
  if (compute.p.value==TRUE){

    cat("Computing p-value...", sep="\n")

    #initiating number of results better than the real psi value
    best.than=0
    psi.reference=unlist(sequences$psi.classic)
    iterations=(nrow(cost)+ncol(cost))*100

    #iterating
    for (i in 1:iterations){

      # #randomize distance matrix
      # random.matrix=cost[sample(1:nrow(cost), replace=FALSE), ]
      # random.matrix=random.matrix[, sample(1:ncol(random.matrix), replace=FALSE)]

      random.matrix=SwapRowCols(reference.matrix=cost, swaps=1)

      #compute least cost path
      random.solution=LeastCost(cost=random.matrix)

      #compute psi
      psi.random=ComputePsi(distance.matrix=random.matrix, least.cost=random.solution, autosum.A=sum.distances.sequence.A, autosum.B=sum.distances.sequence.B)

      #storing result
      if (psi.random[1] < psi.reference){
        best.than = best.than + 1
      }

    }#end of 1000 iterations

    sequences$p.value=best.than/1000

    cat(paste("P-value =", sequences$p.value, sep=" "), sep="\n")

  }



  sequences$p.value=NA

  return(sequences)

}

