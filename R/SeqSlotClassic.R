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

  #SLOTTING
  solution=LeastCost(unlist(sequences$distance.matrix))

  #COMPUTING PSI
  sequences=ComputePsi(sequences=sequences, slotting.solution=solution)

  #printing result
  cat(paste("Psi value =", round(unlist(sequences$psi.classic), 4), sep=" "), sep="\n")


  #COMPUTING P-VALUE
  ##################################
  if (compute.p.value==TRUE){

    cat("Computing p-value...", sep="\n")

    sequences=ComputePvalue(sequences)


  }

  return(sequences)

}

