#' Sequence Slotting
#'
#' It takes as inputs two pollen sequences produced by the function \emph{PrepareInputSequences}, and computes the best slotting.
#'
#' @usage
#'
#' @param sequences A list produced by the function \emph{PrepareInputSequences} or \emph{DistanceMatrix}
#' @param compute.p.value Boolean. Use TRUE to compute p-values via distance matrix randomization.
#' @param method Either "manhattan" computed as \emph{sum(abs(x - y))}, "hellinger" computed as \emph{sqrt(1/2 * sum(sqrt(x)-sqrt(y))^2)}, or "euclidean", computed as \emph{sqrt(sum((x - y)^2))}.
#' @param diagonal Boolean. TRUE to compute the sequence slotting using diagonals (true shortest path). FALSE to compute the sequence slotting without considering diagonals. Using diagonals produces lower psi values but better alignments between sequences.
#' @param silent Boolean. Use TRUE to remove all messages.
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#' #generating random input data
#' results.table=GenerateResultsTable(10)
#' str(results.table)
#' @export
SequenceSlotting=function(sequences=NULL, compute.p.value=NULL, method=NULL, silent=NULL, diagonal=NULL){

  #SILENT?
  if (is.null(silent)){silent=FALSE}

  #initial checks
  if (is.null(sequences)){
    stop("WOOOSH! No input file provided, aborting.")
  }

  #if method is null, set manhattan as default
  if (is.null(method)==TRUE){method="manhattan"}

  #checking if there is a distance matrix in the input object or the user wants to change the method to compute the distance matrix
  if (!is.matrix(sequences$distance.matrix)){
    sequences=DistanceMatrix(sequences=sequences, method=method)
  }

  #default value for compute.p.value
  if (is.null(compute.p.value)){
    compute.p.value=FALSE
  }

  #SLOTTING
  LeastCost.solution=LeastCost(cost=unlist(sequences$distance.matrix), diagonal=diagonal)

  #COMPUTING PSI
  sequences=ComputePsi(sequences=sequences, slotting.solution=LeastCost.solution$cumulative.distance)

  #printing result
  if(silent==FALSE){
  cat(paste("Psi value =", round(unlist(sequences$psi), 4), sep=" "), sep="\n")
  }


  #COMPUTING P-VALUE
  ##################################
  if (compute.p.value==TRUE){

    if(silent==FALSE){
    cat("Computing p-value...", sep="\n")
    }

    sequences=ComputePvalue(sequences, diagonal=diagonal)


  }

  #adding best pairings to sequences
  sequences$pairings<-LeastCost.solution$pairings

  return(sequences)

}

