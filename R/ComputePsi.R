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
ComputePsi=function(sequences=NULL, slotting.solution=NULL){

  if (is.null(sequences)){
    stop("Argument 'sequences' is missing!")
  }

  if (is.null(slotting.solution)){
    stop("Argument 'slotting.solution' is missing!")
  }

  #getting the distance matrix
  distance.matrix=unlist(sequences$distance.matrix)

  #cost of the best solution computed as in the original code
  solution.cost=(slotting.solution*2)+(distance.matrix[1,1]*2)

  #autosum of sequences
  sum.distances.sequences=unlist(sequences$sum.distances.sequence.A) + unlist(sequences$sum.distances.sequence.B)

  #psi.classic
  if (sum.distances.sequences != 0 & solution.cost !=0){
    psi.classic = (solution.cost - sum.distances.sequences) / sum.distances.sequences
    if (psi.classic < 0.0001){psi.classic=0}
  } else {
    psi.classic = NA
  }

  #psi.modern
  psi.modern=slotting.solution/((nrow(distance.matrix)+ncol(distance.matrix))-1)

  #WRITING RESULTS
  sequences$psi.classic=psi.classic
  sequences$psi.modern=psi.modern

  return(sequences)

}
