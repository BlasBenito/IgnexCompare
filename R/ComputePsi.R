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
ComputePsi=function(distance.matrix=NULL, least.cost=NULL, autosum.A=NULL, autosum.B=NULL){

  if (is.null(distance.matrix) | is.null(least.cost) | is.null(autosum.A) | is.null(autosum.B)){
    stop("Arguments missing! You need to define 'distance.matrix', 'least.cost', 'autosum.A', and 'autosum.B'")
  }

  #psi.classic
  solution.cost=(least.cost*2)+(distance.matrix[1,1]*2)
  sum.distances.sequences=autosum.A+autosum.B
  if (sum.distances.sequences != 0 & solution.cost !=0){
    psi.classic = (solution.cost - sum.distances.sequences) / sum.distances.sequences
    if (psi.classic < 0.0001){psi.classic=0}
  } else {
    psi.classic = NA
  }

  #psi.modern
  psi.modern=least.cost/((nrow(distance.matrix)+ncol(distance.matrix))-1)

  #putting results together
  result=c(psi.classic, psi.modern)

  return(result)

}
