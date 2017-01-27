#' Computes a distance matrix
#'
#' Generates a matrix of manhattan or hellinguer distances between the rows of two pollen sequences.
#'
#' @usage
#'
#' @param sequences A list produced by the function \emph{PrepareInputSequences}.
#' @param method Either "manhattan" (default) or "hellinger".
#' @return A matrix with manhattan or hellinger distances among cases of the sequences contained in the list \emph{sequences}. The sequence sequence.A will define the rows, and sequence.B will define the columns.
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#' @export
DistanceMatrix=function(sequences, method=NULL){

  #default value
  if (is.null(method)==TRUE){method="manhattan"}

  #error in method name
  if (method!="manhattan" & method!="hellinger"){
    stop("The 'method' name should be either 'manhattan' or 'hellinger'")
  }

  #getting objects from the list
  sequence.A=sequences$sequence.A
  sequence.B=sequences$sequence.B

  #computing row size
  nrow.sequence.A=nrow(sequence.A)
  nrow.sequence.B=nrow(sequence.B)

  #creating results matrix
  result = matrix(ncol=nrow.sequence.B, nrow=nrow.sequence.A)

  #computing manhattan distance
  if (method=="manhattan"){
    for (i in 1:nrow.sequence.A){
      for (j in 1:nrow.sequence.B){
        result[i,j]=ManhattanDistance(sequence.A[i,], sequence.B[j,])
      }
    }
  }

  #computing hellinger distance
  if (method=="hellinger"){
    for (i in 1:nrow.sequence.A){
      for (j in 1:nrow.sequence.B){
        result[i,j]=HellingerDistance(sequence.A[i,], sequence.B[j,])
      }
    }
  }

  #seting col and row names
  colnames(result)=rownames(sequence.B)
  rownames(result)=rownames(sequence.A)

  return(result)
}
