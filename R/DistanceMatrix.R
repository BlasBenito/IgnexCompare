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
  distance.matrix = matrix(ncol=nrow.sequence.B, nrow=nrow.sequence.A)

  #COMPUTING DISTANCE MATRIX
  ############################################################################################
  #computing manhattan distance
  if (method=="manhattan"){
    for (i in 1:nrow.sequence.A){
      for (j in 1:nrow.sequence.B){
        distance.matrix[i,j]=ManhattanDistance(sequence.A[i,], sequence.B[j,])
      }
    }
  }

  #computing hellinger distance
  if (method=="hellinger"){
    for (i in 1:nrow.sequence.A){
      for (j in 1:nrow.sequence.B){
        distance.matrix[i,j]=HellingerDistance(sequence.A[i,], sequence.B[j,])
      }
    }
  }

  #seting col and row names
  colnames(distance.matrix)=rownames(sequence.B)
  rownames(distance.matrix)=rownames(sequence.A)


  #COMPUTING AUTOSUMS
  ############################################################################################
  distances.sequence.A=vector()
  distances.sequence.B=vector()

  #computing manhattan distance
  if (method=="manhattan"){
    for (i in 1:(nrow.sequence.A-1)){
      distances.sequence.A[i]=ManhattanDistance(sequence.A[i, ], sequence.A[i+1, ])
    }

    for (j in 1:(nrow.sequence.B-1)){
      distances.sequence.B[j]=ManhattanDistance(sequence.B[j, ], sequence.B[j+1, ])
    }
  }


  #computing manhattan distance
  if (method=="hellinger"){
    for (i in 1:nrow.sequence.A-1){
      distances.sequence.A[i]=HellingerDistance(sequence.A[i], sequence.A[i+1])
    }

    for (j in 1:nrow.sequence.B-1){
      distances.sequence.B[j]=HellingerDistance(sequence.B[j], sequence.B[j+1])
    }
  }

  #SUMMING DISTANCES
  sum.distances.sequence.A=sum(distances.sequence.A)
  sum.distances.sequence.B=sum(distances.sequence.B)

  #WRITE RESULTS
  ############################################################################################
  previous.names=names(sequences)

  #writting new elements in the input list
  sequences[[5]]=method
  sequences[[6]]=distance.matrix
  sequences[[7]]=sum.distances.sequence.A
  sequences[[8]]=sum.distances.sequence.B

  #new names
  names(sequences)=c(previous.names, "distance.method", "distance.matrix", "sum.distances.sequence.A", "sum.distances.sequence.B")

  return(sequences)
}
