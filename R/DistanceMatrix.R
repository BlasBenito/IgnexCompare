#' Computes a distance matrix
#'
#' Generates a matrix of manhattan or hellinguer distances between the rows of two pollen sequences.
#'
#' @usage
#'
#' @param sequences A list produced by the function \emph{PrepareInputSequences}.
#' @param method Either "manhattan" computed as \emph{sum(abs(x - y))}, "hellinger" computed as \emph{sqrt(1/2 * sum(sqrt(x)-sqrt(y))^2)}, or "euclidean", computed as \emph{sqrt(sum((x - y)^2))}.
#' @return A matrix with manhattan or hellinger distances among cases of the sequences contained in the list \emph{sequences}. The sequence sequence.A will define the rows, and sequence.B will define the columns.
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#' @export
DistanceMatrix <- function(sequences, method = NULL){

  #default value
  if (is.null(method)==TRUE){method <- "manhattan"}

  #getting objects from the list
  sequence.A <- sequences$sequence.A
  sequence.B <- sequences$sequence.B

  #computing row size
  nrow.sequence.A <- nrow(sequence.A)
  nrow.sequence.B <- nrow(sequence.B)

  #creating results matrix
  distance.matrix <- matrix(ncol = nrow.sequence.B, nrow = nrow.sequence.A)

  #COMPUTING DISTANCE MATRIX
  ############################################################################################
  #computing manhattan distance
  if (method %in% c("manhattan", "Manhattan", "MANHATTAN", "man", "Man", "MAN")){
    for (i in 1:nrow.sequence.A){
      for (j in 1:nrow.sequence.B){
        distance.matrix[i,j] <- .ManhattanDistance(x = sequence.A[i,], y = sequence.B[j,])
      }
    }
  }

  #computing hellinger distance
  if (method %in% c("hellinger", "Hellinger", "HELLINGER", "Hell", "hell", "HELL")){
    for (i in 1:nrow.sequence.A){
      for (j in 1:nrow.sequence.B){
        distance.matrix[i,j] <- .HellingerDistance(x = sequence.A[i,], y = sequence.B[j,])
      }
    }
  }

  #computing chi distance
  if (method %in% c("chi", "Chi", "CHI", "chi.distance", "Chi.distance", "CHI.DISTANCE")){
    for (i in 1:nrow.sequence.A){
      for (j in 1:nrow.sequence.B){
        distance.matrix[i,j] <- .ChiDistance(x = sequence.A[i,], y = sequence.B[j,])
      }
    }
  }

  #computing euclidean distance
  if (method %in% c("euclidean", "Euclidean", "EUCLIDEAN", "euc", "Euc", "EUC")){
    for (i in 1:nrow.sequence.A){
      for (j in 1:nrow.sequence.B){
        distance.matrix[i,j] <- .EuclideanDistance(x = sequence.A[i,], y = sequence.B[j,])
      }
    }
  }

  #seting col and row names
  colnames(distance.matrix) <- rownames(sequence.B)
  rownames(distance.matrix) <- rownames(sequence.A)


  #COMPUTING AUTOSUMS
  ############################################################################################
  distances.sequence.A <- vector()
  distances.sequence.B <- vector()

  #computing manhattan distance
  if (method %in% c("manhattan", "Manhattan", "MANHATTAN", "man", "Man", "MAN")){
    for (i in 1:(nrow.sequence.A-1)){
      distances.sequence.A[i] <- .ManhattanDistance(x = sequence.A[i, ], y = sequence.A[i+1, ])
    }

    for (j in 1:(nrow.sequence.B-1)){
      distances.sequence.B[j] <- .ManhattanDistance(x = sequence.B[j, ], y = sequence.B[j+1, ])
    }
  }


  #computing hellinger distance
  if (method %in% c("hellinger", "Hellinger", "HELLINGER", "Hell", "hell", "HELL")){
    for (i in 1:nrow.sequence.A-1){
      distances.sequence.A[i] <- .HellingerDistance(x = sequence.A[i], y = sequence.A[i+1])
    }

    for (j in 1:nrow.sequence.B-1){
      distances.sequence.B[j] <- .HellingerDistance(x = sequence.B[j], y = sequence.B[j+1])
    }
  }

  #computing euclidean distance
  if (method %in% c("euclidean", "Euclidean", "EUCLIDEAN", "euc", "Euc", "EUC")){
    for (i in 1:nrow.sequence.A-1){
      distances.sequence.A[i] <- .EuclideanDistance(x = sequence.A[i], y = sequence.A[i+1])
    }

    for (j in 1:nrow.sequence.B-1){
      distances.sequence.B[j] <- .EuclideanDistance(x = sequence.B[j], y = sequence.B[j+1])
    }
  }

  #computing chi distance
  if (method %in% c("chi", "Chi", "CHI", "chi.squared", "Chi.squared", "CHI.SQUARED")){
    for (i in 1:nrow.sequence.A-1){
      distances.sequence.A[i] <- .EuclideanDistance(x = sequence.A[i], y = sequence.A[i+1])
    }

    for (j in 1:nrow.sequence.B-1){
      distances.sequence.B[j] <- .EuclideanDistance(x = sequence.B[j], y = sequence.B[j+1])
    }
  }

  #SUMMING DISTANCES
  sum.distances.sequence.A <- sum(distances.sequence.A)
  sum.distances.sequence.B <- sum(distances.sequence.B)

  #WRITE RESULTS
  ############################################################################################
  previous.names <- names(sequences)

  #writting new elements in the input list
  sequences$distance.matrix <- distance.matrix
  sequences$sum.distances.sequence.A <- sum.distances.sequence.A
  sequences$sum.distances.sequence.B <- sum.distances.sequence.B

  return(sequences)
}

#check https://www.rdocumentation.org/packages/analogue/versions/0.17-1/source
#' @export
.ManhattanDistance <- function(x, y){
  sum(abs(x - y))
}

#' @export
.EuclideanDistance <- function(x, y){
  sqrt(sum((x - y)^2))
}

#' @export
.ChordDistance <- function(x, y){
  x <- sqrt(x)
  y <- sqrt(y)
  .EuclideanDistance(x, y)
}

#' @export
.ChiDistance <- function(x, y){
  casewise.sum <- x + y
  y <- y / sum(y)
  x <- x / sum(x)
  sqrt(sum(((x - y)^2) / (casewise.sum / sum(casewise.sum))))
}


#' @export
.HellingerDistance  <-  function(x, y){
  #Legendre: the Hellinger distance is the Chord distance computed on square-rooted species abundance data
  #add 0.0001 so there are no zero entries
  x[x==0] <- 0.00001
  y[y==0] <- 0.00001
  #compute distance
  sqrt(1/2 * sum(sqrt(x) - sqrt(y))^2)
}


