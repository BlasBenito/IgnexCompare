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
ComputePvalue=function(sequences=NULL){

  if (is.null(sequences)){
    stop("Input list not provided!")
  }

  #initiating number of results better than the real psi value
  better.than.reference.value=0
  psi.reference=unlist(sequences$psi.classic)

  #getting the distance matrix
  reference.distance.matrix=unlist(sequences$distance.matrix)
  reference.distance.matrix.nrow=nrow(reference.distance.matrix)
  reference.distance.matrix.ncol=ncol(reference.distance.matrix)

  #number of iterations to compute p-value relative to the size of the matrix
  iterations=(reference.distance.matrix.nrow + reference.distance.matrix.ncol) *100

  #copy of sequences
  sequences.temp=sequences

  #iterating
  for (i in 1:iterations){

    # #randomize distance matrix
    # random.matrix=reference.distance.matrix[sample(1:reference.distance.matrix.nrow, replace=FALSE), ]
    # random.matrix=random.matrix[, sample(1:reference.distance.matrix.ncol, replace=FALSE)]

    #swap two adjacent columns and rows in the matrix (and include it in a fake sequences list)
    random.matrix=.SwapRowCols(reference.distance.matrix)
    sequences.temp$distance.matrix=random.matrix

    #compute least cost path
    random.solution=LeastCost(cost=random.matrix)

    #compute psi
    sequences.temp=ComputePsi(sequences=sequences.temp, slotting.solution=random.solution)

    #storing result
    if (unlist(sequences.temp$psi.classic) < psi.reference){
      better.than.reference.value = better.than.reference.value + 1
    }

  }#end of 1000 iterations

  sequences$p.value=better.than.reference.value/iterations

  cat(paste("p-value =", round(sequences$p.value, 3), sep=" "), sep="\n")

  return(sequences)
}


#' @export
.SwapRowCols=function(reference.distance.matrix=NULL){

  if (is.null(reference.distance.matrix)){stop("No reference matrix provided.")}

  #computing number of swaps (a number between 1 and 1/5th of the minimum number of column or rows of the matrix)
  swaps=max(dim(reference.distance.matrix))
  if (swaps==0){swaps=1}

  #generating the starting matrix
  random.matrix=reference.distance.matrix

  #swapping one column and one row each time
  #indexes of columns to swap (we don't want the first and the last)
  columns.range=2:(ncol(reference.distance.matrix)-1)
  rows.range=2:(nrow(reference.distance.matrix)-1)

  #left or right
  left.or.right=c(-1, 1)

  #swapping rows and columns
  for (i in 1:swaps){

    #choose column to swap
    column.to.swap=sample(columns.range, 1)
    row.to.swap=sample(rows.range, 1)

    #swapping rows
    random.matrix[row.to.swap + 1, ] = reference.distance.matrix[row.to.swap, ]
    random.matrix[row.to.swap, ] = reference.distance.matrix[row.to.swap - 1, ]

    #swapping columns
    random.matrix[, column.to.swap + 1] = reference.distance.matrix[, column.to.swap]
    random.matrix[, column.to.swap] = reference.distance.matrix[, column.to.swap - 1]

  }

  return(random.matrix)
}
