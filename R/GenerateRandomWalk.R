#' Generates a random walk.
#'
#' Generates a random walk based on a distance matrix.
#' Work in progres, put marked code in another function.
#'
#' @usage
#'
#' @param
#' @param
#' @return Data frame, the input random walk with a new column named \emph{distance}
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#' @export
GenerateRandomWalk=function(distance.matrix){

  dm=distance.matrix
  dm.dim=dim(dm)

  #THIS CODE HAS TO GO TO ANOTHER FUNCTION THAT READS THE DISTANCE MATRIX, AND RETURNS A LIST WITH THE DETAILS
  ############################################################################################################
  #NOTE
  #SEQUENCE A IS IN THE COLUMNS
  #SEQUENCE B IS IN THE ROWS

  #bounds of the distance matrix
  a.min=1
  b.min=1
  b.max=dm.dim[1]
  a.max=dm.dim[2]

  #coordinates
  a.coord=a.min:a.max
  b.coord=b.min:b.max

  #base distributions of 1 and 0
  a.base=c(rep(1, a.max-1), rep(0, b.max))
  b.base=c(rep(1, b.max-1), rep(0, a.max))
  ############################################################################

  #THIS SHOULD COME FROM A LIST, TO AVOID READING THE DISTANCE MATRIX EVERY TIME WE CREATE A RANDOM WALK
  #randomized moves (add 0 at the end because we cannot move beyond the last coordinates)
  a.random=c(1, sample(a.base))
  b.random=c(1, sample(b.base))

  #cumulative sum (new coordinates, but with duplicates where a and b are 0 simultaneously)
  a.random.cumsum=cumsum(a.random)
  b.random.cumsum=cumsum(b.random)

  #cumulative sum of the movement vectors
  # random.walk=data.frame(a=a.random.cumsum, b=b.random.cumsum, temp=dm[cbind(a.random.cumsum, b.random.cumsum)])
  random.walk=data.frame(a=a.random.cumsum, b=b.random.cumsum)

  #removing duplicates (where coordinates do not change becaues of consecutive zeros)
  random.walk=random.walk[!duplicated(random.walk), ]

return(random.walk)

}
