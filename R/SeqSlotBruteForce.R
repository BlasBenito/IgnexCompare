#' Sequence Slotting by brute force
#'
#' description
#'
#' @usage
#'
#' @param param Number of rows of the results table.
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#'
#' #generating random input data
#' results.table=GenerateResultsTable(10)
#' str(results.table)
#' @export
SeqSlotBruteForce=function(sequences=NULL, iterations=NULL, compute.p.value=NULL){

  #initial checks
  if (is.null(sequences)){
    stop("WOOOSH! No input file provided, aborting.")
    }

  if (is.null(iterations)){
    cat("Iterations argument was not provided, using default value (500000).", sep="\n")
    iterations=500000
  }

  #initial checks
  if (is.null(compute.p.value)){
    cat("The argument 'compute.p.value' is empty, using default value (FALSE).", sep="\n")
    compute.p.value=FALSE
  }

  #checking if there is a distance matrix in the input object
  if ("distance.matrix" %not-in% names(sequences)){
    message("WARNING: I did not found a distance.matrix object in the input list, but I am very nice, and will compute it for you right away (using the Manhattan method)!")
    sequences=DistanceMatrix(sequences=sequences, method="manhattan")
  }

  #extracting objects from the input list
  distance.matrix=sequences$distance.matrix
  sum.distances.sequence.A=sequences$sum.distances.sequence.A
  sum.distances.sequence.B=sequences$sum.distances.sequence.B

  #results objects
  best.distances=vector()
  best.solution=data.frame()

  #setting initial value for old cost
  starting.random.walk=LeastCostNNRandom(distance.matrix, random.threshold = 0.1)
  old.cost=sum(starting.random.walk$distances)

  #message
  cat(paste("Generating", iterations, "slottings. I'll be roasting your CPU for a while, will be back soon..."), sep="\n")

  #iterating
  for (i in 1:iterations){
    #generating a random walk
    temp.solution=LeastCostNNRandom(distance.matrix, random.threshold = 0.1)

    #computing the new cost
    new.cost=sum(temp.solution$distances)

    #if new.cost is lower or equal than old.cost
    if (new.cost <= old.cost){
      print(paste("Lowest cost =", new.cost, sep=" "))
      old.cost=new.cost
      best.solution = temp.solution
      best.distances=c(best.distances, old.cost)
    }
  }#end of iteration
  cat("Done!", sep="\n")

  #COMPUTING PSI
  best.solution.cost=best.distances[length(best.distances)]
  sum.distances.sequences=sum.distances.sequence.A+sum.distances.sequence.B
  psi = (best.solution.cost - sum.distances.sequences) / sum.distances.sequences

  cat(paste("The psi value is", round(psi, 4), sep=" "), sep="\n")

  #WRITING RESULTS
  #####################################################################
  previous.names=names(sequences)

  sequences[[9]]=iterations
  sequences[[10]]=best.distances
  sequences[[11]]=best.solution
  sequences[[12]]=best.solution.cost
  sequences[[13]]=psi
  sequences[[14]]="Not computed"

  names(sequences)=c(previous.names, "slotting.iterations", "lowest.costs", "best.slotting", "best.slotting.cost", "psi", "p.value")

  #COMPUTING PVALUE
  #################
  if (compute.p.value==TRUE){

    lowest.cost=sequences$best.slotting.cost

    #counting results better than best solution
    best.than=0

    #message
    cat("Computing p-value...", sep="\n")

    #generating random walks
    for (i in 1:iterations){

      random.walk=ComputeDistances(distance.matrix, RandomWalk(distance.matrix))

      if (sum(random.walk$distances) < lowest.cost){
          best.than=best.than + 1
      }

    }

    #p-value
    sequences$p.value=best.than/iterations

    cat(paste("Done! p-value =", sequences$p.value, sep=" "), sep="\n")

  }#end of COMPUTING P VALUE

  #plot
  PlotSlotting(slotting=sequences)

  return(sequences)

}
