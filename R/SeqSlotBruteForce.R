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
SeqSlotBruteForce=function(distance.matrix=NULL, iterations=NULL, compute.p.value=NULL){

  #initial checks
  if (is.null(distance.matrix)){
    stop("WOOOSH! No distance matrix provided, aborting.")
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

  #results objects
  results=list()
  best.distances=vector()
  best.solution=data.frame()

  #setting initial value for old cost
  starting.random.walk=LeastCostNNRandom(distance.matrix, random.threshold = 0.1)
  old.cost=sum(starting.random.walk$distances)

  #message
  cat(paste("Generating", iterations, "search paths. I'll be roasting your CPU for a while, will be back soon..."), sep="\n")

  #iterating
  for (i in 1:iterations){
    #generating a random walk
    temp.solution=LeastCostNNRandom(manhattan.matrix, random.threshold = 0.1)

    #computing the new cost
    new.cost=sum(temp.solution$distances)

    #if new.cost is lower or equal than old.cost
    if (new.cost <= old.cost){
      print(paste("New best solution found, cost =", new.cost, sep=" "))
      old.cost=new.cost
      best.solution = temp.solution
      best.distances=c(best.distances, old.cost)
    }
  }#end of iteration
  cat("Done!", sep="\n")

  #PLOTTING
  par(mfrow=c(2,1))
  PlotDistanceMatrix(distance.matrix, main="Best solution", path=best.solution)
  plot(best.distances, type="l", xlab="Solutions", ylab="Cost", main="Improvement")

  #WRITING RESULTS
  results[[1]]=distance.matrix
  results[[2]]=iterations
  results[[3]]=best.solution
  results[[4]]=best.distances[length(best.distances)]
  results[[5]]=best.distances
  results[[6]]="Not computed"

  names(results)=c("distance.matrix", "iterations", "best.solution", "cost.of.best.solution", "costs.of.all.solutions", "p.value")

  #COMPUTING PVALUE
  #################
  if (compute.p.value==TRUE){

    lowest.cost=results[[4]]=best.distances[length(best.distances)]

    #counting results better than best solution
    best.than=0

    #message
    cat(paste("Generating", iterations, "permutations to compute a p-value. I'll be roasting your CPU for a while, AGAIN!"), sep="\n")

    #generating random walks
    for (i in 1:iterations){

      random.walk=ComputeDistances(distance.matrix, RandomWalk(distance.matrix))

      if (sum(random.walk$distances) < lowest.cost){
          best.than=best.than + 1
      }

    }

    #p-value
    results$p.value=best.than/iterations

    cat(paste("Done! p-value =", results$p.value, sep=" "), sep="\n")

    if (results$p.value==0){
      message("WARNING: the p-value is 0. Maybe the solution is awesome, or maybe the number of iterations was not enough!")
    }

  }#end of COMPUTING P VALUE

  return(results)

}
