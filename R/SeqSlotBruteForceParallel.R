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
SeqSlotBruteForceParallel=function(distance.matrix=NULL, iterations=NULL, compute.p.value=NULL, cores=NULL){

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

  #setting initial value for old cost
  starting.random.walk=LeastCostNNRandom(distance.matrix, random.threshold = 0.1)
  starting.cost=sum(starting.random.walk$distances)

  #libraries
  library(foreach)
  library(parallel)
  library(doParallel)

  #default value for number of cores
  if (is.null(cores)){
    cores=detectCores() - 1
    cat(paste("The argument 'cores' was not provided, I will use", cores, "cores (all but one!).", sep=" "), sep="\n")
  }

  #initiate cluster
  clus=makeCluster(cores)
  registerDoParallel(clus)

  #exporting cluster variables
  clusterExport(cl=clus, varlist=c('distance.matrix', 'compute.p.value', 'SeqSlotBruteForceForParallel', 'LeastCostNNRandom', 'ComputeDistances'), envir=environment())

  #iterations per thread
  iterations.per.thread=round(iterations/cores, 0)

  #parallel execution of SinglePollenTypeToSpecies
  results.foreach=foreach(dummy=1:cores, .combine=cbind) %dopar% SeqSlotBruteForceForParallel(distance.matrix=distance.matrix, iterations=iterations.per.thread, compute.p.value=compute.p.value, starting.cost=starting.cost)

  stopCluster(clus)

  return(results.foreach)

}


#' @export
SeqSlotBruteForceForParallel=function(distance.matrix, iterations, compute.p.value, dummy=dummy, starting.cost=starting.cost){

  #generating a random seed
  set.seed(sample(1:1000, size=1))


  #results objects
  results=list()
  best.distances=vector()
  best.solution=data.frame()
  old.cost=starting.cost

  #iterating
  for (i in 1:iterations){

    #generating a random walk
    temp.solution=LeastCostNNRandom(distance.matrix, random.threshold = 0.1)

    #computing the new cost
    new.cost=sum(temp.solution$distances)

    #if new.cost is lower or equal than old.cost
    if (new.cost <= old.cost){
      old.cost=new.cost
      best.solution=temp.solution
      best.distances=c(best.distances, old.cost)
    }
  }#end of iteration

  #WRITING RESULTS
  results[[1]]=distance.matrix
  results[[2]]=iterations
  results[[3]]=best.solution
  results[[4]]=best.distances[length(best.distances)]
  results[[5]]=best.distances
  results[[6]]="Not computed"

  names(results)=c("distance.matrix", "iterations", "best.solution", "cost.of.best.solution", "costs.of.all.solutions", "random.distances")

  #COMPUTING PVALUE
  #################
  if (compute.p.value==TRUE){

    #results objects
    random.distances=vector()

    lowest.cost=results[[4]]=best.distances[length(best.distances)]

    #counting results better than best solution
    best.than=0

    #generating random walks
    for (i in 1:iterations){

      random.walk=ComputeDistances(distance.matrix, RandomWalk(distance.matrix))

      random.distances[i]=sum(random.walk$distances)

    }

    results$random.distances=random.distances

  }#end of COMPUTING P VALUE

  return(results)

}
