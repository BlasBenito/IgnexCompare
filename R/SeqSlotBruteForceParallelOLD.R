#' Sequence Slotting by brute force
#'
#' It takes as inputs two pollen sequences produced by the function \emph{PrepareInputSequences}, and computes the best slotting.
#'
#' @usage
#'
#' @param sequences A list produced by the function \emph{PrepareInputSequences} or \emph{DistanceMatrix}
#' @param iterations Number of iterations (Default value is 10000).
#' @param compute.p.value Logical value (TRUE or FALSE). Set to TRUE to compute a p-value based on a permutation test.
#' @param max.random.threshold Number in the interval [0, 1], determining the proportion of the steps defined by chance during the search for the best slotting solution.
#' @param cores
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#'
#' #generating random input data
#' results.table=GenerateResultsTable(10)
#' str(results.table)
#' @export
SeqSlotBruteForceParallelOLD=function(sequences=NULL, iterations=NULL, compute.p.value=NULL, cores=NULL, max.random.threshold=NULL, plot.title=NULL){

  #initial checks
  if (is.null(sequences)){
    stop("WOOOSH! No input data provided, aborting.")
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

  #initial checks
  if (is.null(max.random.threshold)){
    cat("The argument max.random.threshold was set to 0.2.", sep="\n")
    max.random.threshold=0.2
  }

  #checking if there is a distance matrix in the input object
  if ("distance.matrix" %not-in% names(sequences)){
    # message("WARNING: I did not found a distance.matrix object in the input list, but I am very nice, and will compute it for you right away (using the Manhattan method)!")
    sequences=DistanceMatrix(sequences=sequences, method="manhattan")
  }

  #plot title
  if (is.null(plot.title)){
    plot.title="Sequence slotting"
  }



  #getting input data
  distance.matrix=sequences$distance.matrix

  #setting initial value for old cost
  starting.random.walk=LeastCostNNRandom(distance.matrix, max.random.threshold = max.random.threshold)
  starting.cost=sum(starting.random.walk$distances)

  #libraries
  library(foreach)
  library(parallel)
  library(doParallel)

  #default value for number of cores
  if (is.null(cores)){
    cores=detectCores() - 1
    # cat(paste("The argument 'cores' was not provided, I will use", cores, "cores (all but one!).", sep=" "), sep="\n")
  }

  #initiate cluster
  clus=makeCluster(cores)
  registerDoParallel(clus)

  #exporting cluster variables
  clusterExport(cl=clus, varlist=c('SeqSlotBruteForceForParallel', 'LeastCostNNRandom', 'ComputeDistances', 'RandomWalk'), envir=environment())

  #iterations per thread
  iterations.per.thread=round(iterations/cores, 0)

  #parallel execution of SinglePollenTypeToSpecies
  results.foreach=foreach(dummy=1:cores, .combine=cbind) %dopar% SeqSlotBruteForceForParallel(sequences=sequences, iterations=iterations.per.thread, compute.p.value=compute.p.value, starting.cost=starting.cost, max.random.threshold=max.random.threshold)

  stopCluster(clus)


  #PROCESSING results.foreach
  ################################
  n.solutions=ncol(results.foreach)

  #final solution (we get the first column of the list of lists)
  parallel.slotting.solution=results.foreach[, 1]

  # slotting.iterations (sum)
  parallel.slotting.solution$slotting.iterations=sum(unlist(results.foreach["slotting.iterations", ]))

  # lowest costs (c and order)
  parallel.slotting.solution$lowest.costs=rev(sort(unlist(results.foreach["lowest.costs", ])))

  #select column with lower best.slotting.cost, and from that column, get
  best.slotting.values=unlist(results.foreach["best.slotting.cost", ])

  #gets the column of the best slotting (or the first one if all values are the same)
  best.slotting.column=which(best.slotting.values==min(best.slotting.values))[1]

  # best.slotting
  parallel.slotting.solution$best.slotting=data.frame(results.foreach["best.slotting", best.slotting.column])

  # best slotting.cost
  parallel.slotting.solution$best.slotting.cost=unlist(results.foreach["best.slotting.cost", best.slotting.column])

  # psi
  if (is.na(results.foreach["psi", 1])){
    parallel.slotting.solution$psi=NA
  } else {
    parallel.slotting.solution$psi=unlist(results.foreach["psi", best.slotting.column])
  }

  # p.value (sum all p.values)
  if (is.na(results.foreach["p.value", 1])){
    parallel.slotting.solution$p.value=NA
    } else {
    parallel.slotting.solution$p.value=sum(unlist(results.foreach["p.value", ]))
    }

  #plot
  PlotSlotting(slotting=parallel.slotting.solution, main=plot.title)

  return(parallel.slotting.solution)

}


#' @export
SeqSlotBruteForceForParallelOLD=function(sequences, iterations, compute.p.value, dummy=dummy, starting.cost=starting.cost, max.random.threshold=max.random.threshold){

  #generating a random seed
  set.seed(sample(1:1000, size=1))

  #extracting objects from the input list
  distance.matrix=sequences$distance.matrix
  sum.distances.sequence.A=sequences$sum.distances.sequence.A
  sum.distances.sequence.B=sequences$sum.distances.sequence.B

  #results objects
  best.distances=vector()
  best.solution=data.frame()
  old.cost=starting.cost

  #iterating
  for (i in 1:iterations){

    #generating a random walk
    temp.solution=LeastCostNNRandom(distance.matrix, max.random.threshold = max.random.threshold)

    #computing the new cost
    new.cost=sum(temp.solution$distances)

    #if new.cost is lower or equal than old.cost
    if (new.cost <= old.cost){
      old.cost=new.cost
      best.solution=temp.solution
      best.distances=c(best.distances, old.cost)
    }
  }#end of iteration

  #COMPUTING PSI
  best.solution.cost=(sum(best.solution$distances)*2)+(best.solution[1, "distances"]*2)
  sum.distances.sequences=sum.distances.sequence.A+sum.distances.sequence.B
  if (sum.distances.sequences != 0 & best.solution.cost !=0){
    psi = (best.solution.cost - sum.distances.sequences) / sum.distances.sequences
  } else {
    psi = NA
  }

  #WRITTING RESULTS
  #############################################################################
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

    #results objects
    random.distances=vector()

    lowest.cost=sequences$best.slotting.cost

    #counting results better than best solution
    best.than=0

    #generating random walks
    #generating random walks
    for (i in 1:iterations){

      random.walk=ComputeDistances(distance.matrix, RandomWalk(distance.matrix))

      if (((sum(random.walk$distances)*2)+(random.walk[1, "distances"]*2)) < lowest.cost){
        best.than=best.than + 1
      }

    }

    #p-value
    sequences$p.value=best.than/iterations

  }#end of COMPUTING P VALUE

  if (compute.p.value==FALSE){
    sequences$p.value=NA
  }

  return(sequences)

}
