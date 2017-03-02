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
#' @author Blas Benito <blasbenito@gmail.com>
#' @examples
#'
#' #generating random input data
#' results.table=GenerateResultsTable(10)
#' str(results.table)
#' @export
SeqSlotBruteForce=function(sequences=NULL, compute.p.value=NULL, max.random.threshold=NULL){

  #initial checks
  if (is.null(sequences)){
    stop("WOOOSH! No input file provided, aborting.")
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
    #message
    cat("Computing distance matrix...", sep="\n")
    sequences=DistanceMatrix(sequences=sequences, method="manhattan")
  }

  #computing convergence criteria (minimum number of times a solution has to be found to be valid).
  convergence.criterion=nrow(sequences$sequence.A)*nrow(sequences$sequence.B)*100

  #extracting objects from the input list
  distance.matrix=sequences$distance.matrix
  sum.distances.sequence.A=sequences$sum.distances.sequence.A
  sum.distances.sequence.B=sequences$sum.distances.sequence.B

  #results objects
  best.costs=vector()
  best.solution=data.frame()

  #setting initial value for best.distances
  starting.random.walk=LeastCostNNRandom(distance.matrix, max.random.threshold = max.random.threshold)
  slotting.steps=nrow(starting.random.walk)
  best.costs=c(best.costs, sum(starting.random.walk$distances))

  #starting convergence
  convergence=0

  #message
  cat("Computing optimal slotting...", sep="\n")

  iterations.to.convergence=0
  iterations=0

  #iterating
  repeat {

    iterations.to.convergence=iterations.to.convergence+1
    iterations=iterations+1

    #generating a random walk
    temp.solution=LeastCostNNRandom(distance.matrix, max.random.threshold = max.random.threshold)

    #computing the new cost
    new.cost=sum(temp.solution$distances)

      #new cost lower than the last cost added to best.distances
      if (new.cost <= best.costs[length(best.costs)]){

        cat(paste("New solution found - cost = ", round(new.cost/slotting.steps, 3), "; iterations = ", iterations, "; random threshold = ", max.random.threshold, sep=" "), sep="\n")

        #resetting convergence
        iterations.to.convergence=0

        #diminishing max.random.threshold
        max.random.threshold = max.random.threshold - 0.01

        #storing solution
        best.solution=temp.solution

        #storing cost
        best.costs=c(best.costs, new.cost)

      }

    if (iterations.to.convergence >= convergence.criterion){
      cat(paste(iterations.to.convergence, "iterations without finding a better solution. I'm done!", sep=" "), sep="\n")
      break
    }

  }#end of while

  #message
  cat("Sequence slotting done!", sep="\n")

  #COMPUTING PSI
  best.solution.cost=(sum(best.solution$distances)*2)+(best.solution[1, "distances"]*2)
  sum.distances.sequences=sum.distances.sequence.A+sum.distances.sequence.B
  if (sum.distances.sequences != 0 & best.solution.cost !=0){
    psi = (best.solution.cost - sum.distances.sequences) / sum.distances.sequences
  } else {
    psi = NA
  }

  cat(paste("Psi value =", round(psi, 4), sep=" "), sep="\n")

  #NORMALIZING BEST DISTANCES BY THE NUMBER OF SLOTTING STEPS
  best.costs=best.costs/slotting.steps

  #WRITING RESULTS
  #####################################################################
  previous.names=names(sequences)

  sequences[[9]]=iterations
  sequences[[10]]=best.costs
  sequences[[11]]=best.solution
  sequences[[12]]=min(best.costs)[1]
  sequences[[13]]=psi
  sequences[[14]]="Not computed"

  names(sequences)=c(previous.names, "iterations", "lowest.costs", "best.slotting", "best.slotting.cost", "psi", "p.value")

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

      if (((sum(random.walk$distances)*2)+(random.walk[1, "distances"]*2)) < lowest.cost){
          best.than=best.than + 1
      }

    }

    #p-value
    sequences$p.value=best.than/iterations

    cat(paste("P-value =", sequences$p.value, sep=" "), sep="\n")

  }#end of COMPUTING P VALUE

  cat("I am done!", sep="\n")

  return(sequences)

}
