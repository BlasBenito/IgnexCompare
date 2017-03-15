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
SeqSlotEqualSamples=function(sequences=NULL, sampling.multiplier=NULL){

  #Default value for sampling.multiplier
  if (is.null(sampling.multiplier)){
    sampling.multiplier=1
  }

  #initial checks
  if (is.null(sequences)){
    stop("WOOOSH! No input file provided, aborting.")
  }

  #checking if there is a distance matrix in the input object
  if (!is.matrix(unlist(sequences$distance.matrix))){
    #message
    sequences=DistanceMatrix(sequences=sequences, method="manhattan")
  }

  #extracting objects from the input list
  cost=unlist(sequences$distance.matrix)

  #what sequence is longer?
  sequence.A.nrow=nrow(sequences[["sequence.A"]])
  sequence.B.nrow=nrow(sequences[["sequence.B"]])

  #defining small and big sequence
  if (sequence.A.nrow > sequence.B.nrow){
    small.sequence.name="sequence.B"
    big.sequence.name="sequence.A"
  }

  if (sequence.A.nrow < sequence.B.nrow){
    small.sequence.name="sequence.A"
    big.sequence.name="sequence.B"
  }

  #if the sequences have the same size
  if (sequence.A.nrow == sequence.B.nrow){
    small.sequence.name="sequence.A"
    big.sequence.name="sequence.B"
    samples=1
  }


  #getting the sequences
  small.sequence.data=sequences[[small.sequence.name]]
  small.sequence.nrow=nrow(small.sequence.data)
  big.sequence.data=sequences[[big.sequence.name]]
  big.sequence.nrow=nrow(big.sequence.data)


  #computing number of samples to draw
  samples = (big.sequence.nrow - small.sequence.nrow) * sampling.multiplier

  #computing autosum of small sequence
  small.sequence.autosum=vector()
  for (i in 1:(small.sequence.nrow-1)){
    small.sequence.autosum[i]=.ManhattanDistance(small.sequence.data[i, ], small.sequence.data[i+1, ])
  }
  small.sequence.autosum=sum(small.sequence.autosum)



  #vector to store results
  psi.classic.values=vector()
  psi.modern.values=vector()

  #iterating through samples
  ################################
  for (sample in 1:samples){

    #generating a sample of the big sequence
    sampling=sort(sample(1:big.sequence.nrow, small.sequence.nrow))

    #subsampling the big sequence
    big.sequence.temp=big.sequence.data[sampling, ]

    #subsetting the cost matrix
    if (big.sequence.name=="sequence.A"){
      cost.temp=cost[sampling, ]
    }

    if (big.sequence.name=="sequence.B"){
      cost.temp=cost[, sampling]
    }

    #computing autosum of big sequence
    big.sequence.autosum=vector()
    for (i in 1:(small.sequence.nrow-1)){
      big.sequence.autosum[i]=.ManhattanDistance(big.sequence.temp[i, ], big.sequence.temp[i+1, ])
    }
    big.sequence.autosum=sum(big.sequence.autosum)


    #COMPUTING SLOTTING
    ######################################################################
    solution=LeastCost(cost.temp)

    #COMPUTING PSI
    #COMPUTING PSI
    psi=ComputePsi(distance.matrix=cost.temp, least.cost=solution, autosum.A=big.sequence.autosum, autosum.B=small.sequence.autosum)

    #writing result
    psi.classic.values[sample]=psi[1]
    psi.modern.values[sample]=psi[2]
  }#END OF ITERATION THROUGH SAMPLES

  #WRITING RESULTS (for compatibility with the other seqslot functions)
  #####################################################################
  previous.names=names(sequences)

  sequences$psi.classic=psi.classic.values
  sequences$psi.modern=psi.modern.values

  cat(paste("Psi average = ", round(mean(psi.classic.values), 3), sep=""), sep="\n")
  cat(paste("Psi standard deviation = ", round(sd(psi.classic.values), 3), sep=""), sep="\n")

  return(sequences)

}
