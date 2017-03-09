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
SeqSlotEqualSampleSize=function(sequences, psi.mode="classic"){

  #initial checks
  if (is.null(sequences)){
    stop("WOOOSH! No input file provided, aborting.")
  }

  #checking if there is a distance matrix in the input object
  if ("distance.matrix" %not-in% names(sequences)){
    #message
    sequences=DistanceMatrix(sequences=sequences, method="manhattan")
  }

  #extracting objects from the input list
  cost=sequences$distance.matrix

  #what sequence is longer?
  sequence.A.nrow=nrow(sequences$sequence.A)
  sequence.B.nrow=nrow(sequences$sequence.B)

  #defining small and big sequence
  if (sequence.A.nrow > sequence.B.nrow){
    small.sequence.name="sequence.B"
    big.sequence.name="sequence.A"
  }

  if (sequence.A.nrow < sequence.B.nrow){
    small.sequence.name="sequence.A"
    big.sequence.name="sequence.B"
  }


  #getting the sequences
  small.sequence.data=sequences[[small.sequence.name]]
  small.sequence.nrow=nrow(small.sequence.data)
  big.sequence.data=sequences[[big.sequence.name]]
  big.sequence.nrow=nrow(big.sequence.data)

  #computing autosum of small sequence
  small.sequence.autosum=vector()
  for (i in 1:(small.sequence.nrow-1)){
    small.sequence.autosum[i]=ManhattanDistance(small.sequence.data[i, ], small.sequence.data[i+1, ])
  }
  small.sequence.autosum=sum(small.sequence.autosum)

  #computing number of samples to draw
  # samples=sqrt(choose(big.sequence.nrow, small.sequence.nrow))
  samples=1000

  #vector to store results
  psi.values=vector()

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
    for (i in 1:(big.sequence.nrow-1)){
      big.sequence.autosum[i]=ManhattanDistance(big.sequence.data[i, ], big.sequence.data[i+1, ])
    }
    big.sequence.autosum=sum(big.sequence.autosum)


    #COMPUTING SLOTTING
    ######################################################################

    #creating matrix to store cumulative cost
    cumulative.cost=matrix(ncol=small.sequence.nrow, nrow=small.sequence.nrow)
    rownames(cumulative.cost)=rownames(cost.temp)
    colnames(cumulative.cost)=colnames(cost.temp)

    #first value
    cumulative.cost[1,1]=cost.temp[1,1]

    #initiating first column
    cumulative.cost[1, ] = cumsum(cost.temp[1, ])

    #initiating the first row
    cumulative.cost[, 1] = cumsum(cost.temp[, 1])

    #rest of the array
    for (column in 1:(small.sequence.nrow-1)){
      for (row in 1:(small.sequence.nrow-1)){

        #just for clarity
        next.row=row+1
        next.column=column+1

        #value of the next cell
        cumulative.cost[next.row, next.column] = min(cumulative.cost[row, next.column], cumulative.cost[next.row, column]) + cost.temp[next.row, next.column]

      }
    }

    #distance
    solution=cumulative.cost[small.sequence.nrow, small.sequence.nrow]

    #COMPUTING PSI
    if (psi.mode=="classic"){
      solution.cost=(solution*2)+(cost.temp[1,1]*2)
      sum.distances.sequences=big.sequence.autosum+small.sequence.autosum

      if (sum.distances.sequences != 0 & solution.cost != 0){
        psi = (solution.cost - sum.distances.sequences) / sum.distances.sequences

      }

      if (sum.distances.sequences == 0 & solution.cost == 0) {
        psi = NA
      }

      #in some cases psi is only close to zero, I have to check why
      if (psi < 0.0001){
        psi=0
      }
    }

    if (psi.mode=="modern"){
      psi=solution/((cost.temp.rows+cost.temp.columns)-1)
    }
    ######################################################################

    #writing result
    psi.values[sample]=psi

  }#END OF ITERATION THROUGH SAMPLES

  #WRITING RESULTS (for compatibility with the other seqslot functions)
  #####################################################################
  previous.names=names(sequences)

  sequences[[9]]=samples
  sequences[[11]]=solution
  sequences[[12]]=NA
  sequences[[13]]=psi.values
  sequences[[14]]=NA

  names(sequences)=c(previous.names, "iterations", "lowest.costs", "best.slotting", "best.slotting.cost", "psi", "p.value")

  return(sequences)

}
