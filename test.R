#development

#update help files
roxygen2::roxygenise()



#test
# install_github("IGNEX/IgnexCompare")
library(IgnexCompare)

#Recreating example of the book "Numerical methods in Quaternary pollen analysis" (Birks and Gordon, 1985)
seqB=SequenceB[c(7,16,25,29,32,36,39), ]
seqA=SequenceA[c(6,11,22,24,29,31,38,39,43,45), ]

#PREPARING INPUT DATA
##################################################
sequences=PrepareSequences(sequence.A=seqA, sequence.B=seqB, sequence.A.name="Abernethy Forest 1970" , sequence.B.name="Abernethy Forest 1974", if.empty.cases="zero", transformation="proportion")
names(sequences)
sequences$metadata
sequences$sequence.A
sequences$sequence.B

#computing distance matrix
sequences=DistanceMatrix(sequences=sequences, method="manhattan")
sequences$distance.matrix

#WITH NO DIAGONALS
#COMPUTING SLOTTING AND p-value THROUGH DISTANCE MATRIX RANDOMIZATION
##########################################################################
sequences.nodiagonal=SequenceSlotting(sequences=sequences, compute.p.value=TRUE, method="manhattan", diagonal=FALSE)
sequences.nodiagonal$psi
sequences.nodiagonal$p.value

#plotting the distance matrix and the best alignment among sequences
PlotSequenceSlotting(sequences.nodiagonal)

#WITH DIAGONALS
#COMPUTING SLOTTING AND p-value THROUGH DISTANCE MATRIX RANDOMIZATION
##########################################################################
sequences.diagonal=SequenceSlotting(sequences=sequences, compute.p.value=TRUE, method="manhattan", diagonal=TRUE)
sequences.diagonal$psi
sequences.diagonal$p.value

#plotting the distance matrix and the best alignment among sequences
PlotSequenceSlotting(sequences.diagonal)


#TESTING CORRELATION BETWEEN PSI VALUES PRODUCED BY DIAGONAL AND NON DIAGONAL ANALYSES
######################################################################################
diagonal.psi<-vector()
nondiagonal.psi<-vector()
rows.A<-vector()
rows.B<-vector()

#repeating analysis
for (iteration in 1:1000){

  #resampling the data (taking a random sample of rows)
  SequenceA.resampled<-SequenceA[sort(sample(nrow(SequenceA), sample(2:nrow(SequenceA), 1))),]
  SequenceB.resampled<-SequenceB[sort(sample(nrow(SequenceB), sample(2:nrow(SequenceB), 1))),]

  #storing number of rows
  rows.A<-c(rows.A, nrow(SequenceA.resampled))
  rows.B<-c(rows.B, nrow(SequenceB.resampled))

  #preparing sequences
  sequences=PrepareInputSequences(sequence.A=SequenceA.resampled, sequence.B=SequenceB.resampled, sequence.A.name="A" , sequence.B.name="B", if.empty.cases="zero", transformation="none")

  #non-diagonal psi
  nondiagonal.psi=c(nondiagonal.psi, SequenceSlotting(sequences=sequences, compute.p.value=FALSE, method="manhattan", diagonal=FALSE, silent=TRUE)$psi)

  #diagonal psi
  diagonal.psi=c(diagonal.psi, SequenceSlotting(sequences=sequences, compute.p.value=FALSE, method="manhattan", diagonal=TRUE, silent=TRUE)$psi)

}#end of iterations

results<-data.frame(diagonal.psi, nondiagonal.psi, rows.A, rows.B, diff.rows=abs(rows.A-rows.B))

pdf("test.pdf", width=10, height=8)
par(mfrow=c(2,2))
plot(density(results$diagonal.psi), col="red", main="Distribution of psi values")
lines(density(results$nondiagonal.psi), col="blue")
legend("topright",c("Diagonal psi", "Non-diagonal psi"), col=c("red", "blue"), lwd=1, cex=0.8)
plot(results$diagonal.psi, results$nondiagonal.psi, xlab="Diagonal psi", ylab="Non-diagonal psi", main="Diagonal vs non-diagonal psi.")
plot(results$diagonal.psi, results$diff.rows, xlab="Diagonal psi", ylab="Rows difference", main="Diagonal psi vs rows difference", col="red")
plot(results$nondiagonal.psi, results$diff.rows, xlab="Non-diagonal psi", ylab="Rows difference", main="Non-diagonal psi vs rows difference", col="blue")
dev.off()





#OTHER OPTIONS
#################################################

#Computing and plotting distance matrices (this step is also done by the function SeqSlotClassic, so it is here as demonstration)
##################################################
#hellinger distance, computed as sqrt(1/2 * sum(sqrt(x)-sqrt(y))^2)
sequences=DistanceMatrix(sequences=sequences, method="hellinger")
PlotDistanceMatrix(sequences$distance.matrix, title="Hellinger distance")

#euclidean distance, computed as sqrt(sum((x - y)^2))
sequences=DistanceMatrix(sequences=sequences, method="euclidean")
PlotDistanceMatrix(sequences$distance.matrix, title="Euclidean distance")


#COMPUTING p-value THROUGH DISTANCE MATRIX RANDOMIZATION USING A DIFFERENT DISTANCE METHOD
slotting.results=SeqSlotClassic(sequences=sequences, compute.p.value=TRUE, method="euclidean")
slotting.results$psi.classic
slotting.results$p.value



########################################################################
########################################################################
#FUNCTION TO APPLY SLOTTING TO MULTIPLE SEQUENCES
########################################################################
########################################################################
#creating example data
sequences.list<-list(SequenceA, SequenceB, SequenceC)
names(sequences.list)<-c("A", "B")

#PREPARING SEQUENCES MULTIPLE
########################################################################
PreparingSequencesMultiple = function(sequences.list=NULL, if.columns.mismatch=NULL, if.empty.cases=NULL, transformation=NULL, silent=NULL){

  #CHECKING ARGUMENTS
  ###################

    #CHECKING sequences.list
    #-----------------------
    #number of sequences
    if (length(sequences.list)==0){
      stop("The list 'sequences.list' seems to be emtpy.")
    }

    if (length(sequences.list)>0 & length(sequences.list)<3){
      stop("There are ")
    }

    #names
    if (is.null(names(sequences.list))){
      stop("The sequences in 'sequences.list' do not have names.")
    }

    if (sum(is.na(names(sequences.list))) > 0){
      stop("Some sequences in 'sequences.list' do not seem to have a name.")
    }

    #structure
    for(dataset in names(sequences.list)){
      dataset.temp<-sequences.list[[dataset]]
      if(is.data.frame(dataset.temp)==FALSE){
        warning(paste("The dataset named ", dataset, " is not a dataframe, I am trying to coerce it into one.", sep=""))
        dataset.temp<-data.frame(dataset.temp)
        if(is.data.frame(dataset.temp)==TRUE){
          cat("It worked out! \n")
          sequences.list[[dataset]]<-dataset.temp
        }
      }
    }


    #CHECKING if.empty.cases
    #-----------------------
    if (is.null(if.empty.cases)){
      message("Warning, argument if.empty.cases is empty, setting its value to 'zero'.")
      if.empty.cases="omit"
      }

    if (if.empty.cases %not-in% c("omit", "zero")){
      stop("Wrong value in the argument 'if.empty.cases'. It can only have the values 'omit' or 'zero'.")
    }


    #CHECKING transformation
    #-----------------------
    if (is.null(transformation)){
      message("Warning, setting the value of 'transformation' to 'none'.")
      transformation="none"
    }

    if (transformation %not-in% c("none", "percentage", "proportion", "hellinger")){
      stop("The 'transformation' argument only accepts the values: 'none', 'percentage', or 'hellinger'.")
    }


    #CHECKING if.columns.mismatch
    #-----------------------
    if (is.null(if.columns.mismatch)){
      message("Warning, setting the value of 'if.columns.mismatch' to 'drop'.")
      if.columns.mismatch="drop"
    }

    if (transformation %not-in% c("drop", "to.zero")){
      stop("The 'if.columns.mismatch' argument only accepts the values: 'drop' and 'to.zero'.")
    }


    #SILENT?
    #-----------------------
    if (is.null(silent)){silent=FALSE}



} #end of PreparingSequencesMultiple

#SEQUENCE SLOTTING
########################################################################
SequenceSlottingMultiple = function(data, groups.column) {

  if (is.null(sampling.multiplier)){sampling.multiplier=10}

  mis.numbers=unique(nMIS)

  psi.matrix=matrix(nrow=length(mis.numbers), ncol=length(mis.numbers))
  psi.sd.matrix=matrix(nrow=length(mis.numbers), ncol=length(mis.numbers))

  for (i in 1:length(mis.numbers)){

    i.mis=mis.numbers[i]

    for (j in  1:length(mis.numbers)){

      if(j>i){next}

      j.mis=mis.numbers[j]

      sequences=PrepareInputSequences(sequence.A=input.data[nMIS==i.mis,], sequence.B=input.data[nMIS == j.mis,], sequence.A.name="sequence.A", sequence.B.name="sequence.B", if.empty.cases="interpolate", output.type="rescaled-proportion", silent=TRUE)


      slotting=SeqSlotEqualSamples(sequences=sequences, sampling.multiplier = sampling.multiplier)

      #write partial result
      #lower part of the matrix
      psi.matrix[i, j]=mean(slotting[["psi.classic"]], na.rm=TRUE)
      psi.sd.matrix[i, j]=sd(slotting[["psi.classic"]])

      #upper part of the matrix
      psi.matrix[j, i]=mean(slotting[["psi.classic"]], na.rm=TRUE)
      psi.sd.matrix[j, i]=sd(slotting[["psi.classic"]])
    }
  }

  #if sd == NA, then 0
  psi.sd.matrix[which(is.na(psi.sd.matrix))] = 0

  rownames(psi.matrix)=paste("MIS", mis.numbers, sep = "")
  colnames(psi.matrix)=paste("MIS", mis.numbers, sep = "")
  rownames(psi.sd.matrix)=paste("MIS", mis.numbers, sep = "")
  colnames(psi.sd.matrix)=paste("MIS", mis.numbers, sep = "")

  results=list(psi.matrix, psi.sd.matrix)
  names(results)=c("average", "deviation")
  return(results)
}

