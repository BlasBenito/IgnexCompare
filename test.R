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
sequences=PrepareInputSequences(sequence.A=seqA, sequence.B=seqB, sequence.A.name="Abernethy Forest 1970" , sequence.B.name="Abernethy Forest 1974", if.empty.cases="zero", transformation="none")
names(sequences)
sequences$metadata
sequences$sequence.A
sequences$sequence.B


#COMPUTING SLOTTING AND p-value THROUGH DISTANCE MATRIX RANDOMIZATION
##########################################################################
sequences=SeqSlotClassic(sequences=sequences, compute.p.value=TRUE, method="manhattan")
sequences$psi.classic
sequences$p.value

#plotting the distance matrix and the best alignment among sequences
PlotDistanceMatrix(sequences$distance.matrix, title="Manhattan distance")
lines(sequences$pairings$A, sequences$pairings$B)


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
