#test
# install_github("IGNEX/IgnexCompare")
library(IgnexCompare)

#Recreating example of the book "Numerical methods in Quaternary pollen analysis" (Birks and Gordon, 1985)
seqB=SequenceB[c(7,16,25,29,32,36,39), ]
seqA=SequenceA[c(6,11,22,24,29,31,38,39,43,45), ]

#Preparing input data
sequences=PrepareInputSequences(sequence.A=seqA, sequence.B=seqB, sequence.A.name="Abernethy Forest 1970" , sequence.B.name="Abernethy Forest 1974", if.empty.cases="interpolate", output.type="proportion")
names(sequences)
sequences$metadata
sequences$sequence.A
sequences$sequence.B

#Computing manhattan distances among samples and plotting distance matrix (this step is also done by the function SeqSlotClassic, so it is here as demonstration)
sequences=DistanceMatrix(sequences=sequences, method="manhattan")
PlotDistanceMatrix(sequences$distance.matrix, main="Manhattan distance")

#Compute sequence slotting
slotting.results=SeqSlotClassic(sequences=sequences)
slotting.results$psi.classic
