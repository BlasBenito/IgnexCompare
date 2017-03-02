# IgnexCompare
Tools to compare pollen sequences (work in progress)

**THIS IS A WORK IN PROGRESS, BE PATIENT!**


##SEQUENCE SLOTTING

```R
#Installing and loading a required package
install.packages("devtools", dep=TRUE)
library(devtools)

#Installing and loading the package IgnexTaxonomy
install_github("IGNEX/IgnexCompare")
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

#Computing manhattan distances among samples and plotting distance matrix
sequences=DistanceMatrix(sequences=sequences, method="manhattan")
PlotDistanceMatrix(sequences$distance.matrix, main="Manhattan distance")

#Compute best slotting using a naive brute force approach
slotting=SeqSlotBruteForce(sequences=sequences, compute.p.value=TRUE, max.random.threshold=0.2)
PlotSlotting(slotting)
slotting$psi #better than the one in the book, the original slotting was wrong.

#Compute best slotting using a parallelized naive brute force approach (for big datasets).
#Requires the pacakges foreach, parallel, and doParallel
slotting.parallel=SeqSlotBruteForceParallel(sequences=sequences, compute.p.value=TRUE, max.random.threshold = 0.2)
PlotSlotting(slotting.parallel)
slotting.parallel$psi
