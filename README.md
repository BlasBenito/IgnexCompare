# IgnexCompare
Tools to compare pollen sequences (work in progress)


##SEQUENCE SLOTTING##

```R
#Installing and loading a required package
install.packages("devtools")
library(devtools)

#Installing and loading the package IgnexTaxonomy
install_github("IGNEX/IgnexCompare")
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

