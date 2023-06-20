from defines import SequenceAlignment

## Define the alignment object using input sequences
alignment = SequenceAlignment("GCATGCG", "GATGA")

###GLOBAL

## Print similarity matrix
alignment.printSimilarityMatrix(type="global")

## Global Reconstruction
alignment.printAlignment(type="local")


###LOCAL
## Print similarity matrix
alignment.printSimilarityMatrix(type="local")

## Global Reconstruction
alignment.printAlignment(type="local")

