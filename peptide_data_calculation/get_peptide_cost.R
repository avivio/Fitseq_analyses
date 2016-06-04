#these are the two libraries needed for this script, you're going to need to install them to run this
require(seqinr)
require(Peptides)
#load the amino acid cost table 
data(aacost)

#read the gene amino acid sequence table
gene.sequences = read.csv('C:\\Users\\dell7\\Documents\\Tzachi\\analysis\\R_project\\peptide amino acid sequences.csv')


#a function to calculate the cost of each sequence by summing each amino acid
get_pep_cost <- function(sequence){
  #split the sequnce to amino acids
  sequence.split <- strsplit(as.character(sequence), "")[[1]]
  #initate cost variable
  tot.cost = 0
  #loop over amino acid sequence, summing up the cost along the way
  for (aa in sequence.split){
    cost = aacost[aacost[,2] == aa,'tot']
    tot.cost = cost + tot.cost
  }
  return(tot.cost)
}
#vectorize the function
get_pep_cost <- Vectorize(get_pep_cost)
#run to get cost for entire set of genes
gene.sequences$cost <- get_pep_cost(gene.sequences$pep.sequence)

# write data
write.csv(aacost,'C://Users//dell7//Documents//Tzachi//workspace//data//peptide_cost.csv', row.names = F)
