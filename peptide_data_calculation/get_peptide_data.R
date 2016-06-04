#these are the two libraries needed for this script, you're going to need to install them to run this
require(seqinr)
require(Peptides)

get.peptide.stats <- function(sequence){
  # a function for calculating the statistics for a specific peptide
  #load tables containing peptide data
  data(aacost)
  aademand  <- read.csv('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\aa_demand.csv')
#initate peptide data list
pep.properties <- list()
#split the peptide into a list of single aa
sequence.split <- strsplit(sequence, "")[[1]]
#initate total cost and total deman
tot.cost = 0
tot.demand = 0
#loop over amino acids get cost/demand from table and add to running total
for (aa in sequence.split){
  cost = aacost[aacost[,2] == aa,'tot']
  tot.cost = cost + tot.cost
  demand = aademand[aademand$aa == aa,'demand']
  tot.demand = demand + tot.demand
  
}
#put total cost and total demand into pep.properties
pep.properties$pep.cost <- tot.cost
pep.properties$pep.demand <- tot.demand

#get amino acid stats from the seqinr package
stats <- AAstat(sequence.split,FALSE)

#get pki from the aastat stats
pep.properties$pep.pi <- stats[[3]]

#get molecular weight, aliphatic index, boman index, peptide charge, hydrophobic moment, hydrophobicity, and instability index
#from the Peptides package and store in pep.properties
pep.properties$pep.mw <- mw(sequence)
pep.properties$pep.aindex <- aindex(sequence)
pep.properties$pep.boman <- boman(sequence)
pep.properties$pep.charge <- charge(sequence)
pep.properties$pep.hmoment<- hmoment(sequence)
pep.properties$pep.hydro <- hydrophobicity(sequence)
pep.properties$pep.instability <- instaindex(sequence)

#turn pep.properties into a data frame
pep.properties <- data.frame(pep.properties)
#get the aa frequency from aastat stats, convert to a data frame and create the names vector
#by converting the amino acid value to full name
aa.freq <- as.list(stats[[1]]/11)
aa.freq <- data.frame(aa.freq)[-1]
names(aa.freq) <- paste0('pep.',aaa(names(aa.freq)))

#get amino acid classes percent from aastat data and convert names to lower class
aa.class.percent <- data.frame(stats[[2]])
names(aa.class.percent) <- paste0('pep.',tolower(names(aa.class.percent)))

#put the frequecny anc classes into pep.properties
pep.properties <- bind_cols(pep.properties,aa.freq)
pep.properties <- bind_cols(pep.properties,aa.class.percent)

return(pep.properties)

}
#vectorize the function
get.peptide.stats <- Vectorize(get.peptide.stats)
#read the amino acid data
gene.sequences = read.csv('C:\\Users\\dell7\\Documents\\Tzachi\\analysis\\R_project\\peptide amino acid sequences.csv')

#get the data, change it to charachters, put into pep stats function and transpose
peptide.properties <- t(get.peptide.stats(as.character(gene.sequences$pep.sequence)))

#bind the sequence data to the pep stats
gene.sequences <- cbind(gene.sequences,peptide.properties)

#get the 2nd amino acid from each peptide
gene.sequences <- gene.sequences %>%
  mutate(pep.2nd.aa = substr(pep.sequence,2,2))
#read fitseq data and left join with peptide data
fitseq.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\goodman_salis_tuller_fitseq_2_mismatch_with_0.csv'
fitseq.data = read.csv(fitseq.data.location)
fitseq.data.with.pep.data <-  left_join(fitseq.data,gene.sequences,by = 'Gene')
#transform all the data to the correct type (character, int, or number)
fitseq.data.with.pep.data <-  transform(fitseq.data.with.pep.data, Prot = as.numeric(as.character(Prot))
                       ,Bin.Pct.1 = as.numeric(as.character((Bin.Pct.1))),
                       Count.RNA = as.numeric(as.character((Count.RNA))),
                       pep.sequence = as.character(pep.sequence),
                       pep.pi=as.numeric(pep.pi),
                       pep.cost=as.numeric(pep.cost),
                       pep.demand=as.numeric(pep.demand),
                       pep.mw=as.numeric(pep.mw),
                       pep.aindex=as.numeric(pep.aindex),
                       pep.boman=as.numeric(pep.boman),
                       pep.charge=as.numeric(pep.charge),
                       pep.hmoment=as.numeric(pep.hmoment),
                       pep.hydro=as.numeric(pep.hydro),
                       pep.instability=as.numeric(pep.instability),
                       pep.Ala=as.numeric(pep.Ala),
                       pep.Cys=as.numeric(pep.Cys),
                       pep.Asp=as.numeric(pep.Asp),
                       pep.Glu=as.numeric(pep.Glu),
                       pep.Phe=as.numeric(pep.Phe),
                       pep.Gly=as.numeric(pep.Gly),
                       pep.His=as.numeric(pep.His),
                       pep.Ile=as.numeric(pep.Ile),
                       pep.Lys=as.numeric(pep.Lys),
                       pep.Leu=as.numeric(pep.Leu),
                       pep.Met=as.numeric(pep.Met),
                       pep.Asn=as.numeric(pep.Asn),
                       pep.Pro=as.numeric(pep.Pro),
                       pep.Gln=as.numeric(pep.Gln),
                       pep.Arg=as.numeric(pep.Arg),
                       pep.Ser=as.numeric(pep.Ser),
                       pep.Thr=as.numeric(pep.Thr),
                       pep.Val=as.numeric(pep.Val),
                       pep.Trp=as.numeric(pep.Trp),
                       pep.Tyr=as.numeric(pep.Tyr),
                       pep.tiny=as.numeric(pep.tiny),
                       pep.small=as.numeric(pep.small),
                       pep.aromatic=as.numeric(pep.aromatic),
                       pep.polar=as.numeric(pep.polar),
                       pep.charged=as.numeric(pep.charged),
                       pep.basic=as.numeric(pep.basic),
                       pep.acidic=as.numeric(pep.acidic),
                       pep.aliphatic=as.numeric(pep.aliphatic),
                       pep.non.polar=as.numeric(pep.non.polar))
#rewrite fitseq data
write.csv(fitseq.data.with.pep.data,
          file = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\pep_data_goodman_salis_tuller_fitseq_2_mismatch_with_0.csv')
  




aggregation <- read.csv('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\aa_aggregation_stats.csv', row.names = 1)



get_agg_data <- function(pep){
  
  pep_df = data.frame(unlist(strsplit(pep,'')))
  names(pep_df) <- 'aa'
  
  pep_df <- pep_df %>% mutate(levy = aggregation[aa,'levy_agg_prop'],aggrescan = aggregation[aa,'aggrescan'],
                              foldamyloid = aggregation[aa,'foldamyloid'],ww_hydrophobicity = aggregation[aa,'ww_hydrophobicity'])
  pep_sum <- pep_df %>% summarise(levy_mean = mean(levy),levy_median = median(levy),levy_stdev = sd(levy),levy_sum = sum(levy),
                                  aggrescan_mean = mean(aggrescan),aggrescan_median = median(aggrescan),aggrescan_stdev = sd(aggrescan),aggrescan_sum = sum(aggrescan),
                                  foldamyloid_mean = mean(foldamyloid),foldamyloid_median = median(foldamyloid),foldamyloid_stdev = sd(foldamyloid),foldamyloid_sum = sum(foldamyloid),
                                  ww_hydrophobicity_mean = mean(ww_hydrophobicity),ww_hydrophobicity_median = median(ww_hydrophobicity),ww_hydrophobicity_stdev = sd(ww_hydrophobicity),ww_hydrophobicity_sum = sum(ww_hydrophobicity))
  
  
  
  # pep_sum$instability <- instaindex(pep)
  return(pep_sum)
}
pep <- 'MSLNFLDFEQP'
get_agg_data(pep)
pep <- 'MSLNFLDFEQP'
peps <- read.csv('C://Users//dell7//Documents//Tzachi//workspace//data//peptide_amino_acid_sequences.csv')
l <- lapply(as.character(peps$sequence),get_agg_data)
df <- data.frame(peps, do.call(rbind, l))
df <- transform(df, Gene = as.character(Gene))

fitseq_genes <- read.csv('C://Users//dell7//Documents//Tzachi//workspace//data//id_gene.csv')
fitseq_genes <- transform(fitseq_genes, Gene = as.character(Gene))

fitseq_genes_agg_data <- left_join(fitseq_genes,df, by = 'Gene' )

write.csv(fitseq_genes_agg_data,'C://Users//dell7//Documents//Tzachi//workspace//data//design_agg_data.csv', row.names = F)
