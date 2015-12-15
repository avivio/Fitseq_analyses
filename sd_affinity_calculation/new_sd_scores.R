# install all packages
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
 install.packages('psych')
# install.packages('dplyr')
# install.packages('ggplot')
# install.packages('tidyr')

#upload packages
require(Biostrings)
require(psych)
require(dplyr)
require(ggplot2)
require(tidyr)

#set constant variables
t = 303.15
R = 8.3144



get.sd.score <- function(seq,position){
  #get the sd score of a sequence at a certain position start from the highest possible positive value
  sd.score <- Inf
  # look at all the positions that are between 8 and 11 bases backward
  for (i in 8:11){
    base <- position-i
    #the site is defined as the position (between 8 and 11 upstream) the 7 bases behind it and 2 in front
    site <- substring(seq,base-7,base+2)
    #get the delta g from the tenmer table
    sd <- tenmers[site,1]
    #if the score is lower than the existing score, update it
    if (sd < sd.score){
      sd.score <- sd
      
    }
  }
  #calculate the Kb of the site by plugging it in to the equation Kb= e^(deltaG/(R*T))
  k.b = exp((sd.score/(R*t)))
  #the velocity is derived from Kb/(Kb+1) and returned
  velocity = k.b/(1+k.b)
  return(velocity)
}



factor.2.rna.char <- function(fact){
  #a technical function that turns a factor value to an RNA string for one of the packages
  return( toupper( as.character(RNAString(DNAString(as.character(fact)))))) 
}
factor.2.rna.char <- Vectorize(factor.2.rna.char)


get.sd.propterties <- function(utr,cds){
  # gets the SD affinity velocity for a coding sequence and an untranslated region
  #get the start location, create the sequence variable, and a vector for the values
  start.location <- nchar(utr)+1
  seq <- paste0(utr,cds)
  sd.vec <- vector(length = nchar(cds), mode = 'double')
  
  index <-  0
  #go over the sequnence and get the sd velocity for each position and return the resulting velocity vector
  for (i in start.location:(nchar(seq))){
    index <-  index +1
    sd <- get.sd.score(seq,i)
    sd.vec[index] <- sd
    
  }
  return(sd.vec)


}
get.sd.propterties <- Vectorize(get.sd.propterties)


#load the tenmers file into memeore
print('load tenmers')


tenmers  <-   read.csv('/Users/aviv.r/Documents/Tzachi/sd_tenmer_data.csv',row.names = 1)

#load the sequence data into memory

print('get sd data for fitseq sequences')
fitseq.data <- read.csv('/Users/aviv.r/Documents/Tzachi/fitseq_cds_utr.csv')

# fitseq.data.test <- fitseq.data[1:5,]
#get the sd properties for the all the sequences
sd.properties <- get.sd.propterties(factor.2.rna.char(fitseq.data$UTR),factor.2.rna.char(fitseq.data$CDS.seq))

#transpose the data because vectorize returns weird values
print('transpose')
sd.properties <- data.frame(t(sd.properties))
print('numeric')

#set the names of the columns of each position
names(sd.properties) <- c('sd_pos_1','sd_pos_2','sd_pos_3','sd_pos_4','sd_pos_5','sd_pos_6','sd_pos_7','sd_pos_8','sd_pos_9','sd_pos_10',
                          'sd_pos_11','sd_pos_12','sd_pos_13','sd_pos_14','sd_pos_15','sd_pos_16','sd_pos_17','sd_pos_18','sd_pos_19','sd_pos_20',
                          'sd_pos_21','sd_pos_22','sd_pos_23','sd_pos_24','sd_pos_25','sd_pos_26','sd_pos_27','sd_pos_28','sd_pos_29','sd_pos_30',
                          'sd_pos_31','sd_pos_32','sd_pos_33')

#summarize the data into means, max, min, count, sum, median
sd.properties.summarized   <- sd.properties %>% 
  rowwise() %>% 
  do(sd_harm_mean = harmonic.mean(as.numeric(.)), sd_art_mean = mean(as.numeric(.)),
     sd_max = max(as.numeric(.)),sd_min = min(as.numeric(.)),
     sd_sum = sum(as.numeric(.)), sd_median = median(as.numeric(.)),sd_count = sum(. >  0.5 ))
     
#convert to numeric
sd.properties.summarized <- transform(sd.properties.summarized,sd_harm_mean = as.numeric(sd_harm_mean), sd_max = as.numeric(sd_max),
                                      sd_min = as.numeric(sd_min), sd_sum = as.numeric(sd_sum),sd_median = as.numeric(sd_median),
                                      sd_art_mean = as.numeric(sd_art_mean), sd_count = as.numeric(sd_count))
#bind to the raw data
sd.properties.name <- bind_cols(fitseq.data %>% select(Name, new.name),sd.properties.summarized, sd.properties)

#print to csv file
print('output data to file')
write.csv(sd.properties.name, row.names=FALSE,
          file = '/Users/aviv.r/Documents/Tzachi/sd_affinity_per_position_new_calc.csv')
print('Done')
