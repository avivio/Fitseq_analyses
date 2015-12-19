#load libraries
require(tidyr)
require(ggplot2)
require(dplyr)

theme_aviv <-     theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
         strip.text.x = element_text(size=16,face="bold"), strip.text.y = element_text(size=16,face="bold"), 
         plot.title = element_text(size = 25,face = "bold"),
         legend.text = element_text(size=16),legend.title = element_text(size=18,face="bold")	)



#load data
#can't upload data to git so I need to upload this data to the right place later
fitseq.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_data_with_meta_no_umi_generations.csv'
fitseq.data = read.csv(fitseq.data.location)

#throw out uneeded goodman data
fitseq.data.tidy <- fitseq.data %>% select(-Count.A.DNA,-Count.A.RNA,-Count.B.DNA,-Count.B.RNA,
                                           -Bin.1,-Bin.2,-Bin.3,-Bin.4,-Bin.5,-Bin.6,-Bin.7,
                                           -Bin.8,-Bin.9,-Bin.10,-Bin.11,-Bin.12,-RNA.A,-RNA.B,
                                           -Bin.Pct.1,-Bin.Pct.2,-Bin.Pct.3,-Bin.Pct.4,-Bin.Pct.5,
                                           -Bin.Pct.6,-Bin.Pct.7,-Bin.Pct.8,-Bin.Pct.9,-Bin.Pct.10,
                                           -Bin.Pct.11,-Bin.Pct.12,-Insuff.Prot,-Insuff.DNA,
                                           -Insuff.RNA,-Fltr.BelowRange,-Fltr.AboveRange ,
                                           -Fltr.SetGood,-full.seq)
#convert from long to wide on sample column, plit sample to generation and lineage
fitseq.data.tidy  <- fitseq.data.tidy %>%   gather(sample, frequency, A_168:F_0)
fitseq.data.tidy  <- fitseq.data.tidy %>% separate(sample, c("lineage", "generation"), '_',convert = T)
#refactor the rbs and the lineages so they look nice in the figure
fitseq.data.tidy$RBS.Display <- factor(fitseq.data.tidy$RBS.Display,
                                       levels = c("Strong", "Mid", "Weak", "WT"))
fitseq.data.tidy$lineage <- factor(fitseq.data.tidy$lineage,
                                   levels = c("Ancestor", "A", "B", "C", "D", "E", "F"))

#create the variable of the fitnss by normalizing the frequency to the sample and then the sample to anc_1 (+1 for log transformation)
fitseq.data.tidy  <- fitseq.data.tidy %>% 
  group_by(generation,lineage) %>%
  mutate(sum.sample= sum(frequency),norm.freq = frequency/sum.sample,
         sum.anc= sum(anc_1),norm.anc = anc_1/sum.anc,freq.norm.anc = norm.freq/norm.anc,
         sum.sample.1= sum(frequency+1),norm.freq.1 = (frequency+1)/sum.sample.1,
         sum.anc.1= sum(anc_1+1),norm.anc.1 = (anc_1+1)/sum.anc.1,freq.norm.anc.1 = norm.freq.1/norm.anc.1)


#historgrams of normalized frequency + 1 log transformed for different day, rbs-promoter groups, and lineages,

line.day.histogram   <- function(lineage_letter,days,fitseq.data.tidy){
  #a funtion the plot the density plots of each day and facet by the rbs and promoter
#start plot, change the label of the figure by the lineage
    png(paste0('days_fill_normalized_frequency_histograms_lineage_',lineage_letter,'_hist.png'),
      type="cairo",    units="in", width=10, height=6, pointsize=12, res=500)    
  #plot the density plots using all the data in the dataset that is in the days vector, is the right lineage and not wt rbs
  #x log normalized frequency, and fill is the day, facet by rbs and promtoere
    p =  ggplot(filter(fitseq.data.tidy,day %in% days,lineage == lineage_letter,RBS.Display != 'WT') ,
              aes(log2(freq.norm.anc.1),fill = factor(day)))  +
    stat_bin(geom='area',position = 'identity',alpha=.5) +
    theme_minimal() + 
    facet_grid(RBS.Display~Promoter.Display)  +
    theme_aviv +
    guides(color = guide_legend("Day"),fill = guide_legend("Day")) +
    xlab('\nLog 2 of normalized frequency + 1') +
    ylab('Count\n') +
      #I expanded the limit by looking at what looks best for all lineages
    expand_limits(y=900)
  
  print(p)
  dev.off()
  
}


line.rbs.histogram   <- function(lineage_letter,generations,fitseq.data.tidy){
  #a funtion the plot the density plots of each rbs and facet by the day and promoter
  #set the levels of the promoter so you can understand the meaning in the facet text
  levels(fitseq.data.tidy$Promoter.Display) <- c("High promoter", "Low promoter")
  #create the figures as png with each figure a different lineage
  png(paste0('rbs_fill_frequency_histograms_lineage_',lineage_letter,'.png'),
      type="cairo",    units="in", width=9, height=11, pointsize=12, res=500)    
  #plot the log normalized freq +1 with the fill being the rbs, using the data that fits the lineage, days vector and no wt rbs
  p =  ggplot(filter(fitseq.data.tidy,generation %in% generations,lineage == lineage_letter,RBS.Display != 'WT') ,
              aes(log2(norm.freq.1),fill = factor(RBS.Display)))  +
    stat_bin(geom='area',position = 'identity',alpha=.5) +
    theme_minimal() + 
    #facet by generation and promoter
    facet_grid(generation~Promoter.Display)  +
    theme_aviv +
    guides(color = guide_legend("Generation"),fill = guide_legend("RBS")) +
    xlab('\nLog 2 frequency') +
    ylab('# designs\n') +
    # I played with the limits by eye so that all the lineages look good
    expand_limits(y=1000) +
    coord_cartesian(xlim = c(-24,-10))
  
  print(p)
  dev.off()
  
}

#go over each lineage and ru figures filled by day and faceted by rbs
for (lineage in c('A','B','C','D','E','F')){
  line.day.histogram(lineage,c(4,12,28),fitseq.data.tidy)
  
}

#go over each lineage and ru figures filled by rbs and faceted by generation
for (lineage in c('A','B','C','D','E','F')){
  line.rbs.histogram(lineage,c(28,56,84,112,140,168,196),fitseq.data.tidy)
  
}
