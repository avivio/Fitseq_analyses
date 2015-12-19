#upload libraries
require(dplyr)
require(tidyr)
require(ggplot2)
require(scales)

#create theme with good text and line sizes
theme_aviv <-     theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
         strip.text.x = element_text(size=16,face="bold"), strip.text.y = element_text(size=16,face="bold"), 
         plot.title = element_text(size = 25,face = "bold"),
         legend.text = element_text(size=16),legend.title = element_text(size=18,face="bold")	)

load.fitseq.data <- function(fitseq.data.location){
  #function for loading fitseq data
  #read from csv
  fitseq.data = read.csv(fitseq.data.location)
#some of the feilds have textual data in goodmans data, I transform to numeric and the text becomes na and is then filtered automatically from figures
  fitseq.data = transform(fitseq.data, Prot = as.numeric(as.character(Prot))
                          ,Bin.Pct.1 = as.numeric(as.character((Bin.Pct.1))),
                          Count.RNA = as.numeric(as.character((Count.RNA))))
#chuck away most of the data, we don't need a lot of the goodman stuff in our analysis
    fitseq.data.tidy <- fitseq.data %>% select(-Count.A.DNA,-Count.A.RNA,-Count.B.DNA,-Count.B.RNA,
                                             -Bin.1,-Bin.2,-Bin.3,-Bin.4,-Bin.5,-Bin.6,-Bin.7,
                                             -Bin.8,-Bin.9,-Bin.10,-Bin.11,-Bin.12,-RNA.A,-RNA.B,
                                             -Bin.Pct.1,-Bin.Pct.2,-Bin.Pct.3,-Bin.Pct.4,-Bin.Pct.5,
                                             -Bin.Pct.6,-Bin.Pct.7,-Bin.Pct.8,-Bin.Pct.9,-Bin.Pct.10,
                                             -Bin.Pct.11,-Bin.Pct.12,-Insuff.Prot,-Insuff.DNA,
                                             -Insuff.RNA,-Fltr.BelowRange,-Fltr.AboveRange ,
                                             -Fltr.SetGood,-full.seq,-pep.sequence,-UTR,-CDS.seq,-Promoter.seq,-RBS.seq,-Promoter,
                                             -variable.seq,-full.peptide,-preRBS,-salis.status)
  #this is the tricky part, I spread the data into wide format so I could use the lineage and day as variables, data get's very heavy on memory 
  #I needed to run on only the data from days 12 because my computer sucks, you can comment the day 12 stuff and return to run on all data points
  # fitseq.data.tidy  <- fitseq.data.tidy %>%   gather(sample, frequency, A_24:F_0)
  fitseq.data.tidy  <- fitseq.data.tidy %>%   gather(sample, frequency, D_12:F_12)
  #I kept the fitseq columns with an underscore seperator becuase r was having trouble parsing periods
  fitseq.data.tidy  <- fitseq.data.tidy %>% separate(sample, c("lineage", "day"), '_',convert = T)
  # I changed the order of the rbs factor so it's presented nicely in figures
    fitseq.data.tidy$RBS.Display <- factor(fitseq.data.tidy$RBS.Display,
                                         levels = c("Strong", "Mid", "Weak", "WT"))
  
  
  #calculate the fitness by comparing day 12 (or any day) to ancestor
  fitseq.data.tidy <- fitseq.data.tidy %>% 
    group_by(day,lineage) %>% 
    mutate(sum.sample= sum(frequency),norm.freq = frequency/sum.sample,
           sum.anc= sum(anc_1),norm.anc = anc_1/sum.anc,freq.norm.anc = norm.freq/norm.anc,
           sum.sample.1= sum(frequency+1),norm.freq.1 = (frequency+1)/sum.sample.1,
           sum.anc.1= sum(anc_1+1),norm.anc.1 = (anc_1+1)/sum.anc.1,freq.norm.anc.1 = norm.freq.1/norm.anc.1)
  
  
  return(fitseq.data.tidy)
}



#load data
#use the main dateset wherever it is on your setup, I put one copy in the home directory of all the scripts
fitesq.location <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_data_with_meta_day_12_no_umi_05_08.csv'

fitseq_data_residuals <- load.fitseq.data(fitesq.location)

#throw away all data except 1 time point since we're looking at goodman data
goodman.data.dna <- fitseq_data_residuals %>% filter(day == 12, lineage == 'A') %>% select(Name, new.name, Promoter.Display, RBS.Display, Count.DNA, Prot, frequency)
#throw away wild type since it has erratic behavior and isn't well defined                                    
goodman.data.dna.no.wt <-goodman.data.dna %>% filter(RBS.Display != 'WT') 

#reorder factor
goodman.data.dna.no.wt$RBS.Display <- factor(goodman.data.dna.no.wt$RBS.Display,
                                             levels = c("Strong", "Mid", "Weak"))
goodman.data.dna.no.wt$Promoter.Display <- factor(goodman.data.dna.no.wt$Promoter.Display,
                                             levels = c("High", "Low"))
#start plotting
png('goodman_dna_vs_promoter_rbs.png',
    type="cairo",    units="in", width=10, height=6, pointsize=12, res=500)    
#using the interaction of rbs and promter creates a factor of all combinations, we show the DNA column + 1 (for 0 data) in log2
ggplot(goodman.data.dna.no.wt, aes(x =interaction(RBS.Display,Promoter.Display), y = log2(Count.DNA+1) )) +
#the rest is astetics, colors, axis labels, and legible tick lables
    geom_boxplot(fill = hue_pal()(2)[2]) +
  ylab('Log 2 DNA coverage') +
  xlab('Promoter and RBS type') +
  scale_x_discrete(labels=c("High Strong", "High Mid", "High Weak","Low Strong", "Low Mid", "Low Weak")) +
  theme_aviv
dev.off()

#this is for doing the wilcox test to get the p value of difference between interesting groups
goodman.data.dna.no.wt.1 <- goodman.data.dna.no.wt %>% filter(Promoter.Display == 'Low',RBS.Display == 'Mid')
goodman.data.dna.no.wt.2 <- goodman.data.dna.no.wt %>% filter(Promoter.Display == 'Low',RBS.Display == 'Weak')

wilcox.test(goodman.data.dna.no.wt.1$Count.DNA,goodman.data.dna.no.wt.2$Count.DNA,alternative = 'less')

# protein level by promoter rbs, we didn't use this figure, log2 prot level instead of dna, as if the clue that led to our experiment
png('goodman_prot_vs_promoter_rbs.png',
    type="cairo",    units="in", width=10, height=6, pointsize=12, res=500) 
ggplot(goodman.data.dna.no.wt, aes(x =interaction(RBS.Display,Promoter.Display), y = log2(Prot + 1) )) +
  geom_boxplot(fill = hue_pal()(2)[2]) +
  ylab('Log 2 expression level') +
  xlab('Promoter and RBS type') +
  scale_x_discrete(labels=c("High Strong", "High Mid", "High Weak","Low Strong", "Low Mid", "Low Weak")) +
  theme_aviv
dev.off()
# differences in dna coverage, fill, showing dna coverage by rbs type
png('dna_coverage_rbs_fill.png',
    type="cairo",    units="in", width=8, height=8, pointsize=12, res=500)    
ggplot(goodman.data.dna.no.wt, aes(log10(Count.DNA), fill = RBS.Display)) +
  stat_bin(geom='area',position = 'identity',alpha=.5) +
  theme_aviv + 
  ylab('# Designs') +
  xlab('Log 10 DNA coverage') +
  guides(fill = guide_legend("RBS"))
dev.off()


# differences in dna coverage, fill, showing dna coverage by rbs type in facet
png('dna_coverage_rbs_facet.png',
    type="cairo",    units="in", width=8, height=8, pointsize=12, res=500)    
ggplot(goodman.data.dna.no.wt, aes(log10(Count.DNA))) +
  stat_bin(geom='area',position = 'identity', fill = hue_pal()(2)[2]) +
  theme_aviv + 
  ylab('# Designs') +
  xlab('Log 10 DNA coverage') +
  facet_grid(RBS.Display~.)
dev.off()