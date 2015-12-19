#load packages
require(dplyr)
require(tidyr)
require(ggplot2)
require(ineq)
require(scales)

#set publication level theme
theme_aviv <-     theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
         strip.text.x = element_text(size=16,face="bold"), strip.text.y = element_text(size=16,face="bold"), 
         plot.title = element_text(size = 25,face = "bold"))

#load data
pre.umi.data.location  <- 'match_count_2_mismatches_generations.csv'
pre.umi.data = read.csv(pre.umi.data.location)
#gather by sample so the data is in wide format, each sample a frequency
pre.umi.data  <- pre.umi.data %>%   gather(sample, frequency, A_168:F_28)
#split the sample name into lineage and generation
pre.umi.data  <- pre.umi.data %>% separate(sample, c("lineage", "generation"), '_',convert = T)
#create a new variable that is the relative frequency of design in the whole sample
pre.umi.data <- pre.umi.data %>%
  group_by(generation,lineage) %>% 
  mutate(sum.sample= sum(frequency),norm.freq = frequency/sum.sample)

#create the gini score for each sample using the ineq package
gini <- pre.umi.data %>% 
  group_by(generation,lineage) %>%
  summarise(
            gini =  ineq(frequency,type="Gini")) 

#chuck the ancestor lineage since it's not nice in the figure, refactor after
gini <- filter(gini,lineage !='Ancestor')
gini$lineage = factor(gini$lineage)


#plot the figure using a barplot with the generation as the x and gini as y, facet on lineage 
png('gini_over_time_lienages.png',
    type="cairo",    units="in", width=22, height=6, pointsize=12, res=500)    
ggplot(gini, aes(y = gini,x = factor(generation))) + 
  geom_bar(stat = 'identity',position = 'dodge', fill = hue_pal()(2)[2]) +
   ylab('Gini index\n') + 
  xlab('\nGeneration') + 
  theme_aviv +
  facet_wrap(~lineage,nrow = 1)
  
dev.off()