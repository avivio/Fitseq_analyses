
require(seqinr)
load_packages()
data(aacost)
names(aacost)[1]  <- 'amino_acid'
aademand  <- read.csv('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\amino_acid_supply_demand.csv')

cost_demand_cor = paste0('\nCorrelation: ',signif( cor(cost_demand$demand,cost_demand$cost, method = 'pearson'),2), '    \np-value: ',
                         signif(cor.test(cost_demand$demand,cost_demand$cost, method = 'pearson')$p.value,4), '    ')
ggplot(cost_demand,aes(x = cost, y = demand)) +
  geom_point(size = 5) +
  geom_text(aes(label = paste(' ', amino_acid)),hjust=0, vjust=0) + geom_smooth(method = 'lm') + 
  theme_aviv +
  annotate ('text', x = Inf, y = Inf,label = cost_demand_cor, hjust = 1, vjust = 1, face = 'bold', color = 'blue'  )


#above n pos.neg


#getting data all
pos_neg_amino_acids <- fitseq_data_residuals_above_14 %>% 
  select(day,Gene,n_pos_neg,pep_Ala:pep_Tyr) %>%
  group_by(day,Gene,n_pos_neg) %>%
  mutate_each(funs(.*11))%>%
  gather(amino_acid,frequency,pep_Ala:pep_Tyr) %>%
  separate(amino_acid,into = c('pep','amino_acid'),'[_]',extra = 'drop') %>%
  select(-pep) %>%
  group_by(day,n_pos_neg,amino_acid) %>%
  summarize(gene_frequency = n(),amino_acid_sum = sum(frequency),
            amino_acid_frequency = amino_acid_sum/(gene_frequency*11)) %>% 
  select(-gene_frequency,-amino_acid_sum) %>% 
  spread(n_pos_neg,amino_acid_frequency) %>%
  left_join(aacost %>% select(amino_acid,tot),by = 'amino_acid') %>% 
  mutate(ratio = Positive/Negative) %>% 
  left_join(aademand ,by = 'amino_acid')
names(pos_neg_amino_acids)[5]  <- 'cost'


pos_neg_amino_acids_rbs <- fitseq_data_residuals_above_14 %>%  
  select(day,Gene,RBS_Display,n_pos_neg,pep_Ala:pep_Tyr) %>%
  group_by(day,Gene,RBS_Display,n_pos_neg) %>%
  mutate_each(funs(.*11),-n_pos_neg)%>%
  gather(amino_acid,frequency,pep_Ala:pep_Tyr) %>%
  separate(amino_acid,into = c('pep','amino_acid'),'[_]',extra = 'drop') %>%
  select(-pep) %>%
  group_by(day,RBS_Display,n_pos_neg,amino_acid) %>%
  summarize(gene_frequency = n(),amino_acid_sum = sum(frequency),
            amino_acid_frequency = amino_acid_sum/(gene_frequency*11)) %>% 
  select(-gene_frequency,-amino_acid_sum) %>% 
  spread(n_pos_neg,amino_acid_frequency) %>%
  filter( RBS_Display != 'WT') %>%
  left_join(aacost %>% select(amino_acid,tot),by = 'amino_acid') %>% 
  mutate(ratio = Positive/Negative) %>% 
  left_join(aademand ,by = 'amino_acid')

names(pos_neg_amino_acids_rbs)[6]  <- 'cost'


#cost plots
pos_neg_lim <- max(c(max(pos_neg_amino_acids$Positive),max(pos_neg_amino_acids$Negative)))


png('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\for_thesis\\amino_acid_pos_neg_cost_col.png',
    type="cairo",    units="in", width=10, height=6, pointsize=12, res=500)  

ggplot(pos_neg_amino_acids %>% filter(day == 12),
       aes(y=Positive,x=Negative,color = cost)) +
  geom_point(size = 5) + 
  scale_color_gradient2(low = hue_pal()(3)[3],
                        high = hue_pal()(3)[1],mid = hue_pal()(3)[2],midpoint =45,name = 'Cost') +
  theme_aviv + 
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label=amino_acid),hjust=0, vjust=0) + 
  coord_equal(xlim =c(0,pos_neg_lim + 0.01),ylim =c(0,pos_neg_lim + 0.01) ) 

dev.off()

cost_cor <- paste0('\nCorrelation: ',signif( cor(pos_neg_amino_acids$ratio,pos_neg_amino_acids$cost, method = 'pearson'),2), '    \np-value: ',
                   signif(cor.test(pos_neg_amino_acids$ratio,pos_neg_amino_acids$cost, method = 'pearson')$p.value,4), '    ')


png('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\for_thesis\\amino_acid_ratio_vs_cost.png',
    type="cairo",    units="in", width=14, height=10, pointsize=12, res=500)  

ggplot(pos_neg_amino_acids %>% filter(day == 12),
       aes(y=ratio,x=cost)) +
  geom_point(size = 5) + 
  theme_aviv + 
  geom_text(aes(label=paste(' ',amino_acid)),hjust=0, vjust=0) +
  geom_smooth(method = 'lm') +
  ylab('Ration of Positive amino acid frequency over Negative') 
  # annotate ('text', x = Inf, y = Inf,label = cost_cor, hjust = 1, vjust = 1, face = 'bold', color = 'blue'  )
print(cost_cor)
dev.off()

pos_neg_amino_acids_rbs$RBS_Display <- factor(pos_neg_amino_acids_rbs$RBS_Display,
                                            levels = c('Strong','Mid','Weak','WT'))

pos_neg_lim_rbs <- max(c(max(pos_neg_amino_acids_rbs$Positive),max(pos_neg_amino_acids_rbs$Negative)))

ggplot(pos_neg_amino_acids_rbs %>% filter(day == 12),
       aes(y=Positive,x=Negative,color = cost)) +
  geom_point() +
  scale_color_gradient2(low = hue_pal()(3)[3],
                        high = hue_pal()(3)[1],mid = hue_pal()(3)[2],midpoint =45)+
  theme_aviv + 
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label=amino_acid),hjust=0, vjust=0) + 
  coord_equal(xlim =c(0,pos_neg_lim_rbs),ylim =c(0,pos_neg_lim_rbs) ) + 
  facet_grid(~RBS_Display)



cost_cor_rbs <- pos_neg_amino_acids_rbs %>% 
  ungroup() %>% 
  group_by(RBS_Display) %>% 
  summarize(cor = cor(ratio,cost, method = 'pearson'),p_value =  cor.test(ratio,cost, method = 'pearson')$p.value) %>% 
  mutate(cor_string = paste0('\nCorrelation: ',signif(cor,2),  '    \np-value: ',signif(p_value,4), '    ')) %>% 
  select(RBS_Display,cor_string)

cost_cor_strings_rbs <- cost_cor_rbs$cor_string
names(cost_cor_strings_rbs) <- cost_cor_rbs$RBS_Display

png('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\for_thesis\\amino_acid_ratio_vs_cost_rbs.png',
    type="cairo",    units="in", width=25, height=10, pointsize=12, res=500)  

ggplot(pos_neg_amino_acids_rbs %>% filter(day == 12),
       aes(y= ratio,x=cost)) +
  geom_point(size = 5) + 
  theme_aviv + 
  geom_text(aes(label=paste(' ',amino_acid)),hjust=0, vjust=0) +
  geom_smooth(method = 'lm')  +
  ylab('Ration of Positive amino acid frequency over Negative') +
  facet_grid(~RBS_Display) 
  # geom_text (aes(label = cost_cor_strings_rbs[RBS_Display]), x = Inf, y = Inf, hjust = 1, vjust = 1, face = 'bold', color = 'blue'  )
dev.off()

#demand plots
pos_neg_lim <- max(c(max(pos_neg_amino_acids$Positive),max(pos_neg_amino_acids$Negative)))

ggplot(pos_neg_amino_acids %>% filter(day == 12),
       aes(y=Positive,x=Negative,color = demand)) +
  geom_point(size = 5) + 
  scale_color_gradient2(low = hue_pal()(3)[3],
                        high = hue_pal()(3)[1],mid = hue_pal()(3)[2],midpoint =1.5e+6,
                        name='Amino acid\ndemand') +
  theme_aviv + 
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label=amino_acid),hjust=0, vjust=0) + 
  coord_equal(xlim =c(0,pos_neg_lim + 0.01),ylim =c(0,pos_neg_lim + 0.01) ) 
  
demand_cor <- paste0('\nCorrelation: ',signif( cor(pos_neg_amino_acids$ratio,pos_neg_amino_acids$demand, method = 'pearson'),2), '    \np-value: ',
                   signif(cor.test(pos_neg_amino_acids$ratio,pos_neg_amino_acids$demand, method = 'pearson')$p.value,4), '    ')

ggplot(pos_neg_amino_acids %>% filter(day == 12),
       aes(y=ratio,x=demand)) +
  geom_point(size = 5) + 
  theme_aviv + 
  geom_text(aes(label=paste(' ',amino_acid)),hjust=0, vjust=0) +
  geom_smooth(method = 'lm')  +
  ylab('Ration of Positive amino acid frequency over Negative') +
  xlab(' acid demand') +
annotate ('text', x = Inf, y = Inf,label = demand_cor, hjust = 1, vjust = 1, face = 'bold', color = 'blue'  )



pos_neg_lim_rbs <- max(c(max(pos_neg_amino_acids_rbs$Positive),max(pos_neg_amino_acids_rbs$Negative)))

ggplot(pos_neg_amino_acids_rbs %>% filter(day == 12),
       aes(y=Positive,x=Negative,color = demand)) +
  geom_point() +
  scale_color_gradient2(low = hue_pal()(3)[3],
                        high = hue_pal()(3)[1],mid = hue_pal()(3)[2],midpoint =1.5e+6,
                        name='Amino acid\ndemand')+
  theme_aviv + 
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label=amino_acid),hjust=0, vjust=0) + 
  coord_equal(xlim =c(0,pos_neg_lim_rbs),ylim =c(0,pos_neg_lim_rbs) ) + 
  facet_grid(~RBS_Display)




demand_cor_rbs <- pos_neg_amino_acids_rbs %>% 
  ungroup() %>% 
  group_by(RBS_Display) %>% 
  summarize(cor = cor(ratio,demand, method = 'pearson'),p_value =  cor.test(ratio,demand, method = 'pearson')$p.value) %>% 
  mutate(cor_string = paste0('\nCorrelation: ',signif(cor,2),  '    \np-value: ',signif(p_value,4), '    ')) %>% 
  select(RBS_Display,cor_string)


demand_cor_strings_rbs <- demand_cor_rbs$cor_string
names(demand_cor_strings_rbs) <- demand_cor_rbs$RBS_Display

ggplot(pos_neg_amino_acids_rbs %>% filter(day == 12),
       aes(y=ratio,x=demand)) +
  geom_point(size = 5) + 
  theme_aviv + 
  geom_text(aes(label=paste(' ',amino_acid)),hjust=0, vjust=0) +
  geom_smooth(method = 'lm')  +
  ylab('Ration of Positive amino acid frequency over Negative') +
  xlab('Amino acid demand') +
  facet_grid(~RBS_Display) +
  geom_text (aes(label = demand_cor_strings_rbs[RBS_Display]), x = Inf, y = Inf, hjust = 1, vjust = 1, face = 'bold', color = 'blue'  )


#demand/supply plots
pos_neg_lim <- max(c(max(pos_neg_amino_acids$Positive),max(pos_neg_amino_acids$Negative)))

ggplot(pos_neg_amino_acids %>% filter(day == 12),
       aes(y=Positive,x=Negative,color = demand_tai)) +
  geom_point(size = 5) + 
  scale_color_gradient2(low = hue_pal()(3)[3],
                        high = hue_pal()(3)[1],mid = hue_pal()(3)[2],midpoint =7e+6,
                        name='Amino acid\ndemand over tAI') +
  theme_aviv + 
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label=amino_acid),hjust=0, vjust=0) + 
  coord_equal(xlim =c(0,pos_neg_lim + 0.01),ylim =c(0,pos_neg_lim + 0.01) ) 

demand_tai_cor <- paste0('\nCorrelation: ',signif( cor(pos_neg_amino_acids$ratio,pos_neg_amino_acids$demand_tai, method = 'pearson'),2), '    \np-value: ',
                     signif(cor.test(pos_neg_amino_acids$ratio,pos_neg_amino_acids$demand_tai, method = 'pearson')$p.value,4), '    ')

ggplot(pos_neg_amino_acids %>% filter(day == 12),
       aes(y=ratio,x=demand_tai)) +
  geom_point(size = 5) + 
  theme_aviv + 
  geom_text(aes(label=paste(' ',amino_acid)),hjust=0, vjust=0) +
  geom_smooth(method = 'lm')  +
  ylab('Ration of Positive amino acid frequency over Negative') +
  xlab('Amnino acid demand over tAI') +
  annotate ('text', x = Inf, y = Inf,label = demand_tai_cor, hjust = 1, vjust = 1, face = 'bold', color = 'blue'  )



pos_neg_lim_rbs <- max(c(max(pos_neg_amino_acids_rbs$Positive),max(pos_neg_amino_acids_rbs$Negative)))

ggplot(pos_neg_amino_acids_rbs %>% filter(day == 12),
       aes(y=Positive,x=Negative,color = demand_tai)) +
  geom_point() +
  scale_color_gradient2(low = hue_pal()(3)[3],
                        high = hue_pal()(3)[1],mid = hue_pal()(3)[2],midpoint =7e+6,
                        name='Amino acid\ndemand over tAI')+
  theme_aviv + 
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label=amino_acid),hjust=0, vjust=0) + 
  coord_equal(xlim =c(0,pos_neg_lim_rbs),ylim =c(0,pos_neg_lim_rbs) ) + 
  facet_grid(~RBS_Display)




demand_tai_cor_rbs <- pos_neg_amino_acids_rbs %>% 
  ungroup() %>% 
  group_by(RBS_Display) %>% 
  summarize(cor = cor(ratio,demand_tai, method = 'pearson'),p_value =  cor.test(ratio,demand_tai, method = 'pearson')$p.value) %>% 
  mutate(cor_string = paste0('\nCorrelation: ',signif(cor,2),  '    \np-value: ',signif(p_value,4), '    ')) %>% 
  select(RBS_Display,cor_string)


demand_tai_cor_strings_rbs <- demand_tai_cor_rbs$cor_string
names(demand_tai_cor_strings_rbs) <- demand_tai_cor_rbs$RBS_Display

ggplot(pos_neg_amino_acids_rbs %>% filter(day == 12),
       aes(y=ratio,x=demand_tai)) +
  geom_point(size = 5) + 
  theme_aviv + 
  geom_text(aes(label=paste(' ',amino_acid)),hjust=0, vjust=0) +
  geom_smooth(method = 'lm')  +
  ylab('Ration of Positive amino acid frequency over Negative') +
  xlab('Amino acid demand over tAI') +
  facet_grid(~RBS_Display) +
  geom_text (aes(label = demand_tai_cor_strings_rbs[RBS_Display]), x = Inf, y = Inf, hjust = 1, vjust = 1, face = 'bold', color = 'blue'  )




#each lineage seperately


#getting data all
pos_neg_amino_acids <- fitseq_data_residuals_above_14 %>% 
  select(day,lineage,Gene,pos_neg,pep_Ala:pep_Tyr) %>%
  group_by(day,lineage,Gene,pos_neg) %>%
  mutate_each(funs(.*11))%>%
  gather(amino_acid,frequency,pep_Ala:pep_Tyr) %>%
  separate(amino_acid,into = c('pep','amino_acid'),'[_]',extra = 'drop') %>%
  select(-pep) %>%
  group_by(day,lineage,pos_neg,amino_acid) %>%
  summarize(gene_frequency = n(),amino_acid_sum = sum(frequency),
            amino_acid_frequency = amino_acid_sum/(gene_frequency*11)) %>% 
  select(-gene_frequency,-amino_acid_sum) %>% 
  spread(pos_neg,amino_acid_frequency) %>%
  left_join(aacost %>% select(amino_acid,tot),by = 'amino_acid') %>% 
  mutate(ratio = Positive/Negative) %>% 
  left_join(aademand ,by = 'amino_acid')
names(pos_neg_amino_acids)[6]  <- 'cost'


pos_neg_amino_acids_rbs <- fitseq_data_residuals_above_14 %>%  
  select(day,lineage,Gene,RBS_Display,pos_neg,pep_Ala:pep_Tyr) %>%
  group_by(day,lineage,Gene,RBS_Display,pos_neg) %>%
  mutate_each(funs(.*11),-pos_neg)%>%
  gather(amino_acid,frequency,pep_Ala:pep_Tyr) %>%
  separate(amino_acid,into = c('pep','amino_acid'),'[_]',extra = 'drop') %>%
  select(-pep) %>%
  group_by(day,lineage,RBS_Display,pos_neg,amino_acid) %>%
  summarize(gene_frequency = n(),amino_acid_sum = sum(frequency),
            amino_acid_frequency = amino_acid_sum/(gene_frequency*11)) %>% 
  select(-gene_frequency,-amino_acid_sum) %>% 
  spread(pos_neg,amino_acid_frequency) %>%
  filter( RBS_Display != 'WT') %>%
  left_join(aacost %>% select(amino_acid,tot),by = 'amino_acid') %>% 
  mutate(ratio = Positive/Negative) %>% 
  left_join(aademand ,by = 'amino_acid')

names(pos_neg_amino_acids_rbs)[7]  <- 'cost'


#cost plots
pos_neg_lim <- max(c(max(pos_neg_amino_acids$Positive),max(pos_neg_amino_acids$Negative)))

ggplot(pos_neg_amino_acids %>% filter(day == 12),
       aes(y=Positive,x=Negative,color = cost)) +
  geom_point(size = 5) + 
  scale_color_gradient2(low = hue_pal()(3)[3],
                        high = hue_pal()(3)[1],mid = hue_pal()(3)[2],midpoint =45) +
  theme_aviv + 
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label=amino_acid),hjust=0, vjust=0) + 
  coord_equal(xlim =c(0,pos_neg_lim + 0.01),ylim =c(0,pos_neg_lim + 0.01) ) + 
  facet_grid(~lineage)


cost_cor_lineage <- pos_neg_amino_acids_rbs %>% 
  ungroup() %>% 
  group_by(lineage) %>% 
  summarize(cor = cor(ratio,cost, method = 'pearson'),p_value =  cor.test(ratio,cost, method = 'pearson')$p.value) %>% 
  mutate(cor_string = paste0('\nCorrelation: ',signif(cor,2),  '    \np-value: ',signif(p_value,4), '    ')) %>% 
  select(lineage,cor_string)
cost_cor_lineage_string <- cost_cor_lineage$cor_string
names(cost_cor_lineage_string) <- cost_cor_lineage$lineage



ggplot(pos_neg_amino_acids %>% filter(day == 12),
       aes(y=ratio,x=cost)) +
  geom_point(size = 5) + 
  theme_aviv + 
  geom_text(aes(label=paste(' ',amino_acid)),hjust=0, vjust=0) +
  geom_smooth(method = 'lm') +
  ylab('Ration of Positive amino acid frequency over Negative') +
  facet_grid(~lineage) +
  geom_text (aes(label = cost_cor_lineage_string[lineage]), x = Inf, y = Inf, hjust = 1, vjust = 1, face = 'bold', color = 'blue'  ) 

pos_neg_lim_rbs <- max(c(max(pos_neg_amino_acids_rbs$Positive),max(pos_neg_amino_acids_rbs$Negative)))

ggplot(pos_neg_amino_acids_rbs %>% filter(day == 12),
       aes(y=Positive,x=Negative,color = cost)) +
  geom_point() +
  scale_color_gradient2(low = hue_pal()(3)[3],
                        high = hue_pal()(3)[1],mid = hue_pal()(3)[2],midpoint =45)+
  theme_aviv + 
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label=amino_acid),hjust=0, vjust=0) + 
  coord_equal(xlim =c(0,pos_neg_lim_rbs),ylim =c(0,pos_neg_lim_rbs) ) + 
  facet_grid(lineage~RBS_Display)




ggplot(pos_neg_amino_acids_rbs %>% filter(day == 12),
       aes(y=ratio,x=cost)) +
  geom_point(size = 5) + 
  theme_aviv + 
  geom_text(aes(label=paste(' ',amino_acid)),hjust=0, vjust=0) +
  geom_smooth(method = 'lm') + 
  ylab('Ration of Positive amino acid frequency over Negative') +
  facet_grid(lineage~RBS_Display) 

#demand plots
pos_neg_lim <- max(c(max(pos_neg_amino_acids$Positive),max(pos_neg_amino_acids$Negative)))

ggplot(pos_neg_amino_acids %>% filter(day == 12),
       aes(y=Positive,x=Negative,color = demand)) +
  geom_point(size = 5) + 
  scale_color_gradient2(low = hue_pal()(3)[3],
                        high = hue_pal()(3)[1],mid = hue_pal()(3)[2],midpoint =1.5e+6,
                        name='Amino acid\ndemand') +
  theme_aviv + 
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label=amino_acid),hjust=0, vjust=0) + 
  coord_equal(xlim =c(0,pos_neg_lim + 0.01),ylim =c(0,pos_neg_lim + 0.01) )  + 
  facet_grid(~lineage)



demand_cor_lineage <- pos_neg_amino_acids_rbs %>% 
  ungroup() %>% 
  group_by(lineage) %>% 
  summarize(cor = cor(ratio,demand, method = 'pearson'),p_value =  cor.test(ratio,demand, method = 'pearson')$p.value) %>% 
  mutate(cor_string = paste0('\nCorrelation: ',signif(cor,2),  '    \np-value: ',signif(p_value,4), '    ')) %>% 
  select(lineage,cor_string)
demand_cor_lineage_string <- demand_cor_lineage$cor_string
names(demand_cor_lineage_string) <- demand_cor_lineage$lineage

ggplot(pos_neg_amino_acids %>% filter(day == 12),
       aes(y=ratio,x=demand)) +
  geom_point(size = 5) + 
  theme_aviv + 
  geom_text(aes(label=paste(' ',amino_acid)),hjust=0, vjust=0) +
  geom_smooth(method = 'lm') +
  facet_grid(~lineage) +
  ylab('Ration of Positive amino acid frequency over Negative') +
  xlab('Amino acid demand') +
  geom_text (aes(label = demand_cor_lineage_string[lineage]), x = Inf, y = Inf, hjust = 1, vjust = 1, face = 'bold', color = 'blue'  ) 

pos_neg_lim_rbs <- max(c(max(pos_neg_amino_acids_rbs$Positive),max(pos_neg_amino_acids_rbs$Negative)))

ggplot(pos_neg_amino_acids_rbs %>% filter(day == 12),
       aes(y=Positive,x=Negative,color = demand)) +
  geom_point() +
  scale_color_gradient2(low = hue_pal()(3)[3],
                        high = hue_pal()(3)[1],mid = hue_pal()(3)[2],midpoint =1.5e+6,
                        name='Amino acid\ndemand')+
  theme_aviv + 
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label=amino_acid),hjust=0, vjust=0) + 
  coord_equal(xlim =c(0,pos_neg_lim_rbs),ylim =c(0,pos_neg_lim_rbs) ) + 
  facet_grid(lineage~RBS_Display)





ggplot(pos_neg_amino_acids_rbs %>% filter(day == 12),
       aes(y=ratio,x=demand)) +
  geom_point(size = 5) + 
  theme_aviv + 
  geom_text(aes(label=paste(' ',amino_acid)),hjust=0, vjust=0) +
  geom_smooth(method = 'lm') + 
  ylab('Ration of Positive amino acid frequency over Negative') +
  xlab('Amino acid demand') +
  facet_grid(lineage~RBS_Display) 
