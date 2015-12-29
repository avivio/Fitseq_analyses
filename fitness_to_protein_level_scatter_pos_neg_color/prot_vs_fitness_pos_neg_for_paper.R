load.packages()


get.data.for.correlation.category <- function(source.data,y.string,y.label,x.string,x.label,col.string,col.label){
  
  
  
  
  fitseq.xy <-  source.data %>%
    mutate_(y = y.string,x = x.string,col = col.string)
  
  fitseq.xy <- fitseq.xy %>%
    filter(!is.na(y) & !is.na(x) & !is.na(col))
  
  fitseq.xy.lm.all <- fitseq.xy %>% 
    ungroup() %>%
    group_by(day,lineage,col) %>% 
    do(mod = lm(y ~ x, data = .)) %>%
    mutate(slope = summary(mod)$coeff[2],
           intercept = summary(mod)$coeff[1],
           rmsd = sqrt(mean(mod$residuals^2)),
           p.value = summary(mod)$coeff[2,4],
           lm.string =paste('\n',  'slope = ', signif(slope,2),'\n',
                            'p-value = ' ,formatC(p.value,format='e',digits = 2),'\n',
                            'RMSD = ', signif(rmsd,2),'\n')) %>% 
    select(-mod)
  
  fitseq.xy <- left_join(fitseq.xy,fitseq.xy.lm.all,by = c('day','lineage','col'))
  
  fitseq.xy.lm.group <- fitseq.xy %>% 
    ungroup() %>%
    group_by(day,lineage,col,Promoter.Display,RBS.Display ) %>%
    do(mod = lm(y ~ x, data = .)) %>%
    mutate(slope = summary(mod)$coeff[2],
           intercept = summary(mod)$coeff[1],
           rmsd = sqrt(mean(mod$residuals^2)),
           p.value = summary(mod)$coeff[2,4],
           lm.string.group = paste('\n',  'slope = ', signif(slope,2),'\n',
                                   'p-value = ' ,formatC(p.value,format='e',digits = 2),'\n',
                                   'RMSD = ', signif(rmsd,2),'\n')) %>% 
    select(day,lineage,col,Promoter.Display,RBS.Display ,lm.string.group)
  
  
  fitseq.xy <- left_join(fitseq.xy,fitseq.xy.lm.group,by = c('day','lineage','col','Promoter.Display','RBS.Display'))
  
  
  fitseq.xy <-  fitseq.xy %>% 
    mutate(text.loc.y = ifelse(above.14,-Inf,Inf),
           text.loc.x= ifelse(above.14,Inf,-Inf),
           hor = ifelse(above.14,1,0),
           ver = ifelse(above.14,0,1))
  
  fitseq.xy <-  fitseq.xy %>% 
    mutate(prediction = intercept + slope*log2(Prot),
           fit.resid = log2(freq.norm.anc.1) -  prediction, 
           pos.neg = ifelse(fit.resid > 0, 'Positive', 'Negative') )
  return(fitseq.xy)
  
}


fitseq.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_data_with_meta_day_12_no_umi_05_08.csv'
fitseq.data.residuals  <- load.fitseq.data(fitseq.data.location)

fitseq.data.residuals <- fitseq.data.residuals %>% 
  mutate(above.14 = log2(Prot)> 14)

#adding something tiny for commit
#define parameters for correlation 

y.string <-  'log2(freq.norm.anc.1)'
y.label <- 'Log 2 of frequency in sample over frequency in ancestor'
x.string  <- 'log2(Prot)'
x.label <- 'Log 2 Protein level'
col.string  <- 'above.14'
# col.label <- 'Number of ribosomes\nper mRNA \n(calulated by TASEP)'
col.label <- 'Protein level\nabove 14'




fitseq.data.residuals <- fitseq.data.residuals %>% 
  filter(     log2(Prot) < 17.5,
              
              Promoter.Display=='High')

# get correlation data into data set
fitseq.data.residuals <- get.data.for.correlation.category(fitseq.data.residuals,y.string,y.label,x.string,x.label,col.string,col.label)



#decide what data set you're working on

fit.summary  <- fitseq.data.residuals %>% 
  group_by(day,lineage,pos.neg) %>% 
  summarise(percentile= quantile(abs(fit.resid),probs = c(0.2))) 

# fitseq.data.residuals <-   left_join(fitseq.data.residuals,fit.summary,by = c('day','lineage','pos.neg'))
# 
# fitseq.data.residuals <- fitseq.data.residuals %>% 
#   filter( 
#     abs(fit.resid) >  percentile
#   )



fitseq.data.residuals.above.14 <- fitseq.data.residuals %>% filter(above.14 == TRUE)



day.number = 12
lineage.letter = 'A'

fitseq.current.sample <-  fitseq.data.residuals %>% filter(day == day.number, lineage == lineage.letter)
fitseq.current.sample  <- fitseq.current.sample  %>% mutate(pos.neg.null = ifelse(above.14==T,pos.neg, ''))
fitseq.current.sample  <- fitseq.current.sample  %>% mutate(pred.null = ifelse(above.14==T,prediction, NA))


png('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\for_thesis\\prot_vs_fit_pos_neg_colors.png',
    type="cairo",    units="in", width=12, height=10, pointsize=12, res=500)    



p <- ggplot(fitseq.current.sample,aes(x= log2(Prot), y = log2(freq.norm.anc.1),color  = pos.neg.null)) +
  geom_point(alpha = 0.3)   +
  # xlim(x.limits) +
  coord_cartesian(ylim = c(-5.5,4)) +
  ylab( 'Log 2 normalized frequency') + 
  xlab('Log 2 Protein level') +
  scale_colour_manual( values = c("grey70",hue_pal()(2)[1],hue_pal()(2)[2]), breaks = c('Positive','Negative')) +
  geom_line(aes(x = x,y= pred.null), color = 'grey50',inherit.aes=FALSE,
            size = 1.2,show_guide  = F)+
  # geom_abline(aes(slope = slope, intercept = intercept), color = 'dodgerblue', size = 1.2) +
  #     geom_smooth(aes_string(linetype = col.string),show_guide  = F,
  #                 color = 'red', size = 1.2, method = 'lm') +
  guides(colour = guide_legend(override.aes = list(alpha = 1),title = 'Fitness\nresidual')) +
  theme_aviv
print(p)

dev.off()