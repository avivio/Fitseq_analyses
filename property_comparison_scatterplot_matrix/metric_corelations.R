
filter.n.fitness.residuals <- function(n,fitseq.data.residuals){
  
  fitseq.data.pos.neg.lineages <- 
    fitseq.data.residuals %>%
    select(Name, day, lineage,pos.neg) %>% 
    spread(lineage,pos.neg)
  
  n.pos <- fitseq.data.pos.neg.lineages %>%
    mutate(A= as.numeric(A=='Positive'),B= as.numeric(B=='Positive'),C= as.numeric(C=='Positive'),
           D= as.numeric(D=='Positive'),E=as.numeric( E=='Positive'),F= as.numeric(F=='Positive')) %>%
    mutate(sum.lineages = A+B+C+D+E+F,n.pos.neg = ifelse(sum.lineages >= n,'Positive',NA) ) %>%
    select(Name,day,n.pos.neg) %>%
    filter(!is.na(n.pos.neg))
  
  n.neg <- fitseq.data.pos.neg.lineages %>%
    mutate(A= as.numeric(A=='Negative'),B= as.numeric(B=='Negative'),C= as.numeric(C=='Negative'),
           D= as.numeric(D=='Negative'),E=as.numeric( E=='Negative'),F= as.numeric(F=='Negative')) %>%
    mutate(sum.lineages = A+B+C+D+E+F,n.pos.neg = ifelse(sum.lineages >= n,'Negative',NA) ) %>%
    select(Name,day,n.pos.neg) %>%
    filter(!is.na(n.pos.neg))
  n.pos.neg <-  bind_rows(n.pos,n.neg)
  
  return(inner_join(fitseq.data.residuals,n.pos.neg))
  
  
}

require(GGally)
load.packages()
fitseq.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_data_with_meta_day_12_no_umi_05_08.csv'

# data_set = 'whole_library'
# fitseq.data = read.csv(fitseq.data.location)
data_set = 'above_5'
fitseq.data  <- load.fitseq.data(fitseq.data.location)
fitseq.data <- fitseq.data %>% 
  mutate(above.14 = log2(Prot)> 14)

y.string <-  'log2(freq.norm.anc.1)'
y.label <- 'Log 2 of frequency in sample over frequency in ancestor'
x.string  <- 'log2(Prot)'
x.label <- 'Log 2 Protein level'
col.string  <- 'above.14'
# col.label <- 'Number of ribosomes\nper mRNA \n(calulated by TASEP)'
col.label <- 'Protein level\nabove 14'


fitseq.data <- fitseq.data %>% 
  filter(     log2(Prot) < 17.5,Promoter.Display=='High')
fitseq.data <- get.data.for.correlation.category(fitseq.data,y.string,y.label,x.string,x.label,col.string,col.label)
fitseq.data <- filter.n.fitness.residuals(5,fitseq.data)
fitseq.data <- fitseq.data %>% filter(above.14 == TRUE) %>% 
  filter(lineage =='A') %>% 
  ungroup() %>% 
  select(-lineage)




fitseq.data <- fitseq.data %>% 
  mutate(salis.log = log2(salis.init)) %>% 
  select(CAI,tAI,CDS.GC,GC, dG, dG.noutr,
         salis.log,contains('TASEP'),
         pep.cost, pep.boman, pep.charge, pep.hydro,pep.polar,pep.basic,pep.acidic ,pep.mw,
         sd.max.score, sd.max.position,  sd.mean,sd.sdev, sd.median, sd.count  ) 
  
column_list = list(1:6,7:13,14:21,22:27,c(6,13,14,17,22,27))
label_list = list(c('','','','','',''),c('','','','','','',''),c('','','','','','','',''),c('','','','','',''),
                  c('','','','','',''))
name_list = list('codon_bias','ribosome_flow_simulation','peptide_properties','sd_affinity','all')
result.dir <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\comparing_metrics\\'
for (i in 1:length(column_list)){
  print(name_list[[i]])
  png(paste0(result.dir,data_set,'_scatter_plot_matrix_',name_list[[i]], '.png'),
      type="cairo",    units="in", width=20, height=16, pointsize=12, res=500)   
  p <- ggpairs(fitseq.data,column_list[[i]],columnLabels = label_list[[i]], upper = list(params = c(size = 10))) +
     # theme_minimal() + 
    theme(  axis.text=element_text(size=18))
  print(p)
  dev.off()

}