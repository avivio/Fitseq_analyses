require(FNN)
get.wilcox.string <-  function(pos,neg){
  pn_max <- max(c(max(pos),max(neg)))
  pn_min <- min(c(min(pos),min(neg)))
  pn_range = pn_max - pn_min
  greater <- wilcox.test(pos,neg,alternative = 'greater')$p.value
  greater.effect <- wilcox.test(pos,neg,greater = 'greater',conf.int=TRUE)$estimate / pn_range
  less <- wilcox.test(pos,neg,alternative = 'less')$p.value
  less.effect <- wilcox.test(pos,neg,alternative = 'less',conf.int=TRUE)$estimate /pn_range
  
  if(greater < less){
    return(paste('\n Positive > Negative\n p value =',formatC(greater,format='e',digits = 2),
                 'estimate = %', signif(greater.effect*100, digits = 3)))
  } else{
    return(paste('\n Positive < Negative\n p value =',formatC(less,format='e',digits = 2),
                 'estimate = %', signif(less.effect*100, digits = 3)))
  }
  
  
}


filter_n_fitness_residuals_csv <- function(n,fitseq_data_residuals){
  
  fitseq_data_pos_neg_lineages <- 
    fitseq_data_residuals %>%
    select(Name, day, lineage, pos_neg)

  fitseq_data_pos_neg_lineages <-  fitseq_data_pos_neg_lineages %>% 
    spread(lineage,pos_neg)
  
  n_pos <- fitseq_data_pos_neg_lineages %>%
    mutate(A= as.numeric(A=='Positive'),B= as.numeric(B=='Positive'),C= as.numeric(C=='Positive'),
           D= as.numeric(D=='Positive'),E=as.numeric( E=='Positive'),F= as.numeric(F=='Positive')) %>%
    mutate(sum_lineages = A+B+C+D+E+F,n_pos_neg = ifelse(sum_lineages >= n,'Positive',NA) ) %>%
    select(Name,day,n_pos_neg) %>%
    filter(!is.na(n_pos_neg))
  
  n_neg <- fitseq_data_pos_neg_lineages %>%
    mutate(A= as.numeric(A=='Negative'),B= as.numeric(B=='Negative'),C= as.numeric(C=='Negative'),
           D= as.numeric(D=='Negative'),E=as.numeric( E=='Negative'),F= as.numeric(F=='Negative')) %>%
    mutate(sum_lineages = A+B+C+D+E+F,n_pos_neg = ifelse(sum_lineages >= n,'Negative',NA) ) %>%
    select(Name,day,n_pos_neg) %>%
    filter(!is.na(n_pos_neg))
  n_pos_neg <-  bind_rows(n_pos,n_neg)
  
  
  
  return(collect(inner_join(fitseq_data_residuals,n_pos_neg, copy = TRUE)))
  
  
}


load_fitseq_data <- function(fitseq_data_location){
  
  fitseq_data = read.csv(fitseq_data_location)
  fitseq_data = transform(fitseq_data, Prot = as.numeric(as.character(Prot))
                          ,Bin_Pct_1 = as.numeric(as.character((Bin_Pct_1))),
                          Count_RNA = as.numeric(as.character((Count_RNA))))
  fitseq_data_tidy <- fitseq_data %>% select(-Count_A_DNA,-Count_A_RNA,-Count_B_DNA,-Count_B_RNA,
                                             -Bin_1,-Bin_2,-Bin_3,-Bin_4,-Bin_5,-Bin_6,-Bin_7,
                                             -Bin_8,-Bin_9,-Bin_10,-Bin_11,-Bin_12,-RNA_A,-RNA_B,
                                             -Bin_Pct_1,-Bin_Pct_2,-Bin_Pct_3,-Bin_Pct_4,-Bin_Pct_5,
                                             -Bin_Pct_6,-Bin_Pct_7,-Bin_Pct_8,-Bin_Pct_9,-Bin_Pct_10,
                                             -Bin_Pct_11,-Bin_Pct_12,-Insuff_Prot,-Insuff_DNA,
                                             -Insuff_RNA,-Fltr_BelowRange,-Fltr_AboveRange ,
                                             -Fltr_SetGood,-full_seq,-pep_sequence,-UTR,-CDS_seq,-Promoter_seq,-RBS_seq,-Promoter,
                                             -variable_seq,-full_peptide,-preRBS,-salis_status)
  fitseq_data_tidy  <- fitseq_data_tidy %>%   gather(sample, frequency, D_12:F_12)
  fitseq_data_tidy  <- fitseq_data_tidy %>% separate(sample, c("lineage", "day"), '_',convert = T)
  fitseq_data_tidy$RBS_Display <- factor(fitseq_data_tidy$RBS_Display,
                                         levels = c("Strong", "Mid", "Weak", "WT"))
  
  
  
  
  fitseq_data_tidy <- fitseq_data_tidy %>% 
    group_by(day,lineage) %>% 
    mutate(sum_sample= sum(frequency),norm_freq = frequency/sum_sample,
           sum_anc= sum(anc_1),norm_anc = anc_1/sum_anc,freq_norm_anc = norm_freq/norm_anc,
           sum_sample_1= sum(frequency+1),norm_freq_1 = (frequency+1)/sum_sample_1,
           sum_anc_1= sum(anc_1+1),norm_anc_1 = (anc_1+1)/sum_anc_1,freq_norm_anc_1 = norm_freq_1/norm_anc_1)
  
  
  return(fitseq_data_tidy)
}

get_data_for_correlation_category <- function(source_data,y_string,y_label,x_string,x_label,col_string,col_label){
  
  
  
  
  fitseq_xy <-  source_data %>%
    mutate_(y = y_string,x = x_string,col = col_string)
  
  fitseq_xy <- fitseq_xy %>%
    filter(!is.na(y) & !is.na(x) & !is.na(col))
  
  fitseq_xy_lm_all <- fitseq_xy %>% 
    ungroup() %>%
    group_by(day,lineage,col) %>% 
    do(mod = lm(y ~ x, data = .)) %>%
    mutate(slope = summary(mod)$coeff[2],
           intercept = summary(mod)$coeff[1],
           rmsd = sqrt(mean(mod$residuals^2)),
           p_value = summary(mod)$coeff[2,4],
           lm_string =paste('\n',  'slope = ', signif(slope,2),'\n',
                            'p-value = ' ,formatC(p_value,format='e',digits = 2),'\n',
                            'RMSD = ', signif(rmsd,2),'\n')) %>% 
    select(-mod)
  
  fitseq_xy <- left_join(fitseq_xy,fitseq_xy_lm_all,by = c('day','lineage','col'))
  
  fitseq_xy_lm_group <- fitseq_xy %>% 
    ungroup() %>%
    group_by(day,lineage,col,Promoter_Display,RBS_Display ) %>%
    do(mod = lm(y ~ x, data = .)) %>%
    mutate(slope = summary(mod)$coeff[2],
           intercept = summary(mod)$coeff[1],
           rmsd = sqrt(mean(mod$residuals^2)),
           p_value = summary(mod)$coeff[2,4],
           lm_string_group = paste('\n',  'slope = ', signif(slope,2),'\n',
                                   'p-value = ' ,formatC(p_value,format='e',digits = 2),'\n',
                                   'RMSD = ', signif(rmsd,2),'\n')) %>% 
    select(day,lineage,col,Promoter_Display,RBS_Display ,lm_string_group)
  
  
  fitseq_xy <- left_join(fitseq_xy,fitseq_xy_lm_group,by = c('day','lineage','col','Promoter_Display','RBS_Display'))
  
  
  fitseq_xy <-  fitseq_xy %>% 
    mutate(text_loc_y = ifelse(above_14,-Inf,Inf),
           text_loc_x= ifelse(above_14,Inf,-Inf),
           hor = ifelse(above_14,1,0),
           ver = ifelse(above_14,0,1))
  
  fitseq_xy <-  fitseq_xy %>% 
    mutate(prediction = intercept + slope*log2(Prot),
           fit_resid = log2(freq_norm_anc_1) -  prediction, 
           pos_neg = ifelse(fit_resid > 0, 'Positive', 'Negative') )
  return(fitseq_xy)
  
}




load_packages()
fiseq_data_location <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_data_with_meta_new_sd_17112015.csv'
fitseq_data_residuals <- load_fitseq_data(fiseq_data_location)
fitseq_data_residuals <-fitseq_data_residuals %>% 
  mutate(above_14 = log2(Prot)> 14) %>% 
  filter(log2(Prot)<17.5
         ,    Promoter_Display=='High'
         )

y_string <-  'log2(freq_norm_anc_1)'
y_label <- 'Log 2 of frequency in sample over frequency in ancestor'
x_string  <- 'log2(Prot)'
x_label <- 'Log 2 Protein level'
col_string  <- 'above_14'
# col_label <- 'Number of ribosomes\nper mRNA \n(calulated by TASEP)'
col_label <- 'Protein level\nabove 14'



# get correlation data into data set
fitseq_data_residuals <- get_data_for_correlation_category(fitseq_data_residuals,y_string,y_label,x_string,x_label,col_string,col_label)

fitseq_data_residuals <- filter_n_fitness_residuals_csv(5,fitseq_data_residuals)

fitseq_data_residuals <- fitseq_data_residuals %>%  
  filter(lineage =='A') %>% ungroup() %>% select(-lineage,-pos_neg,-fit_resid,-prediction,-p_value,-rmsd,-slope,-freq_norm_anc_1,-norm_anc_1,-sum_anc_1,-norm_freq_1,
                                                 -sum_sample_1,-freq_norm_anc,-norm_anc,-sum_anc,-norm_freq,-sum_sample,-frequency,-prediction)


variables <- 
  c('sd_harm_mean',	'sd_art_mean',	'sd_max',	'sd_min',	'sd_sum'	,'sd_median', 'sd_count')
labels <- 
  c('Shine Dalgarno speed harmonic mean','Shine Dalgarno speed arithemetic mean','Shine Dalgarno speed max value in CDS',
    'Shine Dalgarno speed minimum','Sum of Shine Dalgarno speed','Shine Dalgarno affninity speed median',
    'Number of positions with Shine Dalgarno speed above 0.5')
# variables <- 
#   c('TASEP_avgRate',	'TASEP_avgRiboNum',	'TASEP_density_avg',	'TASEP_density_0',	'TASEP_bottle_neck_position'	,'TASEP_bottle_neck_depth')
# labels <- 
#   c(  c('TASEP average translation rate',	'TASEP average # of ribosomes',	'TASEP density average per codon',	'TASEP density at start codon',
#         'TASEP bottleneck position'	,'TASEP bottleneck depth'))
# 

# fitseq_data_residuals_above_14 <- fitseq_data_residuals %>% filter(above_14 == TRUE)
# variables <- 
#   c('levy_mean',	'levy_median',	'levy_stdev',	'levy_sum',	'aggrescan_mean',	'aggrescan_median',	'aggrescan_stdev',	'aggrescan_sum',	'foldamyloid_mean',	
#     'foldamyloid_median', 'foldamyloid_stdev',	'foldamyloid_sum',	'ww_hydrophobicity_mean',	'ww_hydrophobicity_median',	'ww_hydrophobicity_stdev',	
#     'ww_hydrophobicity_sum', 'pep_instability')
# labels <- 
#   c('Mean aggregation propensity by Levy et al. 2012',	'Median aggregation propensity by Levy et al. 2012',
#     'Standard deviation of aggregation propensity by Levy et al. 2012',	'Sum of aggregation propensity values by Levy et al. 2012',
#     'Mean Aggrescan score',	'Median  Aggrescan score',	'Standard deviation of Aggrescan score',	'Sum of  Aggrescan scores',	'Mean FoldAmyloid score',	
#     'Median FoldAmyloid score', 'FoldAmyloid score standard deviation',	'Sum of FoldAmyloid scores',	'Mean Wimley and White hydrophobicity',
#     'Median Wimley and White hydrophobicity',	'Wimley and White hydrophobicity standard deviation',    'Sum of Wimley and White hydrophobicity scores',
#     'Instability index score')

# variables = c('CAI',	'tAI',	'CDS_GC',	'GC',	'dG',	'dG_noutr',	'dG_unif',	'salis_init',	'TASEP_avgRate',	'TASEP_avgRiboNum',	'TASEP_density_avg',	'TASEP_density_0',
# 'TASEP_bottle_neck_position',	'TASEP_bottle_neck_depth',	'sd_rbs',	'sd_max_score',	'sd_max_position',	'sd_mean',	'sd_sdev',	'sd_median',	'sd_count',
# 'sd_gfp_max_score',	'sd_gfp_max_position',	'sd_gfp_mean',	'sd_gfp_sdev',	'sd_gfp_median',	'sd_gfp_count',	'gc_1_cds',	'gc_2_cds',	'gc_3_cds',
# 'codon_TTT',	'codon_TTC',	'codon_TTA',	'codon_TTG',	'codon_TCT',	'codon_TCC',	'codon_TCA',	'codon_TCG',	'codon_TAT',	'codon_TAC',
# 'codon_TGT',	'codon_TGC',		'codon_TGG',	'codon_CTT',	'codon_CTC',	'codon_CTA',	'codon_CTG',	'codon_CCT',	'codon_CCC',
# 'codon_CCA',	'codon_CCG',	'codon_CAT',	'codon_CAC',	'codon_CAA',	'codon_CAG',	'codon_CGT',	'codon_CGC',	'codon_CGA',	'codon_CGG',	'codon_ATT',	
# 'codon_ATC',	'codon_ATA',	'codon_ATG',	'codon_ACT',	'codon_ACC',	'codon_ACA',	'codon_ACG',	'codon_AAT',	'codon_AAC',	'codon_AAA',	'codon_AAG',	'
# codon_AGT',	'codon_AGC',	'codon_AGA',	'codon_AGG',	'codon_GTT',	'codon_GTC',	'codon_GTA',	'codon_GTG',	'codon_GCT',	'codon_GCC',	'codon_GCA',	
# 'codon_GCG',	'codon_GAT',	'codon_GAC',	'codon_GAA',	'codon_GAG',	'codon_GGT',	'codon_GGC',	'codon_GGA',	'codon_GGG ',		'pep_pi',
# 'pep_cost',	'pep_demand',	'pep_mw',	'pep_aindex',	'pep_boman',	'pep_charge',	'pep_hmoment',	'pep_hydro',	'pep_instability',	'pep_tiny',	'pep_small',
# 'pep_aliphatic',	'pep_aromatic',	'pep_non_polar',	'pep_polar',	'pep_charged',	'pep_basic',	'pep_acidic',	'pep_demand',
# 'pep_Ala',	'pep_Cys',	'pep_Asp',	'pep_Glu',	'pep_Phe',	'pep_Gly',	'pep_His',	'pep_Ile',	'pep_Lys',	'pep_Leu',	'pep_Met',	'pep_Asn',	'pep_Pro',
# 'pep_Gln',	'pep_Arg',	'pep_Ser',	'pep_Thr',	'pep_Val',	'pep_Trp',	'pep_Tyr')

# variables <- 
#   c('CAI','tAI','CDS_GC','dG','dG_noutr',	'TASEP_bottle_neck_position',	'TASEP_bottle_neck_depth','pep_cost','pep_hydro','pep_polar',
#     'pep_mw','pep_basic','pep_acidic', 'sd_max_score','sd_mean',  'sd_count')
# labels <- 
#   c('CAI','tAI','Coding sequence GC %','Delta G','Delta G no UTR','Bottle neck position (calculated using TASEP)',
#     'Bottle neck depth (calculated using TASEP)','Peptide cost','Peptide hydrophobicity','Peptide polar amino acid content',
#     'Peptide molecular weight','Peptide basic amino acid content','Peptide acidic amino acid content',
#     'Variable region max SD affinity score ','Variable region mean SD affinity score','Variable region SD affinity score count')

fitseq_data_residuals_above_14 <- fitseq_data_residuals
names(labels) <- variables

for (i in 1:length(labels)){
  variable  <- names(labels)[i]
  label <- labels[i]
  print(label)
# for (variable in variables) {
#   print(variable)
  legend_string_negative  <- fitseq_data_residuals_above_14 %>%
    ungroup %>% 
    filter(n_pos_neg== 'Negative') %>% mutate_(new_variable = variable)  %>% select(new_variable) 
  legend_string_negative <- as.numeric( unlist(legend_string_negative))  
  legend_string_positive  <- fitseq_data_residuals_above_14 %>%
    ungroup %>% 
    filter(n_pos_neg== 'Positive') %>%
    mutate_(new_variable = variable)  %>% select(new_variable) 
  legend_string_positive <- as.numeric( unlist(legend_string_positive)) 
  
  wilcox_string_all <- get.wilcox.string(legend_string_positive,legend_string_negative)
  print(wilcox_string_all)
  
  png(paste0('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\new_sd\\pos_neg_histogram_',variable, '.png'),
      type="cairo",    units="in", width=12, height=10, pointsize=12, res=500)
  
  
  p <- ggplot(fitseq_data_residuals_above_14, aes_string(x = variable,fill ='n_pos_neg' )) +
    geom_density(alpha = 0.4) +
    xlab(label)+
    ylab('Density') +  
    guides(fill=guide_legend(title='Fitness\nresidual')) +
    theme_aviv
  print(p)
  dev.off()
  
}

# }



