
#load packages
require(dplyr)
require(ggplot2)
require(ggvis)
require(scales)
require(tidyr)


#set publication theme
theme_aviv <-     theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
         strip.text.x = element_text(size=16,face="bold"), strip.text.y = element_text(size=16,face="bold"), 
         plot.title = element_text(size = 25,face = "bold"),
         legend.text = element_text(size=16),legend.title = element_text(size=18,face="bold")	)

#load data and transpose it so it can be handaled
summary_data_location = 'final_summary_generations.csv'
summary_data = t(read.csv(summary_data_location,check.names=FALSE,row.names = 1,as.is = T))

#put data into the correct type so we can use it
summary_data <- transform(summary_data,
                          all_reads = as.numeric(as.character(all_reads)),
                          merged_reads= as.numeric(as.character(merged_reads)),
                          trimmed_reads= as.numeric(as.character(trimmed_reads)),
                          found_designs_0_mismatches = as.numeric(as.character(found_designs_0_mismatches)),
                          design_count_pre_umi_0_mismatches = as.numeric(as.character(design_count_pre_umi_0_mismatches)),
                          design_count_post_umi_0_mismatches = as.numeric(as.character(design_count_post_umi_0_mismatches)),
                          found_designs_1_mismatches = as.numeric(as.character(found_designs_1_mismatches)),
                          design_count_pre_umi_1_mismatches = as.numeric(as.character(design_count_pre_umi_1_mismatches)), 
                          design_count_post_umi_1_mismatches = as.numeric(as.character(design_count_post_umi_1_mismatches)),
                          found_designs_2_mismatches = as.numeric(as.character(found_designs_2_mismatches)),   
                          design_count_pre_umi_2_mismatches = as.numeric(as.character(design_count_pre_umi_2_mismatches)),
                          design_count_post_umi_2_mismatches = as.numeric(as.character(design_count_post_umi_2_mismatches)),
                          found_designs_3_mismatches = as.numeric(as.character(found_designs_3_mismatches)),
                          design_count_pre_umi_3_mismatches = as.numeric(as.character(design_count_pre_umi_3_mismatches)), 
                          design_count_post_umi_3_mismatches = as.numeric(as.character(design_count_post_umi_3_mismatches)),
                          found_designs_4_mismatches = as.numeric(as.character(found_designs_4_mismatches)),
                          design_count_pre_umi_4_mismatches = as.numeric(as.character(design_count_pre_umi_4_mismatches)), 
                          design_count_post_umi_4_mismatches = as.numeric(as.character(design_count_post_umi_4_mismatches)),          
                          found_designs_all = as.numeric(as.character(found_designs_all)),   
                          design_count_pre_umi_all = as.numeric(as.character(design_count_pre_umi_all)),
                          design_count_post_umi_all = as.numeric(as.character(design_count_post_umi_all)),
                          lineage = as.character(lineage),
                          generation = as.numeric(as.character(generation))
)

#create a column whih is th percent of found designs
summary_data_design_coverage_2_mismatches <- summary_data %>% mutate(found_design_perc = found_designs_2_mismatches/14234 )  
#order the factor of the lineages so it looks nice in legend
summary_data_design_coverage_2_mismatches$lineage <- factor(summary_data_design_coverage_2_mismatches$lineage,
                                                            levels = c("Ancestor", "A", "B", "C", "D", "E", "F"))

#plot figure

#straight forward x = generation, y percent designs found, color lineage
png('percent_of_library_covered_over_time.png',
    type="cairo",    units="in", width=10, height=6, pointsize=12, res=500)    
ggplot(summary_data_design_coverage_2_mismatches, aes(x = generation,y = found_design_perc, color = lineage)) + 
  geom_point(size = 3) +
  xlab("Generations") +
  ylab("Percent of library covered") + 
  scale_x_continuous(breaks=c(0, 28, 56, 84, 112,140,168,196)) +
  theme_aviv

dev.off()