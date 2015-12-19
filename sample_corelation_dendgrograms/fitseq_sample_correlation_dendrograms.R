# load libraries needed
require(dplyr)
require(ineq)
require(dendextend)
require(scales)

#load data
raw_data_location_2_mismatches = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\match_count_2_mismatches_generations.csv'
fitseq_raw_data_2_mismatches = read.csv(raw_data_location_2_mismatches,check.names=FALSE)
#rename the first column 
colnames(fitseq_raw_data_2_mismatches)[1] = 'design'
# #then throw it out
# fitseq_raw_data_mat_2_mismatches  <- data.matrix(fitseq_raw_data_2_mismatches %>% select(-design))
# fitseq_raw_data_norm_2_mismatches  <- data.frame(sweep(fitseq_raw_data_mat_2_mismatches,2,colSums(fitseq_raw_data_mat_2_mismatches),`/`))
# fitseq_raw_data_norm_log1p_2_mismatches  <- data.frame(log10(sweep(fitseq_raw_data_mat_2_mismatches+1,2,colSums(fitseq_raw_data_mat_2_mismatches+1),`/`)))

#create a color list for each of the lineages
lineage_colors <- hue_pal()(7)
names(lineage_colors) <- c('Ancestor', 'A', 'B', 'C', 'D', 'E', 'F')

#create a color list for each of the generations
generation_colors <- hue_pal()(8)
names(generation_colors) <- seq(0,196,28)
generation_colors["0"] <- lineage_colors['Ancestor']


#clustering 2 mismatches
#spearman
#calculate distance matrix of the spearman corelation using the cor function and the distance function with the distance being 1 - corelation
fitseq_spearman_distance_matrix_2_mismatches  <-  dist(1-cor(fitseq_raw_data_2_mismatches %>% select(-design) , method = "spearman"))
#calculate the upgma clustering using hclust and the average method
fitseq_spearman_clustering_2_mismatches  <- hclust(fitseq_spearman_distance_matrix_2_mismatches,method = 'average')
#turn the clusters into a dendrogram object
spearman.dend <- as.dendrogram(fitseq_spearman_clustering_2_mismatches)

#color the brancesh of the dedrogram by taking the names of the frequency matrix, ordering them by the dendrogram order, spliting on the _ underscore
# then taking the 1st column and then turning into a character, then using these characters as keys in the lineage color list and the resulting vector as the color
branch.col <- lineage_colors[as.character(t(as.data.frame(strsplit(names(fitseq_raw_data_2_mismatches)[-1][order.dendrogram(spearman.dend)],'_'))[1,]))]
# use the color of the vector created column to color the branches
spearman.dend <- color_branches(spearman.dend,col =branch.col )

#do the same thing as befor but take the second column of the split characters which is the generation of each branch, then use them as indices
# for the generation colors list, and use this to color the labels of the dendrograms
# uses the labels color function and then the set function with the branches lwd label
labels_colors(spearman.dend) <- generation_colors[as.character(t(as.data.frame(strsplit(names(fitseq_raw_data_2_mismatches)[-1][order.dendrogram(spearman.dend)],'_'))[2,]))]
spearman.dend <- set(spearman.dend, 'branches_lwd',3)

#plot the dendrogram
png('fitseq_spearman_upgma_clustering_2_mismatches.png',
    type="cairo",    units="in", width=14, height=6, pointsize=12, res=500)

plot(spearman.dend)

dev.off()

#use this to write the spearman corelation matrix
write.csv(cor(fitseq_raw_data %>% select(-design), method = "spearman"),
          paste0(results.dir,'spearman_correlation_matrix_2_mismatches.csv'))

#pearson
#exactly the same as before just use the pearson corelation method
fitseq_pearson_distance_matrix_2_mismatches  <-  dist(1-cor(fitseq_raw_data_2_mismatches %>% select(-design), method = "pearson"))

fitseq_pearson_clustering_2_mismatches  <- hclust(fitseq_pearson_distance_matrix_2_mismatches,method = 'average')

pearson.dend <- as.dendrogram(fitseq_pearson_clustering_2_mismatches)

branch.col <- lineage_colors[as.character(t(as.data.frame(strsplit(names(fitseq_raw_data_2_mismatches)[-1][order.dendrogram(pearson.dend)],'_'))[1,]))]
pearson.dend <- color_branches(pearson.dend,col =branch.col )

labels_colors(pearson.dend) <- generation_colors[as.character(t(as.data.frame(strsplit(names(fitseq_raw_data_2_mismatches)[-1][order.dendrogram(pearson.dend)],'_'))[2,]))]
pearson.dend <- set(pearson.dend, 'branches_lwd',3)

png('fitseq_pearson_upgma_clustering_2_mismatches.png',
    type="cairo",    units="in", width=14, height=6, pointsize=12, res=500)

plot(pearson.dend)

dev.off()

write.csv(cor(fitseq_raw_data %>% select(-design), method = "pearson"),
          paste0(results.dir,'pearson_correlation_matrix_2_mismatches.csv'))


