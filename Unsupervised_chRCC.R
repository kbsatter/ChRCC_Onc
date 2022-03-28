
rm(list = ls())


set.seed(158)
## library
library(tidyverse)
library(magrittr)
library(umap)
library(ggplot2)
library(dbscan)
library(RColorBrewer)
library(pheatmap)
library(gplots)
## load multi-core functions
source("myfunctions.R")


### load data files
genExp = read.csv("TestData1.csv", row.names = 1)
geneID = genExp$GENEID2
genExp %<>% select(-1)
genExp = t(genExp) %>% as.data.frame()
pheno = read.csv("phenoData.csv")
pheno2 = pheno[pheno$Histology2 != "N", ]
genExp2 = genExp[ match(pheno2$geo_accession, rownames(genExp)), ]
allgenes = colnames(genExp)


set.seed(142)
### check UMAP data
u1 = umap(genExp2,
          min_dist = 0.01,
          n_neighbor = 20,
          metric = "cosine")
Histology = pheno2$Histology2

u1$layout %>% as.data.frame() %>% 
  ggplot(aes(V1,V2, color = Histology)) + geom_point() +
  theme_classic() + theme(text = element_text(size = 15)) +
  scale_color_manual(labels = c("chRCC", "Oncocytoma"),
                    values = c("yellowgreen", "violet")) +
  xlab("UMAP_1") + ylab("UMAP_2")
### Grid-Search
### get the best combination
dir.create("./Unsupervised/New")
setwd("./Unsupervised/New")
for (i in c("euclidean", "manhattan", "cosine")){
  #i="euclidean"
  pdf(paste0("01b_Umap_parameters_grid_search_",i,"_",gsub("-","",Sys.Date()),".pdf"))
  umap_parameters_gridsearch(genExp,allgenes,colorlabel)
  dev.off()
}

## final parameters
## final dataset
#2b. Run 1000 umap iterations
umap_position_matrix2 = umapiterations(dataset = genExp2,
                                      genes = allgenes,
                                      no_iters = 1000,
                                      rand_gene_no = 1000,
                                      n_neighbors = 20,
                                      min_dist = 0.01,
                                      metric = "manhattan")

### check plots
par(mfrow = c(3,3))
for(i in 1:9){
  plot_DBUiter(umap_position_matrix2,i)
}

dev.off()
#3c. optimize DBSCAN parameters

dbutest(umap_position_matrix = umap_position_matrix2,
        iter = 50,
        k = 5,
        eps = 2)

#3d. run DBSCAN on all umap iters
db_umap_res<-dbuiters(umap_position_matrix = umap_position_matrix2,
                      k = 5,
                      eps = 2)

############### Consensus clustering ##################
#Find pairwise co-grouping frequencies for each sample
Mfreq <- group_consensus_freq(umap_position_matrix2,db_umap_res)
#make group calls based on clustering pairwise frequencies
k = 2
pheatmap::pheatmap(Mfreq,scale = "none",
                   clustering_distance_cols = "manhattan",
                   clustering_distance_rows = "manhattan",
                   cutree_rows = k,
                   cutree_cols = k)

consens_calls<-makegrpcalls_clust(Mfreq,k)
table(consens_calls$samp_clust6)

consens_calls$final_cluster<-change_grp_names4(consens_calls$samp_clust6)
consens_calls$final_cluster[consens_calls$samp_clus_div_res>7]<-NA
table(consens_calls$final_cluster)

x = 50
plot(umap_position_matrix2[,x:(x+1)], col=consens_calls$final_cluster,xlab="",ylab="",pch=19,cex=0.7,
     main = "UMAP iteration 50")

##
pheno = cbind(pheno, consens_calls)
table(pheno2$Histology2, pheno2$final_cluster)

write.csv(pheno2, file = "PhenoData_unsupervised2.csv")

save.image("Unsupervised2.RData")

######
library(ggalluvial)
Unsupervised = ifelse(pheno2$final_cluster == 4, "DBU1", "DBU2")
pheno2 = cbind(pheno2, Unsupervised)
ggplot(data = pheno2[ , -1],
       aes(axis1 = Batch, axis2 = Histology2, axis3 = Unsupervised)) +
  scale_x_discrete(limits = c("Study", "Diagnosis", "Unsupervised"))+
  xlab("") +
  geom_alluvium(aes(fill = Histology2)) +
  scale_fill_manual(values = c("limegreen","violet")) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() + theme(text = element_text(size = 12))
  ggtitle("chRCC-Onc",
          "Unsupervised and Supervised Model")

load("/Users/khaled/Downloads/chRCC/Unsupervised/Unsupervised2.RData")
  