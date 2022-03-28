

rm(list = ls())
library(tidyverse)
library(umap)
library(pheatmap)
library(fpc)
library(gplots)
library(RColorBrewer)
library(dbscan)
library(magrittr)
set.seed(158)
mlist = readxl::read_xlsx("84GenesDetails.xlsx", sheet = 2)
genExp = read.csv("TestData1.csv", row.names = 1)
pheno = read.csv("Pheno_chRCC_Onc20_30.csv", row.names = 1)
table(pheno$Histology2, pheno$DBUgrp)
pheno2 = pheno %>% filter(Histology2 == "chRCC" | Histology2 == "Oncocytoma")
pheno2 = as.data.frame(pheno2)
genExp2 = genExp[ , colnames(genExp) %in% pheno2$geo_accession]
genExp2 = t(genExp2) %>% as.data.frame()
genExp2 = genExp2[ , colnames(genExp2) %in% tolower(mlist$initial_alias)]
genExp2 = genExp2[ , match(tolower(mlist$initial_alias), colnames(genExp2))]
identical(pheno2$geo_accession, rownames(genExp2)) ## needs to be true

### dotplots

Group = ifelse(pheno2$DBUgrp == 1 & pheno2$Histology2 == "chRCC", "chRCC",
               ifelse(pheno2$DBUgrp == 2 & pheno2$Histology2 == "Oncocytoma", "Oncocytoma", 
                      ifelse(pheno2$DBUgrp == "Ambi", "Ambiguous", "Missclassified")))
table(Group)
identical(rownames(genExp2), pheno2$geo_accession)
qdf = cbind(Group, genExp2) %>% as.data.frame()
genes <- names(qdf)[2:62]

sndata = read.csv("DifferentialExpression.csv")

xvals<-jitter(rep(0,dim(qdf)[1]), factor=3)
xids<-ifelse(qdf$Group=='chRCC', 0.7,1)
gName = sndata[match(colnames(qdf), sndata$Gene), "SYMBOL"]
gName = gName[2:62]
bCuts = sndata[match(colnames(qdf), sndata$Gene), "optimal_cutpoint"]
bCuts = bCuts[2:62]


####### dot plot ggplot for the paper
qdf %>% 
  ggplot(aes(x = Group, y = `218780_at`, color = Group)) +
  geom_jitter() + theme_classic() +
  theme(text = element_text(size = 12, color = "black")) +
  theme(axis.text.x = element_text(angle = 90)) + ylab("HOOK2") +
  xlab("") + theme(legend.box.background = element_rect(linetype = 1))
qdf %>% 
  ggplot(aes(x = Group, y = `225291_at`, color = Group)) +
  geom_jitter() + theme_classic() +
  theme(text = element_text(size = 12, color = "black")) +
  theme(axis.text.x = element_text(angle = 90)) + ylab("PNPT1") +
  xlab("")
################

pdf("Dotplots.pdf", h=7, w=10)
par(oma = c(2,2,2,2))
par(mfrow=c(3,4), mar=c(4,3,0.5,.5))
for(i in 1:length(genes)) {
  ##i=1
  g=genes[i]
  y=qdf[,g]
  x<-xvals+xids
  clrid<-rep('black', length(y))
  clrid[which(qdf$Group=="Oncocytoma")] <- 'red'
  clrid[which(qdf$Group == "Missclassified")] <- "green"
  # table(clrid)
  lgtxt<-c(paste0(g), paste0(gName[i]))
  plot(x,y, xlim=c(0.5,1.2), type='p', xaxt='n',pch=20,font=2, xlab='', ylab='', col=clrid)
  abline(h = bCuts[i], col = "blue")
  title(paste(lgtxt, collapse="  "), '\n')
  axis(1, at=c(.7,1), c('chRCC','Oncocytoma'))
  
  #	xlab= " ", font=2,  las=1, xaxt='n', pch=16, font.lab=2, col="black"))
  
}
dev.off()

mlist2 = readxl::read_xlsx("84GenesDetails.xlsx", sheet = 3)
genExp3 = genExp[rownames(genExp) %in% tolower(mlist2$initial_alias), ]
genExp3 = genExp3[ , colnames(genExp3) %in% pheno2$geo_accession]


### check with unsupervised models
dir.create("SmallMods")
setwd("./SmallMods")
source("myfunctions.r")
genExp3 = t(genExp3) %>% as.data.frame()
allgenes = colnames(genExp3)
for (i in c("euclidean", "manhattan", "cosine")){
  #i="euclidean"
  pdf(paste0("01b_Umap_parameters_grid_search_",i,"_",gsub("-","",Sys.Date()),".pdf"))
  umap_parameters_gridsearch(genExp3,allgenes,colorlabel)
  dev.off()
}


##final parameters##
custom.config = umap.defaults
no_iters=1000
rand_gene_no=20
custom.config$n_neighbors = 30
custom.config$min_dist = 0.1
custom.config$metric = "cosine"

##all data stored in R file. don't have to run code below
# load("TCGAGlioma2015stand_DBU1000iters_20190611_results.R")

########02 1000 iterations of umap, store into data matrix, takes 25 minutes####
umap_position_matrix<-matrix(NA, nrow=dim(genExp3)[1],ncol=0)
Sys.time() -> tm
pdf(paste("01b_Umap no_iters=",no_iters,
          "n_neighbors=",custom.config$n_neighbors, 
          "min_dist=",custom.config$min_dist, 
          "no_genes=",rand_gene_no,
          "metric=",custom.config$metric,
          gsub("-","",Sys.Date()),
          ".pdf"))
par(mfrow = c(3, 3))
for (i in 1:no_iters){
  # i=1
  print(paste("Iteration #",i))
  brain3grp.umap.2 = umap(genExp3[,sample(allgenes,rand_gene_no)],custom.config)
  umap_position_matrix<-cbind.data.frame(umap_position_matrix,brain3grp.umap.2$layout)
  plot(brain3grp.umap.2$layout,pch=16,main=i,xlab = "",ylab = "")
}
dev.off()


Sys.time() - tm


no_points= dim(genExp3)[1]
trythingscount<-2
# 1) DBSCAN clustering algorithm
pdf(paste("01b_DBSCAN UMAP Reps",no_iters,rand_gene_no,"_",trythingscount,".pdf"))
par(mfrow=c(4,3),mar=c(0,0,0,0))
db_umap_res<-c(1:no_points)
no_actual_iters<-(dim(umap_position_matrix)[2]/2)
for (i in 1:no_actual_iters){
  # i=10
  x<-(i*2)-1
  # Compute DBSCAN using fpc package
  df<-umap_position_matrix[,x:(x+1)]
  #Method for determining the optimal eps value
  # dbscan::kNNdistplot(genExp3, k = 2)
  # abline(h = 1.35, lty = 2)
  
  db <- fpc::dbscan(df, eps = 2.6, MinPts = 20)
  plot(db, df, frame = FALSE,xlab="",ylab="")
  legend('topleft',paste(i),bty = 'n',text.col='red',text.font = 2)
  db_umap_res<-cbind.data.frame(db_umap_res,db$cluster)
  print(i)
}
dev.off()

db_umap_res<-db_umap_res[,-1]
colnames(db_umap_res)<-1:no_actual_iters

xx<-heatmap.2(as.matrix(t(db_umap_res)),trace = 'none',
              scale = 'none',
              col = brewer.pal(6,"BuPu"), ylab = "Iterations", xlab = "Samples",
              cexCol = 0.06, cexRow = 0.06)


#run this section if there is clustering inconsistency
hc<-hclust(dist(as.matrix(t(db_umap_res))))
memb <- cutree(hc, k = 2)
sum(memb==2)
xx = which(memb < 2)

xy = db_umap_res[ , xx]
heatmap.2(as.matrix(t(xy)),trace = 'none',
          scale = 'none',
          col = brewer.pal(6,"BuPu"),
          cexCol = 0.06, cexRow = 0.06, lhei = c(1,5), key.title = NA, 
          keysize = 0.8, density.info = "histogram", key.par = list(cex = 0.5),
          ylab = "Iterations", xlab = "Samples")
db_umap_res = xy

rownames(db_umap_res)<- pheno2$geo_accession

# 05 find ambiguous samples ########
dbumapgrps<-db_umap_res#[,-iter_rm]
all_maxgrp<-c()

for (j in 1:dim(dbumapgrps)[1]){
  # j=4
  datab<-table(as.character(dbumapgrps[j,]))
  maxgrp<-names(datab)[(which.max(datab))]
  maxperc<-max(table(as.character(dbumapgrps[j,]))*100/dim(dbumapgrps)[2])
  toadd<-c(maxgrp,maxperc)
  all_maxgrp<-rbind(all_maxgrp,toadd)
  print(j)
}
rownames(all_maxgrp)<-rownames(dbumapgrps)
colnames(all_maxgrp)<-c("Group","Proportion")
cutoff<-70
ambiloc<-which(all_maxgrp[,1]=="0"|as.numeric(all_maxgrp[,2])<cutoff)
all_maxgrp[ambiloc,]->ambig

#define DBUgrp to classify all samples
DBUgrp<-all_maxgrp[,1]
DBUgrp[ambiloc]<-"Ambi"
table(DBUgrp)

all_maxgrp<-cbind.data.frame(all_maxgrp,DBUgrp)


# 06 save everything ##
save(umap_position_matrix,
     db_umap_res,
     dbumapgrps,
     all_maxgrp, file="umapMatData.R")

#adding more to pheno table
newpheno<-merge.data.frame(pheno2,all_maxgrp,by="row.names")
newpheno_ordered<-newpheno[match(rownames(pheno2),newpheno$Row.names),]

table(newpheno_ordered$Histology2,newpheno_ordered$DBUgrp)



# cbind(newpheno_ordered$Row.names,rownames(pheno))
write.csv(newpheno_ordered,"Pheno_chRCC_Onc20_30.csv",row.names = F)





# plot ambiguous samples
pdf(paste0("01b_DBSCAN UMAP Reps",no_iters,rand_gene_no,"ambi_without_control.pdf"))
par(mfrow=c(4,3),mar=c(0,0,0,0))
newvec<-rep(1,1032)
newvec[ambiloc]<-2
for (i in 1:no_iters){
  # i=1
  x<-(i*2)-1
  df<-umap_position_matrix[,x:(x+1)]
  plot(df, col=newvec,xlab="",ylab="")
  legend('topright',paste(i),bty = 'n',text.col='red',text.font = 2)
  
}
dev.off()

# make dot plots

Group = ifelse(newpheno_ordered$final_cluster == 1 & newpheno_ordered$Histology2 == "chRCC" & newpheno_ordered$DBUgrp == 1, "chRCC",
               ifelse(newpheno_ordered$final_cluster == 2 & newpheno_ordered$Histology2 == "Oncocytoma" & newpheno_ordered$DBUgrp == 2, "Oncocytoma", "Missclassified"))


match(rownames(genExp3), newpheno_ordered$geo_accession)
qdf = cbind(Group, genExp3) %>% as.data.frame()
genes <- names(qdf)[2:31]


xvals<-jitter(rep(0,dim(qdf)[1]), factor=3)
xids<-ifelse(qdf$Group=='chRCC', 0.7,1)
gName = sndata[match(colnames(qdf)[-1], sndata$Gene), "SYMBOL"]
bCuts = sndata[match(colnames(qdf)[-1], sndata$Gene), "optimal_cutpoint"]

qdf %>% 
  ggplot(aes(x = Group, y = `225291_at`, color = Group)) +
  geom_jitter()

pdf("Dotplots30.pdf", h=7, w=10)
par(oma = c(2,2,2,2))
par(mfrow=c(3,4), mar=c(4,3,0.5,.5))
for(i in 1:length(genes)) {
  ##i=1
  g=genes[i]
  y=qdf[,g]
  x<-xvals+xids
  clrid<-rep('black', length(y))
  clrid[which(qdf$Group == "chRCC")] <- "black"
  clrid[which(qdf$Group=="Oncocytoma")] <- 'red'
  clrid[which(qdf$Group == "Missclassified")] <- "green"
  # table(clrid)
  lgtxt<-c(paste0(g), paste0(gName[i]))
  plot(x,y, xlim=c(0.5,1.2), type='p', xaxt='n',pch=20,font=2, xlab='', ylab='', col=clrid)
  abline(h = bCuts[i], col = "blue")
  title(paste(lgtxt, collapse="  "), '\n')
  axis(1, at=c(.7,1), c('chRCC','Oncocytoma'))
  
  #	xlab= " ", font=2,  las=1, xaxt='n', pch=16, font.lab=2, col="black"))
  
}
dev.off()

save.image("SmallMods.RData")

#### umap plots
u1 = umap(genExp2)
u1$layout %>% as.data.frame() %>% 
  ggplot(aes(V1, V2, color = as.factor(pheno2$Histology2))) +
  theme_classic() + geom_point()

save.image("SmallMods2.RData")


### figure redrawn
qdf2 = qdf %>% filter(Group == "chRCC" |
                        Group == "Oncocytoma")

qdf2 %>% 
  ggplot(aes(x = Group, y = `218780_at`, fill = Group)) +
  geom_boxplot() + theme_classic() +
  theme(text = element_text(size = 12, color = "black")) +
  scale_fill_manual(labels = c("chRCC", "Oncocytoma"),
                    values = c("yellowgreen", "violet")) +
  theme(axis.text.x = element_text(angle = 90)) + ylab("HOOK2") +
  xlab("") + theme(legend.box.background = element_rect(linetype = 1))
qdf2 %>% 
  ggplot(aes(x = Group, y = `225291_at`, fill = Group)) +
  geom_boxplot() + theme_classic() +
  scale_fill_manual(labels = c("chRCC", "Oncocytoma"),
                    values = c("yellowgreen", "violet")) +
  theme(text = element_text(size = 12, color = "black")) +
  theme(axis.text.x = element_text(angle = 90)) + ylab("PNPT1") +
  xlab("")

######## all gene boxplots

library(hgu133plus2.db)
colnames(qdf)
qdf2 = qdf %<>% filter(Group == "chRCC" | Group == "RO") 
u1 = umap((qdf2[, -1]),
          min_dist = 0.01,
          n_neighbors = 20,
          metric = "cosine")
qdf3 = cbind(u1$layout, qdf2)
colnames(qdf3)[1:2] = c("V1", "V2")

u1$layout %>% as.data.frame %>% 
  ggplot(aes(V1, V2, color = qdf2[ , 1])) +
  geom_point()

geneID = mlist2[match(colnames(qdf3)[4:33], tolower(mlist2$initial_alias)), "Name"]
colnames(qdf)[2:31] = geneID$Name
qdf$Group[qdf$Group == "Oncocytoma"] = "RO"



qdf %>% filter(Group != "Missclassified") %>% 
  dplyr::select(Group, KCNG3:AP1M2) %>%
  gather(Measure, Value, -Group) %>%
  ggplot(aes(x = Group, y = Value, color = Group)) +
  geom_jitter() +
  facet_wrap(~Measure
             , scales = "free_y") +
  theme_classic() +
  scale_color_manual(labels = c("chRCC", "RO"),
                    values = c("limegreen", "violet"))
quartz.save("Figure5.tiff", width = 11, height = 8, type = 'tiff', dpi = 300,
            device = dev.cur(), bg = 'white')


qdf3 %>% dplyr::select(1, 2, 5) %>% 
  ggplot(aes(V1, V2, color = PRDX3)) +
  geom_point() +
  scale_color_gradient(low = "green", high = "red")



for (i in 1:30) {
  assign(paste0("P", i), ggplot(qdf3, aes_string("V1", "V2", color = colnames(qdf3)[i+3])) +
    geom_point() + theme_classic() +
    scale_color_gradient(low = "grey", high = "dodgerblue4"))
  print(i)
}

# library(gridExtra)
grid.arrange(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,
             P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,
             P21,P22,P23,P24,P25,P26,P27,P28,P29,P30,nrow = 6, ncol = 5)
