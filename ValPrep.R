
rm(list = ls())

library(umap)
library(ggplot2)
library(magrittr)
library(tidyverse)
## GS2
GS2 = read.csv("GS2.csv", row.names = 1)
GS1 = read.csv("GS1.csv", row.names = 1)
GS2GID = GS1[match(GS2$x, tolower(GS1$initial_alias)), "Name"]

### GSE12090
valdata = read.csv("GSE12090exp.csv", row.names = 1)
colnames(valdata) = gsub("X", "", colnames(valdata))
valdata = valdata[ , colnames(valdata) %in% GS2$x]
valPheno = read_csv("/Users/khaled/Downloads/final_chrcc/old files/Pheno_chRCC_Onc.csv")

valPheno %<>% filter(batch == "GSE12090")
valPheno %<>% mutate(Group = ifelse(Dx2 == "chRCC", 1,2))

## TCGA
tcga = read.csv("./TCGA validation/TCGA_valData.csv", row.names = 1)
tcga = tcga[ , GS2GID]
GS2PID = tolower(GS1[match(colnames(tcga), GS1$Name), "initial_alias"])
colnames(tcga) = GS2PID
tpheno = read.csv("./TCGA validation/TCGA_valPheno.csv")

### GSE15641
load("GSE15641Data.Rdata")
match(GS2$x, rownames(df2))
df2 = df2[rownames(df2) %in% GS2PID, ]
df2 %<>% t %>% as.data.frame()
pheno41 %<>% filter(Histology != "N")
df2 = df2[rownames(df2) %in% pheno41$`!Sample_geo_accession`, ]
df2 %<>% log2
## tcga, 12090, 15641
dim(df2) ; dim(tcga); dim(valdata)
df2 = df2[ , order(colnames(df2))]
tcga = tcga[ , order(colnames(tcga))]
valdata = valdata[ , order(colnames(valdata))]
identical(colnames(df2), colnames(tcga))
identical(colnames(df2), colnames(valdata))

ValExp = rbind.data.frame(df2, tcga, valdata)

#### valpheno
tpheno = tpheno %>% add_column(Batch = rep("KICH"),
                               Histology = rep("chRCC"))
rownames(tpheno)
tpheno2 = tpheno[ , c("X", "Batch", "Histology")]

valPheno = valPheno %>% dplyr::select(1,4,5)
colnames(valPheno) = c("X", "Batch", "Histology")

pheno41 %<>% add_column(Batch = "GSE15641") %>% dplyr::select(1, 5, 6)
pheno41 = pheno41[ , c(1,3,2)]
pheno41$Histology[which(pheno41$Histology == "RO")] = "Oncocytoma"
colnames(pheno41)[1] = "X"

valPheno2 = rbind.data.frame(valPheno, tpheno2, pheno41)
write_csv(valPheno2, "ValidationPhedata.csv")


### combat
valPheno2 %<>% dplyr::filter(Batch != "GSE15641")
ValExp = ValExp[ rownames(ValExp) %in% valPheno2$X, ]
valPheno2 = valPheno2[match(rownames(ValExp), valPheno2$X), ]
identical(rownames(ValExp), valPheno2$X)

Batch = as.factor(valPheno2$Batch)
Histology = as.factor(valPheno2$Histology)

pca1 = prcomp(ValExp, scale. = T)
score1 = pca1$x %>% as.data.frame()
ggplot(data = score1, aes(x = PC1, y = PC2, label = rownames(score1))) + 
  geom_point(aes(color = Batch)) + theme_classic() +
  theme(axis.text = element_text(size = 15), legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

ggplot(data = score1, aes(x = PC1, y = PC2, label = rownames(score1))) + 
  geom_point(aes(color = Histology)) + theme_classic() +
  theme(axis.text = element_text(size = 15), legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  scale_color_manual(values = c("limegreen", "royalblue3", "violet"))


u1 = umap(ValExp)
u1$layout %>% as.data.frame() %>%  
  ggplot(aes(x = V1, y = V2, color = Batch)) + 
  geom_point() + theme_classic() +
  theme(axis.text = element_text(size = 15), legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

u1$layout %>% as.data.frame() %>%  
  ggplot(aes(x = V1, y = V2, color = Histology)) + 
  geom_point() + theme_classic() +
  theme(axis.text = element_text(size = 15), legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

library(sva)
df1 = ComBat(dat = t(ValExp), mod = Histology, batch = Batch, par.prior = T, mean.only = F) ####ComBat
pca2 = prcomp(t(df1), scale. = T) 
scores2 = pca2$x %>% as.data.frame() 
ggplot(data = scores2, aes(x = PC1, y = PC2, label = rownames(scores2))) + 
  geom_point(aes(color = Batch)) + theme_classic() +
  theme(axis.text = element_text(size = 15), legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

ggplot(data = scores2, aes(x = PC1, y = PC2, label = rownames(scores2))) + 
  geom_point(aes(color = Histology)) + theme_classic() +
  theme(axis.text = element_text(size = 15), legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  scale_color_manual(values = c("limegreen", "violet"))
df_scaled = scale(t(df1), center = T, scale = T)
u2 = umap(t(df1), 
          n_neighbors = 20,
          min_dist = 0.01,
          metric = "cosine",
          input = "data") 
u2$layout %>% as.data.frame() %>%  
  ggplot(aes(V1,V2, color = Batch)) + 
  geom_point() + theme_classic() +
  theme(axis.text = element_text(size = 15), legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))
u2$layout %>% as.data.frame() %>%  
  ggplot(aes(V1,V2, color = Histology)) + 
  geom_point() + theme_classic() +
  theme(axis.text = element_text(size = 15), legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  scale_color_manual(values = c("limegreen", "violet"))
