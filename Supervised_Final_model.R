
rm(list = ls())
library(tidyverse)
library(magrittr)
library(umap)
library(ComplexHeatmap)
library(circlize)
### important file pheno2, genexp3, newpheno_ordered, mlist2
load("./SmallMods/SmallMods.RData")
pheno3 = newpheno_ordered[newpheno_ordered$DBUgrp == 1 &
                            newpheno_ordered$Histology2 == "chRCC" &
                            newpheno_ordered$final_cluster == 1 |
                            newpheno_ordered$DBUgrp == 2 &
                            newpheno_ordered$Histology2 == "Oncocytoma" &
                            newpheno_ordered$final_cluster == 2, ]

Group = pheno3$DBUgrp

genExp3.1 = genExp3[rownames(genExp3) %in% pheno3$geo_accession, ]
genExp3.1 = cbind.data.frame(Group, genExp3.1)

library(caret)
library(randomForest)
library(mlbench)
library(caret)
library(e1071)
index = createDataPartition(genExp3.1$Group, p = 0.7, list = F)
Xtrain = genExp3.1[index, ]
Xtest = genExp3.1[-index, ]
control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=3)
metric = "Accuracy"
set.seed(142)
rf_default <- train(Group ~., 
                    data = Xtrain, 
                    method='rf', 
                    metric='Accuracy', 
                    tuneLength = 20,
                    trControl=control)

mean(predict(rf_default, Xtest) == Xtest$Group)


valdata = read.csv("GSE12090exp.csv", row.names = 1)
colnames(valdata) = gsub("X", "", colnames(valdata))
valdata = valdata[ , colnames(valdata) %in% tolower(mlist2$initial_alias)]
valPheno = read_csv("/Users/khaled/Downloads/final_chrcc/old files/Pheno_chRCC_Onc.csv")

valPheno %<>% filter(batch == "GSE12090")
valPheno %<>% mutate(Group = ifelse(Dx2 == "chRCC", 1,2))
valdata = cbind.data.frame(valPheno$Group, valdata)
mean(predict(rf_default, valdata) == valdata$`valPheno$Group`)
cbind(predict(rf_default, valdata), valdata$`valPheno$Group`)

testres = predict(rf_default, Xtest, reference = Xtest$Group, type = "prob")
valres = predict(rf_default, valdata, reference = valdata$`valPheno$Group`, type = "prob")


library(pROC)
par(pty = "s")
par(mfrow = c(1,2))
plot(roc((as.numeric(Xtest$Group)-1), testres$`1`), lty = 1, lwd = 2, col = "black")

plot(roc((as.numeric(valdata$`valPheno$Group`)-1), valres$`1`), lty = 2, lwd = 2, col = "red")
roc((as.numeric(Xtest$Group)-1), testres$`1`, plot = T, legacy.axis = T,
    percent = T, print.auc = T, main = "Test Data")
roc((as.numeric(valdata$`valPheno$Group`)-1), valres$`1`, plot = T, legacy.axis = T,
    percent = T, print.auc = T, main = "Validation Data")

### 
colnames(valdata)[1] = "Group"
set.seed(142)
Histology = ifelse(valdata$Group == 1, "chRCC", "Oncocytoma")
val_umap = umap(valdata[, -1],
                n_neighbors = 10,
                min_dist = 0.01,
                metric = "cosine")
val_umap$layout %>% as.data.frame() %>% 
  ggplot(aes(V1, V2, color = Histology)) +
  geom_point() + theme_classic() +
  theme(text = element_text(size = 12)) +
  scale_color_manual(labels = c("chRCC", "Oncocytoma"),
                     values = c("yellowgreen", "violet"))


top_anno = HeatmapAnnotation(Histology = pheno2$Histology2,
                             Batch = pheno2$Batch,
                             Unsupervised = pheno2$DBUgrp,
                             col = list(Histology = c("chRCC" = "yellowgreen",
                                                      "Oncocytoma" = "violet")))

tset = genExp3 %>% scale(center = T, scale = T) %>% t() %>%
  Heatmap(cluster_columns = T, cluster_rows = T,
          col = colorRamp2(c(-2,0,2), c("blue", "grey", "red")), 
          row_names_side = "left", row_dend_side = "right", 
          row_names_gp = gpar(cex = .5), 
          row_dend_width = unit(2, "cm"), 
          column_names_gp = gpar(cex = 0.1),
          clustering_distance_rows = "maximum", 
          clustering_method_rows = "ward.D2", top_annotation = top_anno,
          km = 2, column_km = 2,
          name = "Training Set")

val_ano = HeatmapAnnotation(Histology = Histology,
                            col = list(Histology = c("chRCC" = "yellowgreen",
                                                     "Oncocytoma" = "orchid2")))
  
vset = valdata %>% select(-1) %>% scale(center = T, scale = T) %>% t %>% 
  Heatmap(cluster_columns = T, cluster_rows = T,
          col = colorRamp2(c(-2,0,2), c("blue", "grey", "red")), 
          row_names_side = "left", row_dend_side = "right", 
          row_names_gp = gpar(cex = .5), 
          row_dend_width = unit(3, "cm"), 
          column_names_gp = gpar(cex = 0.1),
          clustering_distance_rows = "maximum", 
          clustering_method_rows = "ward.D2", top_annotation = val_ano,
          km = 2, column_km = 2,
          name = "Validation Set")
tset + vset

save.image("ValidationMods.RData")


#### tcga data
tcga_exp = read.csv("./TCGA validation/TCGA-KICH.htseq_fpkm-uq.tsv.gz", sep = "\t", row.names = 1)
tcga_pheno = read.csv("./TCGA validation/TCGA-KICH.GDC_phenotype.tsv.gz", sep = "\t")
