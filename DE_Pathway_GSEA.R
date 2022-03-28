

rm(list = ls())
load("DE_Data2_7152021.RData")
## library
library(tidyverse)
library(ReactomePA)
library(clusterProfiler)
library(ComplexHeatmap)
library(ggplot2)
library(limma)
library(EnhancedVolcano)
library(pathview)
library(caTools)
library(hgu133plus2.db)
library(circlize)

df = read.csv("TestData1.csv", row.names = 1)
allpheno = read.csv("PhenoData_unsupervised.csv")
allpheno = allpheno %>% dplyr::select(-1,-2)
pheno = read.csv("/Users/khaled/Downloads/chRCC/Unsupervised/PhenoData_unsupervised.csv")

## sample proper
pheno2 = pheno %>% filter(Histology2 == "chRCC" & final_cluster == 4 |
                            Histology2 == "Oncocytoma" & final_cluster == 2)

### Differential Expression
df2 = df[ , match(pheno2$geo_accession, colnames(df))]
design = model.matrix(~ pheno2$Histology2)
fit = lmFit(df2, design)
fit = eBayes(fit)
topGene = topTable(fit, n = dim(df2)[1])

EnhancedVolcano(topGene, rownames(topGene), x = "logFC", y = "adj.P.Val",
                title = "chRCC-Oncocytoma comparison", 
                subtitle = "Differential expression with Limma",
                selectLab = "None")

## AUC
AUC = colAUC(t(df2), pheno2$final_cluster, alg = "ROC")
AUC = AUC %>% t %>% as.data.frame()

SYMBOL = mapIds(hgu133plus2.db, keys = rownames(topGene), "SYMBOL", keytype = "PROBEID")
ENTREZ = mapIds(hgu133plus2.db, keys = rownames(topGene), "ENTREZID", keytype = "PROBEID")

topGene2 = cbind(SYMBOL, ENTREZ, topGene)
topGene2 = merge(topGene2, AUC, by = 0)

topGeneX = topGene2 %>% filter(adj.P.Val <= 0.05)
topGeneX %<>% filter(logFC != 0)
topGeneX %<>% filter(abs(logFC) >= 1)
colnames(topGeneX)[10] = "AUC"
topGeneX %<>% filter(AUC > .9)


chRCC = pheno2[pheno2$Histology2 == "chRCC", ]
Oncocytoma = pheno2[!pheno2$Histology2 == "chRCC", ]

mchRCC = df2 %>% dplyr::select(all_of(chRCC$geo_accession)) %>% 
  apply(1, mean)
qchRCC = df2 %>% dplyr::select(all_of(chRCC$geo_accession)) %>% 
  apply(1, quantile, probs = c(.02, .10, .25, .50, .75, .90, .98)) %>% t

mOnc = df2 %>% dplyr::select(all_of(Oncocytoma$geo_accession)) %>% 
  apply(1, mean)
qOnc = df2 %>% dplyr::select(all_of(Oncocytoma$geo_accession)) %>% 
  apply(1, quantile, probs = c(.02, .10, .25, .50, .75, .90, .98)) %>% t

mnq = cbind.data.frame(mchRCC, qchRCC, mOnc, qOnc)

library(cutpointr)
dt2 = cbind.data.frame(pheno2$final_cluster, t(df2)) 
colnames(dt2)[1] = "Group"
totalGenelist = colnames(dt2)[-1]

cp = c()
# i = 3
for (i in seq_along(totalGenelist)) {
  dt3 = dt2 %>% dplyr::select(c(Group, i+1))
  names(dt3)[2] = "gid" 
  cp1 = cutpointr(dt3, gid, Group, method = maximize_metric, metric = sum_sens_spec)
  cp2 = cbind(totalGenelist[i], cp1[[2]],cp1[[3]], cp1[[4]],cp1[[5]],
              cp1[[6]], cp1[[7]],cp1[[8]],cp1[[9]], cp1[[10]],
              cp1[[11]],cp1[[12]], cp1[[13]])
  cp = rbind(cp2, cp)
  print(i)
}
dev.off()

colnames(cp) = c("Gene", names(cp1[2]),names(cp1[3]), 
                 names(cp1[4]),names(cp1[5]),names(cp1[6]), names(cp1[7]),
                 names(cp1[8]),names(cp1[9]), names(cp1[10]),
                 names(cp1[11]),names(cp1[12]), names(cp1[13]))

mData = merge(cp, topGene2, by.x = "Gene", by.y = "Row.names")
mData = merge(mData, mnq, by.x = "Gene", by.y = 0)

write_csv(mData, "DE_mq_cp.csv")


### Differential with N
pheno3 = allpheno %>% filter(Histology2 == "chRCC" & final_cluster == 1 |
                            Histology2 == "Oncocytoma" & final_cluster == 2 |
                            Histology2 == "N" & final_cluster == 3)
df3 = df[ , colnames(df) %in% pheno3$geo_accession]
# nsamp = pheno3$geo_accession[pheno3$final_cluster == 3]
# csamp = pheno3$geo_accession[pheno3$final_cluster == 1]
# osamp = pheno3$geo_accession[pheno3$final_cluster == 2]
# mean_csamp = apply(df3[ , colnames(df3) %in% csamp], 1, mean)
# mean_nsamp = apply(df3[ , colnames(df3) %in% nsamp], 1, mean)
# mean_osamp = apply(df3[ , colnames(df3) %in% osamp], 1, mean)
# mdata2 = data.frame(row.names = names(mean_osamp), mean_csamp, mean_osamp, mean_nsamp)
# topTable_data = merge(TP_chRCC, TP_Onc, by = 0)
# topTable_data = merge(topTable_data, mdata2, by.x = "Row.names", by.y = 0)
# write_csv(topTable_data, "DE2.csv")


df3 = df[ , match(pheno3$geo_accession, colnames(df))]
design2 = model.matrix(~ 0 + Histology2, data = pheno3)
colnames(design2) = c("chRCC", "N", "Oncocytoma")
fit2 = lmFit(df3, design2)
contrast.matrix = makeContrasts(
  "chRCC" = chRCC - N,
  "Oncocytoma" = Oncocytoma - N,
  levels = design2
)
fit2 = contrasts.fit(fit2, contrast.matrix)
fit2 = eBayes(fit2)
TP_chRCC = topTable(fit2, coef = "chRCC", n = dim(df3)[1], adjust.method = "fdr")
TP_Onc = topTable(fit2, coef = "Oncocytoma", n = dim(df3)[1], adjust.method = "fdr")
TP_chRCC = TP_chRCC[order(rownames(TP_chRCC), decreasing = T), ]
TP_Onc = TP_Onc[order(rownames(TP_Onc), decreasing = T), ]

mData = mData[order(mData$Gene, decreasing = T), ]

differential_data = merge(mData, TP_chRCC, by.x = "Gene", by.y = 0)
differential_data = merge(differential_data, TP_Onc, by.x = "Gene", by.y = 0)
write.csv(differential_data, "DifferentialExpression.csv")

save.image("DE_Data.RData")
### fgsea
library(tidyverse)
library(magrittr)
library(org.Hs.eg.db)
library(hgu133plus2.db)
library(tidymodels)
library(clusterProfiler)
library(ReactomePA)
library(fgsea)
library(pathview)
library(circlize)

# TP_chRCC, TP_Onc
TP_chRCC %<>% arrange(desc(logFC))
TP_Onc %<>% arrange(desc(logFC))
fc_chrcc = TP_chRCC$logFC

chRCC_geneID = mapIds(hgu133plus2.db, rownames(TP_chRCC), "ENTREZID", "PROBEID")
names(fc_chrcc) = chRCC_geneID

x1 = enrichPathway(chRCC_geneID, organism = "human", 
                  pAdjustMethod = "BH", minGSSize = 100, readable = T)
p1 = cnetplot(x1, foldChange = fc_chrcc, categorySize = "pvalue", colorEdge = T, 
              node_label = "all")
p1
barplot(x1, showCategory = 10)

chrcc_network = x1@result
chrcc_network %<>% filter(p.adjust <= 0.05)

chrcc_network %>% filter(Count > 100) %>% arrange(p.adjust) %>% dplyr::slice(1:10) %>% 
  ggplot(aes(x = Count, y = Description, color = p.adjust, size = Count)) +
  geom_point(alpha = 0.5) + theme_classic() +
  theme(text = element_text(size = 18))

## kegg and pathview
GID1 = mapIds(hgu133plus2.db, rownames(TP_chRCC), "ENTREZID", "PROBEID")
fchange1 = TP_chRCC$logFC
names(fchange1) = GID1
fchange1 = fchange1[-which(duplicated(names(fchange1)))]
fchange1 = fchange1[order(fchange1, decreasing = T)]
head(fchange1)
summary(fchange1)
kegg1 = gseKEGG(geneList = fchange1,
                organism = "hsa",
                minGSSize = 50,
                # maxGSSize = 800,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                keyType = "ncbi-geneid")
kres1 = kegg1@result
setwd("./KEGG/chRCC")
pathview(gene.data = fchange1, pathway.id = "hsa00190", species = "hsa")
pathview(gene.data = fchange1, pathway.id = "hsa00010", species = "hsa")
pathview(gene.data = fchange1, pathway.id = "hsa00020", species = "hsa")
pathview(gene.data = fchange1, pathway.id = "hsa04115", species = "hsa") ##
pathview(gene.data = fchange1, pathway.id = "hsa05200", species = "hsa") ## pathways in cancer

### onc kegg analysis
GID2 = mapIds(hgu133plus2.db, rownames(TP_Onc), "ENTREZID", "PROBEID")
fchange2 = TP_Onc$logFC
names(fchange2) = GID2
fchange2 = fchange2[!duplicated(names(fchange2))]
fchange2 = fchange2[order(fchange2, decreasing = T)]
kegg2 = gseKEGG(geneList = fchange2,
                            organism = "hsa",
                            minGSSize = 50,
                            maxGSSize = 800,
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            keyType = "ncbi-geneid", eps = 0)
kres2 = kegg2@result

###
kres1 %<>% add_column(Histology = rep("chRCC", dim(kres1)[1]))
kres2 %<>% add_column(Histology = rep("Oncocytoma", dim(kres2)[1]))
kres = rbind.data.frame( kres1, kres2)
ggplot(data = kres, aes(x = Histology, y = Description, color = NES, size = abs(NES))) +
  geom_point() + theme_classic() +
  theme(text = element_text(size = 12))

ggplot(kres2, aes(x = NES, y = Description, size = abs(NES), color = p.adjust)) +
  geom_point() + theme_minimal() + theme(text = element_text(size = 18))
setwd("~/Downloads/chRCC/KEGG/Onco")
pathview(gene.data = fchange2, pathway.id = "hsa00190", species = "hsa")
pathview(gene.data = fchange2, pathway.id = "hsa00010", species = "hsa")
pathview(gene.data = fchange2, pathway.id = "hsa00020", species = "hsa")
pathview(gene.data = fchange2, pathway.id = "hsa04115", species = "hsa") ##
pathview(gene.data = fchange2, pathway.id = "hsa05200", species = "hsa")

## kres 
kresX = kres[c(1,4, 5, 6, 7, 8, 18, 19, 20, 26, 30, 38, 49, 55,59, 60), ]
ggplot(data = kresX, aes(x = Histology, y = Description, color = NES, size = abs(NES))) +
  geom_point() + theme_classic() +
  theme(text = element_text(size = 15),
        axis.text = element_text(color = "black")) + ylab("") + xlab("")

## canonical pathway
## 

###
phenoX = pheno %>% filter(Histology2 == "chRCC" | Histology2 == "Oncocytoma")
TA = HeatmapAnnotation(Histology = phenoX$Histology2,
                       Unsupervised = ifelse(phenoX$final_cluster == 2, "DBU1", "DBU2"),
                       col = list(Histology = c("chRCC" = "limegreen", "Oncocytoma" = "violet")))
draw(TA)
dfX = df[ , colnames(df) %in% phenoX$geo_accession]

dfX %<>% t %>% as.data.frame() %>% dplyr::select(all_of(topGeneX$Row.names))
dfX = scale(t(dfX), center = T, scale = T)
dfX %>% 
  Heatmap(cluster_columns = T, cluster_rows = T, 
          col = colorRamp2(c(-2,0,2), c("blue", "grey", "red")), 
          row_names_side = "left", row_dend_side = "right", 
          row_names_gp = gpar(cex = .1), 
          column_names_gp = gpar(cex = .2), 
          row_dend_width = unit(0.5, "cm"), 
          clustering_distance_rows = "maximum", 
          clustering_method_rows = "ward.D2", 
          top_annotation = TA,
          # right_annotation = row_ha,
          column_km = 2, #km = 3,
          name = "Differential Genes")

save.image("DE_Data2_7152021.RData")

####
x1@result -> chRCC_Onc_reactome
write.csv(chRCC_Onc_reactome, "ReactomeData_chRCC_Onc.csv")

