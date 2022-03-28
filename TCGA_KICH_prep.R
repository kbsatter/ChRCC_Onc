

rm(list = ls())
library(tidyverse)
library(umap)
library(magrittr)
library(ggplot2)
setwd("./TCGA validation")
df = read.csv("TCGA-KICH.htseq_fpkm-uq.tsv.gz", sep = '\t', row.names = 1)
pheno = read.csv("TCGA-KICH.GDC_phenotype.tsv.gz", sep = '\t', row.names = 1)
genecode = read.csv("gencode.gene.info.v22.tsv", sep = '\t')
table(duplicated(genecode$gene_name))
GeneID = genecode[match(rownames(df), genecode$gene_id), "gene_name"]
df2 = cbind.data.frame(GeneID, df)
df2 = df2[!duplicated(df2$GeneID), ]
rownames(df2) = df2$GeneID
df2 %<>% select(-1)
class(df2)
df2 = t(df2) %>% as.data.frame()

pheno2 = pheno %>% filter(sample_type.samples != "Solid Tissue Normal")
rownames(pheno2) = gsub("-", ".", rownames(pheno2))
match(rownames(df2), rownames(pheno2))
pheno2 = pheno2[rownames(pheno2) %in% rownames(df2), ]
df2 = df2[rownames(df2) %in% rownames(pheno2), ]
pheno2 = pheno2[match(rownames(df2), rownames(pheno2)), ]

write.csv(df2, "TCGA_valData.csv")
write.csv(pheno2, "TCGA_valPheno.csv")
###
# get the genelist
mlist = read.csv("/Users/khaled/Downloads/chRCC/GS30_allinformation.csv")
match(mlist$SYMBOL, colnames(df2))
mlist$SYMBOL[6] = "WDR63"

df3 = df2[ , colnames(df2) %in% mlist$SYMBOL]
colnames(df3)[2] = "DNAI3"


## load the model
write.csv(df3, "TCGA_valData.csv")
write.csv(pheno2, "TCGA_valPheno.csv")
