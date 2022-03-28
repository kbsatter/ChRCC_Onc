

rm(list = ls())
library(tidyverse)
library(magrittr)
library(tidymodels)
library(umap)
library(GEOquery)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

### create Dataset
## GSE19982
g1 = getGEO("GSE19982", destdir = getwd())
g2 = getGEO("GSE12090", destdir = getwd())
g3 = getGEO("GSE11151", destdir = getwd())
g4 = getGEO("GSE2109", destdir = getwd())
g5 = getGEO("GSE8271", destdir = getwd())
g6 = getGEO("GSE11024", destdir = getwd())
G19982 = g1$GSE19982_series_matrix.txt.gz

e1 = assayData(G19982)
e1 = e1$exprs
p1 = phenoData(G19982)
p1 = p1@data

G12090 = g2$GSE12090_series_matrix.txt.gz
e2 = assayData(G12090)
e2 = e2$exprs
p2 = phenoData(G12090)
p2 = p2@data

G11151 = g3$GSE11151_series_matrix.txt.gz
e3 = assayData(G11151)
e3 = e3$exprs
p3 = phenoData(G11151)
p3 = p3@data

G2109 = g4$GSE2109_series_matrix.txt.gz
e4 = assayData(G2109)
e4 = e4$exprs
p4 = phenoData(G2109)
p4 = p4@data

G8271 = g5$`GSE8271-GPL4866_series_matrix.txt.gz`
e5 = assayData(G8271)
e5 = e5$exprs
p5 = phenoData(G8271)
p5 = p5@data

G11024 = g6$GSE11024_series_matrix.txt.gz
e6 = assayData(G11024)
e6 = e6$exprs
p6 = phenoData(G11024)
p6 = p6@data

### compile pheodata
p1 %<>% select(geo_accession, `disease state:ch1`)
p2 %<>% select(geo_accession, title)
p3 %<>% select(geo_accession, title)
p4 = p4[grep("Kidney", p4$title), ]

## P4 needs to be formatted
p5 %<>% select(geo_accession, characteristics_ch1)
p6 %<>% select(geo_accession, title, source_name_ch1) 

p4$name <- apply(p4, 1, function(x)as.integer(any(grep("Histology",x))))

p4 %<>% filter(name != 0)
p4.1 = c()
for (i in 1:dim(p4)[1]) {
  gid = p4$geo_accession[i]
  histology = grep("Histology", p4[i, ], value = T)
    mydata = cbind(gid, histology)
  p4.1 = rbind(p4.1, mydata)
  print(i)
}

## fix the column names
colnames(p1)
colnames(p1)[2] = "Histology"
colnames(p2)
colnames(p2)[2] = "Histology"
colnames(p3)
colnames(p3)[2] = "Histology"
colnames(p4.1)
colnames(p4.1) = c("geo_accession", "Histology")
colnames(p5)
colnames(p5)[2] = "Histology"
colnames(p6)
colnames(p6)[2] = "Histology"

p1 %<>% add_column(Batch = "GSE19982")
p2 %<>% add_column(Batch = "GSE12090")
p3 %<>% add_column(Batch = "GSE11151")
p4.1 %<>% as.data.frame %>% add_column(Batch = "GSE2109")
p5 %<>% add_column(Batch = "GSE8271")
p6 %<>% add_column(Batch = "GSE11024")
p6 %<>% select(-3)

phenotype = rbind(p1, p2, p3, p4.1, p5, p6)

table(phenotype$Histology)
Histology2 = vector(mode = "character", length = 495)
Histology2[grep("oncocytoma", phenotype$Histology, ignore.case = T)] = "Oncocytoma"
Histology2[grep("ON_K", phenotype$Histology, ignore.case = T)] = "Oncocytoma"
Histology2[grep("Chromophobe", phenotype$Histology, ignore.case = T)] = "chRCC"
Histology2[grep("CHR_K", phenotype$Histology, ignore.case = T)] = "chRCC"
Histology2[grep("normal", phenotype$Histology, ignore.case = T)] = "N"           
Histology2[grep("NO_K", phenotype$Histology, ignore.case = T)] = "N"           
Histology2[grep("normal", phenotype$Histology, ignore.case = T)] = "N"           

phenotype2 = cbind(phenotype, Histology2)
# write_csv(phenotype2, "P2.csv")

phenotype2 %<>% filter(Histology2 != "")
table(phenotype2$Histology2)
phenotype2 = phenotype2[!duplicated(phenotype2$geo_accession), ]

#### organize the expression data
dim(e1)
e1 = e1[ , colnames(e1) %in% phenotype2$geo_accession]
e2 = e2[ , colnames(e2) %in% phenotype2$geo_accession]
e3 = e3[ , colnames(e3) %in% phenotype2$geo_accession]
e4 = e4[ , colnames(e4) %in% phenotype2$geo_accession]
e5 = e5[ , colnames(e5) %in% phenotype2$geo_accession]
e6 = e6[ , colnames(e6) %in% phenotype2$geo_accession]


### jetset
library(jetset)

### check gene names
gname = function(x)(rownames(x)[1:5])
cbind(gname(e1), gname(e2), gname(e3), gname(e4), gname(e5), gname(e6))


### GSE8271 is loaded separately from local drive
e5 = read.csv("GSE8271.csv", row.names = 1)
cbind(gname(e1), gname(e2), gname(e3), gname(e4), gname(e5), gname(e6))
###
# best probe for 11024
bestProbe = jmap("hgu133plus2", eg = rownames(e5))
bestProbe = bestProbe[!is.na(bestProbe)]
head(bestProbe)
e5 = e5[rownames(e5) %in% names(bestProbe), ]

# best probe for 8271
bestProbe2 = jmap("hgu133plus2", eg = rownames(e6))
bestProbe2 = bestProbe2[!is.na(bestProbe2)]
bestProbe2 = bestProbe2[complete.cases(match(bestProbe2, bestProbe))]
head(bestProbe2)

# best probes for 11151, 19982, 2109, 12090
bestProbe3 = bestProbe2[complete.cases(match(bestProbe2, rownames(e3)))]

"final probe bestProbe3"
e3 = e3[rownames(e3) %in% bestProbe3, ]
e6 = e6[rownames(e6) %in% names(bestProbe3), ]
e5 = e5[ rownames(e5) %in% names(bestProbe3), ]
e1 = e1[ rownames(e1) %in% bestProbe3, ]
e2 = e2[ rownames(e2) %in% bestProbe3, ]
e4 = e4[ rownames(e4) %in% bestProbe3, ]

## transpose
e1 = t(e1) %>% as.data.frame()
e2 = t(e2) %>% as.data.frame()
e3 = t(e3) %>% as.data.frame()
e4 = t(e4) %>% as.data.frame()
e5 = t(e5) %>% as.data.frame()
e6 = t(e6) %>% as.data.frame()

## convert entrez to probid 11024 8271
head(colnames(e6))
vex = bestProbe3[match(colnames(e6), names(bestProbe3))]
head(vex)
colnames(e6) = vex
head(colnames(e6))

vex2 = bestProbe3[match(colnames(e5), names(bestProbe3))]
head(vex2)
colnames(e5) = vex2
###
head(e1)
head(e2)
head(e3) ##
head(e4) ##
head(e5)
head(e6)

e3 %<>% log2 
e4 %<>% log2
df = rbind.data.frame(e1, e3, e4, e5, e6)

write.csv(df, "rawdata.csv")

match(phenotype2$geo_accession, rownames(df))
phenotype3 = phenotype2[match(rownames(df), phenotype2$geo_accession), ]
identical(rownames(df), phenotype3$geo_accession)

library(ggplot2)
library(sva)
library(umap)

Batch = as.factor(phenotype3$Batch)
mod = as.factor(phenotype3$Histology2)

set.seed(158)
identical(rownames(df), phenotype3$geo_accession)
pca1 = prcomp(df, scale. = T)
score1 = pca1$x %>% as.data.frame()
ggplot(data = score1, aes(x = PC1, y = PC2, label = rownames(score1))) + 
  geom_point(aes(color = Batch)) + theme_classic() +
  theme(axis.text = element_text(size = 15), legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

Histology = as.factor(phenotype3$Histology2)
ggplot(data = score1, aes(x = PC1, y = PC2, label = rownames(score1))) + 
  geom_point(aes(color = Histology)) + theme_classic() +
  theme(axis.text = element_text(size = 15), legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  scale_color_manual(values = c("limegreen", "royalblue3", "violet"))
u1 = umap(df)
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

df1 = ComBat(dat = t(df), mod = mod, batch = Batch, par.prior = T, mean.only = F,
             ref.batch = "GSE11024") ####ComBat
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
  scale_color_manual(values = c("limegreen", "royalblue3", "violet"))

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
  scale_color_manual(values = c("limegreen", "royalblue3", "violet"))

##### umap without N
set.seed(142)
phenotype4 = phenotype3 %>% filter(Histology2 != "N")
df3 = df1[ , match(phenotype4$geo_accession, colnames(df1))]
identical(colnames(df3), phenotype4$geo_accession)
u3 = umap(t(df3), 
          n_neighbors = 20,
          min_dist = 0.01,
          metric = "cosine",
          input = "data") 
Histology = phenotype4$Histology2
u3$layout %>% as.data.frame() %>%  
  ggplot(aes(V1,V2, color = Histology)) + 
  geom_point() + theme_classic() +
  theme(axis.text = element_text(size = 15), legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  scale_color_manual(values = c("limegreen", "violet"))

Batch = phenotype4$Batch
u3$layout %>% as.data.frame() %>%  
  ggplot(aes(V1,V2, color = Batch)) + 
  geom_point() + theme_classic() +
  theme(axis.text = element_text(size = 15), legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))
### convert probeID to SYMBOL
library(hgu133plus2.db)
GENEID = select(hgu133plus2.db, rownames(df1), "SYMBOL", "PROBEID")

GENEID2 = GENEID[match(rownames(df1), GENEID$PROBEID), "SYMBOL"]
table(GENEID2)
length(GENEID2)
df2 = cbind.data.frame(GENEID2, df1)

###
write.csv(df2, "TestData1.csv")
write.csv(e2, "GSE12090exp.csv")
write.csv(phenotype3, "PhenoData.csv")
save.image("DataPrep.RData")
