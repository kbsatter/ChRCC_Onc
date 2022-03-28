
rm(list = ls())

library(tidyverse)
library(magrittr)
library(umap)
library(e1071)
library(caret)

genExp = read.csv("TestData1.csv", row.names = 1)
pheno = read.csv("phenoData.csv")
mlist = read.csv("/Users/khaled/Downloads/chRCC/GS30_allinformation.csv")

###
df = genExp %>% select(-1)
df = df[ match(mlist$Gene, rownames(df)), ]
pheno2 = pheno %>% filter(Histology2 != "N")
df = df[ , match(pheno2$geo_accession, colnames(df))]
df = t(df) %>% as.data.frame()
Histology = pheno2$Histology2 %>% as.factor()
df2 = cbind.data.frame(Histology, df)


###
train_control <- trainControl(method="repeatedcv", number=10, repeats=3)
index = createDataPartition(Histology, p = 0.8, list = F)
x_train = df2[index, ]

svmtrain = train(Histology ~ .,
                 data = x_train,
                 method = "svmLinear",
                 trControl = train_control,
                 preProcess c("center", "scale"))

svmtrain$results

### svm with other package
svmfit = svm(Histology ~ ., data = x_train, kernel = "linear", scale = F)
summary(svmfit)
plot(svmfit,  x_train, X.1552897_a_at ~ X.212163_at)
