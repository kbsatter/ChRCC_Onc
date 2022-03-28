

## alluvial for unsupevised data
rm(list = ls())
library(ggalluvial)
fscore = read.csv("/Users/Khaled/Downloads/Thesis/Final Score.csv")


fscore %>% 
  ggplot( aes(axis1 = Histology, axis3 = k_df.cluster,
              axis4 = DBU.clusters, axis5 = nmf_score, axis2 = cc_kmeans,
              axis6 = best_class)) +
  scale_x_discrete(limits = c("Histology", "CC.Kmeans", "Kmeans",
                              "DBU", "NMF", "COnsensus.HC"))+
  xlab("") +
  geom_alluvium(aes(fill = Histology)) +
  scale_fill_manual(values = c("limegreen","violet")) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() + theme(text = element_text(size = 12)) +
  theme(text = element_text(color = "black", size = 12)) +
  ggtitle("chRCC-Onc",
        "Unsupervised Model Comparison")


## alluvial for supevised data

