#Plot the reactome results 

#load library
library(tidyverse)

# read the csv

cluster_1_reactome <- read_csv("Results_csv/IMMERSE paper/SCRIPT 2/Reactome results/cluster 1 reactome.csv")
cluster_2_reactome <- read_csv("Results_csv/IMMERSE paper/SCRIPT 2/Reactome results/cluster 2 reactome.csv")
cluster_3_reactome <- read_csv("Results_csv/IMMERSE paper/SCRIPT 2/Reactome results/cluster 3 reactome.csv")
cluster_4_reactome <- read_csv("Results_csv/IMMERSE paper/SCRIPT 2/Reactome results/cluster 4 reactome.csv")

colnames(cluster_1_reactome)

cluster_1_reactome <- cluster_1_reactome %>% head(10) ##s
colnames(cluster_1_reactome)

cluster1 <- ggplot(cluster_1_reactome, aes(x = -log10(`Entities pValue`), y = reorder(`Pathway name`, -`Entities pValue`)))+
  geom_bar(aes(fill = -log10(`Entities pValue`)),stat = "identity", width=0.6, position=position_dodge(), colour="black")+
  theme_classic()+
  scale_fill_viridis_c()+
  ylab("Pathway")+ theme(text = element_text(size = 16))+
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")+
  theme(legend.position = "bottom")

ggsave(plot = cluster1, "FIGURES/IMMERSE paper/SCRIPT 8/cluster1_reactome_plot.png", dpi = 300, height = 5, width = 9)



########################################################
########################################################
########################################################
cluster_2_reactome <- cluster_2_reactome %>% head(10) ## default is 10
colnames(cluster_1_reactome)

cluster2 <- ggplot(cluster_2_reactome, aes(x = -log10(`Entities pValue`), y = reorder(`Pathway name`, -`Entities pValue`)))+
  geom_bar(aes(fill = -log10(`Entities pValue`)),stat = "identity", width=0.6, position=position_dodge(), colour="black")+
  theme_classic()+
  scale_fill_viridis_c()+
  ylab("Pathway")+ theme(text = element_text(size = 16))+
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")+
  theme(legend.position = "bottom")

ggsave(plot = cluster2, "FIGURES/IMMERSE paper/SCRIPT 8/cluster2_reactome_plot.png", dpi = 300, height = 5, width = 9)


########################################################
########################################################
########################################################
#cluster 3
cluster_3_reactome <- cluster_3_reactome %>% head(10) ## default is 10
colnames(cluster_3_reactome)

cluster3 <- ggplot(cluster_3_reactome, aes(x = -log10(`Entities pValue`), y = reorder(`Pathway name`, -`Entities pValue`)))+
  geom_bar(aes(fill = -log10(`Entities pValue`)),stat = "identity", width=0.6, position=position_dodge(), colour="black")+
  theme_classic()+
  scale_fill_viridis_c()+
  ylab("Pathway")+ theme(text = element_text(size = 16))+
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")+
  theme(legend.position = "bottom")

ggsave(plot = cluster3, "FIGURES/IMMERSE paper/SCRIPT 8/cluster3_reactome_plot.png", dpi = 300, height = 5, width = 9)


########################################################
########################################################
########################################################
#cluster 4
cluster_4_reactome <- cluster_4_reactome %>% head(10) ## default is 10
colnames(cluster_4_reactome)

cluster4 <- ggplot(cluster_4_reactome, aes(x = -log10(`Entities pValue`), y = reorder(`Pathway name`, -`Entities pValue`)))+
  geom_bar(aes(fill = -log10(`Entities pValue`)),stat = "identity", width=0.6, position=position_dodge(), colour="black")+
  theme_classic()+
  scale_fill_viridis_c()+
  ylab("Pathway")+ theme(text = element_text(size = 16))+
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")+
  theme(legend.position = "bottom")

ggsave(plot = cluster4, "FIGURES/IMMERSE paper/SCRIPT 8/cluster4_reactome_plot.png", dpi = 300, height = 5, width = 9)




