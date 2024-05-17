#plot the time course analysis 
# Boxplot of gene expressiong for all comparisons. 

# load packages 
library(plyr)
library(tidyverse)
library(ggbeeswarm)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(scales)


########################################################################################################################################
########################################################################################################################################
# PLot the heatmap of genes that are changing with p = 0.01 and and logfold greater 1 or less than 1 differential expressed genes from LRT analysis 

# make list of the significant genes required for the analysis 
#import the LRT timecourse analysis 

LRT <- read_csv("Results_csv/IMMERSE paper/SCRIPT 1/LRT analysis sepsis timecourse.csv")%>% filter(padj <=0.01)%>% 
  filter(log2FoldChange > 1 |log2FoldChange < -1) %>% 
  mutate(gene_id_2 = coalesce(hgnc_symbol,gene_ID))
colnames(LRT)
#make list
LRT_gene_list <- LRT %>% pull(gene_ID)


# import the transformed normalised VSD df
df <- read_csv("Results_CSV/IMMERSE paper/SCRIPT 1/vsd_df.csv") %>% dplyr::select(gene_ID, 2:259) 
colnames(df)

#filter genes of interest
df <- df %>% 
  filter(gene_ID %in% LRT_gene_list) %>% as.data.frame()


#select variables from LRT to merge so can get real gene names
LRT1 <- LRT %>% dplyr::select(gene_ID, gene_id_2) %>% as.data.frame()
#merge with the df
df <- merge(LRT1,df )
colnames(df)
df <- df %>% dplyr::select(2:260)

#need to make genes variables and patients observations
df1 = setNames(data.frame(t(df[,-1])), df[,1]) #transpose so obvs = patients, variables = genes

#now have a data frame with 5000 of the most variable genes in the dataset.
df1 <- df1 %>% rownames_to_column("sample_id")


metadata <- read_csv("Metadata/metadata.csv") %>% 
  dplyr::select("Cohort", "sample_id", "case_id", 
                "Cohort_time", "Timepoint", "COVID_status", "Library_cohort") %>% #choose the variables I need
  filter(COVID_status == "Negative") %>%  #remove covid patients
  filter(!Library_cohort == "Cardiac_lib2")%>%#remove duplicate samples
  #filter(Cohort == "Sepsis") %>% #keep only sepsis
  filter(!case_id == "SC09", #removed as IVIG
         !case_id == "SC84", #removed as insisted on it 
         !case_id == "SC27") %>%  #removed as screen failure, sofa >2
  as.data.frame()# %>% 


#import the SRS groupings for davenport and extended
srs_dp <- read_csv("Results_csv/IMMERSE paper/SCRIPT 1/srs_groupings_davenport.csv") %>% as.data.frame()
srs_ex <- read_csv("Results_csv/IMMERSE paper/SCRIPT 1/srs_groupings_extended.csv") %>% as.data.frame()
#merge with metadata
metadata <- merge(metadata, srs_dp, by = "sample_id", all.x = TRUE)
metadata <- merge(metadata, srs_ex, by = "sample_id", all.x = TRUE)

#read in the mars groupings
mars <- read_csv("Results_csv/IMMERSE paper/SCRIPT 1/mars_groupings.csv")
metadata <- merge(metadata, mars, by = "sample_id", all.x = TRUE)
colnames(metadata)



#create data frame with all metadata and gene expression
df1 <- merge(metadata, df1, all.x = TRUE)

#Now need to create a data frame with averaged counts by Cohort_time
colnames(df1)

table(df1$Cohort_time)

colnames(df_average)
df_average <- df1 %>% dplyr::select(Cohort_time, 13:730) %>% 
  gather("gene", "expression", 2:719)%>% 
  group_by(Cohort_time, gene)%>% 
  summarise(median = median(expression))%>% 
  ungroup() %>% 
  spread(gene, median) %>% 
  column_to_rownames("Cohort_time") %>% 
  scale() %>%
  as.data.frame() %>%  
  rownames_to_column("Cohort_time") %>% 
  as.matrix()



#######
class(df_average)

df_average1 = setNames(data.frame(t(df_average[,-1])), df_average[,1])

df_average1[1:6] <- sapply(df_average1[1:6],as.numeric)


df_average1 <- df_average1 %>%  as.matrix()


# run clustering

d <- dist(df_average1, method = "euclidean") # Euclidean distance matrix.
H.fit <- hclust(d, method="ward.D")
plot(H.fit) # display dendogram
groups <- cutree(H.fit, k=4) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
plot(H.fit)
rect.hclust(H.fit, k=4, border="red") 
group <- as.data.frame(groups)#extract grouping variable 
table(group)
head(group)

#make annotation
colnames(df_average1)
df_average1.0 <- df_average1 %>% as.data.frame() %>% rownames_to_column("gene_id_2")
colnames(LRT)
colnames(df_average1.0)
df_average1.0 <- merge(df_average1.0, LRT)
df_average1.0 <- cbind(df_average1.0, group)

colnames(df_average1.0)
df_average1.0$gene_id_2
class(df_average1.0)

anno <- df_average1.0  %>% dplyr::select(all_of(c("gene_id_2", "gene_biotype", "groups")))%>% 
  mutate(`Gene cluster` = case_when(groups == "1" ~ "Cluster 3",
                                    groups == "2" ~ "Cluster 2",
                                    groups == "3" ~ "Cluster 1",
                                    groups == "4" ~ "Cluster 4")) %>% 
dplyr::select(-gene_id_2, -groups)
table(anno$`Gene cluster`)
library(ComplexHeatmap)
#create annotation
ha = HeatmapAnnotation(df = anno, 
                       col = list(gene_biotype = c("artifact" = "#E31A1C", 
                                                   "IG_C_gene" = "#6A3D9A",
                                                   "IG_C_pseudogene" = "#FF7F00",
                                                   "IG_V_gene" = "gold1",
                                                   "lncRNA" = "skyblue2", 
                                                   "misc_RNA" = "#FB9A99",
                                                   "processed_pseudogene" = "palegreen2",
                                                   "protein_coding" = "#CAB2D6",
                                                   "TEC" = "#FDBF6F",
                                                   "TR_C_gene" = "blue1",
                                                   "TR_V_gene" = "deeppink1",
                                                   "TR_V_pseudogene" = "maroon",
                                                   "transcribed_processed_pseudogene" = "orchid1",
                                                   "transcribed_unitary_pseudogene" = "khaki2",
                                                   "transcribed_unprocessed_pseudogene" = "gray70",
                                                   "unprocessed_pseudogene" = "steelblue4"),
                                  `Gene cluster`= c("Cluster 1" = "#FFB3B3",
                                                    "Cluster 2" = "#C1EFFF",
                                                    "Cluster 3" = "#FFE9AE",
                                                    "Cluster 4" = "#9ED2C6")),
                       which = "row",
                       gp = gpar(col = "black", lty = 0.5))
table(anno$gene_biotype)



#make a heatmap
#set the split for the columns
column_split = rep("Cardiac", 6)
column_split[3:6] = "Sepsis"
library(ComplexHeatmap)
 hm <- Heatmap(df_average1, name = "Scale",
              cluster_columns = FALSE,
              cluster_rows = H.fit,
              show_heatmap_legend = TRUE,
              #clustering_distance_rows = "euclidean",
              #clustering_method_rows = "ward",
              #bottom_annotation = columnAnnotation(mark=anno),
              # width = ncol(df4)*unit(6, "mm"), 
              # height = nrow(df4)*unit(4, "mm"),
              column_title_gp = gpar(fontsize = 12, fontface = "bold"),
              border = TRUE,
              #row_km = 5,
              #column_km = 3,
              show_column_names = TRUE,
              show_row_names = FALSE,
              #heatmap_legend_param = list(title = "Scale"),
              column_names_gp = grid::gpar(fontsize = 12),
              row_names_gp = grid::gpar(fontsize = 5),
              right_annotation = ha,
              row_dend_width = unit(30, "mm"),
             # top_annotation = ha,
              #right_annotation = ha_row
         row_split = 4,
         column_split = column_split)
        



png("FIGURES/IMMERSE paper/SCRIPT 2/HEATMAP.png",width=8,height=11,units="in",res = 800)
hm
dev.off()


################################################################################
################################################################################
################################################################################

#make line plots to compliment heatmap

group <- group %>% rownames_to_column("gene_id_2")
head(group)
#df_average1.0 <- merge(df_average1.0,group)
colnames(df_average1.0)
df_average1.02 <- df_average1.0 %>% gather(Cohort_time, Expression, 2:7)
df_average1.02$Cohort <-substr(df_average1.02$Cohort_time, 0,1) 
colnames(df_average1.02)
df_average1.03 <- df_average1.02%>% 
  mutate(Cohort = case_when(Cohort == "S" ~ "Sepsis",
                            Cohort == "C" ~ "Cardiac"))%>% 
  mutate(`Gene cluster` = case_when(groups == "1" ~ "Cluster 3",
                            groups == "2" ~ "Cluster 2",
                            groups == "3" ~ "Cluster 1",
                            groups == "4" ~ "Cluster 4"))


# make colour scheme
colours <- c("artifact" = "#E31A1C", 
             "IG_C_gene" = "#6A3D9A",
             "IG_C_pseudogene" = "#FF7F00",
             "IG_V_gene" = "gold1",
             "lncRNA" = "skyblue2", 
             "misc_RNA" = "#FB9A99",
             "processed_pseudogene" = "palegreen2",
             "protein_coding" = "#CAB2D6",
             "TEC" = "#FDBF6F",
             "TR_C_gene" = "blue1",
             "TR_V_gene" = "deeppink1",
             "TR_V_pseudogene" = "maroon",
             "transcribed_processed_pseudogene" = "orchid1",
             "transcribed_unitary_pseudogene" = "khaki2",
             "transcribed_unprocessed_pseudogene" = "gray70",
             "unprocessed_pseudogene" = "steelblue4")

#make info for standard error bars
df.summary <- df_average1.03 %>%
  group_by(Cohort_time, `Gene cluster`) %>%
  summarise(
    sd = sd(Expression, na.rm = TRUE),
    Expression = median(Expression)
  )
df.summary
df.summary$Cohort <-substr(df.summary$Cohort_time, 0,1) 
df.summary <- df.summary%>% 
  mutate(Cohort = case_when(Cohort == "S" ~ "Sepsis",
                            Cohort == "C" ~ "Cardiac"))
#plot

df_average1.03$gene_id_2
c1list <- c("ABCD2","ADAMTS7P1")
p1 <- ggplot(df_average1.03 %>% filter(`Gene cluster` == "Cluster 1"), 
       aes(Cohort_time, Expression)) +
  geom_line(aes(group = gene_id_2, color = gene_biotype), alpha = 0.6, size = 0.2) +
#  geom_label_repel(data = df_average1.03 %>% filter(gene_id_2 %in% c1list),aes(label = gene_id_2),
 #                  nudge_x = 1,
  #                 na.rm = TRUE)+
  geom_line(stat = "summary", fun = "median", colour = "black", size = 1.5, 
            aes(group = 1)) +
  geom_errorbar(data = df.summary %>% filter(`Gene cluster` == "Cluster 1"),
                aes(ymin=Expression-sd, ymax=Expression+sd, group = Cohort ), width=.2,
                position=position_dodge(0.05))+
  geom_point(data = df.summary %>% filter(`Gene cluster` == "Cluster 1"), size=4, shape=21, fill = "#FFB3B3")+
  facet_grid(~Cohort, scales = "free_x")+
  scale_colour_manual(values = colours)+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme(legend.position = "none")

ggsave(plot = p1, "FIGURES/IMMERSE paper/SCRIPT 2/line plot cluster 1.png", width= 10, height = 6, dpi = 300)

p2 <- ggplot(df_average1.03 %>% filter(`Gene cluster` == "Cluster 2"),
       aes(Cohort_time, Expression)) +
  geom_line(aes(group = gene_id_2, color = gene_biotype), alpha = 0.6, size = 0.2) +
  #  geom_label_repel(data = df_average1.03 %>% filter(gene_id_2 %in% c1list),aes(label = gene_id_2),
  #                  nudge_x = 1,
  #                 na.rm = TRUE)+
  geom_line(stat = "summary", fun = "median", colour = "black", size = 1.5, 
            aes(group = 1)) +
  geom_errorbar(data = df.summary %>% filter(`Gene cluster` == "Cluster 2"),
                aes(ymin=Expression-sd, ymax=Expression+sd, group = Cohort ), width=.2,
                position=position_dodge(0.05))+
  geom_point(data = df.summary %>% filter(`Gene cluster` == "Cluster 2"), size=4, shape=21, fill = "#C1EFFF")+
  facet_grid(~Cohort, scales = "free_x")+
  scale_colour_manual(values = colours)+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme(legend.position = "none")

ggsave(plot = p2, "FIGURES/IMMERSE paper/SCRIPT 2/line plot cluster 2.png", width= 10, height = 6, dpi = 300)


p3 <- ggplot(df_average1.03 %>% filter(`Gene cluster` == "Cluster 3"),
       aes(Cohort_time, Expression)) +
  geom_line(aes(group = gene_id_2, color = gene_biotype), alpha = 0.6, size = 0.2) +
  #  geom_label_repel(data = df_average1.03 %>% filter(gene_id_2 %in% c1list),aes(label = gene_id_2),
  #                  nudge_x = 1,
  #                 na.rm = TRUE)+
  geom_line(stat = "summary", fun = "median", colour = "black", size = 1.5, 
            aes(group = 1)) +
  geom_errorbar(data = df.summary %>% filter(`Gene cluster` == "Cluster 3"),
                aes(ymin=Expression-sd, ymax=Expression+sd, group = Cohort ), width=.2,
                position=position_dodge(0.05))+
  geom_point(data = df.summary %>% filter(`Gene cluster` == "Cluster 3"), size=4, shape=21, fill = "#FFE9AE")+
  facet_grid(~Cohort, scales = "free_x")+
  scale_colour_manual(values = colours)+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme(legend.position = "none")


ggsave(plot = p3, "FIGURES/IMMERSE paper/SCRIPT 2/line plot cluster 3.png", width= 10, height = 6, dpi = 300)


p4 <- ggplot(df_average1.03 %>% filter(`Gene cluster` == "Cluster 4"),
       aes(Cohort_time, Expression)) +
  geom_line(aes(group = gene_id_2, color = gene_biotype), alpha = 0.6, size = 0.2) +
  #  geom_label_repel(data = df_average1.03 %>% filter(gene_id_2 %in% c1list),aes(label = gene_id_2),
  #                  nudge_x = 1,
  #                 na.rm = TRUE)+
  geom_line(stat = "summary", fun = "median", colour = "black", size = 1.5, 
            aes(group = 1)) +
  geom_errorbar(data = df.summary %>% filter(`Gene cluster` == "Cluster 4"),
                aes(ymin=Expression-sd, ymax=Expression+sd, group = Cohort ), width=.2,
                position=position_dodge(0.05))+
  geom_point(data = df.summary %>% filter(`Gene cluster` == "Cluster 4"), size=4, shape=21, fill = "#9ED2C6")+
  facet_grid(~Cohort, scales = "free_x")+
  scale_colour_manual(values = colours)+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme(legend.position = "none")

ggsave(plot = p4, "FIGURES/IMMERSE paper/SCRIPT 2/line plot cluster 4.png", width= 10, height = 6, dpi = 300)


p5 <- ggplot(df_average1.03 %>% filter(`Gene cluster` == "Cluster 1"), aes(x = `Gene cluster`, fill = gene_biotype)) +
  geom_bar(position = "fill")+
  scale_fill_manual(values = colours)+
  theme_classic()+
  theme(legend.position = "none")+
  ylab("")+
  xlab("")+
  theme(axis.text=element_text(size=12, face = "bold", colour = "black"))

ggsave(plot = p5, "FIGURES/IMMERSE paper/SCRIPT 2/barplot_gene_biotype_cluster1.png", width= 3, height = 6, dpi = 300)

p6 <- ggplot(df_average1.03 %>% filter(`Gene cluster` == "Cluster 2"), aes(x = `Gene cluster`, fill = gene_biotype)) +
  geom_bar(position = "fill")+
  scale_fill_manual(values = colours)+
  theme_classic()+
  theme(legend.position = "none")+
  ylab("")+
  xlab("")+
  theme(axis.text=element_text(size=12, face = "bold", colour = "black"))

ggsave(plot = p6, "FIGURES/IMMERSE paper/SCRIPT 2/barplot_gene_biotype_cluster2.png", width= 3, height = 6, dpi = 300)

p7 <- ggplot(df_average1.03 %>% filter(`Gene cluster` == "Cluster 3"), aes(x = `Gene cluster`, fill = gene_biotype)) +
  geom_bar(position = "fill")+
  scale_fill_manual(values = colours)+
  theme_classic()+
  theme(legend.position = "none")+
  ylab("")+
  xlab("")+
  theme(axis.text=element_text(size=12, face = "bold", colour = "black"))

ggsave(plot = p7, "FIGURES/IMMERSE paper/SCRIPT 2/barplot_gene_biotype_cluster3.png", width= 3, height = 6, dpi = 300)

p8 <- ggplot(df_average1.03 %>% filter(`Gene cluster` == "Cluster 4"), aes(x = `Gene cluster`, fill = gene_biotype)) +
  geom_bar(position = "fill")+
  scale_fill_manual(values = colours)+
  theme_classic()+
  theme(legend.position = "none")+
  ylab("")+
  xlab("")+
  theme(axis.text=element_text(size=12, face = "bold", colour = "black"))

ggsave(plot = p8, "FIGURES/IMMERSE paper/SCRIPT 2/barplot_gene_biotype_cluster4.png", width= 3, height = 6, dpi = 300)

cluster1 <- df_average1.03 %>% filter(`Gene cluster` == "Cluster 1") %>% filter(Cohort_time == "Sepsis_T0") #use this to get one list
write_csv(cluster1, "Results_csv/IMMERSE paper/SCRIPT 2/longitudinal analysis cluster gene lists/cluster1.csv")
cluster2 <- df_average1.03 %>% filter(`Gene cluster` == "Cluster 2")%>% filter(Cohort_time == "Sepsis_T0") #use this to get one list
write_csv(cluster2, "Results_csv/IMMERSE paper/SCRIPT 2/longitudinal analysis cluster gene lists/cluster2.csv")
cluster3 <- df_average1.03 %>% filter(`Gene cluster` == "Cluster 3")%>% filter(Cohort_time == "Sepsis_T0") #use this to get one list
write_csv(cluster3, "Results_csv/IMMERSE paper/SCRIPT 2/longitudinal analysis cluster gene lists/cluster3.csv")
cluster4 <- df_average1.03 %>% filter(`Gene cluster` == "Cluster 4")%>% filter(Cohort_time == "Sepsis_T0") #use this to get one list
write_csv(cluster4, "Results_csv/IMMERSE paper/SCRIPT 2/longitudinal analysis cluster gene lists/cluster4.csv")


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

# exploring the use of the longitudinal analysis significant gene set in PCA analysis. 

# use df that was created already. 

# 1) PCA
#runPCA
colnames(df1)


pca <- prcomp(df1[13:728], scale. = F)
pca1 <- as.data.frame(pca$x)
pca1

df3.1 <- cbind(df1, pca1[1:4])
df3.1$Cohort_time


mu <- ddply(df3.1, "Cohort_time", summarise, grp.median=median(PC1))
mu$Cohort <- substr(mu$Cohort_time, 0,1)
mu <- mu %>%  mutate(Cohort = case_when(Cohort == "C" ~ "Cardiac",
                                        Cohort == "S" ~ "Sepsis"))
#df4$Cohort <- factor(df1$Cluster,
df3.1$Cohort <- factor(df3.1$Cohort, levels = c("Sepsis", "Cardiac"))
mu$Cohort <- factor(mu$Cohort, levels = c("Sepsis", "Cardiac"))
#levels = c('Phenotype 1','Phenotype 2', "Phenotype 3"),ordered = TRUE)
PC1_histo <- ggplot(df3.1, aes(x=PC1, fill = Cohort_time))+ 
  #geom_histogram(aes(y=..density..), colour="black", fill = "white")+
  geom_density(alpha=.7) +
  scale_fill_manual(values = c("Sepsis_T0" = "#efe350ff",
                               "Sepsis_T1" = "#f68f46ff",
                               "Sepsis_T2" = "#a65c85ff",
                               "Sepsis_T3" = "#403891ff",
                               "Cardiac_T0" ="#56A8CBFF",
                               "Cardiac_T1" = "#DA291CFF"))+
  theme_classic()+
  theme(legend.position = "none")+
  geom_vline(data=mu, aes(xintercept=grp.median, color=Cohort_time, group = Cohort),
             linetype="dashed", size =2)+
  scale_color_manual(values = c("Sepsis_T0" = "#efe350ff",
                                "Sepsis_T1" = "#f68f46ff",
                                "Sepsis_T2" = "#a65c85ff",
                                "Sepsis_T3" = "#403891ff",
                                "Cardiac_T0" ="#56A8CBFF",
                                "Cardiac_T1" = "#DA291CFF"))+
  xlab("")+
  xlim(c(-40,45))+
  theme(legend.position = "none")+
  facet_wrap(~Cohort, ncol = 1)+ theme(text = element_text(size = 12))    
#save the PCA histo
ggsave(plot = PC1_histo, "FIGURES/IMMERSE paper/SCRIPT 2/histo_PCA1.png", width = 5, height = 5, dpi = 300)


mu2 <- ddply(df3.1, "Cohort_time", summarise, grp.median=median(PC2))
mu2$Cohort <- substr(mu2$Cohort_time, 0,1)
mu2 <- mu2 %>%  mutate(Cohort = case_when(Cohort == "C" ~ "Cardiac",
                                          Cohort == "S" ~ "Sepsis"))



PC2_histo <- ggplot(df3.1, aes(x=PC2, fill = Cohort_time)) +
  #geom_histogram(aes(y=..density..), colour="black", fill = "white")+
  geom_density(alpha=.7) +
  scale_fill_manual(values = c("Sepsis_T0" = "#efe350ff",
                               "Sepsis_T1" = "#f68f46ff",
                               "Sepsis_T2" = "#a65c85ff",
                               "Sepsis_T3" = "#403891ff",
                               "Cardiac_T0" ="#56A8CBFF",
                               "Cardiac_T1" = "#DA291CFF"))+
  theme_classic()+
  theme(legend.position = "none")+
  geom_vline(data=mu2, aes(xintercept=grp.median, color=Cohort_time, group = Cohort),
             linetype="dashed", size =2)+
  scale_color_manual(values = c("Sepsis_T0" = "#efe350ff",
                                "Sepsis_T1" = "#f68f46ff",
                                "Sepsis_T2" = "#a65c85ff",
                                "Sepsis_T3" = "#403891ff",
                                "Cardiac_T0" ="#56A8CBFF",
                                "Cardiac_T1" = "#DA291CFF"))+
  coord_flip()+
  xlab("")+
  xlim(c(-20,25))+
  theme(axis.text.y=element_blank(), #remove x axis labels
        axis.ticks.y=element_blank(), #remove x axis ticks
  )+
  theme(legend.position  = "none")+
  facet_wrap(~Cohort)+ theme(text = element_text(size = 12))    
ggsave(plot = PC2_histo, "FIGURES/IMMERSE paper/SCRIPT 2/histo_PCA2.png", width = 5, height = 5, dpi = 300)

#figure out proportion of variation for x and y axis on PCA plot
library(factoextra)
library(scales)
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 50))
devtools::install_github("jtlandis/ggside")
library(ggside)
pca_plot <- ggplot(df3.1, aes(x = PC1, y = PC2))+
  #geom_density_2d(colour = "black", alpha = 0.2)+
  geom_point(aes(fill = Cohort_time, shape = Cohort), alpha = 0.9, size = 4, stroke = 0.7, colour = "black")+
  theme_bw()+ 
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_fill_manual(values = c("Sepsis_T0" = "#efe350ff",
                               "Sepsis_T1" = "#f68f46ff",
                               "Sepsis_T2" = "#a65c85ff",
                               "Sepsis_T3" = "#403891ff",
                               "Cardiac_T0" ="#56A8CBFF",
                               "Cardiac_T1" = "#DA291CFF"))+
  scale_shape_manual(values = c("Sepsis" = 21, "Cardiac" = 24))+
  #theme(aspect.ratio = 1)+
  guides(fill=guide_legend(title="Cohort"))+
  ylab("PC2 (13.4%)")+
  xlab("PC1 (43%)")+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  ylim(c(-20,25))+
  xlim(c(-40,45))+
  theme(legend.position = "none",
        legend.box = "vertical",
        legend.title=element_text(size=10, face = "bold"), 
        legend.text=element_text(size=10),
        axis.title =element_text(size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
       panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(plot = pca_plot, "Figures/IMMERSE paper/SCRIPT 2/pca_plot_nolegend.png", dpi = 300, height = 5, width = 5)


# calculate centroids of Cohort_time in the PCA

colnames(df3.1)
PCA_centroids <- aggregate(df3.1[,731:732], list(Type = df3.1$Cohort_time), median)
PCA_centroids$Cohort <- c("Cardiac", "Cardiac", "Sepsis", "Sepsis", "Sepsis", "Sepsis")

PCA_centroids$Type <- factor(PCA_centroids$Type, levels = c("Cardiac_T1", "Cardiac_T0", "Sepsis_T0", "Sepsis_T1",
                                                            "Sepsis_T2", "Sepsis_T3"))

PCA_centroids_plot <- ggplot(PCA_centroids, aes(x = PC1, y = PC2))+
  geom_density_2d(data = df3.1, aes(x = PC1, y = PC2), colour = "grey", alpha = 0.5)+
  geom_point(aes(fill = Type, shape = Cohort), colour = "black", alpha = 0.8, size = 5, stroke = 0.7)+
  scale_fill_manual(values = c("Sepsis_T0" = "#efe350ff",
                               "Sepsis_T1" = "#f68f46ff",
                               "Sepsis_T2" = "#a65c85ff",
                               "Sepsis_T3" = "#403891ff",
                               "Cardiac_T0" ="#56A8CBFF",
                               "Cardiac_T1" = "#DA291CFF"))+
  scale_colour_manual(values = c("Sepsis_T0" = "#efe350ff",
                                 "Sepsis_T1" = "#f68f46ff",
                                 "Sepsis_T2" = "#a65c85ff",
                                 "Sepsis_T3" = "#403891ff",
                                 "Cardiac_T0" ="#56A8CBFF",
                                 "Cardiac_T1" = "#DA291CFF"))+
  scale_shape_manual(values = c("Sepsis" = 21, "Cardiac" = 24))+
  theme_bw()+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom")+
  ylim(c(-20,25))+
  xlim(c(-40,45))+
  theme(legend.position = "none")+
  geom_line(data = PCA_centroids %>% filter(Cohort == "Cardiac"),aes(group = Cohort),
            color="black",
            arrow = arrow(type = "closed",
                          length=unit(0.15, "inches"),
                          ends = "last"))+
  geom_path(data = PCA_centroids %>% filter(Cohort == "Sepsis"),aes(group = Cohort),
            color="black",
            arrow = arrow(type = "closed",
                          length=unit(0.15, "inches"),
                          ends = "last"))+
  theme(aspect.ratio = 1)+
  ylab("PC2 (13.4%)")+
  xlab("PC1 (43%)")+ theme(text = element_text(size = 12))    

PCA_centroids_plot

ggsave(plot = PCA_centroids_plot, "FIGURES/IMMERSE paper/SCRIPT 2/PCA_centroids_plot.png", dpi=300, height = 5, width = 5, units = "in" )








    