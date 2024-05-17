# script to look at the orgianal microarray data from MARS paper in terms of my states 

#load packages 
library(tidyverse)
library(limma)
library(biomaRt)
library(uwot)
library(phateR)
library(SingleCellExperiment)
library(slingshot)
library(reshape2)
library(ggalt)
library(ggforce)
library(RColorBrewer)
library(ggbeeswarm)
library(tradeSeq)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)
library(ggrepel)
library(clusterExperiment)
library(factoextra)
library(NbClust)
library(stats)
library(stringr)
library(dbscan)



# read in data 

# Sciluna 2016  data set
colnames(Discovery)
Discovery <- read_csv("Public data sets/Mars/GSE65682_series_matrix.csv")

  

#read in gene name and probe id matching given by Brendon
probenames <- read_csv("Public data sets/COMBAT and MARS merged/Annotations_Bowtie.csv") %>% dplyr:: select(-biotype_gene)
colnames(Discovery)

#merge probenames with gene names

df <- merge(probenames, Discovery )

df <- df %>%  mutate(gene_name = coalesce(hgnc_symbol, probeID))

#filter genes for LRT genes only. 

# read in 651 LRT genes

LRT <- read_csv("states gene set/LRT analysis sepsis timecourse.csv")%>% filter(padj <=0.01)%>% 
  filter(log2FoldChange > 1 |log2FoldChange < -1) %>% 
  mutate(gene_id_2 = coalesce(hgnc_symbol,gene_ID))
colnames(LRT)
#make list
LRT_gene_list <- LRT %>% pull(gene_id_2)

#filter the gene set 
df1 <- df%>% 
  filter(gene_name %in% LRT_gene_list) %>% 
  as.data.frame() 

#remove duplicates as these are isoforms
df1 <- df1[!duplicated(df1$hgnc_symbol), ]
rownames(df1) <- NULL
df1 <- df1 %>% 
  column_to_rownames("gene_name") %>% 
  dplyr::select(3:804) %>% 
  t() %>% 
  as.data.frame()%>%
  rownames_to_column("Sample_geo_accession")



#  Sciluna 2016 metadata
discovery_meta <- read_csv("Public data sets/Mars/GSE65682 metadata.csv")

#merge the metadata and the gene data

df2 <- merge(discovery_meta, df1, by = "Sample_geo_accession")  

#############################################################################################################################

# now need to filter for just the sepsis patients and the data frame is ready for downstream analysis


df2 <- df2 %>% filter(Cohort == "discovery" |
                      Cohort == "validation")  
 

####################################################################################################
####################################################################################################
####################################################################################################


# run the analysis 

#library Phate R

colnames(df2)   
set.seed(5648)
tree.phate <- phate(df2[9:255], knn= 30, decay=40, knn.dist.method = "euclidean")
phate<- as.data.frame(tree.phate)
df2 <- cbind(df2,phate)

colnames(df2)

#plot
phate_plot <- ggplot(df2, aes(x = PHATE1, y = PHATE2))+
  #geom_density_2d(colour = "grey", alpha = 0.5,binwidth = 0.01)+
  geom_point( colour = "black", alpha = 1, size = 2)+
  theme_bw()+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(aspect.ratio = 1)#+
#geom_path(aes(group = case_id), size = 0.3, alpha = 0.5,
#         color="black",
#        arrow = arrow(type = "closed",
#                     length=unit(0.07, "inches"),
#                     ends = "last"),alpha = 0.5)

phate_plot
ggsave(plot = phate_plot, "Figures/IMMERSE paper/Scicluna 2015/phatePlot_test.png", dpi = 300, height = 5, width = 5)

##############################################################################################################
#Run the clustering 

#determine number of clusters
# Elbow method
elbow <- fviz_nbclust(df2[9:230], hcut, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method")
elbow
ggsave(plot = elbow, "Figures/IMMERSE paper/Scicluna 2015/elbow.png", height = 5, width = 5, dpi  = 300)

##############################################################################################################################
# run the clustering 
colnames(df2)
res.hc3 <- eclust(df2[9:255],hc_method = "ward.D2",hc_metric = "euclidean",FUNcluster= "hclust", k	=3, stand = TRUE)


group3 <- as.data.frame(res.hc3$cluster)%>% dplyr::rename("hclust3" ="res.hc3$cluster")

#merge into data frame
df2 <- cbind(df2, group3)

colnames(df2)

hclust3_raw_PHATE <- ggplot(df2 , aes(x = PHATE1, y = PHATE2)) +
  geom_point(aes(fill = as.factor(hclust3)), colour = "black", alpha = 1, size = 4, pch = 21)+
  #scale_fill_manual(values = c("SS1" = "#FF9AA2", "SS2" = "#ffffba", "SS3" = "#A7bed3"))+
  theme_bw()+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(aspect.ratio = 1)

ggsave(plot = hclust3_raw_PHATE, "Figures/IMMERSE paper/clustering information/hclust3_raw_PHATE.png", dpi = 300)

##############################################################################################################################
#test
ggplot(df2, aes(x = PHATE1, y = PHATE2))+
  geom_point(aes(fill = CD3D ), pch = 21)
##############################################################################################################################
class(df4)
df2 <- df2 %>% mutate(State = case_when(hclust3 == '2' ~ 'State 2',
                                        hclust3 == '1' ~ 'State 1',
                                        hclust3 == '3' ~ 'State 3',))

##############################################################################################################################
# calculate proportions of each state
colnames(df2)

props <- df2 %>% 
  group_by(State) %>% # Variable to be transformed
  dplyr::count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

pie <- ggplot(props, aes(x = "", y = perc, fill = State)) +
  geom_col(color = "black") +
  geom_label(aes(label = labels), color = c("black"),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE, size = 8) +
  theme(legend.position = "none")+
  guides(fill="none")+
  scale_fill_manual(values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF"))+
  coord_polar(theta = "y") + 
  theme_void()+
  ggtitle("Scicluna 2015")

ggsave(plot = pie, "Figures/IMMERSE paper/Scicluna 2015/pie.png", dpi = 300, height = 5, width = 5)

######################################################################################################################################
# use this one to create a pseudotime based on 3 clusters
######################################################################################################################################
##Create a Slingshot object and add the vst assay information to it.e.g remove the meta data and dimension reduction so this is not in assay section


colnames(df2)
# add in the normalised data from deseq
sce <- SingleCellExperiment(assays = List(norm = as.matrix(t(df2[9:255]))))


#add the reduced dimensions that I want to use here that is the PHATE cordinates. 
reducedDims(sce) <- SimpleList(PHATE = as.matrix(df2 %>%dplyr::select(c("PHATE1", "PHATE2"))))
#here that is the raw data to test the paths- doesnt work as well and cant plot arrows. 
#reducedDims(sce) <- SimpleList(expressiondata = as.matrix(df_sepsis %>% dplyr::select(13:5012)))

#slingshot requires clusters to run analysis. For us this is the clusters as no timepoint in this data

colData(sce)$State <- as.factor(df2$State)

#run slingshot
sce <- slingshot(sce, clusterLabels = 'State', reducedDim = "PHATE", start.clus = 'State 1')



######################################################################################################################################


df2_1 <- cbind(df2, as.data.frame(slingPseudotime(sce, na = FALSE), row.names = FALSE))
df2_1 <- df2_1 %>% rownames_to_column("sample_id")
colnames(df2_1)
#write  the csv for future analysis
#write_csv(vst_1_t_meta_sling, "Results_CSV/vst_1_t_meta_sling.csv")
pseudotimevalues1 <- df2_1%>%dplyr::select(c(PHATE1, PHATE2, Lineage1, sample_id, Sample_geo_accession,
                                             hclust3, State, Cohort, Gender, Age, `Time to event day 28`, 
                                             `MARs group`, `Mortality day 28`, `Pneumonia diagnoses`))

head(pseudotimevalues1)
pseudotimevalues1 <- melt(pseudotimevalues1, id.vars=c('State', 'PHATE1', 'PHATE2', 'Sample_geo_accession',
                                                       'sample_id', 'Cohort', 'Gender', 'Age', 'Time to event day 28', 
                                                       'MARs group', 'Mortality day 28', 'Pneumonia diagnoses',
                                                       'hclust3'), variable.name='lineage', value.name = 'pseudotime')



#exclude any NA's
pseudotimevaluesexclNA1 <- pseudotimevalues1%>%filter(pseudotime!='NA')


#éxtract the curve object for plotting the curve in GGplot 
slingcurves1 <- as.data.frame(slingCurves(sce)[[1]]$s[slingCurves(sce)[[1]]$ord, ])
#slingcurves2 <- as.data.frame(slingCurves(sce)[[2]]$s[slingCurves(sce)[[2]]$ord, ]) # use this if other lineages
#slingcurves3 <- as.data.frame(slingCurves(sce)[[3]]$s[slingCurves(sce)[[3]]$ord, ]) # use this if other lineages
######################################################################################################################################

#plot

#generate colorpalette
colors <- colorRampPalette(rev(brewer.pal(11, 'Spectral'))[-6])(100)
colors1 <- rev(colors)
PHATE_pseudotime <- ggplot(pseudotimevaluesexclNA1%>%arrange(pseudotime), aes(x=PHATE1, y=PHATE2))+
  geom_point(aes(fill=pseudotime),size=3, alpha=1,shape = 21)+
  scale_fill_gradientn(colours=colors1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)+
  geom_path(data = slingcurves1, aes(x=PHATE1, y=PHATE2), 
            arrow = arrow(type = "open", angle = 40,  length = unit(0.2, "inches")), size =2)

ggsave(plot = PHATE_pseudotime, "Figures/IMMERSE paper/Scicluna 2015/PHATE_pseudotime.png", dpi = 300, height = 5, width = 6)

######################################################################################################################################

#PHATE by MARS 
pseudotimevaluesexclNA1$`SRS group`
PHATE_MARs <- ggplot(pseudotimevaluesexclNA1, aes(x=PHATE1, y=PHATE2))+
  geom_point(aes(fill=`MARs group`),size=4, alpha=1,shape = 21)+
  scale_fill_manual(values = c("Mars1" = "#c7deb2",
                               "Mars2" ="#2e64aa",
                               "Mars3" = "#9dacd5",
                               "Mars4" = "#82ba55"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)+
  geom_path(data = slingcurves1, aes(x=PHATE1, y=PHATE2), 
            arrow = arrow(type = "open", angle = 40,  length = unit(0.2, "inches")), size =2)

ggsave(plot = PHATE_MARs, "Figures/IMMERSE paper/Scicluna 2015/PHATE_MARs.png", dpi = 300, height = 5, width = 6)

#PHATE by State
pseudotimevaluesexclNA1
PHATE_state <- ggplot(pseudotimevaluesexclNA1, aes(x=PHATE1, y=PHATE2))+
  geom_point(aes(fill=State ),size=4, alpha=1,shape = 21)+
  scale_fill_manual(values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)+
  geom_path(data = slingcurves1, aes(x=PHATE1, y=PHATE2), 
            arrow = arrow(type = "open", angle = 40,  length = unit(0.2, "inches")), size =2)

ggsave(plot = PHATE_state, "Figures/IMMERSE paper/Scicluna 2015/PHATE_states.png", dpi = 300, height = 5, width = 6)


#Check there is no batch effect between the discovery and validation


PHATE_Cohort <- ggplot(pseudotimevaluesexclNA1, aes(x=PHATE1, y=PHATE2))+
  geom_point(aes(fill=Cohort),size=4, alpha=1,shape = 21)+
  #scale_fill_manual(values = c())+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)+
  geom_path(data = slingcurves1, aes(x=PHATE1, y=PHATE2), 
            arrow = arrow(type = "open", angle = 40,  length = unit(0.2, "inches")), size =2)

ggsave(plot = PHATE_Cohort, "Figures/IMMERSE paper/Scicluna 2015/PHATE_Cohort.png", dpi = 300, height = 5, width = 6)



######################################################################################################################################
######################################################################################################################################
#plot the heatmap

colnames(df4)
df2_1_hm <- df2_1 %>%arrange(Lineage1)

colnames(df2_1_hm)
df2_1_hm_clus <- df2_1_hm %>% 
  dplyr::select(`Sample_geo_accession`,10:256) %>% 
  column_to_rownames("Sample_geo_accession") %>% 
  scale() %>% 
  t()


d <- dist(df2_1_hm_clus, method = "euclidean") # Euclidean distance matrix.
H.fit <- hclust(d, method="ward.D")
plot(H.fit) # display dendogram
groups <- cutree(H.fit, k=4) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
plot(H.fit)
rect.hclust(H.fit, k=4, border="red") 
group <- as.data.frame(groups) %>% rownames_to_column("gene_name")#extract grouping variable 


# merge gene_biotype_df with gene group
colnames(group)
colnames(gene_biotype_df)

group1 <- group %>%  mutate(`Gene cluster` = case_when(groups == "1" ~ "Cluster 1",
                                                       groups == "2" ~ "Cluster 2",
                                                       groups == "3" ~ "Cluster 3",
                                                       groups == "4" ~ "Cluster 4"))
group2 <- group1 %>% dplyr::select(-groups) %>% column_to_rownames("gene_name")

table(group)
head(group2)
## merge cluster group with dataframe

######################################################################################################################################
library("mmand")
colnames(df2_1_hm)
hm_matrix <- df2_1_hm %>% dplyr::select(`Sample_geo_accession`,10:256) %>% column_to_rownames("Sample_geo_accession")%>% scale()
hm_matrix_smooth = gaussianSmooth(hm_matrix, sigma = 10) %>% t()

col_fun = colorRamp2(c(0,0.5, 1), c("blue", "#EEEEEE", "red"))
col_fun(seq(-1, 1))

#make annotation
colnames(df2_1_hm)
df2_1_hm$`Mortality day 28`
anno <- df2_1_hm %>% dplyr::select(Lineage1, State, `MARs group`, `Mortality day 28`)

anno <- anno %>% mutate(`Mortality day 28` = case_when(`Mortality day 28` == "1" ~ "non-survivor",
                                                       `Mortality day 28` == "0" ~ "survivor"))
colnames(anno)
#create colour scheme for pseudotime´
n = nrow(anno)
min_v = min(anno$Lineage1)
max_v = max(anno$Lineage1)
Var = circlize::colorRamp2(seq(min_v, max_v, length = n), hcl.colors(n,"Spectral"))
#create annotation
ha = HeatmapAnnotation(df = anno, 
                       col = list(Lineage1 = Var,
                                  `MARs group` = c("Mars1" = "#c7deb2",
                                                  "Mars2" ="#2e64aa",
                                                  "Mars3" = "#9dacd5",
                                                  "Mars4" = "#82ba55"),
                                  State = c("State 1" = "#FDE725FF", 
                                            "State 2" = "#29AF7FFF", 
                                            "State 3" = "#453781FF"),
                                  `Mortality day 28` = c("survivor" = "#47abd8", 
                                                         "non-survivor" = "#FF4242")),
                       which = "col",
                       gp = gpar(col = "black"))

df4_1_hm$`28 day survival`
# make the column annotation

ha_row = rowAnnotation(df = group2, 
                       col = list(`Gene cluster` = c("Cluster 1" =  "#FF756D",
                                                     "Cluster 2" = "#F6B4C3", 
                                                     "Cluster 3" = "#AAD6FA",
                                                     "Cluster 4" = "#FCE6A9"),
                                  which = "row"))
class(col_fun)

hm <- Heatmap(hm_matrix_smooth, name = "Z-score",
              cluster_columns = FALSE,
              cluster_rows = H.fit,
              show_heatmap_legend = TRUE,
              #clustering_distance_rows = "euclidean",
              #clustering_method_rows = "ward",
              #bottom_annotation = columnAnnotation(mark=anno),
              # width = ncol(df4)*unit(6, "mm"), 
              # height = nrow(df4)*unit(4, "mm"),
              column_title_gp = gpar(fontsize = 3, fontface = "bold"),
              border = TRUE,
              #row_km = 5,
              #column_km = 3,
              show_column_names = TRUE,
              show_row_names = FALSE,
              #heatmap_legend_param = list(title = "Scale"),
              column_names_gp = grid::gpar(fontsize = 4),
              #row_names_gp = grid::gpar(fontsize = 5),
              #right_annotation = ha,
              row_dend_width = unit(30, "mm"),
              top_annotation = ha,
              right_annotation = ha_row)#+
hm

png("Figures/IMMERSE paper/Scicluna 2015/heatmap.png",width=14,height=8,units="in",res = 800)
hm
dev.off()

########################################################################



