# script for analysing the  COMBAT data set 

# load packages 
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
df <-  read_csv("Public data sets/SRS COMABAT/Cano_Gamez_2022_Logcpm_864_20416.csv")

#get gene names using biomaRt
library("biomaRt")
#the useEnsembl line worked as had to redirect to useast server. choose between useEnsembl and useMart option
#ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = 'useast')
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl", host = "useast.ensembl.org")
#mart = useMart("ensembl")
genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                  filters = "ensembl_gene_id",
                  values = df$GeneID,
                  mart = ensembl)
idx <- match(df$GeneID, genemap$ensembl_gene_id)
df$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
#df$gene_biotype <- genemap$gene_biotype[ idx ]

colnames(df)

df1 <- df%>%  mutate(hgnc_symbol = na_if(hgnc_symbol, ""))  %>%  mutate(gene_name = coalesce(hgnc_symbol, GeneID)) 
#need to make genes v


#read in LRT genes 
#need to select the 651 gene set i identified
LRT <- read_csv("states gene set/LRT analysis sepsis timecourse.csv")%>% filter(padj <=0.01)%>% 
  filter(log2FoldChange > 1 |log2FoldChange < -1) %>% 
  mutate(gene_id_2 = coalesce(hgnc_symbol,gene_ID))
colnames(LRT)



#make list
LRT_gene_list <- LRT %>% pull(gene_id_2)

#filter the gene set 
df2 <- df1 %>% 
  filter(gene_name %in% LRT_gene_list) %>% 
  as.data.frame() %>% 
  column_to_rownames("gene_name") %>% 
  dplyr::select(-GeneID, - hgnc_symbol)%>% 
  t() %>% 
  as.data.frame()

# Now need to run analysis  


colnames(df2)
set.seed(5648)
tree.phate <- phate(df2, knn= 30, decay=5, knn.dist.method = "euclidean")
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
ggsave(plot = phate_plot, "Figures/IMMERSE paper/IMMERSE paper/Combat/phatePlot_test.png", dpi = 300, height = 5, width = 5)

##############################################################################################################
#Run the clustering 

#determine number of clusters
# Elbow method
elbow <- fviz_nbclust(df2, hcut, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method")
elbow
ggsave(plot = elbow, "Figures/IMMERSE paper/Combat/elbow.png", height = 5, width = 5, dpi  = 300)

#
##############################################################################################################################
# run the clustering 


res.hc3 <- eclust(df2,hc_method = "ward.D2",hc_metric = "euclidean",FUNcluster= "hclust", k	=3, stand = TRUE)


group3 <- as.data.frame(res.hc3$cluster)%>% dplyr::rename("hclust3" ="res.hc3$cluster")

#merge into data frame
df2 <- cbind(df2, group3)

colnames(df2)
dim(df2)

#quick check to see where clusters lie
hclust3_raw_PHATE <- ggplot(df2 , aes(x = PHATE1, y = PHATE2)) +
  geom_point(aes(fill = as.factor(hclust3)), colour = "black", alpha = 1, size = 4, pch = 21)+
  #scale_fill_manual(values = c("SS1" = "#FF9AA2", "SS2" = "#ffffba", "SS3" = "#A7bed3"))+
  theme_bw()+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(aspect.ratio = 1)

##############################################################################################################################
# quick check to work out cluster assignments 
library(viridis)
colnames(df2)
ggplot(df2, aes(x=PHATE1, y=PHATE2))+
  geom_point(aes(fill=`CD3D`),size=4, alpha=1,shape = 21)+
  scale_fill_viridis(discrete = F)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)

df2$hclust3 <- as.factor("")
df2 <- df2 %>% mutate(State = case_when(hclust3 == '1' ~ 'State 3',
                                                                    hclust3 == '2' ~ 'State 2',
                                                                    hclust3 == '3' ~ 'State 1',))


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
  ggtitle("Combat 2022")

ggsave(plot = pie, "Figures/IMMERSE paper/Combat/pie.png", dpi = 300, height = 5, width = 5)

######################################################################################################################################
# use this one to create a pseudotime based on 3 clusters
######################################################################################################################################
##Create a Slingshot object and add the vst assay information to it.e.g remove the meta data and dimension reduction so this is not in assay section


colnames(df4)
# add in the normalised data from deseq
sce <- SingleCellExperiment(assays = List(norm = as.matrix(t(df2))))


#add the reduced dimensions that I want to use here that is the PHATE cordinates. 
reducedDims(sce) <- SimpleList(PHATE = as.matrix(df2 %>%dplyr::select(c("PHATE1", "PHATE2"))))
#here that is the raw data to test the paths- doesnt work as well and cant plot arrows. 
#reducedDims(sce) <- SimpleList(expressiondata = as.matrix(df_sepsis %>% dplyr::select(13:5012)))

#slingshot requires clusters to run analysis. For us this is the clusters as no timepoint in this data

colData(sce)$State <- as.factor(df2$State)

#run slingshot
sce <- slingshot(sce, clusterLabels = 'State', reducedDim = "PHATE", start.clus = 'State 1')
?slingshot
sce$slingshot

######################################################################################################################################

df2_1 <- cbind(df2, as.data.frame(slingPseudotime(sce, na = FALSE), row.names = FALSE))
df2_1 <- df2_1 %>% rownames_to_column("sample_id")
colnames(df2_1)
#write  the csv for future analysis
#write_csv(vst_1_t_meta_sling, "Results_CSV/vst_1_t_meta_sling.csv")
pseudotimevalues1 <- df2_1%>%dplyr::select(c(PHATE1, PHATE2, Lineage1, sample_id, 
                                             hclust3, State))

head(pseudotimevalues1)
pseudotimevalues1 <- melt(pseudotimevalues1, id.vars=c('State', 'PHATE1', 'PHATE2', 
                                                       'sample_id',
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

ggsave(plot = PHATE_pseudotime, "Figures/IMMERSE paper/Combat/PHATE_pseudotime.png", dpi = 300, height = 5, width = 6)

######################################################################################################################################

#PHATE by state 
pseudotimevaluesexclNA1$`SRS group`
PHATE_State <- ggplot(pseudotimevaluesexclNA1, aes(x=PHATE1, y=PHATE2))+
  geom_point(aes(fill=State),size=2, alpha=1,shape = 21)+
  scale_fill_manual(values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)+
  geom_path(data = slingcurves1, aes(x=PHATE1, y=PHATE2), 
            arrow = arrow(type = "open", angle = 40,  length = unit(0.2, "inches")), size =2)

ggsave(plot = PHATE_State, "Figures/IMMERSE paper/Combat/PHATE_State.png", dpi = 300, height = 5, width = 6)

######################################################################################################################################
######################################################################################################################################
#plot the heatmap

colnames(df2_1_hm)
df2_1_hm <- df2_1 %>%arrange(Lineage1)

colnames(df2_1_hm)
df2_1_hm_clus <- df2_1_hm %>% 
  dplyr::select(`sample_id`,2:645) %>% 
  column_to_rownames("sample_id") %>% 
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
colnames(df4_1_hm)
hm_matrix <- df2_1_hm %>% dplyr::select(sample_id,2:645) %>% column_to_rownames("sample_id")%>% scale()
hm_matrix_smooth = gaussianSmooth(hm_matrix, sigma = 10) %>% t()

col_fun = colorRamp2(c(0,0.5, 1), c("blue", "#EEEEEE", "red"))
col_fun(seq(-1, 1))

#make annotation
colnames(df2_1_hm)
anno <- df2_1_hm %>% dplyr::select(Lineage1, State)
colnames(anno)
#create colour scheme for pseudotime´
n = nrow(anno)
min_v = min(anno$Lineage1)
max_v = max(anno$Lineage1)
Var = circlize::colorRamp2(seq(min_v, max_v, length = n), hcl.colors(n,"Spectral"))
#create annotation
ha = HeatmapAnnotation(df = anno, 
                       col = list(Lineage1 = Var,
                                  State = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF")),
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

png("Figures/IMMERSE paper/Combat/heatmap.png",width=16,height=8,units="in",res = 800)
hm
dev.off()

########################################################################







