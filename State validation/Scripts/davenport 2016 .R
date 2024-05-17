# script to look at the orgianal microarray data from SRS paper in terms of my states 

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

# Davenport 2016 discovery data set

Discovery <- read_csv("Public data sets/Davenport 2016/E-MTAB-4421 discovery/Davenport_sepsis_Jan2016_normalised_265.csv") %>% 
  column_to_rownames("Probe_ID") %>% 
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("Comment[beadchiparray_id]")


# discovery metadata
discovery_meta <- read_csv("Public data sets/Davenport 2016/E-MTAB-4421 discovery/E-MTAB-4421.sdrf.csv")

#merge the metadata and the gene data

Discovery1 <- merge(discovery_meta, Discovery, all.y = T, by = "Comment[beadchiparray_id]")  %>% 
  dplyr::select(-`Comment[beadchiparray_id]`)



################################################################################################    
#### read in the validation dataset 

validation <- read_csv("Public data sets/Davenport 2016/E-MTAB-4451 Validation/Davenport_sepsis_Feb2016_normalised_106.csv")%>% 
  column_to_rownames("Probe_ID") %>% 
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("Source Name")

#read in validation metadata
validation_meta <- read_csv("Public data sets/Davenport 2016/E-MTAB-4451 Validation/E-MTAB-4451.sdrf.csv")

#merge validation

validation1 <- merge(validation_meta, validation, by = "Source Name", all.y = T)

identical(Discovery1, validation1)

# check that data frames have identical columns so can merge the two. 
colnames(validation1[1:34])
colnames(Discovery1[1:34])

my_func <- function(x,y) {
  for (i in names(x)) {
    if (!(i %in% names(y))) {
      print('Warning: Names are not the same')
      break
    }  
    else if(i==tail(names(y),n=1)) {
      print('Names are identical')
    }
  }
}
my_func(validation1,Discovery1)

#they do! 


#combine the data frames by adding rows 

df <- rbind(Discovery1, validation1)

colnames(df[1:20])
################################################################################################################################################

# prepare the dataframe for matrix
# select only the genes and sample idea
# transpose
# allign gene names
# select genes from genes derived by LRT. 
colnames(df1[1:100])
df1 <- df %>% dplyr::select(1,7:26191) %>% 
  column_to_rownames("Source Name") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("probeID")

# read in gene name and probe id matching given by Brendon

probenames <- read_csv("Public data sets/COMBAT and MARS merged/Annotations_Bowtie.csv") %>% dplyr:: select(-biotype_gene)

df1 <- merge(probenames, df1, all.y = T, by = "probeID")%>%  mutate(gene_name = coalesce(hgnc_symbol, probeID))

colnames(df1[1:20])
class(probenames$probeID)
df1$gene_name

# read in 651 LRT genes

LRT <- read_csv("states gene set/LRT analysis sepsis timecourse.csv")%>% filter(padj <=0.01)%>% 
  filter(log2FoldChange > 1 |log2FoldChange < -1) %>% 
  mutate(gene_id_2 = coalesce(hgnc_symbol,gene_ID))
colnames(LRT)
#make list
LRT_gene_list <- LRT %>% pull(gene_id_2)

#filter the gene set 
df2 <- df1 %>% 
  filter(gene_name %in% LRT_gene_list) %>% as.data.frame() 

#remove duplicates as these are isoforms
df3 <- df2[!duplicated(df2$hgnc_symbol), ]

rownames(df3) <- NULL

df4 <- df3 %>% 
  column_to_rownames("gene_name") %>% 
  dplyr::select(-hgnc_symbol, -probeID) %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column("Source Name")

df4 <- merge(df[1:6], df4, by = "Source Name", all.y = T)

# now have a dataframe with the genes in the dataset that relate to LRT analaysis and have combined the discovery
# with the validation set ready for the 

####################################################################################################
####################################################################################################
####################################################################################################


# run the analysis 

#library Phate R

colnames(df4)
set.seed(5648)
tree.phate <- phate(df4[7:266], knn= 30, decay=40, knn.dist.method = "euclidean")
phate<- as.data.frame(tree.phate)
df4 <- cbind(df4,phate)

colnames(df4)


#plot
phate_plot <- ggplot(df4, aes(x = PHATE1, y = PHATE2))+
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
ggsave(plot = phate_plot, "Figures/IMMERSE paper/Davenport 2016/phatePlot_test.png", dpi = 300, height = 5, width = 5)

##############################################################################################################
#Run the clustering 

#determine number of clusters
# Elbow method
elbow <- fviz_nbclust(df4[7:266], hcut, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method")
elbow
ggsave(plot = elbow, "Figures/IMMERSE paper/Davenport 2016/elbow.png", height = 5, width = 5, dpi  = 300)

##############################################################################################################################
# run the clustering 



res.hc3 <- eclust(df4[7:266],hc_method = "ward.D2",hc_metric = "euclidean",FUNcluster= "hclust", k	=3, stand = TRUE)


group3 <- as.data.frame(res.hc3$cluster)%>% dplyr::rename("hclust3" ="res.hc3$cluster")

#merge into data frame
df4 <- cbind(df4, group3)

colnames(df4)

hclust3_raw_PHATE <- ggplot(df4 , aes(x = PHATE1, y = PHATE2)) +
  geom_point(aes(fill = as.factor(hclust3)), colour = "black", alpha = 1, size = 4, pch = 21)+
  #scale_fill_manual(values = c("SS1" = "#FF9AA2", "SS2" = "#ffffba", "SS3" = "#A7bed3"))+
  theme_bw()+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(aspect.ratio = 1)

ggsave(plot = hclust3_raw_PHATE, "Figures/IMMERSE paper/aim 3/clustering information/hclust3_raw_PHATE.png", dpi = 300)

##############################################################################################################################
class(df4)
#df4$hclust3 <- as.factor("")
# check states
colnames(df4)
ggplot(df4, aes(x = PHATE1, y = PHATE2))+
  geom_point(aes(fill = CD3D ), pch = 21)

df4 <- df4 %>% mutate(State = case_when(hclust3 == '1' ~ 'State 3',
                                                                    hclust3 == '2' ~ 'State 2',
                                                                    hclust3 == '3' ~ 'State 1',))

##############################################################################################################################
# calculate proportions of each state
colnames(df7)

props <- df4 %>% 
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
  ggtitle("Davenport 2016")

ggsave(plot = pie, "Figures/IMMERSE paper/Davenport 2016/pie.png", dpi = 300, height = 5, width = 5)

######################################################################################################################################
# use this one to create a pseudotime based on 3 clusters
######################################################################################################################################
##Create a Slingshot object and add the vst assay information to it.e.g remove the meta data and dimension reduction so this is not in assay section


colnames(df4)
# add in the normalised data from deseq
sce <- SingleCellExperiment(assays = List(norm = as.matrix(t(df4[7:266]))))


#add the reduced dimensions that I want to use here that is the PHATE cordinates. 
reducedDims(sce) <- SimpleList(PHATE = as.matrix(df4 %>%dplyr::select(c("PHATE1", "PHATE2"))))
#here that is the raw data to test the paths- doesnt work as well and cant plot arrows. 
#reducedDims(sce) <- SimpleList(expressiondata = as.matrix(df_sepsis %>% dplyr::select(13:5012)))

#slingshot requires clusters to run analysis. For us this is the clusters as no timepoint in this data

colData(sce)$State <- as.factor(df4$State)

#run slingshot
sce <- slingshot(sce, clusterLabels = 'State', reducedDim = "PHATE", start.clus = 'State 1')
?slingshot
sce$slingshot


######################################################################################################################################

df4_1 <- cbind(df4, as.data.frame(slingPseudotime(sce, na = FALSE), row.names = FALSE))
df4_1 <- df4_1 %>% rownames_to_column("sample_id")
colnames(df4_1)
#write  the csv for future analysis
#write_csv(vst_1_t_meta_sling, "Results_CSV/vst_1_t_meta_sling.csv")
pseudotimevalues1 <- df4_1%>%dplyr::select(c(PHATE1, PHATE2, Lineage1, sample_id, 
                                                         hclust3, State, Cohort, Sex, Age, `28 day survival`, 
                                             `SRS group`))

head(pseudotimevalues1)
pseudotimevalues1 <- melt(pseudotimevalues1, id.vars=c('State', 'PHATE1', 'PHATE2', 
                                                       'sample_id', 'Cohort', 'Sex', 'Age', '28 day survival', 
                                                       'SRS group',
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

ggsave(plot = PHATE_pseudotime, "Figures/IMMERSE paper/Davenport 2016/PHATE_pseudotime.png", dpi = 300, height = 5, width = 6)

######################################################################################################################################

#PHATE by state 
pseudotimevaluesexclNA1$`SRS group`
PHATE_SRS <- ggplot(pseudotimevaluesexclNA1, aes(x=PHATE1, y=PHATE2))+
  geom_point(aes(fill=`SRS group`),size=4, alpha=1,shape = 21)+
  scale_fill_manual(values = c("SRS1" = "#810103",
                               "SRS2" ="#4d7faf"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)+
  geom_path(data = slingcurves1, aes(x=PHATE1, y=PHATE2), 
            arrow = arrow(type = "open", angle = 40,  length = unit(0.2, "inches")), size =2)

ggsave(plot = PHATE_SRS, "Figures/IMMERSE paper/Davenport 2016/PHATE_SRS.png", dpi = 300, height = 5, width = 6)


#PHATE by SRS
pseudotimevaluesexclNA1
PHATE_state <- ggplot(pseudotimevaluesexclNA1, aes(x=PHATE1, y=PHATE2))+
  geom_point(aes(fill=State ),size=4, alpha=1,shape = 21)+
  scale_fill_manual(values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)+
  geom_path(data = slingcurves1, aes(x=PHATE1, y=PHATE2), 
            arrow = arrow(type = "open", angle = 40,  length = unit(0.2, "inches")), size =2)

ggsave(plot = PHATE_state, "Figures/IMMERSE paper/Davenport 2016/PHATE_states.png", dpi = 300, height = 5, width = 6)


#Check there is no batch effect between the discovery and validation


  PHATE_Cohort <- ggplot(pseudotimevaluesexclNA1, aes(x=PHATE1, y=PHATE2))+
  geom_point(aes(fill=Cohort),size=4, alpha=1,shape = 21)+
  #scale_fill_manual(values = c())+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)+
  geom_path(data = slingcurves1, aes(x=PHATE1, y=PHATE2), 
            arrow = arrow(type = "open", angle = 40,  length = unit(0.2, "inches")), size =2)

ggsave(plot = PHATE_Cohort, "Figures/IMMERSE paper/Davenport 2016/PHATE_Cohort.png", dpi = 300, height = 5, width = 6)


######################################################################################################################################
######################################################################################################################################
#plot the heatmap

colnames(df4)
df4_1_hm <- df4_1 %>%arrange(Lineage1)

colnames(df4_1)
df4_1_hm_clus <- df4_1_hm %>% 
  dplyr::select(`Source Name`,8:267) %>% 
  column_to_rownames("Source Name") %>% 
  scale() %>% 
  t()


d <- dist(df4_1_hm_clus, method = "euclidean") # Euclidean distance matrix.
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
hm_matrix <- df4_1_hm %>% dplyr::select(sample_id,8:267) %>% column_to_rownames("sample_id")%>% scale()
hm_matrix_smooth = gaussianSmooth(hm_matrix, sigma = 10) %>% t()

col_fun = colorRamp2(c(0,0.5, 1), c("blue", "#EEEEEE", "red"))
col_fun(seq(-1, 1))

#make annotation
colnames(df4_1_hm)
anno <- df4_1_hm %>% dplyr::select(Lineage1, State, `SRS group`, `28 day survival`)
colnames(anno)
#create colour scheme for pseudotime´
n = nrow(anno)
min_v = min(anno$Lineage1)
max_v = max(anno$Lineage1)
Var = circlize::colorRamp2(seq(min_v, max_v, length = n), hcl.colors(n,"Spectral"))
#create annotation
ha = HeatmapAnnotation(df = anno, 
                       col = list(Lineage1 = Var,
                                  `SRS group` = c("SRS1" = "#810103",
                                                "SRS2" ="#4d7faf"),
                                  State = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF"),
                                  `28 day survival` = c("survivor" = "#47abd8", "non-survivor" = "#FF4242")),
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

png("Figures/IMMERSE paper/Davenport 2016/heatmap.png",width=14,height=8,units="in",res = 800)
hm
dev.off()

########################################################################







