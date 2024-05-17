
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


# xiao 2011 data set

df <- read_csv("Public data sets/Xiao 2011/GSE36809.csv")

colnames(df)
# read in gene name and probe id matching given by Brendon

probenames <- read_csv("Public data sets/Xiao 2011/GPL570-55999.csv")

df1 <- merge(probenames, df, all.y = T, by = "probeID")%>%  mutate(gene_name = coalesce(`Gene Symbol`, probeID))

colnames(df1)



#need to select the 651 gene set i identified
LRT <- read_csv("states gene set/LRT analysis sepsis timecourse.csv")%>% filter(padj <=0.01)%>% 
  filter(log2FoldChange > 1 |log2FoldChange < -1) %>% 
  mutate(gene_id_2 = coalesce(hgnc_symbol,gene_ID))
colnames(LRT)

#make list
LRT_gene_list <- LRT %>% pull(gene_id_2)

#filter the gene set 
df2 <- df1 %>% 
  filter(gene_name %in% LRT_gene_list)%>% 
  as.data.frame() %>% filter(!row_number() %in% c(1:89))

#remove duplicates as these are isoforms
df3 <- df2[!duplicated(df2$`Gene Symbol`), ]
rownames(df3) <- NULL
df3 <- df3 %>% column_to_rownames("Gene Symbol")%>% 
  dplyr::select(-gene_name, -probeID)%>% 
  t() %>% 
  as.data.frame() %>% rownames_to_column("sample_title")


#read in metadata and make a time point column. 

metadata <- read_csv("Public data sets/Xiao 2011/metadata.csv")

metadata <- metadata %>% mutate(Days = case_when(`Hours since injury` >= 0 & `Hours since injury` <= 12 ~ 0.5,
                                                      `Hours since injury` > 12 & `Hours since injury` <= 24 ~ 1,
                                                      `Hours since injury` > 24 & `Hours since injury` <= 96 ~ 4,
                                                      `Hours since injury` > 96 & `Hours since injury` <= 168 ~ 7,
                                                      `Hours since injury` > 168 & `Hours since injury` <= 336 ~ 14,
                                                      `Hours since injury` > 336 & `Hours since injury` <= 504 ~ 21,
                                                      `Hours since injury` > 504 & `Hours since injury` <= 1000 ~ 28))

                                                

table(metadata)

#merge metadata and gene list

df4 <- merge(df3, metadata)

df4_1 <- df4
df4_1$sample_id <- substr(df4_1$sample_title, 16,24)
length(unique(df4_1$sample_id))

summarise()
count(df4$sample_id)

# there are some NA's in the genes. Need to remove these genes to run the next section. 
colnames(df4)
df4 <- df4 %>%  dplyr::select(1:381)%>%
  select_if(~ !any(is.na(.)))

#now can perform pseudotime magic. 

df4$Days

colnames(df4)
set.seed(5648)
tree.phate <- phate(df4[2:374], knn= 30, decay=5, knn.dist.method = "euclidean")
phate<- as.data.frame(tree.phate)
df5 <- cbind(df4,phate)

colnames(df5)

df4.1 <- merge(df3, metadata)
colnames(df4.1)

df4.1 <- df4.1 %>% dplyr::select(1, 382:388)

df5 <- merge(df5,df4.1)
colnames(df5)

df5$Days <- as.factor(df5$Days)
df5$Days <- str_replace(df5$Days, "Control", "NA")
df5$Days <- df5$Days %>% replace_na('Control')

df5$Days <- factor(df5$Days, levels=c("Control", "0.5", "1", "4", "7", "14", "21", "28"))


#plot
phate_plot <- ggplot(df5, aes(x = PHATE1, y = PHATE2))+
  #geom_density_2d(colour = "grey", alpha = 0.5,binwidth = 0.01)+
  geom_point(aes(fill = Days), pch =21, alpha = 1, size = 2)+
  theme_bw()+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(aspect.ratio = 1)#+
#geom_path(aes(group = case_id), size = 0.3, alpha = 0.5,
#         color="black",
#        arrow = arrow(type = "closed",
#                     length=unit(0.07, "inches"),
#                     ends = "last"),alpha = 0.5)

phate_plot
ggsave(plot = phate_plot, "Figures/IMMERSE paper/Xiao 2011/phatePlot_timepoints.png", dpi = 300, height = 5, width = 5)



#There is a batch effect based on sample data row count metadata column. 
#It does not look biological variation. 
#Will remove. 

batch_effect <- ggplot(df5, aes(x = PHATE1, y = PHATE2))+
  #geom_density_2d(colour = "grey", alpha = 0.5,binwidth = 0.01)+
  geom_point(aes(fill = as.factor(`!Sample_data_row_count`)), pch =21, alpha = 1, size = 2)+
  theme_bw()+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(aspect.ratio = 1)+
  guides(fill=guide_legend(title="Sample_data_row_count
                           batch?"))

ggsave(plot = batch_effect, "Figures/IMMERSE paper/Xiao 2011/batch_effect.png", dpi = 300, height = 5, width = 5)

# remove samples effected by batch effect. 
#re run PHATE
#it removes 45 samples. 
df6 <- df5 %>% filter(!`!Sample_data_row_count` == 54675) %>%  dplyr::select(-PHATE1, -PHATE2)

#run phate
colnames(df6)
set.seed(5648)
tree.phate <- phate(df6[2:374], knn= 30, decay=10, knn.dist.method = "euclidean")
phate<- as.data.frame(tree.phate)
df6 <- cbind(df6,phate)

df6$Group
#plot phate by days
phate_plot1 <- ggplot(df6, aes(x = PHATE1, y = PHATE2))+
  #geom_density_2d(colour = "grey", alpha = 0.5,binwidth = 0.01)+
  geom_point(aes(fill = Days, shape = Group), alpha = 1, size = 2.5)+
  scale_fill_manual(values = c("0.5" = "#E31A1C", "1" = "#FC6A03", "4" =   "#FB9A99", "7" =  "#33A02C",
                               "14" =  "#B2DF8A" , "21" = "#1F78B4", "28" = "#A6CEE3", "Control" = "#A6CEE3" ))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  scale_shape_manual(values = c(24,21))+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(aspect.ratio = 1)#+

ggsave(plot = phate_plot1, "Figures/IMMERSE paper/Xiao 2011/phate_plot1.png", dpi = 300, height = 5, width = 6)

# can see there is an outlier group. legit or batch? Have not seen this in  other data set. 



##############################################################################################################################
# run the clustering 

df7 <- df6 %>% filter(!Days ==  "Control")
#count number of participants
colnames(df7)[1:10]

subjects <- unique(df7$sample_id) #  755 unique sample_ids so this means 755 patients
colnames(df7)
res.hc3 <- eclust(df7[2:374],hc_method = "ward.D2",hc_metric = "euclidean",FUNcluster= "hclust", k	=3, stand = TRUE)


group3 <- as.data.frame(res.hc3$cluster)%>% dplyr::rename("hclust3" ="res.hc3$cluster")

#merge into data frame
df7 <- cbind(df7, group3)
ggplot

#test
ggplot(df7, aes(x = PHATE1, y = PHATE2))+
  geom_point(aes(fill = CD3D ), pch = 21)

ggplot(df7, aes(x = PHATE1, y = PHATE2))+
  geom_point(aes(fill = as.factor(hclust3 )), pch = 21)

df7 <- df7 %>% mutate(State = case_when(hclust3 == '1' ~ 'State 1',
                             hclust3 == '2' ~ 'State 2',
                             hclust3 == '3' ~ 'State 3',))
colnames(df6)
df6$hclus
#quick check to see where clusters lie
hclust3_raw_PHATE <- ggplot(df7 , aes(x = PHATE1, y = PHATE2)) +
  geom_point(aes(fill = as.factor(State)), colour = "black", alpha = 1, size = 4, pch = 21)+
  scale_fill_manual(values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF"))+
  theme_bw()+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(aspect.ratio = 1)


ggsave(plot = hclust3_raw_PHATE, "Figures/IMMERSE paper/Xiao 2011/phate_plot_states.png", dpi = 300, height = 5, width = 5)


library(viridis)
colnames(df7)
CD3 <- ggplot(df7, aes(x=PHATE1, y=PHATE2))+
  geom_point(aes(fill=`CD3D`),size=4, alpha=1,shape = 21)+
  scale_fill_viridis(discrete = F)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)

ggsave(plot = CD3, "Figures/IMMERSE paper/Xiao 2011/CD3.png", dpi = 300, height = 5, width = 5)


######################################################################################################################################
##############################################################################################################################
# calculate proportions of each state
colnames(df2)

props <- df7%>% filter(Days == 1) %>% # NOTE this is at day one for the pie chart to compare between the other endotypes
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
  ggtitle("Xiou 2011")

ggsave(plot = pie, "Figures/IMMERSE paper/Xiao 2011/pie_Day1 only.png", dpi = 300, height = 5, width = 5)


#use all samples
props <- df7%>% filter(Days == 1) %>% # NOTE this is at day one for the pie chart to compare between the other endotypes
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
  ggtitle("Xiou 2011")

ggsave(plot = pie, "Figures/IMMERSE paper/Xiao 2011/pie_all samples.png", dpi = 300, height = 5, width = 5)


######################################################################################################################################
# use this one to create a pseudotime based on 3 clusters
######################################################################################################################################
##Create a Slingshot object and add the vst assay information to it.e.g remove the meta data and dimension reduction so this is not in assay section


colnames(df7)
# add in the normalised data from deseq
sce <- SingleCellExperiment(assays = List(norm = as.matrix(t(df7[2:374]))))


#add the reduced dimensions that I want to use here that is the PHATE cordinates. 
reducedDims(sce) <- SimpleList(PHATE = as.matrix(df7 %>%dplyr::select(c("PHATE1", "PHATE2"))))
#here that is the raw data to test the paths- doesnt work as well and cant plot arrows. 
#reducedDims(sce) <- SimpleList(expressiondata = as.matrix(df_sepsis %>% dplyr::select(13:5012)))

#slingshot requires clusters to run analysis. For us this is the clusters as no timepoint in this data

colData(sce)$State <- as.factor(df7$State)

#run slingshot
sce <- slingshot(sce, clusterLabels = 'State', reducedDim = "PHATE", start.clus = 'State 1')
?slingshot
sce$slingshot

######################################################################################################################################

df7_1 <- cbind(df7, as.data.frame(slingPseudotime(sce, na = FALSE), row.names = FALSE))
df7_1 <- df7_1 %>% rownames_to_column("sample_id")
colnames(df7_1)
#write  the csv for future analysis
#write_csv(vst_1_t_meta_sling, "Results_CSV/vst_1_t_meta_sling.csv")
pseudotimevalues1 <- df7_1%>%dplyr::select(c(PHATE1, PHATE2, Lineage1, sample_title, 
                                             hclust3, State, Sex, Age, `Hours since injury`, Days))

head(pseudotimevalues1)
pseudotimevalues1 <- melt(pseudotimevalues1, id.vars=c('State', 'PHATE1', 'PHATE2', 
                                                       'sample_title',
                                                       'hclust3',
                                                       "Sex","Age", "Hours since injury","Days"
                                                       ), variable.name='lineage', value.name = 'pseudotime')



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
  geom_point(aes(fill=pseudotime),size=2.5, alpha=1,shape = 21)+
  scale_fill_gradientn(colours=colors1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)+
  geom_path(data = slingcurves1, aes(x=PHATE1, y=PHATE2), 
            arrow = arrow(type = "open", angle = 40,  length = unit(0.2, "inches")), size =2)


ggsave(plot = PHATE_pseudotime, "Figures/IMMERSE paper/Xiao 2011/PHATE_pseudotime.png", dpi = 300, height = 5, width = 5)


#PHATE by state 
pseudotimevaluesexclNA1$`SRS group`
PHATE_State <- ggplot(pseudotimevaluesexclNA1, aes(x=PHATE1, y=PHATE2))+
  geom_point(aes(fill=State),size=2.5, alpha=1,shape = 21)+
  scale_fill_manual(values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)+
  geom_path(data = slingcurves1, aes(x=PHATE1, y=PHATE2), 
            arrow = arrow(type = "open", angle = 40,  length = unit(0.2, "inches")), size =2)

ggsave(plot = PHATE_State, "Figures/IMMERSE paper/Xiao 2011/PHATE_State.png", dpi = 300, height = 5, width = 6)



######################################################################################################################################
######################################################################################################################################
#plot the heatmap

colnames(df7_1_hm)
df7_1_hm <- df7_1 %>%arrange(Lineage1)


df7_1_hm_clus <- df7_1_hm %>% 
  dplyr::select(`sample_title`,2:374) %>% 
  column_to_rownames("sample_title") %>% 
  scale() %>% 
  t()


d <- dist(df7_1_hm_clus, method = "euclidean") # Euclidean distance matrix.
H.fit <- hclust(d, method="ward.D")
plot(H.fit) # display dendogram
groups <- cutree(H.fit, k=3) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
plot(H.fit)
rect.hclust(H.fit, k=3, border="red") 
group <- as.data.frame(groups) %>% rownames_to_column("gene_name")#extract grouping variable 


# merge gene_biotype_df with gene group
colnames(group)
colnames(gene_biotype_df)

group1 <- group %>%  mutate(`Gene cluster` = case_when(groups == "1" ~ "Cluster 1",
                                                       groups == "2" ~ "Cluster 2",
                                                       groups == "3" ~ "Cluster 3"))
group2 <- group1 %>% dplyr::select(-groups) %>% column_to_rownames("gene_name")

table(group)
head(group2)
## merge cluster group with dataframe

######################################################################################################################################
library("mmand")
colnames(df7_1_hm)
hm_matrix <- df7_1_hm %>% dplyr::select(sample_title,2:374) %>% column_to_rownames("sample_title")%>% scale()
hm_matrix_smooth = gaussianSmooth(hm_matrix, sigma = 10) %>% t()

col_fun = colorRamp2(c(0,0.5, 1), c("blue", "#EEEEEE", "red"))
col_fun(seq(-1, 1))

#make annotation
colnames(df7_1_hm)
anno <- df7_1_hm %>% dplyr::select(Lineage1, State, Days)
colnames(anno)
#create colour scheme for pseudotime´
n = nrow(anno)
min_v = min(anno$Lineage1)
max_v = max(anno$Lineage1)
Var = circlize::colorRamp2(seq(min_v, max_v, length = n), hcl.colors(n,"Spectral"))
#create annotation
ha = HeatmapAnnotation(df = anno, 
                       col = list(Lineage1 = Var,
                                  State = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF"),
                                  Days = c("0.5" = "#E31A1C", "1" = "#FC6A03", "4" =   "#FB9A99", "7" =  "#33A02C",
                                           "14" =  "#B2DF8A" , "21" = "#1F78B4", "28" = "#A6CEE3")),
                       which = "col",
                       gp = gpar(col = "black"))


# make the column annotation

ha_row = rowAnnotation(df = group2, 
                       col = list(`Gene cluster` = c("Cluster 1" =  "#FF756D",
                                                     "Cluster 2" = "#F6B4C3", 
                                                     "Cluster 3" = "#AAD6FA"),
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

png("Figures/IMMERSE paper/Xiao 2011/heatmap.png",width=20,height=8,units="in",res = 800)
hm
dev.off()

# I can correlate the hours since admission and pseudotime. 
# I dont think these two variables will correlate 

library(ggpubr)
colnames(pseudotimevaluesexclNA1)
class(pseudotimevaluesexclNA1$lineage)
cor <- ggscatter(pseudotimevaluesexclNA1, x = "pseudotime", y = "Hours since injury",
          color = "black", shape = 21, size = 3, fill = "State", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "spearman", label.x = 0.01, label.sep = "\n")
)+
  scale_fill_manual(values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF"))

ggsave(plot = cor, "Figures/IMMERSE paper/Xiao 2011/pseudo_hours_cor.png", dpi = 300, height = 5, width = 5)

violin_plot <- ggplot(pseudotimevaluesexclNA1, aes(x = pseudotime, y = Days))+
  geom_violin()+
  geom_quasirandom(aes(fill = State), pch = 21, size = 3)+
  scale_fill_manual(values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)

ggsave(plot = violin_plot, "Figures/IMMERSE paper/Xiao 2011/pseudo_violin.png", dpi = 300, height = 5, width = 6)

colnames(pseudotimevaluesexclNA1)
#make a new column so I have the patient ID 16:25

pseudotimevaluesexclNA1$sample_id <- substr(pseudotimevaluesexclNA1$sample_title, 16,24)

ggplot(pseudotimevaluesexclNA1, aes(x = pseudotime, y = Days, group = sample_id)) + 
  geom_point() + geom_line()
pseudotimevaluesexclNA1$sample_title

head(pseudotimevaluesexclNA1)

# Summarise counts 
pseudotimevaluesexclNA1_1 <- pseudotimevaluesexclNA1 %>% 
  dplyr::select(sample_title, sample_id,pseudotime, Days) %>% 
  # for each gene
  group_by(sample_id,Days) %>% 
  # scale the cts column
  mutate(mean_pseudotime = mean(pseudotime))

colnames(pseudotimevaluesexclNA1_1)
pseudotimevaluesexclNA1_1 <- pseudotimevaluesexclNA1 %>% 
  dplyr::select(sample_id ,pseudotime, Days) %>% 
  
  spread(Days,pseudotime) 
  gather
%>% 
  # for each gene, strain and minute
  group_by(gene, strain, minute) %>%
  # calculate the mean (scaled) cts
  summarise(mean_cts_scaled = mean(cts_scaled),
            nrep = n()) %>% 
  ungroup()


ggplot(pseudotimevaluesexclNA1,aes(Days, pseudotime)) +
  geom_line(aes(group = sample_id), alpha = 0.3)


ggplot(pseudotimevaluesexclNA1,aes(`Hours since injury`, pseudotime)) +
  geom_line(aes(group = sample_id), alpha = 0.3)+
  scale_x_log10()

#cluster patient
head(df7)
head(pseudotimevaluesexclNA1)
pseudotimevaluesexclNA1_1 <- pseudotimevaluesexclNA1 %>% dplyr:: select(sample_id, pseudotime)


pseudotimevaluesexclNA1_1 <- setNames(data.frame(t(pseudotimevaluesexclNA1_1[,-1])), pseudotimevaluesexclNA1_1[,1])
pseudotimevaluesexclNA1_1 <- pseudotimevaluesexclNA1_1 %>% as.matrix()
hclust(pseudotimevaluesexclNA1, method = "complete")

clust <- eclust(pseudotimevaluesexclNA1[11],hc_method = "ward.D2",hc_metric = "euclidean",FUNcluster= "hclust", k	=2, stand = TRUE)

clust_dend <- fviz_dend(clust, rect = TRUE, k_colors = c("1" = "#FDE725FF","2" ="#29AF7FFF"),
                          color_labels_by_k = FALSE,
                          cex = 0.2,
                          show_labels = T) # dendrogram

group <- as.data.frame(clust$cluster)%>% dplyr::rename("patient_clust" ="clust$cluster")

#merge into data frame
pseudotimevaluesexclNA1_1 <- cbind(pseudotimevaluesexclNA1, group)

colnames(pseudotimevaluesexclNA1_1)

ggplot(pseudotimevaluesexclNA1_1,aes(Days, pseudotime)) +
  geom_line(aes(group = sample_id,colour = as.factor(patient_clust)), alpha = 0.3)+
  facet_wrap(~as.factor(patient_clust))

head(pseudotimevaluesexclNA1)

cluster <- pseudotimevaluesexclNA1 %>% dplyr::select(sample_id,Days, pseudotime) %>%
  dplyr::slice(-c(70, 111, 115, 251, 300, 364, 373, 384, 391, 396, 421, 463, 491, 502, 533, 566, 577, 587, 615, 620, 636,
         663, 668, 677, 756, 20, 40, 46, 166, 286, 336, 398, 432, 485, 550, 638, 741, 125, 329, 554, 585, 667, 760,
         147, 426, 37, 58, 87, 100, 133, 182, 272, 299, 340, 346, 380, 415, 451, 498, 557, 614, 626, 635, 684, 709,
         716,731) )%>% spread(Days, pseudotime)

colnames(cluster)
cluster1 <- cluster %>% 
  column_to_rownames("sample_id") %>% as.data.frame()

class(cluster1)
library(ClustImpute)
library(tictoc)
nr_iter <- 10 # iterations of procedure
n_end <- 10 # step until convergence of weight function to 1
nr_cluster <- 2 # number of clusters
c_steps <- 50 # numer of cluster steps per iteration
tictoc::tic("Run ClustImpute")
res <- ClustImpute(cluster1,nr_cluster=nr_cluster, nr_iter=nr_iter, c_steps=c_steps, n_end=n_end) 
tictoc::toc()


clusters <- res$clusters

cluster2 <- cbind(cluster1, clusters)

cluster2 <- cluster2 %>% rownames_to_column("sample_id") %>%  gather(Days, pseudotime, 2:8)

cluster2$Days <- factor(cluster2$Days, levels=c('0.5', '1', '4', '7', '14', '21', '28'))



psudotime_cluster <- ggplot(cluster2,aes(Days, pseudotime)) +
  geom_line(aes(group = sample_id,colour = as.factor(clusters)), alpha = 0.1)+
    geom_smooth( size = 1.5, 
              aes(group = clusters, colour = as.factor(clusters)))+
  theme_classic()
  
ggsave(plot = psudotime_cluster, "Figures/IMMERSE paper/Xiao 2011/pseudotime_cluster.png", dpi = 300, height = 5, width = 6)

  ggplot(pseudotimevaluesexclNA1,aes(Days, pseudotime)) +
    geom_line(aes(group = sample_id), alpha = 0.3)


colnames(cluster)
head(cluster)
hclust_matrix <- pseudotimevaluesexclNA1 %>% 
  select(-gene) %>% 
  as.matrix()

# assign rownames
rownames(hclust_matrix) <- trans_cts$gene




########################################################################
# save the new metadata
write_csv(pseudotimevaluesexclNA1, "pseudo_metadata.csv")

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
#run the stats on whole geneset 
#need to run stats based on pseudotime cluster

# clear the environment
rm(list = ls())

#load in data

df <- read_csv("Public data sets/Xiao 2011/GSE36809.csv") %>% na.omit()

#metadata
metadata <- read_csv("CSV output /Xioa 2011/pseudo_metadata.csv")






df_log <- df_1 %>% mutate_at(vars(49:76), log2) %>% select(49:76, sample_id, Cohort_time) %>% as.data.frame() %>% 
  na.omit()
analysis <- df_log %>% select("sample_id", 1:28) 
colnames(analysis)

#make matrix of values samples columns and frequencies as rows
count_comparison_tmydf = setNames(data.frame(t(analysis[,-1])), analysis[,1])%>% as.matrix() 
str(count_comparison_tmydf)

dim(count_comparison_tmydf)
dim(design)

f <- factor(df_log$Cohort_time, levels=c("CC_T0","CC_T1", "SC_T0"))
design <- model.matrix(~0+f)
colnames(design) <- c("CC_T0","CC_T1","SC_T0")
#To make all pair-wise comparisons between the three groups one could proceed
fit <- lmFit(count_comparison_tmydf, design)
?lmFit
contrast.matrix <- makeContrasts(CC_T1-CC_T0, SC_T0-CC_T0, SC_T0-CC_T1,
                                 levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2)
?eBayes
df_stats1 <- topTable(fit2, coef=1, adjust="BH", number = "inf") %>% rownames_to_column("Cytokine")
write_csv(df_stats1, "Stats csv/CC_T0vsCC_T1.csv")




df1 <- merge(probenames, df, all.y = T, by = "probeID")%>%  mutate(gene_name = coalesce(`Gene Symbol`, probeID))

BiocManager::install("GEOquery")
library(GEOquery)
# load series and platform data from GEO

gset <- getGEO("GSE36809", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }

