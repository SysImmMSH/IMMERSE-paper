#Build the clustering modelling based on most significant genes from the LRT test. 
# run the pseudotime analysis 

#load packages
library(tidyverse)
library(ggplot2)
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
library(FlowSOM)
library(cytofkit2)
library(stats)
library(stringr)
library(dbscan)


# make list of the significant genes required for the analysis 
#import the LRT timecourse analysis 
getwd()
LRT <- read_csv("Results_csv/IMMERSE paper/SCRIPT 1/LRT analysis sepsis timecourse.csv")%>% filter(padj <=0.01)%>% 
  filter(log2FoldChange > 1 |log2FoldChange < -1) %>% 
  mutate(gene_id_2 = coalesce(hgnc_symbol,gene_ID))
colnames(LRT)
#make list
LRT_gene_list <- LRT %>% pull(gene_ID)

# import the transformed normalised VSD df
df <- read_csv("Results_CSV/IMMERSE paper/SCRIPT 1/vsd_df.csv") %>% dplyr::select(gene_ID, 2:259) 

#filter genes of interest
df <- df %>% 
  filter(gene_ID %in% LRT_gene_list) %>% as.data.frame()



library("biomaRt")
#the useEnsembl line worked as had to redirect to useast server. choose between useEnsembl and useMart option
#ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = 'useast')
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl", host = "useast.ensembl.org")
#mart = useMart("ensembl")
genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                  filters = "ensembl_gene_id",
                  values = df$gene_ID,
                  mart = ensembl)
idx <- match(df$gene_ID, genemap$ensembl_gene_id)
df$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
df$gene_biotype <- genemap$gene_biotype[ idx ]

colnames(df)

df1 <- df%>%  mutate(hgnc_symbol = na_if(hgnc_symbol, ""))  %>%  mutate(gene_name = coalesce(hgnc_symbol, gene_ID)) %>% dplyr::select(262, 2:259)
colnames(df1)
#need to make genes variables and patients observations
df2 = setNames(data.frame(t(df1[,-1])), df1[,1]) #transpose so obvs = patients, variables = genes

colnames(df2)
#now have a data frame with 5000 of the most variable genes in the dataset.
df2 <- df2 %>% rownames_to_column("sample_id")

colnames(df2)
#read in metadata
colnames(metadata)
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

sweeney <- read_csv("Results_csv/IMMERSE paper/SCRIPT 1/sweeney_groups.csv")
metadata <- merge(metadata, sweeney, by = "sample_id", all.x = TRUE)
colnames(metadata)

#read in clinical
merged_baseline <- read_csv("clinical_csv/merged baseline.csv") %>% dplyr::select(case_id,`Infection source2`)
colnames(merged_baseline)
metadata <- merge(metadata, merged_baseline, by = "case_id", all.x = TRUE)

#read in outcome
merged_baseline <- read_csv("clinical_csv/df_outcome_test for sce.csv") %>% dplyr::select(case_id,`outcome_hospital_survival`)
colnames(merged_baseline)
metadata <- merge(metadata, merged_baseline, by = "case_id", all.x = TRUE)

df3 <- merge(metadata, df2, all.x = TRUE)

df3$outcome_hospital_survival
#Dataframe is now built

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

# Run the unsupervised analysis  
# PHATE dimension reduction analysis


#library Phate R
library(phateR)
colnames(df3)
set.seed(5648)
df4 <- df3 #%>% filter(!Cohort_time == "Cardiac_T1")
tree.phate <- phate(df4[16:733], knn= 50, decay=40, knn.dist.method = "euclidean") #ALWAYS CHECK THE SUBSETTING
phate<- as.data.frame(tree.phate)
df4.1 <- cbind(df4,phate)

#plot
phate <- ggplot(df4.1, aes(x = PHATE1, y = PHATE2))+
  #geom_density_2d(colour = "grey", alpha = 0.5,binwidth = 0.01)+
  geom_point(aes(fill = Cohort_time, shape = Cohort), colour = "black", alpha = 0.9, size = 4)+
  scale_fill_manual(values = c("Sepsis_T0" = "#efe350ff",
                               "Sepsis_T1" = "#f68f46ff",
                               "Sepsis_T2" = "#a65c85ff",
                               "Sepsis_T3" = "#403891ff",
                               "Cardiac_T0" ="#56A8CBFF",
                               "Cardiac_T1" = "#DA291CFF"))+
  scale_shape_manual(values = c("Sepsis" = 21, "Cardiac" = 24))+
  theme_bw()+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(aspect.ratio = 1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(legend.position = "none")+
  theme(text = element_text(size = 20))#+
#geom_path(aes(group = case_id), size = 0.3, alpha = 0.5,
#         color="black",
#        arrow = arrow(type = "closed",
#                     length=unit(0.07, "inches"),
#                     ends = "last"),alpha = 0.5)

phate
ggsave(plot = phate, "FIGURES/IMMERSE paper/SCRIPT 3/phatePlot_nolegend.png", dpi = 300, height = 5, width = 5)



#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#I run this jsut to see other dimension reduction technique
library("destiny")

set.seed(458)
diffusion.map <- DiffusionMap(df4[16:733], verbose = TRUE)
dm<- as.data.frame(diffusion.map@eigenvectors)
df4.1 <- cbind(df4.1,dm[1:4])

colnames(dm)
#view briefly

diiffusionMap <- ggplot(df4.1, aes(x = DC1, y = DC2))+
  #geom_density_2d(colour = "grey", alpha = 0.5,binwidth = 0.01)+
  geom_point(aes(fill = Cohort_time, shape = Cohort), colour = "black", alpha = 1, size = 4)+
  scale_fill_manual(values = c("Sepsis_T0" = "#efe350ff",
                               "Sepsis_T1" = "#f68f46ff",
                               "Sepsis_T2" = "#a65c85ff",
                               "Sepsis_T3" = "#403891ff",
                               "Cardiac_T0" ="#56A8CBFF",
                               "Cardiac_T1" = "#DA291CFF"))+
  scale_shape_manual(values = c("Sepsis" = 21, "Cardiac" = 24))+
  theme_bw()+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(aspect.ratio = 1)


ggsave(plot = diiffusionMap, "Figures/IMMERSE paper/SCRIPT 3/diiffusionMap.png", dpi = 300, height = 10, width = 10)



###################################################################
###### choosing how many clusters of genes#########################
#Save the data of the cardiac patients and sepsis patients (only the phate cordinates, timepoint Cohort_ttime, Cohort, sample_id)
colnames(df4.1)
save <- df4.1 %>% select(case_id, sample_id, Cohort, Cohort_time, PHATE1, PHATE2, Timepoint)
write_csv(df4.1, "Results_csv/IMMERSE paper/SCRIPT 3/metadata of PHATE.csv")
# first create data frame to only include sepsis patients. 

df_sepsis <- df4.1 %>% filter(Cohort == "Sepsis")
table(df_sepsis$Timepoint)
#play with what data to use. 
# evidence to cluster on UMAP https://link.springer.com/chapter/10.1007/978-3-030-51935-3_34#Sec2
colnames(df_sepsis)
# Elbow method
elbow <- fviz_nbclust(df_sepsis[15:732], hcut, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method")
elbow
ggsave(plot = elbow, "FIGURES/IMMERSE paper/SCRIPT 3/elbow.png", height = 5, width = 5, dpi  = 300)


######################################################################################################################################
######################################################################################################################################
#run the clustering
df_sepsis
# RAW data clustering
#run clustering for 3 groups 

df_sepsis <- df_sepsis %>% column_to_rownames("sample_id")
colnames(df_sepsis)
res.hc3 <- eclust(df_sepsis[16:733],hc_method = "ward.D2",hc_metric = "euclidean",FUNcluster= "hclust", k	=3, stand = TRUE)
df_sepsis <- df_sepsis %>% rownames_to_column("sample_id")


group3 <- as.data.frame(res.hc3$cluster)%>% dplyr::rename("hclust3" ="res.hc3$cluster")
clust_dend_3 <- fviz_dend(res.hc3, rect = TRUE, k_colors = c("2" = "#FDE725FF","1" ="#29AF7FFF","3" = "#453781FF"),
                          color_labels_by_k = FALSE,
                          cex = 0.2) # dendrogam
ggsave(plot = clust_dend_3, "FIGURES/IMMERSE paper/SCRIPT 3/dend_hclust3.png", dpi = 300, height = 5, width = 10)

##############################################################################################################
##############################################################################################################
#merge into data frame
df_sepsis <- cbind(df_sepsis, group3)
df_sepsis$hclust3

######################################################################################################################################
######################################################################################################################################
# view clustering on PHATE RAW data
df_sepsis$DC1
# hclust3 raw data 
hclust3_raw_PHATE <- ggplot(df_sepsis , aes(x = PHATE1, y = PHATE2)) +
  geom_point(aes(fill = as.factor(hclust3)), colour = "black", alpha = 1, size = 4, pch = 21)+
  #scale_fill_manual(values = c("SS1" = "#FF9AA2", "SS2" = "#ffffba", "SS3" = "#A7bed3"))+
  theme_bw()+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(aspect.ratio = 1)

ggsave(plot = hclust3_raw_PHATE, "FIGURES/IMMERSE paper/SCRIPT 3/hclust3_raw_PHATE.png", dpi = 300)


# hclust3 raw data 
hclust3_raw_DC <- ggplot(df_sepsis , aes(x = DC1, y = DC2)) +
  geom_point(aes(fill = as.factor(hclust3)), colour = "black", alpha = 1, size = 4, pch = 21)+
  #scale_fill_manual(values = c("SS1" = "#FF9AA2", "SS2" = "#ffffba", "SS3" = "#A7bed3"))+
  theme_bw()+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(aspect.ratio = 1)

ggsave(plot = hclust3_raw_DC, "FIGURES/IMMERSE paper/SCRIPT 3/hclust3_raw_DC.png", dpi = 300)




######################################################################################################################################
# use this one to create a pseudotime based on 3 clusters
######################################################################################################################################
##Create a Slingshot object and add the raw cuont assay information to it.e.g remove the meta data and dimension reduction so this is not in assay section

#get list of genes from vst transformation most variable
# import the raw counts VSD df
samples <- df_sepsis$sample_id
raw_count <- read_csv("Results_CSV/raw count data.csv") %>% dplyr::select(gene_ID, 2:259) %>% 
  column_to_rownames("gene_ID") %>% 
  dplyr::select(samples) %>% 
  as.matrix()
colnames(raw_count)

gene_names <- as.vector(df1$gene_id)
colnames(df_sepsis)
# add in the normalised data from deseq
sce <- SingleCellExperiment(assays = List(counts = raw_count),
                            colData=DataFrame(df_sepsis[1:14]))
sce

#add the reduced dimensions that I want to use here that is the PHATE cordinates. 
reducedDims(sce) <- SimpleList(PHATE = as.matrix(df_sepsis%>%dplyr::select(c("PHATE1", "PHATE2"))))
#here that is the raw data to test the paths- doesnt work as well and cant plot arrows. 
#reducedDims(sce) <- SimpleList(expressiondata = as.matrix(df_sepsis %>% dplyr::select(13:5012)))

#slingshot requires "clusters" to run analysis. For us this is the timepoints

#colData(sce)$Timepoint <- as.factor(df_sepsis$Timepoint)

#run slingshot
sce <- slingshot(sce, clusterLabels = 'Timepoint', reducedDim = "PHATE", start.clus = 'T0')




######################################################################################################################################
######################################################################################################################################
# can be a quick plot in base r
library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PHATE, col = plotcol, pch=16,)
lines(SlingshotDataSet(sce), lwd=2, col='black')

# PLOT in ggplot

df_sepsis_sling <- cbind(df_sepsis, as.data.frame(slingPseudotime(sce, na = FALSE), row.names = FALSE))
colnames(df_sepsis_sling)
#write  the csv for future analysis
#write_csv(vst_1_t_meta_sling, "Results_CSV/vst_1_t_meta_sling.csv")
pseudotimevalues1 <- df_sepsis_sling%>%dplyr::select(c(PHATE1, PHATE2, Timepoint, Lineage1, Cohort, Cohort_time, sample_id, 
                                                   hclust3,SRS_davenport, SRSq_davenport,SRS_extended,SRSq_extended,MARS, `Sweeney 2019`, `Infection source2`,
                                                   outcome_hospital_survival))

head(pseudotimevalues1)
pseudotimevalues1 <- melt(pseudotimevalues1, id.vars=c('Timepoint', 'PHATE1', 'PHATE2', 
                                                       'sample_id', 'Cohort', 'Cohort_time', 'SRS_davenport', 'SRSq_davenport','SRS_extended','SRSq_extended','MARS',
                                                       'Sweeney 2019','Infection source2','outcome_hospital_survival',
                                                        'hclust3'), variable.name='lineage', value.name = 'pseudotime')



#exclude any NA's
pseudotimevaluesexclNA1 <- pseudotimevalues1%>%filter(pseudotime!='NA')

pseudotimevaluesexclNA1 <- pseudotimevaluesexclNA1 %>% mutate(State = case_when(hclust3 == "1" ~ "State 1",
                                                    hclust3 == "2" ~ "State 2",
                                                    hclust3 == "3" ~ "State 3"))

#éxtract the curve object for plotting the curve in GGplot 
slingcurves1 <- as.data.frame(slingCurves(sce)[[1]]$s[slingCurves(sce)[[1]]$ord, ])
#slingcurves2 <- as.data.frame(slingCurves(sce)[[2]]$s[slingCurves(sce)[[2]]$ord, ]) # use this if other lineages
#slingcurves3 <- as.data.frame(slingCurves(sce)[[3]]$s[slingCurves(sce)[[3]]$ord, ]) # use this if other lineages

#make sure the pseudotime is numeric

#generate colorpalette
colors <- colorRampPalette(rev(brewer.pal(11, 'Spectral'))[-6])(100)
colors1 <- rev(colors)
PHATE_pseudotime <- ggplot(pseudotimevaluesexclNA1%>%arrange(pseudotime), aes(x=PHATE1, y=PHATE2))+
  geom_point(aes(fill=pseudotime),size=5, alpha=1,shape = 21)+
  scale_fill_gradientn(colours=colors1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)+
  geom_path(data = slingcurves1, aes(x=PHATE1, y=PHATE2), 
            arrow = arrow(type = "open", angle = 40,  length = unit(0.2, "inches")), size =2)+
  theme(text = element_text(size = 20))+
  theme(legend.position = "bottom",legend.direction="horizontal", legend.text = element_text(size =10))
  

ggsave(plot = PHATE_pseudotime, "FIGURES/IMMERSE paper/SCRIPT 3/PHATE_pseudotime.png", dpi = 300, height = 5, width = 5)

#PHATE by state 
pseudotimevaluesexclNA1$hclust3
PHATE_state <- ggplot(pseudotimevaluesexclNA1, aes(x=PHATE1, y=PHATE2))+
  geom_point(aes(fill=State),size=5, alpha=1,shape = 21)+
  scale_fill_manual(values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)+
  geom_path(data = slingcurves1, aes(x=PHATE1, y=PHATE2), 
            arrow = arrow(type = "open", angle = 40,  length = unit(0.2, "inches")), size =2)+
  theme(legend.position = "none")+
  theme(text = element_text(size = 20))

ggsave(plot = PHATE_state, "FIGURES/IMMERSE paper/SCRIPT 3/PHATE_states_nolegend.png", dpi = 300, height = 5, width = 5)


######################################################################################################################################
# calculate proportions of each state
colnames(pseudotimevaluesexclNA1)
pseudotimevaluesexclNA1
library(dplyr)
props <- pseudotimevaluesexclNA1 %>% 
  group_by(State) %>% # Variable to be transformed
  dplyr::count() %>% 
  ungroup() %>% 
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(perc= perc*100)
  
props$perc <-  signif(props$perc, digits = 3)

props$perc <- as.character(props$perc)
 
props$perc <- paste(props$perc, "%", sep = "")

pie <- ggplot(props, aes(x = "", y = perc, fill = State)) +
  geom_col(color = "black") +
  geom_label(aes(label = perc), color = c("black"),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE, size = 8) +
  theme(legend.position = "none")+
  guides(fill="none")+
  scale_fill_manual(values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF"))+
  coord_polar(theta = "y") + 
  theme_void()+
  ggtitle("Sepsis immune states")

ggsave(plot = pie, "FIGURES/IMMERSE paper/SCRIPT 3/pie.png", dpi = 300, height = 5, width = 5)



#PHATE by davenport 
pseudotimevaluesexclNA1$hclust3
PHATE_SRS_DP <- ggplot(pseudotimevaluesexclNA1, aes(x=PHATE1, y=PHATE2))+
  geom_point(aes(fill=SRS_davenport),size=7, alpha=1,shape = 21)+
  scale_fill_manual(values = c("SRS1" = "#810103",
                               "SRS2" ="#4d7faf",
                               "SRS3" = "#4d7faf"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)+
  geom_path(data = slingcurves1, aes(x=PHATE1, y=PHATE2), 
            arrow = arrow(type = "open", angle = 40,  length = unit(0.2, "inches")), size =2)+
  theme(legend.position = "none")+
  theme(text=element_text(size = 22))

ggsave(plot = PHATE_SRS_DP, "FIGURES/IMMERSE paper/SCRIPT 3/PHATE_SRS_DP.png", dpi = 300, height = 7, width = 7)

#PHATE by extended 
pseudotimevaluesexclNA1$hclust3
PHATE_SRS_EX <- ggplot(pseudotimevaluesexclNA1, aes(x=PHATE1, y=PHATE2))+
  geom_point(aes(fill=SRS_extended),size=7, alpha=1,shape = 21)+
  scale_fill_manual(values = c("SRS1" = "#810103",
                               "SRS2" ="#4d7faf",
                               "SRS3" = "#010089"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)+
  geom_path(data = slingcurves1, aes(x=PHATE1, y=PHATE2), 
            arrow = arrow(type = "open", angle = 40,  length = unit(0.2, "inches")), size =2)+
  theme(legend.position = "none")+
  theme(text=element_text(size = 22))

ggsave(plot = PHATE_SRS_EX, "FIGURES/IMMERSE paper/SCRIPT 3/PHATE_SRS_EX.png", dpi = 300, height = 7, width = 7)

#PHATE by MARS 
pseudotimevaluesexclNA1$MARS
PHATE_SRS_MARS <- ggplot(pseudotimevaluesexclNA1, aes(x=PHATE1, y=PHATE2))+
  geom_point(aes(fill=MARS),size=7, alpha=1,shape = 21)+
  scale_fill_manual(values = c("Mars1" = "#c7deb2",
                               "Mars2" ="#2e64aa",
                               "Mars3" = "#9dacd5",
                               "Mars4" = "#82ba55"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)+
  geom_path(data = slingcurves1, aes(x=PHATE1, y=PHATE2), 
            arrow = arrow(type = "open", angle = 40,  length = unit(0.2, "inches")), size =2)+
  theme(legend.position = "none")+
  theme(text=element_text(size = 22))

ggsave(plot = PHATE_SRS_MARS, "FIGURES/IMMERSE paper/SCRIPT 3/PHATE_MARS.png", dpi = 300, height = 7, width = 7)

#PHATE by Sweeney 2019 
pseudotimevaluesexclNA1$`Sweeney 2019`
PHATE_SRS_Sweeney <- ggplot(pseudotimevaluesexclNA1, aes(x=PHATE1, y=PHATE2))+
  geom_point(aes(fill=`Sweeney 2019`),size=7, alpha=1,shape = 21)+
  scale_fill_manual(values = c("adaptive" = "#6cc0e5",
                               "coagulopathic" ="#fbc93d",
                               "inflammopathic" = "#fb4f4f"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  theme(aspect.ratio = 1)+
  geom_path(data = slingcurves1, aes(x=PHATE1, y=PHATE2), 
            arrow = arrow(type = "open", angle = 40,  length = unit(0.2, "inches")), size =2)+
  theme(legend.position = "none")+
  theme(text=element_text(size = 22))

ggsave(plot = PHATE_SRS_Sweeney, "FIGURES/IMMERSE paper/SCRIPT 3/PHATE_Sweeney.png", dpi = 300, height = 7, width = 7)


# centroid plot for PHATE and states
colnames(pseudotimevaluesexclNA1)
PHATE_centroids <- aggregate(pseudotimevaluesexclNA1[,2:3], list(Type = pseudotimevaluesexclNA1$State), median)


#get the centroids for cardiac
colnames(df4.1)
PHATE_cardiac_centroids <- aggregate(df4.1[,734:735], list(Type = df4.1$Cohort_time), median)

centroids <- rbind(PHATE_centroids,PHATE_cardiac_centroids)

centroids$Cohort <- c("Sepsis", "Sepsis","Sepsis","Cardiac", "Cardiac", "Sepsis", "Sepsis", "Sepsis", "Sepsis")
centroids$State <- c("State 1", "State 2", "State 3", "Pre-Surgery", "Post-surgery", "Admission", "Day 3", "Day 5", "CC discharge")

centroids$Type <- c("State", "State", "State", "Surgery", "Surgery", "Clinical timepoint", "Clinical timepoint", "Clinical timepoint", "Clinical timepoint")

centroids$State <- factor(centroids$State, levels = c("State 1", "State 2", "State 3", "Pre-Surgery", "Post-surgery", "CC discharge", "Day 5", "Day 3", "Admission"))


PCA_centroids_plot <- ggplot(centroids, aes(x = PHATE1, y = PHATE2))+
 # geom_density_2d(data = pseudotimevaluesexclNA1, aes(x = PHATE1, y = PHATE2), colour = "grey", alpha = 0.5)+
  geom_point(data = centroids %>% filter(!Type == "Clinical timepoint", 
                                         !Type == "Surgery"), aes(fill = State, shape = Type), colour = "black", alpha = 1, size = 7, stroke = 0.7)+
  geom_point(data = centroids %>% filter(Type == "Clinical timepoint"|
                                           Type == "Surgery"), aes(fill = State, shape = Type), colour = "black", alpha = 0.5, size = 7, stroke = 0.7)+
  scale_fill_manual(values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF",
                               "CC discharge" = "#403891ff",
                               "Pre-Surgery" ="#56A8CBFF",
                               "Post-surgery" = "#DA291CFF",
                               "Admission" = "#efe350ff",
                               "Day 3" = "#f68f46ff",
                               "Day 5" = "#a65c85ff",
                               "CC discharge" = "#403891ff"))+
  scale_shape_manual(values = c("State" = 21, "Surgery" = 24, "Clinical timepoint" = 23))+
  theme_bw()+ 
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom")+
  theme(legend.position = "none")+
  geom_line(data = centroids %>% filter(Type == "Surgery"),aes(group = Type), alpha = 0.5, linetype = "dotdash",linewidth = 1,
            color="black",
            arrow = arrow(type = "closed",
                          length=unit(0.2, "inches"),
                          ends = "first"))+
  geom_line(data = centroids %>% filter(Type == "Clinical timepoint"),aes(group = Type), alpha = 0.5, linetype = "dashed",linewidth = 1,
            color="black",
            arrow = arrow(type = "closed",
                          length=unit(0.2, "inches"),
                          ends = "last"))+
  geom_path(data = centroids%>% filter(Type == "State"),aes(group = Type), alpha = 0.8, linewidth = 1,
            color="black",
            arrow = arrow(type = "closed",
                          length=unit(0.2, "inches"),
                          ends = "last"))+
  theme(aspect.ratio = 1)+ theme(text = element_text(size = 12))    
?arrow
PCA_centroids_plot

ggsave(plot = PCA_centroids_plot, "FIGURES/IMMERSE paper/SCRIPT 3/PHATE_centroids_plot.png", dpi=300, height = 5, width = 5, units = "in" )



  

# plot boxplot of pseduotime on y axis and cluster id on x
bp_time_pseudo <- ggplot(pseudotimevaluesexclNA1, aes(x=Timepoint, y=pseudotime))+
  geom_hline(yintercept = 0.06, colour= "grey")+
  geom_hline(yintercept = 0.12, colour= "grey")+
  geom_violin()+
    geom_quasirandom(aes(fill=State),size=4, alpha=1,shape = 21)+theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(aspect.ratio = 1)+
    scale_fill_manual(values =c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF"))+
    coord_flip()+
  theme_classic()+
  
  theme(legend.position = "none")+
  theme(text = element_text(size = 20))
  
  
ggsave(plot = bp_time_pseudo, "FIGURES/IMMERSE paper/SCRIPT 3/bp_time_pseudo.png", dpi = 300, height = 5, width = 6)
  

######################################################################################################################################
######################################################################################################################################
#plot heatmap of all genes

df_sepsis_hm <- df_sepsis_sling %>%arrange(Lineage1)%>% mutate(State = case_when(hclust3 == "1" ~ "State 1",
                                                                                 hclust3 == "2" ~ "State 2",
                                                                                 hclust3 == "3" ~ "State 3"))

colnames(df_sepsis_hm)
hm_matrix_clus <- df_sepsis_hm %>% dplyr::select(sample_id,16:733) %>% column_to_rownames("sample_id") %>% scale() %>% t()


d <- dist(hm_matrix_clus, method = "euclidean") # Euclidean distance matrix.
H.fit <- hclust(d, method="ward.D")
plot(H.fit) # display dendogram
groups <- cutree(H.fit, k=3) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
plot(H.fit)
rect.hclust(H.fit, k=3, border="red") 
group <- as.data.frame(groups) %>% rownames_to_column("gene_name")#extract grouping variable 
#merge the gene biotype information to this
gene_biotype_df <- df %>% dplyr::select(hgnc_symbol, gene_biotype,gene_ID)
gene_biotype_df <- gene_biotype_df %>% mutate_all(na_if,"") %>%  mutate(gene_name = coalesce(hgnc_symbol, gene_ID)) %>% 
  dplyr::select(-"hgnc_symbol")

# merge gene_biotype_df with gene group
colnames(group)
colnames(gene_biotype_df)
group1 <- merge(group, gene_biotype_df, all.x = T, sort = F)
group1 <- group1  %>%  mutate(`Gene cluster` = case_when(groups == "1" ~ "Cluster 1",
                                    groups == "2" ~ "Cluster 2",
                                    groups == "3" ~ "Cluster 3"))
#save info
write_csv(group1, "Results_csv/IMMERSE paper/SCRIPT 3/718 gene cluster information.csv")
group2 <- group1 %>% dplyr::select(-groups) %>% column_to_rownames("gene_name")

table(group)
head(group2)
## merge cluster group with dataframe

######################################################################################################################################
library("mmand")

hm_matrix <- df_sepsis_hm %>% dplyr::select(sample_id,16:733) %>% column_to_rownames("sample_id")%>% scale()
hm_matrix_smooth = gaussianSmooth(hm_matrix, sigma = 10) %>% t()

col_fun = colorRamp2(c(0,0.5, 1), c("blue", "#EEEEEE", "red"))
col_fun(seq(-1, 1))

#make annotation
colnames(df_sepsis_hm)
anno <- df_sepsis_hm %>% dplyr::select(Pseudotime=Lineage1, `Clinical timepoint`=Cohort_time, State, SRS_davenport,SRS_extended, MARS, `Sweeney 2019`)
colnames(anno)
#create colour scheme for pseudotime´
n = nrow(anno)
min_v = min(anno$Pseudotime)
max_v = max(anno$Pseudotime)
Var = circlize::colorRamp2(seq(min_v, max_v, length = n), hcl.colors(n,"Spectral"))
#create annotation
ha = HeatmapAnnotation(df = anno, 
                       col = list(Pseudotime = Var,
                                  `Clinical timepoint` = c("Sepsis_T0" = "#efe350ff", "Sepsis_T1" = "#f68f46ff","Sepsis_T2" = "#a65c85ff","Sepsis_T3" = "#403891ff"),
                                  State = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF"),
                                  SRS_davenport = c("SRS1" = "#810103",
                                                    "SRS2" ="#4d7faf",
                                                    "SRS3" = "#4d7faf"),
                                  SRS_extended = c("SRS1" = "#810103",
                                                    "SRS2" ="#4d7faf",
                                                    "SRS3" = "#010089"),
                                  MARS = c("Mars1" = "#c7deb2",
                                           "Mars2" ="#2e64aa",
                                           "Mars3" = "#9dacd5",
                                           "Mars4" = "#82ba55"),
                                  `Sweeney 2019` = c("adaptive" = "#6cc0e5",
                                                    "coagulopathic" ="#fbc93d",
                                                     "inflammopathic" = "#fb4f4f")),
                       
                       which = "col",
                       gp = gpar(col = "black"))


# make the column annotation
group2 <- group2 %>% dplyr::select(-gene_ID)
ha_row = HeatmapAnnotation(df = group2, 
                           col = list(`Gene cluster` = c("Cluster 1" =  "#FF756D",
                                                 "Cluster 2" = "#F6B4C3", 
                                                 "Cluster 3" = "#AAD6FA"),
                                      gene_biotype = c("artifact" = "#E31A1C", 
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
                                                       "unprocessed_pseudogene" = "steelblue4")),
                           which = "row")
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

png("FIGURES/IMMERSE paper/SCRIPT 3/heatmap.png",width=15,height=8,units="in",res = 800)
hm
dev.off()


##########################################################################################
# Show contrtibuting factors to the variation
library(pvca)
library(Biobase)
pct_threshold <- 0.6
#Add metadata to the sepsis patients to include in this analysis. 
merged_baseline <- read_csv("clinical_csv/merged baseline.csv") %>% dplyr::select(case_id, age_years, sex)
colnames(merged_baseline)[1:20]
colnames(df_sepsis_sling)[1:20]
df_sepsis_sling <- merge(merged_baseline, df_sepsis_sling, all.y = T)
colnames(df_sepsis_sling)
test_m <- df_sepsis_sling%>% 
  dplyr::select(18:735) %>% 
  as.matrix()
df_sepsis_sling$case_id

test_m_anno <- df_sepsis_sling %>% 
  dplyr::select("Timepoint","Infection" = "Infection source2","outcome"= "outcome_hospital_survival", "State" = "hclust3",
                "Age" ="age_years", "sex")

metadata_m <- data.frame(LabelDescription = 
                           c("Clinical timepoints",
                             "Infection source",
                             "outcome",
                             "State",
                             "Age",
                             "sex"
                             ),
                            row.names=c("Timepoint","Infection","outcome", "State",
                                        "Age",
                                        "sex"))


all(rownames(test_m_anno)==colnames(t(test_m)))

test_m_annon1 <- new("AnnotatedDataFrame",
                     data = test_m_anno,
                     varMetadata = metadata_m)

batch.factors <- c("Timepoint","Infection","outcome", "State",
                   "Age",
                   "sex")

minimalSet <- ExpressionSet(assayData = t(test_m),
                            phenoData = test_m_annon1)

bpvaObj <- pvcaBatchAssess(minimalSet, batch.factors, pct_threshold)
colnames(barplot_var)
barplot_var <- as.data.frame(bpvaObj$dat)
new <- c(bpvaObj$label)
old <- c(colnames(barplot_var))
barplot_var <- barplot_var %>% rename_with(~new, all_of(old))

colnames(barplot_var)
library(wesanderson)
class(barplot_var$`Timepoint:Infection`)
barplot_variance <- barplot_var %>% dplyr::select(6,14,18,19,20,21) %>% 
  round(4) %>% 
  mutate_if(is.numeric, ~ . *100)%>% 
  gather("Variable", "Proportion", 1:5) %>% 
  ggplot(aes(x = reorder(Variable, - Proportion), y = Proportion, label = Proportion))+
  geom_bar(stat="identity", aes(fill = Variable), colour= "black")+
  geom_text(aes(label = Proportion), vjust=-0.3, size = 6)+
  theme_classic()+
  scale_fill_manual(values = c("#94C9A9", 
                               "#6495ED",
                               "#F6AE2D",
                               "#BC3908",
                               "#D0CFEC"))+
  theme(text = element_text(size = 20))+
  theme(legend.position = "None")+
  xlab("Variable")+
  ylim(c(0,40))

ggsave(plot = barplot_variance, "Figures/IMMERSE paper/SCRIPT 3/variance plot.png", width = 6, height =6, dpi = 300)

##########################################################################################
##########################################################################################



##########################################################################################
##########################################################################################

#name the clusters right in this data frame
df_sepsis_sling <- df_sepsis_sling %>%   mutate(State = case_when(hclust3 == "1" ~ "State 1",
                                               hclust3 == "2" ~ "State 2",
                                               hclust3 == "3" ~ "State 3"))
#########################################################################################################
#########################################################################################################
# Make an alluvial plot
library(scales)
library(ggalluvial)
colnames(df_sepsis_sling)[1:15]


df_sepsis_sling <- df_sepsis_sling %>%   mutate(Timepoint2 = case_when(Timepoint == "T0" ~ "Admission",
                                                                  Timepoint == "T1" ~ "Day 3",
                                                                  Timepoint == "T2" ~ "Day 5",
                                                                  Timepoint == "T3" ~ "CC discharge"))
df_sepsis_sling$Timepoint2 <- factor(df_sepsis_sling$Timepoint2, levels = c('Admission', 'Day 3', 'Day 5', "CC discharge"))


n_id <- length(unique(df_sepsis_sling$case_id))
alluvial2 <- ggplot(df_sepsis_sling, 
       aes(x = Timepoint2, stratum = State, alluvium = case_id,
           fill = State, label = State)) +
  scale_fill_manual(values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF"))+
 # scale_y_continuous(label = scales::percent_format(scale = 100 / n_id)) +
  geom_flow(stat = "flow", knot.pos = 1/4, aes.flow = "forward",
            color = "darkgray", na.rm = FALSE) + 
  geom_stratum() +
  theme(legend.position = "bottom") +
  geom_text(stat = "stratum",
            aes(label = percent(after_stat(prop), accuracy = .1 )), size = 3)+
  theme_classic()+
  ylab("Number of samples")+
  xlab("Timepoint")+
  theme(text = element_text(size = 16))  

ggsave(plot = alluvial2, "FIGURES/IMMERSE paper/SCRIPT 3/patient_alluvial2.png", height = 5, width = 7, dpi = 300)

#plot individual patients over time.
df_sepsis_sling$State
#read in outcomes
df_outcome <- read_csv("clinical_csv/df_outcome.csv")

df_sepsis_sling1 <- merge(df_sepsis_sling, df_outcome)
patient_PT <- ggplot(df_sepsis_sling, aes(x = Lineage1, y = case_id))+ 
  geom_line(aes(group = case_id),color="black")+
  geom_point(aes(fill = State), pch = 21, size = 3) +
  scale_fill_manual(values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF"))+
  theme_classic()+
  ylab("Patient ID")+
  xlab("Sepsis immune trajectory (Pseudotime)")

df_sepsis_sling1$outcome_hospital_survival
ggsave(plot = patient_PT, "FIGURES/IMMERSE paper/SCRIPT 3/patient_dotplot_pseudotime.png", height = 7, width = 7, dpi = 300)



# save the metadata with sates 
#select data 
metadata_states <- df_sepsis_sling %>% dplyr:: select(1:12,Lineage1,State,PHATE1,PHATE2,DC1,DC2)
colnames(df_sepsis_sling)
#save
write_csv(metadata_states, "Metadata/metadata_states_Phate.csv")



