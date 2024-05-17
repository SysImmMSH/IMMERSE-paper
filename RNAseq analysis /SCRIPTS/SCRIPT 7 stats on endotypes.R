# compare the SRS, mars and sweeney to pre surgery. 


#read in packages
library(tximport)
library(tidyverse)
library(GenomicFeatures)
library(DelayedArray) 
library(DESeq2)
library(magrittr)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(apeglm)
library(topGO)
library(pcaExplorer)
library(GeneTonic)

################################################################################################
################################################################################################
################################################################################################
# 1) read in all quant files and with Deseq2. Creates count matrix and stats. 

#load meta data 
samples <- read_csv("Metadata/metadata.csv") %>% #filter(Cohort_time == "Cardiac_T0" |
  #   Cohort_time== "Cardiac_T1" |
  # Cohort_time == "Sepsis_T0") %>% 
  filter(COVID_status == "Negative") %>% 
  filter(!Library_cohort == "Cardiac_lib2")%>%
  filter(!case_id == "SC09", #removed as IVIG
         !case_id == "SC84", #removed as insisted on it 
         !case_id == "SC27") %>%  #removed as screen failure, sofa >2  
  as.data.frame()

#load metadata from states
states_meta <- read_csv("Metadata/metadata_states_Phate.csv") %>% dplyr::select(sample_id, State,Lineage1,SRS_davenport,
                                                                                SRS_extended)

#merge
samples <- merge(samples, states_meta, all.x = TRUE)

#read in sweeney data 
#Note, have manually removed classifcation for the cardiac group and saved csv. 
sweeney <- read_csv("Results_csv/IMMERSE paper/SCRIPT 1/sweeney_groups.csv") %>% filter(!str_detect(sample_id, 'CC')) %>% 
  dplyr::select(sample_id,  Sweeney = `Sweeney 2019`)
#merge
samples <- merge(samples, sweeney, by = "sample_id", all.x = TRUE)

#read in mars groupings 
mars_groupings <- read_csv("Results_csv/IMMERSE paper/SCRIPT 1/mars_groupings.csv")%>% filter(!str_detect(sample_id, 'CC')) %>% 
  dplyr::select(sample_id,  MARS)
#merge
samples <- merge(samples, mars_groupings, by = "sample_id", all.x = TRUE)

#add in cardiac surgery to fill in NAs 
samples <- samples %>%  dplyr::mutate(SRS_davenport = coalesce(SRS_davenport, Cohort_time)) 
colnames(states_meta)

samples <- samples %>%  mutate(SRS_davenport = case_when(SRS_davenport == "SRS3" ~ "SRS2",
                                                         SRS_davenport == "SRS2" ~ "SRS2",
                                                         SRS_davenport == "SRS1" ~ "SRS1",
                                                         SRS_davenport == "Cardiac_T0" ~ "Cardiac_T0"))

samples <- samples %>%  dplyr::mutate(SRS_extended = coalesce(SRS_extended, Cohort_time)) 
colnames(states_meta)

samples <- samples %>%  dplyr::mutate(MARS = coalesce(MARS, Cohort_time)) 
colnames(states_meta)

samples <- samples %>%  dplyr::mutate(`Sweeney` = coalesce(`Sweeney`, Cohort_time)) 
colnames(states_meta)

table(samples$Cohort_time)

samples <- samples %>% filter(!Cohort_time == "Cardiac_T1")

table(samples$`Sweeney`)

####### directory to salmon quant files
#directory of quant.sf files
dir <- "quantfiles"
#list them/check they are there
list.files(file.path(dir))

files <- file.path(dir, samples$name)

names(files) <- samples$sample_id
all(file.exists(files))
?names
#gene annotation imput 
txdb <- makeTxDbFromGFF("gta file/gencode.v38.annotation.gtf.gz")



#select columns required
columns(txdb)
k <- keys(txdb, "GENEID")
tx2gene <- select(txdb, keys = k, keytype = 'GENEID', columns = 'TXNAME')
head(tx2gene)


res <- AnnotationDbi::select(txdb, k, "TXNAME", "GENEID")
tx2gene <- res[,2:1]
head(tx2gene)

#use tximport to read our count files. Since we are using the same gtf file, the versions of the transcripts will be the same, hence the argument ignoreTxVersion = FALSE. If not, we set it to TRUE
txi <- tximport(files = files, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = FALSE, ignoreAfterBar = TRUE)
summary(txi)
head(txi)

#########################################################################################
#########################################################################################
#########################################################################################
# SRS_Davenport first.

#create deseq object
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ library_id + SRS_davenport)



nrow(ddsTxi) #60230 genes measured
ddsTxi <- ddsTxi[rowSums(counts(ddsTxi) >= 10) >= 20]
nrow(ddsTxi)



dds <- DESeq(ddsTxi)



#add symbol to the rowData. Note this adds a dataframe in that section. So when using downstream will need to subset. 
#get the gene name
rownames(dds) <- substr(rownames(dds), 1, 15)
rowData(dds)$SYMBOL1 <- get_annotation_orgdb(dds, "org.Hs.eg.db", "ENSEMBL")
#add to dds 
rowData(dds)$SYMBOL <- as.factor(rowData(dds)$SYMBOL1$gene_name)

########################################################################################################################
########################################################################################################################
library(viridis)
#Gene tonic

#SRS1 vs cardiac T0
SRS1vsCCT0_GT <- results(dds, contrast = c("SRS_davenport", "SRS1", "Cardiac_T0"),lfcThreshold = 1, alpha = 0.05 )
SRS1vsCCT0_GT$SYMBOL <- rowData(dds)$SYMBOL

#make the information for state1vsCCT0_GT
de_symbols_SRS1vsCCT0_GT <- deseqresult2df(SRS1vsCCT0_GT, FDR = 0.05)$SYMBOL
bg_ids <- rowData(dds)$SYMBOL[rowSums(counts(dds)) > 0]


topgoDE_SRS1vsCCT0_GT <-
  pcaExplorer::topGOtable(de_symbols_SRS1vsCCT0_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_SRS1vsCCT0_GT <- shake_topGOtableResult(topgoDE_SRS1vsCCT0_GT)
colnames(topgoDE_Cohort_time_Cardiac_T1_vs_Cardiac_T0)

#create the annotation dataframe
anno_df <- data.frame(
  gene_id = rownames(dds),
  gene_name = mapIds(org.Hs.eg.db, keys = rownames(dds), column = "SYMBOL", keytype = "ENSEMBL"),
  stringsAsFactors = FALSE,
  row.names = rownames(dds)
)
head(anno_df)

#create the aggreegated scores
res_enrich_SRS1vsCCT0_GT <- get_aggrscores(res_enrich = res_enrich_SRS1vsCCT0_GT,
                                             res_de = SRS1vsCCT0_GT,
                                             annotation_obj = anno_df,
                                             aggrfun = mean)
#plot top 10 most significant pathways for SRS2vs cardiac T0
res_enrich_SRS1vsCCT0_GT <- res_enrich_SRS1vsCCT0_GT[1:10,] 

SRS1vsCCT0_BP <- ggplot(res_enrich_SRS1vsCCT0_GT,aes(x=z_score, y =reorder(gs_description, +z_score)))+
  geom_bar(stat = "identity", aes(fill=z_score), colour = "black")+
  scale_fill_viridis(option = "D")+
  theme_classic()+
  theme(text=element_text(size =14))

#save
ggsave(plot = SRS1vsCCT0_BP, "FIGURES/IMMERSE paper/SCRIPT 7/SRS1vsCCT0_BP.png", width = 10, height = 8, dpi = 300 )


##########################################################################################
##########################################################################################

#SRS2 vs cardiac T0

SRS2vsCCT0_GT <- results(dds, contrast = c("SRS_davenport", "SRS2", "Cardiac_T0"),lfcThreshold = 1, alpha = 0.05 )
SRS2vsCCT0_GT$SYMBOL <- rowData(dds)$SYMBOL

#make the information for state1vsCCT0_GT
de_symbols_SRS2vsCCT0_GT <- deseqresult2df(SRS2vsCCT0_GT, FDR = 0.05)$SYMBOL
bg_ids <- rowData(dds)$SYMBOL[rowSums(counts(dds)) > 0]


topgoDE_SRS2vsCCT0_GT <-
  pcaExplorer::topGOtable(de_symbols_SRS2vsCCT0_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_SRS2vsCCT0_GT <- shake_topGOtableResult(topgoDE_SRS2vsCCT0_GT)
colnames(topgoDE_Cohort_time_Cardiac_T1_vs_Cardiac_T0)

#create the annotation dataframe
anno_df <- data.frame(
  gene_id = rownames(dds),
  gene_name = mapIds(org.Hs.eg.db, keys = rownames(dds), column = "SYMBOL", keytype = "ENSEMBL"),
  stringsAsFactors = FALSE,
  row.names = rownames(dds)
)
head(anno_df)

#create the aggreegated scores
res_enrich_SRS2vsCCT0_GT <- get_aggrscores(res_enrich = res_enrich_SRS2vsCCT0_GT,
                                           res_de = SRS2vsCCT0_GT,
                                           annotation_obj = anno_df,
                                           aggrfun = mean)

#plot top 10 most significant pathways for SRS2vs cardiac T0
res_enrich_SRS2vsCCT0_GT <- res_enrich_SRS2vsCCT0_GT[1:10,] 

SRS2vsCCT0_BP <- ggplot(res_enrich_SRS2vsCCT0_GT,aes(x=z_score, y =reorder(gs_description, +z_score)))+
         geom_bar(stat = "identity", aes(fill=z_score), colour = "black")+
  scale_fill_viridis(option = "D")+
  theme_classic()+
  theme(text=element_text(size =14))

#save
ggsave(plot = SRS2vsCCT0_BP, "FIGURES/IMMERSE paper/SCRIPT 7/SRS2vsCCT0_BP.png", width = 10, height = 8, dpi = 300 )







#########################################################################################
#########################################################################################
#########################################################################################
# SRS_Extended

#create deseq object
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ library_id + SRS_extended)



nrow(ddsTxi) #60230 genes measured
ddsTxi <- ddsTxi[rowSums(counts(ddsTxi) >= 10) >= 20]
nrow(ddsTxi)



dds <- DESeq(ddsTxi)



#add symbol to the rowData. Note this adds a dataframe in that section. So when using downstream will need to subset. 
#get the gene name
rownames(dds) <- substr(rownames(dds), 1, 15)
rowData(dds)$SYMBOL1 <- get_annotation_orgdb(dds, "org.Hs.eg.db", "ENSEMBL")
#add to dds 
rowData(dds)$SYMBOL <- as.factor(rowData(dds)$SYMBOL1$gene_name)

########################################################################################################################
########################################################################################################################

#Gene tonic

#SRS1 vs cardiac T0
SRS1vsCCT0_GT <- results(dds, contrast = c("SRS_extended", "SRS1", "Cardiac_T0"),lfcThreshold = 1, alpha = 0.05 )
SRS1vsCCT0_GT$SYMBOL <- rowData(dds)$SYMBOL

#make the information for state1vsCCT0_GT
de_symbols_SRS1vsCCT0_GT <- deseqresult2df(SRS1vsCCT0_GT, FDR = 0.05)$SYMBOL
bg_ids <- rowData(dds)$SYMBOL[rowSums(counts(dds)) > 0]


topgoDE_SRS1vsCCT0_GT <-
  pcaExplorer::topGOtable(de_symbols_SRS1vsCCT0_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_SRS1vsCCT0_GT <- shake_topGOtableResult(topgoDE_SRS1vsCCT0_GT)
colnames(topgoDE_Cohort_time_Cardiac_T1_vs_Cardiac_T0)

#create the annotation dataframe
anno_df <- data.frame(
  gene_id = rownames(dds),
  gene_name = mapIds(org.Hs.eg.db, keys = rownames(dds), column = "SYMBOL", keytype = "ENSEMBL"),
  stringsAsFactors = FALSE,
  row.names = rownames(dds)
)
head(anno_df)

#create the aggreegated scores
res_enrich_SRS1vsCCT0_GT <- get_aggrscores(res_enrich = res_enrich_SRS1vsCCT0_GT,
                                           res_de = SRS1vsCCT0_GT,
                                           annotation_obj = anno_df,
                                           aggrfun = mean)
#plot top 10 most significant pathways for SRS2vs cardiac T0
res_enrich_SRS1vsCCT0_GT <- res_enrich_SRS1vsCCT0_GT[1:10,] 

SRS1vsCCT0_BP <- ggplot(res_enrich_SRS1vsCCT0_GT,aes(x=z_score, y =reorder(gs_description, +z_score)))+
  geom_bar(stat = "identity", aes(fill=z_score), colour = "black")+
  scale_fill_viridis(option = "D")+
  theme_classic()+
  theme(text=element_text(size =14))

#save
ggsave(plot = SRS1vsCCT0_BP, "FIGURES/IMMERSE paper/SCRIPT 7/SRS1_extendedvsCCT0_BP.png", width = 10, height = 8, dpi = 300 )


##########################################################################################
##########################################################################################

#SRS2 vs cardiac T0

SRS2vsCCT0_GT <- results(dds, contrast = c("SRS_extended", "SRS2", "Cardiac_T0"),lfcThreshold = 1, alpha = 0.05 )
SRS2vsCCT0_GT$SYMBOL <- rowData(dds)$SYMBOL

#make the information for state1vsCCT0_GT
de_symbols_SRS2vsCCT0_GT <- deseqresult2df(SRS2vsCCT0_GT, FDR = 0.05)$SYMBOL
bg_ids <- rowData(dds)$SYMBOL[rowSums(counts(dds)) > 0]


topgoDE_SRS2vsCCT0_GT <-
  pcaExplorer::topGOtable(de_symbols_SRS2vsCCT0_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_SRS2vsCCT0_GT <- shake_topGOtableResult(topgoDE_SRS2vsCCT0_GT)
colnames(topgoDE_Cohort_time_Cardiac_T1_vs_Cardiac_T0)

#create the annotation dataframe
anno_df <- data.frame(
  gene_id = rownames(dds),
  gene_name = mapIds(org.Hs.eg.db, keys = rownames(dds), column = "SYMBOL", keytype = "ENSEMBL"),
  stringsAsFactors = FALSE,
  row.names = rownames(dds)
)
head(anno_df)

#create the aggreegated scores
res_enrich_SRS2vsCCT0_GT <- get_aggrscores(res_enrich = res_enrich_SRS2vsCCT0_GT,
                                           res_de = SRS2vsCCT0_GT,
                                           annotation_obj = anno_df,
                                           aggrfun = mean)

#plot top 10 most significant pathways for SRS2vs cardiac T0
res_enrich_SRS2vsCCT0_GT <- res_enrich_SRS2vsCCT0_GT[1:10,] 

SRS2vsCCT0_BP <- ggplot(res_enrich_SRS2vsCCT0_GT,aes(x=z_score, y =reorder(gs_description, +z_score)))+
  geom_bar(stat = "identity", aes(fill=z_score), colour = "black")+
  scale_fill_viridis(option = "D")+
  theme_classic()+
  theme(text=element_text(size =14))

#save
ggsave(plot = SRS2vsCCT0_BP, "FIGURES/IMMERSE paper/SCRIPT 7/SRS2_extendedvsCCT0_BP.png", width = 10, height = 8, dpi = 300 )


##########################################################################################
##########################################################################################

#SRS3 vs cardiac T0

SRS3vsCCT0_GT <- results(dds, contrast = c("SRS_extended", "SRS3", "Cardiac_T0"),lfcThreshold = 1, alpha = 0.05 )
SRS3vsCCT0_GT$SYMBOL <- rowData(dds)$SYMBOL

#make the information for state1vsCCT0_GT
de_symbols_SRS3vsCCT0_GT <- deseqresult2df(SRS3vsCCT0_GT, FDR = 0.05)$SYMBOL
bg_ids <- rowData(dds)$SYMBOL[rowSums(counts(dds)) > 0]


topgoDE_SRS3vsCCT0_GT <-
  pcaExplorer::topGOtable(de_symbols_SRS3vsCCT0_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_SRS3vsCCT0_GT <- shake_topGOtableResult(topgoDE_SRS3vsCCT0_GT)
colnames(topgoDE_Cohort_time_Cardiac_T1_vs_Cardiac_T0)

#create the annotation dataframe
anno_df <- data.frame(
  gene_id = rownames(dds),
  gene_name = mapIds(org.Hs.eg.db, keys = rownames(dds), column = "SYMBOL", keytype = "ENSEMBL"),
  stringsAsFactors = FALSE,
  row.names = rownames(dds)
)
head(anno_df)

#create the aggreegated scores
res_enrich_SRS3vsCCT0_GT <- get_aggrscores(res_enrich = res_enrich_SRS3vsCCT0_GT,
                                           res_de = SRS3vsCCT0_GT,
                                           annotation_obj = anno_df,
                                           aggrfun = mean)

#plot top 10 most significant pathways for SRS3vs cardiac T0
res_enrich_SRS3vsCCT0_GT <- res_enrich_SRS3vsCCT0_GT[1:10,] 

SRS3vsCCT0_BP <- ggplot(res_enrich_SRS3vsCCT0_GT,aes(x=z_score, y =reorder(gs_description, +z_score)))+
  geom_bar(stat = "identity", aes(fill=z_score), colour = "black")+
  scale_fill_viridis(option = "D")+
  theme_classic()+
  theme(text=element_text(size =14))

#save
ggsave(plot = SRS3vsCCT0_BP, "FIGURES/IMMERSE paper/SCRIPT 7/SRS3_extendedvsCCT0_BP.png", width = 10, height = 8, dpi = 300 )



#########################################################################################
#########################################################################################
#########################################################################################
# Sweeney
colnames(samples)
#create deseq object
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ library_id + Sweeney)



nrow(ddsTxi) #60230 genes measured
ddsTxi <- ddsTxi[rowSums(counts(ddsTxi) >= 10) >= 20]
nrow(ddsTxi)



dds <- DESeq(ddsTxi)



#add symbol to the rowData. Note this adds a dataframe in that section. So when using downstream will need to subset. 
#get the gene name
rownames(dds) <- substr(rownames(dds), 1, 15)
rowData(dds)$SYMBOL1 <- get_annotation_orgdb(dds, "org.Hs.eg.db", "ENSEMBL")
#add to dds 
rowData(dds)$SYMBOL <- as.factor(rowData(dds)$SYMBOL1$gene_name)

########################################################################################################################
########################################################################################################################

#Gene tonic
table(samples$Sweeney)
#SRS1 vs cardiac T0
InvsCCT0_GT <- results(dds, contrast = c("Sweeney", "inflammopathic", "Cardiac_T0"),lfcThreshold = 1, alpha = 0.05 )
InvsCCT0_GT$SYMBOL <- rowData(dds)$SYMBOL

#make the information for state1vsCCT0_GT
de_symbols_InvsCCT0_GT <- deseqresult2df(InvsCCT0_GT, FDR = 0.05)$SYMBOL
bg_ids <- rowData(dds)$SYMBOL[rowSums(counts(dds)) > 0]


topgoDE_InvsCCT0_GT <-
  pcaExplorer::topGOtable(de_symbols_InvsCCT0_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_InvsCCT0_GT <- shake_topGOtableResult(topgoDE_InvsCCT0_GT)
colnames(topgoDE_Cohort_time_Cardiac_T1_vs_Cardiac_T0)

#create the annotation dataframe
anno_df <- data.frame(
  gene_id = rownames(dds),
  gene_name = mapIds(org.Hs.eg.db, keys = rownames(dds), column = "SYMBOL", keytype = "ENSEMBL"),
  stringsAsFactors = FALSE,
  row.names = rownames(dds)
)
head(anno_df)

#create the aggreegated scores
res_enrich_InvsCCT0_GT <- get_aggrscores(res_enrich = res_enrich_InvsCCT0_GT,
                                           res_de = InvsCCT0_GT,
                                           annotation_obj = anno_df,
                                           aggrfun = mean)
#plot top 10 most significant pathways for SRS2vs cardiac T0
res_enrich_InvsCCT0_GT <- res_enrich_InvsCCT0_GT[1:10,] 

InvsCCT0_BP <- ggplot(res_enrich_InvsCCT0_GT,aes(x=z_score, y =reorder(gs_description, +z_score)))+
  geom_bar(stat = "identity", aes(fill=z_score), colour = "black")+
  scale_fill_viridis(option = "D")+
  theme_classic()+
  theme(text=element_text(size =14))

#save
ggsave(plot = InvsCCT0_BP, "FIGURES/IMMERSE paper/SCRIPT 7/In_extendedvsCCT0_BP.png", width = 10, height = 8, dpi = 300 )


##########################################################################################
##########################################################################################
colnames(test)
ggplot(test, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05))
#adap vs cardiac T0

adapvsCCT0_GT <- results(dds, contrast = c("Sweeney", "adaptive", "Cardiac_T0"),lfcThreshold = 0.3, alpha = 0.05 )
adapvsCCT0_GT$SYMBOL <- rowData(dds)$SYMBOL


#make the information for state1vsCCT0_GT
de_symbols_adapvsCCT0_GT <- deseqresult2df(adapvsCCT0_GT, FDR = 0.05)$SYMBOL
bg_ids <- rowData(dds)$SYMBOL[rowSums(counts(dds)) > 0]


topgoDE_adapvsCCT0_GT <-
  pcaExplorer::topGOtable(de_symbols_adapvsCCT0_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_adapvsCCT0_GT <- shake_topGOtableResult(topgoDE_adapvsCCT0_GT)
colnames(topgoDE_Cohort_time_Cardiac_T1_vs_Cardiac_T0)

#create the annotation dataframe
anno_df <- data.frame(
  gene_id = rownames(dds),
  gene_name = mapIds(org.Hs.eg.db, keys = rownames(dds), column = "SYMBOL", keytype = "ENSEMBL"),
  stringsAsFactors = FALSE,
  row.names = rownames(dds)
)
head(anno_df)

#create the aggreegated scores
res_enrich_adapvsCCT0_GT <- get_aggrscores(res_enrich = res_enrich_adapvsCCT0_GT,
                                           res_de = adapvsCCT0_GT,
                                           annotation_obj = anno_df,
                                           aggrfun = mean)

#plot top 10 most significant pathways for adapvs cardiac T0
res_enrich_adapvsCCT0_GT <- res_enrich_adapvsCCT0_GT[1:10,] 

adapvsCCT0_BP <- ggplot(res_enrich_adapvsCCT0_GT,aes(x=z_score, y =reorder(gs_description, +z_score)))+
  geom_bar(stat = "identity", aes(fill=z_score), colour = "black")+
  scale_fill_viridis(option = "D")+
  theme_classic()+
  theme(text=element_text(size =14))

#save
ggsave(plot = adapvsCCT0_BP, "FIGURES/IMMERSE paper/SCRIPT 7/adap_vsCCT0_BP.png", width = 10, height = 8, dpi = 300 )


##########################################################################################
##########################################################################################

#coag vs cardiac T0

coagvsCCT0_GT <- results(dds, contrast = c("Sweeney", "coagulopathic", "Cardiac_T0"),lfcThreshold = 1, alpha = 0.05 )
coagvsCCT0_GT$SYMBOL <- rowData(dds)$SYMBOL

#make the information for state1vsCCT0_GT
de_symbols_coagvsCCT0_GT <- deseqresult2df(coagvsCCT0_GT, FDR = 0.05)$SYMBOL
bg_ids <- rowData(dds)$SYMBOL[rowSums(counts(dds)) > 0]


topgoDE_coagvsCCT0_GT <-
  pcaExplorer::topGOtable(de_symbols_coagvsCCT0_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_coagvsCCT0_GT <- shake_topGOtableResult(topgoDE_coagvsCCT0_GT)
colnames(topgoDE_Cohort_time_Cardiac_T1_vs_Cardiac_T0)

#create the annotation dataframe
anno_df <- data.frame(
  gene_id = rownames(dds),
  gene_name = mapIds(org.Hs.eg.db, keys = rownames(dds), column = "SYMBOL", keytype = "ENSEMBL"),
  stringsAsFactors = FALSE,
  row.names = rownames(dds)
)
head(anno_df)

#create the aggreegated scores
res_enrich_coagvsCCT0_GT <- get_aggrscores(res_enrich = res_enrich_coagvsCCT0_GT,
                                           res_de = coagvsCCT0_GT,
                                           annotation_obj = anno_df,
                                           aggrfun = mean)

#plot top 10 most significant pathways for coagvs cardiac T0
res_enrich_coagvsCCT0_GT <- res_enrich_coagvsCCT0_GT[1:10,] 

coagvsCCT0_BP <- ggplot(res_enrich_coagvsCCT0_GT,aes(x=z_score, y =reorder(gs_description, +z_score)))+
  geom_bar(stat = "identity", aes(fill=z_score), colour = "black")+
  scale_fill_viridis(option = "D")+
  theme_classic()+
  theme(text=element_text(size =14))

#save
ggsave(plot = coagvsCCT0_BP, "FIGURES/IMMERSE paper/SCRIPT 7/coag_extendedvsCCT0_BP.png", width = 10, height = 8, dpi = 300 )



#########################################################################################
#########################################################################################
#########################################################################################
# mARS
colnames(samples)
table(samples$MARS)
#create deseq object
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ library_id + MARS)



nrow(ddsTxi) #60230 genes measured
ddsTxi <- ddsTxi[rowSums(counts(ddsTxi) >= 10) >= 20]
nrow(ddsTxi)



dds <- DESeq(ddsTxi)



#add symbol to the rowData. Note this adds a dataframe in that section. So when using downstream will need to subset. 
#get the gene name
rownames(dds) <- substr(rownames(dds), 1, 15)
rowData(dds)$SYMBOL1 <- get_annotation_orgdb(dds, "org.Hs.eg.db", "ENSEMBL")
#add to dds 
rowData(dds)$SYMBOL <- as.factor(rowData(dds)$SYMBOL1$gene_name)

########################################################################################################################
########################################################################################################################

#Gene tonic
table(samples$Sweeney)
#SRS1 vs cardiac T0
Mars1vsCCT0_GT <- results(dds, contrast = c("MARS", "Mars1", "Cardiac_T0"),lfcThreshold = 1, alpha = 0.05 )
Mars1vsCCT0_GT$SYMBOL <- rowData(dds)$SYMBOL

#make the Mars1formation for state1vsCCT0_GT
de_symbols_Mars1vsCCT0_GT <- deseqresult2df(Mars1vsCCT0_GT, FDR = 0.05)$SYMBOL
bg_ids <- rowData(dds)$SYMBOL[rowSums(counts(dds)) > 0]


topgoDE_Mars1vsCCT0_GT <-
  pcaExplorer::topGOtable(de_symbols_Mars1vsCCT0_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use Mars1 genetonic
res_enrich_Mars1vsCCT0_GT <- shake_topGOtableResult(topgoDE_Mars1vsCCT0_GT)
colnames(topgoDE_Cohort_time_Cardiac_T1_vs_Cardiac_T0)

#create the annotation dataframe
anno_df <- data.frame(
  gene_id = rownames(dds),
  gene_name = mapIds(org.Hs.eg.db, keys = rownames(dds), column = "SYMBOL", keytype = "ENSEMBL"),
  strMars1gsAsFactors = FALSE,
  row.names = rownames(dds)
)
head(anno_df)

#create the aggreegated scores
res_enrich_Mars1vsCCT0_GT <- get_aggrscores(res_enrich = res_enrich_Mars1vsCCT0_GT,
                                         res_de = Mars1vsCCT0_GT,
                                         annotation_obj = anno_df,
                                         aggrfun = mean)
#plot top 10 most significant pathways for SRS2vs cardiac T0
res_enrich_Mars1vsCCT0_GT <- res_enrich_Mars1vsCCT0_GT[1:10,] 

Mars1vsCCT0_BP <- ggplot(res_enrich_Mars1vsCCT0_GT,aes(x=z_score, y =reorder(gs_description, +z_score)))+
  geom_bar(stat = "identity", aes(fill=z_score), colour = "black")+
  scale_fill_viridis(option = "D")+
  theme_classic()+
  theme(text=element_text(size =14))

#save
ggsave(plot = Mars1vsCCT0_BP, "FIGURES/IMMERSE paper/SCRIPT 7/Mars1_vsCCT0_BP.png", width = 10, height = 8, dpi = 300 )


##########################################################################################
##########################################################################################
colnames(test)

#adap vs cardiac T0

Mars2vsCCT0_GT <- results(dds, contrast = c("MARS", "Mars2", "Cardiac_T0"),lfcThreshold = 0.3, alpha = 0.05 )
Mars2vsCCT0_GT$SYMBOL <- rowData(dds)$SYMBOL


#make the information for state1vsCCT0_GT
de_symbols_Mars2vsCCT0_GT <- deseqresult2df(Mars2vsCCT0_GT, FDR = 0.05)$SYMBOL
bg_ids <- rowData(dds)$SYMBOL[rowSums(counts(dds)) > 0]


topgoDE_Mars2vsCCT0_GT <-
  pcaExplorer::topGOtable(de_symbols_Mars2vsCCT0_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_Mars2vsCCT0_GT <- shake_topGOtableResult(topgoDE_Mars2vsCCT0_GT)
colnames(topgoDE_Cohort_time_Cardiac_T1_vs_Cardiac_T0)

#create the annotation dataframe
anno_df <- data.frame(
  gene_id = rownames(dds),
  gene_name = mapIds(org.Hs.eg.db, keys = rownames(dds), column = "SYMBOL", keytype = "ENSEMBL"),
  stringsAsFactors = FALSE,
  row.names = rownames(dds)
)
head(anno_df)

#create the aggreegated scores
res_enrich_Mars2vsCCT0_GT <- get_aggrscores(res_enrich = res_enrich_Mars2vsCCT0_GT,
                                           res_de = Mars2vsCCT0_GT,
                                           annotation_obj = anno_df,
                                           aggrfun = mean)

#plot top 10 most significant pathways for Mars2vs cardiac T0
res_enrich_Mars2vsCCT0_GT <- res_enrich_Mars2vsCCT0_GT[1:10,] 

Mars2vsCCT0_BP <- ggplot(res_enrich_Mars2vsCCT0_GT,aes(x=z_score, y =reorder(gs_description, +z_score)))+
  geom_bar(stat = "identity", aes(fill=z_score), colour = "black")+
  scale_fill_viridis(option = "D")+
  theme_classic()+
  theme(text=element_text(size =14))

#save
ggsave(plot = Mars2vsCCT0_BP, "FIGURES/IMMERSE paper/SCRIPT 7/Mars2_vsCCT0_BP.png", width = 10, height = 8, dpi = 300 )


##########################################################################################
##########################################################################################

#Mars3 vs cardiac T0

Mars3vsCCT0_GT <- results(dds, contrast = c("MARS", "Mars3", "Cardiac_T0"),lfcThreshold = 1, alpha = 0.05 )
Mars3vsCCT0_GT$SYMBOL <- rowData(dds)$SYMBOL

#make the information for state1vsCCT0_GT
de_symbols_Mars3vsCCT0_GT <- deseqresult2df(Mars3vsCCT0_GT, FDR = 0.05)$SYMBOL
bg_ids <- rowData(dds)$SYMBOL[rowSums(counts(dds)) > 0]


topgoDE_Mars3vsCCT0_GT <-
  pcaExplorer::topGOtable(de_symbols_Mars3vsCCT0_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_Mars3vsCCT0_GT <- shake_topGOtableResult(topgoDE_Mars3vsCCT0_GT)
colnames(topgoDE_Cohort_time_Cardiac_T1_vs_Cardiac_T0)

#create the annotation dataframe
anno_df <- data.frame(
  gene_id = rownames(dds),
  gene_name = mapIds(org.Hs.eg.db, keys = rownames(dds), column = "SYMBOL", keytype = "ENSEMBL"),
  stringsAsFactors = FALSE,
  row.names = rownames(dds)
)
head(anno_df)

#create the aggreegated scores
res_enrich_Mars3vsCCT0_GT <- get_aggrscores(res_enrich = res_enrich_Mars3vsCCT0_GT,
                                           res_de = Mars3vsCCT0_GT,
                                           annotation_obj = anno_df,
                                           aggrfun = mean)

#plot top 10 most significant pathways for Mars3vs cardiac T0
res_enrich_Mars3vsCCT0_GT <- res_enrich_Mars3vsCCT0_GT[1:10,] 

Mars3vsCCT0_BP <- ggplot(res_enrich_Mars3vsCCT0_GT,aes(x=z_score, y =reorder(gs_description, +z_score)))+
  geom_bar(stat = "identity", aes(fill=z_score), colour = "black")+
  scale_fill_viridis(option = "D")+
  theme_classic()+
  theme(text=element_text(size =14))

#save
ggsave(plot = Mars3vsCCT0_BP, "FIGURES/IMMERSE paper/SCRIPT 7/Mars3_extendedvsCCT0_BP.png", width = 10, height = 8, dpi = 300 )







