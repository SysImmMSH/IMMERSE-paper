
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
library("biomaRt")



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
table(samples$Cohort_time)

#rownames(samples) <- samples$name
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

#create deseq object
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ library_id + Cohort_time)


#ddsTxi$condition <- factor(ddsTxi$Cohort_time, levels = c("Cardiac_T0","Cardiac_T1",
                                             #     "Sepsis_T0", "Sepsis_T1",
                                              #    "Sepsis_T2", "Sepsis_T3" ))
ddsTxi$Cohort_time <- relevel(ddsTxi$Cohort_time, ref = "Cardiac_T0")

table(ddsTxi$sample_id, ddsTxi$name)
#this is a filtering step - removes any counts with less than 10 counts and are expressed in 5 samples
nrow(ddsTxi)
ddsTxi <- ddsTxi[rowSums(counts(ddsTxi) >= 10) >= 5]

nrow(ddsTxi)



#run deseqobject
dds <- DESeq(ddsTxi)


#####################################################################################

#look at how normalisation works
#add this value for adding a pseudo plus one to counts for plotting
epsilon <- 1
colnames(dds)

#plot the normalised data
# Open a pdf file
pdf("FIGURES/IMMERSE PAPER/SCRIPT 1/normalisation.pdf", width = 10, height = 50) 
par(mfrow=c(1,2),cex.lab=0.7)
boxplot(log2(counts(dds)+epsilon),  cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="Raw counts")
boxplot(log2(counts(dds, normalized=TRUE)+epsilon), cex.axis=0.7, 
        las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts") 
# Close the pdf file
dev.off() 




library(affy)
pdf("FIGURES/IMMERSE PAPER/SCRIPT 1/density.pdf", width = 12, height = 10) 
par(mfrow=c(1,2),cex.lab=0.7)
plotDensity(log2(counts(dds)+epsilon), 
            xlab="log2(counts)", cex.lab=0.7, panel.first=grid()) 
plotDensity(log2(counts(dds, normalized=TRUE)+epsilon), 
            xlab="log2(normalized counts)", cex.lab=0.7, panel.first=grid()) 
dev.off() 

###############################################################################################################
###############################################################################################################
# look at the comparisons  
design(dds)
colnames(dds)
resultsNames(dds)


#make results table for cardiacT0vscardiac_T1

Cohort_time_Cardiac_T1_vs_Cardiac_T0 <- results(dds, name = "Cohort_time_Cardiac_T1_vs_Cardiac_T0") %>% as.data.frame() %>% 
  rownames_to_column("GENEID")
Cohort_time_Cardiac_T1_vs_Cardiac_T0$gene_ID <- substr(Cohort_time_Cardiac_T1_vs_Cardiac_T0$GENEID,1,15) #need gene_id to allign with gene name

#######
#####
# THIS IS A GOOD WAY TO MATCH GENE ID's to NAMEs and also entrez number. 
library("biomaRt")
#the useEnsembl line worked as had to redirect to useast server. choose between useEnsembl and useMart option
#ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = 'useast')
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl", host = "useast.ensembl.org")
#mart = useMart("ensembl")
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = Cohort_time_Cardiac_T1_vs_Cardiac_T0$gene_ID,
                  mart = ensembl )
idx <- match( Cohort_time_Cardiac_T1_vs_Cardiac_T0$gene_ID, genemap$ensembl_gene_id )
Cohort_time_Cardiac_T1_vs_Cardiac_T0$entrez <- genemap$entrezgene[ idx ]
Cohort_time_Cardiac_T1_vs_Cardiac_T0$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(Cohort_time_Cardiac_T1_vs_Cardiac_T0)


write_csv(Cohort_time_Cardiac_T1_vs_Cardiac_T0, "Results_CSV/IMMERSE PAPER/SCRIPT 1/Cohort_time_Cardiac_T1_vs_Cardiac_T0.csv")

##############################################################################################################################
#make results table for Cohort_time_Sepsis_T0_vs_Cardiac_T0

res_Cohort_time_Sepsis_T0_vs_Cardiac_T0 <- results(dds, name = "Cohort_time_Sepsis_T0_vs_Cardiac_T0")

Cohort_time_Sepsis_T0_vs_Cardiac_T0 <- results(dds, name = "Cohort_time_Sepsis_T0_vs_Cardiac_T0")%>% as.data.frame() %>% 
  rownames_to_column("GENEID")
Cohort_time_Sepsis_T0_vs_Cardiac_T0$gene_ID <- substr(Cohort_time_Sepsis_T0_vs_Cardiac_T0$GENEID,1,15) #need gene_id to allign with gene name

genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = Cohort_time_Sepsis_T0_vs_Cardiac_T0$gene_ID,
                  mart = ensembl )
idx <- match( Cohort_time_Sepsis_T0_vs_Cardiac_T0$gene_ID, genemap$ensembl_gene_id )
Cohort_time_Sepsis_T0_vs_Cardiac_T0$entrez <- genemap$entrezgene[ idx ]
Cohort_time_Sepsis_T0_vs_Cardiac_T0$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(Cohort_time_Sepsis_T0_vs_Cardiac_T0)
write_csv(Cohort_time_Sepsis_T0_vs_Cardiac_T0, "Results_CSV/IMMERSE PAPER/SCRIPT 1/Cohort_time_Sepsis_T0_vs_Cardiac_T0.csv")

##############################################################################################################################

#make results table for Cohort_time_Sepsis_T1_vs_Cardiac_T0
Cohort_time_Sepsis_T1_vs_Cardiac_T0 <- results(dds, name = "Cohort_time_Sepsis_T1_vs_Cardiac_T0")%>% as.data.frame() %>% 
  rownames_to_column("GENEID")
Cohort_time_Sepsis_T1_vs_Cardiac_T0$gene_ID <- substr(Cohort_time_Sepsis_T1_vs_Cardiac_T0$GENEID,1,15) #need gene_id to allign with gene name


genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = Cohort_time_Sepsis_T1_vs_Cardiac_T0$gene_ID,
                  mart = ensembl )
idx <- match( Cohort_time_Sepsis_T1_vs_Cardiac_T0$gene_ID, genemap$ensembl_gene_id )
Cohort_time_Sepsis_T1_vs_Cardiac_T0$entrez <- genemap$entrezgene[ idx ]
Cohort_time_Sepsis_T1_vs_Cardiac_T0$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(Cohort_time_Sepsis_T1_vs_Cardiac_T0)
write_csv(Cohort_time_Sepsis_T1_vs_Cardiac_T0, "Results_CSV/IMMERSE PAPER/SCRIPT 1/Cohort_time_Sepsis_T1_vs_Cardiac_T0.csv")

##############################################################################################################################
#make results table for Cohort_time_Sepsis_T2_vs_Cardiac_T0
Cohort_time_Sepsis_T2_vs_Cardiac_T0 <- results(dds, name = "Cohort_time_Sepsis_T2_vs_Cardiac_T0")%>% as.data.frame() %>% 
  rownames_to_column("GENEID")

Cohort_time_Sepsis_T2_vs_Cardiac_T0$gene_ID <- substr(Cohort_time_Sepsis_T2_vs_Cardiac_T0$GENEID,1,15) #need gene_id to allign with gene name

genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = Cohort_time_Sepsis_T2_vs_Cardiac_T0$gene_ID,
                  mart = ensembl )
idx <- match( Cohort_time_Sepsis_T2_vs_Cardiac_T0$gene_ID, genemap$ensembl_gene_id )
Cohort_time_Sepsis_T2_vs_Cardiac_T0$entrez <- genemap$entrezgene[ idx ]
Cohort_time_Sepsis_T2_vs_Cardiac_T0$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(Cohort_time_Sepsis_T2_vs_Cardiac_T0)
write_csv(Cohort_time_Sepsis_T2_vs_Cardiac_T0, "Results_CSV/IMMERSE PAPER/SCRIPT 1/Cohort_time_Sepsis_T2_vs_Cardiac_T0.csv")

##############################################################################################################################
#make results table for Cohort_time_Sepsis_T3_vs_Cardiac_T0
Cohort_time_Sepsis_T3_vs_Cardiac_T0 <- results(dds, name = "Cohort_time_Sepsis_T3_vs_Cardiac_T0")%>% as.data.frame() %>% 
  rownames_to_column("GENEID")

Cohort_time_Sepsis_T3_vs_Cardiac_T0$gene_ID <- substr(Cohort_time_Sepsis_T3_vs_Cardiac_T0$GENEID,1,15) #need gene_id to allign with gene name


genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = Cohort_time_Sepsis_T3_vs_Cardiac_T0$gene_ID,
                  mart = ensembl )
idx <- match( Cohort_time_Sepsis_T3_vs_Cardiac_T0$gene_ID, genemap$ensembl_gene_id )
Cohort_time_Sepsis_T3_vs_Cardiac_T0$entrez <- genemap$entrezgene[ idx ]
Cohort_time_Sepsis_T3_vs_Cardiac_T0$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
head(Cohort_time_Sepsis_T3_vs_Cardiac_T0)
write_csv(Cohort_time_Sepsis_T3_vs_Cardiac_T0, "Results_CSV/IMMERSE PAPER/SCRIPT 1/Cohort_time_Sepsis_T3_vs_Cardiac_T0.csv")

##############################################################################################################################
# save the raw count data 
dds_df <- as.data.frame(counts(dds,normalized = FALSE)) %>% rownames_to_column("GENEID")
head(dds_df)
dds_df$gene_ID <- substr(dds_df$GENEID,1,15) #need gene_id to allign with gene name

genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = dds_df$gene_ID,
                  mart = ensembl )
idx <- match( dds_df$gene_ID, genemap$ensembl_gene_id )
dds_df$entrez <- genemap$entrezgene[ idx ]
dds_df$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
head(dds_df)
write_csv(dds_df, "Results_CSV/IMMERSE PAPER/SCRIPT 1/raw count data.csv")


##############################################################################################################################
# scale the data for PCA
vsd <- vst(dds, blind=FALSE)
assay(vsd)
#account for batch effect using limma
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$library_id)


#make data frame for vsd batch corrected count matrix por downstream analysis (making better PCA's/UMAPS)
vsd_df <- as.data.frame(assay(vsd)) %>% rownames_to_column("GENEID")
vsd_df$gene_ID <- substr(vsd_df$GENEID,1,15) #need gene_id to allign with gene name

genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = vsd_df$gene_ID,
                  mart = ensembl )
idx <- match( vsd_df$gene_ID, genemap$ensembl_gene_id )
vsd_df$entrez <- genemap$entrezgene[ idx ]
vsd_df$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
head(vsd_df)

write_csv(vsd_df, "Results_CSV/IMMERSE PAPER/SCRIPT 1/vsd_df.csv")
########################################################################################################################
#make a quick QC PCA to check  data
pcaData <-plotPCA(vsd, intgroup=c("Cohort_time", "library_id"), returnData=TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
batch_corrected_pca <- ggplot(pcaData, aes(PC1, PC2, color=Cohort_time)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
batch_corrected_pca

ggsave(plot = batch_corrected_pca, "FIGURES/IMMERSE PAPER/SCRIPT 1/batch_corrected.png", width = 6, height = 5, dpi = 300)



#####################################################################
#####################################################################
# Do analysis but have cardiac T1 as the comparison
#make new dds 

ddsTxi1 <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ library_id + Cohort_time)
#ddsTxi1 <- ddsTxi

ddsTxi1$Cohort_time <- relevel(ddsTxi1$Cohort_time, ref = "Cardiac_T1")

nrow(ddsTxi1)
ddsTxi1 <- ddsTxi1[rowSums(counts(ddsTxi1) >= 10) >= 5]

nrow(ddsTxi1)

#run deseqobject
dds1 <- DESeq(ddsTxi1)


resultsNames(dds1)

########################################################################################################################
#make results table for Cohort_time_Sepsis_T0_vs_Cardiac_T1
Cohort_time_Sepsis_T0_vs_Cardiac_T1 <- results(dds1, name = "Cohort_time_Sepsis_T0_vs_Cardiac_T1")%>% as.data.frame() %>% 
  rownames_to_column("GENEID")
Cohort_time_Sepsis_T0_vs_Cardiac_T1$gene_ID <- substr(Cohort_time_Sepsis_T0_vs_Cardiac_T1$GENEID,1,15) #need gene_id to allign with gene name

genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = Cohort_time_Sepsis_T0_vs_Cardiac_T1$gene_ID,
                  mart = ensembl )
idx <- match( Cohort_time_Sepsis_T0_vs_Cardiac_T1$gene_ID, genemap$ensembl_gene_id )
Cohort_time_Sepsis_T0_vs_Cardiac_T1$entrez <- genemap$entrezgene[ idx ]
Cohort_time_Sepsis_T0_vs_Cardiac_T1$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(Cohort_time_Sepsis_T0_vs_Cardiac_T1)
write_csv(Cohort_time_Sepsis_T0_vs_Cardiac_T1, "Results_CSV/IMMERSE PAPER/SCRIPT 1/Cohort_time_Sepsis_T0_vs_Cardiac_T1.csv")

########################################################################################################################
#make results table for Cohort_time_Sepsis_T1_vs_Cardiac_T1
Cohort_time_Sepsis_T1_vs_Cardiac_T1 <- results(dds1, name = "Cohort_time_Sepsis_T1_vs_Cardiac_T1")%>% as.data.frame() %>% 
  rownames_to_column("GENEID")
Cohort_time_Sepsis_T1_vs_Cardiac_T1$gene_ID <- substr(Cohort_time_Sepsis_T1_vs_Cardiac_T1$GENEID,1,15) #need gene_id to allign with gene name

genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = Cohort_time_Sepsis_T1_vs_Cardiac_T1$gene_ID,
                  mart = ensembl )
idx <- match( Cohort_time_Sepsis_T1_vs_Cardiac_T1$gene_ID, genemap$ensembl_gene_id )
Cohort_time_Sepsis_T1_vs_Cardiac_T1$entrez <- genemap$entrezgene[ idx ]
Cohort_time_Sepsis_T1_vs_Cardiac_T1$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(Cohort_time_Sepsis_T1_vs_Cardiac_T1)
write_csv(Cohort_time_Sepsis_T1_vs_Cardiac_T1, "Results_CSV/IMMERSE PAPER/SCRIPT 1/Cohort_time_Sepsis_T1_vs_Cardiac_T1.csv")


########################################################################################################################
#make results table for Cohort_time_Sepsis_T2_vs_Cardiac_T1
Cohort_time_Sepsis_T2_vs_Cardiac_T1 <- results(dds1, name = "Cohort_time_Sepsis_T2_vs_Cardiac_T1")%>% as.data.frame() %>% 
  rownames_to_column("GENEID")
Cohort_time_Sepsis_T2_vs_Cardiac_T1$gene_ID <- substr(Cohort_time_Sepsis_T2_vs_Cardiac_T1$GENEID,1,15) #need gene_id to allign with gene name

genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = Cohort_time_Sepsis_T2_vs_Cardiac_T1$gene_ID,
                  mart = ensembl )
idx <- match( Cohort_time_Sepsis_T2_vs_Cardiac_T1$gene_ID, genemap$ensembl_gene_id )
Cohort_time_Sepsis_T2_vs_Cardiac_T1$entrez <- genemap$entrezgene[ idx ]
Cohort_time_Sepsis_T2_vs_Cardiac_T1$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(Cohort_time_Sepsis_T2_vs_Cardiac_T1)
write_csv(Cohort_time_Sepsis_T2_vs_Cardiac_T1, "Results_CSV/IMMERSE PAPER/SCRIPT 1/Cohort_time_Sepsis_T2_vs_Cardiac_T1.csv")

########################################################################################################################
#make results table for Cohort_time_Sepsis_T3_vs_Cardiac_T1
Cohort_time_Sepsis_T3_vs_Cardiac_T1 <- results(dds1, name = "Cohort_time_Sepsis_T3_vs_Cardiac_T1")%>% as.data.frame() %>% 
  rownames_to_column("GENEID")

Cohort_time_Sepsis_T3_vs_Cardiac_T1$gene_ID <- substr(Cohort_time_Sepsis_T3_vs_Cardiac_T1$GENEID,1,15) #need gene_id to allign with gene name


genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = Cohort_time_Sepsis_T3_vs_Cardiac_T1$gene_ID,
                  mart = ensembl )
idx <- match( Cohort_time_Sepsis_T3_vs_Cardiac_T1$gene_ID, genemap$ensembl_gene_id )
Cohort_time_Sepsis_T3_vs_Cardiac_T1$entrez <- genemap$entrezgene[ idx ]
Cohort_time_Sepsis_T3_vs_Cardiac_T1$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(Cohort_time_Sepsis_T3_vs_Cardiac_T1)
write_csv(Cohort_time_Sepsis_T3_vs_Cardiac_T1, "Results_CSV/IMMERSE PAPER/SCRIPT 1/Cohort_time_Sepsis_T3_vs_Cardiac_T1.csv")


########################################################################################################################
#lets look at srs signatures

library(SepstratifieR)
library(data.table)
library(edgeR)


# i only run  this to get the variables for the normalisation downstream after running the cpm.
#it is important to do it in this order. Get counts, calcuate CPM, batch correction.
ntd <- normTransform(dds)
# extract the raw counts. Important it is raw counts here
counts <- counts(dds, normalized = FALSE)

DGEcpm <- as.data.frame(cpm(counts, log = TRUE, prior.count = TRUE), keep.rownames = TRUE)  %>% 
  as.matrix()

DGEcpm <- limma::removeBatchEffect(DGEcpm, ntd$library_id)%>% as.data.frame() %>% rownames_to_column("gene_id")


DGEcpm$gene_id <- substr(DGEcpm$gene_id,1,15) #need gene_id to align with gene name
head(counts)
head(DGEcpm)

# Stratify patients
#select genes for extended srs phenotyping
geneid <- data.frame (gene_id  = c("ENSG00000144659",
                                   "ENSG00000103423",
                                   "ENSG00000135372",
                                   "ENSG00000079134",
                                   "ENSG00000135972",
                                   "ENSG00000087157",
                                   "ENSG00000165006",
                                   "ENSG00000111667",
                                   "ENSG00000182670",
                                   "ENSG00000097033",
                                   "ENSG00000165733",
                                   "ENSG00000103264",
                                   "ENSG00000152219",
                                   "ENSG00000100814",
                                   "ENSG00000127334",
                                   "ENSG00000131355",
                                   "ENSG00000137337",
                                   "ENSG00000156414",
                                   "ENSG00000115085"),
                  gene_name = c("SLC25A38",
                                "DNAJA3",
                                "NAT10",
                                "THOC1",
                                "MRPS9",
                                "PGS1",
                                "UBAP1",
                                "USP5",
                                "TTC3",
                                "SH3GLB1",
                                "BMS1",
                                "FBXO31",
                                "ARL14EP",
                                "CCNB1IP1",
                                "DYRK2",
                                "ADGRE3",
                                "MDC1",
                                "TDRD9",
                                "ZAP70")
)

colnames(DGEcpm)
dds_df_norm_srs <- DGEcpm %>% filter(gene_id  %in% geneid$gene_id) %>% dplyr::select(gene_id, 2:259)


dds_df_norm_srs = setNames(data.frame(t(dds_df_norm_srs[,-1])), dds_df_norm_srs[,1])%>% na.omit()

class(dds_df_norm_srs)
predictions <- stratifyPatients(dds_df_norm_srs, gene_set = "extended")
davenport_SRS_predictions <- stratifyPatients(dds_df_norm_srs, k = 144, gene_set = "davenport")
plotAlignedSamples(predictions)
plotAlignedSamples(davenport_SRS_predictions)


pred <- data.frame(group = predictions@SRS,
                   SRSq = predictions@SRSq) %>% rownames_to_column("sample_id")

write_csv(pred, "Results_CSV/IMMERSE PAPER/SCRIPT 1/srs_groupings_extended.csv")

pred1 <- data.frame(group = davenport_SRS_predictions@SRS,
                   SRSq = davenport_SRS_predictions@SRSq) %>% rownames_to_column("sample_id")# %>% filter(!str_detect(sample_id, 'SC100'))

write_csv(pred1, "Results_CSV/IMMERSE PAPER/SCRIPT 1/srs_groupings_davenport.csv")

pred$cohort <- substr(pred$sample_id, 1,2)
pred$timepoint <- substr(pred$sample_id, 6,7)
pred <- pred %>% unite(Cohort_time, c(cohort, timepoint), sep = "_", remove = FALSE) %>% as.data.frame()


ggplot(pred, aes(x=group, y = SRSq))+
  geom_boxplot()+geom_jitter()

ggplot(pred, aes(x=Cohort_time, y = SRSq))+
  geom_boxplot()+geom_jitter()


pred1$cohort <- substr(pred1$sample_id, 1,2)
pred1$timepoint <- substr(pred1$sample_id, 6,7)
pred1$patient_id <- substr(pred1$sample_id, 1,4)

pred1 <- pred1 %>% unite(Cohort_time, c(cohort, timepoint), sep = "_", remove = FALSE) %>% as.data.frame()


ggplot(pred1, aes(x=group, y = SRSq))+
  geom_boxplot()+geom_jitter()

ggplot(pred1, aes(x=Cohort_time, y = SRSq))+
  geom_boxplot()+geom_jitter()



#######################################################################################################################
# try looking at mars endotypes
### Brendon Scicluna - CEMM AMC
### endotype predictions
### 140 gene classifier in both training and validation/replication sets

library(CMA)
library(e1071)
library(randomForest)

### training set - AMC cohort 140 gene classifier x 306 patients

t <- read_csv("MARS/discovery_140genes_matrix_positiveSil_rename_removed.csv")

g.amc<-read.csv(file="MARS/discovery_140genes_matrix_positiveSil_rename_removed.csv",header =T) %>% column_to_rownames("gene")

g.m.1=data.matrix(g.amc[,-1])
head(g.m.1)
phenos<-read.csv(file="MARS/Positive_silhouette_discoveryAMC_endo_markup.csv",header=T)# phenotype table with endotype markup

#### validation set - UMCU cohort 140 gene classifier x 216 patients

genemap <- getBM( attributes = c( "ensembl_gene_id","hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = DGEcpm$gene_id,
                  mart = ensembl )
idx <- match( DGEcpm$gene_id, genemap$ensembl_gene_id )
DGEcpm$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(DGEcpm)
dds_df_norm_MARS1 <- DGEcpm %>% filter(hgnc_symbol  %in% t$gene) %>% dplyr::select(hgnc_symbol, 2:260) %>% column_to_rownames("hgnc_symbol")
head(t)
check <- t %>% filter(!gene  %in% DGEcpm$hgnc_symbol) 
#If all the genes have alligned (i.e., no updates in gene names, this tibble should be empty)



#g.umc<-read.table(,header=T)
g.m.2=data.matrix(dds_df_norm_MARS1[,-1])


d.1 = sweep(g.m.1,1, apply(g.m.1,1,median,na.rm=T))
d.2 = sweep(g.m.2,1, apply(g.m.2,1,median,na.rm=T))

##### Endotype PREDICTION by random forest method

pred_mars <- prediction(X.tr=t(d.1),y.tr=as.factor(phenos$endotype),X.new=t(d.2), classifier = rfCMA)

table(show(pred_mars))
pred_mars@yhat

#make a data frame with the variables
mars <- cbind(rownames(pred_mars@Xnew), as.data.frame(pred_mars@yhat)) %>% rename("rownames(pred_mars@Xnew)"="sample_id",
                                                                                  "pred_mars@yhat" ="MARS" )
head(mars)

write_csv(mars, "Results_CSV/IMMERSE PAPER/SCRIPT 1/mars_groupings.csv")

##########################################################################################################################################
# Sweeney classifier model 

#my data 

vsd_df_sween <- vsd_df %>% dplyr::select(hgnc_symbol,2:259)
colnames(vsd_df_sween)
  


#gene list
geneid <- c("ARG1",  "LCN2",  "LTF",   "OLFM4","HLA-DMB","YKT6", "PDE4B","POLR1F",  "BTN2A2",   "ZBTB33",
            "PSMB9",    "CAMK4",    "TMEM19",   "SLC12A7",  "TP53BP1", "PLEKHO1",  "SLC25A22", "FRS2",
            "KCNMB4", "CRISP2", "HTRA1",  "PPL","RHBDF2",  "ZCCHC4",    "DDX6",    "SENP5",
            "RAPGEF1", "DTX2",    "RELB","GADD45A", "CD24",    "S100A12", "STX1A" )
out2

### load endotyping file originally derived for Sweeney et al CCM 2018 paper: forwardSearch_SAM_multLR_out.RData
load("/Volumes/Samsung_X5/work/PhD all/RNAseq analysis/RNAseq V1 Current/Sweeney classifier/endotypesShare.Rdata")

out2[[1]] <- c("YKT6",     "PDE4B",   "POLR1F",  "BTN2A2",   "ZBTB33",   "PSMB9",    "CAMK4",
               "TMEM19",   "SLC12A7",  "TP53BP1", "PLEKHO1",  "SLC25A22", "FRS2" )
### how to calculate diff-of-geo-mean scores
calcScore <- function(mtx, pos, neg){
  geomMean <- function (x) return(exp(sum(log(x))/length(x)))
  scores <- unlist(apply(mtx, 2, function(x) geomMean(x[pos]) - length(neg)/length(pos) * geomMean(x[neg]) ))
  scale(scores)
}

#need to get the same genes as sweeney
GenesShare <- genesShare %>% as.data.frame() %>% rownames_to_column("Gene_names") 

colnames(vsd_df_sween)
vsd_df_sween2 <- vsd_df_sween %>% filter(hgnc_symbol  %in% geneid)%>% 
  column_to_rownames("hgnc_symbol")%>% as.matrix()
vsd_df_sween2$hgnc_symbol
out2
endotypeScores <- data.frame(score1 = calcScore(vsd_df_sween2, out1[[1]], out1[[2]]),
                             score2 = calcScore(vsd_df_sween2, out2[[1]], out2[[2]]),
                             score3 = calcScore(vsd_df_sween2, out3[[1]], out3[[2]]) )
### calculate endotype probabilities
### 1 = inflammopathic; 2 = adaptive; 3 = coagulopathic
### 'endotype' is greatest-prob discrete class
### endoProbs are the continuous probabilties for each class
library(nnet)
probs <- data.frame(endoProbs=round(predict(multLR, endotypeScores, type="probs"), 3), 
                    endotype=predict(multLR, endotypeScores, type="class"))
head(probs)
#add in the info and save
probs <- probs %>% mutate(`Sweeney 2019` = case_when(endotype == "1" ~ "inflammopathic",
                                              endotype == "2" ~ "adaptive",
                                              endotype == "3" ~ "coagulopathic")) %>% 
  rownames_to_column("sample_id") %>% 
  dplyr::select(sample_id, `Sweeney 2019`)
write_csv(probs, "Results_CSV/IMMERSE paper/SCRIPT 1/sweeney_groups.csv")


#####################################################################
#####################################################################
# Do analysis but have Sepsis T1 as the comparison. this is  to then compare the ipa analysis from baseline through to discharge
#make new dds 

#ddsTxi <- DESeqDataSetFromTximport(txi,
#                                  colData = samples,
#                                 design = ~ library_id + Cohort_time)
ddsTxi2 <- ddsTxi

ddsTxi2$Cohort_time <- relevel(ddsTxi2$Cohort_time, ref = "Sepsis_T0")


#run deseqobject
dds2 <- DESeq(ddsTxi2)

#look at how normalisation works
#add this value for adding a pseudo plus one to counts for plotting
#epsilon <- 1
#colnames(dds)

resultsNames(dds2)


########################################################################################################################



#make results table for Cohort_time_Sepsis_T1_vs_Sepsis_T0
Cohort_time_Sepsis_T1_vs_Sepsis_T0 <- results(dds2, name = "Cohort_time_Sepsis_T1_vs_Sepsis_T0")%>% as.data.frame() %>% 
  rownames_to_column("GENEID")

Cohort_time_Sepsis_T1_vs_Sepsis_T0$gene_ID <- substr(Cohort_time_Sepsis_T1_vs_Sepsis_T0$GENEID,1,15) #need gene_id to allign with gene name


genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = Cohort_time_Sepsis_T1_vs_Sepsis_T0$gene_ID,
                  mart = ensembl )
idx <- match( Cohort_time_Sepsis_T1_vs_Sepsis_T0$gene_ID, genemap$ensembl_gene_id )
Cohort_time_Sepsis_T1_vs_Sepsis_T0$entrez <- genemap$entrezgene[ idx ]
Cohort_time_Sepsis_T1_vs_Sepsis_T0$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(Cohort_time_Sepsis_T1_vs_Sepsis_T0)
write_csv(Cohort_time_Sepsis_T1_vs_Sepsis_T0, "Results_CSV/IMMERSE paper/SCRIPT 1/Cohort_time_Sepsis_T1_vs_Sepsis_T0.csv")


########################################################################################################################


#make results table for Cohort_time_Sepsis_T2_vs_Sepsis_T0
Cohort_time_Sepsis_T2_vs_Sepsis_T0 <- results(dds2, name = "Cohort_time_Sepsis_T2_vs_Sepsis_T0")%>% as.data.frame() %>% 
  rownames_to_column("GENEID")

Cohort_time_Sepsis_T2_vs_Sepsis_T0$gene_ID <- substr(Cohort_time_Sepsis_T2_vs_Sepsis_T0$GENEID,1,15) #need gene_id to allign with gene name


genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = Cohort_time_Sepsis_T2_vs_Sepsis_T0$gene_ID,
                  mart = ensembl )
idx <- match( Cohort_time_Sepsis_T2_vs_Sepsis_T0$gene_ID, genemap$ensembl_gene_id )
Cohort_time_Sepsis_T2_vs_Sepsis_T0$entrez <- genemap$entrezgene[ idx ]
Cohort_time_Sepsis_T2_vs_Sepsis_T0$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(Cohort_time_Sepsis_T2_vs_Sepsis_T0)
write_csv(Cohort_time_Sepsis_T2_vs_Sepsis_T0, "Results_CSV/Cohort_time_Sepsis_T2_vs_Sepsis_T0.csv")


########################################################################################################################


########################################################################################################################


#make results table for Cohort_time_Sepsis_T3_vs_Sepsis_T0
Cohort_time_Sepsis_T3_vs_Sepsis_T0 <- results(dds2, name = "Cohort_time_Sepsis_T3_vs_Sepsis_T0")%>% as.data.frame() %>% 
  rownames_to_column("GENEID")

Cohort_time_Sepsis_T3_vs_Sepsis_T0$gene_ID <- substr(Cohort_time_Sepsis_T3_vs_Sepsis_T0$GENEID,1,15) #need gene_id to allign with gene name


genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = Cohort_time_Sepsis_T3_vs_Sepsis_T0$gene_ID,
                  mart = ensembl )
idx <- match( Cohort_time_Sepsis_T3_vs_Sepsis_T0$gene_ID, genemap$ensembl_gene_id )
Cohort_time_Sepsis_T3_vs_Sepsis_T0$entrez <- genemap$entrezgene[ idx ]
Cohort_time_Sepsis_T3_vs_Sepsis_T0$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(Cohort_time_Sepsis_T3_vs_Sepsis_T0)
write_csv(Cohort_time_Sepsis_T3_vs_Sepsis_T0, "Results_CSV/Cohort_time_Sepsis_T3_vs_Sepsis_T0.csv")


########################################################################################################################
########################################################################################################################
########################################################################################################################

# Do a time analysis using a likelihood ratio test (LRT) This is from the deseq vignette at the bottom. 
# https://support.bioconductor.org/p/113630/ question on how to run it with one group
#subselect just sepsis patients
dds_sepsis <- dds[,dds$Cohort == "Sepsis"]

# create new desq object with new design matrix
dds_sepsis <- DESeqDataSet(dds_sepsis, ~ library_id + Timepoint)
# run deseq with the LRT model
dds_sepsis <- DESeq(dds_sepsis, test="LRT", reduced = ~ library_id)

# get results 
resSepsis <- results(dds_sepsis)%>% as.data.frame() %>% 
  rownames_to_column("GENEID")
resSepsis$gene_ID <- substr(resSepsis$GENEID,1,15)
x <- listAttributes(ensembl)
x$name
table(x$page)
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol","description","gene_biotype"),
                  filters = "ensembl_gene_id",
                  values = resSepsis$gene_ID,
                  mart = ensembl )
idx <- match( resSepsis$gene_ID, genemap$ensembl_gene_id)
resSepsis$entrez <- genemap$entrezgene[ idx ]
resSepsis$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
resSepsis$description <- genemap$description[ idx ]
resSepsis$gene_biotype <- genemap$gene_biotype[ idx ]

# save the dataframe in results 
write_csv(resSepsis, "Results_csv/IMMERSE paper/SCRIPT 1/LRT analysis sepsis timecourse.csv")

########################################################################################################################
########################################################################################################################
########################################################################################################################

#create a dataframe and save it with gene expression information for the 718 genes. This is for the
#integrated analysis later.
resSepsis1 <- resSepsis %>% filter(padj <=0.01)%>% 
  filter(log2FoldChange > 1 |log2FoldChange < -1) %>% 
  mutate(gene_id_2 = coalesce(hgnc_symbol,gene_ID))

resSepsis1 <- resSepsis1$gene_ID

vsd_df_718 <- vsd_df %>% filter(gene_ID %in% resSepsis1)


#save csv
write_csv(vsd_df_718, "Results_csv/IMMERSE paper/SCRIPT 1/gene_data_718_genes.csv")





