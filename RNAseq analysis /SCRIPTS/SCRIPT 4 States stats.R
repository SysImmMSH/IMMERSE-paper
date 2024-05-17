# need to do stats on states 
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
library("biomaRt")
library(broom)
library(heatmaply)

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
states_meta <- read_csv("Metadata/metadata_states_Phate.csv") %>% dplyr::select(sample_id, State,Lineage1,PHATE1, PHATE2)

samples <- merge(samples, states_meta, all.x = TRUE)

samples <- samples %>%  dplyr::mutate(State = coalesce(State, Cohort_time)) 

####### directory to salmon quant files
#directory of quant.sf files
dir <- "quantfiles"
#list them/check they are there
list.files(file.path(dir))

files <- file.path(dir, samples$name)

names(files) <- samples$sample_id
all(file.exists(files))

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
                                   design = ~ library_id + State)

nrow(ddsTxi) #60230 genes measured
ddsTxi <- ddsTxi[rowSums(counts(ddsTxi) >= 10) >= 20]
nrow(ddsTxi)#25637 left after filtering



dds <- DESeq(ddsTxi)



#########################################################################################
#########################################################################################
#########################################################################################

state1vsCCT0 <- results(dds, contrast = c("State", "State 1", "Cardiac_T0") ) %>% as.data.frame() %>% rownames_to_column("GENEID")
state1vsCCT0$gene_ID <- substr(state1vsCCT0$GENEID,1,15) #need gene_id to allign with gene name


ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl", host = "useast.ensembl.org")
#mart = useMart("ensembl")
genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                  filters = "ensembl_gene_id",
                  values = state1vsCCT0$gene_ID,
                  mart = ensembl )
idx <- match( state1vsCCT0$gene_ID, genemap$ensembl_gene_id )
state1vsCCT0$gene_biotype <- genemap$gene_biotype[ idx ]
state1vsCCT0$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(state1vsCCT0)

state1vsCCT0 <- state1vsCCT0%>%  mutate(gene_name = na_if(hgnc_symbol, ""))   %>%  mutate(gene_name = coalesce(hgnc_symbol,gene_ID )) 

write_csv(state1vsCCT0, "Results_CSV/IMMERSE paper/SCRIPT 4/state1vsCCT0.csv")

#########################################################################################

state2vsCCT0 <- results(dds, contrast = c("State", "State 2", "Cardiac_T0") ) %>% as.data.frame() %>% rownames_to_column("GENEID")
state2vsCCT0$gene_ID <- substr(state2vsCCT0$GENEID,1,15) #need gene_id to allign with gene name

#mart = useMart("ensembl")
genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                  filters = "ensembl_gene_id",
                  values = state2vsCCT0$gene_ID,
                  mart = ensembl )
idx <- match( state2vsCCT0$gene_ID, genemap$ensembl_gene_id )
state2vsCCT0$gene_biotype <- genemap$gene_biotype[ idx ]
state2vsCCT0$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(state2vsCCT0)

state2vsCCT0 <- state2vsCCT0%>%  mutate(gene_name = na_if(hgnc_symbol, ""))   %>%  mutate(gene_name = coalesce(hgnc_symbol, gene_ID)) 


write_csv(state2vsCCT0, "Results_CSV/IMMERSE paper/SCRIPT 4/state2vsCCT0.csv")

#########################################################################################

state3vsCCT0 <- results(dds, contrast = c("State", "State 3", "Cardiac_T0") ) %>% as.data.frame() %>% rownames_to_column("GENEID")
state3vsCCT0$gene_ID <- substr(state3vsCCT0$GENEID,1,15) #need gene_id to allign with gene name


#mart = useMart("ensembl")
genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                  filters = "ensembl_gene_id",
                  values = state3vsCCT0$gene_ID,
                  mart = ensembl )
idx <- match( state3vsCCT0$gene_ID, genemap$ensembl_gene_id )
state3vsCCT0$gene_biotype <- genemap$gene_biotype[ idx ]
state3vsCCT0$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(state3vsCCT0)
state3vsCCT0 <- state3vsCCT0%>%  mutate(gene_name = na_if(hgnc_symbol, ""))  %>%  mutate(gene_name = coalesce(hgnc_symbol, gene_ID)) 


write_csv(state3vsCCT0, "Results_CSV/IMMERSE paper/SCRIPT 4/state3vsCCT0.csv")

#########################################################################################

state1vsCCT1 <- results(dds, contrast = c("State", "State 1", "Cardiac_T1") ) %>% as.data.frame() %>% rownames_to_column("GENEID")
state1vsCCT1$gene_ID <- substr(state1vsCCT1$GENEID,1,15) #need gene_id to allign with gene name


#mart = useMart("ensembl")
genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                  filters = "ensembl_gene_id",
                  values = state1vsCCT1$gene_ID,
                  mart = ensembl )
idx <- match( state1vsCCT1$gene_ID, genemap$ensembl_gene_id )
state1vsCCT1$gene_biotype <- genemap$gene_biotype[ idx ]
state1vsCCT1$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(state1vsCCT1)
state1vsCCT1 <- state1vsCCT1%>%  mutate(gene_name = na_if(hgnc_symbol, ""))  %>%  mutate(gene_name = coalesce(hgnc_symbol, gene_ID)) 


write_csv(state1vsCCT1, "Results_CSV/IMMERSE paper/SCRIPT 4/state1vsCCT1.csv")

#########################################################################################

state2vsCCT1 <- results(dds, contrast = c("State", "State 2", "Cardiac_T1") ) %>% as.data.frame() %>% rownames_to_column("GENEID")
state2vsCCT1$gene_ID <- substr(state2vsCCT1$GENEID,1,15) #need gene_id to allign with gene name


#mart = useMart("ensembl")
genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                  filters = "ensembl_gene_id",
                  values = state2vsCCT1$gene_ID,
                  mart = ensembl )
idx <- match( state2vsCCT1$gene_ID, genemap$ensembl_gene_id )
state2vsCCT1$gene_biotype <- genemap$gene_biotype[ idx ]
state2vsCCT1$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(state2vsCCT1)
state2vsCCT1 <- state2vsCCT1%>%  mutate(gene_name = na_if(hgnc_symbol, ""))  %>%  mutate(gene_name = coalesce(hgnc_symbol, gene_ID)) 


write_csv(state2vsCCT1, "Results_CSV/IMMERSE paper/SCRIPT 4/state2vsCCT1.csv")

#########################################################################################

state3vsCCT1 <- results(dds, contrast = c("State", "State 3", "Cardiac_T1") ) %>% as.data.frame() %>% rownames_to_column("GENEID")
state3vsCCT1$gene_ID <- substr(state3vsCCT1$GENEID,1,15) #need gene_id to allign with gene name


#mart = useMart("ensembl")
genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                  filters = "ensembl_gene_id",
                  values = state3vsCCT1$gene_ID,
                  mart = ensembl )
idx <- match( state3vsCCT1$gene_ID, genemap$ensembl_gene_id )
state3vsCCT1$gene_biotype <- genemap$gene_biotype[ idx ]
state3vsCCT1$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(state3vsCCT1)
state3vsCCT1 <- state3vsCCT1%>%  mutate(gene_name = na_if(hgnc_symbol, ""))  %>%  mutate(gene_name = coalesce(hgnc_symbol, gene_ID)) 


write_csv(state3vsCCT1, "Results_CSV/IMMERSE paper/SCRIPT 4/state3vsCCT1.csv")

#########################################################################################

state2vsstate1 <- results(dds, contrast = c("State", "State 2", "State 1") ) %>% as.data.frame() %>% rownames_to_column("GENEID")
state2vsstate1$gene_ID <- substr(state2vsstate1$GENEID,1,15) #need gene_id to align with gene name


#mart = useMart("ensembl")
genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                  filters = "ensembl_gene_id",
                  values = state2vsstate1$gene_ID,
                  mart = ensembl )
idx <- match( state2vsstate1$gene_ID, genemap$ensembl_gene_id )
state2vsstate1$gene_biotype <- genemap$gene_biotype[ idx ]
state2vsstate1$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(state3vsCCT1)
state2vsstate1 <- state2vsstate1 %>%  mutate(gene_name = na_if(hgnc_symbol, "")) %>%  mutate(gene_name = coalesce(hgnc_symbol, gene_ID)) 


write_csv(state2vsstate1, "Results_CSV/IMMERSE paper/SCRIPT 4/state2vsstate1.csv")

#########################################################################################

state3vsstate1 <- results(dds, contrast = c("State", "State 3", "State 1") ) %>% as.data.frame() %>% rownames_to_column("GENEID")
state3vsstate1$gene_ID <- substr(state3vsstate1$GENEID,1,15) #need gene_id to align with gene name


#mart = useMart("ensembl")
genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                  filters = "ensembl_gene_id",
                  values = state3vsstate1$gene_ID,
                  mart = ensembl )
idx <- match( state3vsstate1$gene_ID, genemap$ensembl_gene_id )
state3vsstate1$gene_biotype <- genemap$gene_biotype[ idx ]
state3vsstate1$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

head(state3vsstate1)
state3vsstate1 <- state3vsstate1 %>%  mutate(gene_name = na_if(hgnc_symbol, "")) %>%  mutate(gene_name = coalesce(hgnc_symbol, gene_ID)) 


write_csv(state3vsstate1, "Results_CSV/IMMERSE paper/SCRIPT 4/state3vsstate1.csv")



##############################################################################################################
##############################################################################################################
##############################################################################################################

########################################################################################################################################
########################################################################################################################################

# Want to plot the differntial gene expression of states vs pre-surgery and states 2 & 3 vs state 1

########################################################################################################################################

# State 1 vs pre surgery
#state1vsCCT0 <- read_csv("Results_csv/States/Successfully Uploaded 2022-09-08 14-55-46/state1vsCCT0.csv")



state1vsCCT0 <- state1vsCCT0 %>% dplyr::select(log2FoldChange, gene_ID, hgnc_symbol, padj) %>% 
  dplyr::rename("State 1 vs pre_log2FoldChange"="log2FoldChange",
                "State 1 vs pre_padj" = "padj")


########################################################################################################################################

# State 2 vs pre surgery
#state2vsCCT0 <- read_csv("Results_csv/States/Successfully Uploaded 2022-09-08 14-55-46/state2vsCCT0.csv")



state2vsCCT0 <- state2vsCCT0 %>% dplyr::select(log2FoldChange, gene_ID, hgnc_symbol, padj) %>% 
  dplyr::rename("State 2 vs pre_log2FoldChange"="log2FoldChange",
                "State 2 vs pre_padj" = "padj")



########################################################################################################################################

# State 3 vs pre surgery
#state3vsCCT0 <- read_csv("Results_csv/States/Successfully Uploaded 2022-09-08 14-55-46/state3vsCCT0.csv")



state3vsCCT0 <- state3vsCCT0 %>% dplyr::select(log2FoldChange, gene_ID, hgnc_symbol, padj) %>% 
  dplyr::rename("State 3 vs pre_log2FoldChange"="log2FoldChange",
                "State 3 vs pre_padj" = "padj")

########################################################################################################################################

# State 1 vs post surgery
#state1vsCCT1 <- read_csv("Results_csv/States/Successfully Uploaded 2022-09-08 14-55-46/state1vsCCT1.csv")



state1vsCCT1 <- state1vsCCT1 %>% dplyr::select(log2FoldChange, gene_ID, hgnc_symbol, padj) %>% 
  dplyr::rename("State 1 vs post_log2FoldChange"="log2FoldChange",
                "State 1 vs post_padj" = "padj")

########################################################################################################################################

# State 2 vs post surgery
#state2vsCCT1 <- read_csv("Results_csv/States/Successfully Uploaded 2022-09-08 14-55-46/state2vsCCT1.csv")



state2vsCCT1 <- state2vsCCT1 %>% dplyr::select(log2FoldChange, gene_ID, hgnc_symbol, padj) %>% 
  dplyr::rename("State 2 vs post_log2FoldChange"="log2FoldChange",
                "State 2 vs post_padj" = "padj")

########################################################################################################################################

# State 3 vs post surgery
#state3vsCCT1 <- read_csv("Results_csv/States/Successfully Uploaded 2022-09-08 14-55-46/state3vsCCT1.csv")



state3vsCCT1 <- state3vsCCT1 %>% dplyr::select(log2FoldChange, gene_ID, hgnc_symbol, padj) %>% 
  dplyr::rename("State 3 vs post_log2FoldChange"="log2FoldChange",
                "State 3 vs post_padj" = "padj")


########################################################################################################################################

# State 2 vs state 1
#state2vsstate1 <-read_csv("Results_csv/States/state2vsstate1.csv") 



state2vsstate1 <- state2vsstate1 %>% dplyr::select(log2FoldChange, gene_ID, hgnc_symbol, padj) %>% 
  dplyr::rename("State 2 vs 1_log2FoldChange"="log2FoldChange",
                "State 2 vs 1_padj" = "padj") %>% 
  dplyr::select(-hgnc_symbol)

########################################################################################################################################

# State 3 vs state 1
#state3vsstate1 <-read_csv("Results_csv/States/state3vsstate1.csv") 



state3vsstate1 <- state3vsstate1 %>% dplyr::select(log2FoldChange, gene_ID, hgnc_symbol, padj) %>% 
  dplyr::rename("State 3 vs 1_log2FoldChange"="log2FoldChange",
                "State 3 vs 1_padj" = "padj") %>% 
  dplyr::select(-hgnc_symbol)


identical(state2vsstate1[['hgnc_symbol']],state2vsCCT1[['hgnc_symbol']])
########################################################################################################################################
########################################################################################################################################
# Combine and pivot
df <- merge(state1vsCCT0, state2vsCCT0, all = T)
df <- merge(df,state3vsCCT0, all = T)
df <- merge(df,state1vsCCT1, all = T)
df <- merge(df,state2vsCCT1, all = T)
df <- merge(df,state3vsCCT1, all = T)
df <- merge(df,state2vsstate1, all = T)
df <- merge(df,state3vsstate1, all = T)

colnames(df)

df1 <- df %>%  pivot_longer(cols = c(3:18),
                            names_to = c("Comparison", "statistic"),
                            names_sep = "_",
                            values_to = "value") %>% 
  spread(statistic,value) %>% 
  filter(`padj` <= 0.01) %>% 
  filter(log2FoldChange > 1 |log2FoldChange < -1) %>% 
  filter(log2FoldChange < 15 & !log2FoldChange < -15)

head(df1)
# add a column of NAs
df1$diffexpressed <- "NA"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df1$diffexpressed[df1$log2FoldChange > 1 & df1$padj < 0.01] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df1$diffexpressed[df1$log2FoldChange < -1 & df1$padj < 0.01] <- "DOWN"




df1$diffexpressed1 <- "NA"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df1$diffexpressed1[df1$log2FoldChange > 3 & df1$padj < 0.01] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df1$diffexpressed1[df1$log2FoldChange < -3 & df1$padj < 0.01] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
df1$delabel <- NA
df1$delabel[df1$diffexpressed1 != "NA"] <- df1$hgnc_symbol[df1$diffexpressed1 != "NA"]



########################################################################################################################################
########################################################################################################################################

# Plot the graph to as a dot plot in columns 
df1$Comparison
df1$Comparison <- factor(df1$Comparison, levels = c("State 1 vs pre", "State 2 vs pre", "State 3 vs pre",
                                                    "State 1 vs post","State 2 vs post","State 3 vs post",
                                                    "State 2 vs 1", "State 3 vs 1"))
table(df1$Comparison)
mycolors <- c("#47abd8", "#D01B1B")
names(mycolors) <- c("DOWN", "UP")
library(ggbeeswarm)

df1 <- df1 %>% mutate(
  `-log10 padj` = -log10(padj)
)
df1$delabel


plot <- ggplot(df1, aes(x= Comparison, y = log2FoldChange, label = delabel))+ 
  geom_quasirandom(method = "tukey", aes(fill =diffexpressed, size =`-log10 padj`), pch = 21, colour = "black", alpha = 0.5)+
  scale_fill_manual(values = mycolors)+
  scale_size_continuous(range = c(1, 5),
             breaks =  c(5, 10, 50, 100)) +
  #geom_text_repel(max.overlaps =  20, position=position_quasirandom(method = "tukey"))+
  theme(panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white"))+
  geom_hline(yintercept = 1, color = "black")+
  geom_hline(yintercept = -1, color = "black")+
  ylim(c(-16,16))+
  xlab("Comparison")+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))+ 
 # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_label_repel(max.overlaps = 20, force = 10)
plot


ggsave(plot = plot, "FIGURES/IMMERSE paper/SCRIPT 4/States_gene_expression_barplot_1.png",width=25, height=12,units="in",dpi = 300)


#want to show the number of genes changing in this comparison 

count <- table(df1$Comparison, df1$diffexpressed) %>% as.data.frame() 

count1 <- count %>% filter(Var1 ==  "State 1 vs pre")

p1 <- ggplot(count1, aes(y=Freq, x = Var2))+
  geom_bar(stat= "identity", aes(fill = Var2), colour = "black", alpha = 0.8)+
  scale_fill_manual(values = mycolors)+
  theme_classic()+
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=2, size = 16)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(size = 24),axis.title = element_text(size = 24))+
  ylab("Count")+
  xlab("")+        # Remove background elements manually
  theme(legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA))

ggsave(plot = p1, "FIGURES/IMMERSE paper/SCRIPT 4/S1 vs pre bar.png",width=5, height=6,units="in",dpi = 300)

count2 <- count %>% filter(Var1 ==  "State 2 vs pre")

p2 <- ggplot(count2, aes(y=Freq, x = Var2))+
  geom_bar(stat= "identity", aes(fill = Var2), colour = "black", alpha = 0.8)+
  scale_fill_manual(values = mycolors)+
  theme_classic()+
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=2, size = 16)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(size = 24),axis.title = element_text(size = 24))+
  ylab("Count")+
  xlab("")+        # Remove background elements manually
  theme(legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA))

ggsave(plot = p2, "FIGURES/IMMERSE paper/SCRIPT 4/S2 vs pre bar.png",width=5, height=6,units="in",dpi = 300)

count3 <- count %>% filter(Var1 ==  "State 3 vs pre")

p3 <- ggplot(count3, aes(y=Freq, x = Var2))+
  geom_bar(stat= "identity", aes(fill = Var2), colour = "black", alpha = 0.8)+
  scale_fill_manual(values = mycolors)+
  theme_classic()+
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=2, size = 16)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(size = 24),axis.title = element_text(size = 24))+
  ylab("Count")+
  xlab("")+        # Remove background elements manually
  theme(legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA))

ggsave(plot = p3, "FIGURES/IMMERSE paper/SCRIPT 4/S3 vs pre bar.png",width=5, height=6,units="in",dpi = 300)


count4<- count %>% filter(Var1 ==  "State 1 vs post")

p4 <- ggplot(count4, aes(y=Freq, x = Var2))+
  geom_bar(stat= "identity", aes(fill = Var2), colour = "black", alpha = 0.8)+
  scale_fill_manual(values = mycolors)+
  theme_classic()+
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=2, size = 16)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(size = 24),axis.title = element_text(size = 24))+
  ylab("Count")+
  xlab("")+        # Remove background elements manually
  theme(legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA))

ggsave(plot = p4, "FIGURES/IMMERSE paper/SCRIPT 4/State 1 vs postbar.png",width=5, height=6,units="in",dpi = 300)

count5<- count %>% filter(Var1 ==  "State 2 vs post")

p5 <- ggplot(count5, aes(y=Freq, x = Var2))+
  geom_bar(stat= "identity", aes(fill = Var2), colour = "black", alpha = 0.8)+
  scale_fill_manual(values = mycolors)+
  theme_classic()+
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=2, size = 16)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(size = 24),axis.title = element_text(size = 24))+
  ylab("Count")+
  xlab("")+        # Remove background elements manually
  theme(legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA))

ggsave(plot = p5, "FIGURES/IMMERSE paper/SCRIPT 4/State 2 vs postbar.png",width=5, height=6,units="in",dpi = 300)

count6<- count %>% filter(Var1 ==  "State 3 vs post")

p6 <- ggplot(count6, aes(y=Freq, x = Var2))+
  geom_bar(stat= "identity", aes(fill = Var2), colour = "black", alpha = 0.8)+
  scale_fill_manual(values = mycolors)+
  theme_classic()+
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=2, size = 16)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(size = 24),axis.title = element_text(size = 24))+
  ylab("Count")+
  xlab("")+        # Remove background elements manually
  theme(legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA))

ggsave(plot = p6, "FIGURES/IMMERSE paper/SCRIPT 4/State 3 vs postbar.png",width=5, height=6,units="in",dpi = 300)

count7<- count %>% filter(Var1 ==  "State 2 vs 1")

p7 <- ggplot(count7, aes(y=Freq, x = Var2))+
  geom_bar(stat= "identity", aes(fill = Var2), colour = "black", alpha = 0.8)+
  scale_fill_manual(values = mycolors)+
  theme_classic()+
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=2, size = 16)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(size = 24),axis.title = element_text(size = 24))+
  ylab("Count")+
  xlab("")+        # Remove background elements manually
  theme(legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA))

ggsave(plot = p7, "FIGURES/IMMERSE paper/SCRIPT 4/State 2 vs 1 bar.png",width=5, height=6,units="in",dpi = 300)

count8<- count %>% filter(Var1 ==  "State 3 vs 1")

p8 <- ggplot(count8, aes(y=Freq, x = Var2))+
  geom_bar(stat= "identity", aes(fill = Var2), colour = "black", alpha = 0.8)+
  scale_fill_manual(values = mycolors)+
  theme_classic()+
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=2, size = 16)+
  theme(legend.position = "none")+
  theme(axis.text = element_text(size = 24),axis.title = element_text(size = 24))+
  ylab("Count")+
  xlab("")+        # Remove background elements manually
  theme(legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA))

ggsave(plot = p8, "FIGURES/IMMERSE paper/SCRIPT 4/State 3 vs 1 bar.png",width=5, height=6,units="in",dpi = 300)



####################################################################################################################################
#####################################################################################################################################################################
####################################################################################################################################
#use genetonic for pathway analysis 
library(pcaExplorer)
library("AnnotationDbi")
library("topGO")
#add symbol to the rowData. Note this adds a dataframe in that section. So when using downstream will need to subset. 
#get the gene name
rownames(dds) <- substr(rownames(dds), 1, 15)
rowData(dds)$SYMBOL1 <- get_annotation_orgdb(dds, "org.Hs.eg.db", "ENSEMBL")
#add to dds 
rowData(dds)$SYMBOL <- as.factor(rowData(dds)$SYMBOL1$gene_name)

dds
state1vsCCT0_GT <- results(dds, contrast = c("State", "State 1", "Cardiac_T0"),lfcThreshold = 1, alpha = 0.05 )
state1vsCCT0_GT$SYMBOL <- rowData(dds)$SYMBOL
state1vsCCT1_GT <- results(dds, contrast = c("State", "State 1", "Cardiac_T1"),lfcThreshold = 1, alpha = 0.05 )
state1vsCCT1_GT$SYMBOL <- rowData(dds)$SYMBOL
CCT1vsCCT0_GT <- results(dds, contrast = c("State", "Cardiac_T1", "Cardiac_T0"),lfcThreshold = 1, alpha = 0.05 )
CCT1vsCCT0_GT$SYMBOL <- rowData(dds)$SYMBOL

state2vsCCT0_GT <- results(dds, contrast = c("State", "State 2", "Cardiac_T0"),lfcThreshold = 1, alpha = 0.05 )
state2vsCCT0_GT$SYMBOL <- rowData(dds)$SYMBOL

state2vsCCT1_GT <- results(dds, contrast = c("State", "State 2", "Cardiac_T1"),lfcThreshold = 1, alpha = 0.05 )
state2vsCCT1_GT$SYMBOL <- rowData(dds)$SYMBOL


state3vsCCT0_GT <- results(dds, contrast = c("State", "State 3", "Cardiac_T0"),lfcThreshold = 1, alpha = 0.05 )
state3vsCCT0_GT$SYMBOL <- rowData(dds)$SYMBOL

state3vsCCT1_GT <- results(dds, contrast = c("State", "State 3", "Cardiac_T1"),lfcThreshold = 1, alpha = 0.05 )
state3vsCCT1_GT$SYMBOL <- rowData(dds)$SYMBOL

state2vsstate1_GT <- results(dds, contrast = c("State", "State 2", "State 1"),lfcThreshold = 1, alpha = 0.05 )
state2vsstate1_GT$SYMBOL <- rowData(dds)$SYMBOL

state3vsstate1_GT <- results(dds, contrast = c("State", "State 3", "State 1"),lfcThreshold = 1, alpha = 0.05 )
state3vsstate1_GT$SYMBOL <- rowData(dds)$SYMBOL



#make the information for state1vsCCT0_GT
de_symbols_state1vsCCT0_GT <- deseqresult2df(state1vsCCT0_GT, FDR = 0.05)$SYMBOL
bg_ids <- rowData(dds)$SYMBOL[rowSums(counts(dds)) > 0]


topgoDE_state1vsCCT0_GT <-
  pcaExplorer::topGOtable(de_symbols_state1vsCCT0_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_state1vsCCT0_GT <- shake_topGOtableResult(topgoDE_state1vsCCT0_GT)
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
res_enrich_state1vsCCT0_GT <- get_aggrscores(res_enrich = res_enrich_state1vsCCT0_GT,
                                        res_de = state1vsCCT0_GT,
                                        annotation_obj = anno_df,
                                        aggrfun = mean)

write_csv(res_enrich_state1vsCCT0_GT, "Results_csv/IMMERSE paper/SCRIPT 4/genetonic/res_enrich_state1vsCCT0_GT.csv" )

#make the information for state1vsCCT1_GT
de_symbols_state1vsCCT1_GT <- deseqresult2df(state1vsCCT1_GT, FDR = 0.05)$SYMBOL
bg_ids <- rowData(dds)$SYMBOL[rowSums(counts(dds)) > 0]


topgoDE_state1vsCCT1_GT <-
  pcaExplorer::topGOtable(de_symbols_state1vsCCT1_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_state1vsCCT1_GT <- shake_topGOtableResult(topgoDE_state1vsCCT1_GT)
colnames(res_enrich_state1vsCCT1_GT)


#create the aggreegated scores
res_enrich_state1vsCCT1_GT <- get_aggrscores(res_enrich = res_enrich_state1vsCCT1_GT,
                                             res_de = state1vsCCT1_GT,
                                             annotation_obj = anno_df,
                                             aggrfun = mean)

write_csv(res_enrich_state1vsCCT1_GT, "Results_csv/IMMERSE paper/SCRIPT 4/genetonic/res_enrich_state1vsCCT1_GT.csv" )

##########################################################3
#make the information for CCT1vsCCT0_GT
de_symbols_CCT1vsCCT0_GT <- deseqresult2df(CCT1vsCCT0_GT, FDR = 0.05)$SYMBOL
bg_ids <- rowData(dds)$SYMBOL[rowSums(counts(dds)) > 0]


topgoDE_CCT1vsCCT0_GT <-
  pcaExplorer::topGOtable(de_symbols_CCT1vsCCT0_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_CCT1vsCCT0_GT <- shake_topGOtableResult(topgoDE_CCT1vsCCT0_GT)
colnames(res_enrich_state1vsCCT1_GT)


#create the aggreegated scores
res_enrich_CCT1vsCCT0_GT <- get_aggrscores(res_enrich = res_enrich_CCT1vsCCT0_GT,
                                             res_de = CCT1vsCCT0_GT,
                                             annotation_obj = anno_df,
                                             aggrfun = mean)

write_csv(res_enrich_CCT1vsCCT0_GT, "Results_csv/IMMERSE paper/SCRIPT 4/genetonic/res_enrich_CCT1vsCCT0_GT.csv" )

##########################################################3

#make the information for State2vsCCT0_GT
de_symbols_state2vsCCT0_GT <- deseqresult2df(state2vsCCT0_GT, FDR = 0.05)$SYMBOL

topgoDE_state2vsCCT0_GT <-
  pcaExplorer::topGOtable(de_symbols_state2vsCCT0_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_state2vsCCT0_GT <- shake_topGOtableResult(topgoDE_state2vsCCT0_GT)
colnames(res_enrich_state2vsCCT0_GT)

#create the aggreegated scores
res_enrich_state2vsCCT0_GT <- get_aggrscores(res_enrich = res_enrich_state2vsCCT0_GT,
                                           res_de = state2vsCCT0_GT,
                                           annotation_obj = anno_df,
                                           aggrfun = mean)
write_csv(res_enrich_state2vsCCT0_GT, "Results_csv/IMMERSE paper/SCRIPT 4/genetonic/res_enrich_state2vsCCT0_GT.csv" )

##########################################################3

#make the information for State2vsCCT1_GT
de_symbols_state2vsCCT1_GT <- deseqresult2df(state2vsCCT1_GT, FDR = 0.05)$SYMBOL

topgoDE_state2vsCCT1_GT <-
  pcaExplorer::topGOtable(de_symbols_state2vsCCT1_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_state2vsCCT1_GT <- shake_topGOtableResult(topgoDE_state2vsCCT1_GT)
colnames(res_enrich_state2vsCCT1_GT)

#create the aggreegated scores
res_enrich_state2vsCCT1_GT <- get_aggrscores(res_enrich = res_enrich_state2vsCCT1_GT,
                                             res_de = state2vsCCT1_GT,
                                             annotation_obj = anno_df,
                                             aggrfun = mean)
write_csv(res_enrich_state2vsCCT1_GT, "Results_csv/IMMERSE paper/SCRIPT 4/genetonic/res_enrich_state2vsCCT1_GT.csv" )

##########################################################3

#make the information for State3vsCCT0_GT
de_symbols_state3vsCCT0_GT <- deseqresult2df(state3vsCCT0_GT, FDR = 0.05)$SYMBOL

topgoDE_state3vsCCT0_GT <-
  pcaExplorer::topGOtable(de_symbols_state3vsCCT0_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_state3vsCCT0_GT <- shake_topGOtableResult(topgoDE_state3vsCCT0_GT)
colnames(res_enrich_state3vsCCT0_GT)


#create the aggreegated scores
res_enrich_state3vsCCT0_GT <- get_aggrscores(res_enrich = res_enrich_state3vsCCT0_GT,
                                           res_de = state3vsCCT0_GT,
                                           annotation_obj = anno_df,
                                           aggrfun = mean)
write_csv(res_enrich_state3vsCCT0_GT, "Results_csv/IMMERSE paper/SCRIPT 4/genetonic/res_enrich_state3vsCCT0_GT.csv" )

##########################################################3

#make the information for State3vsCCT1_GT
de_symbols_state3vsCCT1_GT <- deseqresult2df(state3vsCCT1_GT, FDR = 0.05)$SYMBOL

topgoDE_state3vsCCT1_GT <-
  pcaExplorer::topGOtable(de_symbols_state3vsCCT1_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_state3vsCCT1_GT <- shake_topGOtableResult(topgoDE_state3vsCCT1_GT)
colnames(res_enrich_state3vsCCT1_GT)


#create the aggreegated scores
res_enrich_state3vsCCT1_GT <- get_aggrscores(res_enrich = res_enrich_state3vsCCT1_GT,
                                             res_de = state3vsCCT1_GT,
                                             annotation_obj = anno_df,
                                             aggrfun = mean)
write_csv(res_enrich_state3vsCCT1_GT, "Results_csv/IMMERSE paper/SCRIPT 4/genetonic/res_enrich_state3vsCCT1_GT.csv" )

##########################################################3

#make the information for state2vsstate1_GT
de_symbols_state2vsstate1_GT <- deseqresult2df(state2vsstate1_GT, FDR = 0.05)$SYMBOL

topgoDE_state2vsstate1_GT <-
  pcaExplorer::topGOtable(de_symbols_state2vsstate1_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_state2vsstate1_GT <- shake_topGOtableResult(topgoDE_state2vsstate1_GT)
colnames(res_enrich_state3vsCCT1_GT)


#create the aggreegated scores
res_enrich_state2vsstate1_GT <- get_aggrscores(res_enrich = res_enrich_state2vsstate1_GT,
                                             res_de = state2vsstate1_GT,
                                             annotation_obj = anno_df,
                                             aggrfun = mean)
write_csv(res_enrich_state2vsstate1_GT, "Results_csv/IMMERSE paper/SCRIPT 4/genetonic/res_enrich_state2vsstate1_GT.csv" )

##########################################################3

#make the information for state3vsstate1_GT
de_symbols_state3vsstate1_GT <- deseqresult2df(state3vsstate1_GT, FDR = 0.05)$SYMBOL

topgoDE_state3vsstate1_GT <-
  pcaExplorer::topGOtable(de_symbols_state3vsstate1_GT,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Hs.eg.db",
                          geneID = "SYMBOL",
                          topTablerows = 500)

#convert for use in genetonic
res_enrich_state3vsstate1_GT <- shake_topGOtableResult(topgoDE_state3vsstate1_GT)
colnames(res_enrich_state3vsstate1_GT)


#create the aggreegated scores
res_enrich_state3vsstate1_GT <- get_aggrscores(res_enrich = res_enrich_state3vsstate1_GT,
                                               res_de = state3vsstate1_GT,
                                               annotation_obj = anno_df,
                                               aggrfun = mean)

write_csv(res_enrich_state3vsstate1_GT, "Results_csv/IMMERSE paper/SCRIPT 4/genetonic/res_enrich_state3vsstate1_GT.csv" )
#################################################################################
####################################################################################################

#plot heatmap by patient for state1 comparison

# get the genes into a dataframe so I can then add the annotation info
vst_dds <- vst(dds)
assay(vst_dds) <- limma::removeBatchEffect(assay(vst_dds), vst_dds$library_id) #remove batch effect

scores_mat <- gs_scores(
  se = vst_dds,
  res_de = state1vsCCT0_GT,
  res_enrich = res_enrich_state1vsCCT0_GT,
  annotation_obj = anno_df
)
scores_mat <- scores_mat[1:30,]
scores_df <- scores_mat %>% t() %>%  as.data.frame() %>% 
  rownames_to_column("sample_id")
#save this for integrated analysis
write_csv(scores_df, "Results_csv/IMMERSE paper/SCRIPT 4/integrated/pathways_top_30.csv")
#merge with the information from metadata
scores_df <- merge(samples, scores_df)

colnames(scores_df)

scores_df$State <- factor(scores_df$State, levels = c("Cardiac_T0", "Cardiac_T1", "State 1", "State 2","State 3" ))

scores_df <- scores_df %>% arrange(State)
colnames(scores_df)
hm_matrix <- scores_df %>% 
  dplyr::select(sample_id, 19:48) %>% 
  column_to_rownames("sample_id") %>%
  as.matrix() %>% t()


#draw label annotations for column (patients)
colnames(scores_df)
scores_df$Cohort_time
library(ComplexHeatmap)
anno <- scores_df %>% dplyr::select(sample_id,Cohort_time, State) %>% column_to_rownames("sample_id")
ha = HeatmapAnnotation(df = anno, 
                       col = list( Cohort_time = c("Sepsis_T0" = "#efe350ff", "Sepsis_T1" = "#f68f46ff","Sepsis_T2" = "#a65c85ff",
                                                  "Sepsis_T3" = "#403891ff","Cardiac_T0"= "#56A8CBFF", "Cardiac_T1" ="#DA291CFF"),
                                  State = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF",
                                            "Cardiac_T0"= "#56A8CBFF", "Cardiac_T1" ="#DA291CFF" )),
                       which = "col",
                       gp = gpar(col = "black"))

#make annotation rows (patients)



#ha = rowAnnotation(df = anno_df, 
 #                  col = list(Mygroup = c("Sepsis_T0" = "#efe350ff",
  #                                        "Sepsis_T1" = "#f68f46ff",
   #                                       "Sepsis_T2" = "#a65c85ff",
    #                                      "Sepsis_T3" = "#403891ff",
     #                                     "Cardiac_T0" ="#56A8CBFF",
      #                                    "Cardiac_T1" = "#DA291CFF", lwd = 2),
       #                       which = "row"))
ha
draw(ha, 1:95)

table(scores_df$State)
column_split = rep("Pre-Surgery", 258)

column_split[31:60] = "Post-surgery"
column_split[61:104] = "State 1"
column_split[105:200] = "State 2"
column_split[201:258] = "State 3"

column_split <- factor(column_split, levels = c("Pre-Surgery", "Post-surgery","State 1","State 2", "State 3" ))



library(circlize)

#sort out colour scheme#CC0001
col_fun = colorRamp2(c(-2, 0, 2), c("#0F4392", "white", "#CC0001"))
col_fun(seq(-2, 2))
# do heatmap for this comparison
hm <- Heatmap(hm_matrix, name = "Unknown",
              col = col_fun,
              cluster_columns = FALSE,
              cluster_rows = TRUE,
              show_heatmap_legend = TRUE,
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "ward",
              #bottom_annotation = columnAnnotation(mark=anno),
              #width = ncol(df3.1.1.1)*unit(5, "mm"), 
              #height = nrow(df3.1.1.1)*unit(1, "mm"),
              border = TRUE,
              row_km = 2,
              #column_km = 3,
              show_column_names = TRUE,
              heatmap_legend_param = list(title = "Z-score",legend_direction = "horizontal"),
              column_names_gp = grid::gpar(fontsize = 0),
              row_names_gp = grid::gpar(fontsize = 12),
              top_annotation = ha,
              row_names_side = c( "left"),
              row_dend_width = unit(20, "mm"),
              row_dend_side = c("right"),
              column_dend_height = unit(10, "mm"),
              column_split = column_split,
              column_title_gp = gpar(fontsize = 12, fontface = "bold"),
             width = ncol(hm_matrix)*unit(0.5, "mm"), 
             height = nrow(hm_matrix)*unit(6, "mm"))#+

#rowAnnotation(mark=row_anno)
hm


png("FIGURES/IMMERSE paper/SCRIPT 4/genetonic/HEATMAP_GO_sepsis states vs cardiac.png",width=15, height=12,units="in",res=500)
draw(hm, heatmap_legend_side="bottom", 
     legend_grouping = "original",
     merge_legends = T,padding = unit(c(0, 100, 0, 10), "mm"))
dev.off()

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

#Make a list of the pathways I want 
#The top 20 from each comparison. 
#CCT1vsCCT0
CCT1vsCCT0list <- res_enrich_CCT1vsCCT0_GT %>% slice(1:20)
CCT1vsCCT0list <- CCT1vsCCT0list$gs_description

#S1vsCCT0
S1vsCCT0list <- res_enrich_state1vsCCT0_GT %>% slice(1:20)
S1vsCCT0list <- S1vsCCT0list$gs_description

#S2vsCCT0
S2vsCCT0list <- res_enrich_state2vsCCT0_GT %>% slice(1:20)
S2vsCCT0list <- S2vsCCT0list$gs_description

#S3vsCCT0
S3vsCCT0list <- res_enrich_state3vsCCT0_GT %>% slice(1:20)
S3vsCCT0list <- S3vsCCT0list$gs_description

#S1vsCCT1
S1vsCCT1list <- res_enrich_state1vsCCT1_GT %>% slice(1:20)
S1vsCCT1list <- S1vsCCT1list$gs_description

#S2vsCCT1
S2vsCCT1list <- res_enrich_state2vsCCT1_GT %>% slice(1:20)
S2vsCCT1list <- S2vsCCT1list$gs_description

#S3vsCCT1
S3vsCCT1list <- res_enrich_state3vsCCT1_GT %>% slice(1:20)
S3vsCCT1list <- S3vsCCT1list$gs_description

#S2vsS1
S2vsS1list <- res_enrich_state2vsstate1_GT %>% slice(1:20)
S2vsS1list <- S2vsS1list$gs_description

#S3vsS1
S3vsS1list <- res_enrich_state3vsstate1_GT %>% slice(1:20)
S3vsS1list <- S3vsS1list$gs_description


pathwaylist <- c(CCT1vsCCT0list,S1vsCCT0list,S2vsCCT0list,S3vsCCT0list,S1vsCCT1list,S2vsCCT1list,S3vsCCT1list,S2vsS1list,S3vsS1list)
pathwaylist <- unique(pathwaylist)

# make summary heatmap for each group
colnames(res_enrich_CCT1vsCCT0_GT)
CCT1vsCCT0_hm <- res_enrich_CCT1vsCCT0_GT %>%
  dplyr::select(gs_description,gs_pvalue,gs_de_count,gs_bg_count,z_score) %>% 
  mutate(expression_ratio = gs_de_count/gs_bg_count) %>% 
  rename("pvalue_CCT1vsCCT0" = "gs_pvalue",
         'bgcount_CCT1vsCCT0' = "gs_bg_count",
         "decount_CCT1vsCCT0" = "gs_de_count",
         'Zscore_CCT1vsCCT0' = "z_score",
         'expressionratio_CCT1vsCCT0' = "expression_ratio")
  
S1vsCCT0_hm <- res_enrich_state1vsCCT0_GT %>%
  dplyr::select(gs_description,gs_pvalue,gs_de_count,gs_bg_count,z_score) %>% 
  mutate(expression_ratio = gs_de_count/gs_bg_count) %>% 
  rename("pvalue_S1vsCCT0" = "gs_pvalue",
         'bgcount_S1vsCCT0' = "gs_bg_count",
         "decount_S1vsCCT0" = "gs_de_count",
         'Zscore_S1vsCCT0' = "z_score",
         'expressionratio_S1vsCCT0' ="expression_ratio")

S2vsCCT0_hm <- res_enrich_state2vsCCT0_GT %>%
  dplyr::select(gs_description,gs_pvalue,gs_de_count,gs_bg_count,z_score) %>% 
  mutate(expression_ratio = gs_de_count/gs_bg_count) %>% 
  rename("pvalue_S2vsCCT0" = "gs_pvalue",
         'bgcount_S2vsCCT0' = "gs_bg_count",
         "decount_S2vsCCT0" = "gs_de_count",
         'Zscore_S2vsCCT0' = "z_score",
         'expressionratio_S2vsCCT0' ="expression_ratio")

S3vsCCT0_hm <- res_enrich_state3vsCCT0_GT %>%
  dplyr::select(gs_description,gs_pvalue,gs_de_count,gs_bg_count,z_score) %>% 
  mutate(expression_ratio = gs_de_count/gs_bg_count) %>% 
  rename("pvalue_S3vsCCT0" = "gs_pvalue",
         'bgcount_S3vsCCT0' = "gs_bg_count",
         "decount_S3vsCCT0" = "gs_de_count",
         'Zscore_S3vsCCT0' = "z_score",
         'expressionratio_S3vsCCT0' ="expression_ratio")

S1vsCCT1_hm <- res_enrich_state1vsCCT1_GT %>%
  dplyr::select(gs_description,gs_pvalue,gs_de_count,gs_bg_count,z_score) %>% 
  mutate(expression_ratio = gs_de_count/gs_bg_count) %>% 
  rename("pvalue_S1vsCCT1" = "gs_pvalue",
         'bgcount_S1vsCCT1' = "gs_bg_count",
         "decount_S1vsCCT1" = "gs_de_count",
         'Zscore_S1vsCCT1' = "z_score",
         'expressionratio_S1vsCCT1' ="expression_ratio")

S2vsCCT1_hm <- res_enrich_state2vsCCT1_GT %>%
  dplyr::select(gs_description,gs_pvalue,gs_de_count,gs_bg_count,z_score) %>% 
  mutate(expression_ratio = gs_de_count/gs_bg_count) %>% 
  rename("pvalue_S2vsCCT1" = "gs_pvalue",
         'bgcount_S2vsCCT1' = "gs_bg_count",
         "decount_S2vsCCT1" = "gs_de_count",
         'Zscore_S2vsCCT1' = "z_score",
         'expressionratio_S2vsCCT1' ="expression_ratio")

S3vsCCT1_hm <- res_enrich_state3vsCCT1_GT %>%
  dplyr::select(gs_description,gs_pvalue,gs_de_count,gs_bg_count,z_score) %>% 
  mutate(expression_ratio = gs_de_count/gs_bg_count) %>% 
  rename("pvalue_S3vsCCT1" = "gs_pvalue",
         'bgcount_S3vsCCT1' = "gs_bg_count",
         "decount_S3vsCCT1" = "gs_de_count",
         'Zscore_S3vsCCT1' = "z_score",
         'expressionratio_S3vsCCT1' ="expression_ratio")

S2vsS1_hm <- res_enrich_state2vsstate1_GT %>%
  dplyr::select(gs_description,gs_pvalue,gs_de_count,gs_bg_count,z_score) %>% 
  mutate(expression_ratio = gs_de_count/gs_bg_count) %>% 
  rename("pvalue_S2vsS1" = "gs_pvalue",
         'bgcount_S2vsS1' = "gs_bg_count",
         "decount_S2vsS1" = "gs_de_count",
         'Zscore_S2vsS1' = "z_score",
         'expressionratio_S2vsS1' ="expression_ratio")

S3vsS1_hm <- res_enrich_state3vsstate1_GT %>%
  dplyr::select(gs_description,gs_pvalue,gs_de_count,gs_bg_count,z_score) %>% 
  mutate(expression_ratio = gs_de_count/gs_bg_count) %>% 
  rename("pvalue_S3vsS1" = "gs_pvalue",
         'bgcount_S3vsS1' = "gs_bg_count",
         "decount_S3vsS1" = "gs_de_count",
         'Zscore_S3vsS1' = "z_score",
         'expressionratio_S3vsS1' ="expression_ratio")



colnames(hm_merge)
hm_merge <- merge(CCT1vsCCT0_hm,S1vsCCT0_hm, all = T)
hm_merge <- merge(hm_merge,S2vsCCT0_hm, all = T)
hm_merge <- merge(hm_merge,S3vsCCT0_hm, all = T)
hm_merge <- merge(hm_merge,S1vsCCT1_hm, all = T)
hm_merge <- merge(hm_merge,S2vsCCT1_hm, all = T)
hm_merge <- merge(hm_merge,S3vsCCT1_hm, all = T)
hm_merge <- merge(hm_merge,S2vsS1_hm, all = T)
hm_merge <- merge(hm_merge,S3vsS1_hm, all = T)%>% 
  mutate(Zscore_CCT1vsCCT0 = coalesce(Zscore_CCT1vsCCT0, 0),
         Zscore_S1vsCCT0 = coalesce(Zscore_S1vsCCT0, 0),
         Zscore_S2vsCCT0 = coalesce(Zscore_S2vsCCT0, 0),
         Zscore_S3vsCCT0 = coalesce(Zscore_S3vsCCT0, 0),
         Zscore_S1vsCCT1 = coalesce(Zscore_S1vsCCT1, 0),
         Zscore_S2vsCCT1 = coalesce(Zscore_S2vsCCT1, 0),
         Zscore_S3vsCCT1 = coalesce(Zscore_S3vsCCT1, 0),
         Zscore_S2vsS1 = coalesce(Zscore_S2vsS1, 0),
         Zscore_S3vsS1 = coalesce(Zscore_S3vsS1, 0),
         
         pvalue_CCT1vsCCT0 = coalesce(pvalue_CCT1vsCCT0, 1),
         pvalue_S1vsCCT0 = coalesce(pvalue_S1vsCCT0, 1),
         pvalue_S2vsCCT0 = coalesce(pvalue_S2vsCCT0, 1),
         pvalue_S3vsCCT0 = coalesce(pvalue_S3vsCCT0, 1),
         pvalue_S1vsCCT1 = coalesce(pvalue_S1vsCCT1, 1),
         pvalue_S2vsCCT1 = coalesce(pvalue_S2vsCCT1, 1),
         pvalue_S3vsCCT1 = coalesce(pvalue_S3vsCCT1, 1),
         pvalue_S2vsS1 = coalesce(pvalue_S2vsS1, 1),
         pvalue_S3vsS1 = coalesce(pvalue_S3vsS1, 1)) %>% 
  filter(gs_description %in% pathwaylist)

colnames(hm_merge)
hm_merge1 <- hm_merge %>%  pivot_longer(cols = c(2:46),
                            names_to = c("statistic", "comparison"),
                            names_sep = "_",
                            values_to = "value") %>% 
  spread(statistic,value)


# make the matrix for the Z-score 
colnames(hm_merge1)
# make the matrix for the log fold change 
matrixZ <- hm_merge1 %>% 
  dplyr::select(1,2,7) %>%
  spread(comparison, Zscore) %>% 
  column_to_rownames("gs_description") %>% 
  dplyr::select(CCT1vsCCT0, S1vsCCT0, S2vsCCT0,S3vsCCT0,S1vsCCT1,S2vsCCT1,S3vsCCT1,S2vsS1,S3vsS1) %>% 
   as.matrix()

# make the matrix for the log fold change 
colnames(hm_merge1)
matrix_pvalue1 <- hm_merge1 %>% 
  dplyr::select(1,2,6)  %>% 
  spread(comparison, pvalue)%>% 
  column_to_rownames("gs_description")%>% 
  dplyr::select(CCT1vsCCT0, S1vsCCT0, S2vsCCT0,S3vsCCT0,S1vsCCT1,S2vsCCT1,S3vsCCT1,S2vsS1,S3vsS1) %>% 
  as.matrix()
class(matrix_pvalue1)


col_fun = colorRamp2(c(-3, 0, 3), c("#4575B4", "white", "#D73027"))
col_fun(seq(-3, 3))


column_split = rep("Surgery", 9)
column_split[2:4] = "States vs. pre-surgery"
column_split[5:7] = "States vs. post-surgery"
column_split[8:9] = "vs. State 1"



col.subsections <- c(1, 3, 3, 2)
col_split = data.frame(rep(c("Surgery", "Vs. Admission", "Vs. Post-surgery", "Vs. Pre-surgery"), col.subsections))

hm <- Heatmap(matrixZ,name = "Z-score", cluster_columns=FALSE,cluster_rows=T,border = F, col = col_fun,
        column_split = col_split,
        column_gap=unit(.02, "npc"),
        row_names_side = c( "left"),
        row_dend_width = unit(20, "mm"),
        row_dend_side = c("right"), 
        column_dend_height = unit(10, "mm"),
        heatmap_legend_param = list(title = "Z-score",legend_direction = "horizontal"),
        column_names_gp = grid::gpar(fontsize = 10),
        row_names_gp = grid::gpar(fontsize = 14),
        width = ncol(matrixZ)*unit(6, "mm"), 
        height = nrow(matrixZ)*unit(5, "mm"),
        row_gap=unit(.02, "npc"),
        #row_split = row_split,
        cell_fun = function(j, i, x, y, w, h, fill) {
          grid.rect(x = x, y = y, width = w, height = h, 
                    gp = gpar(col = "black", fill = NA))
          if(matrix_pvalue1[i, j] < 0.001) {
            grid.text("***", x, y, just = "centre", vjust = 0.8)
          } else if(matrix_pvalue1[i, j] < 0.01) {
            grid.text("**", x, y, just = "centre", vjust = 0.8)
          } else if(matrix_pvalue1[i, j] < 0.05) {
            grid.text("*", x, y, just = "centre", vjust = 0.8)
          }
        })
hm


png("FIGURES/IMMERSE paper/SCRIPT 4/genetonic/genetonic summary heatmap.png",width=15, height=23,units="in",res=500)
draw(hm,heatmap_legend_side="bottom", 
     legend_grouping = "original",
     merge_legends = T,padding = unit(c(0, 120, 0, -70), "mm"))
dev.off()
####################################################################################################

####################################################################################################
####################################################################################################
# Make line plots to show the changes in pathways over the course of pseudotime
#make plots to show resolution vs tolerance 
#first use the scores mmatrix from SIS1 vs Cardiac T0. This should bring up pathways that I need here

colnames(res_enrich_state1vsCCT0_GT)

scores_mat_S1vsCCT0 <- gs_scores(
  se = vst_dds,
  res_de = State1vsCCT0_GT,
  res_enrich = res_enrich_state1vsCCT0_GT[1:500,],
  annotation_obj = anno_df
)


scores_mat_S1vsCCT0 <- scores_mat_S1vsCCT0 %>% t() %>%  as.data.frame() %>% 
  rownames_to_column("sample_id")

T_pathways <- scores_mat_S1vsCCT0 %>% 
  dplyr::select("sample_id",
                "response to wounding|GO:0009611",
                "extracellular matrix organization|GO:0030198",
                "extracellular matrix disassembly|GO:0022617",
                "cellular response to oxygen-containing compound|GO:1901701", 
                "regeneration|GO:0031099",
                "animal organ regeneration|GO:0031100")




T_pathways <- merge(T_pathways, samples, all.x =TRUE)


colnames(T_pathways)
T_pathways <- T_pathways %>% gather(pathway, value, 2:7)
T_pathways <- T_pathways%>% tidyr::separate(pathway, into = c('Pathway', 'GO number'), sep = "\\|GO") %>% dplyr::select(-`GO number`) %>%  spread(Pathway, value)

colnames(T_pathways)


T_pathways <- T_pathways %>% gather(Pathway, Z_score, 19:24)



tolerance <- T_pathways %>% group_by(Pathway) %>% filter(State == "State 1"  |
                                                           State == "State 2"|
                                                           State == "State 3") %>% 
  do(augment(loess(Z_score~Lineage1, .))) %>% 
  {ggplot(., aes(x = Lineage1, y= Z_score, label = Pathway)) + 
      geom_vline(xintercept = 0.03, colour= "grey")+
      geom_vline(xintercept = 0.15, colour= "grey")+
      stat_smooth(aes(group = Pathway, linetype = Pathway, label = Pathway),method = "loess", se = FALSE, colour =  "lightblue", size =2)+
      stat_smooth(colour = "black", size = 2)+
      guides(color = F) +
      geom_label_repel(data = filter(., Lineage1 == max(Lineage1)),
                       aes(Lineage1, .fitted), size = 9, nudge_x = 0.1,nudge_y = 0.6, max.overlaps = 20)+
      theme_classic()+
      theme(legend.position = "none")+
      theme(text = element_text(size = 24))+
      xlab("Pseudotime")}


ggsave(plot = tolerance,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/Tolerance.png", width = 10, height = 7, dpi = 300 )

#plot by the PHATE plot
library(heatmaply)
colnames(T_pathways)
phate_woundhealing <- T_pathways  %>% filter(State == "State 1"  |
                                              State == "State 2"|
                                              State == "State 3") %>% 
  spread(Pathway, Z_score) %>% 
  ggplot( aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = `response to wounding`), pch = 21, size = 12, colour = "black")+
  scale_fill_viridis_c()+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = 24))
ggsave(plot = phate_woundhealing,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/phate/response to wounding.png", width =7, height = 7, dpi = 300 )

#plot the boxplot of each condition
woundhealing_boxpot <- T_pathways %>% 
  spread(Pathway, Z_score) %>% 
  ggplot( aes(x = State, y= `response to wounding`))+
  geom_boxplot()+
  geom_quasirandom(aes(fill = State), pch = 21, size = 8, colour = "black")+
  scale_fill_manual( values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF",
                                "Cardiac_T0"= "#56A8CBFF", "Cardiac_T1" ="#DA291CFF"))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(text = element_text(size = 24))

ggsave(plot = woundhealing_boxpot,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/BP/response to wounding.png", width = 6, height = 6, dpi = 300 )


matrix <- assay(vst_dds) %>% t() %>% as.data.frame() %>% rownames_to_column("sample_id") 
matrix <- merge(samples, matrix)

#PDGFB inhibits MMPs, anti apoptotic and promote proliferation
PDGFB <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000100311)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())
ggsave(plot = PDGFB,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/wound healing/PDGFB.png", width = 5, height = 5, dpi = 300 )

#TGFA can promote proliferation or differrntiation of cells
IL10 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000136634)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())
ggsave(plot = IL10,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/wound healing/IL10.png", width = 5, height = 5, dpi = 300 )

#IGF1
IGF1 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000017427)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())
ggsave(plot = IGF1,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/wound healing/IGF1.png", width = 5, height = 5, dpi = 300 )
#GATA2 Regulates developments and proliferation of hemeopoeitc cells and endocrine cells
GATA2 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000179348)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())
ggsave(plot = GATA2,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/wound healing/GATA2.png", width = 5, height = 5, dpi = 300 )
####################################################################################################
####################################################################################################
####################################################################################################

R_pathways <- scores_mat_S1vsCCT0 %>% 
  dplyr::select("sample_id",
                "pyroptosis|GO:0070269",
                "acute inflammatory response|GO:0002526",
                "positive regulation of interferon-gamma production|GO:0032729",
                "cellular response to interferon-gamma|GO:0071346")




R_pathways <- merge(R_pathways, samples, all.x =TRUE)


colnames(R_pathways)
R_pathways <- R_pathways %>% gather(pathway, value, 2:5)
R_pathways <- R_pathways%>% tidyr::separate(pathway, into = c('Pathway', 'GO number'), sep = "\\|GO") %>% dplyr::select(-`GO number`) %>%  spread(Pathway, value)

colnames(R_pathways)


R_pathways <- R_pathways %>% gather(Pathway, Z_score, 19:22)


resistance <- R_pathways %>% group_by(Pathway)%>% filter(State == "State 1"  |
                                                           State == "State 2"|
                                                           State == "State 3") %>% 
  do(augment(loess(Z_score~Lineage1, .))) %>% 
  {ggplot(., aes(x = Lineage1, y= Z_score, label = Pathway)) + 
      geom_vline(xintercept = 0.03, colour= "grey")+
      geom_vline(xintercept = 0.15, colour= "grey")+
      
      stat_smooth(aes(group = Pathway, linetype = Pathway, label = Pathway),method = "loess", se = FALSE, colour =  "#FF4F4B", size =2)+
      stat_smooth(colour = "black", size = 2)+
      guides(color = F) +
      geom_label_repel(data = filter(., Lineage1 == max(Lineage1)),
                       aes(Lineage1, .fitted), size = 9, nudge_x = 0.06, max.overlaps = 20)+
      theme_classic()+
      theme(legend.position = "none")+
      theme(text = element_text(size = 24))+
      xlab("Pseudotime")}



ggsave(plot = resistance,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/Resistance.png", width = 10, height = 7, dpi = 300 )

#plot by the PHATE plot
library(heatmaply)
colnames(cellkilling_pathways)
phate_acuteinflam <- R_pathways  %>% filter(State == "State 1"  |
                                                        State == "State 2"|
                                                        State == "State 3") %>% 
  spread(Pathway, Z_score) %>% 
  ggplot( aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = `acute inflammatory response`), pch = 21, size = 12, colour = "black")+
  scale_fill_viridis_c()+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = 24))
ggsave(plot = phate_acuteinflam,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/phate/Acute inflammatory response.png", width =7, height = 7, dpi = 300 )

#plot the boxplot of each condition
acuteinflam_boxpot <- R_pathways %>% 
  spread(Pathway, Z_score) %>% 
  ggplot( aes(x = State, y= `acute inflammatory response`))+
  geom_boxplot()+
  geom_quasirandom(aes(fill = State), pch = 21, size = 8, colour = "black")+
  scale_fill_manual( values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF",
                                "Cardiac_T0"= "#56A8CBFF", "Cardiac_T1" ="#DA291CFF"))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(text = element_text(size = 24))

ggsave(plot = acuteinflam_boxpot,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/BP/acute inflammatory response.png", width = 6, height = 6, dpi = 300 )


#ALOX5AP ENSG00000163221
ALOX5AP <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000132965)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())
ggsave(plot = ALOX5AP,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/acute inflammatory response/ALOX5AP.png", width = 5, height = 5, dpi = 300 )
#C2CD4B ENSG00000205502
C2CD4B <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000205502)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())
ggsave(plot = C2CD4B,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/acute inflammatory response/C2CD4B.png", width = 5, height = 5, dpi = 300 )


#CD163 ENSG00000177575
#acute phase-regulated receptor involved in the clearance and endocytosis of hemoglobin/haptoglobin complexes by macrophages
# immune sensor for bacteria and inducer of inflammation
CD163 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000177575)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())
ggsave(plot = CD163,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/acute inflammatory response/CD163.png", width = 5, height = 5, dpi = 300 )
#CEBPB ENSG00000172216
CEBPB <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000172216)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())
ggsave(plot = CEBPB,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/acute inflammatory response/CEBPB.png", width = 5, height = 5, dpi = 300 )

#ORM1 ENSG00000229314 
ORM1 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000229314)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())
ggsave(plot = ORM1,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/acute inflammatory response/ORM1.png", width = 5, height = 5, dpi = 300 )


#SERPINA1 ENSG00000197249  
SERPINA1 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000197249)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())
ggsave(plot = SERPINA1,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/acute inflammatory response/SERPINA1.png", width = 5, height = 5, dpi = 300 )


#OSM ENSG00000099985  
OSM <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000099985)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())
ggsave(plot = OSM,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/acute inflammatory response/OSM.png", width = 5, height = 5, dpi = 300 )

############################################################################################################################################################
############################################################################################################################################################
# Do bigger celldeath plot
celldeath_pathways <- scores_mat_S1vsCCT0 %>% 
  dplyr::select("sample_id",
                "positive regulation of endothelial cell apoptotic process|GO:2000353",
                "apoptotic process|GO:0006915",
                "negative regulation of T cell apoptotic process|GO:0070233",
                "regulation of apoptotic process|GO:0042981")




celldeath_pathways <- merge(celldeath_pathways, samples, all.x =TRUE)

celldeath_pathways <- celldeath_pathways
colnames(celldeath_pathways)
celldeath_pathways <- celldeath_pathways %>% gather(pathway, value, 2:5)
celldeath_pathways <- celldeath_pathways%>% tidyr::separate(pathway, into = c('Pathway', 'GO number'), sep = "\\|GO") %>% dplyr::select(-`GO number`) %>%  spread(Pathway, value)

colnames(celldeath_pathways)


celldeath_pathways <- celldeath_pathways %>% gather(Pathway, Z_score, 19:22)


celldeath <- celldeath_pathways %>% group_by(Pathway) %>% filter(State == "State 1"  |
                                                                   State == "State 2"|
                                                                   State == "State 3") %>% 
  do(augment(loess(Z_score~Lineage1, .))) %>% 
  {ggplot(., aes(x = Lineage1, y= Z_score, label = Pathway)) + 
      geom_vline(xintercept = 0.03, colour= "grey")+
      geom_vline(xintercept = 0.15, colour= "grey")+
      stat_smooth(aes(group = Pathway, linetype = Pathway, label = Pathway),method = "loess", se = FALSE, colour =  "coral", size =2)+
      stat_smooth(colour = "black", size = 2)+
      guides(color = F) +
      geom_label_repel(data = filter(., Lineage1 == max(Lineage1)),
                       aes(Lineage1, .fitted), size = 9, nudge_x = 0.06, max.overlaps = 20)+
      theme_classic()+
      theme(legend.position = "none")+
      theme(text = element_text(size = 24))+
      xlab("Pseudotime")}


ggsave(plot = celldeath,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/celldeath.png", width = 10, height = 7, dpi = 300 )


#phate plot
phate_celldeath <- celldeath_pathways  %>% filter(State == "State 1"  |
                                                        State == "State 2"|
                                                        State == "State 3") %>% 
  spread(Pathway, Z_score) %>% 
  ggplot( aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = `apoptotic process`), pch = 21, size = 12, colour = "black")+
  scale_fill_viridis_c()+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = 24))
ggsave(plot = phate_celldeath,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/phate/apoptotic process.png", width =7, height = 7, dpi = 300 )

#plot the boxplot of each condition

celldeath_boxplot <- celldeath_pathways %>% 
  spread(Pathway, Z_score) %>% 
  ggplot( aes(x = State, y= `apoptotic process`))+
  geom_boxplot()+
  geom_quasirandom(aes(fill = State), pch = 21, size = 8, colour = "black")+
  scale_fill_manual( values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF",
                                "Cardiac_T0"= "#56A8CBFF", "Cardiac_T1" ="#DA291CFF"))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(text = element_text(size = 24))

ggsave(plot = celldeath_boxplot,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/BP/apoptotic process.png", width = 6, height = 6, dpi = 300 )

#PDCD1 ENSG00000188389
PDCD1 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000188389)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())
ggsave(plot = PDCD1,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/positive regulation of apoptotic process/PDCD1.png", width = 5, height = 5, dpi = 300 )

#NLRC4 ENSG00000091106
NLRC4 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000091106)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())
ggsave(plot = NLRC4,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/positive regulation of apoptotic process/NLRC4.png", width = 5, height = 5, dpi = 300 )

#LDHA  ENSG00000134333
LDHA <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize( ENSG00000134333)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())
ggsave(plot = LDHA,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/positive regulation of apoptotic process/LDHA.png", width = 5, height = 5, dpi = 300 )

#FASLG  ENSG00000117560
FASLG <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize( ENSG00000117560)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())
ggsave(plot = FASLG,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/positive regulation of apoptotic process/FASLG.png", width = 5, height = 5, dpi = 300 )

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

# Do bigger cellkilling plot
cellkilling_pathways <- scores_mat_S1vsCCT0 %>% 
  dplyr::select("sample_id",
                "killing of cells of other organism|GO:0031640",
                "killing by host of symbiont cells|GO:0051873",
                "neutrophil-mediated killing of symbiont cell|GO:0070943")




cellkilling_pathways <- merge(cellkilling_pathways, samples, all.x =TRUE)

cellkilling_pathways <- cellkilling_pathways
colnames(cellkilling_pathways)
cellkilling_pathways <- cellkilling_pathways %>% gather(pathway, value, 2:4)
cellkilling_pathways <- cellkilling_pathways%>% tidyr::separate(pathway, into = c('Pathway', 'GO number'), sep = "\\|GO") %>% dplyr::select(-`GO number`) %>%  spread(Pathway, value)

colnames(cellkilling_pathways)


cellkilling_pathways <- cellkilling_pathways %>% gather(Pathway, Z_score, 19:21)

library(broom)

cellkilling <- cellkilling_pathways %>% filter(State == "State 1"  |
                                                 State == "State 2"|
                                                 State == "State 3") %>% group_by(Pathway) %>% 
  do(augment(loess(Z_score~Lineage1, .))) %>% 
  {ggplot(., aes(x = Lineage1, y= Z_score, label = Pathway)) + 
      geom_vline(xintercept = 0.03, colour= "grey")+
      geom_vline(xintercept = 0.15, colour= "grey")+
      stat_smooth(aes(group = Pathway, linetype = Pathway, label = Pathway),method = "loess", se = FALSE, colour =  "coral", size =2)+
      stat_smooth(colour = "black", size = 2)+
      guides(color = F) +
      geom_label_repel(data = filter(., Lineage1 == max(Lineage1)),
                       aes(Lineage1, .fitted), size = 9, nudge_x = 0.06, max.overlaps = 20)+
      theme_classic()+
      theme(legend.position = "none")+
      theme(text = element_text(size = 24))+
      xlab("Pseudotime")}


ggsave(plot = cellkilling,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/cellkilling.png", width = 10, height = 7, dpi = 300 )

#plot by the PHATE plot
library(heatmaply)
colnames(cellkilling_pathways)
phate_cellkilling <- cellkilling_pathways  %>% filter(State == "State 1"  |
                                       State == "State 2"|
                                       State == "State 3") %>% 
  spread(Pathway, Z_score) %>% 
  ggplot( aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = `killing of cells of other organism`), pch = 21, size = 12, colour = "black")+
  scale_fill_viridis_c()+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = 24))
ggsave(plot = phate_cellkilling,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/phate/killing of cells of other organism.png", width =7, height = 7, dpi = 300 )

library(ggbeeswarm)
#plot the boxplot of each condition
KC1_boxpot <- cellkilling_pathways %>% 
  spread(Pathway, Z_score) %>% 
  ggplot( aes(x = State, y= `killing of cells of other organism`))+
  geom_boxplot()+
  geom_quasirandom(aes(fill = State), pch = 21, size = 8, colour = "black")+
  scale_fill_manual( values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF",
"Cardiac_T0"= "#56A8CBFF", "Cardiac_T1" ="#DA291CFF"))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(text = element_text(size = 24))

ggsave(plot = KC1_boxpot,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/BP/killing of cells of other organism.png", width = 6, height = 6, dpi = 300 )


#take two genes. 
samples
matrix <- assay(vst_dds) %>% t() %>% as.data.frame() %>% rownames_to_column("sample_id") 
matrix <- merge(samples, matrix)

#S100A12 ENSG00000163221
S100A12 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000163221)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())
ggsave(plot = S100A12,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/killing of cells of other organism/S100A12.png", width = 5, height = 5, dpi = 300 )
  
#PRF1 ENSG00000180644
PRF1 <-   ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
    geom_point(aes(fill = normalize(ENSG00000180644)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())
ggsave(plot = PRF1,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/killing of cells of other organism/PRF1.png", width = 5, height = 5, dpi = 300 )

#GZMA
GZMA <-   ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
    geom_point(aes(fill = normalize(ENSG00000145649)), pch = 21, size = 12, colour = "black")+
    theme_minimal()+
    scale_fill_viridis_c(option = "A")+
    theme(legend.position = "none")+
    theme(text = element_text(size = 26))+
    theme( axis.line=element_blank(),
           axis.text.x=element_blank(),
           axis.text.y=element_blank(),
           axis.ticks=element_blank(),
           axis.title.x=element_blank(),
           axis.title.y=element_blank(),
           legend.position="none",
           panel.background=element_blank(),
           panel.border=element_blank(),
           panel.grid.major=element_blank(),
           panel.grid.minor=element_blank(),
           plot.background=element_blank())
ggsave(plot = GZMA,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/killing of cells of other organism/GZMA.png", width = 5, height = 5, dpi = 300 )
  
# GNLY ENSG00000115523 released upon antigen stimulation.
#This protein is present in cytotoxic granules of cytotoxic T lymphocytes and natural killer cells, and it has antimicrobial

GNLY <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
    geom_point(aes(fill = normalize(ENSG00000115523)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())

 ggsave(plot = GNLY,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/killing of cells of other organism/GNLY.png", width = 5, height = 5, dpi = 300 )
    
#ENSG00000197561 ELANE hydrolyzes proteins within specialized neutrophil lysosomes, called azurophil granules. Neutropenia. 
ELANE <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
      geom_point(aes(fill = normalize(ENSG00000197561)), pch = 21, size = 12, colour = "black")+
      theme_minimal()+
      scale_fill_viridis_c(option = "A")+
    theme(legend.position = "none")+
    theme(text = element_text(size = 26))+
 theme( axis.line=element_blank(),
axis.text.x=element_blank(),
axis.text.y=element_blank(),
axis.ticks=element_blank(),
axis.title.x=element_blank(),
axis.title.y=element_blank(),
legend.position="none",
panel.background=element_blank(),
panel.border=element_blank(),
panel.grid.major=element_blank(),
panel.grid.minor=element_blank(),
plot.background=element_blank())

ggsave(plot = ELANE,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/killing of cells of other organism/ELANE.png", width = 5, height = 5, dpi = 300 )
    
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

# Do bigger innate plot
innate_pathways <- scores_mat_S1vsCCT0 %>% 
  dplyr::select("sample_id",
                "neutrophil activation|GO:0042119",
                "neutrophil degranulation|GO:0043312",
                "positive regulation of monocyte differentiation|GO:0045657",
                "dendritic cell antigen processing and presentation|GO:0002468",
                "antigen processing and presentation of endogenous antigen|GO:0019883")




innate_pathways <- merge(innate_pathways, samples, all.x =TRUE)


colnames(innate_pathways)
innate_pathways <- innate_pathways %>% gather(pathway, value, 2:6)
innate_pathways <- innate_pathways%>% tidyr::separate(pathway, into = c('Pathway', 'GO number'), sep = "\\|GO") %>% dplyr::select(-`GO number`) %>%  spread(Pathway, value)

colnames(innate_pathways)


innate_pathways <- innate_pathways %>% gather(Pathway, Z_score, 19:23)

library(broom)

innate <- innate_pathways %>% filter(State == "State 1"  |
                                                 State == "State 2"|
                                                 State == "State 3") %>% group_by(Pathway) %>% 
  do(augment(loess(Z_score~Lineage1, .))) %>% 
  {ggplot(., aes(x = Lineage1, y= Z_score, label = Pathway)) + 
      geom_vline(xintercept = 0.03, colour= "grey")+
      geom_vline(xintercept = 0.15, colour= "grey")+
      stat_smooth(aes(group = Pathway, linetype = Pathway, label = Pathway),method = "loess", se = FALSE, colour =  "#FDDA0D", size =2)+
      stat_smooth(colour = "black", size = 2)+
      guides(color = F) +
      geom_label_repel(data = filter(., Lineage1 == max(Lineage1)),
                       aes(Lineage1, .fitted), size = 9, nudge_x = 0.09,nudge_y =0.6, max.overlaps = 20)+
      theme_classic()+
      theme(legend.position = "none")+
      theme(text = element_text(size = 24))+
      xlab("Pseudotime")}


ggsave(plot = innate,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/neut activation.png", width = 10, height = 7, dpi = 300 )

#plot by the PHATE plot
library(heatmaply)
colnames(cellkilling_pathways)
phate_cellinnate <- innate_pathways  %>% filter(State == "State 1"  |
                                                        State == "State 2"|
                                                        State == "State 3") %>% 
  spread(Pathway, Z_score) %>% 
  ggplot( aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = `neutrophil activation`), pch = 21, size = 12, colour = "black")+
  scale_fill_viridis_c()+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = 24))
ggsave(plot = phate_cellinnate,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/phate/neutrophil activation.png", width =7, height = 7, dpi = 300 )

library(ggbeeswarm)
#plot the boxplot of each condition
innate_boxpot <- innate_pathways %>% 
  spread(Pathway, Z_score) %>% 
  ggplot( aes(x = State, y= `neutrophil activation`))+
  geom_boxplot()+
  geom_quasirandom(aes(fill = State), pch = 21, size = 8, colour = "black")+
  scale_fill_manual( values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF",
                                "Cardiac_T0"= "#56A8CBFF", "Cardiac_T1" ="#DA291CFF"))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(text = element_text(size = 24))

ggsave(plot = innate_boxpot,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/BP/neutrophil activation.png", width = 6, height = 6, dpi = 300 )


CD177 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000204936)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())

ggsave(plot = CD177,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/neutrophil activation/CD177.png", width = 5, height = 5, dpi = 300 )

ANXA3 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000138772)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())

ggsave(plot = ANXA3,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/neutrophil activation/ANXA3.png", width = 5, height = 5, dpi = 300 )

ITGAM <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000169896)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())

ggsave(plot = ITGAM,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/neutrophil activation/ITGAM.png", width = 5, height = 5, dpi = 300 )


IL18RAP <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000115607)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())

ggsave(plot = IL18RAP,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/neutrophil activation/IL18RAP.png", width = 5, height = 5, dpi = 300 )

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

# Do bigger T cell plot
Tcell_pathways <- scores_mat_S1vsCCT0 %>% 
  dplyr::select("sample_id",
                "T cell receptor signaling pathway|GO:0050852",
                "T cell costimulation|GO:0031295",
                "positive regulation of T cell activation|GO:0050870",
                "T-helper 2 cell differentiation|GO:0045064",
                "T-helper 17 cell differentiation|GO:0072539",
                "T-helper 1 cell differentiation|GO:0045063")




Tcell_pathways <- merge(Tcell_pathways, samples, all.x =TRUE)


colnames(Tcell_pathways)
Tcell_pathways <- Tcell_pathways %>% gather(pathway, value, 2:7)
Tcell_pathways <- Tcell_pathways%>% tidyr::separate(pathway, into = c('Pathway', 'GO number'), sep = "\\|GO") %>% dplyr::select(-`GO number`) %>%  spread(Pathway, value)

colnames(Tcell_pathways)


Tcell_pathways <- Tcell_pathways %>% gather(Pathway, Z_score, 19:24)

library(broom)

Tcell <- Tcell_pathways %>% filter(State == "State 1"  |
                                       State == "State 2"|
                                       State == "State 3") %>% group_by(Pathway) %>% 
  do(augment(loess(Z_score~Lineage1, .))) %>% 
  {ggplot(., aes(x = Lineage1, y= Z_score, label = Pathway)) + 
      geom_vline(xintercept = 0.03, colour= "grey")+
      geom_vline(xintercept = 0.15, colour= "grey")+
      stat_smooth(aes(group = Pathway, linetype = Pathway, label = Pathway),method = "loess", se = FALSE, colour =  "darkblue", size =2)+
      stat_smooth(colour = "black", size = 2)+
      guides(color = F) +
      geom_label_repel(data = filter(., Lineage1 == max(Lineage1)),
                       aes(Lineage1, .fitted), size = 9, nudge_x = 0.06, max.overlaps = 20)+
      theme_classic()+
      theme(legend.position = "none")+
      theme(text = element_text(size = 24))+
      xlab("Pseudotime")}


ggsave(plot = Tcell,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/Tcell.png", width = 10, height = 7, dpi = 300 )

#plot by the PHATE plot
library(heatmaply)
colnames(Tcell_pathways)
phate_Tcell <- Tcell_pathways  %>% filter(State == "State 1"  |
                                                  State == "State 2"|
                                                  State == "State 3") %>% 
  spread(Pathway, Z_score) %>% 
  ggplot( aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = `T-helper 1 cell differentiation`), pch = 21, size = 12, colour = "black")+
  scale_fill_viridis_c()+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = 24))
ggsave(plot = phate_Tcell,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/phate/T helper 1.png", width =7, height = 7, dpi = 300 )

library(ggbeeswarm)
#plot the boxplot of each condition
T_boxpot <- Tcell_pathways %>% 
  spread(Pathway, Z_score) %>% 
  ggplot( aes(x = State, y= `T-helper 1 cell differentiation`))+
  geom_boxplot()+
  geom_quasirandom(aes(fill = State), pch = 21, size = 8, colour = "black")+
  scale_fill_manual( values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF",
                                "Cardiac_T0"= "#56A8CBFF", "Cardiac_T1" ="#DA291CFF"))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(text = element_text(size = 24))
ggsave(plot = T_boxpot,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/BP/T_boxpot_1.png", width = 6, height = 6, dpi = 300 )


Tcell_pathways$State <- factor(Tcell_pathways$State, levels = c("Cardiac_T0", "Cardiac_T1",
                                                                "State 1", "State 2", "State 3"))
TcellBP <- ggplot(Tcell_pathways %>% filter(Pathway == "T-helper 2 cell differentiation"|
                                   Pathway == "T-helper 1 cell differentiation"|
                                   Pathway == "T-helper 17 cell differentiation"), aes(x = State, y= Z_score, fill = Pathway))+
  geom_boxplot()+
  #geom_quasirandom(aes(fill = State), pch = 21, size = 8, colour = "black")+
  scale_fill_manual( values = c("T-helper 2 cell differentiation" = "#008080",
                                "T-helper 1 cell differentiation" = "#F7BBDB", 
                                "T-helper 17 cell differentiation"= "gold"))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(text = element_text(size = 24))
  

ggsave(plot = TcellBP,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/BP/TcellBP.png", width = 6, height = 6, dpi = 300 )
library(heatmaply)
GATA3 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000107485)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())

ggsave(plot = GATA3,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/Tcell/GATA3.png", width = 5, height = 5, dpi = 300 )

TBX21 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000073861)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())

ggsave(plot = TBX21,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/Tcell/TBX21.png", width = 5, height = 5, dpi = 300 )

RORC <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000143365)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())

ggsave(plot = RORC,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/Tcell/RORC.png", width = 5, height = 5, dpi = 300 )

JAK3 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000105639)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())

ggsave(plot = JAK3,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/Tcell/JAK3.png", width = 5, height = 5, dpi = 300 )
####################################################################################################################
####################################################################################################################
####################################################################################################################

# Do bigger B cel plot
bcell_pathways <- scores_mat_S1vsCCT0 %>% 
  dplyr::select("sample_id",
                "negative regulation of B cell activation|GO:0050869",
                "positive regulation of B cell proliferation|GO:0030890",
                "B cell differentiation|GO:0030183",
                "mature B cell differentiation involved in immune response|GO:0002313",
                "antimicrobial humoral response|GO:0019730")




bcell_pathways <- merge(bcell_pathways, samples, all.x =TRUE)


colnames(bcell_pathways)
bcell_pathways <- bcell_pathways %>% gather(pathway, value, 2:6)
bcell_pathways <- bcell_pathways%>% tidyr::separate(pathway, into = c('Pathway', 'GO number'), sep = "\\|GO") %>% dplyr::select(-`GO number`) %>%  spread(Pathway, value)

colnames(bcell_pathways)


bcell_pathways <- bcell_pathways %>% gather(Pathway, Z_score, 19:23)



bcell <- bcell_pathways %>% filter(State == "State 1"  |
                                     State == "State 2"|
                                     State == "State 3") %>% group_by(Pathway) %>% 
  do(augment(loess(Z_score~Lineage1, .))) %>% 
  {ggplot(., aes(x = Lineage1, y= Z_score, label = Pathway)) + 
      geom_vline(xintercept = 0.03, colour= "grey")+
      geom_vline(xintercept = 0.15, colour= "grey")+
      stat_smooth(aes(group = Pathway, linetype = Pathway, label = Pathway),method = "loess", se = FALSE, colour =  "darkgreen", size =2)+
      stat_smooth(colour = "black", size = 2)+
      guides(color = F) +
      geom_label_repel(data = filter(., Lineage1 == max(Lineage1)),
                       aes(Lineage1, .fitted), size = 9, nudge_x = 0.06, max.overlaps = 20)+
      theme_classic()+
      theme(legend.position = "none")+
      theme(text = element_text(size = 24))+
      xlab("Pseudotime")}


ggsave(plot = bcell,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/bcell.png", width = 10, height = 7, dpi = 300 )

#plot by the PHATE plot

colnames(bcell_pathways)
phate_bcell <- bcell_pathways  %>% filter(State == "State 1"  |
                                            State == "State 2"|
                                            State == "State 3") %>% 
  spread(Pathway, Z_score) %>% 
  ggplot( aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = `B cell differentiation`), pch = 21, size = 12, colour = "black")+
  scale_fill_viridis_c()+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = 24))
ggsave(plot = phate_bcell,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/phate/B cell differentiation.png", width =7, height = 7, dpi = 300 )

library(ggbeeswarm)
#plot the boxplot of each condition
b_boxpot <- bcell_pathways %>% 
  spread(Pathway, Z_score) %>% 
  ggplot( aes(x = State, y= `B cell differentiation`))+
  geom_boxplot()+
  geom_quasirandom(aes(fill = State), pch = 21, size = 8, colour = "black")+
  scale_fill_manual( values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF",
                                "Cardiac_T0"= "#56A8CBFF", "Cardiac_T1" ="#DA291CFF"))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(text = element_text(size = 24))

ggsave(plot = b_boxpot,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/BP/BcellBP.png", width = 6, height = 6, dpi = 300 )

CD27 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000139193)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())

ggsave(plot = CD27,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/Bcell/CD27.png", width = 5, height = 5, dpi = 300 )

CD40LG <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000102245)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())

ggsave(plot = CD40LG,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/Bcell/CD40LG.png", width = 5, height = 5, dpi = 300 )



DOCK10 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000135905)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())

ggsave(plot = DOCK10,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/Bcell/DOCK10.png", width = 5, height = 5, dpi = 300 )

BCL6 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000113916)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())

ggsave(plot = BCL6,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/Bcell/BCL6.png", width = 5, height = 5, dpi = 300 )
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

# Do bigger metabolism plot
metabolism_pathways <- scores_mat_S1vsCCT0 %>% 
  dplyr::select("sample_id",
                "nitric oxide metabolic process|GO:0046209",
                "cellular modified amino acid metabolic process|GO:0006575",
                "positive regulation of lipid metabolic process|GO:0045834",
                "energy reserve metabolic process|GO:0006112",
                "arginine metabolic process|GO:0006525",
                "cell division|GO:0051301",
                "positive regulation of mitotic cell cycle phase transition|GO:1901992")




metabolism_pathways <- merge(metabolism_pathways, samples, all.x =TRUE)

colnames(metabolism_pathways)
metabolism_pathways <- metabolism_pathways %>% gather(pathway, value, 2:8)
metabolism_pathways <- metabolism_pathways%>% tidyr::separate(pathway, into = c('Pathway', 'GO number'), sep = "\\|GO") %>% dplyr::select(-`GO number`) %>%  spread(Pathway, value)

colnames(metabolism_pathways)


metabolism_pathways <- metabolism_pathways %>% gather(Pathway, Z_score, 19:25)

library(broom)

metabolism <- metabolism_pathways %>% group_by(Pathway) %>% 
  do(augment(loess(Z_score~Lineage1, .))) %>% 
  {ggplot(., aes(x = Lineage1, y= Z_score, label = Pathway)) + 
      geom_vline(xintercept = 0.03, colour= "grey")+
      geom_vline(xintercept = 0.15, colour= "grey")+
      stat_smooth(aes(group = Pathway, linetype = Pathway, label = Pathway),method = "loess", se = FALSE, colour =  "salmon", size =2)+
      stat_smooth(colour = "black", size = 2)+
      guides(color = F) +
      geom_label_repel(data = filter(., Lineage1 == max(Lineage1)),
                       aes(Lineage1, .fitted), size = 9, nudge_x = 0.1,nudge_y = 1.3, max.overlaps = 20)+
      theme_classic()+
      theme(legend.position = "none")+
      theme(text = element_text(size = 24))+
      xlab("Pseudotime")}


ggsave(plot = metabolism,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/metabolism_CC.png", width = 10, height = 7, dpi = 300 )

#plot by the PHATE plot
library(heatmaply)
colnames(Tcell_pathways)
phate_met <- metabolism_pathways  %>% filter(State == "State 1"  |
                                            State == "State 2"|
                                            State == "State 3") %>% 
  spread(Pathway, Z_score) %>% 
  ggplot( aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = `cell division`), pch = 21, size = 12, colour = "black")+
  scale_fill_viridis_c()+
  theme_classic()+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = 24))
ggsave(plot = phate_met,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/phate/cell division.png", width =7, height = 7, dpi = 300 )

library(ggbeeswarm)
#plot the boxplot of each condition
met_boxpot <- metabolism_pathways %>% 
  spread(Pathway, Z_score) %>% 
  ggplot( aes(x = State, y= `cell division`))+
  geom_boxplot()+
  geom_quasirandom(aes(fill = State), pch = 21, size = 8, colour = "black")+
  scale_fill_manual( values = c("State 1" = "#FDE725FF", "State 2" = "#29AF7FFF", "State 3" = "#453781FF",
                                "Cardiac_T0"= "#56A8CBFF", "Cardiac_T1" ="#DA291CFF"))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(text = element_text(size = 24))

ggsave(plot = met_boxpot,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Lineplots/BP/celldivisionBP.png", width = 6, height = 6, dpi = 300 )

CDK1 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000170312)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())

ggsave(plot = CDK1,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/cell division/CDK1.png", width = 5, height = 5, dpi = 300 )

CDC20 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000117399)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())

ggsave(plot = CDC20,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/cell division/CDC20.png", width = 5, height = 5, dpi = 300 )

SH3GLB1 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000097033)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())

ggsave(plot = SH3GLB1,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/cell division/SH3GLB1.png", width = 5, height = 5, dpi = 300 )

map9 <- ggplot(matrix, aes(x = PHATE1, y= PHATE2))+
  geom_point(aes(fill = normalize(ENSG00000164114)), pch = 21, size = 12, colour = "black")+
  theme_minimal()+
  scale_fill_viridis_c(option = "A")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 26))+
  theme( axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.position="none",
         panel.background=element_blank(),
         panel.border=element_blank(),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         plot.background=element_blank())

ggsave(plot = map9,"FIGURES/IMMERSE paper/SCRIPT 4/genetonic/Genes/cell division/MAP9.png", width = 5, height = 5, dpi = 300 )



####################################################################################################
