
# IMMERSE paper 

Analysis script for creation of an trajectory sepsis immune defined by sepsis immune states (SIS) and validation in previously published datasets.

This is the first upload for reviewers to see and use. 

Files for analysis are stored in the folder called quant files. 

Load scripts for creation of analysis in paper. 

## Workflow - States analysis

Script 1 = Analysis of transcriptomic samples by clinical time points to identify the 718 genes used for analyses. Predicts previous published endotypes on our data. 

Script 2 = Clustering and PCA of data by clinical timepoints. 

Script 3 = Identification of sepsis trajectory and immune states

Script 4 = Statistics and genetonic (Uses GO pathways to conduct a gene set enrichment analysis) on sepsis immune states compared to pre-surgery and post-surgery

Script 5 = Using statistics in script 4 we correlate the log2fold change of genes in the comparisons of Post- vs pre-surgery compared to the states vs pre-surgery.

Script 6 = Plot the proportions of previously published endotypes in states.

Script 7 = Runs statistics and uses genetonic to identify pathways based upon the previously reported endotypes

Script 8 = Plotting of reactome analysis. 


## Workflow - states validation

Validating states in sepsis and non infectious insult (trauma)using four previously published data sets. 

1) Davenport 2016 (SIS)
2) Scicluna 2017 (MARS)
3) Sweeney 2019
4) Xioa 2011 (trauma)
