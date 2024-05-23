rm(list=ls())
getwd()
setwd("/Users/jake/Documents/R_analysis/R_Analysis/ENIGMA_Analysis/FBR_DESeq/Final_Code/")
getwd()

Packages<-c("dplyr", "FactoMineR", "factoextra", "ggplot2", "graphics", "grDevices", "matrixStats", "pheatmap", 
            "RColorBrewer","readr", "reshape2", "stats", "SummarizedExperiment", "tidyr","textshape",
            "readr", "tidyverse", "corrplot")
#install.packages(Packages)
lapply(Packages, library, character.only =T)

#BiocManager packages
BiocManager_Packages<-c("Biostrings", "DESeq2")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install(BiocManager_Packages)
lapply(BiocManager_Packages, library, character.only =T)

##############################################################################################
#### DVH Info from microbes online
DvH_info<-read.delim(file = "./Final_Data_Files/DvHgenomeInfo.txt", header = T)
Mmp_info<-read.delim(file = "./Final_Data_Files/MmpgenomeInfo.txt", header = T)
#### importing gene info ####
gene_info<-readDNAStringSet("./Final_Data_Files/dvh_mmp_merged.cdna.all.fa") ### Merged.FA file
## split data and return relevant info
gene_names = data.frame(ID=names(gene_info))
gene_names_split<-colsplit(gene_names$ID, " ", names=c("EMBL_Annotation", "Source", "Chromosome", "Gene_ID", 
                                                       "Gene_Biotype","Transcript_Biotype", "Gene_Symbol+Description"))
####preparing gene name info
rm(gene_names_split2)
gene_names_split2<-cbind(gene_names_split,colsplit(gene_names_split$Chromosome, ":", names = c("a","chromosome","b", "loci")))
gene_names_split2<-gene_names_split2[,-c(3,5,6,8,10)]
gene_names_split2$Gene_ID<-gsub("gene:","", gene_names_split2$Gene_ID)
gene_names_split2$Gene_ID<-gsub("[_]","", gene_names_split2$Gene_ID, ignore.case = T)
gene_names_index<-gene_names_split2
gene_names_index$`Gene_Symbol+Description`<-gsub("description:", "", gene_names_index$`Gene_Symbol+Description`)
gene_names_index$`Gene_Symbol+Description`<-gsub("gene_symbol:", "", gene_names_index$`Gene_Symbol+Description`)
gene_names_index$`Gene_Symbol+Description`<-gsub("^$", NA, gene_names_index$`Gene_Symbol+Description`)
gene_names_index$`Gene_Symbol+Description`<-ifelse(is.na(gene_names_index$`Gene_Symbol+Description`), gene_names_index$Gene_ID, gene_names_index$`Gene_Symbol+Description`)
DvH_gene_names_index<-gene_names_index[1:3531,]
Mmp_gene_names_index<-gene_names_index[3532:nrow(gene_names_index),]

### Count matrix import for DvH ###############################################################################
DvH_cnts<-as.matrix(read.table("./Final_Data_Files/dv_kallisto_TPM_norm_matrix.txt",
                          header = TRUE,  check.names = FALSE))
head(DvH_cnts,2) #sediment A_5-12 was not sequenced
Mmp_cnts<-as.matrix(read.table("./Final_Data_Files/mm_kallisto_TPM_norm_matrix.txt",
                           header = TRUE,  check.names = FALSE))
head(Mmp_cnts,2)
#convert EMBL ides to Gene names
DvH_cts<-data.frame(genes = DvH_gene_names_index$Gene_ID[ match(rownames(DvH_cnts), DvH_gene_names_index$EMBL_Annotation) ],DvH_cnts, check.names = FALSE)
rownames(DvH_cts)<-DvH_cts[,1]
DvH_cts<-DvH_cts[,-1]
DvH_cts<-as.matrix(DvH_cts)
is.numeric(DvH_cts)

#convert EMBL ides to Gene names
Mmp_cts<-data.frame(genes = Mmp_gene_names_index$Gene_ID[ match(rownames(Mmp_cnts), Mmp_gene_names_index$EMBL_Annotation) ],Mmp_cnts, check.names = FALSE)
rownames(Mmp_cts)<-Mmp_cts[,1]
Mmp_cts<-Mmp_cts[,-1]
Mmp_cts<-as.matrix(Mmp_cts)
is.numeric(Mmp_cts)

## check if cnts and cts is true
table(DvH_cts==DvH_cnts); table(Mmp_cts==Mmp_cnts)

### Import metadata file ###
coldata <- read.csv(file = "./Final_Data_Files/merged-coldata.csv", row.names = 1,header=TRUE)
coldata<-coldata %>% tibble::column_to_rownames(var ="Sample_Name")

### Add Day Column
coldata$Day<-coldata$Day-7

### Add new column data as factor
coldata$Sample_type_Day <- factor(paste(coldata$Sample_type, coldata$Day,sep = "_"))
##
coldata$EPD <- factor(coldata$EPD)
coldata$condition <- factor(coldata$condition)
coldata$Type <- factor(coldata$Type)
coldata$Day <- factor(coldata$Day)
coldata$Sample_type <- factor(coldata$Sample_type)

## From DESeq2
# It is absolutely critical that the columns of the count matrix and the rows of the 
#column data (information about samples) are in the same order. If TRUE continue
all(rownames(coldata) %in% colnames(DvH_cts)); all(rownames(coldata) %in% colnames(Mmp_cts)) #are the rownames the same as the colnames
all(rownames(coldata) == colnames(DvH_cts)); all(rownames(coldata) == colnames(Mmp_cts))  #are the rownames in the same order as colnames
DvH_cts <- DvH_cts[, rownames(coldata)]; Mmp_cts <- Mmp_cts[, rownames(coldata)] #organize the colnames of cts to be in order of rownames of the metadata
all(rownames(coldata) == colnames(DvH_cts)); all(rownames(coldata) == colnames(Mmp_cts)) #are the rownames the same as the colnames, MUST BE TRUE

####Running DESeq for DvH all planktonic against all sediment
DvH_dds_all <- DESeqDataSetFromMatrix(countData = round(DvH_cts),
                              colData = coldata,
                              design = ~ Sample_type) ### for design between all P/S using Sample_type instead of the condition
                              #design = ~ condition) ####  for P/S for each day use condition
DvH_dds_all

DvH_dds_all <- collapseReplicates(DvH_dds_all, groupby=DvH_dds_all$condition) ### groups by condition removing Replicate information
####
#  pre-filter low count genes before running the DESeq2 functions
Dvh_keep_all <- rowSums(counts(DvH_dds_all)) >= 100
DvH_dds_all <- DvH_dds_all[Dvh_keep_all,]
# leveling factors (not done for DVH Time Series) but this is how to set the reference
DvH_dds_all$Sample_type <- relevel(DvH_dds_all$Sample_type, ref = "Sediment") ### Planktonic/Sediment
DvH_dds_all$Sample_type <- droplevels(DvH_dds_all$Sample_type)

DvH_dds_all<-DESeq(DvH_dds_all)

DvH_res_all<-as.data.frame(results(DvH_dds_all))
DvH_res_all<-cbind("GeneID" = rownames(DvH_res_all),DvH_res_all)
write.csv(as.data.frame(DvH_res_all), "./Final_Output/DEG/DvH//DvH_Results_Planktonic_vs_Sediment.csv", row.names = F)
#order our results table by the smallest p value:
DvH_resOrdered_all <- as.data.frame(DvH_res_all[order(DvH_res_all$baseMean, decreasing = TRUE),])
#write.csv(as.data.frame(resOrdered), file="Results/condition_HC_vs_LC.csv")
DvH_resSig_all <- as.data.frame(subset(DvH_resOrdered_all, padj < 0.01 & abs(log2FoldChange) > 1))

DvH_vsd <- vst(DvH_dds_all, blind=FALSE)
DvH_rld <- rlog(DvH_dds_all, blind = FALSE)
DvH_norm_counts<-counts(DvH_dds_all, normalized=T)

##reorder norm_counts matrix
DvH_norm_counts<-DvH_norm_counts[,c("EPD.P.5.8","EPD.P.5.9","EPD.P.5.10","EPD.P.5.11","EPD.P.5.12","EPD.P.5.13",
                                    "EPD.S.5.8","EPD.S.5.9","EPD.S.5.10","EPD.S.5.11","EPD.S.5.12","EPD.S.5.13")]

## make metadata dataframe to match norm_counts matrix
rm(Sample_metadata)
Sample_metadata<-cbind("Sample_Label" = colnames(DvH_norm_counts), colsplit(colnames(DvH_norm_counts), "\\.", names=c("EPD", "Phase", "Month", "Date")))
Sample_metadata<-Sample_metadata %>% tibble::column_to_rownames(var ="Sample_Label")
Sample_metadata$Day<-Sample_metadata$Date
Sample_metadata$Day<-gsub("8","1", Sample_metadata$Day)
Sample_metadata$Day<-gsub("9","2", Sample_metadata$Day)
Sample_metadata$Day<-gsub("10","3", Sample_metadata$Day)
Sample_metadata$Day<-gsub("11","4", Sample_metadata$Day)
Sample_metadata$Day<-gsub("12","5", Sample_metadata$Day)
Sample_metadata$Day<-gsub("13","6", Sample_metadata$Day)
Sample_metadata[sapply(Sample_metadata, is.character)]<-lapply(Sample_metadata[sapply(Sample_metadata,is.character)], as.factor)

colnames(DvH_norm_counts)==rownames(Sample_metadata)

####Build Long dataframes for Normalized counts and z-scores
DvH_long_df_counts<-melt(DvH_norm_counts, value.name = "norm_counts", varnames = c("Gene_ID", "Sample_Label"))
DvH_long_df_counts<-cbind(DvH_long_df_counts, colsplit(DvH_long_df_counts$Sample_Label, "\\.", names=c("EPD", "Phase", "Month", "Date")))
DvH_long_df_counts$Day<-DvH_long_df_counts$Date
DvH_long_df_counts$Day<-gsub("8","1", DvH_long_df_counts$Day) #### can be loopified
DvH_long_df_counts$Day<-gsub("9","2", DvH_long_df_counts$Day)
DvH_long_df_counts$Day<-gsub("10","3", DvH_long_df_counts$Day)
DvH_long_df_counts$Day<-gsub("11","4", DvH_long_df_counts$Day)
DvH_long_df_counts$Day<-gsub("12","5", DvH_long_df_counts$Day)
DvH_long_df_counts$Day<-gsub("13","6", DvH_long_df_counts$Day)

### Creating z-scores for plots https://www.r-bloggers.com/2012/03/r-tutorial-series-centering-variables-and-generating-z-scores-with-the-scale-function/
DvH_z_counts <- t(scale(t(DvH_norm_counts),center = TRUE, scale = TRUE)) #
DvH_long_df_zscore<-melt(DvH_z_counts, value.name = "zScore", varnames = c("Gene_ID", "Sample_Label"))
DvH_long_df_zscore<-cbind(DvH_long_df_zscore, colsplit(DvH_long_df_zscore$Sample_Label, "\\.", names=c("EPD", "Phase", "Month", "Date")))
DvH_long_df_zscore$Day<-DvH_long_df_zscore$Date
DvH_long_df_zscore$Day<-gsub("8","1", DvH_long_df_zscore$Day) #### can be loopified
DvH_long_df_zscore$Day<-gsub("9","2", DvH_long_df_zscore$Day)
DvH_long_df_zscore$Day<-gsub("10","3", DvH_long_df_zscore$Day)
DvH_long_df_zscore$Day<-gsub("11","4", DvH_long_df_zscore$Day)
DvH_long_df_zscore$Day<-gsub("12","5", DvH_long_df_zscore$Day)
DvH_long_df_zscore$Day<-gsub("13","6", DvH_long_df_zscore$Day)

####################################################################
####Running DESeq for DvH and each Days
DvH_dds_days <- DESeqDataSetFromMatrix(countData = round(DvH_cts),
                                      colData = coldata,
                                      design = ~ condition) ####  for P/S for each day use condition
DvH_dds_days

####
#  pre-filter low count genes before running the DESeq2 functions
DvH_keep_days <- rowSums(counts(DvH_dds_days)) >= 100
DvH_dds_days <- DvH_dds_days[DvH_keep_days,]

DvH_dds_days<-DESeq(DvH_dds_days)

###create consecutive differentially expressed genes over the days
resultsNames(DvH_dds_days)
#using contrast to compare different sets of conditions: example:  results(dds_Acer, contrast = c("Treatment_TP", "Control_T4", "Control_T1"))
comparisons<-read.csv(file = "./Final_Data_Files/DEG_comparisons.csv", header = T, check.names = F)

results(DvH_dds_days, contrast = c("condition", "EPD.P.5.8", "EPD.S.5.8")) ## Check on DEseq design

for (i in rownames(comparisons)){
  print(i)
  DEG_comp<-comparisons[i,"Comparison"]
  DEG_design<-comparisons[i,"Design"]
  ref<-comparisons[i,"Reference"]
  cond<-comparisons[i,"Condition"]
  print(DEG_comp)
  print(ref)
  print(cond)
  results_comp<-results(DvH_dds_days, contrast = c(DEG_design, cond, ref))
  results_comp$GeneID=row.names(results_comp)
  results_comp<-results_comp[,c(7,1:6)]
  write.csv(data.frame(results_comp), 
            file = paste("./Final_Output/DEG/DvH/DvH_Results_", DEG_comp, ".csv", sep = ""), row.names = F)
}

#### read files for DEG analysis
setwd("Final_Output/DEG/DvH")
DvH_Results_list<-gsub("\\.csv$", "", list.files(pattern ="*.csv", full.names=F))
DvH_DEG_list<-gsub("\\.csv$", "", list.files(pattern ="*.csv", full.names=F))

for (i in DvH_Results_list){
  assign(i, as.data.frame(read.csv(paste(i, ".csv", sep=""))))
}

for (i in DvH_DEG_list){
  assign(i, as.data.frame(subset(read.csv(paste(i, ".csv", sep="")), padj < 0.01 & abs(log2FoldChange) > 1)))
}

setwd("/Users/jake/Documents/R_analysis/R_Analysis/ENIGMA_Analysis/FBR_DESeq/Final_Code/")
getwd()

#### Create list of dataframes for all combinations
DvH_DEG_DF_list<-setNames(lapply(DvH_DEG_list, get), DvH_DEG_list)

###number of genes
DvH_DEGs<-t(data.frame("DvH_Results_Pairwise_Days_EPD.P.5.8_vs_EPD.S.5.8" = nrow(DvH_DEG_DF_list$DvH_Results_Pairwise_Days_EPD.P.5.8_vs_EPD.S.5.8),
                      "DvH_Results_Pairwise_Days_EPD.P.5.9_vs_EPD.S.5.9" = nrow(DvH_DEG_DF_list$DvH_Results_Pairwise_Days_EPD.P.5.9_vs_EPD.S.5.9),
                      "DvH_Results_Pairwise_Days_EPD.P.5.10_vs_EPD.S.5.10" = nrow(DvH_DEG_DF_list$DvH_Results_Pairwise_Days_EPD.P.5.10_vs_EPD.S.5.10),
                      "DvH_Results_Pairwise_Days_EPD.P.5.11_vs_EPD.S.5.11" = nrow(DvH_DEG_DF_list$DvH_Results_Pairwise_Days_EPD.P.5.11_vs_EPD.S.5.11),
                      "DvH_Results_Pairwise_Days_EPD.P.5.12_vs_EPD.S.5.12" = nrow(DvH_DEG_DF_list$DvH_Results_Pairwise_Days_EPD.P.5.12_vs_EPD.S.5.12),
                      "DvH_Results_Pairwise_Days_EPD.P.5.13_vs_EPD.S.5.13" = nrow(DvH_DEG_DF_list$DvH_Results_Pairwise_Days_EPD.P.5.13_vs_EPD.S.5.13),
                      "DvH_Results_Planktonic_vs_Sediment" = nrow(DvH_DEG_DF_list$DvH_Results_Planktonic_vs_Sediment)
))
                      
colnames(DvH_DEGs)[1]<-"Sig_DEG"
DvH_DEGs
### creating DEG list for DvH
DvH_long_DEG<-NULL ### create empty DF
for (i in DvH_Results_list){
  print(i)
  temp1<-as.data.frame(DvH_DEG_DF_list[i], col.names = NULL)
  temp1$Design<-i
  temp1$Org<-"DvH"
  DvH_long_DEG<-rbind(DvH_long_DEG,temp1)
}

DvH_long_DEG$Day<-DvH_long_DEG$Design
DvH_long_DEG$Day<-gsub("DvH_Results_Pairwise_Days_EPD.P.5.8_vs_EPD.S.5.8","1", DvH_long_DEG$Day)
DvH_long_DEG$Day<-gsub("DvH_Results_Pairwise_Days_EPD.P.5.9_vs_EPD.S.5.9","2", DvH_long_DEG$Day)
DvH_long_DEG$Day<-gsub("DvH_Results_Pairwise_Days_EPD.P.5.10_vs_EPD.S.5.10","3", DvH_long_DEG$Day)
DvH_long_DEG$Day<-gsub("DvH_Results_Pairwise_Days_EPD.P.5.11_vs_EPD.S.5.11","4", DvH_long_DEG$Day)
DvH_long_DEG$Day<-gsub("DvH_Results_Pairwise_Days_EPD.P.5.12_vs_EPD.S.5.12","5", DvH_long_DEG$Day)
DvH_long_DEG$Day<-gsub("DvH_Results_Pairwise_Days_EPD.P.5.13_vs_EPD.S.5.13","6", DvH_long_DEG$Day)
DvH_long_DEG$Day<-gsub("DvH_Results_Planktonic_vs_Sediment","All", DvH_long_DEG$Day)#### can be loopified

#########################
#######DESeq for Mmp##########
########################

####Running DESeq for Mmp all planktonic against all sediment
Mmp_dds_all <- DESeqDataSetFromMatrix(countData = round(Mmp_cts),
                                      colData = coldata,
                                      design = ~ Sample_type) ### for design between all P/S using Sample_type instead of the condition
#design = ~ condition) ####  for P/S for each day use condition
Mmp_dds_all

Mmp_dds_all <- collapseReplicates(Mmp_dds_all, groupby=Mmp_dds_all$condition) ### groups by condition removing Replicate information
####
#  pre-filter low count genes before running the DESeq2 functions
Mmp_keep_all <- rowSums(counts(Mmp_dds_all)) >= 100
Mmp_dds_all <- Mmp_dds_all[Mmp_keep_all,]
# leveling factors (not done for Mmp Time Series) but this is how to set the reference
Mmp_dds_all$Sample_type <- relevel(Mmp_dds_all$Sample_type, ref = "Sediment") ### Planktonic/Sediment
Mmp_dds_all$Sample_type <- droplevels(Mmp_dds_all$Sample_type)

Mmp_dds_all<-DESeq(Mmp_dds_all)

Mmp_res_all<-as.data.frame(results(Mmp_dds_all))
Mmp_res_all<-cbind("GeneID" = rownames(Mmp_res_all),Mmp_res_all)
write.csv(as.data.frame(Mmp_res_all), "./Final_Output/DEG/Mmp/Mmp_Results_Planktonic_vs_Sediment.csv", row.names = F)
#order our results table by the smallest p value:
Mmp_resOrdered_all <- as.data.frame(Mmp_res_all[order(Mmp_res_all$baseMean, decreasing = TRUE),])
#write.csv(as.data.frame(resOrdered), file="Results/condition_HC_vs_LC.csv")
Mmp_resSig_all <- as.data.frame(subset(Mmp_resOrdered_all, padj < 0.01 & abs(log2FoldChange) > 1))

Mmp_vsd <- vst(Mmp_dds_all, blind=FALSE)
Mmp_rld <- rlog(Mmp_dds_all, blind = FALSE)
Mmp_norm_counts<-counts(Mmp_dds_all, normalized=T)

##reorder norm_counts matrix
Mmp_norm_counts<-Mmp_norm_counts[,c("EPD.P.5.8","EPD.P.5.9","EPD.P.5.10","EPD.P.5.11","EPD.P.5.12","EPD.P.5.13",
                                    "EPD.S.5.8","EPD.S.5.9","EPD.S.5.10","EPD.S.5.11","EPD.S.5.12","EPD.S.5.13")]

## make metadata dataframe to match norm_counts matrix
rm(Sample_metadata)
Sample_metadata<-cbind("Sample_Label" = colnames(Mmp_norm_counts), colsplit(colnames(Mmp_norm_counts), "\\.", names=c("EPD", "Phase", "Month", "Date")))
Sample_metadata<-Sample_metadata %>% tibble::column_to_rownames(var ="Sample_Label")
Sample_metadata$Day<-Sample_metadata$Date
Sample_metadata$Day<-gsub("8","1", Sample_metadata$Day)
Sample_metadata$Day<-gsub("9","2", Sample_metadata$Day)
Sample_metadata$Day<-gsub("10","3", Sample_metadata$Day)
Sample_metadata$Day<-gsub("11","4", Sample_metadata$Day)
Sample_metadata$Day<-gsub("12","5", Sample_metadata$Day)
Sample_metadata$Day<-gsub("13","6", Sample_metadata$Day)
Sample_metadata[sapply(Sample_metadata, is.character)]<-lapply(Sample_metadata[sapply(Sample_metadata,is.character)], as.factor)

colnames(Mmp_norm_counts)==rownames(Sample_metadata)

####Build Long dataframes for Normalized counts and z-scores
Mmp_long_df_counts<-melt(Mmp_norm_counts, value.name = "norm_counts", varnames = c("Gene_ID", "Sample_Label"))
Mmp_long_df_counts<-cbind(Mmp_long_df_counts, colsplit(Mmp_long_df_counts$Sample_Label, "\\.", names=c("EPD", "Phase", "Month", "Date")))
Mmp_long_df_counts$Day<-Mmp_long_df_counts$Date
Mmp_long_df_counts$Day<-gsub("8","1", Mmp_long_df_counts$Day) #### can be loopified
Mmp_long_df_counts$Day<-gsub("9","2", Mmp_long_df_counts$Day)
Mmp_long_df_counts$Day<-gsub("10","3", Mmp_long_df_counts$Day)
Mmp_long_df_counts$Day<-gsub("11","4", Mmp_long_df_counts$Day)
Mmp_long_df_counts$Day<-gsub("12","5", Mmp_long_df_counts$Day)
Mmp_long_df_counts$Day<-gsub("13","6", Mmp_long_df_counts$Day)

### Creating z-scores for plots https://www.r-bloggers.com/2012/03/r-tutorial-series-centering-variables-and-generating-z-scores-with-the-scale-function/
Mmp_z_counts <- t(scale(t(Mmp_norm_counts),center = TRUE, scale = TRUE)) #
Mmp_long_df_zscore<-melt(Mmp_z_counts, value.name = "zScore", varnames = c("Gene_ID", "Sample_Label"))
Mmp_long_df_zscore<-cbind(Mmp_long_df_zscore, colsplit(Mmp_long_df_zscore$Sample_Label, "\\.", names=c("EPD", "Phase", "Month", "Date")))
Mmp_long_df_zscore$Day<-Mmp_long_df_zscore$Date
Mmp_long_df_zscore$Day<-gsub("8","1", Mmp_long_df_zscore$Day) #### can be loopified
Mmp_long_df_zscore$Day<-gsub("9","2", Mmp_long_df_zscore$Day)
Mmp_long_df_zscore$Day<-gsub("10","3", Mmp_long_df_zscore$Day)
Mmp_long_df_zscore$Day<-gsub("11","4", Mmp_long_df_zscore$Day)
Mmp_long_df_zscore$Day<-gsub("12","5", Mmp_long_df_zscore$Day)
Mmp_long_df_zscore$Day<-gsub("13","6", Mmp_long_df_zscore$Day)

####################################################################
####Running DESeq for Mmp and each Days
Mmp_dds_days <- DESeqDataSetFromMatrix(countData = round(Mmp_cts),
                                       colData = coldata,
                                       design = ~ condition) ####  for P/S for each day use condition
Mmp_dds_days

####
#  pre-filter low count genes before running the DESeq2 functions
Mmp_keep_days <- rowSums(counts(Mmp_dds_days)) >= 100
Mmp_dds_days <- Mmp_dds_days[Mmp_keep_days,]

Mmp_dds_days<-DESeq(Mmp_dds_days)

###create consecutive differentially expressed genes over the days
resultsNames(Mmp_dds_days)
#using contrast to compare different sets of conditions: example:  results(dds_Acer, contrast = c("Treatment_TP", "Control_T4", "Control_T1"))
#comparisons<-read.csv(file = "DEG_comparisons_DvH_Mmp.csv", header = T, check.names = F)

results(Mmp_dds_days, contrast = c("condition", "EPD.P.5.8", "EPD.S.5.8")) ## Check on DEseq design

for (i in rownames(comparisons)){
  print(i)
  DEG_comp<-comparisons[i,"Comparison"]
  DEG_design<-comparisons[i,"Design"]
  ref<-comparisons[i,"Reference"]
  cond<-comparisons[i,"Condition"]
  print(DEG_comp)
  print(ref)
  print(cond)
  results_comp<-results(Mmp_dds_days, contrast = c(DEG_design, cond, ref))
  results_comp$GeneID=row.names(results_comp)
  results_comp<-results_comp[,c(7,1:6)]
  write.csv(data.frame(results_comp), 
            file = paste("./Final_Output/DEG/Mmp/Mmp_Results_", DEG_comp, ".csv", sep = ""), row.names = F)
}

#spotcheck differential expression with norm_counts (double check with added samples)
norm_counts_Mmp_test<-as.data.frame(counts(Mmp_dds_days, normalized=T))
results_test<-as.data.frame(results(Mmp_dds_days, contrast = c("condition", "EPD.P.5.8", "EPD.S.5.8")))
results_test<-subset(results_test, padj < 0.001 & abs(log2FoldChange) > 1)
results_test2<-as.data.frame(read.csv("./Final_Output/DEG/Mmp/Mmp_Results_Pairwise_Days_EPD.P.5.8_vs_EPD.S.5.8.csv", header = TRUE))
results_test2<-subset(results_test2, padj < 0.001 & abs(log2FoldChange) > 1)
round(results_test2$pvalue, 3) == round(results_test$pvalue, 3)
rowSums(norm_counts_Mmp_test["MMP0369",])/length(norm_counts_Mmp_test)
results_test["MMP0369", "baseMean"]
M_T4<-rowSums(norm_counts_Mmp_test["MMP0369",c("EPD_S_A_5-8", "EPD_S_B_5-8", "EPD_S_C_5-8")])/3 ### Sediment
M_T5<-rowSums(norm_counts_Mmp_test["MMP0369",c("EPD_P_A_5-8", "EPD_P_B_5-8", "EPD_P_C_5-8")])/3 ##Planktonic

log2(M_T5/M_T4)
results_test["MMP0369", "log2FoldChange"]

#### read files for DEG analysis
setwd("./Final_Output/DEG/Mmp")
Mmp_Results_list<-gsub("\\.csv$", "", list.files(pattern ="*.csv", full.names=F))
Mmp_DEG_list<-gsub("\\.csv$", "", list.files(pattern ="*.csv", full.names=F))

for (i in Mmp_Results_list){
  assign(i, as.data.frame(read.csv(paste(i, ".csv", sep=""))))
}

for (i in Mmp_DEG_list){
  assign(i, as.data.frame(subset(read.csv(paste(i, ".csv", sep="")), padj < 0.01 & abs(log2FoldChange) > 1)))
}

setwd("/Users/jake/Documents/R_analysis/R_Analysis/ENIGMA_Analysis/FBR_DESeq/Final_Code/")
getwd()

#### Create list of dataframes for all combinations
Mmp_DEG_DF_list<-setNames(lapply(Mmp_DEG_list, get), Mmp_DEG_list)

###number of genes
Mmp_DEGs<-t(data.frame("Mmp_Results_Pairwise_Days_EPD.P.5.8_vs_EPD.S.5.8" = nrow(Mmp_DEG_DF_list$Mmp_Results_Pairwise_Days_EPD.P.5.8_vs_EPD.S.5.8),
                       "Mmp_Results_Pairwise_Days_EPD.P.5.9_vs_EPD.S.5.9" = nrow(Mmp_DEG_DF_list$Mmp_Results_Pairwise_Days_EPD.P.5.9_vs_EPD.S.5.9),
                       "Mmp_Results_Pairwise_Days_EPD.P.5.10_vs_EPD.S.5.10" = nrow(Mmp_DEG_DF_list$Mmp_Results_Pairwise_Days_EPD.P.5.10_vs_EPD.S.5.10),
                       "Mmp_Results_Pairwise_Days_EPD.P.5.11_vs_EPD.S.5.11" = nrow(Mmp_DEG_DF_list$Mmp_Results_Pairwise_Days_EPD.P.5.11_vs_EPD.S.5.11),
                       "Mmp_Results_Pairwise_Days_EPD.P.5.12_vs_EPD.S.5.12" = nrow(Mmp_DEG_DF_list$Mmp_Results_Pairwise_Days_EPD.P.5.12_vs_EPD.S.5.12),
                       "Mmp_Results_Pairwise_Days_EPD.P.5.13_vs_EPD.S.5.13" = nrow(Mmp_DEG_DF_list$Mmp_Results_Pairwise_Days_EPD.P.5.13_vs_EPD.S.5.13),
                       "Mmp_Results_Planktonic_vs_Sediment" = nrow(Mmp_DEG_DF_list$Mmp_Results_Planktonic_vs_Sediment)
))

colnames(Mmp_DEGs)[1]<-"Sig_DEG"
Mmp_DEGs
### creating DEG list for Mmp
Mmp_long_DEG<-NULL ### create empty DF
for (i in Mmp_Results_list){
  print(i)
  
  temp1<-as.data.frame(Mmp_DEG_DF_list[i], col.names = NULL)
  temp1$Design<-i
  temp1$Org<-"Mmp"
  Mmp_long_DEG<-rbind(Mmp_long_DEG,temp1)
}

Mmp_long_DEG$Day<-Mmp_long_DEG$Design
Mmp_long_DEG$Day<-gsub("Mmp_Results_Pairwise_Days_EPD.P.5.8_vs_EPD.S.5.8","1", Mmp_long_DEG$Day)
Mmp_long_DEG$Day<-gsub("Mmp_Results_Pairwise_Days_EPD.P.5.9_vs_EPD.S.5.9","2", Mmp_long_DEG$Day)
Mmp_long_DEG$Day<-gsub("Mmp_Results_Pairwise_Days_EPD.P.5.10_vs_EPD.S.5.10","3", Mmp_long_DEG$Day)
Mmp_long_DEG$Day<-gsub("Mmp_Results_Pairwise_Days_EPD.P.5.11_vs_EPD.S.5.11","4", Mmp_long_DEG$Day)
Mmp_long_DEG$Day<-gsub("Mmp_Results_Pairwise_Days_EPD.P.5.12_vs_EPD.S.5.12","5", Mmp_long_DEG$Day)
Mmp_long_DEG$Day<-gsub("Mmp_Results_Pairwise_Days_EPD.P.5.13_vs_EPD.S.5.13","6", Mmp_long_DEG$Day)
Mmp_long_DEG$Day<-gsub("Mmp_Results_Planktonic_vs_Sediment","All", Mmp_long_DEG$Day)#### can be loopified

DvH_Mmp_all_DEGs<-rbind(DvH_long_DEG,Mmp_long_DEG)

