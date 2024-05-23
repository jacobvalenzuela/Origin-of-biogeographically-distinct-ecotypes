rm(list=ls())
getwd()
setwd("/Users/jake/Documents/R_analysis/R_Analysis/ENIGMA_Analysis/FBR_DESeq/Final_Code/")
getwd()

Packages<-c("dplyr", "FactoMineR", "factoextra", "ggplot2", "graphics", "grDevices", "matrixStats", "pheatmap", 
            "RColorBrewer","readr", "reshape2", "stats", "SummarizedExperiment", "tidyr","textshape",
            "readr", "tidyverse", "corrplot", "ggrepel")
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
DvH_cnts<-data.frame(genes = DvH_gene_names_index$Gene_ID[ match(rownames(DvH_cnts), DvH_gene_names_index$EMBL_Annotation) ],DvH_cnts, check.names = FALSE)
rownames(DvH_cnts)<-DvH_cnts[,1]
DvH_cnts<-DvH_cnts[,-1]
DvH_cnts<-as.matrix(DvH_cnts)
is.numeric(DvH_cnts)

#convert EMBL ides to Gene names
Mmp_cnts<-data.frame(genes = Mmp_gene_names_index$Gene_ID[ match(rownames(Mmp_cnts), Mmp_gene_names_index$EMBL_Annotation) ],Mmp_cnts, check.names = FALSE)
rownames(Mmp_cnts)<-Mmp_cnts[,1]
Mmp_cnts<-Mmp_cnts[,-1]
Mmp_cnts<-as.matrix(Mmp_cnts)
is.numeric(Mmp_cnts)

colnames(DvH_cnts) == colnames(Mmp_cnts)

SynCom_cnts<-rbind(DvH_cnts, Mmp_cnts)
colnames(SynCom_cnts)<-gsub("EPD", "1K-EPD", colnames(SynCom_cnts))

###### Build WT Matrix
WT_SynCom_cnts<-merge(read.table("./Final_Data_Files//01-L1-WT-2.genes.results",
                              header = TRUE,  check.names = FALSE)[,c(1,6)],
                   read.table("./Final_Data_Files/02-L1-WT-4.genes.results",
                              header = TRUE,  check.names = FALSE)[,c(1,6)], 
                   by = "gene_id")

colnames(WT_SynCom_cnts)[2:3]<-c("WT_P_A_2017", "WT_P_B_2017")
head(WT_SynCom_cnts,8); tail(WT_SynCom_cnts,8)
WT_SynCom_cnts$gene_id<-gsub("DVU_", "DVU", WT_SynCom_cnts$gene_id)
WT_SynCom_cnts<-cbind(colsplit(WT_SynCom_cnts$gene_id, pattern = "_", names = c("gene_ID", "gene_symbol")), WT_SynCom_cnts)
head(WT_SynCom_cnts,8); tail(WT_SynCom_cnts,8)
WT_SynCom_cnts<-WT_SynCom_cnts[,-c(2,3)]

rownames(WT_SynCom_cnts)<-WT_SynCom_cnts[,1]
WT_SynCom_cnts<-WT_SynCom_cnts[,-1]
WT_SynCom_cnts<-as.matrix(WT_SynCom_cnts)
is.numeric(WT_SynCom_cnts)

nrow(SynCom_cnts); nrow(WT_SynCom_cnts)
head(SynCom_cnts,8); head(WT_SynCom_cnts,8); 
### Merge all data sets
TPMs_all<-merge(WT_SynCom_cnts, SynCom_cnts, by = "row.names", all.x = T)
rownames(TPMs_all)<-TPMs_all[,1]
TPMs_all<-TPMs_all[,-1]
TPMs_all<-as.matrix(TPMs_all)
nrow(TPMs_all)
TPMs_all<-as.matrix(na.omit(TPMs_all))
nrow(TPMs_all)

#### remove non-overlapping genes
rm(list=setdiff(ls(), c("TPMs_all", "DvH_gene_names_index", "Mmp_gene_names_index")))
colnames(TPMs_all)

TPMs_long<-melt(TPMs_all, varnames =c("gene_id", "Sample"), value.name = "TPMs")

theme_tpms<-theme(plot.title = element_text(angle = 0, size = 7),
                        axis.text = element_text(angle = 0, size = 6),
                        axis.title= element_text(angle = 0, size = 7),
                        strip.text = element_text(angle = 0, size = 5, color = "black", margin = margin(c(.07,.07,.07,.07), unit='cm')),
                        legend.box.background = element_rect(fill= "white", color = "grey", size =.25), 
                        legend.title = element_text(size=5),
                        legend.position = "none", legend.text = element_text(size = 5), # c(.09,.89)
                        legend.key.size = unit(.05, "in"), legend.key.width = unit(.2, "in"),
                        legend.justification = 'center', legend.title.align=0.5, legend.margin=margin(c(.05,.05,.05,.05), unit='cm'))
## check expression patterns
ggplot(data = TPMs_long)+
  geom_boxplot(aes(x = Sample, y = log10(TPMs)), outlier.size = .2, linewidth = .05)+
  theme_tpms+theme(axis.text.x = element_text(angle = 90, size = 6))+
  ggtitle("Raw Data (TPMs)")

### quantile normalization
#BiocManager::install("preprocessCore")
library("preprocessCore")
TPMs_all_quantnorm<-normalize.quantiles(TPMs_all,copy=TRUE, keep.names = T)
colnames(TPMs_all_quantnorm)<-colnames(TPMs_all)

TPMs_norm_long<-melt(TPMs_all_quantnorm, varnames =c("gene_id", "Sample"), value.name = "TPMs")

ggplot(data = TPMs_norm_long)+
  geom_boxplot(aes(x = Sample, y = log10(TPMs)), outlier.size = .2, linewidth = .05)+
  theme_tpms+theme(axis.text.x = element_text(angle = 90, size = 6))+
  ggtitle("Quantile Normalization (TPMs)")

########################
## create coldata ##
coldata<-data.frame("Sample_type" = colnames(TPMs_all_quantnorm), colsplit(colnames(TPMs_all_quantnorm), pattern = "_", names = c("Strain", "Phase", "Rep", "Date")))
rownames(coldata)<-coldata[,1]
#change date to day
coldata$Date<-gsub("2017","1", coldata$Date)
coldata$Date<-gsub("5-8","1", coldata$Date)
coldata$Date<-gsub("5-9","2", coldata$Date)
coldata$Date<-gsub("5-10","3", coldata$Date)
coldata$Date<-gsub("5-11","4", coldata$Date)
coldata$Date<-gsub("5-12","5", coldata$Date)
coldata$Date<-gsub("5-13","6", coldata$Date)
#add condition
coldata$condition<-paste(coldata$Strain, coldata$Phase, coldata$Date, sep = "_")
##
coldata[sapply(coldata, is.character)] <- lapply(coldata[sapply(coldata, is.character)], 
                                       as.factor)
# It is absolutely critical that the columns of the count matrix and the rows of the 
#column data (information about samples) are in the same order. If TRUE continue
all(rownames(coldata) %in% colnames(TPMs_all_quantnorm)) #are the rownames the same as the colnames
all(rownames(coldata) == colnames(TPMs_all_quantnorm)) #are the rownames in the same order as colnames

####Running DESeq for DvH all planktonic against all sediment
SynCom_dds_all <- DESeqDataSetFromMatrix(countData = round(TPMs_all_quantnorm),
                              colData = coldata,
                              design = ~ condition) ### for design between all P/S using Sample_type instead of the condition
                              #design = ~ condition) ####  for P/S for each day use condition
SynCom_dds_all

SynCom_dds_all <- estimateSizeFactors(SynCom_dds_all)
####

SynCom_vsd <- vst(SynCom_dds_all, blind = T)
pcaData <-plotPCA(SynCom_vsd, intgroup=c("Strain", "Phase", "Date"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData$Date<-ifelse(pcaData$Strain == "WT", yes = "WT", no = pcaData$Date)
pcaData$Date<-as.factor(pcaData$Date)
pcaData$Date<-factor(pcaData$Date, levels = c("WT", "1", "2", "3", "4", "5", "6"))

ggplot(pcaData, aes(PC1, PC2, fill=Phase, shape = Date, size=Strain, color = group)) +
  geom_point(stroke = .2) +
  scale_color_manual(values = c("WT:P:1" = "black", "1K-EPD:P:1" ="#1f78b4", "1K-EPD:P:2"  ="black", "1K-EPD:P:3"  ="black", "1K-EPD:P:4"  ="black", "1K-EPD:P:5"  ="black", "1K-EPD:P:6" = "black",
                                "1K-EPD:S:1" ="#a6cee3", "1K-EPD:S:2"  ="black", "1K-EPD:S:3"  ="black", "1K-EPD:S:4"  ="black", "1K-EPD:S:5"  ="black", "1K-EPD:S:6" = "black"), 
                     guide = NULL)+
  scale_fill_manual(values = c("P" =  "#1f78b4", "S" = "#a6cee3"), guide =NULL, name = "Phase")+
  scale_size_manual(values = c("WT" = 3.4, "1K-EPD" = 1.7))+
  scale_shape_manual(values = c("WT" = 13, "1" =  8, "2" =  21, "3" =  22, "4" =  23, "5" =  24, "6" = 25), guide =NULL, name = "Day")+
  
  #scale_shape_manual(values = c("P" = 21, "S" = 22))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme_bw()+theme_tpms+theme(legend.position = "right")+
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "black",size = 1.5, stroke = .2,
                                                 fill = c("P" =  "#1f78b4", "S" = "#a6cee3"))),
         shape =guide_legend(override.aes = list( colour = "black", size = 1.5, stroke = .2,
                                                 fill = c("WT" = "white", "1" =  "white", "2" =  "white", "3" =  "white", "4" =  "white", "5" =  "white", "6" =  "white"))))
ggsave("./Final_Plots/vsd_PCA_All_TPMs.pdf", width = 3.5, height = 3.75, dpi = "retina")
write.csv(assay(SynCom_vsd), "./Source_Data/Figure_3B.csv")

