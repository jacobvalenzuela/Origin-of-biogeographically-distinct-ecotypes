Heatmap_Packages<-c("dplyr", "FactoMineR", "factoextra", "ggplot2", "graphics", "grDevices", "matrixStats", "pheatmap", 
            "RColorBrewer","readr", "reshape2", "stats", "SummarizedExperiment", "tidyr","textshape",
            "readr", "tidyverse", "corrplot", "ggpubr", "pals", "ComplexHeatmap", "circlize", "paletteer")
#install.packages(Heatmap_Packages)
lapply(Heatmap_Packages, library, character.only =T)

#####
'%notin%' <- Negate('%in%')
##### Building Dvh Mmp Heatmaps for transcriptional states (use complexHeatmap)
## Run DvH_Mmp_DESeq_Final.R to generate DvH and Mmp results
#create dataframe for all interactions for DvH
df1_DvH<-data.frame("GeneID" = DvH_DEG_DF_list$DvH_Results_Pairwise_Days_EPD.P.5.8_vs_EPD.S.5.8$GeneID, "Day1" = 1)
df2_DvH<-data.frame("GeneID" = DvH_DEG_DF_list$DvH_Results_Pairwise_Days_EPD.P.5.9_vs_EPD.S.5.9$GeneID, "Day2" = 1)
df3_DvH<-data.frame("GeneID" = DvH_DEG_DF_list$DvH_Results_Pairwise_Days_EPD.P.5.10_vs_EPD.S.5.10$GeneID, "Day3" = 1)
df4_DvH<-data.frame("GeneID" = DvH_DEG_DF_list$DvH_Results_Pairwise_Days_EPD.P.5.11_vs_EPD.S.5.11$GeneID, "Day4" = 1)
df5_DvH<-data.frame("GeneID" = DvH_DEG_DF_list$DvH_Results_Pairwise_Days_EPD.P.5.12_vs_EPD.S.5.12$GeneID, "Day5" = 1)
df6_DvH<-data.frame("GeneID" = DvH_DEG_DF_list$DvH_Results_Pairwise_Days_EPD.P.5.13_vs_EPD.S.5.13$GeneID, "Day6" = 1)
df7_DvH<-data.frame("GeneID" = DvH_DEG_DF_list$DvH_Results_Planktonic_vs_Sediment$GeneID, "All" = 1)

df_list_DvH <- list(df1_DvH, df2_DvH, df3_DvH, df4_DvH, df5_DvH, df6_DvH, df7_DvH)
DvH_DEG_Overlap<-Reduce(function(x, y) merge(x, y, by = "GeneID", all=TRUE), df_list_DvH)  
DvH_DEG_Overlap[is.na(DvH_DEG_Overlap)] <- 0   
DvH_DEG_Overlap<-tibble::column_to_rownames(DvH_DEG_Overlap, var = "GeneID")

#create dataframe for all interactions for Mmp
df1_Mmp<-data.frame("GeneID" = Mmp_DEG_DF_list$Mmp_Results_Pairwise_Days_EPD.P.5.8_vs_EPD.S.5.8$GeneID, "Day1" = 1)
df2_Mmp<-data.frame("GeneID" = Mmp_DEG_DF_list$Mmp_Results_Pairwise_Days_EPD.P.5.9_vs_EPD.S.5.9$GeneID, "Day2" = 1)
df3_Mmp<-data.frame("GeneID" = Mmp_DEG_DF_list$Mmp_Results_Pairwise_Days_EPD.P.5.10_vs_EPD.S.5.10$GeneID, "Day3" = 1)
df4_Mmp<-data.frame("GeneID" = Mmp_DEG_DF_list$Mmp_Results_Pairwise_Days_EPD.P.5.11_vs_EPD.S.5.11$GeneID, "Day4" = 1)
df5_Mmp<-data.frame("GeneID" = Mmp_DEG_DF_list$Mmp_Results_Pairwise_Days_EPD.P.5.12_vs_EPD.S.5.12$GeneID, "Day5" = 1)
df6_Mmp<-data.frame("GeneID" = Mmp_DEG_DF_list$Mmp_Results_Pairwise_Days_EPD.P.5.13_vs_EPD.S.5.13$GeneID, "Day6" = 1)
df7_Mmp<-data.frame("GeneID" = Mmp_DEG_DF_list$Mmp_Results_Planktonic_vs_Sediment$GeneID, "All" = 1)

df_list_Mmp <- list(df1_Mmp, df2_Mmp, df3_Mmp, df4_Mmp, df5_Mmp, df6_Mmp, df7_Mmp)
Mmp_DEG_Overlap<-Reduce(function(x, y) merge(x, y, by = "GeneID", all=TRUE), df_list_Mmp)  
Mmp_DEG_Overlap[is.na(Mmp_DEG_Overlap)] <- 0   
Mmp_DEG_Overlap<-tibble::column_to_rownames(Mmp_DEG_Overlap, var = "GeneID")

#extract z-score counts from DEG lists
DvH_Sig_mat_All<-(as.matrix(subset(DvH_z_counts, rownames(DvH_z_counts) %in% DvH_Results_Planktonic_vs_Sediment$GeneID))) #padj < 0.001 & abs(log2FoldChange) > 1
Mmp_Sig_mat_All<-(as.matrix(subset(Mmp_z_counts, rownames(Mmp_z_counts) %in% Mmp_Results_Planktonic_vs_Sediment$GeneID))) #padj < 0.001 & abs(log2FoldChange) > 1

rownames(DvH_Results_Planktonic_vs_Sediment)<-DvH_Results_Planktonic_vs_Sediment$GeneID; rownames(Mmp_Results_Planktonic_vs_Sediment)<-Mmp_Results_Planktonic_vs_Sediment$GeneID;
#DvH_Results_Planktonic_vs_Sediment_mat<-as.matrix(DvH_Results_Planktonic_vs_Sediment)
#Mmp_Results_Planktonic_vs_Sediment_mat<-as.matrix(Mmp_Results_Planktonic_vs_Sediment)

membrane_transport<-c("DVU2683","DVU2617","DVU1818", "DVU0095", "DVU3090", "DVU2373", "DVU0943","DVU1671", "DVU0251", "DVU2342","DVU1548")
Flagella<-c("DVU0517","DVU1441","DVU1443","DVU2444") 

### build an ordered dataframe of day 1-6 (split each day by positive and negative) DvH
DvH_Day_1<-DvH_Results_Pairwise_Days_EPD.P.5.8_vs_EPD.S.5.8[order(DvH_Results_Pairwise_Days_EPD.P.5.8_vs_EPD.S.5.8$log2FoldChange, decreasing = T),]
DvH_Day_2<-DvH_Results_Pairwise_Days_EPD.P.5.9_vs_EPD.S.5.9[order(DvH_Results_Pairwise_Days_EPD.P.5.9_vs_EPD.S.5.9$log2FoldChange, decreasing = T),]
DvH_Day_3<-DvH_Results_Pairwise_Days_EPD.P.5.10_vs_EPD.S.5.10[order(DvH_Results_Pairwise_Days_EPD.P.5.10_vs_EPD.S.5.10$log2FoldChange, decreasing = T),]
DvH_Day_4<-DvH_Results_Pairwise_Days_EPD.P.5.11_vs_EPD.S.5.11[order(DvH_Results_Pairwise_Days_EPD.P.5.11_vs_EPD.S.5.11$log2FoldChange, decreasing = T),]
DvH_Day_5<-DvH_Results_Pairwise_Days_EPD.P.5.12_vs_EPD.S.5.12[order(DvH_Results_Pairwise_Days_EPD.P.5.12_vs_EPD.S.5.12$log2FoldChange, decreasing = T),]
DvH_Day_6<-DvH_Results_Pairwise_Days_EPD.P.5.13_vs_EPD.S.5.13[order(DvH_Results_Pairwise_Days_EPD.P.5.13_vs_EPD.S.5.13$log2FoldChange, decreasing = T),]
DvH_Day_all<-DvH_Results_Planktonic_vs_Sediment[order(DvH_Results_Planktonic_vs_Sediment$log2FoldChange, decreasing = T),]

DvH_Day_1$Day<-"Day1"; DvH_Day_2$Day<-"Day2"; DvH_Day_3$Day<-"Day3"; DvH_Day_4$Day<-"Day4";DvH_Day_5$Day<-"Day5"; DvH_Day_6$Day<-"Day6"; DvH_Day_all$Day<-"All"

DvH_Ord_DEG<-rbind(subset(DvH_Day_1, DvH_Day_1$log2FoldChange > 0), subset(DvH_Day_2, DvH_Day_2$log2FoldChange > 0), subset(DvH_Day_3, DvH_Day_3$log2FoldChange > 0), 
                   subset(DvH_Day_4, DvH_Day_4$log2FoldChange > 0), subset(DvH_Day_5, DvH_Day_5$log2FoldChange > 0), subset(DvH_Day_6, DvH_Day_6$log2FoldChange > 0),
                   subset(DvH_Day_all, DvH_Day_all$log2FoldChange > 0),
                   subset(DvH_Day_1, DvH_Day_1$log2FoldChange < 0), subset(DvH_Day_2, DvH_Day_2$log2FoldChange < 0), subset(DvH_Day_3, DvH_Day_3$log2FoldChange < 0), 
                   subset(DvH_Day_4, DvH_Day_4$log2FoldChange < 0), subset(DvH_Day_5, DvH_Day_5$log2FoldChange < 0), subset(DvH_Day_6, DvH_Day_6$log2FoldChange < 0),
                   subset(DvH_Day_all, DvH_Day_all$log2FoldChange < 0))

### generate a unique dataframe and matrix
DvH_Ord_DEG_IDs<-unique(DvH_Ord_DEG$GeneID)
DvH_Ord_DEG_unique <- DvH_Ord_DEG[!duplicated(DvH_Ord_DEG$GeneID),]
rownames(DvH_Ord_DEG_unique)<-DvH_Ord_DEG_unique$GeneID

DvH_Sig_mat_Ord<-as.data.frame(subset(DvH_z_counts, rownames(DvH_z_counts) %in% DvH_Ord_DEG_IDs))
DvH_Sig_mat_Ord<-rownames_to_column(DvH_Sig_mat_Ord, var = "rowname")

### order unique matrix the same as the order of concatenated dataframe in sequential order (day 1 up, day 2 up, etc.)
merge(DvH_Sig_mat_Ord, data.frame("GeneID"= DvH_Ord_DEG_IDs), by.x = "rowname", by.y = "GeneID", sort=FALSE)
DvH_Sig_mat_Ordered<-merge(data.frame("GeneID"= DvH_Ord_DEG_IDs), DvH_Sig_mat_Ord,  by.y = "rowname", by.x = "GeneID", sort=FALSE)
DvH_Sig_mat_Ordered<-tibble::column_to_rownames(DvH_Sig_mat_Ordered, var = "GeneID")

#Split data on breaks between up and down
DvH_splits<-data.frame("GeneID" =rownames(DvH_Sig_mat_Ordered), "Group" = c(rep("DvH UpReg (P)", 172), rep("DvH UpReg (S)", nrow(DvH_Sig_mat_Ordered)-172)))
DvH_splits<-tibble::column_to_rownames(DvH_splits, var = "GeneID")

### build an ordered dataframe of day 1-6 (split each day by positive and negative) Mmp
Mmp_Day_1<-Mmp_Results_Pairwise_Days_EPD.P.5.8_vs_EPD.S.5.8[order(Mmp_Results_Pairwise_Days_EPD.P.5.8_vs_EPD.S.5.8$log2FoldChange, decreasing = T),]
Mmp_Day_2<-Mmp_Results_Pairwise_Days_EPD.P.5.9_vs_EPD.S.5.9[order(Mmp_Results_Pairwise_Days_EPD.P.5.9_vs_EPD.S.5.9$log2FoldChange, decreasing = T),]
Mmp_Day_3<-Mmp_Results_Pairwise_Days_EPD.P.5.10_vs_EPD.S.5.10[order(Mmp_Results_Pairwise_Days_EPD.P.5.10_vs_EPD.S.5.10$log2FoldChange, decreasing = T),]
Mmp_Day_4<-Mmp_Results_Pairwise_Days_EPD.P.5.11_vs_EPD.S.5.11[order(Mmp_Results_Pairwise_Days_EPD.P.5.11_vs_EPD.S.5.11$log2FoldChange, decreasing = T),]
Mmp_Day_5<-Mmp_Results_Pairwise_Days_EPD.P.5.12_vs_EPD.S.5.12[order(Mmp_Results_Pairwise_Days_EPD.P.5.12_vs_EPD.S.5.12$log2FoldChange, decreasing = T),]
Mmp_Day_6<-Mmp_Results_Pairwise_Days_EPD.P.5.13_vs_EPD.S.5.13[order(Mmp_Results_Pairwise_Days_EPD.P.5.13_vs_EPD.S.5.13$log2FoldChange, decreasing = T),]
Mmp_Day_all<-Mmp_Results_Planktonic_vs_Sediment[order(Mmp_Results_Planktonic_vs_Sediment$log2FoldChange, decreasing = T),]

Mmp_Day_1$Day<-"Day1"; Mmp_Day_2$Day<-"Day2"; Mmp_Day_3$Day<-"Day3"; Mmp_Day_4$Day<-"Day4";Mmp_Day_5$Day<-"Day5"; Mmp_Day_6$Day<-"Day6"; Mmp_Day_all$Day<-"All"

Mmp_Ord_DEG<-rbind(subset(Mmp_Day_1, Mmp_Day_1$log2FoldChange > 0), subset(Mmp_Day_2, Mmp_Day_2$log2FoldChange > 0), subset(Mmp_Day_3, Mmp_Day_3$log2FoldChange > 0), 
                   subset(Mmp_Day_4, Mmp_Day_4$log2FoldChange > 0), subset(Mmp_Day_5, Mmp_Day_5$log2FoldChange > 0), subset(Mmp_Day_6, Mmp_Day_6$log2FoldChange > 0),
                   subset(Mmp_Day_all, Mmp_Day_all$log2FoldChange > 0),
                   subset(Mmp_Day_1, Mmp_Day_1$log2FoldChange < 0), subset(Mmp_Day_2, Mmp_Day_2$log2FoldChange < 0), subset(Mmp_Day_3, Mmp_Day_3$log2FoldChange < 0), 
                   subset(Mmp_Day_4, Mmp_Day_4$log2FoldChange < 0), subset(Mmp_Day_5, Mmp_Day_5$log2FoldChange < 0), subset(Mmp_Day_6, Mmp_Day_6$log2FoldChange < 0),
                   subset(Mmp_Day_all, Mmp_Day_all$log2FoldChange < 0))

### generate a unique dataframe and matrix
Mmp_Ord_DEG_IDs<-unique(Mmp_Ord_DEG$GeneID)
Mmp_Ord_DEG_unique <- Mmp_Ord_DEG[!duplicated(Mmp_Ord_DEG$GeneID),]
rownames(Mmp_Ord_DEG_unique)<-Mmp_Ord_DEG_unique$GeneID

Mmp_Sig_mat_Ord<-as.data.frame(subset(Mmp_z_counts, rownames(Mmp_z_counts) %in% Mmp_Ord_DEG_IDs))
Mmp_Sig_mat_Ord<-rownames_to_column(Mmp_Sig_mat_Ord, var = "rowname")

### order unique matrix the same as the order of concatenated dataframe in sequential order (day 1 up, day 2 up, etc.)
merge(Mmp_Sig_mat_Ord, data.frame("GeneID"= Mmp_Ord_DEG_IDs), by.x = "rowname", by.y = "GeneID", sort=FALSE)
Mmp_Sig_mat_Ordered<-merge(data.frame("GeneID"= Mmp_Ord_DEG_IDs), Mmp_Sig_mat_Ord,  by.y = "rowname", by.x = "GeneID", sort=FALSE)
Mmp_Sig_mat_Ordered<-tibble::column_to_rownames(Mmp_Sig_mat_Ordered, var = "GeneID")

#Split data on breaks between up and down
Mmp_splits<-data.frame("GeneID" = rownames(Mmp_Sig_mat_Ordered), "Group" = c(rep("Mmp UpReg (P)", 150), rep("Mmp UpReg (S)", nrow(Mmp_Sig_mat_Ordered)-150)))
Mmp_splits<-tibble::column_to_rownames(Mmp_splits, var = "GeneID")

#
#DvH_Sig_mat_Ordered$rowid <- seq(1, nrow(DvH_Sig_mat_Ordered)*2, by = 2) # or simply 1:nrow(x)
#Mmp_Sig_mat_Ordered$rowid <- seq(2, nrow(Mmp_Sig_mat_Ordered)*2, by = 2)


cbind(DvH_Sig_mat_Ordered, "order" = seq(1, nrow(DvH_Sig_mat_Ordered), by =1))
SynCom_Ordered<-rbind(cbind(DvH_Sig_mat_Ordered, "order" = seq(1, nrow(DvH_Sig_mat_Ordered), by =1)), 
                      cbind(Mmp_Sig_mat_Ordered, "order" = seq(nrow(DvH_Sig_mat_Ordered)+1, nrow(DvH_Sig_mat_Ordered)+nrow(Mmp_Sig_mat_Ordered), by =1)))
SynCom_Ordered<-tibble::rownames_to_column(SynCom_Ordered, var = "GeneID")
SynCom_Ordered_index<-subset(gene_names_index, gene_names_index$Gene_ID %in% SynCom_Ordered$GeneID)
SynCom_Ordered_index<-merge(SynCom_Ordered, SynCom_Ordered_index, by.x ="GeneID", by.y ="Gene_ID")
SynCom_Ordered_index<-SynCom_Ordered_index[order(SynCom_Ordered_index$order), ]

write.csv(SynCom_Ordered_index[,c(14,1,17)], "./Final_Output/DEG/SynCom/SynCom_Ordered_index.csv")

#creating lists of GO enrichments
match(c("DVU1922","DVU1921","DVU1769","DVU2329"),rownames(DvH_Sig_mat_Ordered))

SynCom_Ordered_index[,c(14,1,17)]

### build metadata data frame
LacAceMeth<-read.csv(file = "./Final_Data_Files/FBR_metadata.csv", nrows = 45, check.names = F)
colnames(LacAceMeth)
summs_meta<-LacAceMeth %>%
  dplyr::group_by(Sample_type, Days) %>%
  summarise(dplyr::across(c(OD_Balch, pH, DO, mM_Lactate, mM_Acetate, mM_Methane), #ng_DNA_per_ml, ng_DNA_per_gram, ng_RNA_per_ml, ng_RNA_per_gram, ug_mL_Protein, ng_Protein_per_mg
                   list("mean" = mean, "sd" = sd, "n" = length))) 
summs_meta$Days2<-summs_meta$Days-1

lac_df<-rbind(summs_meta[3:8,"mM_Lactate_mean"],summs_meta[3:8,"mM_Lactate_mean"])
ace_df<-rbind(summs_meta[3:8,"mM_Acetate_mean"],summs_meta[3:8,"mM_Acetate_mean"])
ch4_df<-rbind(summs_meta[3:8,"mM_Methane_mean"],summs_meta[3:8,"mM_Methane_mean"])

#building heatmaps
ht_opt(RESET = T)                                                                            
ht_opt(heatmap_row_title_gp = gpar(fontsize = 6),
       legend_border = "black",
       heatmap_border = TRUE,
       annotation_border = TRUE)  

### color schemes
col_rnorm = colorRamp2(c(-3, 0, 3), c("#1b7837","white","#762a83"))
col_padj = colorRamp2(c(0, 40), c("white","magenta"))
col_basemean = colorRamp2(c(0,2, 5), c("white","lightgrey", "#3D1778"))
col_FC = colorRamp2(c(-1.5,0,.5), c("#3D1778","white", "Darkgreen"))
col_base = colorRamp2(c(1, 4), c("grey","#800000"))
col_mM_Lac = colorRamp2(c(0, 35), c("white","#006666"))
col_mM_Ace = colorRamp2(c(0, 35), c("white","#009999"))
col_mM_CH4 = colorRamp2(c(0, 1), c("white","#662700"))
col_Zscore_DEG = colorRamp2(c(-1, 0, 1), c("Blue", "white", "Red"), space = "RGB", )
#col_Zscore_FBA = colorRamp2(c(-1.5, 0, 1.5), c("Blue", "white", "Red"), space = "RGB") Original
col_Zscore_FBA = colorRamp2(c(-1.5, 0, 1.5), c("Blue", "white", "Red"), space = "RGB") ## Revised
col_FBA_Mat = colorRamp2(c(-1000, 0, 1000), c("Blue", "white", "Red"), space = "RGB")


ha1 = HeatmapAnnotation(Phase = c("P","P","P","P","P","P","S","S","S","S","S","S"), col = list(Phase = c("P" =  "#1f78b4", "S" = "#a6cee3")),
                        annotation_legend_param = list(
                          Phase = list(title = "Phase",
                                       at = c("P", "S"),
                                       labels = c("P" ="Planktonic", "S"="Sediment"))),
                        gp = gpar(col = "white"),
                        simple_anno_size = unit(.2, "cm"),
                        annotation_name_side = "left"
                        #show_annotation_name = c("Phase" = FALSE)
)

ha_meta = HeatmapAnnotation("mM Lactate" = lac_df$mM_Lactate_mean,
                           "mM Acetate" = ace_df$mM_Acetate_mean,
                           "mM Methane" = ch4_df$mM_Methane_mean,
                           col = list("mM Lactate" = col_mM_Lac,
                                      "mM Acetate" = col_mM_Ace,
                                      "mM Methane" = col_mM_CH4),
                           show_legend = c(TRUE,TRUE,TRUE),
                           annotation_legend_param = list(
                             "mM Lactate" = list(title = "mM Lactate"),
                             "mM Acetate" = list(title = "mM Acetate"),
                             "mM Methane" = list(title = "mM Methane")),
                        gp = gpar(col = "white"),
                        annotation_name_gp = gpar(fontsize =6),
                        simple_anno_size = unit(.2, "cm")
                        #annotation_name_side = "left"
                        #show_annotation_name = c("Phase" = FALSE)
)

#### build a dataframe for specific labels
gene_ID_labels<-read.csv("./Final_Data_Files/SynCom_Ordered_index_processed.csv")
#Hydrogenases<-subset(gene_ID_labels, gene_ID_labels$Note %in% c("Hydrogenase", "Hydrogenase(LA)", "Lactate", "Methanogenesis", "TF", "Biofilm"))
Hydrogenases<-subset(gene_ID_labels, gene_ID_labels$Note %in% c("Flagella", "Hydrogenase", "Hydrogenase(LA)", "Lactate"))
Hydrogenases[which(Hydrogenases$Org == "Mmp"), "GeneID"]

####### 
ht5 = Heatmap(DvH_Sig_mat_Ordered, name = "Z-score",  row_names_gp = gpar(fontsize = 0), col = col_Zscore_DEG,
              top_annotation = ha_meta,
              cluster_rows = F, cluster_columns = F, clustering_distance_rows = "binary",
              column_names_gp = gpar(fontsize = 6), width = unit(2, "in"),  
              row_split = DvH_splits, 
              right_annotation = rowAnnotation("log10(BaseMean)" = log10(DvH_Ord_DEG_unique[,"baseMean"]), 
                                                                                #gp = gpar(fill = "black", lwd =0), axis_param = list(side = "top"), width = unit(.2, "in")),
                                               show_annotation_name = c("log10(BaseMean)" = FALSE), annotation_name_rot = c(90, 0),
                                               annotation_name_gp = gpar(fontsize = 6),col = list("log10(BaseMean)" = col_basemean),
                                               
                                               "GoEnrich" = anno_mark(at = match(Hydrogenases[which(Hydrogenases$Org == "DvH"), "GeneID"],rownames(DvH_Sig_mat_Ordered)), 
                                                                      labels= Hydrogenases[which(Hydrogenases$Org == "DvH"), "Gene_name_edited"], labels_gp = gpar(fontsize =3.7)
                                               )
                                               
                                               
              )
)
ht5

ht6 = Heatmap(Mmp_Sig_mat_Ordered, name = "Z-score",  row_names_gp = gpar(fontsize = 0), col = col_Zscore_DEG,
              cluster_rows = F, cluster_columns = F, clustering_distance_rows = "binary",
              bottom_annotation = ha1, column_labels = rep(c("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6"), 2), 
              column_names_gp = gpar(fontsize = 6), width = unit(2, "in"),  
              row_split = Mmp_splits,
              right_annotation = rowAnnotation("log10(BaseMean)" = log10(Mmp_Ord_DEG_unique[,"baseMean"]), 
                                                                                #gp = gpar(fill = "black", lwd = 0 ), axis_param = list(side = "bottom"), width = unit(.305, "in")),
                                               show_annotation_name = c("log10(BaseMean)" = TRUE), annotation_name_rot = c(90, 0),
                                               annotation_name_gp = gpar(fontsize = 4),col = list("log10(BaseMean)" = col_basemean),
                                               "GoEnrich" = anno_mark(at = match(Hydrogenases[which(Hydrogenases$Org == "Mmp"), "GeneID"],rownames(Mmp_Sig_mat_Ordered)),
                                                                      labels= Hydrogenases[which(Hydrogenases$Org == "Mmp"), "Gene_name_edited"], labels_gp = gpar(fontsize =3.7)
                                               )
                                               
              )
)
ht6

pdf("./Final_Plots/SynCom_DEG_Heatmap_TimeCourse_RGB.pdf", height = 6.7)
ht_list3 = ht5 %v% ht6
draw(ht_list3, column_title_gp = gpar(fontsize = 0), row_title_gp = gpar(fontsize = 3), merge_legend = TRUE, column_km =2, 
     heatmap_legend_side = "left")
dev.off()

write.csv(rbind(DvH_Sig_mat_Ordered, Mmp_Sig_mat_Ordered), "./Source_Data/Figure_3A.csv")

          