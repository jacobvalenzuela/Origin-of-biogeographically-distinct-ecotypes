rm(list=ls())
getwd()
setwd("/Users/jake/Documents/R_analysis/R_Analysis/ENIGMA_Analysis/FBR_DESeq/Final_Code/")
getwd()

FBA_Packages<-c("dplyr", "FactoMineR", "factoextra", "ggplot2", "graphics", "grDevices", "matrixStats", "pheatmap", 
            "RColorBrewer","readr", "reshape2", "stats", "SummarizedExperiment", "tidyr","textshape",
             "readr", "tidyverse", "corrplot", "ggpubr", "data.table", "ggExtra", "ggrepel", "ggpmisc", "ggbreak", "colorRamp2")
#install.packages(FBA_Packages)
lapply(FBA_Packages, library, character.only =T)
'%notin%' <- Negate('%in%')
#####
DvH_info<-read.delim(file = "./Final_Data_Files/DvHgenomeInfo.txt", header = T)
###GPR 
gpr<-read.csv(file = "./Final_Data_Files/SynComModelGenes-GPRdelineated-Rxn-Pathways.csv", header = T) ##  GPR

### Differential reactions
FBA_methane<-read.csv(file = "./Final_Data_Files/FluxTable_methaneAsObj.csv", 
                      header = T, check.names = F, skip = 1) ### New Revised methane Spreadsheet
                      
colnames(FBA_methane)[c(1:2,5,6)]<-c("Org", "Rxn", "Rxn_name", "Subsystem_Pathway") ### rename columns
### order rows by column alphabetically
FBA_methane<-FBA_methane[order(FBA_methane$Metabolic_Overview),]

FBA_methane_raw<-FBA_methane
rownames(FBA_methane_raw)<-FBA_methane_raw[,"Rxn"] 
FBA_methane_raw[sapply(FBA_methane_raw, is.character)]<-lapply(FBA_methane_raw[sapply(FBA_methane_raw,is.character)], as.factor) ## convert chr col to factors
FBA_methane_mat_raw<-FBA_methane_raw[,c(13:ncol(FBA_methane_raw))] #### Create the raw matrix for rxn fluxes
### remove spaces from colnames
colnames(FBA_methane_mat_raw)<-gsub(" ", "_", colnames(FBA_methane_mat_raw))
#### removed non finite rows
is.finite(as.matrix(FBA_methane_mat_raw)) #### if all true continue, if some false remove non-finite rows
#FBA_methane_mat_raw <- FBA_methane_mat_raw[is.finite(rowSums(FBA_methane_mat_raw)),] ### did not change

### calculate row.variance across samples
rowVariance<-data.frame("rowVar" = apply(FBA_methane_mat_raw, 1, var),
                        "rowMax" = apply(FBA_methane_mat_raw, 1, max), 
                        "rowMin" = apply(FBA_methane_mat_raw, 1, min),
                        "rowMaxAbs" = apply(abs(FBA_methane_mat_raw), 1, max),
                        "rowMean" = apply(FBA_methane_mat_raw, 1, mean),
                        "rowSD" = apply(FBA_methane_mat_raw, 1, sd),
                        "rowMedian" = apply(FBA_methane_mat_raw, 1, median))

TwoSTDEV_Rxns<- subset(rowVariance, rowVariance$rowSD < 2) #95% of values will be within 2 standard deviations of the mean
var_threshold<-(max(TwoSTDEV_Rxns$rowVar))

Theme_variance<-theme(plot.title = element_text(angle = 0, size = 8, face = "plain", hjust = .5, vjust = -2),
                    axis.text = element_text(angle = 0, size = 6),
                    axis.title= element_text(angle = 0, size = 6.8),
                    strip.text = element_text(angle = 0, size = 5, color = "black", margin = margin(c(.07,.07,.07,.07), unit='cm')),
                    legend.box.background = element_rect(fill= "white", color = "grey", size =.25), 
                    legend.title = element_text(size=5),
                    legend.position = "none", legend.text = element_text(size = 5), # c(.09,.89)
                    legend.key.size = unit(.02, "in"), legend.key.width = unit(.02, "in"),
                    legend.justification = 'center', legend.title.align=0.5, legend.margin=margin(c(.05,.05,.05,.05), unit='cm'))
### Plotting Reaction Flux Variability
p<- ggplot(data = rowVariance, aes(log10(rowVar), rowSD)) +
  geom_point()+
  ylab("Reaction Standard Deviation")+ xlab ("log10(Reaction Flux Variance)")+
  geom_vline(xintercept = log10(max(TwoSTDEV_Rxns$rowVar)), linetype = 2)+
  annotate('text', x = -8, y = 12, label = paste("n =", nrow(rowVariance)))+
  theme_bw()+Theme_variance
p

# Bin width customization
p2<-ggMarginal(p, type = "histogram", 
               bins = 100, margins = "x")
p2

p3<- ggplot(TwoSTDEV_Rxns, 
            aes(log10(rowVar), rowSD)) +
  ylab("Reaction Standard Deviation")+ xlab ("log10(Reaction Flux Variance)")+
  geom_point()+ylim(0,2)+
  annotate('text', x = -8, y = 2, label = paste("n =", nrow(TwoSTDEV_Rxns)))+
  scale_x_continuous(breaks = seq(-10,5,5))+
  theme_bw()+Theme_variance
p3

# Bin width customization
p4<-ggMarginal(p3, type = "histogram", 
               bins = 100, margins = "x")
p4
rxnVar<-ggarrange(p2, p4,  ncol = 2)
rxnVar
#####
#####################################
## Analyzing core reaction fluxes

#remove rows with no variation
FBA_methane_mat_raw_VarHigh<-FBA_methane_mat_raw[apply(FBA_methane_mat_raw, 1, var) > var_threshold,]
FBA_methane_mat<-FBA_methane_mat_raw[apply(FBA_methane_mat_raw, 1, var) < var_threshold,]

#remove rows with no variation
FBA_methane_mat_noVar<-FBA_methane_mat[apply(FBA_methane_mat, 1, var) == 0, ] ##### rows with zero variation 

###subsetting rows with negative and positive values (flux flipping)
FBA_methane_mat_NegPos<-subset(FBA_methane_mat,(rowSums(sign(FBA_methane_mat)<0)>0) & (rowSums(sign(FBA_methane_mat)>0)>0))

FBA_methane_mat<-subset(FBA_methane_mat, rownames(FBA_methane_mat) %notin% rownames(FBA_methane_mat_NegPos))
### Calculate fold change
FBA_FC_Meth<-(data.frame(Day1=(FBA_methane_mat[,"Day_1_p"]/FBA_methane_mat[,"Day_1_s"]),
                        Day2=(FBA_methane_mat[,"Day_2_p"]/FBA_methane_mat[,"Day_2_s"]),
                        Day3=(FBA_methane_mat[,"Day_3_p"]/FBA_methane_mat[,"Day_3_s"]),
                        Day4=(FBA_methane_mat[,"Day_4_p"]/FBA_methane_mat[,"Day_4_s"]),
                        Day5=(FBA_methane_mat[,"Day_5_p"]/FBA_methane_mat[,"Day_5_s"]),
                        Day6=(FBA_methane_mat[,"Day_6_p"]/FBA_methane_mat[,"Day_6_s"])
                             ))
#### Positive is a flux through Planktonic and Negative is a Flux through Sediment
rownames(FBA_FC_Meth)<-rownames(FBA_methane_mat)
#### removed non finite rows
FBA_FC_Meth_NAN <- rbind(FBA_FC_Meth[is.infinite(rowSums(FBA_FC_Meth)),], FBA_FC_Meth[is.nan(rowSums(FBA_FC_Meth)),]) ### remove nan numbers
FBA_FC_Meth_INF <- FBA_FC_Meth[is.infinite(rowSums(FBA_FC_Meth)),] ### remove nan numbers
FBA_FC_Meth2 <- FBA_FC_Meth[is.finite(rowSums(FBA_FC_Meth)),] ### remove infinite numbers

#FBA_FC_Meth2<- FBA_FC_Meth2 %>% mutate(na_if,"") ### replace blanks with NA (may error just bypass)
#FBA_FC_Meth2<- na_if(FBA_FC_Meth2,"") ### replace blanks with NA (may error)
### round small values
FBA_FC_Meth2<-round(FBA_FC_Meth2, digits = 5) ### round values

#remove rows with no variation
FBA_FC_Meth3_noVar<-FBA_FC_Meth2[apply(FBA_FC_Meth2, 1, var) == 0, ]
FBA_FC_Meth3<-FBA_FC_Meth2[apply(FBA_FC_Meth2, 1, var) != 0, ]

summary(log(FBA_FC_Meth3,2) > 0)
#remove rows in which fold change is less than 1
###duplicated rows
FBA_FC_Meth_dups<-FBA_FC_Meth3[duplicated(FBA_FC_Meth3), ] ### duplicated rows
FBA_FC_Meth4<-FBA_FC_Meth3[!duplicated(FBA_FC_Meth3), ] ### non-duplicated rows

##
FBA_Summary_Table<-rbind(data.frame("Flux_Matrices" = "FBA_methane_mat_raw", "no. rxns" = nrow(FBA_methane_mat_raw)),
                 data.frame("Flux_Matrices" = "FBA_methane_mat", "no. rxns" = nrow(FBA_methane_mat)),
                 data.frame("Flux_Matrices" = "FBA_methane_mat_raw_VarHigh", "no. rxns" = nrow(FBA_methane_mat_raw_VarHigh)),
                 
                 data.frame("Flux_Matrices" = "FBA_methane_mat_noVar", "no. rxns" = nrow(FBA_methane_mat_noVar)),
                 data.frame("Flux_Matrices" = "FBA_methane_mat_NegPos", "no. rxns" = nrow(FBA_methane_mat_NegPos)),
                 data.frame("Flux_Matrices" = "FBA_FC_Meth_NAN", "no. rxns" = nrow(FBA_FC_Meth_NAN)),
                 data.frame("Flux_Matrices" = "FBA_FC_Meth_INF", "no. rxns" = nrow(FBA_FC_Meth_INF)),
                 data.frame("Flux_Matrices" = "FBA_FC_Meth2", "no. rxns" = nrow(FBA_FC_Meth2)),
                 data.frame("Flux_Matrices" = "FBA_FC_Meth3_noVar", "no. rxns" = nrow(FBA_FC_Meth3_noVar)),
                 data.frame("Flux_Matrices" = "FBA_FC_Meth_dups", "no. rxns" = nrow(FBA_FC_Meth_dups)),
                 data.frame("Flux_Matrices" = "FBA_FC_Meth3_Final", "no. rxns" = nrow(FBA_FC_Meth3)))
FBA_Summary_Table

###############################################
###Building a tnse
#BiocManager::install("Rtsne")
library(Rtsne)
head(FBA_methane_raw)
legend.position<-c(.1,.68)
theme_tsne<-theme(plot.title = element_text(angle = 0, size = 8, face = "plain", hjust = .5, vjust = -2),
                  axis.text = element_text(angle = 0, size = 6),
                  axis.title= element_text(angle = 0, size = 6.8),
                  strip.text = element_text(angle = 0, size = 5, color = "black", margin = margin(c(.07,.07,.07,.07), unit='cm')),
                  legend.box.background = element_rect(fill= "white", color = "grey", size =.25), 
                  legend.title = element_text(size=5),
                  legend.position = legend.position, legend.text = element_text(size = 5), # c(.09,.89)
                  legend.key.size = unit(.02, "in"), legend.key.width = unit(.02, "in"),
                  legend.justification = 'center', legend.title.align=0.5, legend.margin=margin(c(.05,.05,.05,.05), unit='cm'))

tsne_df<-subset(FBA_methane_mat, rownames(FBA_methane_mat) %in% rownames(FBA_FC_Meth3)) #### use matrix of core fluxes (FBA_FC_Meth3)
# Set a seed if you want reproducible results
set.seed(1)
tsne_out <-Rtsne::Rtsne(t(tsne_df),
                        perplex =3, check_duplicates = T, theta = 0.1, pca_center = F, pca_scale = T) #RUN tSNE
### build output dataframe for ggplot
tsne_out_df<-cbind(tsne_out$Y,colnames(tsne_df))
colnames(tsne_out_df)<-c("xdim", "ydim", "sample")
tsne_out_df<-as.data.frame(cbind(tsne_out_df, colsplit(tsne_out_df[,"sample"], pattern = "_", names = c("Day", "Days", "Phase"))))
tsne_out_df$xdim<-as.numeric(tsne_out_df$xdim)
tsne_out_df$ydim<-as.numeric(tsne_out_df$ydim)
head(tsne_out_df)

tsne_plt<-ggplot(data =tsne_out_df, aes(x= xdim, y= ydim))+
  geom_point(aes(color = factor(sample), fill = factor(Phase), shape = factor(Days)), size = 2.5)+
  scale_color_manual(values = c("Day_1_p" =  "#1f78b4", "Day_1_s" = "#a6cee3",
                                "Day_2_p" =  "black", "Day_2_s" = "black",
                                "Day_3_p" =  "black", "Day_3_s" = "black",
                                "Day_4_p" =  "black", "Day_4_s" = "black",
                                "Day_5_p" =  "black", "Day_5_s" = "black",
                                "Day_6_p" =  "black", "Day_6_s" = "black"), guide =NULL, name = "Phase")+
  scale_fill_manual(values = c("p" =  "#1f78b4", "s" = "#a6cee3"), guide =NULL, name = "Phase")+
  scale_shape_manual(values = c(8,21:25), guide =NULL, name = "Day")+
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "black",size = 1.5, stroke = .2,
                                                 fill = c("P" =  "#1f78b4", "S" = "#a6cee3"))),
         shape =guide_legend(override.aes = list(shape = c(8,21:25), colour = "black", size = 1.5, stroke = .2,
                                                 fill = c("1" =  "white", "2" =  "white", "3" =  "white", "4" =  "white", "5" =  "white", "6" =  "white"))))+
  scale_x_continuous(expand = c(.1,.1), name = "tSNE dimension 1")+
  scale_y_continuous(expand = c(.1,.1), name = "tSNE dimension 2")+
  theme_bw()+theme_tsne
tsne_plt
ggsave("./Final_Plots/tsne_FBA_methane.pdf", width = 3, height = 2.75, units = "in", dpi = "retina")
write.csv(tsne_df, "./Source_Data/Figure_4C.csv")

#### Bar plots for Methanogenesis activity (Flux)
#remotes::install_github("coolbutuseless/ggpattern")
library("ggpattern")
### Load FBA Flux outputs
methano_flux<-read.csv("./Final_Data_Files/compiled_methanogenesis_flux_revised.csv", header = T, check.names = T)
methano_flux$Type_Desc<-paste(methano_flux$Type, methano_flux$Desc, sep = "_")

theme_barplot<-theme(plot.title = element_text(angle = 0, size = 4.5, face = "italic", hjust = .8, vjust = -4),
                  axis.text = element_text(angle = 0, size = 5),
                  axis.title= element_text(angle = 0, size = 5.5),
                  strip.text = element_text(angle = 0, size = 4.5, color = "black", margin = margin(c(.02,.02,.02,.02), unit='cm'), face = "italic"),
                  #strip.background = element_rect(fill= "white", color = "black", size =.25),
                  strip.background = element_blank(),
                  legend.box.background = element_rect(fill= "white", color = "grey", size =.25), 
                  legend.title = element_text(size=5),
                  legend.position = "none", legend.text = element_text(size = 5), # c(.09,.89)
                  legend.key.size = unit(.02, "in"), legend.key.width = unit(.02, "in"),
                  legend.justification = 'center', legend.title.align=0.5, legend.margin=margin(c(.05,.05,.05,.05), unit='cm'),
                  axis.ticks = element_line(linewidth = .25),
                  panel.spacing.y = unit(.3, "lines"))

Lac_flux<-ggplot(data = NULL)+
  geom_col(data =subset(methano_flux, methano_flux$Metabolite %in% "Lactate"),
           aes(x = factor(Days),y= Flux, fill = factor(Type_Desc)), 
           position = position_dodge(preserve = "single"),  color = "Black", linewidth = .2)+
  scale_fill_manual(values = c("Planktonic_Total" =  alpha("#1f78b4", 1), "Sediment_Total" = alpha("#a6cee3", 1),
                               "Planktonic_Used" =  "#1f78b4", "Sediment_Used" = "#a6cee3"), guide =NULL, name = "Phase")+
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "black",
                                                 fill = c("P" =  "#1f78b4", "S" = "#a6cee3"))))+
  scale_x_discrete(name = "Day (Fluidization)")+
  scale_y_continuous(name = "Lactate Utilization (mMol/gDCW/h)", limits = c(0, 33),)+
  ggtitle("Dv")+
  theme_bw()+theme_barplot
Lac_flux

### Significance test
Mean_Pct_Lactate<-t.test(subset(methano_flux, methano_flux$Metabolite %in% "Lactate" & methano_flux$Desc %in% "Used")$Flux[1:6],
                          subset(methano_flux, methano_flux$Metabolite %in% "Lactate" & methano_flux$Desc %in% "Used")$Flux[7:12])
Mean_Pct_Lactate<-cbind(t(Mean_Pct_Lactate$estimate), "p.value" = Mean_Pct_Lactate$p.value)
colnames(Mean_Pct_Lactate)[1:2]<-c("P_pct_Lactate", "S_pct_Lactate")
Mean_Pct_Lactate<-data.frame(Mean_Pct_Lactate)
Mean_Pct_Lactate

H2_flux<-ggplot(data = NULL)+
  geom_col_pattern(data =subset(methano_flux, methano_flux$Metabolite %in% "Hydrogen" & methano_flux$Desc %in% "Total"),
           aes(x = factor(Days), y= Flux, fill = factor(Type_Desc)),
           position = position_dodge(preserve = "single"), linewidth = .2,
           pattern = "stripe", #Pattern name string e.g. 'stripe' (default), 'crosshatch', 'point', 'circle', 'none'                 
           color = "black", 
           pattern_fill = "white",
           pattern_angle = 45, 
           pattern_density = 0.01, #Approximate fill fraction of the pattern. Usually in range [0, 1], but can be higher. default: 0.2
           pattern_spacing = 0.04, #Spacing of the pattern as a fraction of the plot size. default: 0.05
           pattern_size =.2) + 
  geom_col(data =subset(methano_flux, methano_flux$Metabolite %in% "Hydrogen" & methano_flux$Desc %in% "Used"),
           aes(x = factor(Days),y= Flux, fill = factor(Type_Desc)), 
           position = position_dodge(preserve = "single"),  color = "Black", linewidth = .2)+
  scale_fill_manual(values = c("Planktonic_Total" =  alpha("#1f78b4", .9), "Sediment_Total" = alpha("#a6cee3",.9),
                               "Planktonic_Used" =  "#1f78b4", "Sediment_Used" = "#a6cee3"), guide =NULL, name = "Phase")+
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "black",
                                                fill = c("P" =  "#1f78b4", "S" = "#a6cee3"))))+
  scale_x_discrete(name = "Day (Fluidization)")+
  scale_y_continuous(name = "Hydrogen (mMol/gDCW/h)", limits = c(0, 65), breaks = seq(0,60, by = 15))+
  ggtitle("Dv & Mm")+
  theme_bw()+theme_barplot
H2_flux

### Significance test
Mean_Pct_Unsed_H2<-t.test(round((subset(methano_flux, methano_flux$Metabolite %in% "Hydrogen" & methano_flux$Desc %in% "Unused")$Flux/
                     subset(methano_flux, methano_flux$Metabolite %in% "Hydrogen" & methano_flux$Desc %in% "Total")$Flux)*100, 2)[1:6],
       round((subset(methano_flux, methano_flux$Metabolite %in% "Hydrogen" & methano_flux$Desc %in% "Unused")$Flux/
                subset(methano_flux, methano_flux$Metabolite %in% "Hydrogen" & methano_flux$Desc %in% "Total")$Flux)*100, 2)[7:12]
       )
Mean_Pct_Unsed_H2<-cbind(t(Mean_Pct_Unsed_H2$estimate), "p.value" = Mean_Pct_Unsed_H2$p.value)
colnames(Mean_Pct_Unsed_H2)[1:2]<-c("P_pct_unused", "S_pct_unused")  
Mean_Pct_Unsed_H2<-data.frame(Mean_Pct_Unsed_H2)
Mean_Pct_Unsed_H2
  
Formate_flux<-ggplot(data = NULL)+
  geom_col_pattern(data =subset(methano_flux, methano_flux$Metabolite %in% "Formate" & methano_flux$Desc %in% "Total"),
                   aes(x = factor(Days), y= Flux, fill = factor(Type_Desc)),
                   position = position_dodge(preserve = "single"), linewidth = .2,
                   pattern = "stripe", #Pattern name string e.g. 'stripe' (default), 'crosshatch', 'point', 'circle', 'none'                 
                   color = "black", 
                   pattern_fill = "white",
                   pattern_angle = 45, 
                   pattern_density = 0.01, #Approximate fill fraction of the pattern. Usually in range [0, 1], but can be higher. default: 0.2
                   pattern_spacing = 0.04, #Spacing of the pattern as a fraction of the plot size. default: 0.05
                   pattern_size =.2) + 
  geom_col(data =subset(methano_flux, methano_flux$Metabolite %in% "Formate" & methano_flux$Desc %in% "Used"),
           aes(x = factor(Days),y= Flux, fill = factor(Type_Desc)), 
           position = position_dodge(preserve = "single"),  color = "Black", linewidth = .2)+
  scale_fill_manual(values = c("Planktonic_Total" =  alpha("#1f78b4", .9), "Sediment_Total" = alpha("#a6cee3",.9),
                               "Planktonic_Used" =  "#1f78b4", "Sediment_Used" = "#a6cee3"), guide =NULL, name = "Phase")+
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "black",
                                                 fill = c("P" =  "#1f78b4", "S" = "#a6cee3"))))+
  scale_x_discrete(name = "Day (Fluidization)")+
  scale_y_continuous(name = "Formate (mMol/gDCW/h)", limits = c(0, .4))+
  ggtitle("Dv & Mm")+
  theme_bw()+theme_barplot
Formate_flux

Acetate_flux<-ggplot(data = NULL)+
  geom_col_pattern(data =subset(methano_flux, methano_flux$Metabolite %in% "Acetate" & methano_flux$Desc %in% "Total"),
                   aes(x = factor(Days), y= Flux, fill = factor(Type_Desc)),
                   position = position_dodge(preserve = "single"), linewidth = .2,
                   pattern = "stripe", #Pattern name string e.g. 'stripe' (default), 'crosshatch', 'point', 'circle', 'none'                 
                   color = "black", 
                   pattern_fill = "white",
                   pattern_angle = 45, 
                   pattern_density = 0.01, #Approximate fill fraction of the pattern. Usually in range [0, 1], but can be higher. default: 0.2
                   pattern_spacing = 0.04, #Spacing of the pattern as a fraction of the plot size. default: 0.05
                   pattern_size =.2) + 
  geom_col(data =subset(methano_flux, methano_flux$Metabolite %in% "Acetate" & methano_flux$Desc %in% "Used"),
           aes(x = factor(Days),y= Flux, fill = factor(Type_Desc)), 
           position = position_dodge(preserve = "single"),  color = "Black", linewidth = .2)+
  scale_fill_manual(values = c("Planktonic_Total" =  alpha("#1f78b4", .9), "Sediment_Total" = alpha("#a6cee3",.9),
                               "Planktonic_Used" =  "#1f78b4", "Sediment_Used" = "#a6cee3"), guide =NULL, name = "Phase")+
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "black",
                                                 fill = c("P" =  "#1f78b4", "S" = "#a6cee3"))))+
  scale_x_discrete(name = "Day (Fluidization)")+
  scale_y_continuous(name = "Acetate (mMol/gDCW/h)", limits = c(0, 33))+
  ggtitle("Dv & Mm")+
  theme_bw()+theme_barplot
Acetate_flux

CH4_flux<-ggplot(data = NULL)+
  geom_col(data =subset(methano_flux, methano_flux$Metabolite %in% "Methane"),
           aes(x = factor(Days),y= Flux, fill = factor(Type_Desc)), 
           position = position_dodge(preserve = "single"),  color = "Black", linewidth = .2)+
  scale_fill_manual(values = c("Planktonic_Total" =  alpha("#1f78b4", 1), "Sediment_Total" = alpha("#a6cee3",1),
                               "Planktonic_Used" =  "#1f78b4", "Sediment_Used" = "#a6cee3"), guide =NULL, name = "Phase")+
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "black",
                                                 fill = c("P" =  "#1f78b4", "S" = "#a6cee3"))))+
  scale_x_discrete(name = "Day (Fluidization)")+
  scale_y_continuous(name = "CH4 (mMol/gDCW/h)", limits = c(0, 12), breaks = seq(0,12, by = 3))+
  ggtitle("Mm")+
  theme_bw()+theme_barplot
CH4_flux

Biomass_flux<-ggplot(data = NULL)+
  geom_col(data =subset(methano_flux, methano_flux$Metabolite %in% c("DvH_Biomass", "Mmp_Biomass")),
           aes(x = factor(Days),y= Flux, fill = factor(Type_Desc)), 
           position = position_dodge(preserve = "single"),  color = "Black", linewidth = .2)+
  scale_fill_manual(values = c("Planktonic_Total" =  alpha("#1f78b4", 1), "Sediment_Total" = alpha("#a6cee3",1),
                               "Planktonic_Used" =  "#1f78b4", "Sediment_Used" = "#a6cee3"), guide =NULL, name = "Phase")+
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "black",
                                                 fill = c("P" =  "#1f78b4", "S" = "#a6cee3"))))+
  scale_x_discrete(name = "Days")+
  scale_y_continuous(name = "Biomass (h-1)")+
  facet_wrap(~Metabolite, ncol = 2, scales = "free_y")+
  theme_bw()+theme_barplot
Biomass_flux
ggsave(plot =  Biomass_flux, filename =  "./Final_Plots/SynCom_Biomass_barplot.pdf", width = 3.2, height = 2, units = "in", dpi = "retina")

Methano_barplot<-ggarrange(Lac_flux,H2_flux, CH4_flux,  ncol = 3) # vertical
ggsave(plot =  Methano_barplot, filename =  "./Final_Plots/Syntrophy_barplot.pdf", width = 3.5, height = 2, units = "in", dpi = "retina")
write.csv(methano_flux, "./Source_Data/Figure4DEFGH.csv")
#### CORRELATION PLOTS BI-LEVEL OPTIMIZATION
#Bio_Constrain_df<-read.csv("./FBA_data/CorrelationPlot-BilevelOptimization.csv", header = T, check.names = T) ##old
Bio_Constrain_df<-read.csv("./Final_Data_Files/CorrelationPlot-Pre_Post_BiomassConstraints.csv", header = T, check.names = T) ##old
#methano_flux$Flux2<-methano_flux$Flux
Bio_Constrain_df<-cbind(Bio_Constrain_df, colsplit(Bio_Constrain_df$label, pattern = "_", names = c("Sample", "Org")))
Bio_Constrain_df$Day<-gsub("Day", "", Bio_Constrain_df$Sample)
Bio_Constrain_df$Days<-gsub("p","",Bio_Constrain_df$Day)
Bio_Constrain_df$Days<-gsub("s","",Bio_Constrain_df$Days)
Bio_Constrain_df$Phase<-gsub("[[:digit:]]","",Bio_Constrain_df$Day)
Bio_Constrain_df$Phase<-gsub("p","P",Bio_Constrain_df$Phase); Bio_Constrain_df$Phase<-gsub("s","S",Bio_Constrain_df$Phase)

Bio_Constrain_plt<-ggplot(data = Bio_Constrain_df)+
  geom_abline(intercept = 0, slope = .5/.5, color ="darkgrey", lwd =.5,  linetype =3, alpha = .6)+

  stat_poly_line(
                 aes(x=Prediction1, y = Experiment), color ="#A09CC0", se =T, linetype =2, fill = "#A09CC0", alpha =.2, linewidth =.5) +
  stat_poly_eq(
    aes(x=Prediction1, y = Experiment, label = after_stat(eq.label)), label.y = 0.14, label.x = .93, size = 2.5) +
  stat_poly_eq(
    aes(x=Prediction1, y = Experiment), label.y = 0.09, label.x = .93, size = 2.5) +
  geom_point(
    aes(x = Prediction1,y= Experiment, fill = factor(Phase), shape = factor(Days), color = factor(Sample)), size = 1.5, stroke =.2, alpha =.3)+
  stat_poly_line(
                 aes(x=Prediction2, y = Experiment), color ="Black", se =T, linetype =1, fill = "grey", alpha =.2, linewidth =.5) +
  stat_poly_eq(
    aes(x=Prediction2, y = Experiment, label = after_stat(eq.label)), label.y = .9, label.x = .3, size = 2.5) +
  stat_poly_eq(
    aes(x=Prediction2, y = Experiment), label.y = .85, label.x = .3, size = 2.5) +###
  geom_point(
    aes(x = Prediction2,y= Experiment, fill = factor(Phase), shape = factor(Days), color = factor(Sample)), size = 2.5, stroke =.5)+ 
  scale_fill_manual(values = c("P" =  "#1f78b4", "S" = "#a6cee3"), guide =NULL, name = "Phase")+
  scale_color_manual(values = c("Day1p" =  "#1f78b4", "Day1s" = "#a6cee3",
                                "Day2p" =  "black", "Day2s" = "black",
                                "Day3p" =  "black", "Day3s" = "black",
                                "Day4p" =  "black", "Day4s" = "black",
                                "Day5p" =  "black", "Day5s" = "black",
                                "Day6p" =  "black", "Day6s" = "black"), guide =NULL, name = "Phase")+
  scale_shape_manual(values = c(8,21:25), guide =NULL, name = "Day")+
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "black",size = 1.5, stroke = .2,
                                                 fill = c("P" =  "#1f78b4", "S" = "#a6cee3"))),
         shape =guide_legend(override.aes = list(shape = c(8,21:25), colour = "black", size = 1.5, stroke = .2,
                                                 fill = c("1" =  "white", "2" =  "white", "3" =  "white", "4" =  "white", "5" =  "white", "6" =  "white"))))+
  scale_x_continuous(limits = c(.25,1.05), breaks = seq(.25,1, by =.25), name = "DvH Model Prediction Biomass Fraction ")+
  scale_y_continuous(limits = c(.25,1.05), breaks = seq(.25,1, by =.25), name = "DvH Experimental Biomass Fraction")+
  #scale_y_continuous(expand = c(.1,.1), name = "tSNE dimension 2")+
  theme_bw()+theme_tsne
Bio_Constrain_plt
ggsave(plot =  Bio_Constrain_plt, filename =  "./Final_Plots/Bio_Constrain_correlation_plt_combined_Pre_Post.pdf", width = 3, height = 2.75, units = "in", dpi = "retina")
write.csv(Bio_Constrain_df, "./Source_Data/Figure_4B.csv")

