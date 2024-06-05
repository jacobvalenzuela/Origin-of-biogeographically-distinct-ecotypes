#rm(list=ls())
getwd()
setwd("/Users/jake/Documents/R_analysis/R_Analysis/ENIGMA_Analysis/FBR_DESeq/Final_Code/")
getwd()

Lineplot_Packages<-c("dplyr", "FactoMineR", "factoextra", "ggplot2", "graphics", "grDevices", "matrixStats", "pheatmap", 
            "RColorBrewer","readr", "reshape2", "stats", "SummarizedExperiment", "tidyr","textshape",
            "readr", "tidyverse", "corrplot", "ggpubr", "pals", "ComplexHeatmap", "circlize", "paletteer")
#install.packages(Lineplot_Packages)
lapply(Lineplot_Packages, library, character.only =T)

#####
"%notin%" <- Negate("%in%")
##### functions and objects used for line plots
`%notin%` <- Negate(`%in%`)

# Function for median and IQR
IQR <- function(x) {
  data.frame(
    ymin = quantile(x)[2], # 1st quartile
    ymax = quantile(x)[4])  # 3rd quartile
}

# Function for mean 
mean_fun <- function(x) {
  data.frame(y = mean(x))  #mean
}

# Function for sd 
sd_range <- function(x) {
  data.frame(
    ymin = mean(x) - sd(x) ,
    ymax = mean(x) + sd(x) ) 
}
# Function for min, max values
range <- function(x) {
  data.frame(ymin=min(x),
             ymax=max(x))
}
# Themes for enrichment  plots
theme_enrichment<-theme(plot.title = element_text(angle = 0, size = 8),
                        axis.text = element_text(angle = 0, size = 6),
                        axis.title= element_text(angle = 0, size = 8),
                        strip.text = element_text(angle = 0, size = 5, color = "black", margin = margin(c(.07,.07,.07,.07), unit='cm')),
                        legend.box.background = element_rect(fill= "white", color = "grey", size =.25), 
                        legend.title = element_text(size=5),
                        legend.position = "none", legend.text = element_text(size = 5), # c(.09,.89)
                        legend.key.size = unit(.02, "in"), legend.key.width = unit(.02, "in"),
                        legend.justification = 'center', legend.title.align=0.5, legend.margin=margin(c(.05,.05,.05,.05), unit='cm'))
#### Run DvH_Mmp_DESeq_Final.R first to generate TPM files
#### making line plots for subsets of genes
#### build a dataframe for specific labels
gene_ID_labels<-read.csv("./Final_Data_Files/SynCom_Ordered_index_processed.csv")
Flagella<-subset(gene_ID_labels, gene_ID_labels$Note %in% c("Flagella"))
Lactate_transport<-subset(gene_ID_labels, gene_ID_labels$Note %in% c("Lactate"))
Hydrogenases_LA_mmp<-subset(gene_ID_labels, gene_ID_labels$Note %in% c("Hydrogenase(LA)") & gene_ID_labels$Org %in% c("Mmp"))
Hydrogenases_HA_mmp<-subset(gene_ID_labels, gene_ID_labels$Note %in% c("Hydrogenase") & gene_ID_labels$Org %in% c("Mmp"))
Hydrogenases_DvH<-subset(gene_ID_labels, gene_ID_labels$Note %in% c("Hydrogenase") & gene_ID_labels$Org %in% c("DvH"))

#### Enriched Gene Plots (z-scores)
### build data.frame for enrichment genes
Enrich_df<-rbind(cbind(subset(DvH_long_df_zscore, DvH_long_df_zscore$Gene_ID %in% Flagella$GeneID), data.frame("Enrichment" = "Flagella", "Org" = "DvH")),
cbind(subset(Mmp_long_df_zscore, Mmp_long_df_zscore$Gene_ID %in% Hydrogenases_HA_mmp$GeneID), data.frame("Enrichment" = "Hydrogenases (Mmp)", "Org" = "Mmp")),
cbind(subset(DvH_long_df_zscore, DvH_long_df_zscore$Gene_ID %in% Hydrogenases_DvH$GeneID), data.frame("Enrichment" = "Hydrogenases (DvH)", "Org" = "DvH")),
cbind(subset(Mmp_long_df_zscore, Mmp_long_df_zscore$Gene_ID %in% Hydrogenases_LA_mmp$GeneID), data.frame("Enrichment" = "Low Affinity Hydrogenase (Mmp)", "Org" = "Mmp")),
cbind(subset(DvH_long_df_zscore, DvH_long_df_zscore$Gene_ID %in% Lactate_transport$GeneID), data.frame("Enrichment" = "Lactate Permeases", "Org" = "DvH")))
Enrich_df$Enrichment<-factor(Enrich_df$Enrichment, levels = c("Lactate Permeases","Flagella", "Hydrogenases (DvH)", "Hydrogenases (Mmp)", "Low Affinity Hydrogenase (Mmp)" ))

pos.d<-position_dodge(.6)
Enrich_plt<-ggplot(Enrich_df, 
                  aes(x=factor(Day), y=zScore, fill=Phase, color=Phase, group = Phase))+
  stat_summary(geom = "path",
               fun.data = mean_fun, 
               position = pos.d,
               size=.3, color = "black")+
  stat_summary(geom = "linerange",
               fun.data = sd_range, 
               position = pos.d,
               size=.2, color = "black")+
  geom_point(position=pos.d, 
              alpha = 1, size = .5, shape = 20)+
  stat_summary(geom = "point",
               fun.data = mean_fun, 
               position = pos.d,
               size=1.8, shape = 21, color = "black")+
  scale_color_manual(values= c(P = "#1f78b4", S = "#a6cee3"))+
  scale_fill_manual(values= c(P = "#1f78b4", S = "#a6cee3"),guide=NULL)+
  ylab("Expression (z-score)")+xlab("Days")+
  #ggtitle("Flagellum associated genes (DvH)")+
  facet_grid(rows = vars(Enrichment))+
  theme_light()+theme_enrichment
Enrich_plt
### color in the strip backgrounds https://stackoverflow.com/questions/19440069/ggplot2-facet-wrap-strip-color-based-on-variable-in-data-set
g <- ggplot_gtable(ggplot_build(Enrich_plt))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("#D94801", "#D94801", "#D94801", "#BEBEBE", "#BEBEBE")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
pdf("./Final_Plots/SynCom_Enrich_Zscore_2.pdf", width = 2.25, height = 6.6)
grid::grid.draw(g)
dev.off()

summs_emrichment<-Enrich_df %>%
  group_by(Enrichment) %>%
  summarise(across(c(Gene_ID),
                   list("n" = length))) 
summs_emrichment$Gene_ID_n<-summs_emrichment$Gene_ID_n/12
summs_emrichment
write.csv(Enrich_df, "./Source_Data/Figure_3C.csv")
#######
### build metadata data frame
LacAceMeth<-read.csv(file = "./Final_Data_Files/FBR_metadata.csv")
colnames(LacAceMeth)
summs_meta<-LacAceMeth %>%
  group_by(Sample_type, Days, Hours) %>%
  summarise(across(c(OD_Balch, pH, DO, mM_Lactate, mM_Acetate, mM_Methane, ng_DNA_per_ml, ng_DNA_per_gram, ng_RNA_per_ml, ng_RNA_per_gram, ug_mL_Protein, ng_Protein_per_mg),
                   list("mean" = mean, "sd" = sd, "n" = length))) 
summs_meta$Days2<-summs_meta$Days-1
colnames(summs_meta)

RNA_times<-summs_meta$Hours[3:8]
point_size<-2.5
line.size<-.3
point_size_RNA<-1

theme_meta<-theme(plot.title = element_text(angle = 0, size = 8, face = "plain", hjust = .5, vjust = -2),
                        axis.text = element_text(angle = 0, size = 6),
                        axis.title= element_text(angle = 0, size = 8),
                        strip.text = element_text(angle = 0, size = 5, color = "black", margin = margin(c(.07,.07,.07,.07), unit='cm')),
                        legend.box.background = element_rect(fill= "white", color = "grey", size =.25), 
                        legend.title = element_text(size=5),
                        legend.position = "none", legend.text = element_text(size = 5), # c(.09,.89)
                        legend.key.size = unit(.02, "in"), legend.key.width = unit(.02, "in"),
                        legend.justification = 'center', legend.title.align=0.5, legend.margin=margin(c(.05,.05,.05,.05), unit='cm'))

OD_plt<-ggplot(summs_meta, 
                   aes(x=Hours, fill=Sample_type, group = Sample_type))+
  geom_rect(aes(xmin = 0, xmax = 46.5, ymin = -Inf, ymax = Inf), fill = "grey90", alpha =.05)+
  geom_linerange(aes(ymin = OD_Balch_mean - OD_Balch_sd ,
                 ymax = OD_Balch_mean + OD_Balch_sd), color = "Black", size = line.size)+
  geom_path(aes(y = OD_Balch_mean), color="black", size = line.size)+
  geom_point(aes(y = OD_Balch_mean), size = point_size, shape = 21)+
  annotate(geom = "point", x = RNA_times, y = .63, shape = 25, fill = "grey", size = point_size_RNA, stroke = .2)+
  scale_color_manual(values= c(Planktonic = "#1f78b4", Sediment = "#a6cee3"), guide=NULL)+
  scale_fill_manual(values= c(Planktonic = "#1f78b4", Sediment = "#a6cee3"), guide = NULL)+
  scale_x_continuous(limits = c(0,168), breaks = seq(0,168, by = 24),
                     labels = c(1,2,1:6))+
  scale_y_continuous(limits = c(0,.63), breaks = seq(0,.6, by = .1))+
  ylab("O.D. 600 nm")+xlab(NULL)+
  ggtitle("Planktonic")+
  theme_light()+theme_meta+theme(axis.text.x=element_blank())
OD_plt

protein_plt<-ggplot(summs_meta, 
               aes(x=Hours, fill=Sample_type, group = Sample_type))+
  geom_rect(aes(xmin = 0, xmax = 46.5, ymin = -Inf, ymax = Inf), fill = "grey90", alpha =.05)+
  geom_linerange(aes(ymin = ng_Protein_per_mg_mean - ng_Protein_per_mg_sd,
                     ymax = ng_Protein_per_mg_mean + ng_Protein_per_mg_sd), color = "Black", size = line.size)+
  geom_path(aes(y = ng_Protein_per_mg_mean), color="black", size = line.size)+
  geom_point(aes(y = ng_Protein_per_mg_mean), size = point_size, shape = 21)+
  annotate(geom = "point", x = RNA_times, y = 23, shape = 25, fill = "grey", size = point_size_RNA, stroke = .2)+
  scale_color_manual(values= c(Planktonic = "#1f78b4", Sediment = "#a6cee3"), guide=NULL)+
  scale_fill_manual(values= c(Planktonic = "#1f78b4", Sediment = "#a6cee3"), guide = NULL)+
  scale_x_continuous(limits = c(0,168), breaks = seq(0,168, by = 24),
                     labels = c(1,2,1:6))+
  scale_y_continuous(limits = c(0,23), breaks = seq(0,20, by = 5))+
  ylab("ng Protein / mg Sediment")+xlab(NULL)+
  ggtitle("Sediment (attached)")+
  theme_light()+theme_meta+theme(axis.text.x=element_blank())
protein_plt

organics_plt<-ggplot(summs_meta[!is.na(summs_meta$mM_Lactate_mean),], aes(x=Hours))+
  geom_rect(aes(xmin = 0, xmax = 46.5, ymin = -Inf, ymax = Inf), fill = "grey90", alpha =.12)+
  geom_linerange(aes(ymin = mM_Lactate_mean - mM_Lactate_sd,
                     ymax = mM_Lactate_mean + mM_Lactate_sd), color = "Black", size = line.size)+
  geom_linerange(aes(ymin = mM_Acetate_mean - mM_Acetate_sd,
                     ymax = mM_Acetate_mean + mM_Acetate_sd), color = "Black", size = line.size)+
  geom_path(aes(y = mM_Lactate_mean), color="black", size = line.size)+
  geom_point(aes(y = mM_Lactate_mean), fill = "#006666", size = point_size, shape = 22)+
  geom_path(aes(y = mM_Acetate_mean), color="black", size = line.size)+
  geom_point(aes(y = mM_Acetate_mean), fill = "#009999", size = point_size, shape = 22)+
  annotate(geom = "point", x = RNA_times, y = 45, shape = 25, fill = "grey", size = point_size_RNA, stroke = .2)+
  scale_x_continuous(limits = c(0,168), breaks = seq(0,168, by = 24),
                     labels = c(1,2,1:6), name = "Days (Fluidization)")+
  scale_y_continuous(limits = c(0,45), breaks = seq(0,45, by = 10))+
  ylab("mM Lactate & Acetate")+xlab("Hours")+
  annotate(geom= "rect", xmin = c(.75), xmax=c(32), ymin = c(40.75), ymax = c(43.5), 
           color = "black", fill = c("white"), size = .2)+
  annotate(geom= "point", x= c(5,5), y = c(42.75, 41.5), size = point_size_RNA, stroke = .2, 
           color = "black", fill = c("#006666", "#009999"), shape =c(22,22))+
  annotate(geom= "text", x= c(20,20), y = c(42.75, 41.5), label = c("Lactate", "Acetate"), 
           color = "black",  size = 1)+
  ggtitle("Reactor Organics")+
  theme_light()+theme_meta+theme(legend.position = c(.1,.9))
organics_plt

meth_plt<-ggplot(summs_meta, aes(x=Hours))+
  geom_rect(aes(xmin = 0, xmax = 46.5, ymin = -Inf, ymax = Inf), fill = "grey90", alpha =.05)+
  geom_linerange(aes(ymin = mM_Methane_mean - mM_Methane_sd,
                     ymax = mM_Methane_mean + mM_Methane_sd), color = "Black", size = line.size)+
  geom_path(aes(y = mM_Methane_mean), color="black", size = line.size)+
  geom_point(aes(y = mM_Methane_mean), fill = "#662700", size = point_size, shape = 24)+
  annotate(geom = "point", x = RNA_times, y = 2, shape = 25, fill = "grey", size = point_size_RNA, stroke =.2)+
  scale_x_continuous(limits = c(0,168), breaks = seq(0,168, by = 24),
                     labels = c(1,2,1:6), name = "Days (Fluidization)")+
  scale_y_continuous(limits = c(0,2), breaks = seq(0,2, by = .5))+
  ylab("mM Methane")+xlab("Hours")+
  ggtitle("Reactor Methane")+
  theme_light()+theme_meta+theme(legend.position = c(.1,.9))
meth_plt

fig_1_meta<-ggarrange(OD_plt, protein_plt, organics_plt,meth_plt, 
                      labels = c("B)", "C)", "D)", "E)"),  font.label = list(size = 9, color = "black", face = "plain"),
                      nrow =2, ncol = 2, heights = c(2,2.05))
fig_1_meta
ggsave(fig_1_meta, file ="./Final_Plots/SynCom_metadata_lineplots.pdf", width = 3.75, height = 4.5, units = "in", dpi = "retina")
dev.off()
write.csv(summs_meta, "./Source_Data/Figure_1BCDE.csv")
### calculating growth rate
### build the model
#### testing growthrates package 
#y0 = initial abundance
#mumax = maximum growth rate (1/time), mass of cells produced/ original mass of cells * time
#K = carrying capacity (max. abundance)
#y_shift = (y_shift) h0 parameter specifying the initial physiological state of organisms (e.g. cells) and in consequence the lag phase (h0 = max growth rate * lag phase)
library("growthrates")
grow_logistic_yshift <- function(time, parms) {
  with(as.list(parms), {
    y <- (K * y0) / (y0 + (K - y0) * exp(-mumax * time)) + y_shift
    as.matrix(data.frame(time = time, y = y))
  })
}
grow_logistic_yshift <- growthmodel(grow_logistic_yshift,
                                    c("y0", "mumax", "K", "y_shift"))
#construct your fit model based on your x and y inputs and estimate your parameters
y = LacAceMeth$OD_Balch
x = LacAceMeth$Hours/24

y2 = LacAceMeth$ng_Protein_per_mg
x2 = LacAceMeth$Hours/24

fit <- fit_growthmodel(grow_logistic_yshift,
                       p = c(y0 = min(y, na.rm = T), 
                             mumax = .7,  
                             K = max(y, na.rm = T), 
                             y_shift = 0),
                       time = x, y = y)

fit2 <- fit_growthmodel(grow_logistic_yshift,
                        p = c(y0 = min(y2, na.rm = T), 
                              mumax = .7,  
                              K = max(y2, na.rm = T), 
                              y_shift = 0),
                        time = x2, y = y2)
plot(fit2)

mu_max<-fit@fit$par["mumax"]
doubling_time<-log(2)/mu_max  

growth_stats1<-data.frame("mumax" =fit@fit$par["mumax"], 
                          "doubling_time" = log(2)/fit@fit$par["mumax"]) 
growth_stats2<-data.frame("mumax" =fit2@fit$par["mumax"], 
                          "doubling_time" = log(2)/fit2@fit$par["mumax"]) 
growth_stats<-rbind(growth_stats1, growth_stats2)
rownames(growth_stats)<-c("Planktonic (OD)", "Sediment (ng/mg)")
growth_stats

#individual plots
pH_DO_plt<-ggplot(summs_meta, aes(x=Hours))+
  geom_rect(aes(xmin = 0, xmax = 46.5, ymin = -Inf, ymax = Inf), fill = "grey90", alpha =.05)+
  geom_linerange(aes(ymin = pH_mean - pH_sd,
                     ymax = pH_mean + pH_sd), color = "Black", size = line.size)+
  geom_linerange(aes(ymin = DO_mean - DO_sd,
                     ymax = DO_mean + DO_sd), color = "Black", size = line.size)+
  geom_path(aes(y = pH_mean), color="black", size = line.size)+
  geom_path(aes(y = DO_mean), color="black", size = line.size)+
  geom_point(aes(y = pH_mean), fill = "#F28E2BFF", size = point_size, shape = 23)+
  geom_point(aes(y = DO_mean), fill = "#59A14FFF", size = point_size, shape = 23)+
  annotate(geom = "point", x = RNA_times, y = 8.1, shape = 25, fill = "grey", size = point_size_RNA, stroke =.2)+
  scale_x_continuous(limits = c(0,168), breaks = seq(0,168, by = 24),
                     labels = c(1,2,1:6), name = "Days (Fluidization)")+
  scale_y_continuous(limits = c(0,8.5), breaks = seq(0,8.5, by = 1))+
  ylab("Dissolved Oxygen (mg/L) & pH")+xlab("Hours")+
  annotate(geom= "rect", xmin = c(.75), xmax=c(32), ymin = c(.35), ymax = c(.9), 
           color = "black", fill = c("white"), size = .2)+
  annotate(geom= "point", x= c(5,5), y = c(.75, .5), size = point_size_RNA, stroke = .2, 
           color = "black", fill = c("#F28E2BFF", "#59A14FFF"), shape =c(23,23))+
  annotate(geom= "text", x= c(13,20), y = c(.75, .5), label = c("pH", "DO (mg/L)"), 
           color = "black",  size = 1)+
  theme_light()+theme_meta+theme(legend.position = c(.1,.9))
pH_DO_plt
ggsave(pH_DO_plt, file ="./Final_Plots/SynCom_pH_DO_plt.pdf", width = 2.25, height = 2.5,units = "in", dpi = "retina")

##### MMP archaellum USE TPM count matrix 
Mmp_cts_1<-Mmp_cts+1
Mmp_cts_z<-scale(Mmp_cts+1) ### From TPM expression matrix
Mmp_flaB_cts<-melt(subset(Mmp_cts_1, rownames(Mmp_cts_1) %in% "MMP1667"), varnames = c("GeneID", "Sample"), value.name = "counts")
#Mmp_flaB_cts<-melt(subset(Mmp_cts_z, rownames(Mmp_cts_z) %in% "MMP1667"), varnames = c("GeneID", "Sample"), value.name = "counts")
Mmp_flaB_cts<-cbind(Mmp_flaB_cts, colsplit(Mmp_flaB_cts$Sample, pattern = "_", names = c("EPD", "Phase", "Reactor", "Date")))
Mmp_flaB_cts$Day<-Mmp_flaB_cts$Date
Mmp_flaB_cts$Day<-gsub("5-8","1", Mmp_flaB_cts$Day)
Mmp_flaB_cts$Day<-gsub("5-9","2", Mmp_flaB_cts$Day)
Mmp_flaB_cts$Day<-gsub("5-10","3", Mmp_flaB_cts$Day)
Mmp_flaB_cts$Day<-gsub("5-11","4", Mmp_flaB_cts$Day)
Mmp_flaB_cts$Day<-gsub("5-12","5", Mmp_flaB_cts$Day)
Mmp_flaB_cts$Day<-gsub("5-13","6", Mmp_flaB_cts$Day)

### Box plots Day 3-6 
Day_Range<-c(3:6)
flab_last_days<- t.test(Mmp_flaB_cts$counts[Mmp_flaB_cts$Phase == "P" & Mmp_flaB_cts$Day %in% Day_Range],
                        Mmp_flaB_cts$counts[Mmp_flaB_cts$Phase == "S" & Mmp_flaB_cts$Day %in% Day_Range],
                          var.equal = TRUE )
flab_last_days$p.value
flab_last_day_FC<-log2(mean(Mmp_flaB_cts$counts[Mmp_flaB_cts$Phase == "P" & Mmp_flaB_cts$Day %in% Day_Range])/
       mean(Mmp_flaB_cts$counts[Mmp_flaB_cts$Phase == "S" & Mmp_flaB_cts$Day %in% Day_Range]))
flab_last_day_FC; flab_last_days$p.value

MMP1667_cts<-ggplot(subset(Mmp_flaB_cts, Mmp_flaB_cts$Day %in% Day_Range), 
                    aes(x=factor(Phase), y=counts, fill=Phase, group = Phase))+
  geom_boxplot(color ="black", size =.25)+
  geom_jitter(aes(shape = Reactor), color ="black", width = .1,  stroke = .15, size = 1.2)+
  #geom_text(aes(x =2., y=145, label = paste("pvalue =", round(Mmp_res_all["MMP1667",]$pvalue,2), sep = " ")), check_overlap = T, size =1.7)+
  geom_text(aes(x =1.9, y=109, label = paste("pvalue =", round(flab_last_days$p.value,4), sep = " ")), check_overlap = T, size =1.9)+
  geom_text(aes(x =1.9, y=117, label = paste("log2(FC: P/S) =", round(flab_last_day_FC,4), sep = " ")), check_overlap = T, size =1.9)+
  scale_color_manual(values= c(P = "#1f78b4", S = "#a6cee3"))+
  scale_fill_manual(values= c(P = "#1f78b4", S = "#a6cee3"),guide=NULL)+
  scale_x_discrete(labels = c("P" = "Planktonic", "S" = "Sediment"))+
  scale_y_continuous(name = "Expression (TPMs + 1)")+
  scale_shape_manual(values = c(21:23))+
  ylab("Expression (TPMs + 1)")+xlab(NULL)+
  #facet_wrap(~Reactor)+
  ggtitle("MMP1667 (archaellin flaB)")+
  #facet_grid(rows = vars(Enrichment))+
  theme_light()+theme_enrichment+theme(plot.title = element_text(face = "italic", vjust = -2, size= 7))
MMP1667_cts
ggsave(MMP1667_cts, file ="./Final_Plots/flaB_TPMs_Reps_Days3_6.pdf", width = 1.8, height = 2, units = "in", dpi = "retina")
write.csv(subset(Mmp_flaB_cts, Mmp_flaB_cts$Day %in% Day_Range), "./Source_Data/Figure_2F.csv")
