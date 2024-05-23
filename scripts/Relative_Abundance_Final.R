rm(list=ls())
getwd()
setwd("/Users/jake/Documents/R_analysis/R_Analysis/ENIGMA_Analysis/FBR_DESeq/Final_Code/")
getwd()

Packages<-c("dplyr", "FactoMineR", "factoextra", "ggplot2", "graphics", "grDevices", "matrixStats", "pheatmap", 
            "RColorBrewer","readr", "reshape2", "stats", "SummarizedExperiment", "tidyr","textshape",
            "readr", "tidyverse", "corrplot", "ggpubr", "pals", "ComplexHeatmap", "circlize", "paletteer", "growthrates")
#install.packages(Packages)
lapply(Packages, library, character.only =T)

##############################################################################################
'%notin%' <- Negate('%in%')
### define color scheme
syncom_col<-c("Dv" = "#D94801", "Mm" = "#BEBEBE","vulgaris" = "#D94801", "maripaludis" = "#BEBEBE" )

### set theme
theme_meta<-theme(plot.title = element_text(angle = 0, size = 8, face = "plain", hjust = .5, vjust = -2),
                  axis.text = element_text(angle = 0, size = 8),
                  axis.title= element_text(angle = 0, size = 8),
                  #axis.title.x = element_blank(),
                  strip.text = element_blank(),#,angle = 0, size = 8, color = "black", margin = margin(c(.07,.07,.07,.07), unit='cm')),
                  legend.box.background = element_rect(fill= "white", color = "grey", size =.25), 
                  legend.title = element_text(size=6),
                  legend.position = "bottom", legend.text = element_text(size = 6), # c(.09,.89)
                  legend.key.size = unit(.02, "in"), legend.key.width = unit(.02, "in"),
                  legend.justification = 'center', legend.title.align=0.5, legend.margin=margin(c(.05,.05,.05,.05), unit='cm'))

######## Analyse Percent GC from MultiQC. output
GC_perc<-read.table("./Final_Data_Files/mqc_fastqc_per_sequence_gc_content_plot_Percentages.txt",
                    header = TRUE,  check.names = FALSE)

#### format data for analysis and plotting
GC_df<-melt(GC_perc)
summary(GC_df)
colnames(GC_df)[2:3]<-c("Percent_GC", "percent_of_reads")

GC_df$Percent_GC<-as.numeric(GC_df$Percent_GC)
summary(GC_df)
GC_df<-cbind(GC_df, colsplit(GC_df$Sample, pattern = "_", names = c("Shortname", "Lane")))
GC_df<-cbind(GC_df, colsplit(GC_df$Shortname, pattern = "-", names = c("EPD", "Reactor", "Month", "Date", "Phase")))
GC_df$Day<-GC_df$Date-7
#### Calculate Percent Reads at GC Peaks
GC_peaks<-subset(GC_df, Percent_GC %in% c(30:40,60:70))
GC_peaks$Org<-if_else(condition = GC_peaks$Percent_GC <  41, true = "Mm", false = "Dv")
GC_peaks$Day<-as.factor(GC_peaks$Day)
ggplot(data = GC_peaks, aes(Percent_GC, percent_of_reads, fill = Org))+
  geom_point(shape = 21)+
  scale_fill_manual(values = syncom_col)

## Calculate mean and SD for relatvie abundance of reps  summarise_each: 
GC_peaks_summs<-GC_peaks %>% ### calulcates percentages by read abundance for each reactor, org, day, phase
  group_by(Org, Day, Phase, Reactor) %>%
  summarise(across(c(percent_of_reads),
                   list("mean" = ~ mean(.x, na.rm = TRUE), "sd" = ~ sd(.x, na.rm = TRUE), "sem"= ~ sd(.x)/sqrt(length(.x)),
                        "n" = ~ sum(!is.na(.x)))))

GC_Percent_abundance<-GC_peaks_summs %>%
  group_by(Day, Phase, Reactor) %>%
  dplyr::summarise(percent_of_reads_mean, Org) %>%
  dplyr::mutate("percent" = (percent_of_reads_mean/ sum(percent_of_reads_mean))*100)

GC_Percent_Stats<-GC_Percent_abundance %>% ### calulcates percentages by read abundance for each reactor, org, day, phase
  group_by(Org, Day, Phase) %>%
  summarise(across(c(percent_of_reads_mean, percent),
                   list("mean" = ~ mean(.x, na.rm = TRUE), "sd" = ~ sd(.x, na.rm = TRUE), "sem"= ~ sd(.x)/sqrt(length(.x)),
                        "n" = ~ sum(!is.na(.x)))))
GC_Percent_Stats <- transform(GC_Percent_Stats, lower=percent_mean-percent_sem, upper=percent_mean+percent_sem)
GC_Percent_Stats2 <- transform(GC_Percent_Stats, lower=percent_mean-percent_sd, upper=percent_mean+percent_sd)

ggplot()+
  geom_col(data = GC_Percent_Stats2, aes(y = percent_mean, x= Day, fill = Org), color ="black", size =.15, alpha=.95)+
  geom_errorbar(data=GC_Percent_Stats2[GC_Percent_Stats2$Org == "Mm",], aes(x = Day, ymax=upper,  ymin=percent_mean), width=0.15, size =.15)+
  #geom_point(data = GC_Percent_abundance,
  #           aes(x= Day, y = percent, shape = Reactor, fill = Org), position = position_dodge(width = .5), stroke = .15, size = .75)+
  scale_fill_manual(values = syncom_col)+xlab("Days (Fluidization)")+
  scale_shape_manual(values = c(21:23))+
  facet_wrap(~Phase)+scale_y_continuous(expand = c(.01,.01), name = str_wrap("Relative Abundance (% reads per GC content)", 21))+
  theme_bw()+theme_meta
ggsave("./Final_Plots/GC_content_Relative_Abundance_Mean_SD_revised.pdf", width = 2.75, height = 2.2, units = "in", dpi = "retina")
write.csv(GC_Percent_Stats2, "./Source_Data/FIgure_2A.csv", row.names = F )


