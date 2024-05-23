rm(list=ls())
getwd()
setwd("/Users/jake/Documents/R_analysis/R_Analysis/ENIGMA_Analysis/FBR_DESeq/Final_Code")
getwd()

##Install the required packages
Packages<-c("reshape2", "ggplot2", "RColorBrewer", "scales", "gridExtra", "dplyr", "ggpubr", 
            "tidyverse", "ggridges", "hrbrthemes", "vegan", "pheatmap", "plyr", "paletteer",
            "ggbeeswarm", "gghalves","ggpmisc" ,"ComplexHeatmap", "ggrepel", "grDevices", 'circlize')
#install.packages(Packages)
lapply(Packages, library, character.only =T)

##############################################################################################
'%notin%' <- Negate('%in%')
##########################################################################################################
###### Remove all lists ########
#rm(list=ls())
getwd()

####load Raw Ancestor DNA output get a list of Ancestral mutations
DNA_Anc_all_DvH <- cbind(read.delim("./Final_Data_Files/dvh_mutations_allsamples_attributes_3112020.txt", 
                      header = TRUE, na.strings = FALSE, check.names = TRUE), "Org" = "DvH")
DNA_Anc_all_Mmp <- cbind(read.delim("./Final_Data_Files/Mmp_mutations_allsamples_attributes_3112020.txt", 
                              header = TRUE, na.strings = FALSE, check.names = TRUE), "Org" = "Mmp")

DNA_Anc_all_Mmp<-DNA_Anc_all_Mmp[,-(length(DNA_Anc_all_Mmp)-1)] ### remove X column

table(colnames(DNA_Anc_all_DvH) == colnames(DNA_Anc_all_Mmp))

DNA_Anc_all<-rbind(DNA_Anc_all_DvH, DNA_Anc_all_Mmp)
#replace blanks with NA
DNA_Anc_all[DNA_Anc_all == ""] <- NA

##Convert frequencies into numeric
DNA_Anc_all$freq<-as.numeric(sub("%","", DNA_Anc_all$freq))

### Create a list of Ancestral mutations that can be used to subset EPD_FBR data
Anc_Mutations<-subset(DNA_Anc_all, DNA_Anc_all$experiment %in% "Ancestors")
Anc_position<-as.numeric(Anc_Mutations$position)

####load Raw FBR_DNA for DvH ####
DNA_FBR_DVH <- cbind(read.delim("./Final_Data_Files/Dvh_mutations_allsamples_attributes_3202020.txt", 
                         header = TRUE, na.strings = FALSE, check.names = TRUE, row.names = NULL), "Org" = "DvH")
colnames(DNA_FBR_DVH)<-colnames(DNA_Anc_all) #re-label colnmames

DNA_FBR_Mmp <- cbind(read.delim("./Final_Data_Files/Mmp_mutations_allsamples_attributes_3202020.txt", 
                          header = TRUE, na.strings = FALSE, check.names = TRUE, row.names = NULL), "Org" = "Mmp")
colnames(DNA_FBR_Mmp)<-colnames(DNA_Anc_all) #re-label colnmames
### check column names
table(colnames(DNA_FBR_DVH) == colnames(DNA_FBR_Mmp))

DNA_FBR<-rbind(DNA_FBR_DVH, DNA_FBR_Mmp)
#replace blanks with NA
DNA_FBR[DNA_FBR == ""] <- NA
### re-label experiment as FBR not syntrophy
DNA_FBR$experiment<-gsub("syntrophy", "FBR", DNA_FBR$experiment)

###Check colnames for symetry in order to rbind two dataframes
(colnames(DNA_Anc_all))==(colnames(DNA_FBR))
### Combine ancestral mutations with mutations in EPD-FBR treatments
DNA_All<-rbind(DNA_Anc_all,DNA_FBR)
### re-factor datafraome
DNA_All <- as.data.frame(lapply(DNA_All, function (x) if (is.factor(x)) factor(x) else x))

####Working with Dataframe of just FBR mutations ####
### re-factor datafraome
DNA_FBR <- as.data.frame(lapply(DNA_FBR, function (x) if (is.factor(x)) factor(x) else x))

#### Adding categorical variables
FBR_split<-colsplit(DNA_FBR$sample, '-', names =  c('Strain',"Rep","Month","Day", "Phase"))
FBR_split<-FBR_split %>% unite(Date, c("Month", "Day"))
### Add Day column
FBR_split$Day<-NA
for(i in FBR_split$Date){
  FBR_split$Day[FBR_split$Date %in% c("5_8")] <- 1
  FBR_split$Day[FBR_split$Date %in% c("5_9")] <- 2
  FBR_split$Day[FBR_split$Date %in% c("5_10")] <- 3
  FBR_split$Day[FBR_split$Date %in% c("5_11")] <- 4
  FBR_split$Day[FBR_split$Date %in% c("5_12")] <- 5
  FBR_split$Day[FBR_split$Date %in% c("5_13")] <- 6
}

### re-combine dataframe and new labels
DNA_FBR<-cbind(FBR_split, DNA_FBR)
### replace Phase definition "P" with Planktonic and S with Attached
DNA_FBR$Phase<-sub("P","Planktonic", DNA_FBR$Phase)
DNA_FBR$Phase<-sub("S","Attached", DNA_FBR$Phase)
##Convert frequencies into numeric
DNA_FBR$freq<-as.numeric(sub("%","", DNA_FBR$freq))
DNA_FBR$freq1<-as.numeric(sub("%","", DNA_FBR$freq1))
DNA_FBR$freq2<-as.numeric(sub("%","", DNA_FBR$freq2))
###Removing Anc_Mutations from FBR Mutations by position
FBR_Mutations<-subset(DNA_FBR, DNA_FBR$position %notin% Anc_position)

############################################
### Combine strain and phase ###
FBR_Mutations$Phase_Rep<-as.factor(paste(FBR_Mutations$Phase, FBR_Mutations$Rep, sep = "_"))
FBR_Mutations$Phase_Day<-as.factor(paste(FBR_Mutations$Phase, FBR_Mutations$Day, sep = "_"))
FBR_Mutations$Phase_Day_Rep<-as.factor(paste(FBR_Mutations$Phase, FBR_Mutations$Day,  FBR_Mutations$Rep, sep = "_"))
FBR_Mutations$Org_Phase_Day_Rep<-as.factor(paste(FBR_Mutations$Org, FBR_Mutations$Phase, FBR_Mutations$Day,  FBR_Mutations$Rep, sep = "_"))
FBR_Mutations$Org_Phase<-as.factor(paste(FBR_Mutations$Org, FBR_Mutations$Phase, sep = "_"))
FBR_Mutations$Org_Phase_Day<-as.factor(paste(FBR_Mutations$Org, FBR_Mutations$Phase, FBR_Mutations$Day, sep = "_"))

FBR_Mutations$mean <- (round(ave(FBR_Mutations$freq, as.factor(interaction(FBR_Mutations$Day,FBR_Mutations$Phase)), FUN=mean), digits = 2))
FBR_Mutations$Ghostpointy<-45
FBR_Mutations$Ghostpointy2<-10
FBR_Mutations$Ghostpointx<-60
FBR_Mutations$freq3 <- 100-FBR_Mutations$freq
#convert remove underscore in gene_id
FBR_Mutations$gene_id <- ifelse(is.na(FBR_Mutations$gene_id), FBR_Mutations$variant_id, FBR_Mutations$gene_id)
FBR_Mutations$gene_id<-gsub("_","", FBR_Mutations$gene_id)

###DvH mutation that may be outliers
spikes_DvH<-c("Chr-IG-495942","Chr-IG-1623932", "Chr-DVU3172-3329290", "Chr-DVU2395-2498730", "Chr-DVU1632-1716114", 
          "Chr-DVU1571-1653920", "Chr-DVU1214-1305196", "Chr-DVU0030-37335")
###Mmp mutation that may be outliers
spikes_Mmp<-c("MMP1325","MMP0815", "MMP0762", "MMP0111", "MMP0117")

spikes<-c(spikes_DvH, spikes_Mmp)
#######
FBR_Mutations2<-subset(FBR_Mutations, FBR_Mutations$variant_id %notin% spikes)
FBR_Mutations2<-subset(FBR_Mutations2, FBR_Mutations2$gene_id %notin% spikes)### FBR_mutations2<- removes gene variants that flash
FBR_Mutations2$gene_id <- ifelse(is.na(FBR_Mutations2$gene_id), FBR_Mutations2$variant_id, FBR_Mutations2$gene_id)
#write.csv(x = FBR_Mutations2, file = "~/Documents/ENIGMA/FBR_Manuscript/Manuscript/Supplemental/FBR_Mutations2.csv")

## Calculate mean frequency of reps dplyr::summarise and grouping by reactor ##
FBR_mutations_mean<-group_by(FBR_Mutations2, Org, Day, Phase, gene_id, variant_id) %>%
  dplyr::summarise(
    count = as.factor(n()),
    mean_freq = mean(freq),
    sd_freq=sd(freq)
  )
#write.csv(FBR_mutations_mean, "~/Documents/ENIGMA/FBR_Manuscript/Manuscript/Supplemental/mean_mut_frequencies.csv")

#### Calculate Averaged for Ridgeline plots 
avg<-FBR_Mutations2 %>%
  group_by(Phase, Org) %>%
  dplyr::summarise(across(freq, list("mean" = mean, "sd" = sd, "n" = length))) 
avg$nvalue<-paste("n = ", avg$freq_n, sep ="")

avg1<-subset(FBR_Mutations2, FBR_Mutations2$freq != 100) %>%
  group_by(Phase, Org) %>%
  dplyr::summarise(across(freq, list("mean" = mean, "sd" = sd, "n" = length))) 
avg1$nvalue<-paste("n = ", avg1$freq_n, sep ="")

avg2<-FBR_Mutations2 %>%
  group_by(Phase, Org, Day) %>%
  dplyr::summarise(across(freq, list("mean" = mean, "sd" = sd, "n" = length))) 
avg2$nvalue<-paste("n = ", avg2$freq_n, sep ="")

avg3<-subset(FBR_Mutations2, FBR_Mutations2$freq != 100) %>%
  group_by(Phase, Org, Day) %>%
  dplyr::summarise(across(freq, list("mean" = mean, "sd" = sd, "n" = length))) 
avg3$nvalue<-paste("n = ", avg3$freq_n, sep ="")

avg4<-subset(FBR_Mutations2, FBR_Mutations2$freq == 100) %>%
  group_by(Phase, Org, Day) %>%
  dplyr::summarise(across(freq, list("mean" = mean, "sd" = sd, "n" = length))) 
avg4$nvalue<-paste("n = ", avg4$freq_n, sep ="")


FBR_Mutations2$Phase<-factor(FBR_Mutations2$Phase, levels = c("Planktonic", "Attached"))
FBR_Mutations2 <- as.data.frame(lapply(FBR_Mutations2, function (x) if (is.factor(x)) factor(x) else x))
colnames(FBR_Mutations2)

summary(FBR_Mutations2$Phase_Day_Rep)
summary(FBR_Mutations2$Org_Phase_Day_Rep)

######## Reshape FBR_mutations into matrix
Var_Freq_check<-FBR_Mutations2[,c("Phase_Day_Rep", "Phase", "freq", "gene_id")]

Var_Freq_check$Phase<-as.factor(Var_Freq_check$Phase);Var_Freq_check$gene_id<-as.factor(Var_Freq_check$gene_id)

Var_Freq_wide_check<-reshape2::dcast(Var_Freq_check, Phase_Day_Rep + Phase ~ gene_id, value.var = "freq" )

### replace NAs with 0
Var_Freq_wide_check <- Var_Freq_wide_check %>% replace(is.na(.), 0)
#### spotcheck mutations
Var_Freq_wide_Mat_check<-Var_Freq_wide_check[,3:length(Var_Freq_wide_check)]
rownames(Var_Freq_wide_Mat_check)<-Var_Freq_wide_check[,1]

#quick check of the the mutation landscape
col_Freq = colorRamp2(c(0, 50, 100), c("Blue", "white", "Red"), space = "RGB", )
ComplexHeatmap::Heatmap(t(Var_Freq_wide_Mat_check), cluster_columns = F, col = col_Freq,
                        row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 8))
dev.off()
######## 

############ Data Vis #########
histo_theme<-theme(strip.background    = element_blank(),
                   axis.text.x         = element_text(size = 6.5),
                   axis.text.y         = element_text(size = 6.5),
                   axis.title.y        = element_text(size = 7),
                   axis.title.x        = element_text(size = 7),
                   axis.line.x         = element_blank(),
                   axis.ticks.x        = element_line(colour = "gray90"),
                   #axis.ticks.length.x = unit(30, "points"),
                   strip.text   = element_text(angle = 0, size = 7, color = "Black", vjust = -1.5),
                   panel.grid.major.x  = element_line(colour = "gray90"),
                   legend.position = "none",
                   legend.background = element_rect(color ="black", fill = "white", linewidth =.2),
                   legend.title = element_blank(), 
                   legend.text=element_text(size=5), 
                   legend.key.size = unit(.2, "cm"),
                   legend.margin = margin(-5,1,1,1, unit = "points"))
# basic example
DvH_muts2<-subset(FBR_Mutations2, FBR_Mutations2$Org %in% "DvH" & FBR_Mutations2$freq != 100)
Mmp_muts2<-subset(FBR_Mutations2, FBR_Mutations2$Org %in% "Mmp" & FBR_Mutations2$freq != 100) # & FBR_Mutations2$freq != 100)

nvalues_DvH<-data.frame(summary(DvH_muts2$Org_Phase))
nvalues_Mmp<-data.frame(summary(Mmp_muts2$Org_Phase))
nvalues2<-data.frame(summary(FBR_Mutations2$Org_Phase_Day))
##### DvH Distributions
freq_day_density<-ggplot(data = NULL) +
  geom_density_ridges(data =DvH_muts2,
                      aes(x = freq, y = factor(Day), fill = Phase),
                      scale =3, alpha = .9, color ="black")+
  scale_x_continuous(breaks = c(40,50,60,70,80,90,100), limits = c(45,102), 
                     expand = c(.05,.05),name = "Mutation Frequency") +
  scale_y_discrete( breaks = c(1:6), expand = c(.07,-.1), name = "Days", position = "right", limits = rev)+
  scale_fill_manual(values = c("Planktonic" =  "#1f78b4", "Attached" = "#a6cee3"), 
                    name = "Phase", label = c("Planktonic" = "Planktonic", "Attached" = "Sediment"))+
  facet_wrap(~factor(Phase, levels = c("Planktonic", "Attached")), nrow = 2, strip.position = "left")+
  theme_bw()+histo_theme+theme(legend.position = "none")
freq_day_density

freq_halves <- ggplot(data =DvH_muts2, aes(x = Phase, y = freq, fill = factor(Phase, levels =rev(c("Planktonic", "Attached"))))) +
  geom_half_point(alpha = .9, color ="black", size = 0.05,
                  transformation = position_quasirandom(width = 0.15),
                  side = "l") +
  geom_half_violin(aes(fill = Phase), side = "r", color = "black", linewidth =.25) + 
  geom_segment(data = subset(avg1, avg1$Org %in% "DvH" ), 
               aes(y=freq_mean, yend=freq_mean,
                   x = rev(c(2,1)), xend = rev(c(2.4, 1.4))),
               color="black", size = .25, linetype =2)+
  scale_fill_manual(values = c("Planktonic" =  "#1f78b4", "Attached" = "#a6cee3"), 
                    name = "Phase", label = c("Planktonic" = "Planktonic", "Attached" = "Sediment"))+
  scale_x_discrete(limits= rev, position = "top")+
  scale_y_continuous(limits =c(45,102), breaks = c(40,50,60,70,80,90,100), 
                     expand = c(.05,.05),name = "Mutation Frequency" )+
  guides(color = "none", fill = "none") + xlab(NULL)+
  coord_flip()+
  theme_bw()+histo_theme+theme(axis.text.y = element_text(size = 6.5, angle = 90))
freq_halves

density_plt_DvH<-ggarrange(freq_halves,freq_day_density, ncol = 2, heights = c(2), widths = c(2,2.5), 
                           labels = c("B)", "C)"),vjust =  0,font.label = list(size = 8, color = "black", face = "plain"))
density_plt_DvH<-annotate_figure(density_plt_DvH, top = text_grob("Desulfovibrio vulgaris", color = "black", 
                                                                  face = "italic", size = 7, vjust = 1.7, hjust = 1.6,))
density_plt_DvH

ggsave(density_plt_DvH, file ="./Final_Plots/density_plt_DvH_freq_horizontal_revised.pdf", width = 3.5, height = 2.5, units = "in", dpi = "retina")
dev.off()
write.csv(DvH_muts2, "./Source_Data/Figure_2BC.csv")
### Calculate p. values
#Overall mean distribution using Wilcox text
#Parametric tests are those that make assumptions about the parameters of the population distribution from which the sample is drawn. 
#This is often the assumption that the population data are normally distributed. Non-parametric tests are “distribution-free” and, as such, 
#can be used for non-Normal variables.
DvH_muts2$Phase
t_test_Phase<-t.test(subset(DvH_muts2, Phase == "Attached")$freq,
                          subset(DvH_muts2, Phase == "Planktonic")$freq, alternative= "two.sided")
t_test_Phase$p.value

pvalue_phase_total<-as.data.frame(list(t_test_Phase$estimate[1],t_test_Phase$estimate[2],
                                       t_test_Phase$p.value), col.names= c("mean_x", "mean_y", "pvalue"), row.names = "EPD_Phase")
#### T-test mean distribution by of each day compared to day 1+
Day1pvalue_sed<-NULL
for (i in unique(DvH_muts2$Day)) {
  print(i)
  temp_pvalue<-t.test(subset(DvH_muts2, Phase %in% "Attached" & Day %in% "1")$freq,
         subset(DvH_muts2, Phase %in% "Attached" & Day %in% i)$freq, alternative= "two.sided")
  Day1pvalue_temp<-as.data.frame(list(temp_pvalue$estimate[1],temp_pvalue$estimate[2],
                                 temp_pvalue$p.value), col.names= c("mean_Day_1","mean_Day_y" , "pvalue"), 
                                 row.names = paste("Attached", i, sep = " "))
  Day1pvalue_temp$Phase<-"Attached"
  Day1pvalue_sed<-rbind(Day1pvalue_sed, Day1pvalue_temp)
}

Day1pvalue_plank<-NULL
for (i in unique(DvH_muts2$Day)) {
  print(i)
  temp_pvalue<-t.test(subset(DvH_muts2, Phase %in% "Planktonic" & Day %in% "1")$freq,
                      subset(DvH_muts2, Phase %in% "Planktonic" & Day %in% i)$freq, alternative= "two.sided")
  Day1pvalue_temp<-as.data.frame(list(temp_pvalue$estimate[1],temp_pvalue$estimate[2],
                                      temp_pvalue$p.value), col.names= c("mean_Day_1","mean_Day_y" , "pvalue"), 
                                 row.names = paste("Planktonic", i, sep = " "))
  Day1pvalue_temp$Phase<-"Planktonic"
  Day1pvalue_plank<-rbind(Day1pvalue_plank, Day1pvalue_temp)
}
Day1pvalue<-rbind(Day1pvalue_sed, Day1pvalue_plank)
Day1pvalue$pvalue_round<-round(Day1pvalue$pvalue, 5)
Day1pvalue
##########
#### T-test mean distribution by of each day between Attached and planktonic
PhaseDaypvalue<-NULL
for (i in unique(DvH_muts2$Day)) {
  print(i)
  temp_pvalue<-t.test(subset(DvH_muts2, Phase %in% "Planktonic" & Day %in% i)$freq,
                      subset(DvH_muts2, Phase %in% "Attached" & Day %in% i)$freq, alternative= "two.sided")
  pvalue_temp<-as.data.frame(list(temp_pvalue$estimate[1],temp_pvalue$estimate[2],
                                      temp_pvalue$p.value), col.names= c("mean_Day_P","mean_Day_S" , "pvalue"), 
                                 row.names = paste("Day_", i, sep = ""))
  pvalue_temp$Day<-i
  PhaseDaypvalue<-rbind(PhaseDaypvalue, pvalue_temp)
}
PhaseDaypvalue$pvalue_round<-round(PhaseDaypvalue$pvalue, 5)
PhaseDaypvalue
#########
#### Distribution analysis for Mm
freq_day_density<-ggplot(data = NULL) +
  geom_density_ridges(data =Mmp_muts2,
                      aes(x = freq, y = factor(Day), fill = Phase),
                      scale =3, alpha = .9, color ="black", size = 0.25)+
  scale_x_continuous(breaks = c(40,50,60,70,80,90,100), limits = c(45,102), 
                     expand = c(.05,.05),name = "Mutation Frequency") +
  scale_y_discrete( breaks = c(1:6), expand = c(.07,-.1), name = "Days", position = "right", limits = rev)+
  scale_fill_manual(values = c("Planktonic" =  "#1f78b4", "Attached" = "#a6cee3"), 
                    name = "Phase", label = c("Planktonic" = "Planktonic", "Attached" = "Sediment"))+
  facet_wrap(~factor(Phase, levels = c("Planktonic", "Attached")), nrow = 2, strip.position = "left")+
  theme_bw()+histo_theme+theme(legend.position = "none")
freq_day_density

freq_halves <- ggplot(data =Mmp_muts2, aes(x = Phase, y = freq, fill = factor(Phase, levels =rev(c("Planktonic", "Attached"))))) +
  geom_half_point(alpha = .9, color ="black", size = 0.05,
                  transformation = position_quasirandom(width = 0.15),
                  side = "l") +
  geom_half_violin(aes(fill = Phase), side = "r", color = "black", size =.25) + 
  geom_segment(data = subset(avg1, avg1$Org %in% "Mmp" ), 
               aes(y=freq_mean, yend=freq_mean,
                   x = rev(c(2,1)), xend = rev(c(2.4, 1.4))),
               color="black", size = .25, linetype =2)+
  scale_fill_manual(values = c("Planktonic" =  "#1f78b4", "Attached" = "#a6cee3"), 
                    name = "Phase", label = c("Planktonic" = "Planktonic", "Attached" = "Sediment"))+
  scale_x_discrete(limits= rev, position = "top")+
  scale_y_continuous(limits =c(45,102), breaks = c(40,50,60,70,80,90,100), 
                     expand = c(.05,.05),name = "Mutation Frequency" )+
  #scale_y_continuous(breaks = rev(c(40,50,60,70,80,90,100)), limits = c(102,45),
  #                   expand = c(.05,.05),name = "Mutation Frequency" ) +
  guides(color = "none", fill = "none") + xlab(NULL)+
  coord_flip()+
  theme_bw()+histo_theme+theme(axis.text.y = element_text(size = 6.5, angle = 90))
freq_halves

density_plt_Mmp<-ggarrange(freq_halves,freq_day_density, ncol = 2, heights = c(2), widths = c(2,2.5), 
                           labels = c("B)", "C)"),vjust =  0,font.label = list(size = 8, color = "black", face = "plain"))
density_plt_Mmp<-annotate_figure(density_plt_Mmp, top = text_grob("Methanococcus maripaludis", color = "black", 
                                                                  face = "italic", size = 7, vjust = 1.7, hjust = 1.6,))
density_plt_Mmp

#### T-test mean distribution by of each day compared to day 1+
# Day1pvalue_sed<-NULL
# for (i in unique(Mmp_muts2$Day)) {
#   print(i)
#   temp_pvalue<-t.test(subset(Mmp_muts2, Phase %in% "Attached" & Day %in% "1")$freq,
#                       subset(Mmp_muts2, Phase %in% "Attached" & Day %in% i)$freq, alternative= "two.sided")
#   Day1pvalue_temp<-as.data.frame(list(temp_pvalue$estimate[1],temp_pvalue$estimate[2],
#                                       temp_pvalue$p.value), col.names= c("mean_Day_1","mean_Day_y" , "pvalue"), 
#                                  row.names = paste("Attached", i, sep = " "))
#   Day1pvalue_temp$Phase<-"Attached"
#   Day1pvalue_sed<-rbind(Day1pvalue_sed, Day1pvalue_temp)
# }
# 
# Day1pvalue_plank<-NULL
# for (i in unique(Mmp_muts2$Day)) {
#   print(i)
#   temp_pvalue<-t.test(subset(Mmp_muts2, Phase %in% "Planktonic" & Day %in% "1")$freq,
#                       subset(Mmp_muts2, Phase %in% "Planktonic" & Day %in% i)$freq, alternative= "two.sided")
#   Day1pvalue_temp<-as.data.frame(list(temp_pvalue$estimate[1],temp_pvalue$estimate[2],
#                                       temp_pvalue$p.value), col.names= c("mean_Day_1","mean_Day_y" , "pvalue"), 
#                                  row.names = paste("Planktonic", i, sep = " "))
#   Day1pvalue_temp$Phase<-"Planktonic"
#   Day1pvalue_plank<-rbind(Day1pvalue_plank, Day1pvalue_temp)
# }
# Day1pvalue<-rbind(Day1pvalue_sed, Day1pvalue_plank)
# Day1pvalue$pvalue_round<-round(Day1pvalue$pvalue, 5)
# Day1pvalue
# ##########
# #### T-test mean distribution by of each day between Attached and planktonic
# PhaseDaypvalue<-NULL
# for (i in unique(Mmp_muts2$Day)) {
#   print(i)
#   temp_pvalue<-t.test(subset(Mmp_muts2, Phase %in% "Planktonic" & Day %in% i)$freq,
#                       subset(Mmp_muts2, Phase %in% "Attached" & Day %in% i)$freq, alternative= "two.sided")
#   pvalue_temp<-as.data.frame(list(temp_pvalue$estimate[1],temp_pvalue$estimate[2],
#                                   temp_pvalue$p.value), col.names= c("mean_Day_P","mean_Day_S" , "pvalue"), 
#                              row.names = paste("Day_", i, sep = ""))
#   pvalue_temp$Day<-i
#   PhaseDaypvalue<-rbind(PhaseDaypvalue, pvalue_temp)
# }
# PhaseDaypvalue$pvalue_round<-round(PhaseDaypvalue$pvalue, 5)
# PhaseDaypvalue

####
#### Line plots for mutations
##### line plots for selected DvH Genes
############
lineplot_theme<-theme(#axis.text = element_text(angle = 0, size=14),
  plot.title=element_text(size =5, vjust = -2),
  axis.title = element_text(size = 7), 
  legend.box.background = element_rect(color = "grey"), 
  legend.key.size = unit(.01, "cm"), legend.key.width = unit(.25, "cm"),
  axis.text.x = element_text(size=7, color = "Black", angle=0),
  axis.text.y = element_text(size=7, color = "Black"),
  strip.background = element_blank(),
  strip.text  = element_blank(),#color = "black", fill = "white", size = 1.0, linetype = "solid"),
  #strip.text.x = element_text(size=8, color = "Black", angle=0),
  #strip.text.y = element_text(size=8, color = "Black", angle=0),
  legend.title = element_text(size = 6, color="Black"),legend.text = element_text(size = 5, color="Black"),
  legend.position = c("none"))
######  Replace Nan with 0 frequency (Use Var_Freq_wide)
###################################
mutation_freq<-melt(Var_Freq_wide_check) ## or use Var_Freq_wide
mutation_freq<-mutation_freq[,-2]
mutation_freq<-cbind(mutation_freq, colsplit(mutation_freq$Phase_Day_Rep, pattern = "_", names = c("Phase", "Day", "Rep")))
colnames(mutation_freq)[c(2,3)]<-c("gene_id", "freq")
mutation_freq$gene_id<-as.character(mutation_freq$gene_id)
mutation_freq$Org <- ifelse(startsWith(mutation_freq$gene_id, "MM"), "Mmp", "DvH")

## Calculate mean frequency of reps dplyr::summarise and grouping by reactor ##
mutation_freq_mean<-group_by(mutation_freq, Org, Day, Phase, gene_id) %>%
  dplyr::summarise(
    count = as.factor(n()),
    mean_freq = mean(freq),
    sd_freq=sd(freq)
  )

#plot both Mm and Dv
mutation_freq_mean
mutation_freq_mean$Time_Days<-mutation_freq_mean$Day

mutation_freq_mean$Time_Days<-if_else(mutation_freq_mean$Time_Days == 1, true = 47.5, false = mutation_freq_mean$Time_Days)
mutation_freq_mean$Time_Days<-if_else(mutation_freq_mean$Time_Days == 2, true = 71, false = mutation_freq_mean$Time_Days)
mutation_freq_mean$Time_Days<-if_else(mutation_freq_mean$Time_Days == 3, true = 95, false = mutation_freq_mean$Time_Days)
mutation_freq_mean$Time_Days<-if_else(mutation_freq_mean$Time_Days == 4, true = 118, false = mutation_freq_mean$Time_Days)
mutation_freq_mean$Time_Days<-if_else(mutation_freq_mean$Time_Days == 5, true = 138, false = mutation_freq_mean$Time_Days)
mutation_freq_mean$Time_Days<-if_else(mutation_freq_mean$Time_Days == 6, true = 166, false = mutation_freq_mean$Time_Days)
#mutation_freq_mean$Time_Days<-as.numeric(mutation_freq_mean$Time_Days)
mutation_freq_mean$Time_Days<-(mutation_freq_mean$Time_Days/24)-1

#DvH_mutation_freq_mean<-subset(mutation_freq_mean, mutation_freq_mean$Org %in% "DvH")
SynCom_mutation_freq_mean<-mutation_freq_mean
SynCom_mutation_freq_mean_Planktonic_mean<-subset(SynCom_mutation_freq_mean, SynCom_mutation_freq_mean$Phase %notin% "Attached")# & FBR_Mutations2$gene_id %notin% "DVU1283")
SynCom_mutation_freq_mean_Sediment_mean<-subset(SynCom_mutation_freq_mean, SynCom_mutation_freq_mean$Phase %in% "Attached")# & FBR_Mutations2$gene_id %notin% "DVU1283")
SynCom_sediment_genes_mean<-unique(factor(SynCom_mutation_freq_mean_Sediment_mean$gene_id))
SynCom_regression_df_mean<-NULL

for (i in SynCom_sediment_genes_mean) {
  print (i)
  temp_model_S <- lm(mean_freq ~ Time_Days, data = SynCom_mutation_freq_mean_Sediment_mean[SynCom_mutation_freq_mean_Sediment_mean$gene_id == i,])
  temp_model_P <- lm(mean_freq ~ Time_Days, data = SynCom_mutation_freq_mean_Planktonic_mean[SynCom_mutation_freq_mean_Planktonic_mean$gene_id == i,])
  temp_lm_S<-cbind(data.frame("gene_ID" = i), t(temp_model_S$coefficients), data.frame("R2" = summary(temp_model_S)$r.squared),"Pvalue" = summary(temp_model_S)$coefficients[,4][2], data.frame("Phase" = "Attached"))
  temp_lm_P<-cbind(data.frame("gene_ID" = i), t(temp_model_P$coefficients), data.frame("R2" = summary(temp_model_P)$r.squared), "Pvalue" = summary(temp_model_P)$coefficients[,4][2], data.frame("Phase" = "Planktonic"))
  SynCom_regression_df_mean<-rbind(SynCom_regression_df_mean, temp_lm_S, temp_lm_P)
  
}
colnames(SynCom_regression_df_mean)[3]<-"Slope"
SynCom_regression_df_mean$Pvalue<-round(SynCom_regression_df_mean$Pvalue,4)
SynCom_regression_df_mean
#reorder regression_Df_mean
SynCom_regression_df_mean<-SynCom_regression_df_mean[order(SynCom_regression_df_mean$R2, decreasing = T),]
SynCom_mutation_freq_mean$gene_id<-factor(SynCom_mutation_freq_mean$gene_id, levels = unique(SynCom_regression_df_mean$gene_ID))


SynCom_mutation_freq_mean_Planktonic_mean<-subset(SynCom_mutation_freq_mean, SynCom_mutation_freq_mean$Phase %notin% "Attached")# & FBR_Mutations2$gene_id %notin% "DVU1283")
SynCom_mutation_freq_mean_Sediment_mean<-subset(SynCom_mutation_freq_mean, SynCom_mutation_freq_mean$Phase %in% "Attached")# & FBR_Mutations2$gene_id %notin% "DVU1283")

SynCom_regression_df_mean
SynCom_regression_df_mean$Org<-if_else(startsWith(SynCom_regression_df_mean$gene_ID, "M") == TRUE, true = "Mm", false = "Dv")
SynCom_regression_df_mean$Org<-as.factor(SynCom_regression_df_mean$Org)

DvH_Reg_rev_plt_mean<-ggplot(data = SynCom_regression_df_mean[SynCom_regression_df_mean$Org %in% "Dv",],
                            aes(y = R2, x = Slope))+
  geom_hline(yintercept = .5, linetype = 2, color = "black", size = .25)+
  geom_vline(xintercept = c(-2.5, 2.5), linetype = 2, color = "black", size = .25)+
  geom_point(aes(fill = Phase), size = 2.5, shape = 21)+ 
  geom_label_repel(data = na.omit(SynCom_regression_df_mean[SynCom_regression_df_mean$R2 > .5 & SynCom_regression_df_mean$Org == "Dv" | SynCom_regression_df_mean$gene_ID == "DVU1283" ,]),
                   aes(label = gene_ID),label.size = .2, size = 1.5, label.padding = .1, ## https://ggrepel.slowkow.com/articles/examples.html
                   min.segment.length = 0, ##Use min.segment.length = 0 to draw all line segments, no matter how short they are
                   box.padding   = .5,max.overlaps = Inf,
                   point.padding = .5,
                   segment.color = 'grey50')+
  scale_fill_manual(values=c("Planktonic" ="#1f78b4","Attached" ="#a6cee3"), name = 'Phase')+
  scale_x_continuous(name = "Slope (Percent Frequency/Day)")+
  theme_light()+ guides(fill=guide_legend(nrow=2))+
  lineplot_theme
DvH_Reg_rev_plt_mean

p_regression_Dv_rev<-ggplot()+ 
  stat_poly_line(data = SynCom_mutation_freq_mean_Planktonic_mean[SynCom_mutation_freq_mean_Planktonic_mean$Org %in% "DvH",], 
                 aes(x= Time_Days, y=mean_freq, fill=factor(Phase)), color = "Black",
                 se=T, size=.3, linetype =3, level = .95, alpha = .85, show.legend = T) +
  stat_poly_line(data = SynCom_mutation_freq_mean_Sediment_mean[SynCom_mutation_freq_mean_Sediment_mean$Org %in% "DvH",], 
                 aes(x= Time_Days, y=mean_freq, fill=factor(Phase)), color = "Black",
                 se=T, size=.3, linetype =2, level = .95, alpha = .85) +
  stat_poly_eq(data = SynCom_mutation_freq_mean_Planktonic_mean[SynCom_mutation_freq_mean_Planktonic_mean$Org %in% "DvH",],
               aes(x= Time_Days, y=mean_freq, fill=factor(Phase),
                   label = paste(after_stat(rr.label), sep = "*\", \"*"), fill = Phase), 
               label.x = c(.07), label.y = .98, color = "white",
               geom = "label_npc", label.size = 0, label.padding = unit(0.04, "cm"),
               size =1.3, alpha =.95) +
  stat_poly_eq(data = SynCom_mutation_freq_mean_Sediment_mean[SynCom_mutation_freq_mean_Sediment_mean$Org %in% "DvH",],
               aes(x= Time_Days, y=mean_freq, fill=factor(Phase),
                   label = paste(after_stat(rr.label), sep = "*\", \"*"), fill = Phase), 
               label.x = c(.95), label.y = .98, color = "Black",
               geom = "label_npc", label.size = 0, label.padding = unit(0.04, "cm"),
               size =1.3, alpha =.95) +
  geom_point(data = SynCom_mutation_freq_mean[SynCom_mutation_freq_mean$Org %in% "DvH",], 
             aes(x= Time_Days, y=mean_freq, fill=factor(Phase), group = Phase),
             color ="black", size =1, shape = 21, stroke = .25)+
  scale_fill_manual(values=c("Planktonic" ="#1f78b4","Attached" ="#a6cee3"), name = 'Phase', guide = NULL)+
  scale_y_continuous(limits = c(-15,120), breaks = seq(0,100, by =20), name = "Mutation Frequency")+
  scale_x_continuous(limits = c(.5, 6.2), breaks = c(1:6.), name = "Day (Fluidization)", expand = c(0,0)) +
  facet_wrap(~gene_id, ncol =5, scales = "fixed")+
  theme_light()+ guides(fill=guide_legend(nrow=2))+#ggtitle("Dv Mutations")+
  lineplot_theme+theme(strip.background = element_rect(fill = "#D94801", color = "#D94801", size = .5, linetype = "solid"),
                     strip.text.x = element_text(size=5, color = "White", angle=0, margin = margin(c(.05,.05,.05,.05), unit='cm')))
p_regression_Dv_rev

Dv_revision_supplementary_LN_regression<-ggarrange(DvH_Reg_rev_plt_mean, p_regression_Dv_rev, ncol = 2, widths = c(2.75,4.25))
ggsave(plot = Dv_revision_supplementary_LN_regression, filename = "./Final_Plots/Dv_supplementary_LN_regression.pdf", width = 6.5, height = 3,
       units = "in", dpi = "retina")

FBR_Mutations_galU_Revision<-subset(SynCom_mutation_freq_mean, SynCom_mutation_freq_mean$gene_id %in% "DVU1283")
FBR_Mutations_galU_Sediment_Revision<-subset(SynCom_mutation_freq_mean, SynCom_mutation_freq_mean$Phase %in% "Attached" & SynCom_mutation_freq_mean$gene_id %in% "DVU1283")
FBR_Mutations_galU_Planktonic_Revision<-subset(SynCom_mutation_freq_mean, SynCom_mutation_freq_mean$Phase %notin% "Attached" & SynCom_mutation_freq_mean$gene_id %in% "DVU1283")

p_GalU_revision2<-ggplot()+ 
  stat_poly_line(data = FBR_Mutations_galU_Planktonic_Revision, 
                 aes(x= Time_Days, y=mean_freq, fill=factor(Phase)), color = "Black",
                 se=T, size=.25, linetype =3, level = .95, alpha = .85, show.legend = T) +
  stat_poly_line(data = FBR_Mutations_galU_Sediment_Revision, 
                 aes(x= Time_Days, y=mean_freq, fill=factor(Phase)), color = "Black",
                 se=T, size=.25, linetype =2, level = .95, alpha = .85) +
  stat_poly_eq(data = FBR_Mutations_galU_Planktonic_Revision,
               aes(x= Time_Days, y=mean_freq, fill=factor(Phase),
                   label = paste(after_stat(rr.label), sep = "*\", \"*"), fill = Phase), 
               label.x = c(.1), label.y = .98, color = "White",
               geom = "label_npc", label.size = 0, label.padding = unit(0.04, "cm"),
               size =2.25, alpha =.85) +
  stat_poly_eq(data = FBR_Mutations_galU_Sediment_Revision,
               aes(x= Time_Days, y=mean_freq, fill=factor(Phase),
                   label = paste(after_stat(rr.label), sep = "*\", \"*"), fill = Phase), 
               label.x = c(.75), label.y = .98, color = "Black",
               geom = "label_npc", label.size = 0, label.padding = unit(0.04, "cm"),
               size =2.25, alpha =.85) +
  geom_point(data = FBR_Mutations_galU_Revision, 
             aes(x= Time_Days, y=mean_freq, fill=factor(Phase), color = Phase, group = Phase),
             color ="black", size =1.5, shape = 21, stroke = .25)+
  scale_fill_manual(values=c("Planktonic" ="#1f78b4","Attached" ="#a6cee3"), name = 'Phase', guide = NULL)+
  scale_y_continuous(limits = c(0,115), breaks = seq(0,100, by =20), name = "Mutation Frequency")+
  scale_x_continuous( breaks = c(1:6.), name = "Day (Fluidization)", expand = c(.025, .025)) +
  theme_light()+ guides(fill=guide_legend(nrow=2))+
  lineplot_theme
p_GalU_revision2
ggsave("./Final_Plots/galU_Regression_plts_revision.pdf", width = 2, height = 2,
       units = "in", dpi = "retina")
write.csv(FBR_Mutations_galU_Revision, "./Source_Data/Figure_2G.csv")
#######
######## NDMS analysis
######## Reshape FBR_mutations into matrix
Var_Freq<-FBR_Mutations2[,c("Phase_Day_Rep", "Phase", "freq", "gene_id")]

Var_Freq$Phase<-as.factor(Var_Freq$Phase);Var_Freq$gene_id<-as.factor(Var_Freq$gene_id)

Var_Freq_wide<-reshape2::dcast(Var_Freq, Phase_Day_Rep + Phase ~ gene_id, value.var = "freq" )

### replace NAs with 0
Var_Freq_wide <- Var_Freq_wide %>% replace(is.na(.), 0)
#### spotcheck mutations
Var_Freq_wide_Mat<-Var_Freq_wide[,3:length(Var_Freq_wide)]
rownames(Var_Freq_wide_Mat)<-Var_Freq_wide[,1]

data_1<-Var_Freq_wide[,3:length(Var_Freq_wide)]
data_2<-Var_Freq_wide[,1:2]

NMDS<-metaMDS(data_1, distance = "bray", k =2)

#### Bootstrapping and testing for differences between groups
fit<-adonis2(data_1 ~ Phase, data = data_2, permutations = 999, method = "bray")
fit

#### Check assumption of homogeneity of multivariate dispersion
distances_data<-vegdist(data_1)
anova(betadisper(distances_data, data_2$Phase, bias.adjust = T, sqrt.dist = T  )) #### w

###Data visualiaton
NMDS_2<-cbind(data_2, NMDS$points)
### quick look
ggplot(data = NMDS_2)+
  #geom_point(aes(x = MDS1, y = MDS2, fill = Phase), size = 4, shape = 21, alpha = .3)
  geom_jitter(aes(x = MDS1, y = MDS2, fill = Phase), size = 4, shape = 21, alpha = .7, width = 0.07, height = 0.07)+
  ggtitle("DvH-MMP Variants")

#Data visualisation (THIS IS AN UPDATED VERSION OF THE SCRIPT, NOW USING GGPLOT)
#Extract the axes scores
datascores <- as.data.frame(scores(NMDS)$sites)  #extract the site scores
datascores_sp <- as.data.frame(scores(NMDS)$species)  #extract the variant scores
datascores_sp<-cbind(datascores_sp, "ABS_Dist" = abs(datascores_sp$NMDS1) + abs(datascores_sp$NMDS2))
datascores_sp_ord <- datascores_sp[order(datascores_sp$ABS_Dist, decreasing = T),][1:6,]

#Add/calculate spider diagram
scores <- cbind(as.data.frame(datascores), Phase = data_2$Phase)
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Phase, data = scores, FUN = mean)
seg <- merge(scores, setNames(centroids, c('Phase','oNMDS1','oNMDS2')),
             by = 'Phase', sort = FALSE)

theme_NDMS<-theme(plot.title = element_text(angle = 0, size = 8, face = "italic", hjust = .0, vjust = -2),
                  axis.text = element_text(angle = 0, size = 6),
                  axis.title= element_text(angle = 0, size = 6.8),
                  strip.text = element_text(angle = 0, size = 5, color = "black", margin = margin(c(.07,.07,.07,.07), unit='cm')),
                  legend.box.background = element_rect(fill= "white", color = "grey", size =.25), 
                  legend.title = element_text(size=7),
                  legend.position = "right", legend.text = element_text(size = 7), # c(.09,.89)
                  legend.key.size = unit(.02, "in"), legend.key.width = unit(.02, "in"),
                  legend.justification = 'center', legend.title.align=0.5, legend.margin=margin(c(.05,.05,.05,.05), unit='cm'))

#plot
NMDS_plt<-ggplot(scores, aes(x = NMDS1, y = NMDS2)) +
  geom_vline(aes(xintercept = 0), linetype = 2, alpha =.5, size =.5)+
  geom_hline(aes(yintercept = 0), linetype = 2, alpha =.5, size =.5)+
  #geom_segment(data = seg,
  #             mapping = aes(xend = oNMDS1, yend = oNMDS2)) + # add spiders                  
  geom_segment(data = datascores_sp_ord, aes(x= 0, y=0, xend = NMDS1, yend = NMDS2), 
               arrow = arrow(length = unit(0.005, "npc"), type = "closed"), lwd =.2)+
  geom_jitter(aes(fill = Phase), shape = 23, size = 1.5, alpha = .95, width = 0.07, height = 0.07, stroke =.25) +  
  geom_point(data = centroids,aes(fill = Phase), size = 2.5, shape = 21, alpha =1) +   # add centroids
  geom_text(data = datascores_sp_ord, aes(x = NMDS1*1.2, y = NMDS2*1.2, label = rownames(datascores_sp_ord)), size = 1.4)+
  scale_fill_manual(values = c("Planktonic" =  "#1f78b4", "Attached" = "#a6cee3"), name = "Phase")+
  ggtitle("D. vulgaris + M. maripaludis")+
  theme_bw()+ theme_NDMS+
  theme(legend.position="none",legend.direction='horizontal')
NMDS_plt
ggsave(NMDS_plt, file ="./Final_Plots/NMDS_plt_DvH_Mmp.pdf", width = 2.5, height = 2.5, units = "in", dpi = "retina")
write.csv(scores, "./Source_Data/Figure_2D.csv")
write.csv(datascores_sp_ord, "./Source_Data/Figure_2D_.csv") ## combine in excel
###############################
## Heatmap. ####
Var_Freq<-FBR_Mutations2[,c("Phase_Day_Rep", "Phase", "Day", "Rep", "mutation", "effect", "freq", "gene_id")]
Var_Freq <- as.data.frame(lapply(Var_Freq, function (x) if (is.character(x)) factor(x) else x))

Var_Freq_wide_Mat<-reshape2::dcast(Var_Freq, Phase_Day_Rep + Phase + Day + Rep ~ gene_id, value.var = "freq" )
Var_Freq_wide_Mat<-Var_Freq_wide_Mat[c(19:nrow(Var_Freq_wide_Mat),1:18),]

### replace NAs with 0
Var_Freq_wide_Mat <- Var_Freq_wide_Mat %>% replace(is.na(.), 0)
rownames(Var_Freq_wide_Mat)<-Var_Freq_wide_Mat$Phase_Day_Rep
Var_Freq_wide_Matrix <-t(Var_Freq_wide_Mat[,-c(1:4)])

### build metadata data frame
Var_Freq_Sample_Annotations<-Var_Freq_wide_Mat[,1:4]
Var_Freq_Sample_Annotations$Day<-as.factor(Var_Freq_Sample_Annotations$Day)

#remove duplcated rows
Var_Freq_Annotations<-Var_Freq[,c("mutation", "effect", "gene_id")]
Var_Freq_Annotations<-Var_Freq_Annotations[!duplicated(Var_Freq_Annotations), ]
rownames(Var_Freq_Annotations)<-NULL
rownames(Var_Freq_Annotations)<-Var_Freq_Annotations$gene_id


### use hclust to get the clustering order for row annotations
hc = hclust(dist(Var_Freq_wide_Matrix))

#testing the order of right annotations #df[match(target, df$name),]
## match annotations with #df[match(target, df$name),]
Var_Freq_Annotations[match(hc$labels, Var_Freq_Annotations$gene_id),] #### Test the order of the annotations
Var_Freq_Annotations[match(hc$labels, Var_Freq_Annotations$gene_id),]$effect #### Test the order of the annotations

## define annotation colors paletteer_d("ggthemes::Tableau_20")
col_heat = colorRamp2(c(0,50, 100), c("Blue","white","Red"), space = "RGB")
effect_col<-setNames(as.character(paletteer_d("dichromat::BluetoDarkOrange_18")[c(1,3,2)]),unique(Var_Freq_Annotations[,"effect"]))                                                                        
mutation_col<-setNames(as.character(paletteer_d("ggthemes::Tableau_20")[c(5:8)]),unique(Var_Freq_Annotations[,"mutation"]))                                                                        
Phase_col<- c("Planktonic" =  "#1f78b4", "Attached" = "#a6cee3")
Day_col<-setNames(as.character(paletteer_d("ggsci::deep_purple_material")[c(2:7)]),unique( Var_Freq_Sample_Annotations$Day))                                                                        
Rep_col<-setNames(as.character(paletteer_d("dichromat::GreentoMagenta_16")[c(1:3)]),unique( Var_Freq_Sample_Annotations$Rep))                                                                        

###selecting splits
DvH_splits<-data.frame("gene_id" =rownames(Var_Freq_wide_Matrix), "Group" = c(rep("DvH", 15), rep("Mmp", nrow(Var_Freq_wide_Matrix)-15)))
DvH_splits<-tibble::column_to_rownames(DvH_splits, var = "gene_id")

row_meta =rowAnnotation("Effect" = Var_Freq_Annotations[match(hc$labels, Var_Freq_Annotations$gene_id),]$effect,
                        annotation_name_rot = c(90, 90),
                        annotation_name_gp = gpar(fontsize = 6),
                        col = list("Effect" = effect_col)
)

col_meta =HeatmapAnnotation("Phase" = Var_Freq_Sample_Annotations$Phase,
                            "Day" = Var_Freq_Sample_Annotations$Day,
                            "Rep" = Var_Freq_Sample_Annotations$Rep,
                            #gp = gpar(col = "white"),
                            annotation_name_rot = c(0,0,0),
                            annotation_name_gp = gpar(fontsize = 6), annotation_name_side = "left",
                            col = list("Phase" = Phase_col, "Day" = Day_col, "Rep" = Rep_col)
)

ht1 = Heatmap(Var_Freq_wide_Matrix, name = "mutation frequency",  
              cluster_rows = T, cluster_columns = F, col = col_heat,
              column_names_gp = gpar(fontsize = 0), width = unit(2, "in"),  
              row_names_gp = gpar(fontsize = 6), height = unit(2, "in"),
              row_split = DvH_splits, show_parent_dend_line =F,
              right_annotation = row_meta,
              bottom_annotation = col_meta
)
ht1
pdf("./Final_Plots/Mutation_Freq_Heatmap_DvH_Mmp_RGB.pdf", width = 7, height = 3)
draw(ht1, merge_legend = TRUE, heatmap_legend_side = "right")
dev.off()
write.csv(Var_Freq_wide_Matrix, "./Source_Data/Figure_2E.csv")




