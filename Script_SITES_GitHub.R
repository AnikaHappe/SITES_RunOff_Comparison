#_______________________________________________________________________________
#
# R SCRIPT FOR ANALYSIS OF THE RUN-OFF MESOCOSM EXPERIMENT
# Author: Anika Happe
# Last updated: 31st March 2026
#
#_______________________________________________________________________________

rm(list=ls()) #Empty the environment
setwd("/Users/Anika/Library/Mobile Documents/com~apple~CloudDocs/PhD/AquaCosm/SITES_RunOff_Comparison")

#Load packages
library(readxl)
library(tidyverse)
library(PKNCA) #For AUC function
library(readr)
library(tidyr)
library(dplyr)
library(drc)
library(DescTools)
library(bayestestR)
library(writexl)
library(lubridate)
library(vegan)
library(vegdist)

#Important Papers for Methodology: 
#Hillebrand et al. 2018 (DOI: 10.1111/ele.12867), 
#Urrutia-Cordero et al. 2021 (DOI: 10.1111/1365-2745.13804)

#_______________________________________________________________________________

#### FIG1: Nutrient Additions ####

Table_Additions <- read_excel("R Scripts/Data for R/Table_Additions.xlsx")
Table_Additions <- pivot_longer(Table_Additions, cols = starts_with("Treatment_"), names_to = "treatment", values_to = "value")

plot1 <- ggplot(Table_Additions, aes(x = Day, y = value, fill = treatment)) +
  annotate("rect", xmin = 21, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey87")+
  geom_bar(stat = "identity", position = "dodge", width = 1.0) +
  labs(x = "Day", y = "Value", fill = "Treatment") +
  scale_fill_manual(values = c("Treatment_E" = "indianred", "Treatment_I" = "skyblue3", "Treatment_D" = "palegreen4"))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + 
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20))
plot1

#ggsave("Fig1_Additions.eps", plot1 , unit = "cm", device = "eps", width = 35, height = 13, dpi = 300)
#ggsave("SITES2_Fig1_Additions.png", plot1 , unit = "cm", device = "png", width = 20, height = 13, dpi = 300)

#_______________________________________________________________________________

#Read data and preparation of main data set
data_phy <- read_excel("R Scripts/Data for R/Phytoplankton_Stoichiometry_PANGAEA_Data.xlsx")

str(data_phy) #Check the data structure and transform is necessary
data_phy$Site <- as.factor(data_phy$Site) 
data_phy$Season <- as.factor(data_phy$Season) 
data_phy$Mesocosm <- as.numeric(data_phy$Mesocosm) 
data_phy$Level <- as.factor(data_phy$Level)
data_phy$Treatment <- as.factor(data_phy$Treatment)
data_phy$C_µmol_L <- as.numeric(data_phy$C_µmol_L)
data_phy$N_µmol_L <- as.numeric(data_phy$N_µmol_L)
data_phy$P_µmol_L <- as.numeric(data_phy$P_µmol_L)
data_phy$Si_µmol_L <- as.numeric(data_phy$Si_µmol_L)

#data_phy1 <- filter(data, Level == "Smaller105µm") #We only need this level in this step
data_phy1 <- filter(data_phy, Treatment != "Lake") #We exclude the lake samples
data_phy1 <- filter(data_phy, Treatment != "L") #We exclude the lake samples
data_phy1 <- subset(data_phy1, !(Site == "Erken" & Season == "Summer" & Mesocosm == 2))

#Data from mesocosm #2 is excluded since it received one wrong nutrient addition
#early on and behaved quite differently leading to large standard deviations.

################################################################################

#Calculate stoichiometric molar ratios
data_phy1 <- mutate(data_phy1, CN = C_µmol_L/N_µmol_L, CP = C_µmol_L/P_µmol_L, NP = N_µmol_L/P_µmol_L,
                    CS = C_µmol_L/Si_µmol_L, SP = Si_µmol_L/P_µmol_L, SN = Si_µmol_L/N_µmol_L)

# Define a function to create the plot
template_plot_timeseries <- function(data, response_variable) {
  plot <- ggplot(data, aes_string(x = "Day - 1", y = response_variable, color = "Treatment")) +
    annotate("rect", xmin = 21, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey87") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1) +
    scale_colour_manual(values = c("#B9B9B9", "#0D7B5D", "#54B2E6", "#E69F00")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +  
    theme(axis.title.x = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13)) +
    facet_grid(Season ~ Site)
  
  return(plot)
}

#Check C data for further analysis + create time series plot

data_C <- data_phy1
data_C <- data_C[complete.cases(data_C[, "C_µmol_L"]), ] #Delete NAs (otherwise difficulties later on)
percentile <- quantile(data_C$C_µmol_L, 0.995, na.rm = TRUE) #Calculate 99.5th percentile
boxplot(data_C$C_µmol_L) #Create boxplot
print(percentile) #Get upper whisker value
outliers <- data_C$C_µmol_L[data_C$C_µmol_L > percentile] #Values above 99.5th percentile
cat("Outliers:", length(outliers), "\n") #Count of values above 99.5th percentile
print(outliers) #Get outliers

data_C1 <- filter(data_C, C_µmol_L < 1000) #One outlier is excluded with a C value of < 1000
boxplot(data_C1$C_µmol_L)

data_C1 <- group_by(data_C1, Day, Treatment, Site, Season)
data_C1 <- mutate(data_C1, Mean_C = mean(C_µmol_L), Sd_C = sd(C_µmol_L))

data_C1$Treatment <- factor(data_C1$Treatment, levels = c("C", "D", "I", "E"))
data_C1$Site <- factor(data_C1$Site, levels = c("Erken", "Bolmen"))
plot_C <- template_plot_timeseries(data_C1, "Mean_C")
plot_C <- plot_C + geom_errorbar(aes(ymin = Mean_C - Sd_C, ymax = Mean_C + Sd_C), width = 0.4, size = 0.8)
plot_C

#ggsave("SITES2_Timeseries_C.eps", plot_C , unit = "cm", device = "eps", width = 20, height = 12, dpi = 300)
#ggsave("SITES2_Timeseries_C.png", plot_C , unit = "cm", device = "png", width = 17, height = 10, dpi = 300)

#All next to each other
plot_C2 <- template_plot_timeseries(data_C1, "Mean_C")
plot_C2 <- plot_C2 + geom_errorbar(aes(ymin = Mean_C - Sd_C, ymax = Mean_C + Sd_C), width = 0.4, size = 0.8) +
  theme(text = element_text(family = "Times New Roman"))+
  facet_wrap(Site~Season)
plot_C2

#ggsave("SITES2_Timeseries_C_1row.eps", plot_C2 , unit = "cm", device = "eps", width = 25, height = 8, dpi = 300)
#ggsave("SITES2_Timeseries_C_June25.png", plot_C2 , unit = "cm", device = "png", width = 25, height = 8, dpi = 300)

#Check CN data for further analysis + create time series plot

data_CN <- data_phy1
data_CN <- data_CN[complete.cases(data_CN[, "CN"]), ]
percentile <- quantile(data_CN$CN, 0.995, na.rm = TRUE) #Calculate 99.5th percentile
boxplot(data_CN$CN) #Create boxplot
print(percentile) #Get upper whisker value
outliers <- data_CN$CN[data_CN$CN > percentile] #Values above 99.5th percentile
cat("Outliers:", length(outliers), "\n") #Count of values above 99.5th percentile
print(outliers) #Get outliers

data_CN1 <- filter(data_CN, C_µmol_L < 1000)
boxplot(data_CN1$CN)

data_CN1 <- group_by(data_CN1, Day, Treatment, Site, Season)
data_CN1 <- mutate(data_CN1, Mean_CN = mean(CN), Sd_CN = sd(CN))

data_CN1$Treatment <- factor(data_CN1$Treatment, levels = c("C", "D", "I", "E"))
data_C1$Site <- factor(data_C1$Site, levels = c("Erken", "Bolmen"))
plot_CN <- template_plot_timeseries(data_CN1, "Mean_CN")
plot_CN <- plot_CN + geom_errorbar(aes(ymin = Mean_CN - Sd_CN, ymax = Mean_CN + Sd_CN), width = 0.4, size = 0.9)
plot_CN

#ggsave("SITES2_Timeseries_CN.eps", plot_CN , unit = "cm", device = "eps", width = 23, height = 15, dpi = 300)
#ggsave("SITES2_Timeseries_CN.png", plot_CN , unit = "cm", device = "png", width = 17, height = 10, dpi = 300)

#All next to each other
plot_CN2 <- template_plot_timeseries(data_CN1, "Mean_CN")
plot_CN2 <- plot_CN2 + geom_errorbar(aes(ymin = Mean_CN - Sd_CN, ymax = Mean_CN + Sd_CN), width = 0.4, size = 0.8) +
  theme(text = element_text(family = "Times New Roman"))+
  facet_wrap(Site~Season)
plot_CN2

#ggsave("SITES2_Timeseries_C_1row.eps", plot_C2 , unit = "cm", device = "eps", width = 25, height = 8, dpi = 300)
#ggsave("SITES2_Timeseries_CN_June25.png", plot_CN2 , unit = "cm", device = "png", width = 25, height = 8, dpi = 300)

#Check CP data for further analysis + create time series plot

data_CP <- data_phy1
data_CP <- data_CP[complete.cases(data_CP[, "CP"]), ]
percentile <- quantile(data_CP$CP, 0.995, na.rm = TRUE) #Calculate 99.5th percentile
boxplot(data_CP$CP) #Create boxplot
print(percentile) #Get upper whisker value
outliers <- data_CP$CP[data_CP$CP > percentile] #Values above 99.5th percentile
cat("Outliers:", length(outliers), "\n") #Count of values above 99.5th percentile
print(outliers) #Get outliers

data_CP1 <- filter(data_CP, CP < 1500)
boxplot(data_CP1$CP)

data_CP1 <- group_by(data_CP1, Day, Treatment, Site, Season)
data_CP1 <- mutate(data_CP1, Mean_CP = mean(CP), Sd_CP = sd(CP))

data_CP1$Treatment <- factor(data_CP1$Treatment, levels = c("C", "D", "I", "E"))
data_CP1$Site <- factor(data_CP1$Site, levels = c("Erken", "Bolmen"))
plot_CP <- template_plot_timeseries(data_CP1, "Mean_CP")
plot_CP <- plot_CP + geom_errorbar(aes(ymin = Mean_CP - Sd_CP, ymax = Mean_CP + Sd_CP), width = 0.4, size = 0.9)
plot_CP

#ggsave("SITES2_Timeseries_CP.eps", plot_CP , unit = "cm", device = "eps", width = 23, height = 15, dpi = 300)
#ggsave("SITES2_Timeseries_CP.png", plot_CP , unit = "cm", device = "png", width = 17, height = 10, dpi = 300)

#All next to each other
plot_CP2 <- template_plot_timeseries(data_CP1, "Mean_CP")
plot_CP2 <- plot_CP2 + geom_errorbar(aes(ymin = Mean_CP - Sd_CP, ymax = Mean_CP + Sd_CP), width = 0.4, size = 0.8) +
  theme(text = element_text(family = "Times New Roman"))+
  facet_wrap(Site~Season)
plot_CP2

#ggsave("SITES2_Timeseries_C_1row.eps", plot_C2 , unit = "cm", device = "eps", width = 25, height = 8, dpi = 300)
#ggsave("SITES2_Timeseries_CP_June25.png", plot_CP2 , unit = "cm", device = "png", width = 25, height = 8, dpi = 300)

#Check NP data for further analysis + create time series plot

data_NP <- data_phy1
data_NP <- data_NP[complete.cases(data_NP[, "NP"]), ]
percentile <- quantile(data_NP$NP, 0.995, na.rm = TRUE) #Calculate 99.5th percentile
boxplot(data_NP$NP) #Create boxplot
print(percentile) #Get upper whisker value
outliers <- data_NP$NP[data_NP$NP > percentile] #Values above 99.5th percentile
cat("Outliers:", length(outliers), "\n") #Count of values above 99.5th percentile
print(outliers) #Get outliers

data_NP1 <- filter(data_NP, NP < 100)

data_NP1 <- group_by(data_NP1, Day, Treatment, Site, Season)
data_NP1 <- mutate(data_NP1, Mean_NP = mean(NP), Sd_NP = sd(NP))

data_NP1$Treatment <- factor(data_NP1$Treatment, levels = c("C", "D", "I", "E"))
data_NP1$Site <- factor(data_NP1$Site, levels = c("Erken", "Bolmen"))
plot_NP <- template_plot_timeseries(data_NP1, "Mean_NP")
plot_NP <- plot_NP + geom_errorbar(aes(ymin = Mean_NP - Sd_NP, ymax = Mean_NP + Sd_NP), width = 0.4, size = 0.9)
plot_NP

#ggsave("SITES2_Timeseries_NP.eps", plot_NP , unit = "cm", device = "eps", width = 23, height = 15, dpi = 300)
#ggsave("SITES2_Timeseries_NP.png", plot_NP , unit = "cm", device = "png", width = 17, height = 10, dpi = 300)

#Check CS data for further analysis + create time series plot

data_CS <- data_phy1
data_CS <- data_CS[complete.cases(data_CS[, "CS"]), ]
percentile <- quantile(data_CS$CS, 0.995, na.rm = TRUE) #Calculate 99.5th percentile
boxplot(data_CS$CS) #Create boxplot
print(percentile) #Get upper whisker value
outliers <- data_CS$CS[data_CS$CS > percentile] #Values above 99.5th percentile
cat("Outliers:", length(outliers), "\n") #Count of values above 99.5th percentile
print(outliers) #Get outliers

#The outliers are all the same treatment on the same day, so it is not an outlier
#But from the time series plots, there was an outlier visible for Bolmen, summer day 0, mesocosm 15 due to a very high C value
data_CS <- filter(data_CS, C_µmol_L < 1000) #One outlier is excluded with a C value of < 1000

data_CS1 <- group_by(data_CS, Day, Treatment, Site, Season)
data_CS1 <- mutate(data_CS1, Mean_CS = mean(CS), Sd_CS = sd(CS))

data_CS1$Treatment <- factor(data_CS1$Treatment, levels = c("C", "D", "I", "E"))
data_CS1$Site <- factor(data_CS1$Site, levels = c("Erken", "Bolmen"))
plot_CS <- template_plot_timeseries(data_CS1, "Mean_CS")
plot_CS <- plot_CS + geom_errorbar(aes(ymin = Mean_CS - Sd_CS, ymax = Mean_CS + Sd_CS), width = 0.4, size = 0.9)+facet_grid(Site ~ Season, scales = "free")
plot_CS

#ggsave("SITES2_Timeseries_CS.eps", plot_CS , unit = "cm", device = "eps", width = 23, height = 15, dpi = 300)

#All next to each other
plot_CS2 <- template_plot_timeseries(data_CS1, "Mean_CS")
plot_CS2 <- plot_CS2 + geom_errorbar(aes(ymin = Mean_CS - Sd_CS, ymax = Mean_CS + Sd_CS), width = 0.4, size = 0.8) +
  theme(text = element_text(family = "Times New Roman"))+
  facet_wrap(Site~Season)
plot_CS2

#ggsave("SITES2_Timeseries_C_1row.eps", plot_C2 , unit = "cm", device = "eps", width = 25, height = 8, dpi = 300)
#ggsave("SITES2_Timeseries_CS_June25.png", plot_CS2 , unit = "cm", device = "png", width = 25, height = 8, dpi = 300)

#Check SN data for further analysis + create time series plot

data_SN <- data_phy1
data_SN <- data_SN[complete.cases(data_SN[, "SN"]), ]
percentile <- quantile(data_SN$SN, 0.995, na.rm = TRUE) #Calculate 99.5th percentile
boxplot(data_SN$SN) #Create boxplot
print(percentile) #Get upper whisker value
outliers <- data_SN$SN[data_SN$SN > percentile] #Values above 99.5th percentile
cat("Outliers:", length(outliers), "\n") #Count of values above 99.5th percentile
print(outliers) #Get outliers

#The outliers are all the same treatment on the same day, so it is not an outlier
#But from the time series plots, there was an outlier visible for Bolmen, summer day 0, mesocosm 15 due to a very high C value
#data_SN <- filter(data_SN, C_µmol_L < 1000) #One outlier is excluded with a C value of < 1000

data_SN1 <- group_by(data_SN, Day, Treatment, Site, Season)
data_SN1 <- mutate(data_SN1, Mean_SN = mean(SN), Sd_SN = sd(SN))

data_SN1$Treatment <- factor(data_SN1$Treatment, levels = c("C", "D", "I", "E"))
data_SN1$Site <- factor(data_SN1$Site, levels = c("Erken", "Bolmen"))
plot_SN <- template_plot_timeseries(data_SN1, "Mean_SN")
plot_SN <- plot_SN + geom_errorbar(aes(ymin = Mean_SN - Sd_SN, ymax = Mean_SN + Sd_SN), width = 0.4, size = 0.9)+facet_grid(Season ~ Site, scales = "free_y")
plot_SN

#ggsave("SITES2_Timeseries_SN.eps", plot_SN , unit = "cm", device = "eps", width = 23, height = 15, dpi = 300)

#Check SP data for further analysis + create time series plot

data_SP <- data_phy1
data_SP <- data_SP[complete.cases(data_SP[, "SP"]), ]
percentile <- quantile(data_SP$SP, 0.995, na.rm = TRUE) #Calculate 99.5th percentile
boxplot(data_SP$SP) #Create boxplot
print(percentile) #Get upper whisker value
outliers <- data_SP$SP[data_SP$SP > percentile] #Values above 99.5th percentile
cat("Outliers:", length(outliers), "\n") #Count of values above 99.5th percentile
print(outliers) #Get outliers

data_SP1 <- filter(data_SP, SP < 200) 

data_SP1 <- group_by(data_SP1, Day, Treatment, Site, Season)
data_SP1 <- mutate(data_SP1, Mean_SP = mean(SP), Sd_SP = sd(SP))

data_SP1$Treatment <- factor(data_SP1$Treatment, levels = c("C", "D", "I", "E"))
data_SP1$Site <- factor(data_SP1$Site, levels = c("Erken", "Bolmen"))
plot_SP <- template_plot_timeseries(data_SP1, "Mean_SP")
plot_SP <- plot_SP + geom_errorbar(aes(ymin = Mean_SP - Sd_SP, ymax = Mean_SP + Sd_SP), width = 0.4, size = 0.9)+facet_grid(Season ~ Site, scales = "free_y")
plot_SP

ggsave("SITES2_Timeseries_SP.eps", plot_SP , unit = "cm", device = "eps", width = 23, height = 15, dpi = 300)

#_______________________________________________________________________________

#### FIG: LRRs, OEV and recovery of phytoplankton ####

#### RECOVERY ####

#To include the variance, we calculate the LRR for each replicate. So,
#ln(disturbed from replicate 1 / mean control) and so on. With this
#we can calculate a LRR in the end.

#1) Paste the value of the mean undisturbed (control) treatment in a new column
#2) Calculate the LRR for Recovery = ln(Disturbed/Control) for each response

#C
data_phy_C_recovery <- (data_C1 %>% dplyr::filter(Sno == 10) %>% dplyr::group_by(Season, Site, Treatment) %>% dplyr::mutate(Control_C = mean(C_µmol_L)))
data_phy_C_recovery$Control_C[data_phy_C_recovery$Site == "Erken" & data_phy_C_recovery$Season == "Summer"] <- data_phy_C_recovery$Control_C[data_phy_C_recovery$ID == 150] 
data_phy_C_recovery$Control_C[data_phy_C_recovery$Site == "Erken" & data_phy_C_recovery$Season == "Spring"] <- data_phy_C_recovery$Control_C[data_phy_C_recovery$ID == 320] 
data_phy_C_recovery$Control_C[data_phy_C_recovery$Site == "Bolmen" & data_phy_C_recovery$Season == "Summer"] <- data_phy_C_recovery$Control_C[data_phy_C_recovery$ID == 462] 
data_phy_C_recovery <- filter (data_phy_C_recovery, Treatment != "C") #Not needed anymore!
data_phy_C_recovery <- mutate(data_phy_C_recovery, Recovery = log(C_µmol_L/Control_C), Stab_Para = "C")
data_phy_C_recovery <- group_by(data_phy_C_recovery, Site, Season, Treatment)
data_phy_C_recovery_mean <- mutate(data_phy_C_recovery, mean_Recovery = mean(Recovery), sd_Recovery = sd(Recovery))
data_phy_C_recovery_mean <- filter(data_phy_C_recovery_mean, Replicate == 4) #They should all be the same now!

#CN
data_phy_CN_recovery <- (data_CN1 %>% dplyr::filter(Sno == 10) %>% dplyr::group_by(Season, Site, Treatment) %>% dplyr::mutate(Control_CN = mean(CN)))
data_phy_CN_recovery$Control_CN[data_phy_CN_recovery$Site == "Erken" & data_phy_CN_recovery$Season == "Summer"] <- data_phy_CN_recovery$Control_CN[data_phy_CN_recovery$ID == 150] 
data_phy_CN_recovery$Control_CN[data_phy_CN_recovery$Site == "Erken" & data_phy_CN_recovery$Season == "Spring"] <- data_phy_CN_recovery$Control_CN[data_phy_CN_recovery$ID == 320] 
data_phy_CN_recovery$Control_CN[data_phy_CN_recovery$Site == "Bolmen" & data_phy_CN_recovery$Season == "Summer"] <- data_phy_CN_recovery$Control_CN[data_phy_CN_recovery$ID == 462] 
data_phy_CN_recovery <- filter (data_phy_CN_recovery, Treatment != "C") #Not needed anymore!
data_phy_CN_recovery <- mutate(data_phy_CN_recovery, Recovery = log(CN/Control_CN), Stab_Para = "CN")
data_phy_CN_recovery <- group_by(data_phy_CN_recovery, Site, Season, Treatment)
data_phy_CN_recovery_mean <- mutate(data_phy_CN_recovery, mean_Recovery = mean(Recovery), sd_Recovery = sd(Recovery))
data_phy_CN_recovery_mean <- filter(data_phy_CN_recovery_mean, Replicate == 4) #They should all be the same now!

#CP
data_phy_CP_recovery <- (data_CP1 %>% dplyr::filter(Sno == 10) %>% dplyr::group_by(Season, Site, Treatment) %>% dplyr::mutate(Control_CP = mean(CP)))
data_phy_CP_recovery$Control_CP[data_phy_CP_recovery$Site == "Erken" & data_phy_CP_recovery$Season == "Summer"] <- data_phy_CP_recovery$Control_CP[data_phy_CP_recovery$ID == 150] 
data_phy_CP_recovery$Control_CP[data_phy_CP_recovery$Site == "Erken" & data_phy_CP_recovery$Season == "Spring"] <- data_phy_CP_recovery$Control_CP[data_phy_CP_recovery$ID == 320] 
data_phy_CP_recovery$Control_CP[data_phy_CP_recovery$Site == "Bolmen" & data_phy_CP_recovery$Season == "Summer"] <- data_phy_CP_recovery$Control_CP[data_phy_CP_recovery$ID == 462] 
data_phy_CP_recovery <- filter (data_phy_CP_recovery, Treatment != "C") #Not needed anymore!
data_phy_CP_recovery <- mutate(data_phy_CP_recovery, Recovery = log(CP/Control_CP), Stab_Para = "CP")
data_phy_CP_recovery <- group_by(data_phy_CP_recovery, Site, Season, Treatment)
data_phy_CP_recovery_mean <- mutate(data_phy_CP_recovery, mean_Recovery = mean(Recovery), sd_Recovery = sd(Recovery))
data_phy_CP_recovery_mean <- filter(data_phy_CP_recovery_mean, Replicate == 4) #They should all be the same now!

#NP
data_phy_NP_recovery <- (data_NP1 %>% dplyr::filter(Sno == 10) %>% dplyr::group_by(Season, Site, Treatment) %>% dplyr::mutate(Control_NP = mean(NP)))
data_phy_NP_recovery$Control_NP[data_phy_NP_recovery$Site == "Erken" & data_phy_NP_recovery$Season == "Summer"] <- data_phy_NP_recovery$Control_NP[data_phy_NP_recovery$ID == 150] 
data_phy_NP_recovery$Control_NP[data_phy_NP_recovery$Site == "Erken" & data_phy_NP_recovery$Season == "Spring"] <- data_phy_NP_recovery$Control_NP[data_phy_NP_recovery$ID == 320] 
data_phy_NP_recovery$Control_NP[data_phy_NP_recovery$Site == "Bolmen" & data_phy_NP_recovery$Season == "Summer"] <- data_phy_NP_recovery$Control_NP[data_phy_NP_recovery$ID == 462] 
data_phy_NP_recovery <- filter (data_phy_NP_recovery, Treatment != "C") #Not needed anymore!
data_phy_NP_recovery <- mutate(data_phy_NP_recovery, Recovery = log(NP/Control_NP), Stab_Para = "NP")
data_phy_NP_recovery <- group_by(data_phy_NP_recovery, Site, Season, Treatment)
data_phy_NP_recovery_mean <- mutate(data_phy_NP_recovery, mean_Recovery = mean(Recovery), sd_Recovery = sd(Recovery))
data_phy_NP_recovery_mean <- filter(data_phy_NP_recovery_mean, Replicate == 4) #They should all be the same now!

#CS
data_phy_CS_recovery <- (data_CS1 %>% dplyr::filter(Sno == 10) %>% dplyr::group_by(Season, Site, Treatment) %>% dplyr::mutate(Control_CS = mean(CS)))
data_phy_CS_recovery$Control_CS[data_phy_CS_recovery$Site == "Erken" & data_phy_CS_recovery$Season == "Summer"] <- data_phy_CS_recovery$Control_CS[data_phy_CS_recovery$ID == 150] 
data_phy_CS_recovery$Control_CS[data_phy_CS_recovery$Site == "Erken" & data_phy_CS_recovery$Season == "Spring"] <- data_phy_CS_recovery$Control_CS[data_phy_CS_recovery$ID == 320] 
data_phy_CS_recovery$Control_CS[data_phy_CS_recovery$Site == "Bolmen" & data_phy_CS_recovery$Season == "Summer"] <- data_phy_CS_recovery$Control_CS[data_phy_CS_recovery$ID == 462] 
data_phy_CS_recovery <- filter (data_phy_CS_recovery, Treatment != "C") #Not needed anymore!
data_phy_CS_recovery <- mutate(data_phy_CS_recovery, Recovery = log(CS/Control_CS), Stab_Para = "CS")
data_phy_CS_recovery <- group_by(data_phy_CS_recovery, Site, Season, Treatment)
data_phy_CS_recovery_mean <- mutate(data_phy_CS_recovery, mean_Recovery = mean(Recovery), sd_Recovery = sd(Recovery))
data_phy_CS_recovery_mean <- filter(data_phy_CS_recovery_mean, Replicate == 4) #They should all be the same now!

#SN
data_phy_SN_recovery <- (data_SN1 %>% dplyr::filter(Sno == 10) %>% dplyr::group_by(Season, Site, Treatment) %>% dplyr::mutate(Control_SN = mean(SN)))
data_phy_SN_recovery$Control_SN[data_phy_SN_recovery$Site == "Erken" & data_phy_SN_recovery$Season == "Summer"] <- data_phy_SN_recovery$Control_SN[data_phy_SN_recovery$ID == 150] 
data_phy_SN_recovery$Control_SN[data_phy_SN_recovery$Site == "Erken" & data_phy_SN_recovery$Season == "Spring"] <- data_phy_SN_recovery$Control_SN[data_phy_SN_recovery$ID == 320] 
data_phy_SN_recovery$Control_SN[data_phy_SN_recovery$Site == "Bolmen" & data_phy_SN_recovery$Season == "Summer"] <- data_phy_SN_recovery$Control_SN[data_phy_SN_recovery$ID == 462] 
data_phy_SN_recovery <- filter (data_phy_SN_recovery, Treatment != "C") #Not needed anymore!
data_phy_SN_recovery <- mutate(data_phy_SN_recovery, Recovery = log(SN/Control_SN), Stab_Para = "SN")
data_phy_SN_recovery <- group_by(data_phy_SN_recovery, Site, Season, Treatment)
data_phy_SN_recovery_mean <- mutate(data_phy_SN_recovery, mean_Recovery = mean(Recovery), sd_Recovery = sd(Recovery))
data_phy_SN_recovery_mean <- filter(data_phy_SN_recovery_mean, Replicate == 4) #They should all be the same now!

#SP
data_phy_SP_recovery <- (data_SP1 %>% dplyr::filter(Sno == 10) %>% dplyr::group_by(Season, Site, Treatment) %>% dplyr::mutate(Control_SP = mean(SP)))
data_phy_SP_recovery$Control_SP[data_phy_SP_recovery$Site == "Erken" & data_phy_SP_recovery$Season == "Summer"] <- data_phy_SP_recovery$Control_SP[data_phy_SP_recovery$ID == 150] 
data_phy_SP_recovery$Control_SP[data_phy_SP_recovery$Site == "Erken" & data_phy_SP_recovery$Season == "Spring"] <- data_phy_SP_recovery$Control_SP[data_phy_SP_recovery$ID == 320] 
data_phy_SP_recovery$Control_SP[data_phy_SP_recovery$Site == "Bolmen" & data_phy_SP_recovery$Season == "Summer"] <- data_phy_SP_recovery$Control_SP[data_phy_SP_recovery$ID == 462] 
data_phy_SP_recovery <- filter (data_phy_SP_recovery, Treatment != "C") #Not needed anymore!
data_phy_SP_recovery <- mutate(data_phy_SP_recovery, Recovery = log(SP/Control_SP), Stab_Para = "SP")
data_phy_SP_recovery <- group_by(data_phy_SP_recovery, Site, Season, Treatment)
data_phy_SP_recovery_mean <- mutate(data_phy_SP_recovery, mean_Recovery = mean(Recovery), sd_Recovery = sd(Recovery))
data_phy_SP_recovery_mean <- filter(data_phy_SP_recovery_mean, Replicate == 4) #They should all be the same now!

#Put all the data sets together
data_recovery_mean <- rbind(data_phy_C_recovery_mean, data_phy_CN_recovery_mean, data_phy_CP_recovery_mean, data_phy_NP_recovery_mean,
                            data_phy_CS_recovery_mean, data_phy_SN_recovery_mean, data_phy_SP_recovery_mean)
data_recovery_mean <- dplyr::select(data_recovery_mean, Site, Season, Mesocosm, Treatment, Stab_Para, mean_Recovery, sd_Recovery)
data_recovery <- rbind(data_phy_C_recovery, data_phy_CN_recovery, data_phy_CP_recovery, data_phy_NP_recovery,
                       data_phy_CS_recovery, data_phy_SN_recovery, data_phy_SP_recovery)
data_recovery <- dplyr::select(data_recovery, Site, Season, Mesocosm, Treatment, Replicate, Stab_Para, Recovery)

#write_xlsx(data_recovery, "SITES_2022_SmallPhyto_Recovery.xlsx")

#### PREP FOR AUC: CALCULATE LOG-RESPONSE RATIOS (LRR) ####

# Define a function to create the plot
template_plot_LRR <- function(data, response_variable) {
  plot <- ggplot(data, aes_string(x = "Day-1", y = response_variable, color = "Treatment")) +
    annotate("rect", xmin = 21, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "grey87") +
    geom_line(linewidth = 1.2) +
    geom_point(size = 1.7) +
    scale_colour_manual(values = c("#0D7B5D", "#54B2E6", "#E69F00")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +  
    theme(axis.title.x = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13)) +
    facet_grid(Season ~ Site)
  
  return(plot)
}

#Needs to be done for each experiment individually (I did not find a smoother solution).
#1) Put the mean value of the control for each sampling day in a new column (with the for-loop)
#2) Put data sets (of individual experiments) back together

#C (all three experiments)
data_C2 <- data_C1
data_C2$Sno <- as.numeric(data_C2$Sno)
data_C2 <-mutate(data_C2, Mean_C_Control = NA)
data_phy_C_LRR_Erken_Summer <- filter(data_C2, Site == "Erken" & Season == "Summer")
data_phy_C_LRR_Erken_Spring <- filter(data_C2, Site == "Erken" & Season == "Spring")
data_phy_C_LRR_Bolmen_Summer <- filter(data_C2, Site == "Bolmen" & Season == "Summer")
for (i in unique(data_phy_C_LRR_Erken_Summer$Sno)) {
  data_phy_C_LRR_Erken_Summer$Mean_C_Control[data_phy_C_LRR_Erken_Summer$Sno == i] <-
    data_phy_C_LRR_Erken_Summer$Mean_C[data_phy_C_LRR_Erken_Summer$Sno == i & data_phy_C_LRR_Erken_Summer$Treatment == "C" & data_phy_C_LRR_Erken_Summer$Replicate == 4]
}
for (i in unique(data_phy_C_LRR_Erken_Spring$Sno)) {
  data_phy_C_LRR_Erken_Spring$Mean_C_Control[data_phy_C_LRR_Erken_Spring$Sno == i] <-
    data_phy_C_LRR_Erken_Spring$Mean_C[data_phy_C_LRR_Erken_Spring$Sno == i & data_phy_C_LRR_Erken_Spring$Treatment == "C" & data_phy_C_LRR_Erken_Spring$Replicate == 4]
}
for (i in unique(data_phy_C_LRR_Bolmen_Summer$Sno)) {
  data_phy_C_LRR_Bolmen_Summer$Mean_C_Control[data_phy_C_LRR_Bolmen_Summer$Sno == i] <-
    data_phy_C_LRR_Bolmen_Summer$Mean_C[data_phy_C_LRR_Bolmen_Summer$Sno == i & data_phy_C_LRR_Bolmen_Summer$Treatment == "C" & data_phy_C_LRR_Bolmen_Summer$Replicate == 1]
}
data_C_LRR <- rbind(data_phy_C_LRR_Erken_Summer, data_phy_C_LRR_Erken_Spring, data_phy_C_LRR_Bolmen_Summer) #Put data sets back together
data_C_LRR <- mutate(data_C_LRR, LRR_C = log(C_µmol_L/Mean_C_Control))
data_C_LRR_mean <- group_by(data_C_LRR, Site, Season, Treatment, Sno)
data_C_LRR_mean <- mutate(data_C_LRR_mean, mean_LRR_C = mean(LRR_C), sd_LRR_C = sd(LRR_C))
data_C_LRR_mean <- data_C_LRR_mean %>% filter(ifelse(Site == "Erken", Replicate == 4, Replicate == 1))
data_C_LRR_mean <- filter(data_C_LRR_mean, Treatment != "C")
plot_LRR_C <- template_plot_LRR(data_C_LRR_mean, "mean_LRR_C")+geom_errorbar(aes(ymin = mean_LRR_C - sd_LRR_C, ymax = mean_LRR_C + sd_LRR_C), width = 0.4, size = 0.9)
plot_LRR_C
#ggsave("SITES2_LRR_C.eps", plot_LRR_C , unit = "cm", device = "eps", width = 23, height = 15, dpi = 300)

#CN (all three experiments)
data_CN2 <- data_CN1
data_CN2$Sno <- as.numeric(data_CN2$Sno)
data_CN2 <-mutate(data_CN2, Mean_CN_Control = NA)
data_phy_CN_LRR_Erken_Summer <- filter(data_CN2, Site == "Erken" & Season == "Summer")
data_phy_CN_LRR_Erken_Spring <- filter(data_CN2, Site == "Erken" & Season == "Spring")
data_phy_CN_LRR_Bolmen_Summer <- filter(data_CN2, Site == "Bolmen" & Season == "Summer")
for (i in unique(data_phy_CN_LRR_Erken_Summer$Sno)) {
  data_phy_CN_LRR_Erken_Summer$Mean_CN_Control[data_phy_CN_LRR_Erken_Summer$Sno == i] <-
    data_phy_CN_LRR_Erken_Summer$Mean_CN[data_phy_CN_LRR_Erken_Summer$Sno == i & data_phy_CN_LRR_Erken_Summer$Treatment == "C" & data_phy_CN_LRR_Erken_Summer$Replicate == 4]
}
for (i in unique(data_phy_CN_LRR_Erken_Spring$Sno)) {
  data_phy_CN_LRR_Erken_Spring$Mean_CN_Control[data_phy_CN_LRR_Erken_Spring$Sno == i] <-
    data_phy_CN_LRR_Erken_Spring$Mean_CN[data_phy_CN_LRR_Erken_Spring$Sno == i & data_phy_CN_LRR_Erken_Spring$Treatment == "C" & data_phy_CN_LRR_Erken_Spring$Replicate == 4]
}
for (i in unique(data_phy_CN_LRR_Bolmen_Summer$Sno)) {
  data_phy_CN_LRR_Bolmen_Summer$Mean_CN_Control[data_phy_CN_LRR_Bolmen_Summer$Sno == i] <-
    data_phy_CN_LRR_Bolmen_Summer$Mean_CN[data_phy_CN_LRR_Bolmen_Summer$Sno == i & data_phy_CN_LRR_Bolmen_Summer$Treatment == "C" & data_phy_CN_LRR_Bolmen_Summer$Replicate == 1]
}
data_CN_LRR <- rbind(data_phy_CN_LRR_Erken_Summer, data_phy_CN_LRR_Erken_Spring, data_phy_CN_LRR_Bolmen_Summer) #Put data sets back together
data_CN_LRR <- mutate(data_CN_LRR, LRR_CN = log(CN/Mean_CN_Control))
data_CN_LRR_mean <- group_by(data_CN_LRR, Site, Season, Treatment, Sno)
data_CN_LRR_mean <- mutate(data_CN_LRR_mean, mean_LRR_CN = mean(LRR_CN), sd_LRR_CN = sd(LRR_CN))
data_CN_LRR_mean <- data_CN_LRR_mean %>% filter(ifelse(Site == "Erken", Replicate == 4, Replicate == 1))
data_CN_LRR_mean <- filter(data_CN_LRR_mean, Treatment != "C")
plot_LRR_CN <- template_plot_LRR(data_CN_LRR_mean, "mean_LRR_CN")+geom_errorbar(aes(ymin = mean_LRR_CN - sd_LRR_CN, ymax = mean_LRR_CN + sd_LRR_CN), width = 0.4, size = 0.9)
plot_LRR_CN
#ggsave("SITES2_LRR_CN.eps", plot_LRR_CN , unit = "cm", device = "eps", width = 23, height = 15, dpi = 300)

#CP (all three experiments)
data_CP2 <- data_CP1
data_CP2$Sno <- as.numeric(data_CP2$Sno)
data_CP2 <-mutate(data_CP2, Mean_CP_Control = NA)
data_phy_CP_LRR_Erken_Summer <- filter(data_CP2, Site == "Erken" & Season == "Summer")
data_phy_CP_LRR_Erken_Spring <- filter(data_CP2, Site == "Erken" & Season == "Spring")
data_phy_CP_LRR_Bolmen_Summer <- filter(data_CP2, Site == "Bolmen" & Season == "Summer")
for (i in unique(data_phy_CP_LRR_Erken_Summer$Sno)) {
  data_phy_CP_LRR_Erken_Summer$Mean_CP_Control[data_phy_CP_LRR_Erken_Summer$Sno == i] <-
    data_phy_CP_LRR_Erken_Summer$Mean_CP[data_phy_CP_LRR_Erken_Summer$Sno == i & data_phy_CP_LRR_Erken_Summer$Treatment == "C" & data_phy_CP_LRR_Erken_Summer$Replicate == 3]
}
for (i in unique(data_phy_CP_LRR_Erken_Spring$Sno)) {
  data_phy_CP_LRR_Erken_Spring$Mean_CP_Control[data_phy_CP_LRR_Erken_Spring$Sno == i] <-
    data_phy_CP_LRR_Erken_Spring$Mean_CP[data_phy_CP_LRR_Erken_Spring$Sno == i & data_phy_CP_LRR_Erken_Spring$Treatment == "C" & data_phy_CP_LRR_Erken_Spring$Replicate == 4]
}
for (i in unique(data_phy_CP_LRR_Bolmen_Summer$Sno)) {
  data_phy_CP_LRR_Bolmen_Summer$Mean_CP_Control[data_phy_CP_LRR_Bolmen_Summer$Sno == i] <-
    data_phy_CP_LRR_Bolmen_Summer$Mean_CP[data_phy_CP_LRR_Bolmen_Summer$Sno == i & data_phy_CP_LRR_Bolmen_Summer$Treatment == "C" & data_phy_CP_LRR_Bolmen_Summer$Replicate == 3]
}
data_CP_LRR <- rbind(data_phy_CP_LRR_Erken_Summer, data_phy_CP_LRR_Erken_Spring, data_phy_CP_LRR_Bolmen_Summer) #Put data sets back together
data_CP_LRR <- mutate(data_CP_LRR, LRR_CP = log(CP/Mean_CP_Control))
data_CP_LRR_mean <- group_by(data_CP_LRR, Site, Season, Treatment, Sno)
data_CP_LRR_mean <- mutate(data_CP_LRR_mean, mean_LRR_CP = mean(LRR_CP), sd_LRR_CP = sd(LRR_CP))
data_CP_LRR_mean <- data_CP_LRR_mean %>% filter(ifelse(Site == "Erken", Replicate == 4, Replicate == 1))
data_CP_LRR_mean <- filter(data_CP_LRR_mean, Treatment != "C")
plot_LRR_CP <- template_plot_LRR(data_CP_LRR_mean, "mean_LRR_CP")+geom_errorbar(aes(ymin = mean_LRR_CP - sd_LRR_CP, ymax = mean_LRR_CP + sd_LRR_CP), width = 0.4, size = 0.9)
plot_LRR_CP
#ggsave("SITES2_LRR_CP.eps", plot_LRR_CP , unit = "cm", device = "eps", width = 23, height = 15, dpi = 300)

#NP (all three experiments)
data_NP2 <- data_NP1
data_NP2$Sno <- as.numeric(data_NP2$Sno)
data_NP2 <-mutate(data_NP2, Mean_NP_Control = NA)
data_phy_NP_LRR_Erken_Summer <- filter(data_NP2, Site == "Erken" & Season == "Summer")
data_phy_NP_LRR_Erken_Spring <- filter(data_NP2, Site == "Erken" & Season == "Spring")
data_phy_NP_LRR_Bolmen_Summer <- filter(data_NP2, Site == "Bolmen" & Season == "Summer")
for (i in unique(data_phy_NP_LRR_Erken_Summer$Sno)) {
  data_phy_NP_LRR_Erken_Summer$Mean_NP_Control[data_phy_NP_LRR_Erken_Summer$Sno == i] <-
    data_phy_NP_LRR_Erken_Summer$Mean_NP[data_phy_NP_LRR_Erken_Summer$Sno == i & data_phy_NP_LRR_Erken_Summer$Treatment == "C" & data_phy_NP_LRR_Erken_Summer$Replicate == 3]
}
for (i in unique(data_phy_NP_LRR_Erken_Spring$Sno)) {
  data_phy_NP_LRR_Erken_Spring$Mean_NP_Control[data_phy_NP_LRR_Erken_Spring$Sno == i] <-
    data_phy_NP_LRR_Erken_Spring$Mean_NP[data_phy_NP_LRR_Erken_Spring$Sno == i & data_phy_NP_LRR_Erken_Spring$Treatment == "C" & data_phy_NP_LRR_Erken_Spring$Replicate == 4]
}
for (i in unique(data_phy_NP_LRR_Bolmen_Summer$Sno)) {
  data_phy_NP_LRR_Bolmen_Summer$Mean_NP_Control[data_phy_NP_LRR_Bolmen_Summer$Sno == i] <-
    data_phy_NP_LRR_Bolmen_Summer$Mean_NP[data_phy_NP_LRR_Bolmen_Summer$Sno == i & data_phy_NP_LRR_Bolmen_Summer$Treatment == "C" & data_phy_NP_LRR_Bolmen_Summer$Replicate == 3]
}

#It somehow did not work for S10 of Bolmen Summer NP
data_phy_NP_LRR_Bolmen_Summer$Mean_NP_Control[data_phy_NP_LRR_Bolmen_Summer$Sno == 10] <- 27.564286

data_NP_LRR <- rbind(data_phy_NP_LRR_Erken_Summer, data_phy_NP_LRR_Erken_Spring, data_phy_NP_LRR_Bolmen_Summer) #Put data sets back together
data_NP_LRR <- mutate(data_NP_LRR, LRR_NP = log(NP/Mean_NP_Control))
data_NP_LRR_mean <- group_by(data_NP_LRR, Site, Season, Treatment, Sno)
data_NP_LRR_mean <- mutate(data_NP_LRR_mean, mean_LRR_NP = mean(LRR_NP), sd_LRR_NP = sd(LRR_NP))
data_NP_LRR_mean <- data_NP_LRR_mean %>% filter(ifelse(Site == "Erken", Replicate == 4, Replicate == 1))
data_NP_LRR_mean <- filter(data_NP_LRR_mean, Treatment != "C")
plot_LRR_NP <- template_plot_LRR(data_NP_LRR_mean, "mean_LRR_NP")+geom_errorbar(aes(ymin = mean_LRR_NP - sd_LRR_NP, ymax = mean_LRR_NP + sd_LRR_NP), width = 0.4, size = 0.9)
plot_LRR_NP
#ggsave("SITES2_LRR_NP.eps", plot_LRR_NP , unit = "cm", device = "eps", width = 23, height = 15, dpi = 300)

#CS (all three experiments)
data_CS2 <- data_CS1
data_CS2$Sno <- as.numeric(data_CS2$Sno)
data_CS2 <-mutate(data_CS2, Mean_CS_Control = NA)
data_phy_CS_LRR_Erken_Summer <- filter(data_CS2, Site == "Erken" & Season == "Summer")
data_phy_CS_LRR_Erken_Spring <- filter(data_CS2, Site == "Erken" & Season == "Spring")
data_phy_CS_LRR_Bolmen_Summer <- filter(data_CS2, Site == "Bolmen" & Season == "Summer")
for (i in unique(data_phy_CS_LRR_Erken_Summer$Sno)) {
  data_phy_CS_LRR_Erken_Summer$Mean_CS_Control[data_phy_CS_LRR_Erken_Summer$Sno == i] <-
    data_phy_CS_LRR_Erken_Summer$Mean_CS[data_phy_CS_LRR_Erken_Summer$Sno == i & data_phy_CS_LRR_Erken_Summer$Treatment == "C" & data_phy_CS_LRR_Erken_Summer$Replicate == 3]
}
for (i in unique(data_phy_CS_LRR_Erken_Spring$Sno)) {
  data_phy_CS_LRR_Erken_Spring$Mean_CS_Control[data_phy_CS_LRR_Erken_Spring$Sno == i] <-
    data_phy_CS_LRR_Erken_Spring$Mean_CS[data_phy_CS_LRR_Erken_Spring$Sno == i & data_phy_CS_LRR_Erken_Spring$Treatment == "C" & data_phy_CS_LRR_Erken_Spring$Replicate == 4]
}
for (i in unique(data_phy_CS_LRR_Bolmen_Summer$Sno)) {
  data_phy_CS_LRR_Bolmen_Summer$Mean_CS_Control[data_phy_CS_LRR_Bolmen_Summer$Sno == i] <-
    data_phy_CS_LRR_Bolmen_Summer$Mean_CS[data_phy_CS_LRR_Bolmen_Summer$Sno == i & data_phy_CS_LRR_Bolmen_Summer$Treatment == "C" & data_phy_CS_LRR_Bolmen_Summer$Replicate == 3]
}
data_CS_LRR <- rbind(data_phy_CS_LRR_Erken_Summer, data_phy_CS_LRR_Erken_Spring, data_phy_CS_LRR_Bolmen_Summer) #Put data sets back together
data_CS_LRR <- mutate(data_CS_LRR, LRR_CS = log(CS/Mean_CS_Control))
data_CS_LRR_mean <- group_by(data_CS_LRR, Site, Season, Treatment, Sno)
data_CS_LRR_mean <- mutate(data_CS_LRR_mean, mean_LRR_CS = mean(LRR_CS), sd_LRR_CS = sd(LRR_CS))
data_CS_LRR_mean <- data_CS_LRR_mean %>% filter(ifelse(Site == "Erken", Replicate == 4, Replicate == 1))
data_CS_LRR_mean <- filter(data_CS_LRR_mean, Treatment != "C")
plot_LRR_CS <- template_plot_LRR(data_CS_LRR_mean, "mean_LRR_CS")+geom_errorbar(aes(ymin = mean_LRR_CS - sd_LRR_CS, ymax = mean_LRR_CS + sd_LRR_CS), width = 0.4, size = 0.9)
plot_LRR_CS
#ggsave("SITES2_LRR_CS.eps", plot_LRR_CS , unit = "cm", device = "eps", width = 23, height = 15, dpi = 300)

#SN (all three experiments)
data_SN2 <- data_SN1
data_SN2$Sno <- as.numeric(data_SN2$Sno)
data_SN2 <-mutate(data_SN2, Mean_SN_Control = NA)
data_phy_SN_LRR_Erken_Summer <- filter(data_SN2, Site == "Erken" & Season == "Summer")
data_phy_SN_LRR_Erken_Spring <- filter(data_SN2, Site == "Erken" & Season == "Spring")
data_phy_SN_LRR_Bolmen_Summer <- filter(data_SN2, Site == "Bolmen" & Season == "Summer")
for (i in unique(data_phy_SN_LRR_Erken_Summer$Sno)) {
  data_phy_SN_LRR_Erken_Summer$Mean_SN_Control[data_phy_SN_LRR_Erken_Summer$Sno == i] <-
    data_phy_SN_LRR_Erken_Summer$Mean_SN[data_phy_SN_LRR_Erken_Summer$Sno == i & data_phy_SN_LRR_Erken_Summer$Treatment == "C" & data_phy_SN_LRR_Erken_Summer$Replicate == 3]
}
for (i in unique(data_phy_SN_LRR_Erken_Spring$Sno)) {
  data_phy_SN_LRR_Erken_Spring$Mean_SN_Control[data_phy_SN_LRR_Erken_Spring$Sno == i] <-
    data_phy_SN_LRR_Erken_Spring$Mean_SN[data_phy_SN_LRR_Erken_Spring$Sno == i & data_phy_SN_LRR_Erken_Spring$Treatment == "C" & data_phy_SN_LRR_Erken_Spring$Replicate == 4]
}
for (i in unique(data_phy_SN_LRR_Bolmen_Summer$Sno)) {
  data_phy_SN_LRR_Bolmen_Summer$Mean_SN_Control[data_phy_SN_LRR_Bolmen_Summer$Sno == i] <-
    data_phy_SN_LRR_Bolmen_Summer$Mean_SN[data_phy_SN_LRR_Bolmen_Summer$Sno == i & data_phy_SN_LRR_Bolmen_Summer$Treatment == "C" & data_phy_SN_LRR_Bolmen_Summer$Replicate == 3]
}

#It did not work for S10 of Bolmen Summer SN
data_phy_SN_LRR_Bolmen_Summer$Mean_SN_Control[data_phy_SN_LRR_Bolmen_Summer$Sno == 10] <- 1.3127124

data_SN_LRR <- rbind(data_phy_SN_LRR_Erken_Summer, data_phy_SN_LRR_Erken_Spring, data_phy_SN_LRR_Bolmen_Summer) #Put data sets back together
data_SN_LRR <- mutate(data_SN_LRR, LRR_SN = log(SN/Mean_SN_Control))
data_SN_LRR_mean <- group_by(data_SN_LRR, Site, Season, Treatment, Sno)
data_SN_LRR_mean <- mutate(data_SN_LRR_mean, mean_LRR_SN = mean(LRR_SN), sd_LRR_SN = sd(LRR_SN))
data_SN_LRR_mean <- data_SN_LRR_mean %>% filter(ifelse(Site == "Erken", Replicate == 4, Replicate == 1))
data_SN_LRR_mean <- filter(data_SN_LRR_mean, Treatment != "C")
plot_LRR_SN <- template_plot_LRR(data_SN_LRR_mean, "mean_LRR_SN")+geom_errorbar(aes(ymin = mean_LRR_SN - sd_LRR_SN, ymax = mean_LRR_SN + sd_LRR_SN), width = 0.4, size = 0.9)
plot_LRR_SN
#ggsave("SITES2_LRR_SN.eps", plot_LRR_SN , unit = "cm", device = "eps", width = 23, height = 15, dpi = 300)

#SP (all three experiments)
data_SP2 <- data_SP1
data_SP2$Sno <- as.numeric(data_SP2$Sno)
data_SP2 <-mutate(data_SP2, Mean_SP_Control = NA)
data_phy_SP_LRR_Erken_Summer <- filter(data_SP2, Site == "Erken" & Season == "Summer")
data_phy_SP_LRR_Erken_Spring <- filter(data_SP2, Site == "Erken" & Season == "Spring")
data_phy_SP_LRR_Bolmen_Summer <- filter(data_SP2, Site == "Bolmen" & Season == "Summer")
for (i in unique(data_phy_SP_LRR_Erken_Summer$Sno)) {
  data_phy_SP_LRR_Erken_Summer$Mean_SP_Control[data_phy_SP_LRR_Erken_Summer$Sno == i] <-
    data_phy_SP_LRR_Erken_Summer$Mean_SP[data_phy_SP_LRR_Erken_Summer$Sno == i & data_phy_SP_LRR_Erken_Summer$Treatment == "C" & data_phy_SP_LRR_Erken_Summer$Replicate == 3]
}
for (i in unique(data_phy_SP_LRR_Erken_Spring$Sno)) {
  data_phy_SP_LRR_Erken_Spring$Mean_SP_Control[data_phy_SP_LRR_Erken_Spring$Sno == i] <-
    data_phy_SP_LRR_Erken_Spring$Mean_SP[data_phy_SP_LRR_Erken_Spring$Sno == i & data_phy_SP_LRR_Erken_Spring$Treatment == "C" & data_phy_SP_LRR_Erken_Spring$Replicate == 4]
}
for (i in unique(data_phy_SP_LRR_Bolmen_Summer$Sno)) {
  data_phy_SP_LRR_Bolmen_Summer$Mean_SP_Control[data_phy_SP_LRR_Bolmen_Summer$Sno == i] <-
    data_phy_SP_LRR_Bolmen_Summer$Mean_SP[data_phy_SP_LRR_Bolmen_Summer$Sno == i & data_phy_SP_LRR_Bolmen_Summer$Treatment == "C" & data_phy_SP_LRR_Bolmen_Summer$Replicate == 3]
}
data_SP_LRR <- rbind(data_phy_SP_LRR_Erken_Summer, data_phy_SP_LRR_Erken_Spring, data_phy_SP_LRR_Bolmen_Summer) #Put data sets back together
data_SP_LRR <- mutate(data_SP_LRR, LRR_SP = log(SP/Mean_SP_Control))
data_SP_LRR_mean <- group_by(data_SP_LRR, Site, Season, Treatment, Sno)
data_SP_LRR_mean <- mutate(data_SP_LRR_mean, mean_LRR_SP = mean(LRR_SP), sd_LRR_SP = sd(LRR_SP))
data_SP_LRR_mean <- data_SP_LRR_mean %>% filter(ifelse(Site == "Erken", Replicate == 4, Replicate == 1))
data_SP_LRR_mean <- filter(data_SP_LRR_mean, Treatment != "C")
plot_LRR_SP <- template_plot_LRR(data_SP_LRR_mean, "mean_LRR_SP")+geom_errorbar(aes(ymin = mean_LRR_SP - sd_LRR_SP, ymax = mean_LRR_SP + sd_LRR_SP), width = 0.4, size = 0.9)
plot_LRR_SP
#ggsave("SITES2_LRR_SP.eps", plot_LRR_SP , unit = "cm", device = "eps", width = 23, height = 15, dpi = 300)

#### OVERALL ECOLOGICAL VULNERABILITY (OEV) ####

#Take the LRR of the parameters that I want to look at (over entire time series)
AUC_C <- data_C_LRR #From the calculations above
AUC_CN <- data_CN_LRR #From the calculations above
AUC_CP <- data_CP_LRR #From the calculations above
AUC_NP <- data_NP_LRR #From the calculations above
AUC_CS <- data_CS_LRR #From the calculations above
AUC_SN <- data_SN_LRR #From the calculations above
AUC_SP <- data_SP_LRR #From the calculations above

#Fill the gaps of missing data, otherwise the function does not work!
#I manually took the mean between the point before and after

#C
M6_R2_TreatE_S1 <- data.frame(ID = 1001, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 1, LRR_C = 0.072154238)
M7_R2_TreatC_S1 <- data.frame(ID = 1002, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 1, LRR_C = 0.072154238)
M8_R2_TreatD_S1 <- data.frame(ID = 1003, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 1, LRR_C = 0.072154238)
M9_R3_TreatE_S1 <- data.frame(ID = 1004, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 1, LRR_C = 0.072154238)
M11_R3_TreatD_S1 <- data.frame(ID = 1005, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 1, LRR_C = 0.072154238)
M12_R3_TreatI_S1 <- data.frame(ID = 1006, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 1, LRR_C = 0.072154238)
M14_R4_TreatI_S1 <- data.frame(ID = 1007, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 1, LRR_C = 0.072154238)
M15_R4_TreatE_S1 <- data.frame(ID = 1008, Site = "Bolmen", Season = "Summer", Mesocosm = 15, Replicate = 4, Treatment = "E", Sno = 1, LRR_C = 0.072154238)
M16_R4_TreatC_S1 <- data.frame(ID = 1009, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 1, LRR_C = 0.072154238)
M6_R2_TreatE_S2 <- data.frame(ID = 1010, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 2, LRR_C = 0.270104702)
M7_R2_TreatC_S2 <- data.frame(ID = 1011, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 2, LRR_C = -0.10509017)
M8_R2_TreatD_S2 <- data.frame(ID = 1012, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 2, LRR_C = 0.330845897)
M9_R3_TreatE_S2 <- data.frame(ID = 1013, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 2, LRR_C = 0.270104702)
M11_R3_TreatD_S2 <- data.frame(ID = 1014, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 2, LRR_C = 0.330845897)
M12_R3_TreatI_S2 <- data.frame(ID = 1015, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 2, LRR_C = 0.134064944)
M14_R4_TreatI_S2 <- data.frame(ID = 1016, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 2, LRR_C = 0.13406494)
M16_R4_TreatC_S2 <- data.frame(ID = 1017, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 2, LRR_C = -0.10509017)
M6_R2_TreatE_S3 <- data.frame(ID = 1018, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 3, LRR_C = (0.270104702+0.349286756193823)/2)
M7_R2_TreatC_S3 <- data.frame(ID = 1019, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 3, LRR_C = (-0.10509017+(-0.108487898234566))/2)
M8_R2_TreatD_S3 <- data.frame(ID = 1020, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 3, LRR_C = (0.330845897+0.0765154498570228)/2)
M9_R3_TreatE_S3 <- data.frame(ID = 1021, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 3, LRR_C = (0.270104702+(-0.0900347340944318))/2)
M11_R3_TreatD_S3 <- data.frame(ID = 1022, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 3, LRR_C = (0.330845897+0.249646377840533)/2)
M12_R3_TreatI_S3 <- data.frame(ID = 1023, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 3, LRR_C = (0.134064944+0.446619611148291)/2)
M14_R4_TreatI_S3 <- data.frame(ID = 1024, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 3, LRR_C = (0.13406494+0.631775713381486)/2)
M16_R4_TreatC_S3 <- data.frame(ID = 1025, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 3, LRR_C = (-0.10509017+0.376274955697022)/2)
AUC_C_full <- rbind(AUC_C, M6_R2_TreatE_S1, M7_R2_TreatC_S1, M8_R2_TreatD_S1, M9_R3_TreatE_S1,
                    M11_R3_TreatD_S1, M12_R3_TreatI_S1, M14_R4_TreatI_S1, M15_R4_TreatE_S1,
                    M16_R4_TreatC_S1, M6_R2_TreatE_S2, M7_R2_TreatC_S2, M8_R2_TreatD_S2, 
                    M9_R3_TreatE_S2, M11_R3_TreatD_S2, M12_R3_TreatI_S2, M14_R4_TreatI_S2,
                    M16_R4_TreatC_S2, M6_R2_TreatE_S3, M7_R2_TreatC_S3, M8_R2_TreatD_S3,
                    M9_R3_TreatE_S3, M11_R3_TreatD_S3, M12_R3_TreatI_S3, M14_R4_TreatI_S3,
                    M16_R4_TreatC_S3)

#CN
M6_R2_TreatE_S1 <- data.frame(ID = 1001, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 1, LRR_CN = 0.005749165)
M7_R2_TreatC_S1 <- data.frame(ID = 1002, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 1, LRR_CN = 0.005749165)
M8_R2_TreatD_S1 <- data.frame(ID = 1003, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 1, LRR_CN = 0.005749165)
M9_R3_TreatE_S1 <- data.frame(ID = 1004, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 1, LRR_CN = 0.005749165)
M11_R3_TreatD_S1 <- data.frame(ID = 1005, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 1, LRR_CN = 0.005749165)
M12_R3_TreatI_S1 <- data.frame(ID = 1006, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 1, LRR_CN = 0.005749165)
M14_R4_TreatI_S1 <- data.frame(ID = 1007, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 1, LRR_CN = 0.005749165)
M15_R4_TreatE_S1 <- data.frame(ID = 1008, Site = "Bolmen", Season = "Summer", Mesocosm = 15, Replicate = 4, Treatment = "E", Sno = 1, LRR_CN = 0.005749165)
M16_R4_TreatC_S1 <- data.frame(ID = 1009, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 1, LRR_CN = 0.005749165)
M6_R2_TreatE_S2 <- data.frame(ID = 1010, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 2, LRR_CN = -0.258668635)
M7_R2_TreatC_S2 <- data.frame(ID = 1011, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 2, LRR_CN = -0.027680527)
M8_R2_TreatD_S2 <- data.frame(ID = 1012, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 2, LRR_CN = -0.236585891)
M9_R3_TreatE_S2 <- data.frame(ID = 1013, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 2, LRR_CN = -0.258668635)
M11_R3_TreatD_S2 <- data.frame(ID = 1014, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 2, LRR_CN = -0.236585891)
M12_R3_TreatI_S2 <- data.frame(ID = 1015, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 2, LRR_CN = -0.225550556)
M14_R4_TreatI_S2 <- data.frame(ID = 1016, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 2, LRR_CN = -0.225550556)
M16_R4_TreatC_S2 <- data.frame(ID = 1017, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 2, LRR_CN = -0.027680527)
M6_R2_TreatE_S3 <- data.frame(ID = 1018, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 3, LRR_CN = (-0.2586686352+(-0.0940307529898581))/2)
M7_R2_TreatC_S3 <- data.frame(ID = 1019, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 3, LRR_CN = (-0.027680527+(-0.0217138795311652))/2)
M8_R2_TreatD_S3 <- data.frame(ID = 1020, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 3, LRR_CN = (-0.236585891+(-0.139270869756275))/2)
M9_R3_TreatE_S3 <- data.frame(ID = 1021, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 3, LRR_CN = (-0.258668635+(0.0871274736784992))/2)
M11_R3_TreatD_S3 <- data.frame(ID = 1022, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 3, LRR_CN = (-0.236585891+(-0.179758584760324))/2)
M12_R3_TreatI_S3 <- data.frame(ID = 1023, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 3, LRR_CN = (-0.225550556+(-0.081606027903596))/2)
M14_R4_TreatI_S3 <- data.frame(ID = 1024, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 3, LRR_CN = (-0.225550556+(-0.268196136529427))/2)
M16_R4_TreatC_S3 <- data.frame(ID = 1025, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 3, LRR_CN = (-0.027680527+(-0.0830617953504457))/2)
M15_R4_TreatE_S5_ES <- data.frame(ID = 1026, Site = "Erken", Season = "Spring", Mesocosm = 15, Replicate = 4, Treatment = "E", Sno = 5, LRR_CN = -0.105851717)
M15_R4_TreatE_S8_BS <- data.frame(ID = 1027, Site = "Bolmen", Season = "Summer", Mesocosm = 15, Replicate = 4, Treatment = "E", Sno = 8, LRR_CN = (0.0143920788+(-0.2186996953))/2)
M10_R3_TreatC_S10_BS <- data.frame(ID = 1028, Site = "Bolmen", Season = "Summer", Mesocosm = 10, Replicate = 3, Treatment = "C", Sno = 10, LRR_CN = 0.2774227959)

#M15_R4_TreatE_S9_BS <- data.frame(ID = 1027, Site = "Bolmen", Season = "Summer", Mesocosm = 15, Replicate = 4, Treatment = "E", Sno = 9, LRR_CN = (0.0143920788+0.1380045285)/2)
#M5_R2_TreatI_S8_ES <- data.frame(ID = 1027, Site = "Erken", Season = "Spring", Mesocosm = 5, Replicate = 2, Treatment = "I", Sno = 8, LRR_CN = -0.050961016)
#M10_R3_TreatC_S8_ES <- data.frame(ID = 1028, Site = "Erken", Season = "Spring", Mesocosm = 10, Replicate = 3, Treatment = "C", Sno = 8, LRR_CN = -0.005821125)

AUC_CN_full <- rbind(AUC_CN, M6_R2_TreatE_S1, M7_R2_TreatC_S1, M8_R2_TreatD_S1, M9_R3_TreatE_S1,
                    M11_R3_TreatD_S1, M12_R3_TreatI_S1, M14_R4_TreatI_S1, M15_R4_TreatE_S1,
                    M16_R4_TreatC_S1, M6_R2_TreatE_S2, M7_R2_TreatC_S2, M8_R2_TreatD_S2, 
                    M9_R3_TreatE_S2, M11_R3_TreatD_S2, M12_R3_TreatI_S2, M14_R4_TreatI_S2,
                    M16_R4_TreatC_S2, M6_R2_TreatE_S3, M7_R2_TreatC_S3, M8_R2_TreatD_S3,
                    M9_R3_TreatE_S3, M11_R3_TreatD_S3, M12_R3_TreatI_S3, M14_R4_TreatI_S3,
                    M16_R4_TreatC_S3)

AUC_CN_full <- rbind(AUC_CN_full, M15_R4_TreatE_S5_ES, M15_R4_TreatE_S8_BS, M10_R3_TreatC_S10_BS)

#CP
M6_R2_TreatE_S1 <- data.frame(ID = 1001, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 1, LRR_CP = -0.20692206)
M7_R2_TreatC_S1 <- data.frame(ID = 1002, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 1, LRR_CP = -0.20692206)
M8_R2_TreatD_S1 <- data.frame(ID = 1003, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 1, LRR_CP = -0.20692206)
M9_R3_TreatE_S1 <- data.frame(ID = 1004, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 1, LRR_CP = -0.20692206)
M11_R3_TreatD_S1 <- data.frame(ID = 1005, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 1, LRR_CP = -0.20692206)
M12_R3_TreatI_S1 <- data.frame(ID = 1006, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 1, LRR_CP = -0.20692206)
M14_R4_TreatI_S1 <- data.frame(ID = 1007, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 1, LRR_CP = -0.20692206)
M15_R4_TreatE_S1 <- data.frame(ID = 1008, Site = "Bolmen", Season = "Summer", Mesocosm = 15, Replicate = 4, Treatment = "E", Sno = 1, LRR_CP = -0.20692206)
M16_R4_TreatC_S1 <- data.frame(ID = 1009, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 1, LRR_CP = -0.20692206)
M6_R2_TreatE_S2 <- data.frame(ID = 1010, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 2, LRR_CP = -0.03492518)
M7_R2_TreatC_S2 <- data.frame(ID = 1011, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 2, LRR_CP = -0.004471603)
M8_R2_TreatD_S2 <- data.frame(ID = 1012, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 2, LRR_CP = -0.449930592)
M9_R3_TreatE_S2 <- data.frame(ID = 1013, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 2, LRR_CP = -0.03492518)
M11_R3_TreatD_S2 <- data.frame(ID = 1014, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 2, LRR_CP = -0.449930592)
M12_R3_TreatI_S2 <- data.frame(ID = 1015, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 2, LRR_CP = -0.430550321)
M14_R4_TreatI_S2 <- data.frame(ID = 1016, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 2, LRR_CP = -0.430550321)
M16_R4_TreatC_S2 <- data.frame(ID = 1017, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 2, LRR_CP = -0.004471603)
M6_R2_TreatE_S3 <- data.frame(ID = 1018, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 3, LRR_CP = (-0.03492518+(-0.880412123697145))/2)
M7_R2_TreatC_S3 <- data.frame(ID = 1019, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 3, LRR_CP = (-0.004471603+(-0.121867149760085))/2)
M8_R2_TreatD_S3 <- data.frame(ID = 1020, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 3, LRR_CP = (-0.449930592+(-0.467145558541313))/2)
M9_R3_TreatE_S3 <- data.frame(ID = 1021, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 3, LRR_CP = (-0.03492518+(-0.604561895807426))/2)
M11_R3_TreatD_S3 <- data.frame(ID = 1022, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 3, LRR_CP = (-0.449930592+(-0.53897000075762))/2)
M12_R3_TreatI_S3 <- data.frame(ID = 1023, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 3, LRR_CP = (-0.430550321+(-0.615569656831044))/2)
M14_R4_TreatI_S3 <- data.frame(ID = 1024, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 3, LRR_CP = (-0.430550321+(-0.786267048082095))/2)
M16_R4_TreatC_S3 <- data.frame(ID = 1025, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 3, LRR_CP = (-0.004471603+(0.330882879521312))/2)
M15_R4_TreatE_S5_ES <- data.frame(ID = 1026, Site = "Erken", Season = "Spring", Mesocosm = 15, Replicate = 4, Treatment = "E", Sno = 5, LRR_CP = -0.276761392)
M15_R2_TreatE_S3_ESu <- data.frame(ID = 1029, Site = "Erken", Season = "Summer", Mesocosm = 15, Replicate = 4, Treatment = "E", Sno = 3, LRR_CP = -0.279057431)
M16_R3_TreatC_S3_ESu <- data.frame(ID = 1030, Site = "Erken", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 3, LRR_CP = -0.11894558)
M12_R3_TreatI_S8_ESu <- data.frame(ID = 1031, Site = "Erken", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 8, LRR_CP = -0.413968508)
M1_R1_TreatC_S1 <- data.frame(ID = 1033, Site = "Bolmen", Season = "Summer", Mesocosm = 1, Replicate = 1, Treatment = "C", Sno = 1, LRR_CP = -0.20692206)

AUC_CP_full <- rbind(AUC_CP, M6_R2_TreatE_S1, M7_R2_TreatC_S1, M8_R2_TreatD_S1, M9_R3_TreatE_S1,
                     M11_R3_TreatD_S1, M12_R3_TreatI_S1, M14_R4_TreatI_S1, M15_R4_TreatE_S1,
                     M16_R4_TreatC_S1, M6_R2_TreatE_S2, M7_R2_TreatC_S2, M8_R2_TreatD_S2, 
                     M9_R3_TreatE_S2, M11_R3_TreatD_S2, M12_R3_TreatI_S2, M14_R4_TreatI_S2,
                     M16_R4_TreatC_S2, M6_R2_TreatE_S3, M7_R2_TreatC_S3, M8_R2_TreatD_S3,
                     M9_R3_TreatE_S3, M11_R3_TreatD_S3, M12_R3_TreatI_S3, M14_R4_TreatI_S3,
                     M16_R4_TreatC_S3, M15_R4_TreatE_S5_ES, M15_R2_TreatE_S3_ESu, M16_R3_TreatC_S3_ESu,
                     M12_R3_TreatI_S8_ESu, M1_R1_TreatC_S1)

#NP
M6_R2_TreatE_S1 <- data.frame(ID = 1001, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 1, LRR_NP = -0.171550938)
M7_R2_TreatC_S1 <- data.frame(ID = 1002, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 1, LRR_NP = -0.171550938)
M8_R2_TreatD_S1 <- data.frame(ID = 1003, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 1, LRR_NP = -0.171550938)
M9_R3_TreatE_S1 <- data.frame(ID = 1004, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 1, LRR_NP = -0.171550938)
M11_R3_TreatD_S1 <- data.frame(ID = 1005, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 1, LRR_NP = -0.171550938)
M12_R3_TreatI_S1 <- data.frame(ID = 1006, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 1, LRR_NP = -0.171550938)
M14_R4_TreatI_S1 <- data.frame(ID = 1007, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 1, LRR_NP = -0.171550938)
#M15_R4_TreatE_S1 <- data.frame(ID = 1008, Site = "Bolmen", Season = "Summer", Mesocosm = 15, Replicate = 4, Treatment = "E", Sno = 1, LRR_NP = -0.171550938)
M16_R4_TreatC_S1 <- data.frame(ID = 1009, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 1, LRR_NP = -0.171550938)
M6_R2_TreatE_S2 <- data.frame(ID = 1010, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 2, LRR_NP = 0.211776422)
M7_R2_TreatC_S2 <- data.frame(ID = 1011, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 2, LRR_NP = -0.002598373)
M8_R2_TreatD_S2 <- data.frame(ID = 1012, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 2, LRR_NP = -0.225311735)
M9_R3_TreatE_S2 <- data.frame(ID = 1013, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 2, LRR_NP = 0.211776422)
M11_R3_TreatD_S2 <- data.frame(ID = 1014, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 2, LRR_NP = -0.225311735)
M12_R3_TreatI_S2 <- data.frame(ID = 1015, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 2, LRR_NP = -0.2169668)
M14_R4_TreatI_S2 <- data.frame(ID = 1016, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 2, LRR_NP = -0.2169668)
M16_R4_TreatC_S2 <- data.frame(ID = 1017, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 2, LRR_NP = -0.002598373)
M6_R2_TreatE_S3 <- data.frame(ID = 1018, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 3, LRR_NP = (0.211776422+(-0.801336545970006))/2)
M7_R2_TreatC_S3 <- data.frame(ID = 1019, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 3, LRR_NP = (-0.002598373+(-0.115108445491639))/2)
M8_R2_TreatD_S3 <- data.frame(ID = 1020, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 3, LRR_NP = (-0.225311735+(-0.342829864047757))/2)
M9_R3_TreatE_S3 <- data.frame(ID = 1021, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 3, LRR_NP = (0.211776422+(-0.706644544748644))/2)
M11_R3_TreatD_S3 <- data.frame(ID = 1022, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 3, LRR_NP = (-0.225311735+(-0.374166591260015))/2)
M12_R3_TreatI_S3 <- data.frame(ID = 1023, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 3, LRR_NP = (-0.2169668+(-0.548918804190168))/2)
M14_R4_TreatI_S3 <- data.frame(ID = 1024, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 3, LRR_NP = (-0.2169668+(-0.533026086815388))/2)
M16_R4_TreatC_S3 <- data.frame(ID = 1025, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 3, LRR_NP = (-0.002598373+(0.398989499609038))/2)
M6_R2_TreatE_S6 <- data.frame(ID = 1032, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 6, LRR_NP = -0.648283351)
M15_R2_TreatE_S8 <- data.frame(ID = 1033, Site = "Bolmen", Season = "Summer", Mesocosm = 15, Replicate = 4, Treatment = "E", Sno = 8, LRR_NP = -0.423952027)
M15_R2_TreatE_S3_ESu <- data.frame(ID = 1029, Site = "Erken", Season = "Summer", Mesocosm = 15, Replicate = 4, Treatment = "E", Sno = 3, LRR_NP = -0.267884766)
M16_R3_TreatC_S3_ESu <- data.frame(ID = 1030, Site = "Erken", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 3, LRR_NP = -0.133706248)
M12_R3_TreatI_S8_ESu <- data.frame(ID = 1031, Site = "Erken", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 8, LRR_NP = -0.288062497)
M15_R4_TreatE_S5_ES <- data.frame(ID = 1026, Site = "Erken", Season = "Spring", Mesocosm = 15, Replicate = 4, Treatment = "E", Sno = 5, LRR_NP = -0.185608958)
M1_R1_TreatC_S1 <- data.frame(ID = 1033, Site = "Bolmen", Season = "Summer", Mesocosm = 1, Replicate = 1, Treatment = "C", Sno = 1, LRR_NP = -0.171550938)

AUC_NP_full <- rbind(AUC_NP, M6_R2_TreatE_S1, M7_R2_TreatC_S1, M8_R2_TreatD_S1, M9_R3_TreatE_S1,
                     M11_R3_TreatD_S1, M12_R3_TreatI_S1, M14_R4_TreatI_S1,
                     M16_R4_TreatC_S1, M6_R2_TreatE_S2, M7_R2_TreatC_S2, M8_R2_TreatD_S2, 
                     M9_R3_TreatE_S2, M11_R3_TreatD_S2, M12_R3_TreatI_S2, M14_R4_TreatI_S2,
                     M16_R4_TreatC_S2, M6_R2_TreatE_S3, M7_R2_TreatC_S3, M8_R2_TreatD_S3,
                     M9_R3_TreatE_S3, M11_R3_TreatD_S3, M12_R3_TreatI_S3, M14_R4_TreatI_S3,
                     M16_R4_TreatC_S3, M6_R2_TreatE_S6, M15_R2_TreatE_S8, M15_R2_TreatE_S3_ESu, 
                     M16_R3_TreatC_S3_ESu, M12_R3_TreatI_S8_ESu, M15_R4_TreatE_S5_ES, M1_R1_TreatC_S1)

#CS 

M15_R4_TreatE_S5_ES <- data.frame(ID = 1026, Site = "Erken", Season = "Spring", Mesocosm = 15, Replicate = 4, Treatment = "E", Sno = 5, LRR_CS = (0.141155696+0.401531323)/2)
M6_R2_TreatE_S1 <- data.frame(ID = 1001, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 1, LRR_CS = 0.07068658857)
M7_R2_TreatC_S1 <- data.frame(ID = 1002, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 1, LRR_CS = 0.07068658857)
M8_R2_TreatD_S1 <- data.frame(ID = 1003, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 1, LRR_CS = 0.07068658857)
M9_R3_TreatE_S1 <- data.frame(ID = 1004, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 1, LRR_CS = 0.07068658857)
M11_R3_TreatD_S1 <- data.frame(ID = 1005, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 1, LRR_CS = 0.07068658857)
M12_R3_TreatI_S1 <- data.frame(ID = 1006, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 1, LRR_CS = 0.07068658857)
M14_R4_TreatI_S1 <- data.frame(ID = 1007, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 1, LRR_CS = 0.07068658857)
M15_R4_TreatE_S1 <- data.frame(ID = 1008, Site = "Bolmen", Season = "Summer", Mesocosm = 15, Replicate = 4, Treatment = "E", Sno = 1, LRR_CS = 0.07068658857)
M16_R4_TreatC_S1 <- data.frame(ID = 1009, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 1, LRR_CS = 0.07068658857)

M6_R2_TreatE_S2 <- data.frame(ID = 1010, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 2, LRR_CS = (0.31419189+0.19517600)/2)
M7_R2_TreatC_S2 <- data.frame(ID = 1011, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 2, LRR_CS = (0.17381668+(-0.21052033))/2)
M8_R2_TreatD_S2 <- data.frame(ID = 1012, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 2, LRR_CS = (0.34352857+0.00765856)/2)
M9_R3_TreatE_S2 <- data.frame(ID = 1013, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 2, LRR_CS = (0.31419189+0.19517600)/2)
M11_R3_TreatD_S2 <- data.frame(ID = 1014, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 2, LRR_CS = (0.34352857+0.00765856)/2)
M12_R3_TreatI_S2 <- data.frame(ID = 1015, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 2, LRR_CS = (0.09730095+(-0.01137597))/2)
M14_R4_TreatI_S2 <- data.frame(ID = 1016, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 2, LRR_CS = (0.09730095+(-0.01137597))/2)
M16_R4_TreatC_S2 <- data.frame(ID = 1017, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 2, LRR_CS = (0.17381668+(-0.21052033))/2)

M6_R2_TreatE_S3 <- data.frame(ID = 1018, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 3, LRR_CS = (0.2546839+(-0.06786466))/2)
M7_R2_TreatC_S3 <- data.frame(ID = 1019, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 3, LRR_CS = ((-0.26489552)+(-0.01835183))/2)
M8_R2_TreatD_S3 <- data.frame(ID = 1020, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 3, LRR_CS = ((-0.21394647)+0.1755936)/2)
M9_R3_TreatE_S3 <- data.frame(ID = 1021, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 3, LRR_CS = (0.01246597+0.2546839)/2)
M11_R3_TreatD_S3 <- data.frame(ID = 1022, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 3, LRR_CS = ((-0.04682543)+0.1755936)/2)
M12_R3_TreatI_S3 <- data.frame(ID = 1023, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 3, LRR_CS = ((-0.18989850)+0.04296249)/2)
M14_R4_TreatI_S3 <- data.frame(ID = 1024, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 3, LRR_CS = (0.2037305+0.04296249)/2)
M16_R4_TreatC_S3 <- data.frame(ID = 1025, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 3, LRR_CS = (0.35012469+(-0.01835183))/2)

M15_R4_TreatE_S8_BS <- data.frame(ID = 443, Site = "Bolmen", Season = "Summer", Mesocosm = 15, Replicate = 4, Treatment = "E", Sno = 8, LRR_CS = ((-0.26104931)+(-0.39697053))/2)

AUC_CS_full <- rbind(AUC_CS, M6_R2_TreatE_S1, M7_R2_TreatC_S1, M8_R2_TreatD_S1, M9_R3_TreatE_S1,
                     M11_R3_TreatD_S1, M12_R3_TreatI_S1, M14_R4_TreatI_S1, M15_R4_TreatE_S1,
                     M16_R4_TreatC_S1, M6_R2_TreatE_S2, M7_R2_TreatC_S2, M8_R2_TreatD_S2, 
                     M9_R3_TreatE_S2, M11_R3_TreatD_S2, M12_R3_TreatI_S2, M14_R4_TreatI_S2,
                     M16_R4_TreatC_S2, M6_R2_TreatE_S3, M7_R2_TreatC_S3, M8_R2_TreatD_S3,
                     M9_R3_TreatE_S3, M11_R3_TreatD_S3, M12_R3_TreatI_S3, M14_R4_TreatI_S3,
                     M16_R4_TreatC_S3, M15_R4_TreatE_S8_BS, M15_R4_TreatE_S5_ES)

#SN

M15_R4_TreatE_S5_ES <- data.frame(ID = 1026, Site = "Erken", Season = "Spring", Mesocosm = 15, Replicate = 4, Treatment = "E", Sno = 5, LRR_SN = ((-0.2157070)+(-0.5877709))/2)

M6_R2_TreatE_S1 <- data.frame(ID = 1001, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 1, LRR_SN = -0.0855397)
M7_R2_TreatC_S1 <- data.frame(ID = 1002, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 1, LRR_SN = -0.0855397)
M8_R2_TreatD_S1 <- data.frame(ID = 1003, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 1, LRR_SN = -0.0855397)
M9_R3_TreatE_S1 <- data.frame(ID = 1004, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 1, LRR_SN = -0.0855397)
M11_R3_TreatD_S1 <- data.frame(ID = 1005, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 1, LRR_SN = -0.0855397)
M12_R3_TreatI_S1 <- data.frame(ID = 1006, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 1, LRR_SN = -0.0855397)
M14_R4_TreatI_S1 <- data.frame(ID = 1007, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 1, LRR_SN = -0.0855397)
M16_R4_TreatC_S1 <- data.frame(ID = 1009, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 1, LRR_SN = -0.0855397)

M6_R2_TreatE_S2 <- data.frame(ID = 1010, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 2, LRR_SN = -0.58094118)
M7_R2_TreatC_S2 <- data.frame(ID = 1011, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 2, LRR_SN = -0.06307703)
M8_R2_TreatD_S2 <- data.frame(ID = 1012, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 2, LRR_SN = -0.47976805)
M9_R3_TreatE_S2 <- data.frame(ID = 1013, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 2, LRR_SN = -0.58094118)
M11_R3_TreatD_S2 <- data.frame(ID = 1014, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 2, LRR_SN = -0.47976805)
M12_R3_TreatI_S2 <- data.frame(ID = 1015, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 2, LRR_SN = -0.33610164)
M14_R4_TreatI_S2 <- data.frame(ID = 1016, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 2, LRR_SN = -0.33610164)
M16_R4_TreatC_S2 <- data.frame(ID = 1017, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 2, LRR_SN = -0.06307703)

M6_R2_TreatE_S3 <- data.frame(ID = 1018, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 3, LRR_SN = ((-0.58094118)+(-0.08896122))/2)
M7_R2_TreatC_S3 <- data.frame(ID = 1019, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 3, LRR_SN = ((-0.06307703)+(0.18038652))/2)
M8_R2_TreatD_S3 <- data.frame(ID = 1020, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 3, LRR_SN = ((-0.47976805)+0.01188048)/2)
M9_R3_TreatE_S3 <- data.frame(ID = 1021, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 3, LRR_SN = ((-0.58094118)+0.01186638)/2)
M11_R3_TreatD_S3 <- data.frame(ID = 1022, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 3, LRR_SN = ((-0.47976805)+(-0.19572828))/2)
M12_R3_TreatI_S3 <- data.frame(ID = 1023, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 3, LRR_SN = ((-0.33610164)+0.045497349)/2)
M14_R4_TreatI_S3 <- data.frame(ID = 1024, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 3, LRR_SN = ((0.33610164)+(-0.5347218))/2)
M16_R4_TreatC_S3 <- data.frame(ID = 1025, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 3, LRR_SN = ((-0.06307703)+(-0.49598161))/2)

M15_R4_TreatE_S8_BS <- data.frame(ID = 443, Site = "Bolmen", Season = "Summer", Mesocosm = 15, Replicate = 4, Treatment = "E", Sno = 8, LRR_SN = ((0.25394362)+(-0.02442311))/2)
M10_R4_TreatE_S10_BS <- data.frame(ID = 471, Site = "Bolmen", Season = "Summer", Mesocosm = 10, Replicate = 3, Treatment = "C", Sno = 10, LRR_SN = -0.15394054)


AUC_SN_full <- rbind(AUC_SN, M6_R2_TreatE_S1, M7_R2_TreatC_S1, M8_R2_TreatD_S1, M9_R3_TreatE_S1,
                     M11_R3_TreatD_S1, M12_R3_TreatI_S1, M14_R4_TreatI_S1,
                     M16_R4_TreatC_S1, M6_R2_TreatE_S2, M7_R2_TreatC_S2, M8_R2_TreatD_S2, 
                     M9_R3_TreatE_S2, M11_R3_TreatD_S2, M12_R3_TreatI_S2, M14_R4_TreatI_S2,
                     M16_R4_TreatC_S2, M6_R2_TreatE_S3, M7_R2_TreatC_S3, M8_R2_TreatD_S3,
                     M9_R3_TreatE_S3, M11_R3_TreatD_S3, M12_R3_TreatI_S3, M14_R4_TreatI_S3,
                     M16_R4_TreatC_S3, M15_R4_TreatE_S8_BS, M10_R4_TreatE_S10_BS, M15_R4_TreatE_S5_ES)

#SP

M16_R4_TreatC_S3_ES <- data.frame(ID = 49, Site = "Erken", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 3, LRR_SP = ((-0.061122558)+(-0.111085221))/2)
M12_R3_TreatI_S8_ES <- data.frame(ID = 128, Site = "Erken", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 8, LRR_SP = ((-0.98554725)+(-0.31173503))/2)

M6_R2_TreatE_S1 <- data.frame(ID = 1001, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 1, LRR_SP = -0.9401808)
M7_R2_TreatC_S1 <- data.frame(ID = 1002, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 1, LRR_SP = -0.9401808)
M8_R2_TreatD_S1 <- data.frame(ID = 1003, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 1, LRR_SP = -0.9401808)
M9_R3_TreatE_S1 <- data.frame(ID = 1004, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 1, LRR_SP = -0.9401808)
M11_R3_TreatD_S1 <- data.frame(ID = 1005, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 1, LRR_SP = -0.9401808)
M12_R3_TreatI_S1 <- data.frame(ID = 1006, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 1, LRR_SP = -0.9401808)
M14_R4_TreatI_S1 <- data.frame(ID = 1007, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 1, LRR_SP = -0.9401808)
M16_R4_TreatC_S1 <- data.frame(ID = 1009, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 1, LRR_SP = -0.9401808)

M6_R2_TreatE_S2 <- data.frame(ID = 1010, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 2, LRR_SP = -0.34406682)
M7_R2_TreatC_S2 <- data.frame(ID = 1011, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 2, LRR_SP = -0.04057747)
M8_R2_TreatD_S2 <- data.frame(ID = 1012, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 2, LRR_SP = -0.67998185)
M9_R3_TreatE_S2 <- data.frame(ID = 1013, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 2, LRR_SP = -0.34406682)
M11_R3_TreatD_S2 <- data.frame(ID = 1014, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 2, LRR_SP = -0.67998185)
M12_R3_TreatI_S2 <- data.frame(ID = 1015, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 2, LRR_SP = -0.52797051)
M14_R4_TreatI_S2 <- data.frame(ID = 1016, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 2, LRR_SP = -0.52797051)
M16_R4_TreatC_S2 <- data.frame(ID = 1017, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 2, LRR_SP = -0.04057747)

M6_R2_TreatE_S3 <- data.frame(ID = 1018, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 3, LRR_SP = ((-0.82824649)+(-0.34406682))/2)
M7_R2_TreatC_S3 <- data.frame(ID = 1019, Site = "Bolmen", Season = "Summer", Mesocosm = 7, Replicate = 2, Treatment = "C", Sno = 3, LRR_SP = ((0.12732935)+(-0.04057747))/2)
M8_R2_TreatD_S3 <- data.frame(ID = 1020, Site = "Bolmen", Season = "Summer", Mesocosm = 8, Replicate = 2, Treatment = "D", Sno = 3, LRR_SP = ((-0.26889812)+(-0.67998185))/2)
M9_R3_TreatE_S3 <- data.frame(ID = 1021, Site = "Bolmen", Season = "Summer", Mesocosm = 9, Replicate = 3, Treatment = "E", Sno = 3, LRR_SP = ((-0.63272689)+(-0.34406682))/2)
M11_R3_TreatD_S3 <- data.frame(ID = 1022, Site = "Bolmen", Season = "Summer", Mesocosm = 11, Replicate = 3, Treatment = "D", Sno = 3, LRR_SP = ((-0.50784360)+(-0.67998185))/2)
M12_R3_TreatI_S3 <- data.frame(ID = 1023, Site = "Bolmen", Season = "Summer", Mesocosm = 12, Replicate = 3, Treatment = "I", Sno = 3, LRR_SP = ((-0.44137018)+(-0.52797051))/2)
M14_R4_TreatI_S3 <- data.frame(ID = 1024, Site = "Bolmen", Season = "Summer", Mesocosm = 14, Replicate = 4, Treatment = "I", Sno = 3, LRR_SP = ((-1.00569661)+(-0.52797051))/2)
M16_R4_TreatC_S3 <- data.frame(ID = 1025, Site = "Bolmen", Season = "Summer", Mesocosm = 16, Replicate = 4, Treatment = "C", Sno = 3, LRR_SP = ((-1.16317140)+(-0.04057747))/2)

M6_R2_TreatE_S6_BS <- data.frame(ID = 401, Site = "Bolmen", Season = "Summer", Mesocosm = 6, Replicate = 2, Treatment = "E", Sno = 6, LRR_SP = ((-0.5437692)+(-0.8635829))/2)

AUC_SP_full <- rbind(AUC_SP, M6_R2_TreatE_S1, M7_R2_TreatC_S1, M8_R2_TreatD_S1, M9_R3_TreatE_S1,
                     M11_R3_TreatD_S1, M12_R3_TreatI_S1, M14_R4_TreatI_S1,
                     M16_R4_TreatC_S1, M6_R2_TreatE_S2, M7_R2_TreatC_S2, M8_R2_TreatD_S2, 
                     M9_R3_TreatE_S2, M11_R3_TreatD_S2, M12_R3_TreatI_S2, M14_R4_TreatI_S2,
                     M16_R4_TreatC_S2, M6_R2_TreatE_S3, M7_R2_TreatC_S3, M8_R2_TreatD_S3,
                     M9_R3_TreatE_S3, M11_R3_TreatD_S3, M12_R3_TreatI_S3, M14_R4_TreatI_S3,
                     M16_R4_TreatC_S3, M16_R4_TreatC_S3_ES, M12_R3_TreatI_S8_ES, M6_R2_TreatE_S6_BS)

AUC_C_full1 <- AUC_C_full %>% mutate(mesocosmID = paste(Site, Season, Mesocosm, Treatment, Replicate, sep = "_"))
AUC_C_full1 <- AUC_C_full1[with(AUC_C_full1, order(mesocosmID, Sno)),]

AUC_CN_full1 <- AUC_CN_full %>% mutate(mesocosmID = paste(Site, Season, Mesocosm, Treatment, Replicate, sep = "_"))
AUC_CN_full1 <- AUC_CN_full1[with(AUC_CN_full1, order(mesocosmID, Sno)),]

AUC_CP_full1 <- AUC_CP_full %>% mutate(mesocosmID = paste(Site, Season, Mesocosm, Treatment, Replicate, sep = "_"))
AUC_CP_full1 <- AUC_CP_full1[with(AUC_CP_full1, order(mesocosmID, Sno)),]

AUC_NP_full1 <- AUC_NP_full %>% mutate(mesocosmID = paste(Site, Season, Mesocosm, Treatment, Replicate, sep = "_"))
AUC_NP_full1 <- AUC_NP_full1[with(AUC_NP_full1, order(mesocosmID, Sno)),]

AUC_CS_full1 <- AUC_CS_full %>% mutate(mesocosmID = paste(Site, Season, Mesocosm, Treatment, Replicate, sep = "_"))
AUC_CS_full1 <- AUC_CS_full1[with(AUC_CS_full1, order(mesocosmID, Sno)),]

AUC_SN_full1 <- AUC_SN_full %>% mutate(mesocosmID = paste(Site, Season, Mesocosm, Treatment, Replicate, sep = "_"))
AUC_SN_full1 <- AUC_SN_full1[with(AUC_SN_full1, order(mesocosmID, Sno)),]

AUC_SP_full1 <- AUC_SP_full %>% mutate(mesocosmID = paste(Site, Season, Mesocosm, Treatment, Replicate, sep = "_"))
AUC_SP_full1 <- AUC_SP_full1[with(AUC_SP_full1, order(mesocosmID, Sno)),]

#The first sampling point needs to be 0 in order to use the AUC function
AUC_C_full1$Sno <- AUC_C_full1$Sno-1
AUC_CN_full1$Sno <- AUC_CN_full1$Sno-1
AUC_CP_full1$Sno <- AUC_CP_full1$Sno-1
AUC_NP_full1$Sno <- AUC_NP_full1$Sno-1
AUC_CS_full1$Sno <- AUC_CS_full1$Sno-1
AUC_SN_full1$Sno <- AUC_SN_full1$Sno-1
AUC_SP_full1$Sno <- AUC_SP_full1$Sno-1

#Create an USI (unique identifier for the experimental unit)
USI_C <- unique(AUC_C_full1$mesocosmID)
USI_CN <- unique(AUC_CN_full1$mesocosmID)
USI_CP <- unique(AUC_CP_full1$mesocosmID)
USI_NP <- unique(AUC_NP_full1$mesocosmID)
USI_CS <- unique(AUC_CS_full1$mesocosmID)
USI_SN <- unique(AUC_SN_full1$mesocosmID)
USI_SP <- unique(AUC_SP_full1$mesocosmID)

#Create empty data frame for AUC
area_C <- data.frame()
area_CN <- data.frame()
area_CP <- data.frame()
area_NP <- data.frame()
area_CS <- data.frame()
area_SN <- data.frame()
area_SP <- data.frame()

#Conduct a loop for calculating AUC
for(i in 1:length(USI_C)){
  temp <- AUC_C_full1[AUC_C_full1$mesocosmID==USI_C[i], ] #temp is a temporary data frame
  area_value <- pk.calc.auc(abs(temp$LRR_C),temp$Sno)
  area_C <- rbind(area_C,data.frame(temp, area_value))
} #For C
for(i in 1:length(USI_CN)){
  temp <- AUC_CN_full1[AUC_CN_full1$mesocosmID==USI_CN[i], ] #temp is a temporary data frame
  area_value <- pk.calc.auc(abs(temp$LRR_CN),temp$Sno)
  area_CN <- rbind(area_CN,data.frame(temp, area_value))
} #For CN
for(i in 1:length(USI_CP)){
  temp <- AUC_CP_full1[AUC_CP_full1$mesocosmID==USI_CP[i], ] #temp is a temporary data frame
  area_value <- pk.calc.auc(abs(temp$LRR_CP),temp$Sno)
  area_CP <- rbind(area_CP,data.frame(temp, area_value))
} #For CP
for(i in 1:length(USI_NP)){
  temp <- AUC_NP_full1[AUC_NP_full1$mesocosmID==USI_NP[i], ] #temp is a temporary data frame
  area_value <- pk.calc.auc(abs(temp$LRR_NP),temp$Sno)
  area_NP <- rbind(area_NP,data.frame(temp, area_value))
} #For NP
for(i in 1:length(USI_CS)){
  temp <- AUC_CS_full1[AUC_CS_full1$mesocosmID==USI_CS[i], ] #temp is a temporary data frame
  area_value <- pk.calc.auc(abs(temp$LRR_CS),temp$Sno)
  area_CS <- rbind(area_CS,data.frame(temp, area_value))
} #For CS
for(i in 1:length(USI_SN)){
  temp <- AUC_SN_full1[AUC_SN_full1$mesocosmID==USI_SN[i], ] #temp is a temporary data frame
  area_value <- pk.calc.auc(abs(temp$LRR_SN),temp$Sno)
  area_SN <- rbind(area_SN,data.frame(temp, area_value))
} #For SN
for(i in 1:length(USI_SP)){
  temp <- AUC_SP_full1[AUC_SP_full1$mesocosmID==USI_SP[i], ] #temp is a temporary data frame
  area_value <- pk.calc.auc(abs(temp$LRR_SP),temp$Sno)
  area_SP <- rbind(area_SP,data.frame(temp, area_value))
} #For SP

#Remove duplicated rows
area_C <- area_C %>% distinct(mesocosmID, .keep_all = TRUE) %>% mutate(Stab_Para = "C") %>% dplyr::select(Site, Season, Mesocosm, Treatment, Replicate, Stab_Para, AUC = area_value)
area_CN <- area_CN %>% distinct(mesocosmID, .keep_all = TRUE) %>% mutate(Stab_Para = "CN") %>% dplyr::select(Site, Season, Mesocosm, Treatment, Replicate, Stab_Para, AUC = area_value)
area_CP <- area_CP %>% distinct(mesocosmID, .keep_all = TRUE) %>% mutate(Stab_Para = "CP") %>% dplyr::select(Site, Season, Mesocosm, Treatment, Replicate, Stab_Para, AUC = area_value)
area_NP <- area_NP %>% distinct(mesocosmID, .keep_all = TRUE) %>% mutate(Stab_Para = "NP") %>% dplyr::select(Site, Season, Mesocosm, Treatment, Replicate, Stab_Para, AUC = area_value)
area_CS <- area_CS %>% distinct(mesocosmID, .keep_all = TRUE) %>% mutate(Stab_Para = "CS") %>% dplyr::select(Site, Season, Mesocosm, Treatment, Replicate, Stab_Para, AUC = area_value)
area_SN <- area_SN %>% distinct(mesocosmID, .keep_all = TRUE) %>% mutate(Stab_Para = "SN") %>% dplyr::select(Site, Season, Mesocosm, Treatment, Replicate, Stab_Para, AUC = area_value)
area_SP <- area_SP %>% distinct(mesocosmID, .keep_all = TRUE) %>% mutate(Stab_Para = "SP") %>% dplyr::select(Site, Season, Mesocosm, Treatment, Replicate, Stab_Para, AUC = area_value)

#Put it into form to merge with final stability table
AUC_full <- rbind(area_C, area_CN, area_CP, area_NP, area_CS, area_SN, area_SP)
AUC_full <- filter(AUC_full, Treatment != "C")
final_stability <- left_join(AUC_full, data_recovery, by = c("Site", "Season", "Mesocosm", "Treatment", "Replicate", "Stab_Para"))

#write_xlsx(final_stability, "SITES2_SmallPhyto_Stability_Nov2024.xlsx")

#PREPARE FOR PLOTTING

final_stability$Stab_Para <- as.factor(final_stability$Stab_Para)
final_stability$Rep <- as.numeric(final_stability$Replicate)
final_stability_mean1 <- group_by(final_stability, Site, Season, Treatment, Stab_Para)
final_stability_mean1 <- mutate(final_stability_mean1, Mean_AUC = mean(AUC), Sd_AUC = sd(AUC), Mean_Recovery = mean(Recovery), Sd_Recovery = sd(Recovery))
final_stability_mean1 <- dplyr::select(final_stability_mean1, -Mesocosm, -Replicate, -AUC, -Recovery, - Rep)
final_stability_mean1 <- unique(final_stability_mean1)
final_stability_mean1 <- pivot_longer(final_stability_mean1, cols = starts_with(c("Mean_", "Sd_")), names_to = c(".value", "variable"), names_pattern = "^(Mean|Sd)_([A-Za-z_]+)")

final_stability_mean1$Treatment <- factor(final_stability_mean1$Treatment, levels = c("D", "I", "E"))
#final_stability_mean1_basic <- filter(final_stability_mean1, Stab_Para == "C" | Stab_Para == "CN" | Stab_Para == "CP" | Stab_Para == "NP")
final_stability_mean1_basic <- final_stability_mean1
final_stability_mean1_basic$Stab_Para <- factor(final_stability_mean1_basic$Stab_Para, levels = c("C", "CN", "CP", "CS", "NP", "SN", "SP"))
final_stability_mean1_basic <- mutate(final_stability_mean1_basic, Unique_ID = paste(Site, Season, sep = "_"))
Fig3 <- ggplot(final_stability_mean1_basic, aes(x = Mean, y = Unique_ID, color = Treatment)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_point(position = position_dodge(width = -0.5), size = 4) +
  geom_errorbar(aes(xmin = Mean - Sd, xmax = Mean + Sd), width = 0.3, position = position_dodge(width = -0.5))+
  scale_colour_manual(values = c("palegreen4", "skyblue3", "indianred"))+
  xlab("") +
  ylab("") +
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+  
  theme(axis.title.x=element_blank())+
  theme(plot.title=element_text(hjust = 0.5))+
  facet_grid(Stab_Para~variable, scales = "free")
Fig3

Fig3_new <- ggplot(final_stability_mean1_basic, aes(x = Mean, y = Treatment, shape = Unique_ID, color = Treatment)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_point(position = position_dodge(width = -0.5), size = 3) +
  geom_errorbar(aes(xmin = Mean - Sd, xmax = Mean + Sd), width = 0.3, position = position_dodge(width = -0.5))+
  scale_colour_manual(values = c("palegreen4", "skyblue3", "indianred"))+
  scale_shape_manual(values=c(2, 19, 17))+
  xlab("") +
  ylab("") +
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+  
  theme(axis.title.x=element_blank())+
  theme(plot.title=element_text(hjust = 0.5))+
  facet_grid(Stab_Para~variable, scales = "free")
Fig3_new

#ggsave("Fig3_new.eps", Fig3_new , unit = "cm", device = "eps", width = 15, height = 20, dpi = 300)
#ggsave("Stability_SITES2.png", Fig3 , unit = "cm", device = "png", width = 15, height = 20, dpi = 300)

#_______________________________________________________________________________

#Plots for silicate ratios only

final_stability_mean1$Treatment <- factor(final_stability_mean1$Treatment, levels = c("D", "I", "E"))
final_stability_mean1_Si <- filter(final_stability_mean1, Stab_Para == "CS" | Stab_Para == "SN" | Stab_Para == "SP" | Stab_Para == "C")
final_stability_mean1_Si$Stab_Para <- factor(final_stability_mean1_Si$Stab_Para, levels = c("CS", "SN", "SP"))
final_stability_mean1_Si <- mutate(final_stability_mean1_Si, Unique_ID = paste(Site, Season, sep = "_"))
Fig3_new <- ggplot(final_stability_mean1_Si, aes(x = Mean, y = Treatment, shape = Unique_ID, color = Treatment)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_point(position = position_dodge(width = -0.5), size = 3) +
  geom_errorbar(aes(xmin = Mean - Sd, xmax = Mean + Sd), width = 0.3, position = position_dodge(width = -0.5))+
  scale_colour_manual(values = c("palegreen4", "skyblue3", "indianred"))+
  scale_shape_manual(values=c(2, 19, 17))+
  xlab("") +
  ylab("") +
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+  
  theme(axis.title.x=element_blank())+
  theme(plot.title=element_text(hjust = 0.5))+
  facet_grid(Stab_Para~variable, scales = "free")
Fig3_new

#_______________________________________________________________________________

#Histograms for delta AUC and delta recovery

library(writexl)
data_histo <- read_xlsx("stability_mean_SITES2_deltas.xlsx") #It was quicker in Excel
data_histo <- data_histo %>% mutate(StabID = paste(Stab_Para, variable, sep = "_"))

data_histo_AUC <- filter(data_histo, variable == "AUC")
data_histo_Recovery <- filter(data_histo, variable == "Recovery")

Histo_Deltas_AUC <- ggplot(data_histo_AUC, aes(x = abs(Value_Delta), y = Delta)) +
  geom_bar(stat = "identity") +
  facet_wrap(Treatment~StabID)+
  xlab("") +
  ylab("") +
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+  
  theme(axis.title.x=element_blank())+
  theme(plot.title=element_text(hjust = 0.5))
Histo_Deltas_AUC

Histo_Deltas_Recovery <- ggplot(data_histo_Recovery, aes(x = abs(Value_Delta), y = Delta)) +
  geom_bar(stat = "identity") +
  facet_wrap(Treatment~StabID)+
  xlab("") +
  ylab("") +
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+  
  theme(axis.title.x=element_blank())+
  theme(plot.title=element_text(hjust = 0.5))
Histo_Deltas_Recovery

#ggsave("Histo_Deltas_Recovery.eps", Histo_Deltas_Recovery , unit = "cm", device = "eps", width = 15, height = 20, dpi = 300)

#_______________________________________________________________________________

#Stats: Small size fraction (AUC)

#https://www.appsilon.com/post/manova-in-r

library(stats)
library(mvnormtest)
library(effectsize)
library(MASS) #For LDA post-hoc test
library(car)
library(MVN) #For multi-variate normality test
library(biotools)
library(readxl)
library(dplyr)
library(tidyverse)

setwd("/Users/Anika/Library/Mobile Documents/com~apple~CloudDocs/PhD/AquaCosm/SITES_RunOff_Comparison")
data_stats <- read_excel("SITES2_SmallPhyto_Stability_Nov2024.xlsx")
data_manova_Cratios <- filter(data_stats, Stab_Para == "C" | Stab_Para == "CN" | Stab_Para == "CP" | Stab_Para == "CS")

#Prepare the data set in a wide format
AUC_stats <- dplyr::select(data_manova_Cratios, -Recovery, -Replicate)
AUC_stats <- AUC_stats %>% pivot_wider(names_from = Stab_Para, values_from = AUC)
AUC_stats$Treatment <- as.factor(AUC_stats$Treatment)
AUC_stats$Site <- as.factor(AUC_stats$Site)
AUC_stats$Season <- as.factor(AUC_stats$Season)

#First look at the differences across seasons
AUC_stats_season <- dplyr::filter(AUC_stats, Site == "Erken")
AUC_stats_site <- dplyr::filter(AUC_stats, Season == "Summer")
AUC_stats_season <- dplyr::select(AUC_stats_season, Site, Mesocosm, Season, Treatment, C, CN, CP, CS)
AUC_stats_site <- dplyr::select(AUC_stats_site, Site, Mesocosm, Season, Treatment, C, CN, CP, CS)

#AUC Season

#Shapiro Test: Test assumptions of a MANOVA (Multivariate normality)
mardia_test <- mvn(AUC_stats_season[, c("C", "CN", "CP", "CS")], mvnTest = "mardia")
print(mardia_test) #Mardia’s skewness and kurtosis for multivariate normality
#Multivariate Normality: The data do not show significant deviations from multivariate normality based on Mardia’s skewness and kurtosis tests.
#Univariate Normality: C and CS shows significant deviation from normality.

#QQPlot: Test assumptions of a MANOVA (Multivariate normality)
response_data <- AUC_stats_season[5:8]
maha_dist <- mahalanobis(response_data, colMeans(response_data), cov(response_data)) #Calculate Mahalanobis distances, accounts for the correlations between variables and the variance of each variable
qqplot(qchisq(ppoints(nrow(response_data)), df = ncol(response_data)), maha_dist,
       main = "Chi-Square Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
abline(0, 1, col = "red")

#Homogeneity of Variance-Covariance Matrices
leveneTest(C ~ Treatment * Season, data = AUC_stats_season)
leveneTest(CN ~ Treatment * Season, data = AUC_stats_season)
leveneTest(CP ~ Treatment * Season, data = AUC_stats_season)
leveneTest(CS ~ Treatment * Season, data = AUC_stats_season)
#All values are above 0.05.

#Multicollinearity
cor_matrix <- cor(AUC_stats_season[, c("C", "CN", "CP", "CS")])
print(cor_matrix)
high_correlations <- which(abs(cor_matrix) > 0.9 & abs(cor_matrix) < 1, arr.ind = TRUE)
print(high_correlations)
#There are no high correlations (above 0.9) among the dependent variables.

#Conduct a MANOVA with Pillai’s Trace (most robust to departures from assumptions)
res.man <- manova(cbind(C, CN, CP, CS) ~ Treatment * Season, data = AUC_stats_season)
summary(res.man)
summary.aov(res.man)
eta_squared(res.man)

#Normality of residuals (Mardia's test for multivariate normality)
residuals_manova <- residuals(res.man)
mvn(residuals_manova, mvnTest = "mardia")
#Multivariate Normality: The residuals from the MANOVA model meet the assumption 
#of multivariate normality based on Mardia's test (both skewness and kurtosis).

#AUC Site

#AUC_stats_site <- mutate(AUC_stats_site, log_CP = log(CP))
#AUC_stats_site <- mutate(AUC_stats_site, log_NP = log(NP))

#Test assumptions of a MANOVA (Multivariate normality)
mardia_test <- mvn(AUC_stats_site[, c("C", "CN", "CP", "CS")], mvnTest = "mardia")
print(mardia_test) #Mardia’s skewness and kurtosis for multivariate normality
#The multivariate normality is accepted, although univariate normality is not given for CP and CSi.

#QQPlot: Test assumptions of a MANOVA (Multivariate normality)
response_data <- AUC_stats_site[5:8]
maha_dist <- mahalanobis(response_data, colMeans(response_data), cov(response_data)) #Calculate Mahalanobis distances, accounts for the correlations between variables and the variance of each variable
qqplot(qchisq(ppoints(nrow(response_data)), df = ncol(response_data)), maha_dist,
       main = "Chi-Square Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
abline(0, 1, col = "red")
#The QQ plot looks fine

#Homogeneity of Variance-Covariance Matrices
leveneTest(C ~ Treatment * Site, data = AUC_stats_site)
leveneTest(CN ~ Treatment * Site, data = AUC_stats_site)
leveneTest(CP ~ Treatment * Site, data = AUC_stats_site)
leveneTest(CS ~ Treatment * Site, data = AUC_stats_site)
#All good, all above 0.05.

#Multicollinearity
cor_matrix <- cor(AUC_stats_site[, c("C", "CN", "CP", "CS")])
print(cor_matrix)
high_correlations <- which(abs(cor_matrix) > 0.9 & abs(cor_matrix) < 1, arr.ind = TRUE)
print(high_correlations)
#There are no high correlations (above 0.9) among the dependent variables.

#Conduct a MANOVA with Pillai’s Trace (most robust to departures from assumptions)
res.man_site <- manova(cbind(C, CN, CP, CS) ~ Treatment * Site, data = AUC_stats_site)
summary(res.man_site)
summary.aov(res.man_site)
eta_squared(res.man_site)

#Normality of residuals (Mardia's test for multivariate normality)
residuals_manova <- residuals(res.man_site)
mvn(residuals_manova, mvnTest = "mardia")
#Multivariate Normality: The residuals from the MANOVA model meet the assumption 
#of multivariate normality based on Mardia's test (both skewness and kurtosis).

#Recovery Season

#data_stats <- read_excel("stability_for_stats_SITES2.xlsx")

#Prepare the data set in a wide format
Rec_stats <- dplyr::select(data_stats, -AUC, -Replicate)
Rec_stats <- Rec_stats %>% pivot_wider(names_from = Stab_Para, values_from = Recovery)
Rec_stats$Treatment <- as.factor(Rec_stats$Treatment)
Rec_stats$Site <- as.factor(Rec_stats$Site)
Rec_stats$Season <- as.factor(Rec_stats$Season)

#First look at the differences across seasons
Rec_stats_season <- dplyr::filter(Rec_stats, Site == "Erken")
Rec_stats_site <- dplyr::filter(Rec_stats, Season == "Summer")
Rec_stats_season <- dplyr::select(Rec_stats_season, Site, Mesocosm, Season, Treatment, C, CN, CP, CS)
Rec_stats_site <- dplyr::select(Rec_stats_site, Site, Mesocosm, Season, Treatment, C, CN, CP, CS)

#Rec_stats_season <- mutate(Rec_stats_season, log_C = (C^2-1)/2) #Leads to violation of normality of residuals

#Test assumptions of a MANOVA (Multivariate normality)
mardia_test <- mvn(Rec_stats_season[, c("C", "CN", "CP", "CS")], mvnTest = "mardia")
print(mardia_test) #Mardia’s skewness and kurtosis for multivariate normality
#All good in multivariate normality

#QQPlot: Test assumptions of a MANOVA (Multivariate normality)
response_data <- Rec_stats_season[5:8]
maha_dist <- mahalanobis(response_data, colMeans(response_data), cov(response_data)) #Calculate Mahalanobis distances, accounts for the correlations between variables and the variance of each variable
qqplot(qchisq(ppoints(nrow(response_data)), df = ncol(response_data)), maha_dist,
       main = "Chi-Square Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
abline(0, 1, col = "red")
#The QQ plot looks fine

#Homogeneity of Variance-Covariance Matrices
leveneTest(C ~ Treatment * Season, data = Rec_stats_season)
leveneTest(CN ~ Treatment * Season, data = Rec_stats_season)
leveneTest(CP ~ Treatment * Season, data = Rec_stats_season)
leveneTest(CS ~ Treatment * Season, data = Rec_stats_season)
#All good, all above 0.05.

#Multicollinearity
cor_matrix <- cor(Rec_stats_season[, c("C", "CN", "CP", "CS")])
print(cor_matrix)
high_correlations <- which(abs(cor_matrix) > 0.9 & abs(cor_matrix) < 1, arr.ind = TRUE)
print(high_correlations)

#Conduct a MANOVA with Pillai’s Trace (most robust to departures from assumptions)
res.man_seasonrec <- manova(cbind(C, CN, CP, CS) ~ Treatment * Season, data = Rec_stats_season)
#summary(res.man_seasonrec)
summary.aov(res.man_seasonrec)
eta_squared(res.man_seasonrec)

#Normality of residuals (Mardia's test for multivariate normality)
residuals_manova <- residuals(res.man_seasonrec)
mvn(residuals_manova, mvnTest = "mardia")

#Recovery Site

#Rec_stats_site <- mutate(Rec_stats_site, log_C = 1/C) #Leads to violation of normality of residuals
#Rec_stats_site <- mutate(Rec_stats_site, log_C = (C^2-1)/2) #Leads to violation of normality of residuals

#Test assumptions of a MANOVA (Multivariate normality)
mardia_test <- mvn(Rec_stats_site[, c("log_C", "CN", "CP", "CS")], mvnTest = "mardia")
print(mardia_test) #Mardia’s skewness and kurtosis for multivariate normality
#All good!

#QQPlot: Test assumptions of a MANOVA (Multivariate normality)
response_data <- Rec_stats_site[5:8]
maha_dist <- mahalanobis(response_data, colMeans(response_data), cov(response_data)) #Calculate Mahalanobis distances, accounts for the correlations between variables and the variance of each variable
qqplot(qchisq(ppoints(nrow(response_data)), df = ncol(response_data)), maha_dist,
       main = "Chi-Square Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
abline(0, 1, col = "red")
#The QQ plot looks fine

#Homogeneity of Variance-Covariance Matrices
leveneTest(C ~ Treatment * Site, data = Rec_stats_site)
leveneTest(CN ~ Treatment * Site, data = Rec_stats_site)
leveneTest(CP ~ Treatment * Site, data = Rec_stats_site)
leveneTest(CS ~ Treatment * Site, data = Rec_stats_site)
#All good, none above 0.05.

#Multicollinearity
cor_matrix <- cor(Rec_stats_site[, c("C", "CN", "CP", "CS")])
print(cor_matrix)
high_correlations <- which(abs(cor_matrix) > 0.9 & abs(cor_matrix) < 1, arr.ind = TRUE)
print(high_correlations)

#Conduct a MANOVA with Pillai’s Trace (most robust to departures from assumptions)
res.man_siterec <- manova(cbind(C, CN, CP, CS) ~ Treatment * Site, data = Rec_stats_site)
#summary(res.man_siterec)
summary.aov(res.man_siterec)
eta_squared(res.man_siterec)

#Normality of residuals (Mardia's test for multivariate normality)
residuals_manova <- residuals(res.man_siterec)
mvn(residuals_manova, mvnTest = "mardia")
#This multivariate normality is all good.

#QQPlot: Test assumptions of a MANOVA (Multivariate normality)
lapply(1:ncol(residuals_manova), function(i) {
  ggplot(data.frame(residuals = residuals_manova[, i]), aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line() +
    theme_minimal() +
    labs(title = paste("Q-Q Plot for Residuals of Variable", colnames(residuals_manova)[i]), x = "Theoretical Quantiles", y = "Sample Quantiles")
})

#_______________________________________________________________________________

#Subsequent ANOVAs and Tukey HSD tests for significant site/season effects

library(multcomp)
library(emmeans)

#AUC C for season

AUC_stats_season
model <- aov(C ~ Treatment * Season, data = AUC_stats_season)
summary(model)
shapiro.test(residuals(model)) #Normality of residuals
qqnorm(residuals(model)) #Q-Q plot for visual inspection of normality
qqline(residuals(model)) #Q-Q plot for visual inspection of normality
leveneTest(C ~ Treatment * Season, data = AUC_stats_season) #Levene's Test for homogeneity of variances
tukey <- TukeyHSD(model) #Conduct Tukey HSD test
print(tukey)

# Apply Bonferroni correction for Treatment:Season interaction
emmeans_interaction <- emmeans(model, ~ Treatment * Season)
bonferroni_interaction <- pairs(emmeans_interaction, adjust = "bonferroni")
print(bonferroni_interaction)

posthoc_C_season_AUC <- as.data.frame(bonferroni_interaction)

# Export the results to an Excel file
write_xlsx(posthoc_C_season_AUC, "posthoc_C_season_AUC.xlsx")

#_______________________________________________________________________________

#AUC C for site

AUC_stats_site
model <- aov(log(C) ~ Treatment * Site, data = AUC_stats_site)
summary(model)
shapiro.test(residuals(model)) #Normality of residuals
qqnorm(residuals(model)) #Q-Q plot for visual inspection of normality
qqline(residuals(model)) #Q-Q plot for visual inspection of normality
leveneTest(log(C) ~ Treatment * Season, data = AUC_stats_site) #Levene's Test for homogeneity of variances
tukey <- TukeyHSD(model) #Conduct Tukey HSD test
print(tukey)

# Apply Bonferroni correction for Treatment:Season interaction
emmeans_interaction <- emmeans(model, ~ Treatment * Site)
bonferroni_interaction <- pairs(emmeans_interaction, adjust = "bonferroni")
print(bonferroni_interaction)

posthoc_C_site_AUC <- as.data.frame(bonferroni_interaction)

# Export the results to an Excel file
write_xlsx(posthoc_C_site_AUC, "posthoc_C_site_AUC.xlsx")

#_______________________________________________________________________________

#AUC CN for season

AUC_stats_season
model <- aov(CN ~ Treatment * Season, data = AUC_stats_season)
summary(model)
shapiro.test(residuals(model)) #Normality of residuals
qqnorm(residuals(model)) #Q-Q plot for visual inspection of normality
qqline(residuals(model)) #Q-Q plot for visual inspection of normality
leveneTest(CN ~ Treatment * Season, data = AUC_stats_season) #Levene's Test for homogeneity of variances
tukey <- TukeyHSD(model) #Conduct Tukey HSD test
print(tukey)

# Apply Bonferroni correction for Treatment:Season interaction
emmeans_interaction <- emmeans(model, ~ Treatment * Season)
bonferroni_interaction <- pairs(emmeans_interaction, adjust = "bonferroni")
print(bonferroni_interaction)

posthoc_CN_season_AUC <- as.data.frame(bonferroni_interaction)

# Export the results to an Excel file
write_xlsx(posthoc_CN_season_AUC, "posthoc_CN_season_AUC.xlsx")

#_______________________________________________________________________________

#AUC CN for site

AUC_stats_site
model <- aov(CN ~ Treatment * Site, data = AUC_stats_site)
summary(model)
shapiro.test(residuals(model)) #Normality of residuals
qqnorm(residuals(model)) #Q-Q plot for visual inspection of normality
qqline(residuals(model)) #Q-Q plot for visual inspection of normality
leveneTest(CN ~ Treatment * Season, data = AUC_stats_site) #Levene's Test for homogeneity of variances
tukey <- TukeyHSD(model) #Conduct Tukey HSD test
print(tukey)

# Apply Bonferroni correction for Treatment:Season interaction
emmeans_interaction <- emmeans(model, ~ Treatment * Site)
bonferroni_interaction <- pairs(emmeans_interaction, adjust = "bonferroni")
print(bonferroni_interaction)

posthoc_CN_site_AUC <- as.data.frame(bonferroni_interaction)

# Export the results to an Excel file
write_xlsx(posthoc_CN_site_AUC, "posthoc_CN_site_AUC.xlsx")

#_______________________________________________________________________________

#AUC CP for season

AUC_stats_season
model <- aov(CP ~ Treatment * Season, data = AUC_stats_season)
summary(model)
shapiro.test(residuals(model)) #Normality of residuals
qqnorm(residuals(model)) #Q-Q plot for visual inspection of normality
qqline(residuals(model)) #Q-Q plot for visual inspection of normality
leveneTest(CP ~ Treatment * Season, data = AUC_stats_season) #Levene's Test for homogeneity of variances
tukey <- TukeyHSD(model) #Conduct Tukey HSD test
print(tukey)

# Apply Bonferroni correction for Treatment:Season interaction
emmeans_interaction <- emmeans(model, ~ Treatment * Season)
bonferroni_interaction <- pairs(emmeans_interaction, adjust = "bonferroni")
print(bonferroni_interaction)

posthoc_CP_season_AUC <- as.data.frame(bonferroni_interaction)

# Export the results to an Excel file
write_xlsx(posthoc_CP_season_AUC, "posthoc_CP_season_AUC.xlsx")

#_______________________________________________________________________________

#AUC CN for site

AUC_stats_site
model <- aov(log(CP) ~ Treatment * Site, data = AUC_stats_site)
summary(model)
shapiro.test(residuals(model)) #Normality of residuals
qqnorm(residuals(model)) #Q-Q plot for visual inspection of normality
qqline(residuals(model)) #Q-Q plot for visual inspection of normality
leveneTest(log(CP) ~ Treatment * Season, data = AUC_stats_site) #Levene's Test for homogeneity of variances
tukey <- TukeyHSD(model) #Conduct Tukey HSD test
print(tukey)

# Apply Bonferroni correction for Treatment:Season interaction
emmeans_interaction <- emmeans(model, ~ Treatment * Site)
bonferroni_interaction <- pairs(emmeans_interaction, adjust = "bonferroni")
print(bonferroni_interaction)

posthoc_CP_site_AUC <- as.data.frame(bonferroni_interaction)

# Export the results to an Excel file
write_xlsx(posthoc_CP_site_AUC, "posthoc_CP_site_AUC.xlsx")

#_______________________________________________________________________________

#AUC CS for season

AUC_stats_season
model <- aov(CS ~ Treatment * Season, data = AUC_stats_season)
summary(model)
shapiro.test(residuals(model)) #Normality of residuals
qqnorm(residuals(model)) #Q-Q plot for visual inspection of normality
qqline(residuals(model)) #Q-Q plot for visual inspection of normality
leveneTest(CS ~ Treatment * Season, data = AUC_stats_season) #Levene's Test for homogeneity of variances
tukey <- TukeyHSD(model) #Conduct Tukey HSD test
print(tukey)

# Apply Bonferroni correction for Treatment:Season interaction
emmeans_interaction <- emmeans(model, ~ Treatment * Season)
bonferroni_interaction <- pairs(emmeans_interaction, adjust = "bonferroni")
print(bonferroni_interaction)

posthoc_CS_season_AUC <- as.data.frame(bonferroni_interaction)

# Export the results to an Excel file
write_xlsx(posthoc_CS_season_AUC, "posthoc_CS_season_AUC.xlsx")

#_______________________________________________________________________________

#AUC CS for site

AUC_stats_site
model <- aov(CS ~ Treatment * Site, data = AUC_stats_site)
summary(model)
shapiro.test(residuals(model)) #Normality of residuals
qqnorm(residuals(model)) #Q-Q plot for visual inspection of normality
qqline(residuals(model)) #Q-Q plot for visual inspection of normality
leveneTest(CS ~ Treatment * Season, data = AUC_stats_site) #Levene's Test for homogeneity of variances
tukey <- TukeyHSD(model) #Conduct Tukey HSD test
print(tukey)

# Apply Bonferroni correction for Treatment:Season interaction
emmeans_interaction <- emmeans(model, ~ Treatment * Site)
bonferroni_interaction <- pairs(emmeans_interaction, adjust = "bonferroni")
print(bonferroni_interaction)

posthoc_CS_site_AUC <- as.data.frame(bonferroni_interaction)

# Export the results to an Excel file
write_xlsx(posthoc_CS_site_AUC, "posthoc_CS_site_AUC.xlsx")

#_______________________________________________________________________________

#C Recovery Season

Rec_stats_season

model <- aov(log(C) ~ Treatment * Season, data = Rec_stats_season)
summary(model)
shapiro.test(residuals(model)) #Normality of residuals
qqnorm(residuals(model)) #Q-Q plot for visual inspection of normality
qqline(residuals(model)) #Q-Q plot for visual inspection of normality
leveneTest(C ~ Treatment * Season, data = AUC_stats_season) #Levene's Test for homogeneity of variances
tukey <- TukeyHSD(model) #Conduct Tukey HSD test
print(tukey)

# Apply Bonferroni correction for Treatment:Season interaction
emmeans_interaction <- emmeans(model, ~ Treatment * Season)
bonferroni_interaction <- pairs(emmeans_interaction, adjust = "bonferroni")
print(bonferroni_interaction)

posthoc_C_season_Rec <- as.data.frame(bonferroni_interaction)

# Export the results to an Excel file
write_xlsx(posthoc_C_season_Rec, "posthoc_C_season_Rec.xlsx")

#_______________________________________________________________________________

#CN Recovery Season

Rec_stats_season

model <- aov(CN ~ Treatment * Season, data = Rec_stats_season)
summary(model)
shapiro.test(residuals(model)) #Normality of residuals
qqnorm(residuals(model)) #Q-Q plot for visual inspection of normality
qqline(residuals(model)) #Q-Q plot for visual inspection of normality
leveneTest(CN ~ Treatment * Season, data = AUC_stats_season) #Levene's Test for homogeneity of variances
tukey <- TukeyHSD(model) #Conduct Tukey HSD test
print(tukey)

# Apply Bonferroni correction for Treatment:Season interaction
emmeans_interaction <- emmeans(model, ~ Treatment * Season)
bonferroni_interaction <- pairs(emmeans_interaction, adjust = "bonferroni")
print(bonferroni_interaction)

posthoc_CN_season_Rec <- as.data.frame(bonferroni_interaction)

# Export the results to an Excel file
write_xlsx(posthoc_CN_season_Rec, "posthoc_CN_season_Rec.xlsx")

#_______________________________________________________________________________

#CP Recovery Season

Rec_stats_season

model <- aov(CP ~ Treatment * Season, data = Rec_stats_season)
summary(model)
shapiro.test(residuals(model)) #Normality of residuals
qqnorm(residuals(model)) #Q-Q plot for visual inspection of normality
qqline(residuals(model)) #Q-Q plot for visual inspection of normality
leveneTest(CP ~ Treatment * Season, data = AUC_stats_season) #Levene's Test for homogeneity of variances
tukey <- TukeyHSD(model) #Conduct Tukey HSD test
print(tukey)

# Apply Bonferroni correction for Treatment:Season interaction
emmeans_interaction <- emmeans(model, ~ Treatment * Season)
bonferroni_interaction <- pairs(emmeans_interaction, adjust = "bonferroni")
print(bonferroni_interaction)

posthoc_CP_season_Rec <- as.data.frame(bonferroni_interaction)

# Export the results to an Excel file
write_xlsx(posthoc_CP_season_Rec, "posthoc_CP_season_Rec.xlsx")

#_______________________________________________________________________________

#CP Recovery Season

Rec_stats_season

model <- aov(CS ~ Treatment * Season, data = Rec_stats_season)
summary(model)
shapiro.test(residuals(model)) #Normality of residuals
qqnorm(residuals(model)) #Q-Q plot for visual inspection of normality
qqline(residuals(model)) #Q-Q plot for visual inspection of normality
leveneTest(CS ~ Treatment * Season, data = AUC_stats_season) #Levene's Test for homogeneity of variances
tukey <- TukeyHSD(model) #Conduct Tukey HSD test
print(tukey)

# Apply Bonferroni correction for Treatment:Season interaction
emmeans_interaction <- emmeans(model, ~ Treatment * Season)
bonferroni_interaction <- pairs(emmeans_interaction, adjust = "bonferroni")
print(bonferroni_interaction)

posthoc_CS_season_Rec <- as.data.frame(bonferroni_interaction)

# Export the results to an Excel file
write_xlsx(posthoc_CS_season_Rec, "posthoc_CS_season_Rec.xlsx")

#_______________________________________________________________________________

#CP Recovery Site

model <- aov(CP ~ Treatment * Site, data = Rec_stats_site)
summary(model)
shapiro.test(residuals(model)) #Normality of residuals
qqnorm(residuals(model)) #Q-Q plot for visual inspection of normality
qqline(residuals(model)) #Q-Q plot for visual inspection of normality
leveneTest(CP ~ Treatment * Site, data = AUC_stats_site) #Levene's Test for homogeneity of variances
tukey <- TukeyHSD(model) #Conduct Tukey HSD test
print(tukey)

# Apply Bonferroni correction for Treatment:Season interaction
emmeans_interaction <- emmeans(model, ~ Treatment * Site)
bonferroni_interaction <- pairs(emmeans_interaction, adjust = "bonferroni")
print(bonferroni_interaction)

posthoc_CP_site_Rec <- as.data.frame(bonferroni_interaction)

# Export the results to an Excel file
write_xlsx(posthoc_CP_site_Rec, "posthoc_CP_site_Rec.xlsx")

#_______________________________________________________________________________

#CS Recovery Site

model <- aov(CS ~ Treatment * Site, data = Rec_stats_site)
summary(model)
shapiro.test(residuals(model)) #Normality of residuals
qqnorm(residuals(model)) #Q-Q plot for visual inspection of normality
qqline(residuals(model)) #Q-Q plot for visual inspection of normality
leveneTest(CS ~ Treatment * Site, data = AUC_stats_site) #Levene's Test for homogeneity of variances
tukey <- TukeyHSD(model) #Conduct Tukey HSD test
print(tukey)

# Apply Bonferroni correction for Treatment:Season interaction
emmeans_interaction <- emmeans(model, ~ Treatment * Site)
bonferroni_interaction <- pairs(emmeans_interaction, adjust = "bonferroni")
print(bonferroni_interaction)

posthoc_CS_site_Rec <- as.data.frame(bonferroni_interaction)

# Export the results to an Excel file
write_xlsx(posthoc_CS_site_Rec, "posthoc_CS_site_Rec.xlsx")

#_______________________________________________________________________________







