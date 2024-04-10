# In this script, I do simulations with my heavy datasets, utilising the HPC
# resources :D

# Call libraries ----------------------------------------------------------

library(styler)
library(dplyr)
library(tidyr)
library(ggplot2)
library(conflicted)
library(mice) #for 2l.pan imputation
library(pan) #for 2l.pan imputation
library(miceadds) #for 2l.pmm imputation
library(viridis) #for using viridis color palette
library(grid) #to put annotation
library(gridExtra) #to make plot grids
library(cowplot) #to make plot grids
library(pROC) #for making roc curves
library(readr) #part of tidyr #data wrangling
library(purrr) #to sample from a ecdf
library(xpose) #for making vpcs
library(xpose4) #for making vpcs #having issues with lattice package
library(vpc) #for making vpcs
library(PsNR) #for making vpcs
library(ggtext) #create latex axis labels - does not work
library(latex2exp) #to create latex label
library(table1) # create descriptive statistic table1
library(ggpubr) # for ggarrange function 
library(extrafont)


# Analyzing simulation datasets -------------------------------------------

## Firstly, CKDEPI effect -------------------------------------------

### Creating datasets for PTA plots -------------------------------------------

#### Standard dosing PTA -------------------------------------------

# Load the dataset
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/CKDEPI_Sim/'
simulations_tablefile <- paste0(working.directory, '/sim_CKDEPI_01.1.npctab.dta') 
dataframe_simulations <- read_nonmem_table(simulations_tablefile)

# Extract only the trough events
CKDEPI_dos_01 <- dataframe_simulations[dataframe_simulations$TROUGH == 1,]

# Convert ID into numeric
CKDEPI_dos_01$ID <- as.integer(CKDEPI_dos_01$ID)

# Create a new variable AUC24 for each ID
CKDEPI_dos_01 <- CKDEPI_dos_01 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

CKDEPI_dos_01 <- CKDEPI_dos_01 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
CKDEPI_dos_01$fAUC <- CKDEPI_dos_01$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and CKDEPI columns
PTA_01 <- subset(CKDEPI_dos_01, select = c("ID", "DV", "fAUC", "DAY", "CKDEPI"))

# Calculate percentage of DV >= 7.5 per DAY per CKDEPI (PTA_Cmin_75)
PTA_01$PTA_Cmin_75 <- 100 * ave(PTA_01$DV >= 7.5, PTA_01$DAY, PTA_01$CKDEPI, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per CKDEPI (PTA_Cmin_80)
PTA_01$PTA_Cmin_80 <- 100 * ave(PTA_01$DV >= 80, PTA_01$DAY, PTA_01$CKDEPI, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per CKDEPI (PTA_fAUC_200)
PTA_01$PTA_fAUC_200 <- 100 * ave(PTA_01$fAUC >= 200, PTA_01$DAY, PTA_01$CKDEPI, FUN = mean)

# Create PTA_01_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and CKDEPI
PTA_01_overall <- unique(PTA_01[c("PTA_Cmin_75", "PTA_Cmin_80", "PTA_fAUC_200", "DAY", "CKDEPI")])

# Export PTA_dos_01 dataset for standard dosing
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/CKDEPI_Sim/")
write.csv(PTA_01_overall, "PTA_dos_01.csv", quote = F, row.names = FALSE)

#### 1000 LD 400 MD PTA -------------------------------------------

# Load the dataset
working.directory0 <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/CKDEPI_Sim/'
directory_CKDEPI_02 <- paste0(working.directory, '/sim_CKDEPI_02.1.npctab.dta') 
CKDEPI_02 <- read_nonmem_table(directory_CKDEPI_02)

# Extract only the trough events
CKDEPI_dos_02 <- CKDEPI_02[CKDEPI_02$TROUGH == 1,]

# Convert ID into numeric
CKDEPI_dos_02$ID <- as.integer(CKDEPI_dos_02$ID)

# Create a new variable AUC24 for each ID
CKDEPI_dos_02 <- CKDEPI_dos_02 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

CKDEPI_dos_02 <- CKDEPI_dos_02 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
CKDEPI_dos_02$fAUC <- CKDEPI_dos_02$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and CKDEPI columns
PTA_02 <- subset(CKDEPI_dos_02, select = c("ID", "DV", "fAUC", "DAY", "CKDEPI"))

# Calculate percentage of DV >= 7.5 per DAY per CKDEPI (PTA_Cmin_75)
PTA_02$PTA_Cmin_75 <- 100 * ave(PTA_02$DV >= 7.5, PTA_02$DAY, PTA_02$CKDEPI, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per CKDEPI (PTA_Cmin_80)
PTA_02$PTA_Cmin_80 <- 100 * ave(PTA_02$DV >= 80, PTA_02$DAY, PTA_02$CKDEPI, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per CKDEPI (PTA_fAUC_200)
PTA_02$PTA_fAUC_200 <- 100 * ave(PTA_02$fAUC >= 200, PTA_02$DAY, PTA_02$CKDEPI, FUN = mean)

# Create PTA_02_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and CKDEPI
PTA_02_overall <- unique(PTA_02[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200", "DAY", "CKDEPI")])

# Export PTA_dos_02 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/CKDEPI_Sim/")
write.csv(PTA_02_overall, "PTA_dos_02.csv",quote = F,row.names = FALSE)

#### 1200 LD 400 MD PTA -------------------------------------------

# Load the dataset
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/CKDEPI_Sim/'
directory_CKDEPI_03 <- paste0(working.directory, '/sim_CKDEPI_03.1.npctab.dta') 
CKDEPI_03 <- read_nonmem_table(directory_CKDEPI_03)

# Extract only the trough events
CKDEPI_dos_03 <- CKDEPI_03[CKDEPI_03$TROUGH == 1,]

# Convert ID into numeric
CKDEPI_dos_03$ID <- as.integer(CKDEPI_dos_03$ID)

# Create a new variable AUC24 for each ID
CKDEPI_dos_03 <- CKDEPI_dos_03 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

CKDEPI_dos_03 <- CKDEPI_dos_03 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2, AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
CKDEPI_dos_03$fAUC <- CKDEPI_dos_03$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and CKDEPI columns
PTA_03 <- subset(CKDEPI_dos_03, select = c("ID", "DV", "fAUC", "DAY", "CKDEPI"))

# Calculate percentage of DV >= 7.5 per DAY per CKDEPI (PTA_Cmin_75)
PTA_03$PTA_Cmin_75 <- 100 * ave(PTA_03$DV >= 7.5, PTA_03$DAY, PTA_03$CKDEPI, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per CKDEPI (PTA_Cmin_80)
PTA_03$PTA_Cmin_80 <- 100 * ave(PTA_03$DV >= 80, PTA_03$DAY, PTA_03$CKDEPI, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per CKDEPI (PTA_fAUC_200)
PTA_03$PTA_fAUC_200 <- 100 * ave(PTA_03$fAUC >= 200, PTA_03$DAY, PTA_03$CKDEPI, FUN = mean)

# Create PTA_03_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and CKDEPI
PTA_03_overall <- unique(PTA_03[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200", "DAY", "CKDEPI")])

# Export PTA_dos_03 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/CKDEPI_Sim/")
write.csv(PTA_03_overall, "PTA_dos_03.csv",quote = F,row.names = FALSE)

#### 12 mg/kg LD 6 mg/kg MD PTA -------------------------------------------

# Load the dataset
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/CKDEPI_Sim/'
directory_CKDEPI_04 <- paste0(working.directory, '/sim_CKDEPI_04.1.npctab.dta') 
CKDEPI_04 <- read_nonmem_table(directory_CKDEPI_04)

# Extract only the trough events
CKDEPI_dos_04 <- CKDEPI_04[CKDEPI_04$TROUGH == 1,]

# Convert ID into numeric
CKDEPI_dos_04$ID <- as.integer(CKDEPI_dos_04$ID)

# Create a new variable AUC24 for each ID
CKDEPI_dos_04 <- CKDEPI_dos_04 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

CKDEPI_dos_04 <- CKDEPI_dos_04 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2, AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
CKDEPI_dos_04$fAUC <- CKDEPI_dos_04$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and CKDEPI columns
PTA_04 <- subset(CKDEPI_dos_04, select = c("ID", "DV", "fAUC", "DAY", "CKDEPI"))

# Calculate percentage of DV >= 7.5 per DAY per CKDEPI (PTA_Cmin_75)
PTA_04$PTA_Cmin_75 <- 100 * ave(PTA_04$DV >= 7.5, PTA_04$DAY, PTA_04$CKDEPI, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per CKDEPI (PTA_Cmin_80)
PTA_04$PTA_Cmin_80 <- 100 * ave(PTA_04$DV >= 80, PTA_04$DAY, PTA_04$CKDEPI, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per CKDEPI (PTA_fAUC_200)
PTA_04$PTA_fAUC_200 <- 100 * ave(PTA_04$fAUC >= 200, PTA_04$DAY, PTA_04$CKDEPI, FUN = mean)

# Create PTA_04_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and CKDEPI
PTA_04_overall <- unique(PTA_04[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200", "DAY", "CKDEPI")])

# Export PTA_dos_04 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/CKDEPI_Sim/")
write.csv(PTA_04_overall, "PTA_dos_04.csv",quote = F, row.names = FALSE)

#### 1200 LD 600 MD PTA -------------------------------------------

# Load the dataset
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/CKDEPI_Sim/'
directory_CKDEPI_05 <- paste0(working.directory, '/sim_CKDEPI_05.1.npctab.dta') 
CKDEPI_05 <- read_nonmem_table(directory_CKDEPI_05)

# Extract only the trough events
CKDEPI_dos_05 <- CKDEPI_05[CKDEPI_05$TROUGH == 1,]

# Convert ID into numeric
CKDEPI_dos_05$ID <- as.integer(CKDEPI_dos_05$ID)

# Create a new variable AUC24 for each ID
CKDEPI_dos_05 <- CKDEPI_dos_05 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

CKDEPI_dos_05 <- CKDEPI_dos_05 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2, AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
CKDEPI_dos_05$fAUC <- CKDEPI_dos_05$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and CKDEPI columns
PTA_05 <- subset(CKDEPI_dos_05, select = c("ID", "DV", "fAUC", "DAY", "CKDEPI"))

# Calculate percentage of DV >= 7.5 per DAY per CKDEPI (PTA_Cmin_75)
PTA_05$PTA_Cmin_75 <- 100 * ave(PTA_05$DV >= 7.5, PTA_05$DAY, PTA_05$CKDEPI, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per CKDEPI (PTA_Cmin_80)
PTA_05$PTA_Cmin_80 <- 100 * ave(PTA_05$DV >= 80, PTA_05$DAY, PTA_05$CKDEPI, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per CKDEPI (PTA_fAUC_200)
PTA_05$PTA_fAUC_200 <- 100 * ave(PTA_05$fAUC >= 200, PTA_05$DAY, PTA_05$CKDEPI, FUN = mean)

# Create PTA_05_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and CKDEPI
PTA_05_overall <- unique(PTA_05[c("PTA_Cmin_75", "PTA_Cmin_80", "PTA_fAUC_200", "DAY", "CKDEPI")])

# Export PTA_dos_05 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/CKDEPI_Sim/")
write.csv(PTA_05_overall, "PTA_dos_05.csv", quote = F, row.names = FALSE)

#### 800 LD 200 MD PTA -------------------------------------------

# Load the dataset
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/CKDEPI_Sim/'
directory_CKDEPI_06 <- paste0(working.directory, '/sim_CKDEPI_06.1.npctab.dta') 
CKDEPI_06 <- read_nonmem_table(directory_CKDEPI_06)

# Extract only the trough events
CKDEPI_dos_06 <- CKDEPI_06[CKDEPI_06$TROUGH == 1,]

# Convert ID into numeric
CKDEPI_dos_06$ID <- as.integer(CKDEPI_dos_06$ID)

# Create a new variable AUC24 for each ID
CKDEPI_dos_06 <- CKDEPI_dos_06 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

CKDEPI_dos_06 <- CKDEPI_dos_06 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2, AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
CKDEPI_dos_06$fAUC <- CKDEPI_dos_06$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and CKDEPI columns
PTA_06 <- subset(CKDEPI_dos_06, select = c("ID", "DV", "fAUC", "DAY", "CKDEPI"))

# Calculate percentage of DV >= 7.5 per DAY per CKDEPI (PTA_Cmin_75)
PTA_06$PTA_Cmin_75 <- 100 * ave(PTA_06$DV >= 7.5, PTA_06$DAY, PTA_06$CKDEPI, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per CKDEPI (PTA_Cmin_80)
PTA_06$PTA_Cmin_80 <- 100 * ave(PTA_06$DV >= 80, PTA_06$DAY, PTA_06$CKDEPI, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per CKDEPI (PTA_fAUC_200)
PTA_06$PTA_fAUC_200 <- 100 * ave(PTA_06$fAUC >= 200, PTA_06$DAY, PTA_06$CKDEPI, FUN = mean)

# Create PTA_06_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and CKDEPI
PTA_06_overall <- unique(PTA_06[c("PTA_Cmin_75", "PTA_Cmin_80", "PTA_fAUC_200", "DAY", "CKDEPI")])

# Export PTA_dos_06 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/CKDEPI_Sim/")
write.csv(PTA_06_overall, "PTA_dos_06.csv", quote = F, row.names = FALSE)

### Making PTA plots -------------------------------------------

#### Some preps -------------------------------------------

# Import PTA_dos_ datasets
PTA_dos_01 <- read.csv("PTA_CKDEPI/PTA_dos_01.csv")
PTA_dos_02 <- read.csv("PTA_CKDEPI/PTA_dos_02.csv")
PTA_dos_03 <- read.csv("PTA_CKDEPI/PTA_dos_03.csv")
PTA_dos_04 <- read.csv("PTA_CKDEPI/PTA_dos_04.csv")
PTA_dos_05 <- read.csv("PTA_CKDEPI/PTA_dos_05.csv")
PTA_dos_06 <- read.csv("PTA_CKDEPI/PTA_dos_06.csv")

# Add regimens
PTA_dos_01$Regimen <- 1
PTA_dos_02$Regimen <- 2
PTA_dos_03$Regimen <- 3
PTA_dos_04$Regimen <- 4
PTA_dos_05$Regimen <- 5
PTA_dos_06$Regimen <- 6

# Combine 6 PTA datasets into 1 - PTA_CKDEPI
PTA_CKDEPI <- rbind(PTA_dos_01, PTA_dos_02, PTA_dos_03, 
                    PTA_dos_04, PTA_dos_05, PTA_dos_06)

# Define the labels and colors
dosing_labels_day12 <- c("800 mg q24h LD",
                        "1000 mg q24h LD",
                        "1200 mg q24h LD",
                        "12 mg/kg q24h LD")
# Regimens 6, 2, 5, 4, respectively

dosing_labels_day714 <- c("200 mg q24h MD", 
                         "400 mg q24h MD", 
                         "600 mg q24h MD",
                         "6 mg/kg q24h MD")
# Regimens 6, 2, 5, 4, respectively

dosing_colors_day12714 <- c("#5ec962",
                           "#21918c",
                           "#3b528b",
                           "#440154")

# Reorder the levels of the dosing regimens
PTA_CKDEPI$Regimen <- factor(PTA_CKDEPI$Regimen,
                             levels = c(6, 2, 5, 4, 3, 1))

#### Day 1 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_CKDEPI200_DAY1 <-  ggplot(PTA_CKDEPI[PTA_CKDEPI$DAY == 1 & PTA_CKDEPI$Regimen %in% c(6, 2, 5, 4), ], 
                              aes(x = CKDEPI, 
                                  y = PTA_fAUC_200, 
                                  group = factor(Regimen), 
                                  color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(5, 215), breaks = seq(5, 215, by = 10), expand = c(0,0)) +
  xlab(TeX(r"($\eGFR_{CKD-EPI}\ (ml/min/1.73m^2)$)")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) +
  ylab(TeX(r"($Probability\ of\ \textit{f}AUC_{0-24}\ \geq 200 \, \mg \times h/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 1") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day12714, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day12) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "gray50", size = 0.6) 
PTA_CKDEPI200_DAY1

#### Day 14 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_CKDEPI200_DAY14 <-  ggplot(PTA_CKDEPI[PTA_CKDEPI$DAY == 14 & PTA_CKDEPI$Regimen %in% c(6, 2, 5, 4), ], 
                              aes(x = CKDEPI, 
                                  y = PTA_fAUC_200, 
                                  group = factor(Regimen), 
                                  color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(5, 215), breaks = seq(5, 215, by = 10), expand = c(0,0)) +
  xlab(TeX(r"($\eGFR_{CKD-EPI}\ (ml/min/1.73m^2)$)")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) +
  ylab(TeX(r"($Probability\ of\ \textit{f}AUC_{0-24}\ \geq 200 \, \mg \times h/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 14") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day12714, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day714) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "gray50", size = 0.6) 
PTA_CKDEPI200_DAY14

#### Day 2 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_CKDEPI200_DAY2 <-  ggplot(PTA_CKDEPI[PTA_CKDEPI$DAY == 2 & PTA_CKDEPI$Regimen %in% c(6, 2, 5, 4), ], 
                              aes(x = CKDEPI, 
                                  y = PTA_fAUC_200, 
                                  group = factor(Regimen), 
                                  color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(5, 215), breaks = seq(5, 215, by = 10), expand = c(0,0)) +
  xlab(TeX(r"($\eGFR_{CKD-EPI}\ (ml/min/1.73m^2)$)")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) +
  ylab(TeX(r"($Probability\ of\ \textit{f}AUC_{0-24}\ \geq 200 \, \mg \times h/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 2") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day12714, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day12) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "gray50", size = 0.6) 
PTA_CKDEPI200_DAY2

#### Day 7 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_CKDEPI200_DAY7 <-  ggplot(PTA_CKDEPI[PTA_CKDEPI$DAY == 7 & PTA_CKDEPI$Regimen %in% c(6, 2, 5, 4), ], 
                               aes(x = CKDEPI, 
                                   y = PTA_fAUC_200, 
                                   group = factor(Regimen), 
                                   color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(5, 215), breaks = seq(5, 215, by = 10), expand = c(0,0)) +
  xlab(TeX(r"($\eGFR_{CKD-EPI}\ (ml/min/1.73m^2)$)")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) +
  ylab(TeX(r"($Probability\ of\ \textit{f}AUC_{0-24}\ \geq 200 \, \mg \times h/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 7") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day12714, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day714) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "gray50", size = 0.6) 
PTA_CKDEPI200_DAY7

#### Combining these 4 plots into 1 -------------------------------------------

## Combining DAY 1 & 2 into 1

# Define common legend for both plots
common_legend_CKDEPI_200_DAY12 <- ggpubr::get_legend(PTA_CKDEPI200_DAY1)

# Remove legends from individual plots
PTA_CKDEPI200_DAY1_no_legend <- PTA_CKDEPI200_DAY1 + 
  theme(legend.position = "none")

PTA_CKDEPI200_DAY2_no_legend <- PTA_CKDEPI200_DAY2 +
  theme(legend.position = "none")

# Combine the two plots
PTA_CKDEPI200_DAY12_no_legend <- ggarrange(PTA_CKDEPI200_DAY1_no_legend, 
                                           PTA_CKDEPI200_DAY2_no_legend,
                                           labels = c("a", "b"),
                                           ncol = 2, nrow = 1,
                                           font.label = list(size = 8, face = "bold"))

# Add common legend

grid_CKDEPI_200_DAY12 <- plot_grid(PTA_CKDEPI200_DAY12_no_legend, 
                                   common_legend_CKDEPI_200_DAY12,
                                   rel_widths = c(1, 0.13),
                                   nrow = 1)

grid_CKDEPI_200_DAY12

## Combining these DAY 7 & 14 into 1

# Define common legend for both plots
common_legend_CKDEPI_200_DAY714 <- ggpubr::get_legend(PTA_CKDEPI200_DAY7)

# Remove legends from individual plots
PTA_CKDEPI200_DAY7_no_legend <- PTA_CKDEPI200_DAY7 + 
  theme(legend.position = "none")

PTA_CKDEPI200_DAY14_no_legend <- PTA_CKDEPI200_DAY14 +
  theme(legend.position = "none")

# Combine the two plots
PTA_CKDEPI200_DAY714_no_legend <- ggarrange(PTA_CKDEPI200_DAY7_no_legend, 
                                            PTA_CKDEPI200_DAY14_no_legend,
                                           labels = c("c", "d"),
                                           ncol = 2, nrow = 1,
                                           font.label = list(size = 8, face = "bold"))

# Add common legend

grid_CKDEPI_200_DAY714 <- plot_grid(PTA_CKDEPI200_DAY714_no_legend, 
                                   common_legend_CKDEPI_200_DAY714,
                                   rel_widths = c(1, 0.13),
                                   nrow = 1)

grid_CKDEPI_200_DAY714

## Then combining these 4 plots
grid_CKDEPI_200 <- plot_grid(grid_CKDEPI_200_DAY12, 
                             grid_CKDEPI_200_DAY714, 
                             nrow = 2, 
                             align = "v")

grid_CKDEPI_200

## Combining DAY 1 & 14 into 1

# Combine the two plots
PTA_CKDEPI200_DAY114 <- ggarrange(PTA_CKDEPI200_DAY1,
                                  PTA_CKDEPI200_DAY14,
                                     labels = c("a", "b"),
                                     ncol = 1, 
                                     nrow = 2,
                                     font.label = list(size = 8, face = "bold"))

# Export the combined plots
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Plots/Dose_finding_simulations/CKDEPI")

# SVG
ggsave("grid_CKDEPI_200.svg", 
       grid_CKDEPI_200, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")
# This is the part when I realised that PTA of day 2 is affected by not just
# the loading dose, but also the maintenance dose. The same for PTA of day 7,
# which is not just influenced by the MD alone. So I decided to plot PTA of 
# day 1 and 14 only as below. Doing so will make sure that the effects of LD and MD are
# reflected to the fullest.

ggsave("PTA_CKDEPI200_DAY114.svg", 
       PTA_CKDEPI200_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

# JPEG
ggsave("PTA_CKDEPI200_DAY114.JPEG", 
       PTA_CKDEPI200_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

# EPS
ggsave("PTA_CKDEPI200_DAY114.EPS", 
       PTA_CKDEPI200_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

### Assessing toxicity PTA Cmin80 --------------------------------------------

# PTA_CKDEPI dataset was created above

#### Day 1 -------------------------------------------

# Check the title if it's correct
plot(TeX(r"($Probability\ of\ \C_{min}\ \geq 80 \, \mg/L\ (\%)$)"))

# Here we use Colorblind-Friendly Palette Viridis
PTA_CKDEPI80_DAY1 <-  ggplot(PTA_CKDEPI[PTA_CKDEPI$DAY == 1 & PTA_CKDEPI$Regimen %in% c(6, 2, 5, 4), ], 
                              aes(x = CKDEPI, 
                                  y = PTA_Cmin_80, 
                                  group = factor(Regimen), 
                                  color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(5, 215), breaks = seq(5, 215, by = 10), expand = c(0,0)) +
  xlab(TeX(r"($\eGFR_{CKD-EPI}\ (ml/min/1.73m^2)$)")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) +
  ylab(TeX(r"($Probability\ of\ \C_{min}\ \geq 80 \, \mg/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 1") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day12714, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day12) 
PTA_CKDEPI80_DAY1

#### Day 14 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_CKDEPI80_DAY14 <-  ggplot(PTA_CKDEPI[PTA_CKDEPI$DAY == 14 & PTA_CKDEPI$Regimen %in% c(6, 2, 5, 4), ], 
                               aes(x = CKDEPI, 
                                   y = PTA_Cmin_80, 
                                   group = factor(Regimen), 
                                   color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(5, 215), breaks = seq(5, 215, by = 10), expand = c(0,0)) +
  xlab(TeX(r"($\eGFR_{CKD-EPI}\ (ml/min/1.73m^2)$)")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) +
  ylab(TeX(r"($Probability\ of\ \C_{min}\ \geq 80 \, \mg/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 14") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day12714, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day714)  
PTA_CKDEPI80_DAY14

#### Combining these 4 plots into 1 -------------------------------------------

## Combining DAY 1 & 14 into 1

# Combine the two plots
PTA_CKDEPI80_DAY114 <- ggarrange(PTA_CKDEPI80_DAY1,
                                 PTA_CKDEPI80_DAY14,
                                 labels = c("a", "b"),
                                 ncol = 1,
                                 nrow = 2,
                                 font.label = list(size = 8, face = "bold"))

# Export the combined plots
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Plots/Dose_finding_simulations/CKDEPI")

# SVG
ggsave("PTA_CKDEPI80_DAY114.svg", 
       PTA_CKDEPI80_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

# JPEG
ggsave("PTA_CKDEPI80_DAY114.JPEG", 
       PTA_CKDEPI80_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

# EPS
ggsave("PTA_CKDEPI80_DAY114.EPS", 
       PTA_CKDEPI80_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

## Secondly, BW non-CRRT effect -------------------------------------------

### Creating datasets for PTA plots -------------------------------------------

#### Standard dosing PTA -------------------------------------------

# Use read_nonmem_table function to import data
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/'
directory_BW_noCRRT_01 <- paste0(working.directory, '/sim_BW_noCRRT_01.1.npctab.dta') 
BW_noCRRT_01 <- read_nonmem_table(directory_BW_noCRRT_01)

# Extract only the trough events
BW_noCRRT_dos_01 <- BW_noCRRT_01[BW_noCRRT_01$TROUGH == 1,]

# Convert ID into numeric
BW_noCRRT_dos_01$ID <- as.integer(BW_noCRRT_dos_01$ID)

# Create a new variable AUC24 for each ID
BW_noCRRT_dos_01 <- BW_noCRRT_dos_01 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

BW_noCRRT_dos_01 <- BW_noCRRT_dos_01 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
BW_noCRRT_dos_01$fAUC <- BW_noCRRT_dos_01$AUC24 * 0.89
 
# Extract ID, DV, fAUC, DAY, and BW columns
PTA_noCRRT_01 <- subset(BW_noCRRT_dos_01, select = c("ID", "DV", "fAUC", "DAY", "BW"))

# Calculate percentage of DV >= 7.5 per DAY per BW (PTA_Cmin_75)
PTA_noCRRT_01$PTA_Cmin_75 <- 100 * ave(PTA_noCRRT_01$DV >= 7.5, PTA_noCRRT_01$DAY, PTA_noCRRT_01$BW, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per BW (PTA_Cmin_80)
PTA_noCRRT_01$PTA_Cmin_80 <- 100 * ave(PTA_noCRRT_01$DV >= 80, PTA_noCRRT_01$DAY, PTA_noCRRT_01$BW, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per BW (PTA_fAUC_200)
PTA_noCRRT_01$PTA_fAUC_200 <- 100 * ave(PTA_noCRRT_01$fAUC >= 200, PTA_noCRRT_01$DAY, PTA_noCRRT_01$BW, FUN = mean)

# Calculate percentage of fAUC >= 400 per DAY per BW (PTA_fAUC_400) #Isabel'comment on simulating for MIC of 4 mg/L
PTA_noCRRT_01$PTA_fAUC_400 <- 100 * ave(PTA_noCRRT_01$fAUC >= 400, PTA_noCRRT_01$DAY, PTA_noCRRT_01$BW, FUN = mean)

# Create PTA_noCRRT_01_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and BW
PTA_noCRRT_01_overall <- unique(PTA_noCRRT_01[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200","PTA_fAUC_400", "DAY", "BW")])

# Export PTA_noCRRT_dos_01 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/")
write.csv(PTA_noCRRT_01_overall, "PTA_noCRRT_dos_01.csv", quote = F, row.names = FALSE)

# Export BW_noCRRT_dos_01.csv for future use
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/")
write.csv(BW_noCRRT_dos_01, "BW_noCRRT_dos_01.csv", quote = F, row.names = FALSE)

#### 1000 LD 400 MD PTA -------------------------------------------

# Use read_nonmem_table function to import data
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/'
directory_BW_noCRRT_02 <- paste0(working.directory, '/sim_BW_noCRRT_02.1.npctab.dta') 
BW_noCRRT_02 <- read_nonmem_table(directory_BW_noCRRT_02)

# Extract only the trough events
BW_noCRRT_dos_02 <- BW_noCRRT_02[BW_noCRRT_02$TROUGH == 1,]

# Convert ID into numeric
BW_noCRRT_dos_02$ID <- as.integer(BW_noCRRT_dos_02$ID)

# Create a new variable AUC24 for each ID
BW_noCRRT_dos_02 <- BW_noCRRT_dos_02 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

BW_noCRRT_dos_02 <- BW_noCRRT_dos_02 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
BW_noCRRT_dos_02$fAUC <- BW_noCRRT_dos_02$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and BW columns
PTA_noCRRT_02 <- subset(BW_noCRRT_dos_02, select = c("ID", "DV", "fAUC", "DAY", "BW"))

# Calculate percentage of DV >= 7.5 per DAY per BW (PTA_Cmin_75)
PTA_noCRRT_02$PTA_Cmin_75 <- 100 * ave(PTA_noCRRT_02$DV >= 7.5, PTA_noCRRT_02$DAY, PTA_noCRRT_02$BW, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per BW (PTA_Cmin_80)
PTA_noCRRT_02$PTA_Cmin_80 <- 100 * ave(PTA_noCRRT_02$DV >= 80, PTA_noCRRT_02$DAY, PTA_noCRRT_02$BW, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per BW (PTA_fAUC_200)
PTA_noCRRT_02$PTA_fAUC_200 <- 100 * ave(PTA_noCRRT_02$fAUC >= 200, PTA_noCRRT_02$DAY, PTA_noCRRT_02$BW, FUN = mean)

# Calculate percentage of fAUC >= 400 per DAY per BW (PTA_fAUC_400) #Isabel'comment on simulating for MIC of 4 mg/L
PTA_noCRRT_02$PTA_fAUC_400 <- 100 * ave(PTA_noCRRT_02$fAUC >= 400, PTA_noCRRT_02$DAY, PTA_noCRRT_02$BW, FUN = mean)

# Create PTA_noCRRT_02_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and BW
PTA_noCRRT_02_overall <- unique(PTA_noCRRT_02[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200","PTA_fAUC_400", "DAY", "BW")])

# Export PTA_noCRRT_dos_02 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/")
write.csv(PTA_noCRRT_02_overall, "PTA_noCRRT_dos_02.csv", quote = F, row.names = FALSE)

# Export BW_noCRRT_dos_02.csv for future use
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/")
write.csv(BW_noCRRT_dos_02, "BW_noCRRT_dos_02.csv", quote = F, row.names = FALSE)

#### 1200 LD 400 MD PTA -------------------------------------------

# Use read_nonmem_table function to import data
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/'
directory_BW_noCRRT_03 <- paste0(working.directory, '/sim_BW_noCRRT_03.1.npctab.dta') 
BW_noCRRT_03 <- read_nonmem_table(directory_BW_noCRRT_03)

# Extract only the trough events
BW_noCRRT_dos_03 <- BW_noCRRT_03[BW_noCRRT_03$TROUGH == 1,]

# Convert ID into numeric
BW_noCRRT_dos_03$ID <- as.integer(BW_noCRRT_dos_03$ID)

# Create a new variable AUC24 for each ID
BW_noCRRT_dos_03 <- BW_noCRRT_dos_03 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

BW_noCRRT_dos_03 <- BW_noCRRT_dos_03 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
BW_noCRRT_dos_03$fAUC <- BW_noCRRT_dos_03$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and BW columns
PTA_noCRRT_03 <- subset(BW_noCRRT_dos_03, select = c("ID", "DV", "fAUC", "DAY", "BW"))

# Calculate percentage of DV >= 7.5 per DAY per BW (PTA_Cmin_75)
PTA_noCRRT_03$PTA_Cmin_75 <- 100 * ave(PTA_noCRRT_03$DV >= 7.5, PTA_noCRRT_03$DAY, PTA_noCRRT_03$BW, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per BW (PTA_Cmin_80)
PTA_noCRRT_03$PTA_Cmin_80 <- 100 * ave(PTA_noCRRT_03$DV >= 80, PTA_noCRRT_03$DAY, PTA_noCRRT_03$BW, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per BW (PTA_fAUC_200)
PTA_noCRRT_03$PTA_fAUC_200 <- 100 * ave(PTA_noCRRT_03$fAUC >= 200, PTA_noCRRT_03$DAY, PTA_noCRRT_03$BW, FUN = mean)

# Calculate percentage of fAUC >= 400 per DAY per BW (PTA_fAUC_400) #Isabel'comment on simulating for MIC of 4 mg/L
PTA_noCRRT_03$PTA_fAUC_400 <- 100 * ave(PTA_noCRRT_03$fAUC >= 400, PTA_noCRRT_03$DAY, PTA_noCRRT_03$BW, FUN = mean)

# Create PTA_noCRRT_03_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and BW
PTA_noCRRT_03_overall <- unique(PTA_noCRRT_03[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200","PTA_fAUC_400", "DAY", "BW")])

# Export PTA_noCRRT_dos_03 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/")
write.csv(PTA_noCRRT_03_overall, "PTA_noCRRT_dos_03.csv", quote = F, row.names = FALSE)

# Export BW_noCRRT_dos_03.csv for future use
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/")
write.csv(BW_noCRRT_dos_03, "BW_noCRRT_dos_03.csv", quote = F, row.names = FALSE)

#### 1400 LD 400 MD PTA -------------------------------------------

# Use read_nonmem_table function to import data
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/'
directory_BW_noCRRT_04 <- paste0(working.directory, '/sim_BW_noCRRT_04.1.npctab.dta') 
BW_noCRRT_04 <- read_nonmem_table(directory_BW_noCRRT_04)

# Extract only the trough events
BW_noCRRT_dos_04 <- BW_noCRRT_04[BW_noCRRT_04$TROUGH == 1,]

# Convert ID into numeric
BW_noCRRT_dos_04$ID <- as.integer(BW_noCRRT_dos_04$ID)

# Create a new variable AUC24 for each ID
BW_noCRRT_dos_04 <- BW_noCRRT_dos_04 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

BW_noCRRT_dos_04 <- BW_noCRRT_dos_04 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
BW_noCRRT_dos_04$fAUC <- BW_noCRRT_dos_04$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and BW columns
PTA_noCRRT_04 <- subset(BW_noCRRT_dos_04, select = c("ID", "DV", "fAUC", "DAY", "BW"))

# Calculate percentage of DV >= 7.5 per DAY per BW (PTA_Cmin_75)
PTA_noCRRT_04$PTA_Cmin_75 <- 100 * ave(PTA_noCRRT_04$DV >= 7.5, PTA_noCRRT_04$DAY, PTA_noCRRT_04$BW, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per BW (PTA_Cmin_80)
PTA_noCRRT_04$PTA_Cmin_80 <- 100 * ave(PTA_noCRRT_04$DV >= 80, PTA_noCRRT_04$DAY, PTA_noCRRT_04$BW, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per BW (PTA_fAUC_200)
PTA_noCRRT_04$PTA_fAUC_200 <- 100 * ave(PTA_noCRRT_04$fAUC >= 200, PTA_noCRRT_04$DAY, PTA_noCRRT_04$BW, FUN = mean)

# Calculate percentage of fAUC >= 400 per DAY per BW (PTA_fAUC_400) #Isabel'comment on simulating for MIC of 4 mg/L
PTA_noCRRT_04$PTA_fAUC_400 <- 100 * ave(PTA_noCRRT_04$fAUC >= 400, PTA_noCRRT_04$DAY, PTA_noCRRT_04$BW, FUN = mean)

# Create PTA_noCRRT_04_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and BW
PTA_noCRRT_04_overall <- unique(PTA_noCRRT_04[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200","PTA_fAUC_400", "DAY", "BW")])

# Export PTA_noCRRT_dos_04 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/")
write.csv(PTA_noCRRT_04_overall, "PTA_noCRRT_dos_04.csv", quote = F, row.names = FALSE)

# Export BW_noCRRT_dos_04.csv for future use
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/")
write.csv(BW_noCRRT_dos_04, "BW_noCRRT_dos_04.csv", quote = F, row.names = FALSE)

#### 1600 LD 400 MD PTA -------------------------------------------

# Use read_nonmem_table function to import data
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/'
directory_BW_noCRRT_05 <- paste0(working.directory, '/sim_BW_noCRRT_05.1.npctab.dta') 
BW_noCRRT_05 <- read_nonmem_table(directory_BW_noCRRT_05)

# Extract only the trough events
BW_noCRRT_dos_05 <- BW_noCRRT_05[BW_noCRRT_05$TROUGH == 1,]

# Convert ID into numeric
BW_noCRRT_dos_05$ID <- as.integer(BW_noCRRT_dos_05$ID)

# Create a new variable AUC24 for each ID
BW_noCRRT_dos_05 <- BW_noCRRT_dos_05 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

BW_noCRRT_dos_05 <- BW_noCRRT_dos_05 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
BW_noCRRT_dos_05$fAUC <- BW_noCRRT_dos_05$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and BW columns
PTA_noCRRT_05 <- subset(BW_noCRRT_dos_05, select = c("ID", "DV", "fAUC", "DAY", "BW"))

# Calculate percentage of DV >= 7.5 per DAY per BW (PTA_Cmin_75)
PTA_noCRRT_05$PTA_Cmin_75 <- 100 * ave(PTA_noCRRT_05$DV >= 7.5, PTA_noCRRT_05$DAY, PTA_noCRRT_05$BW, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per BW (PTA_Cmin_80)
PTA_noCRRT_05$PTA_Cmin_80 <- 100 * ave(PTA_noCRRT_05$DV >= 80, PTA_noCRRT_05$DAY, PTA_noCRRT_05$BW, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per BW (PTA_fAUC_200)
PTA_noCRRT_05$PTA_fAUC_200 <- 100 * ave(PTA_noCRRT_05$fAUC >= 200, PTA_noCRRT_05$DAY, PTA_noCRRT_05$BW, FUN = mean)

# Calculate percentage of fAUC >= 400 per DAY per BW (PTA_fAUC_400) #Isabel'comment on simulating for MIC of 4 mg/L
PTA_noCRRT_05$PTA_fAUC_400 <- 100 * ave(PTA_noCRRT_05$fAUC >= 400, PTA_noCRRT_05$DAY, PTA_noCRRT_05$BW, FUN = mean)

# Create PTA_noCRRT_05_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and BW
PTA_noCRRT_05_overall <- unique(PTA_noCRRT_05[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200","PTA_fAUC_400", "DAY", "BW")])

# Export PTA_noCRRT_dos_05 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/")
write.csv(PTA_noCRRT_05_overall, "PTA_noCRRT_dos_05.csv", quote = F, row.names = FALSE)

# Export BW_noCRRT_dos_05.csv for future use
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/")
write.csv(BW_noCRRT_dos_05, "BW_noCRRT_dos_05.csv", quote = F, row.names = FALSE)

#### 1800 LD 400 MD PTA -------------------------------------------

# Use read_nonmem_table function to import data
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/'
directory_BW_noCRRT_06 <- paste0(working.directory, '/sim_BW_noCRRT_06.1.npctab.dta') 
BW_noCRRT_06 <- read_nonmem_table(directory_BW_noCRRT_06)

# Extract only the trough events
BW_noCRRT_dos_06 <- BW_noCRRT_06[BW_noCRRT_06$TROUGH == 1,]

# Convert ID into numeric
BW_noCRRT_dos_06$ID <- as.integer(BW_noCRRT_dos_06$ID)
 
# Create a new variable AUC24 for each ID
BW_noCRRT_dos_06 <- BW_noCRRT_dos_06 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

BW_noCRRT_dos_06 <- BW_noCRRT_dos_06 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
BW_noCRRT_dos_06$fAUC <- BW_noCRRT_dos_06$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and BW columns
PTA_noCRRT_06 <- subset(BW_noCRRT_dos_06, select = c("ID", "DV", "fAUC", "DAY", "BW"))

# Calculate percentage of DV >= 7.5 per DAY per BW (PTA_Cmin_75)
PTA_noCRRT_06$PTA_Cmin_75 <- 100 * ave(PTA_noCRRT_06$DV >= 7.5, PTA_noCRRT_06$DAY, PTA_noCRRT_06$BW, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per BW (PTA_Cmin_80)
PTA_noCRRT_06$PTA_Cmin_80 <- 100 * ave(PTA_noCRRT_06$DV >= 80, PTA_noCRRT_06$DAY, PTA_noCRRT_06$BW, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per BW (PTA_fAUC_200)
PTA_noCRRT_06$PTA_fAUC_200 <- 100 * ave(PTA_noCRRT_06$fAUC >= 200, PTA_noCRRT_06$DAY, PTA_noCRRT_06$BW, FUN = mean)

# Calculate percentage of fAUC >= 400 per DAY per BW (PTA_fAUC_400) #Isabel'comment on simulating for MIC of 4 mg/L
PTA_noCRRT_06$PTA_fAUC_400 <- 100 * ave(PTA_noCRRT_06$fAUC >= 400, PTA_noCRRT_06$DAY, PTA_noCRRT_06$BW, FUN = mean)

# Create PTA_noCRRT_06_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and BW
PTA_noCRRT_06_overall <- unique(PTA_noCRRT_06[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200","PTA_fAUC_400", "DAY", "BW")])

# Export PTA_noCRRT_dos_06 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/")
write.csv(PTA_noCRRT_06_overall, "PTA_noCRRT_dos_06.csv", quote = F, row.names = FALSE)

# Export BW_noCRRT_dos_06.csv for future use
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/")
write.csv(BW_noCRRT_dos_06, "BW_noCRRT_dos_06.csv", quote = F, row.names = FALSE)

#### 12 mg/kg LD 6 mg/kg MD PTA -------------------------------------------

# Use read_nonmem_table function to import data
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/'
directory_BW_noCRRT_07 <- paste0(working.directory, '/sim_BW_noCRRT_07.1.npctab.dta') 
BW_noCRRT_07 <- read_nonmem_table(directory_BW_noCRRT_07)

# Extract only the trough events
BW_noCRRT_dos_07 <- BW_noCRRT_07[BW_noCRRT_07$TROUGH == 1,]

# Convert ID into numeric
BW_noCRRT_dos_07$ID <- as.integer(BW_noCRRT_dos_07$ID)

# Create a new variable AUC24 for each ID
BW_noCRRT_dos_07 <- BW_noCRRT_dos_07 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

BW_noCRRT_dos_07 <- BW_noCRRT_dos_07 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
BW_noCRRT_dos_07$fAUC <- BW_noCRRT_dos_07$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and BW columns
PTA_noCRRT_07 <- subset(BW_noCRRT_dos_07, select = c("ID", "DV", "fAUC", "DAY", "BW"))

# Calculate percentage of DV >= 7.5 per DAY per BW (PTA_Cmin_75)
PTA_noCRRT_07$PTA_Cmin_75 <- 100 * ave(PTA_noCRRT_07$DV >= 7.5, PTA_noCRRT_07$DAY, PTA_noCRRT_07$BW, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per BW (PTA_Cmin_80)
PTA_noCRRT_07$PTA_Cmin_80 <- 100 * ave(PTA_noCRRT_07$DV >= 80, PTA_noCRRT_07$DAY, PTA_noCRRT_07$BW, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per BW (PTA_fAUC_200)
PTA_noCRRT_07$PTA_fAUC_200 <- 100 * ave(PTA_noCRRT_07$fAUC >= 200, PTA_noCRRT_07$DAY, PTA_noCRRT_07$BW, FUN = mean)

# Calculate percentage of fAUC >= 400 per DAY per BW (PTA_fAUC_400) #Isabel'comment on simulating for MIC of 4 mg/L
PTA_noCRRT_07$PTA_fAUC_400 <- 100 * ave(PTA_noCRRT_07$fAUC >= 400, PTA_noCRRT_07$DAY, PTA_noCRRT_07$BW, FUN = mean)

# Create PTA_noCRRT_07_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and BW
PTA_noCRRT_07_overall <- unique(PTA_noCRRT_07[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200","PTA_fAUC_400", "DAY", "BW")])

# Export PTA_noCRRT_dos_07 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/")
write.csv(PTA_noCRRT_07_overall, "PTA_noCRRT_dos_07.csv", quote = F, row.names = FALSE)

# Export BW_noCRRT_dos_07.csv for future use
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/")
write.csv(BW_noCRRT_dos_07, "BW_noCRRT_dos_07.csv", quote = F, row.names = FALSE)

#### 800 LD 200 MD PTA -------------------------------------------

# Use read_nonmem_table function to import data
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/'
directory_BW_noCRRT_08 <- paste0(working.directory, '/sim_BW_noCRRT_08.1.npctab.dta') 
BW_noCRRT_08 <- read_nonmem_table(directory_BW_noCRRT_08)

# Extract only the trough events
BW_noCRRT_dos_08 <- BW_noCRRT_08[BW_noCRRT_08$TROUGH == 1,]

# Convert ID into numeric
BW_noCRRT_dos_08$ID <- as.integer(BW_noCRRT_dos_08$ID)

# Create a new variable AUC24 for each ID
BW_noCRRT_dos_08 <- BW_noCRRT_dos_08 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

BW_noCRRT_dos_08 <- BW_noCRRT_dos_08 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
BW_noCRRT_dos_08$fAUC <- BW_noCRRT_dos_08$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and BW columns
PTA_noCRRT_08 <- subset(BW_noCRRT_dos_08, select = c("ID", "DV", "fAUC", "DAY", "BW"))

# Calculate percentage of DV >= 7.5 per DAY per BW (PTA_Cmin_75)
PTA_noCRRT_08$PTA_Cmin_75 <- 100 * ave(PTA_noCRRT_08$DV >= 7.5, PTA_noCRRT_08$DAY, PTA_noCRRT_08$BW, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per BW (PTA_Cmin_80)
PTA_noCRRT_08$PTA_Cmin_80 <- 100 * ave(PTA_noCRRT_08$DV >= 80, PTA_noCRRT_08$DAY, PTA_noCRRT_08$BW, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per BW (PTA_fAUC_200)
PTA_noCRRT_08$PTA_fAUC_200 <- 100 * ave(PTA_noCRRT_08$fAUC >= 200, PTA_noCRRT_08$DAY, PTA_noCRRT_08$BW, FUN = mean)

# Calculate percentage of fAUC >= 400 per DAY per BW (PTA_fAUC_400) #Isabel'comment on simulating for MIC of 4 mg/L
PTA_noCRRT_08$PTA_fAUC_400 <- 100 * ave(PTA_noCRRT_08$fAUC >= 400, PTA_noCRRT_08$DAY, PTA_noCRRT_08$BW, FUN = mean)

# Create PTA_noCRRT_08_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and BW
PTA_noCRRT_08_overall <- unique(PTA_noCRRT_08[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200","PTA_fAUC_400", "DAY", "BW")])

# Export PTA_noCRRT_dos_08 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/")
write.csv(PTA_noCRRT_08_overall, "PTA_noCRRT_dos_08.csv", quote = F, row.names = FALSE)

# Export BW_noCRRT_dos_08.csv for future use
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_noCRRT/")
write.csv(BW_noCRRT_dos_08, "BW_noCRRT_dos_08.csv", quote = F, row.names = FALSE)

### Making PTA plots -------------------------------------------

#### Some preps -------------------------------------------

# Import PTA_dos_ datasets
PTA_noCRRT_dos_01 <- read.csv("PTA_BW_noCRRT/PTA_noCRRT_dos_01.csv")
PTA_noCRRT_dos_02 <- read.csv("PTA_BW_noCRRT/PTA_noCRRT_dos_02.csv")
PTA_noCRRT_dos_03 <- read.csv("PTA_BW_noCRRT/PTA_noCRRT_dos_03.csv")
PTA_noCRRT_dos_04 <- read.csv("PTA_BW_noCRRT/PTA_noCRRT_dos_04.csv")
PTA_noCRRT_dos_05 <- read.csv("PTA_BW_noCRRT/PTA_noCRRT_dos_05.csv")
PTA_noCRRT_dos_06 <- read.csv("PTA_BW_noCRRT/PTA_noCRRT_dos_06.csv")
PTA_noCRRT_dos_07 <- read.csv("PTA_BW_noCRRT/PTA_noCRRT_dos_07.csv")
PTA_noCRRT_dos_08 <- read.csv("PTA_BW_noCRRT/PTA_noCRRT_dos_08.csv")

# Add regimens
PTA_noCRRT_dos_01$Regimen <- 1
PTA_noCRRT_dos_02$Regimen <- 2
PTA_noCRRT_dos_03$Regimen <- 3
PTA_noCRRT_dos_04$Regimen <- 4
PTA_noCRRT_dos_05$Regimen <- 5
PTA_noCRRT_dos_06$Regimen <- 6
PTA_noCRRT_dos_07$Regimen <- 7
PTA_noCRRT_dos_08$Regimen <- 8

# Combine 6 PTA datasets into 1 - PTA_BW_noCRRT
PTA_BW_noCRRT <- rbind(PTA_noCRRT_dos_01, PTA_noCRRT_dos_02, PTA_noCRRT_dos_03, 
                    PTA_noCRRT_dos_04, PTA_noCRRT_dos_05, PTA_noCRRT_dos_06,
                    PTA_noCRRT_dos_07, PTA_noCRRT_dos_08)

# Define the labels and colors
dosing_labels_day12 <- c("800 mg q24h LD",
                         "1000 mg q24h LD",
                         "1200 mg q24h LD",
                         "1400 mg q24h LD",
                         "1600 mg q24h LD",
                         "12 mg/kg q24h LD")
# Regimens 8, 2, 3, 4, 5, 6, 7, respectively

dosing_labels_day714 <- c("200 mg q24h MD", 
                          "400 mg q24h MD", 
                          "6 mg/kg q24h MD")
# Regimens 8, 2, 7, respectively

dosing_colors_day12 <- c("#a0da39",
                         "#4ac16d",
                         "#1fa187",
                         "#277f8e",
                         "#365c8d",
                         "#46327e",
                         "#440154")

dosing_colors_day714 <- c("#a0da39",
                          "#4ac16d",
                          "#440154")

# Reorder the levels of the dosing regimens
PTA_BW_noCRRT$Regimen <- factor(PTA_BW_noCRRT$Regimen,
                             levels = c(8, 2, 3, 4, 5, 6, 7, 1))

# It's not necessary to have regimen 6 after I did some explorations

# I also enforce limits on vertical lines here

# Firstly, create dataset with BW, Regimen and PTA on day 1 
BW_Reg <- data.frame( #BW and Regimen dataset
  BW = c(60, 80, 100, 120),
  Regimen = c(8, 2, 3, 4)
)

BW_Reg$Regimen <- factor(BW_Reg$Regimen,
                         levels = c(8, 2, 3, 4))

# Then, make a dataset of PTA on Day 1 per BW per Regimen
BW_PTA_DaY1 <- BW_Reg %>%
  left_join(
    PTA_BW_noCRRT %>% 
      dplyr::filter(DAY == 1),
    by = c("BW", "Regimen")
  )

# Similarly for day 2
BW_PTA_DaY2 <- BW_Reg %>%
  left_join(
    PTA_BW_noCRRT %>% 
      dplyr::filter(DAY == 2),
    by = c("BW", "Regimen")
  )

#### Day 1 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_BW_noCRRT200_DAY1 <-  ggplot(PTA_BW_noCRRT[PTA_BW_noCRRT$DAY == 1 & !PTA_BW_noCRRT$Regimen %in% c(1, 6), ], 
                              aes(x = BW, 
                                  y = PTA_fAUC_200, 
                                  group = factor(Regimen), 
                                  color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(30, 150), breaks = seq(30, 150, by = 5), expand = c(0,0)) +
  xlab(TeX(r"(Total body weight (kg))")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) +
  ylab(TeX(r"($Probability\ of\ \textit{f}AUC_{0-24}\ \geq 200 \, \mg \times h/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 1") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day12, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day12) +
  geom_hline(yintercept = 90, 
             linetype = "dashed", 
             color = "gray50", 
             size = 0.6) +
  geom_segment(
    data = BW_PTA_DaY1,
    aes(
      x = BW,
      xend = BW,
      y = 0,  # Starting from the x-axis
      yend = PTA_fAUC_200  # Ending at the corresponding PTA_fAUC_200 value
    ),
    linetype = "dashed",
    color = "gray50",
    size = 0.6
  )
PTA_BW_noCRRT200_DAY1

#### Day 2 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_BW_noCRRT200_DAY2 <-  ggplot(PTA_BW_noCRRT[PTA_BW_noCRRT$DAY == 2 & !PTA_BW_noCRRT$Regimen %in% c(1, 6), ], 
                                 aes(x = BW, 
                                     y = PTA_fAUC_200, 
                                     group = factor(Regimen), 
                                     color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(30, 150), breaks = seq(30, 150, by = 5), expand = c(0,0)) +
  xlab(TeX(r"(Total body weight (kg))")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) +
  ylab(TeX(r"($Probability\ of\ \textit{f}AUC_{0-24}\ \geq 200 \, \mg \times h/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 2") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day12, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day12) +
  geom_hline(yintercept = 90, 
             linetype = "dashed", 
             color = "gray50", 
             size = 0.6) +
  geom_segment(
    data = BW_PTA_DaY2,
    aes(
      x = BW,
      xend = BW,
      y = 0,  # Starting from the x-axis
      yend = PTA_fAUC_200  # Ending at the corresponding PTA_fAUC_200 value
    ),
    linetype = "dashed",
    color = "gray50",
    size = 0.6
  )
PTA_BW_noCRRT200_DAY2

#### Day 14 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_BW_noCRRT200_DAY14 <-  ggplot(PTA_BW_noCRRT[PTA_BW_noCRRT$DAY == 14 & PTA_BW_noCRRT$Regimen %in% c(8, 2, 7), ], 
                                 aes(x = BW, 
                                     y = PTA_fAUC_200, 
                                     group = factor(Regimen), 
                                     color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(30, 150), breaks = seq(30, 150, by = 5), expand = c(0,0)) +
  xlab(TeX(r"(Total body weight (kg))")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) +
  ylab(TeX(r"($Probability\ of\ \textit{f}AUC_{0-24}\ \geq 200 \, \mg \times h/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 14") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day714, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day714) +
  geom_hline(yintercept = 90, 
             linetype = "dashed", 
             color = "gray50", 
             size = 0.6) 
PTA_BW_noCRRT200_DAY14

#### Day 7 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_BW_noCRRT200_DAY7 <-  ggplot(PTA_BW_noCRRT[PTA_BW_noCRRT$DAY == 7 & PTA_BW_noCRRT$Regimen %in% c(8, 2, 7), ], 
                                  aes(x = BW, 
                                      y = PTA_fAUC_200, 
                                      group = factor(Regimen), 
                                      color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(30, 150), breaks = seq(30, 150, by = 5), expand = c(0,0)) +
  xlab(TeX(r"(Total body weight (kg))")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) +
  ylab(TeX(r"($Probability\ of\ \textit{f}AUC_{0-24}\ \geq 200 \, \mg \times h/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 7") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day714, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day714) +
  geom_hline(yintercept = 90, 
             linetype = "dashed", 
             color = "gray50", 
             size = 0.6) 
PTA_BW_noCRRT200_DAY7

#### Combining these 4 plots into 1 -------------------------------------------

## Combining DAY 1 & 14 into 1

# Combine the two plots
PTA_BW_noCRRT200_DAY114 <- ggarrange(PTA_BW_noCRRT200_DAY1, 
                                               PTA_BW_noCRRT200_DAY14,
                                               labels = c("a", "b"),
                                               ncol = 1, 
                                               nrow = 2,
                                               font.label = list(size = 8, face = "bold"))

# Export the combined plots
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Plots/Dose_finding_simulations/BW_noCRRT")

# SVG
ggsave("PTA_BW_noCRRT200_DAY114.svg", 
       PTA_BW_noCRRT200_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

# JPEG
ggsave("PTA_BW_noCRRT200_DAY114.JPEG", 
       PTA_BW_noCRRT200_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

# EPS
ggsave("PTA_BW_noCRRT200_DAY114.EPS", 
       PTA_BW_noCRRT200_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

### Assessing toxicity PTA Cmin80 -------------------------------------------

# PTA_CKDEPI dataset was created above

#### Day 1 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_BW_noCRRT80_DAY1 <-  ggplot(PTA_BW_noCRRT[PTA_BW_noCRRT$DAY == 1 & !PTA_BW_noCRRT$Regimen %in% c(1, 6), ], 
                             aes(x = BW, 
                                 y = PTA_Cmin_80, 
                                 group = factor(Regimen), 
                                 color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(30, 150), breaks = seq(30, 150, by = 5), expand = c(0,0)) +
  xlab(TeX(r"(Total body weight (kg))")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) +
  ylab(TeX(r"($Probability\ of\ \C_{min}\ \geq 80 \, \mg/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 1") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day12, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day12) +
  geom_segment(
    data = BW_PTA_DaY1,
    aes(
      x = BW,
      xend = BW,
      y = 0,  # Starting from the x-axis
      yend = PTA_Cmin_80  # Ending at the corresponding PTA_Cmin_80 value
    ),
    linetype = "dashed",
    color = "gray50",
    size = 0.6
  )
PTA_BW_noCRRT80_DAY1

#### Day 14 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_BW_noCRRT80_DAY14 <-  ggplot(PTA_BW_noCRRT[PTA_BW_noCRRT$DAY == 14 & PTA_BW_noCRRT$Regimen %in% c(8, 2, 7), ], 
                                aes(x = BW, 
                                    y = PTA_Cmin_80, 
                                    group = factor(Regimen), 
                                    color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(30, 150), breaks = seq(30, 150, by = 5), expand = c(0,0)) +
  xlab(TeX(r"(Total body weight (kg))")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) +
  ylab(TeX(r"($Probability\ of\ \C_{min}\ \geq 80 \, \mg/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 14") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day714, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day714) 
PTA_BW_noCRRT80_DAY14

#### Combining these 4 plots into 1 -------------------------------------------

## Combining DAY 1 & 14 into 1

# Combine the two plots
PTA_BW_noCRRT80_DAY114 <- ggarrange(PTA_BW_noCRRT80_DAY1,
                                    PTA_BW_noCRRT80_DAY14,
                                    labels = c("a", "b"),
                                    ncol = 1,
                                    nrow = 2,
                                    font.label = list(size = 8, face = "bold"))

# Export the combined plots
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Plots/Dose_finding_simulations/BW_noCRRT")

# SVG
ggsave("PTA_BW_noCRRT80_DAY114.svg", 
       PTA_BW_noCRRT80_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

# JPEG
ggsave("PTA_BW_noCRRT80_DAY114.JPEG", 
       PTA_BW_noCRRT80_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

# EPS
ggsave("PTA_BW_noCRRT80_DAY114.EPS", 
       PTA_BW_noCRRT80_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

## Third, BW CRRT effect -------------------------------------------

### Creating datasets for PTA plots -------------------------------------------

#### Standard dosing PTA -------------------------------------------

# Use read_nonmem_table function to import data
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/'
directory_BW_CRRT_01 <- paste0(working.directory, '/sim_BW_CRRT_01.1.npctab.dta') 
BW_CRRT_01 <- read_nonmem_table(directory_BW_CRRT_01)

# Extract only the trough events
BW_CRRT_dos_01 <- BW_CRRT_01[BW_CRRT_01$TROUGH == 1,]

# Convert ID into numeric
BW_CRRT_dos_01$ID <- as.integer(BW_CRRT_dos_01$ID)

# Create a new variable AUC24 for each ID
BW_CRRT_dos_01 <- BW_CRRT_dos_01 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

BW_CRRT_dos_01 <- BW_CRRT_dos_01 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
BW_CRRT_dos_01$fAUC <- BW_CRRT_dos_01$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and BW columns
PTA_CRRT_01 <- subset(BW_CRRT_dos_01, select = c("ID", "DV", "fAUC", "DAY", "BW"))

# Calculate percentage of DV >= 7.5 per DAY per BW (PTA_Cmin_75)
PTA_CRRT_01$PTA_Cmin_75 <- 100 * ave(PTA_CRRT_01$DV >= 7.5, PTA_CRRT_01$DAY, PTA_CRRT_01$BW, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per BW (PTA_Cmin_80)
PTA_CRRT_01$PTA_Cmin_80 <- 100 * ave(PTA_CRRT_01$DV >= 80, PTA_CRRT_01$DAY, PTA_CRRT_01$BW, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per BW (PTA_fAUC_200)
PTA_CRRT_01$PTA_fAUC_200 <- 100 * ave(PTA_CRRT_01$fAUC >= 200, PTA_CRRT_01$DAY, PTA_CRRT_01$BW, FUN = mean)

# Calculate percentage of fAUC >= 400 per DAY per BW (PTA_fAUC_400) #Isabel'comment on simulating for MIC of 4 mg/L
PTA_CRRT_01$PTA_fAUC_400 <- 100 * ave(PTA_CRRT_01$fAUC >= 400, PTA_CRRT_01$DAY, PTA_CRRT_01$BW, FUN = mean)

# Create PTA_CRRT_01_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and BW
PTA_CRRT_01_overall <- unique(PTA_CRRT_01[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200","PTA_fAUC_400", "DAY", "BW")])

# Export PTA_CRRT_dos_01 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/")
write.csv(PTA_CRRT_01_overall, "PTA_CRRT_dos_01.csv", quote = F, row.names = FALSE)

# Export BW_CRRT_dos_01.csv for future use
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/")
write.csv(BW_CRRT_dos_01, "BW_CRRT_dos_01.csv", quote = F, row.names = FALSE)

#### 1000 LD 800 MD PTA -------------------------------------------

# Use read_nonmem_table function to import data
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/'
directory_BW_CRRT_02 <- paste0(working.directory, '/sim_BW_CRRT_02.1.npctab.dta') 
BW_CRRT_02 <- read_nonmem_table(directory_BW_CRRT_02)

# Extract only the trough events
BW_CRRT_dos_02 <- BW_CRRT_02[BW_CRRT_02$TROUGH == 1,]

# Convert ID into numeric
BW_CRRT_dos_02$ID <- as.integer(BW_CRRT_dos_02$ID)

# Create a new variable AUC24 for each ID
BW_CRRT_dos_02 <- BW_CRRT_dos_02 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

BW_CRRT_dos_02 <- BW_CRRT_dos_02 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
BW_CRRT_dos_02$fAUC <- BW_CRRT_dos_02$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and BW columns
PTA_CRRT_02 <- subset(BW_CRRT_dos_02, select = c("ID", "DV", "fAUC", "DAY", "BW"))

# Calculate percentage of DV >= 7.5 per DAY per BW (PTA_Cmin_75)
PTA_CRRT_02$PTA_Cmin_75 <- 100 * ave(PTA_CRRT_02$DV >= 7.5, PTA_CRRT_02$DAY, PTA_CRRT_02$BW, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per BW (PTA_Cmin_80)
PTA_CRRT_02$PTA_Cmin_80 <- 100 * ave(PTA_CRRT_02$DV >= 80, PTA_CRRT_02$DAY, PTA_CRRT_02$BW, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per BW (PTA_fAUC_200)
PTA_CRRT_02$PTA_fAUC_200 <- 100 * ave(PTA_CRRT_02$fAUC >= 200, PTA_CRRT_02$DAY, PTA_CRRT_02$BW, FUN = mean)

# Calculate percentage of fAUC >= 400 per DAY per BW (PTA_fAUC_400) #Isabel'comment on simulating for MIC of 4 mg/L
PTA_CRRT_02$PTA_fAUC_400 <- 100 * ave(PTA_CRRT_02$fAUC >= 400, PTA_CRRT_02$DAY, PTA_CRRT_02$BW, FUN = mean)

# Create PTA_CRRT_02_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and BW
PTA_CRRT_02_overall <- unique(PTA_CRRT_02[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200","PTA_fAUC_400", "DAY", "BW")])

# Export PTA_CRRT_dos_02 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/")
write.csv(PTA_CRRT_02_overall, "PTA_CRRT_dos_02.csv", quote = F, row.names = FALSE)

# Export BW_CRRT_dos_02.csv for future use
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/")
write.csv(BW_CRRT_dos_02, "BW_CRRT_dos_02.csv", quote = F, row.names = FALSE)

#### 1200 LD 800 MD PTA -------------------------------------------

# Use read_nonmem_table function to import data
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/'
directory_BW_CRRT_03 <- paste0(working.directory, '/sim_BW_CRRT_03.1.npctab.dta') 
BW_CRRT_03 <- read_nonmem_table(directory_BW_CRRT_03)

# Extract only the trough events
BW_CRRT_dos_03 <- BW_CRRT_03[BW_CRRT_03$TROUGH == 1,]

# Convert ID into numeric
BW_CRRT_dos_03$ID <- as.integer(BW_CRRT_dos_03$ID)

# Create a new variable AUC24 for each ID
BW_CRRT_dos_03 <- BW_CRRT_dos_03 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

BW_CRRT_dos_03 <- BW_CRRT_dos_03 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
BW_CRRT_dos_03$fAUC <- BW_CRRT_dos_03$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and BW columns
PTA_CRRT_03 <- subset(BW_CRRT_dos_03, select = c("ID", "DV", "fAUC", "DAY", "BW"))

# Calculate percentage of DV >= 7.5 per DAY per BW (PTA_Cmin_75)
PTA_CRRT_03$PTA_Cmin_75 <- 100 * ave(PTA_CRRT_03$DV >= 7.5, PTA_CRRT_03$DAY, PTA_CRRT_03$BW, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per BW (PTA_Cmin_80)
PTA_CRRT_03$PTA_Cmin_80 <- 100 * ave(PTA_CRRT_03$DV >= 80, PTA_CRRT_03$DAY, PTA_CRRT_03$BW, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per BW (PTA_fAUC_200)
PTA_CRRT_03$PTA_fAUC_200 <- 100 * ave(PTA_CRRT_03$fAUC >= 200, PTA_CRRT_03$DAY, PTA_CRRT_03$BW, FUN = mean)

# Calculate percentage of fAUC >= 400 per DAY per BW (PTA_fAUC_400) #Isabel'comment on simulating for MIC of 4 mg/L
PTA_CRRT_03$PTA_fAUC_400 <- 100 * ave(PTA_CRRT_03$fAUC >= 400, PTA_CRRT_03$DAY, PTA_CRRT_03$BW, FUN = mean)

# Create PTA_CRRT_03_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and BW
PTA_CRRT_03_overall <- unique(PTA_CRRT_03[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200","PTA_fAUC_400", "DAY", "BW")])

# Export PTA_CRRT_dos_03 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/")
write.csv(PTA_CRRT_03_overall, "PTA_CRRT_dos_03.csv", quote = F, row.names = FALSE)

# Export BW_CRRT_dos_03.csv for future use
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/")
write.csv(BW_CRRT_dos_03, "BW_CRRT_dos_03.csv", quote = F, row.names = FALSE)

#### 1400 LD 800 MD PTA -------------------------------------------

# Use read_nonmem_table function to import data
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/'
directory_BW_CRRT_04 <- paste0(working.directory, '/sim_BW_CRRT_04.1.npctab.dta') 
BW_CRRT_04 <- read_nonmem_table(directory_BW_CRRT_04)

# Extract only the trough events
BW_CRRT_dos_04 <- BW_CRRT_04[BW_CRRT_04$TROUGH == 1,]

# Convert ID into numeric
BW_CRRT_dos_04$ID <- as.integer(BW_CRRT_dos_04$ID)

# Create a new variable AUC24 for each ID
BW_CRRT_dos_04 <- BW_CRRT_dos_04 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

BW_CRRT_dos_04 <- BW_CRRT_dos_04 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
BW_CRRT_dos_04$fAUC <- BW_CRRT_dos_04$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and BW columns
PTA_CRRT_04 <- subset(BW_CRRT_dos_04, select = c("ID", "DV", "fAUC", "DAY", "BW"))

# Calculate percentage of DV >= 7.5 per DAY per BW (PTA_Cmin_75)
PTA_CRRT_04$PTA_Cmin_75 <- 100 * ave(PTA_CRRT_04$DV >= 7.5, PTA_CRRT_04$DAY, PTA_CRRT_04$BW, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per BW (PTA_Cmin_80)
PTA_CRRT_04$PTA_Cmin_80 <- 100 * ave(PTA_CRRT_04$DV >= 80, PTA_CRRT_04$DAY, PTA_CRRT_04$BW, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per BW (PTA_fAUC_200)
PTA_CRRT_04$PTA_fAUC_200 <- 100 * ave(PTA_CRRT_04$fAUC >= 200, PTA_CRRT_04$DAY, PTA_CRRT_04$BW, FUN = mean)

# Calculate percentage of fAUC >= 400 per DAY per BW (PTA_fAUC_400) #Isabel'comment on simulating for MIC of 4 mg/L
PTA_CRRT_04$PTA_fAUC_400 <- 100 * ave(PTA_CRRT_04$fAUC >= 400, PTA_CRRT_04$DAY, PTA_CRRT_04$BW, FUN = mean)

# Create PTA_CRRT_04_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and BW
PTA_CRRT_04_overall <- unique(PTA_CRRT_04[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200","PTA_fAUC_400", "DAY", "BW")])

# Export PTA_CRRT_dos_04 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/")
write.csv(PTA_CRRT_04_overall, "PTA_CRRT_dos_04.csv", quote = F, row.names = FALSE)

# Export BW_CRRT_dos_04.csv for future use
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/")
write.csv(BW_CRRT_dos_04, "BW_CRRT_dos_04.csv", quote = F, row.names = FALSE)

#### 1600 LD 800 MD PTA -------------------------------------------

# Use read_nonmem_table function to import data
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/'
directory_BW_CRRT_05 <- paste0(working.directory, '/sim_BW_CRRT_05.1.npctab.dta') 
BW_CRRT_05 <- read_nonmem_table(directory_BW_CRRT_05)

# Extract only the trough events
BW_CRRT_dos_05 <- BW_CRRT_05[BW_CRRT_05$TROUGH == 1,]

# Convert ID into numeric
BW_CRRT_dos_05$ID <- as.integer(BW_CRRT_dos_05$ID)

# Create a new variable AUC24 for each ID
BW_CRRT_dos_05 <- BW_CRRT_dos_05 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

BW_CRRT_dos_05 <- BW_CRRT_dos_05 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
BW_CRRT_dos_05$fAUC <- BW_CRRT_dos_05$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and BW columns
PTA_CRRT_05 <- subset(BW_CRRT_dos_05, select = c("ID", "DV", "fAUC", "DAY", "BW"))

# Calculate percentage of DV >= 7.5 per DAY per BW (PTA_Cmin_75)
PTA_CRRT_05$PTA_Cmin_75 <- 100 * ave(PTA_CRRT_05$DV >= 7.5, PTA_CRRT_05$DAY, PTA_CRRT_05$BW, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per BW (PTA_Cmin_80)
PTA_CRRT_05$PTA_Cmin_80 <- 100 * ave(PTA_CRRT_05$DV >= 80, PTA_CRRT_05$DAY, PTA_CRRT_05$BW, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per BW (PTA_fAUC_200)
PTA_CRRT_05$PTA_fAUC_200 <- 100 * ave(PTA_CRRT_05$fAUC >= 200, PTA_CRRT_05$DAY, PTA_CRRT_05$BW, FUN = mean)

# Calculate percentage of fAUC >= 400 per DAY per BW (PTA_fAUC_400) #Isabel'comment on simulating for MIC of 4 mg/L
PTA_CRRT_05$PTA_fAUC_400 <- 100 * ave(PTA_CRRT_05$fAUC >= 400, PTA_CRRT_05$DAY, PTA_CRRT_05$BW, FUN = mean)

# Create PTA_CRRT_05_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and BW
PTA_CRRT_05_overall <- unique(PTA_CRRT_05[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200","PTA_fAUC_400", "DAY", "BW")])

# Export PTA_CRRT_dos_05 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/")
write.csv(PTA_CRRT_05_overall, "PTA_CRRT_dos_05.csv", quote = F, row.names = FALSE)

# Export BW_CRRT_dos_05.csv for future use
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/")
write.csv(BW_CRRT_dos_05, "BW_CRRT_dos_05.csv", quote = F, row.names = FALSE)

#### 1800 LD 800 MD PTA -------------------------------------------

# Use read_nonmem_table function to import data
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/'
directory_BW_CRRT_06 <- paste0(working.directory, '/sim_BW_CRRT_06.1.npctab.dta') 
BW_CRRT_06 <- read_nonmem_table(directory_BW_CRRT_06)

# Extract only the trough events
BW_CRRT_dos_06 <- BW_CRRT_06[BW_CRRT_06$TROUGH == 1,]

# Convert ID into numeric
BW_CRRT_dos_06$ID <- as.integer(BW_CRRT_dos_06$ID)

# Create a new variable AUC24 for each ID
BW_CRRT_dos_06 <- BW_CRRT_dos_06 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

BW_CRRT_dos_06 <- BW_CRRT_dos_06 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
BW_CRRT_dos_06$fAUC <- BW_CRRT_dos_06$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and BW columns
PTA_CRRT_06 <- subset(BW_CRRT_dos_06, select = c("ID", "DV", "fAUC", "DAY", "BW"))

# Calculate percentage of DV >= 7.5 per DAY per BW (PTA_Cmin_75)
PTA_CRRT_06$PTA_Cmin_75 <- 100 * ave(PTA_CRRT_06$DV >= 7.5, PTA_CRRT_06$DAY, PTA_CRRT_06$BW, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per BW (PTA_Cmin_80)
PTA_CRRT_06$PTA_Cmin_80 <- 100 * ave(PTA_CRRT_06$DV >= 80, PTA_CRRT_06$DAY, PTA_CRRT_06$BW, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per BW (PTA_fAUC_200)
PTA_CRRT_06$PTA_fAUC_200 <- 100 * ave(PTA_CRRT_06$fAUC >= 200, PTA_CRRT_06$DAY, PTA_CRRT_06$BW, FUN = mean)

# Calculate percentage of fAUC >= 400 per DAY per BW (PTA_fAUC_400) #Isabel'comment on simulating for MIC of 4 mg/L
PTA_CRRT_06$PTA_fAUC_400 <- 100 * ave(PTA_CRRT_06$fAUC >= 400, PTA_CRRT_06$DAY, PTA_CRRT_06$BW, FUN = mean)

# Create PTA_CRRT_06_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and BW
PTA_CRRT_06_overall <- unique(PTA_CRRT_06[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200","PTA_fAUC_400", "DAY", "BW")])

# Export PTA_CRRT_dos_06 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/")
write.csv(PTA_CRRT_06_overall, "PTA_CRRT_dos_06.csv", quote = F, row.names = FALSE)

# Export BW_CRRT_dos_06.csv for future use
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/")
write.csv(BW_CRRT_dos_06, "BW_CRRT_dos_06.csv", quote = F, row.names = FALSE)

#### 12 mg/kg LD 6 mg/kg MD PTA -------------------------------------------

# Use read_nonmem_table function to import data
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/'
directory_BW_CRRT_07 <- paste0(working.directory, '/sim_BW_CRRT_07.1.npctab.dta') 
BW_CRRT_07 <- read_nonmem_table(directory_BW_CRRT_07)

# Extract only the trough events
BW_CRRT_dos_07 <- BW_CRRT_07[BW_CRRT_07$TROUGH == 1,]

# Convert ID into numeric
BW_CRRT_dos_07$ID <- as.integer(BW_CRRT_dos_07$ID)

# Create a new variable AUC24 for each ID
BW_CRRT_dos_07 <- BW_CRRT_dos_07 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

BW_CRRT_dos_07 <- BW_CRRT_dos_07 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
BW_CRRT_dos_07$fAUC <- BW_CRRT_dos_07$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and BW columns
PTA_CRRT_07 <- subset(BW_CRRT_dos_07, select = c("ID", "DV", "fAUC", "DAY", "BW"))

# Calculate percentage of DV >= 7.5 per DAY per BW (PTA_Cmin_75)
PTA_CRRT_07$PTA_Cmin_75 <- 100 * ave(PTA_CRRT_07$DV >= 7.5, PTA_CRRT_07$DAY, PTA_CRRT_07$BW, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per BW (PTA_Cmin_80)
PTA_CRRT_07$PTA_Cmin_80 <- 100 * ave(PTA_CRRT_07$DV >= 80, PTA_CRRT_07$DAY, PTA_CRRT_07$BW, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per BW (PTA_fAUC_200)
PTA_CRRT_07$PTA_fAUC_200 <- 100 * ave(PTA_CRRT_07$fAUC >= 200, PTA_CRRT_07$DAY, PTA_CRRT_07$BW, FUN = mean)

# Calculate percentage of fAUC >= 400 per DAY per BW (PTA_fAUC_400) #Isabel'comment on simulating for MIC of 4 mg/L
PTA_CRRT_07$PTA_fAUC_400 <- 100 * ave(PTA_CRRT_07$fAUC >= 400, PTA_CRRT_07$DAY, PTA_CRRT_07$BW, FUN = mean)

# Create PTA_CRRT_07_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and BW
PTA_CRRT_07_overall <- unique(PTA_CRRT_07[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200","PTA_fAUC_400", "DAY", "BW")])

# Export PTA_CRRT_dos_07 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/")
write.csv(PTA_CRRT_07_overall, "PTA_CRRT_dos_07.csv", quote = F, row.names = FALSE)

# Export BW_CRRT_dos_07.csv for future use
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/")
write.csv(BW_CRRT_dos_07, "BW_CRRT_dos_07.csv", quote = F, row.names = FALSE)

#### 800 mg LD 600 mg MD PTA -------------------------------------------

# Use read_nonmem_table function to import data
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/'
directory_BW_CRRT_08 <- paste0(working.directory, '/sim_BW_CRRT_08.1.npctab.dta') 
BW_CRRT_08 <- read_nonmem_table(directory_BW_CRRT_08)

# Extract only the trough events
BW_CRRT_dos_08 <- BW_CRRT_08[BW_CRRT_08$TROUGH == 1,]

# Convert ID into numeric
BW_CRRT_dos_08$ID <- as.integer(BW_CRRT_dos_08$ID)

# Create a new variable AUC24 for each ID
BW_CRRT_dos_08 <- BW_CRRT_dos_08 |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

BW_CRRT_dos_08 <- BW_CRRT_dos_08 |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
BW_CRRT_dos_08$fAUC <- BW_CRRT_dos_08$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, and BW columns
PTA_CRRT_08 <- subset(BW_CRRT_dos_08, select = c("ID", "DV", "fAUC", "DAY", "BW"))

# Calculate percentage of DV >= 7.5 per DAY per BW (PTA_Cmin_75)
PTA_CRRT_08$PTA_Cmin_75 <- 100 * ave(PTA_CRRT_08$DV >= 7.5, PTA_CRRT_08$DAY, PTA_CRRT_08$BW, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per BW (PTA_Cmin_80)
PTA_CRRT_08$PTA_Cmin_80 <- 100 * ave(PTA_CRRT_08$DV >= 80, PTA_CRRT_08$DAY, PTA_CRRT_08$BW, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per BW (PTA_fAUC_200)
PTA_CRRT_08$PTA_fAUC_200 <- 100 * ave(PTA_CRRT_08$fAUC >= 200, PTA_CRRT_08$DAY, PTA_CRRT_08$BW, FUN = mean)

# Calculate percentage of fAUC >= 400 per DAY per BW (PTA_fAUC_400) #Isabel'comment on simulating for MIC of 4 mg/L
PTA_CRRT_08$PTA_fAUC_400 <- 100 * ave(PTA_CRRT_08$fAUC >= 400, PTA_CRRT_08$DAY, PTA_CRRT_08$BW, FUN = mean)

# Create PTA_CRRT_08_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, and BW
PTA_CRRT_08_overall <- unique(PTA_CRRT_08[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200","PTA_fAUC_400", "DAY", "BW")])

# Export PTA_CRRT_dos_08 dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/")
write.csv(PTA_CRRT_08_overall, "PTA_CRRT_dos_08.csv", quote = F, row.names = FALSE)

# Export BW_CRRT_dos_08.csv for future use
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/BW_CRRT/")
write.csv(BW_CRRT_dos_08, "BW_CRRT_dos_08.csv", quote = F, row.names = FALSE)

### Making PTA plots -------------------------------------------

#### Some preps -------------------------------------------

# Import PTA_dos_ datasets
PTA_CRRT_dos_01 <- read.csv("PTA_BW_CRRT/PTA_CRRT_dos_01.csv")
PTA_CRRT_dos_02 <- read.csv("PTA_BW_CRRT/PTA_CRRT_dos_02.csv")
PTA_CRRT_dos_03 <- read.csv("PTA_BW_CRRT/PTA_CRRT_dos_03.csv")
PTA_CRRT_dos_04 <- read.csv("PTA_BW_CRRT/PTA_CRRT_dos_04.csv")
PTA_CRRT_dos_05 <- read.csv("PTA_BW_CRRT/PTA_CRRT_dos_05.csv")
PTA_CRRT_dos_06 <- read.csv("PTA_BW_CRRT/PTA_CRRT_dos_06.csv")
PTA_CRRT_dos_07 <- read.csv("PTA_BW_CRRT/PTA_CRRT_dos_07.csv")
PTA_CRRT_dos_08 <- read.csv("PTA_BW_CRRT/PTA_CRRT_dos_08.csv")

# Add regimens
PTA_CRRT_dos_01$Regimen <- 1
PTA_CRRT_dos_02$Regimen <- 2
PTA_CRRT_dos_03$Regimen <- 3
PTA_CRRT_dos_04$Regimen <- 4
PTA_CRRT_dos_05$Regimen <- 5
PTA_CRRT_dos_06$Regimen <- 6
PTA_CRRT_dos_07$Regimen <- 7
PTA_CRRT_dos_08$Regimen <- 8

# Combine 8 PTA datasets into 1 - PTA_BW_CRRT
PTA_BW_CRRT <- rbind(PTA_CRRT_dos_01, PTA_CRRT_dos_02, PTA_CRRT_dos_03, 
                       PTA_CRRT_dos_04, PTA_CRRT_dos_05, PTA_CRRT_dos_06,
                       PTA_CRRT_dos_07, PTA_CRRT_dos_08)

# Define the labels and colors
dosing_labels_day12 <- c("800  mg  q24h LD",
                         "1000 mg  q24h LD",
                         "1200 mg  q24h LD",
                         "1400 mg  q24h LD",
                         "1600 mg  q24h LD",
                         "1800 mg  q24h LD",
                         "12 mg/kg q24h LD")
# Regimens 1, 2, 3, 4, 5, 6, 7, respectively

dosing_labels_day714 <- c("400 mg q24h MD", 
                          "600 mg q24h MD",
                          "800 mg q24h MD",
                          "6 mg/kg q24h MD")
# Regimens 1, 8, 3, 7, respectively

dosing_colors_day12 <- c("#addc30",
                         "#5ec962",
                         "#28ae80",
                         "#21918c",
                         "#2c728e",
                         "#3b528b",
                         "#472d7b")

dosing_colors_day714 <- c("#addc30",
                          "#440154",
                          "#28ae80",
                          "#3b528b")

# Reorder the levels of the dosing regimens
PTA_BW_CRRT$Regimen <- factor(PTA_BW_CRRT$Regimen,
                                levels = c(1, 8, 2, 3, 4, 5, 6, 7))


# I also enforce limits on vertical lines here

# Firstly, create dataset with BW, Regimen and PTA on day 1 
BW_Reg <- data.frame( #BW and Regimen dataset
  BW = c(60, 80, 100, 120, 140),
  Regimen = c(2, 3, 4, 5, 6)
)

BW_Reg$Regimen <- factor(BW_Reg$Regimen,
                         levels = c(2, 3, 4, 5, 6))

# Then, make a dataset of PTA on Day 1 per BW per Regimen
BW_PTA_DaY1 <- BW_Reg %>%
  left_join(
    PTA_BW_CRRT %>% 
      dplyr::filter(DAY == 1),
    by = c("BW", "Regimen")
  )

# Similarly for day 2
BW_PTA_DaY2 <- BW_Reg %>%
  left_join(
    PTA_BW_CRRT %>% 
      dplyr::filter(DAY == 2),
    by = c("BW", "Regimen")
  )

#### Day 1 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_BW_CRRT200_DAY1 <-  ggplot(PTA_BW_CRRT[PTA_BW_CRRT$DAY == 1 & !PTA_BW_CRRT$Regimen == 8, ], 
                                 aes(x = BW, 
                                     y = PTA_fAUC_200, 
                                     group = factor(Regimen), 
                                     color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(30, 150), breaks = seq(30, 150, by = 5), expand = c(0,0)) +
  xlab(TeX(r"(Total body weight (kg))")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) +
  ylab(TeX(r"($Probability\ of\ \textit{f}AUC_{0-24}\ \geq 200 \, \mg \times h/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 1") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day12, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day12) +
  geom_hline(yintercept = 90, 
             linetype = "dashed", 
             color = "gray50", 
             size = 0.6) +
  geom_segment(
    data = BW_PTA_DaY1,
    aes(
      x = BW,
      xend = BW,
      y = 0,  # Starting from the x-axis
      yend = PTA_fAUC_200  # Ending at the corresponding PTA_fAUC_200 value
    ),
    linetype = "dashed",
    color = "gray50",
    size = 0.6
  )
PTA_BW_CRRT200_DAY1

#### Day 14 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_BW_CRRT200_DAY14 <-  ggplot(PTA_BW_CRRT[PTA_BW_CRRT$DAY == 14 & PTA_BW_CRRT$Regimen %in% c(1, 8, 3, 7), ], 
                                  aes(x = BW, 
                                      y = PTA_fAUC_200, 
                                      group = factor(Regimen), 
                                      color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(30, 150), breaks = seq(30, 150, by = 5), expand = c(0,0)) +
  xlab(TeX(r"(Total body weight (kg))")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) +
  ylab(TeX(r"($Probability\ of\ \textit{f}AUC_{0-24}\ \geq 200 \, \mg \times h/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 14") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day714, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day714) +
  geom_hline(yintercept = 90, 
             linetype = "dashed", 
             color = "gray50", 
             size = 0.6) 
PTA_BW_CRRT200_DAY14

#### Combining these 4 plots into 1 -------------------------------------------

## Combining DAY 1 & 14 into 1

# Combine the two plots
PTA_BW_CRRT200_DAY114 <- ggarrange(PTA_BW_CRRT200_DAY1, 
                                     PTA_BW_CRRT200_DAY14,
                                     labels = c("a", "b"),
                                     ncol = 1, 
                                     nrow = 2,
                                     font.label = list(size = 8, face = "bold"))

# Export the combined plots
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Plots/Dose_finding_simulations/BW_CRRT")

# SVG
ggsave("PTA_BW_CRRT200_DAY114.svg", 
       PTA_BW_CRRT200_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

# JPEG
ggsave("PTA_BW_CRRT200_DAY114.JPEG", 
       PTA_BW_CRRT200_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

# EPS
ggsave("PTA_BW_CRRT200_DAY114.EPS", 
       PTA_BW_CRRT200_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

### Assessing toxicity PTA Cmin80 -------------------------------------------

# PTA_CKDEPI dataset was created above

#### Day 1 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_BW_CRRT80_DAY1 <-  ggplot(PTA_BW_CRRT[PTA_BW_CRRT$DAY == 1 & !PTA_BW_CRRT$Regimen == 8, ], 
                                aes(x = BW, 
                                    y = PTA_Cmin_80, 
                                    group = factor(Regimen), 
                                    color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(30, 150), breaks = seq(30, 150, by = 5), expand = c(0,0)) +
  xlab(TeX(r"(Total body weight (kg))")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) +
  ylab(TeX(r"($Probability\ of\ \C_{min}\ \geq 80 \, \mg/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 1") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day12, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day12) +
  geom_segment(
    data = BW_PTA_DaY1,
    aes(
      x = BW,
      xend = BW,
      y = 0,  # Starting from the x-axis
      yend = PTA_Cmin_80  # Ending at the corresponding PTA_Cmin_80 value
    ),
    linetype = "dashed",
    color = "gray50",
    size = 0.6
  )
PTA_BW_CRRT80_DAY1

#### Day 14 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_BW_CRRT80_DAY14 <-  ggplot(PTA_BW_CRRT[PTA_BW_CRRT$DAY == 14 & PTA_BW_CRRT$Regimen %in% c(1, 8, 3, 7), ], 
                                 aes(x = BW, 
                                     y = PTA_Cmin_80, 
                                     group = factor(Regimen), 
                                     color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(30, 150), breaks = seq(30, 150, by = 5), expand = c(0,0)) +
  xlab(TeX(r"(Total body weight (kg))")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) +
  ylab(TeX(r"($Probability\ of\ \C_{min}\ \geq 80 \, \mg/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 14") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day714, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day714) 
PTA_BW_CRRT80_DAY14

#### Combining these 4 plots into 1 -------------------------------------------

## Combining DAY 1 & 14 into 1

# Combine the two plots
PTA_BW_CRRT80_DAY114 <- ggarrange(PTA_BW_CRRT80_DAY1,
                                  PTA_BW_CRRT80_DAY14,
                                  labels = c("a", "b"),
                                  ncol = 1,
                                  nrow = 2,
                                  font.label = list(size = 8, face = "bold"))

# Export the combined plots
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Plots/Dose_finding_simulations/BW_CRRT")

# SVG
ggsave("PTA_BW_CRRT80_DAY114.svg", 
       PTA_BW_CRRT80_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

# JPEG
ggsave("PTA_BW_CRRT80_DAY114.JPEG", 
       PTA_BW_CRRT80_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

# EPS
ggsave("PTA_BW_CRRT80_DAY114.EPS", 
       PTA_BW_CRRT80_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

## Fourth, CKDEPI optimal dose -------------------------------------------

### Creating datasets for PTA plots -------------------------------------------

# Use read_nonmem_table function to import data
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/CKDEPI_Opt_dose/'
directory_CKDEPI_Opt <- paste0(working.directory, '/BW_stra_LD_CKDEPI.1.npctab.dta') 
CKDEPI_Opt <- read_nonmem_table(directory_CKDEPI_Opt)

# Extract only the trough events
CKDEPI_Opt_dos <- CKDEPI_Opt[CKDEPI_Opt$TROUGH == 1,]

# Convert ID into numeric
CKDEPI_Opt_dos$ID <- as.integer(CKDEPI_Opt_dos$ID)

# Create a new variable AUC24 for each ID
CKDEPI_Opt_dos <- CKDEPI_Opt_dos |> 
  group_by(ID) |> 
  mutate(AUC24 = AUC2 - dplyr::lag(AUC2, default = AUC2[1])) |> 
  ungroup()

CKDEPI_Opt_dos <- CKDEPI_Opt_dos |> 
  mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
CKDEPI_Opt_dos$fAUC <- CKDEPI_Opt_dos$AUC24 * 0.89

# Extract ID, DV, fAUC, DAY, CKDEPI, and Regimen columns
PTA_CKDEPI_Opt <- subset(CKDEPI_Opt_dos, select = c("ID", "DV", "fAUC", "DAY", "CKDEPI", "Regimen"))

# Calculate percentage of DV >= 7.5 per DAY per CKDEPI per Regimen (PTA_Cmin_75)
PTA_CKDEPI_Opt$PTA_Cmin_75 <- 100 * ave(PTA_CKDEPI_Opt$DV >= 7.5, PTA_CKDEPI_Opt$DAY, PTA_CKDEPI_Opt$CKDEPI, PTA_CKDEPI_Opt$Regimen, FUN = mean)

# Calculate percentage of DV >= 80 per DAY per BW (PTA_Cmin_80)
PTA_CKDEPI_Opt$PTA_Cmin_80 <- 100 * ave(PTA_CKDEPI_Opt$DV >= 80, PTA_CKDEPI_Opt$DAY, PTA_CKDEPI_Opt$CKDEPI, PTA_CKDEPI_Opt$Regimen, FUN = mean)

# Calculate percentage of fAUC >= 200 per DAY per BW (PTA_fAUC_200)
PTA_CKDEPI_Opt$PTA_fAUC_200 <- 100 * ave(PTA_CKDEPI_Opt$fAUC >= 200, PTA_CKDEPI_Opt$DAY, PTA_CKDEPI_Opt$CKDEPI, PTA_CKDEPI_Opt$Regimen, FUN = mean)

# Calculate percentage of fAUC >= 400 per DAY per BW (PTA_fAUC_400) #Isabel'comment on simulating for MIC of 4 mg/L
PTA_CKDEPI_Opt$PTA_fAUC_400 <- 100 * ave(PTA_CKDEPI_Opt$fAUC >= 400, PTA_CKDEPI_Opt$DAY, PTA_CKDEPI_Opt$CKDEPI, PTA_CKDEPI_Opt$Regimen, FUN = mean)

# Create PTA_CKDEPI_Opt_overall dataset with unique values of PTA_Cmin, PTA_fAUC, DAY, CKDEPI, and Regimen
PTA_CKDEPI_Opt_overall <- unique(PTA_CKDEPI_Opt[c("PTA_Cmin_75","PTA_Cmin_80", "PTA_fAUC_200","PTA_fAUC_400", "DAY", "CKDEPI", "Regimen")])

# Export PTA_CKDEPI_Opt_dos dataset
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/CKDEPI_Opt_dose/")
write.csv(PTA_CKDEPI_Opt_overall, "PTA_CKDEPI_Opt_dos.csv", quote = F, row.names = FALSE)

# Export CKDEPI_Opt_dos.csv for future use
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/CKDEPI_Opt_dose/")
write.csv(CKDEPI_Opt_dos, "CKDEPI_Opt_dos.csv", quote = F, row.names = FALSE)

### Making PTA plots -------------------------------------------

#### Some preps -------------------------------------------

# Import PTA_dos_ datasets
PTA_CKDEPI_Opt <- read.csv("PTA_CKDEPI_Opt/PTA_CKDEPI_Opt_dos.csv")

# Define the labels and colors
dosing_labels_day1 <-  c("800  mg  q24h LD",
                         "1000 mg  q24h LD",
                         "1200 mg  q24h LD",
                         "1400 mg  q24h LD",
                         "1600 mg  q24h LD")

dosing_labels_day27 <- c("800  mg  q24h LD / 400 mg q24h MD",
                         "1000 mg  q24h LD / 400 mg q24h MD",
                         "1200 mg  q24h LD / 400 mg q24h MD",
                         "1400 mg  q24h LD / 400 mg q24h MD",
                         "1600 mg  q24h LD / 400 mg q24h MD")

dosing_labels_day14 <- c("400 mg q24h MD")

dosing_colors_day127 <- c("#7ad151",
                          "#22a884",
                          "#2a788e",
                          "#414487",
                          "#440154")

dosing_colors_day14  <- c("#440154")

#### Day 1 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_CKDEPI_Opt200_DAY1 <-  ggplot(PTA_CKDEPI_Opt[PTA_CKDEPI_Opt$DAY == 1, ], 
                               aes(x = CKDEPI, 
                                   y = PTA_fAUC_200, 
                                   group = factor(Regimen), 
                                   color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(5, 215), breaks = seq(5, 215, by = 10), expand = c(0,0)) +
  xlab(TeX(r"($\eGFR_{CKD-EPI}\ (ml/min/1.73m^2)$)")) + 
  scale_y_continuous(limits = c(80, 100), breaks = seq(80, 100, by = 5), expand = c(0,0)) + 
  # The y-axis ranging from 80 to 100 to make the distinction clearer
  ylab(TeX(r"($Probability\ of\ \textit{f}AUC_{0-24}\ \geq 200 \, \mg \times h/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 1") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day127, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day1) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "gray50", size = 0.6) 
PTA_CKDEPI_Opt200_DAY1

# Tomorrow work from here

#### Day 14 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_CKDEPI_Opt200_DAY14 <-  ggplot(PTA_CKDEPI_Opt[PTA_CKDEPI_Opt$DAY == 14 & PTA_CKDEPI_Opt$Regimen %in% c(1), ], 
                                aes(x = CKDEPI, 
                                    y = PTA_fAUC_200, 
                                    group = factor(Regimen), 
                                    color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(5, 215), breaks = seq(5, 215, by = 10), expand = c(0,0)) +
  xlab(TeX(r"($\eGFR_{CKD-EPI}\ (ml/min/1.73m^2)$)")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) + 
  # The y-axis ranging from 80 to 100 to make the distinction clearer
  ylab(TeX(r"($Probability\ of\ \textit{f}AUC_{0-24}\ \geq 200 \, \mg \times h/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 14") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day14, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day14) +
  geom_hline(yintercept = 90, 
             linetype = "dashed", 
             color = "gray50", 
             size = 0.6) 
PTA_CKDEPI_Opt200_DAY14

#### Combining these 4 plots into 1 -------------------------------------------

## Combining DAY 1 & 14 into 1

# Combine the two plots
PTA_CKDEPI_Opt200_DAY114 <- ggarrange(PTA_CKDEPI_Opt200_DAY1, 
                                   PTA_CKDEPI_Opt200_DAY14,
                                   labels = c("a", "b"),
                                   ncol = 1, 
                                   nrow = 2,
                                   font.label = list(size = 8, face = "bold"))

# Export the combined plots
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Plots/Dose_finding_simulations/CKDEPI_Opt")

# SVG
ggsave("PTA_CKDEPI_Opt200_DAY114.svg", 
       PTA_CKDEPI_Opt200_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

# JPEG
ggsave("PTA_CKDEPI_Opt200_DAY114.JPEG", 
       PTA_CKDEPI_Opt200_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

# EPS
ggsave("PTA_CKDEPI_Opt200_DAY114.EPS", 
       PTA_CKDEPI_Opt200_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

### Assessing toxicity PTA Cmin80 -------------------------------------------

# PTA_CKDEPI dataset was created above

#### Day 1 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_CKDEPI_Opt80_DAY1 <-  ggplot(PTA_CKDEPI_Opt[PTA_CKDEPI_Opt$DAY == 1, ], 
                              aes(x = CKDEPI, 
                                  y = PTA_Cmin_80, 
                                  group = factor(Regimen), 
                                  color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(5, 215), breaks = seq(5, 215, by = 10), expand = c(0,0)) +
  xlab(TeX(r"($\eGFR_{CKD-EPI}\ (ml/min/1.73m^2)$)")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) + 
  ylab(TeX(r"($Probability\ of\ \C_{min}\ \geq 80 \, \mg/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 1") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day127, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day1) 
PTA_CKDEPI_Opt80_DAY1

#### Day 14 -------------------------------------------

# Here we use Colorblind-Friendly Palette Viridis
PTA_CKDEPI_Opt80_DAY14 <-  ggplot(PTA_CKDEPI_Opt[PTA_CKDEPI_Opt$DAY == 14 & PTA_CKDEPI_Opt$Regimen %in% c(1), ], 
                               aes(x = CKDEPI, 
                                   y = PTA_Cmin_80, 
                                   group = factor(Regimen), 
                                   color = factor(Regimen))) +
  geom_line(size = 1) +
  scale_x_continuous(limits = c(5, 215), breaks = seq(5, 215, by = 10), expand = c(0,0)) +
  xlab(TeX(r"($\eGFR_{CKD-EPI}\ (ml/min/1.73m^2)$)")) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0,0)) + 
  ylab(TeX(r"($Probability\ of\ \C_{min}\ \geq 80 \, \mg/L\ (\%)$)")) +
  theme_minimal() +
  ggtitle("Day 14") +
  theme(legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 6, family = "Helvetica"),
        axis.text = element_text(size = 5.5, family = "Helvetica"),
        legend.title = element_text(size = 5, family = "Helvetica"),
        legend.text = element_text(size = 5, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + # remove all the clutters, i.e., grid minor & grid major
  labs(color = "Regimen") +
  scale_color_manual(values = dosing_colors_day14, 
                     name = "Dosing regimen",
                     labels = dosing_labels_day14) +
  geom_hline(yintercept = 90, 
             linetype = "dashed", 
             color = "gray50", 
             size = 0.6)
PTA_CKDEPI_Opt80_DAY14

#### Combining these 4 plots into 1 -------------------------------------------

## Combining DAY 1 & 14 into 1

# Combine the two plots
PTA_BW_CRRT80_DAY114 <- ggarrange(PTA_BW_CRRT80_DAY1,
                                  PTA_BW_CRRT80_DAY14,
                                  labels = c("a", "b"),
                                  ncol = 1,
                                  nrow = 2,
                                  font.label = list(size = 8, face = "bold"))

# Export the combined plots
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Plots/Dose_finding_simulations/BW_CRRT")

# SVG
ggsave("PTA_BW_CRRT80_DAY114.svg", 
       PTA_BW_CRRT80_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

# JPEG
ggsave("PTA_BW_CRRT80_DAY114.JPEG", 
       PTA_BW_CRRT80_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

# EPS
ggsave("PTA_BW_CRRT80_DAY114.EPS", 
       PTA_BW_CRRT80_DAY114, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

# Pop_dosing analysis ----------------------------------------------------------

## Conc Time plot -----------------------------------------

### Standard dosing -------------------------------------------

# Load the dataset
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Pop_sim/'
std_tablefile <- paste0(working.directory, '/pop_std_conc.1.npctab.dta') 
pop_dosing_std <- read_nonmem_table(std_tablefile)

# Extract only TIME and DV columns
DV_TIME_std <- pop_dosing_std %>% select(DV, TIME)
setwd('/lustre1/scratch/357/vsc35700/Fluco_revised/Pop_sim/')
# write.csv(DV_TIME_std, "DV_TIME_std.csv",quote = F,row.names = FALSE)

# Remove the rows where DV is 0 and TIME is not 0
DV_TIME_std <- DV_TIME_std[DV_TIME_std$DV != 0 | DV_TIME_std$TIME == 0, ]

# Calculate 5th and 95th percentiles
data_summary_std <- DV_TIME_std %>%
  group_by(TIME) %>%
  summarize(median = median(DV),
            min = quantile(DV, 0.05),
            max = quantile(DV, 0.95),
            .groups = 'drop')

# Export data to the project folder
setwd('/data/leuven/357/vsc35700/Fluco_revised/R_code/Pop_sim/')
write.csv(data_summary_std, "data_summary_std.csv", quote = F, row.names = F)

# Round the values to 2 decimal places
data_summary_std <- data_summary_std %>%
  mutate(across(everything(), ~ round(., 2)))

# ribbon plot std dosing
pop_sim_std <- ggplot(data = data_summary_std, aes(x = TIME/24, y = median)) +
  geom_ribbon(aes(ymin = min, ymax = max), fill = "#21918c", alpha = 0.5) +
  geom_line(color = "#440154", size = 0.8) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "gray50",size = 0.6) +
  xlab("Day") +
  ylab("Fluconazole total concentration (mg/L)") +
  scale_x_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 1),expand = c(0,0.01)) +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 10),expand = c(0,0.01)) +
  ggtitle("Standard dosing regimen") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), 
        axis.title = element_text(size = 7, family = "Helvetica"),
        axis.text = element_text(size = 7, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pop_sim_std

### Optimised dosing -------------------------------------------

# Load the dataset
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Pop_sim/'
opt_tablefile <- paste0(working.directory, '/pop_opt_conc.1.npctab.dta') 
pop_dosing_opt <- read_nonmem_table(opt_tablefile)

# Extract only TIME and DV columns
DV_TIME_opt <- pop_dosing_opt %>% select(DV, TIME)
setwd('/lustre1/scratch/357/vsc35700/Fluco_revised/Pop_sim/')
# write.csv(DV_TIME_opt, "DV_TIME_opt.csv",quote = F,row.names = FALSE)

# Remove the rows where DV is 0 and TIME is not 0
DV_TIME_opt <- DV_TIME_opt[DV_TIME_opt$DV != 0 | DV_TIME_opt$TIME == 0, ]

# Calculate 5th and 95th percentiles
data_summary_opt <- DV_TIME_opt %>%
  group_by(TIME) %>%
  summarize(median = median(DV),
            min = quantile(DV, 0.05),
            max = quantile(DV, 0.95),
            .groups = 'drop')

# Round the values to 2 decimal places
data_summary_opt <- data_summary_opt %>%
  mutate(across(everything(), ~ round(., 2)))

# Export data to the project folder
setwd('/data/leuven/357/vsc35700/Fluco_revised/R_code/Pop_sim/')
write.csv(data_summary_opt, "data_summary_opt.csv", quote = F, row.names = F)

# ribbon plot opt dosing
pop_sim_opt <- ggplot(data = data_summary_opt, aes(x = TIME/24, y = median)) +
  geom_ribbon(aes(ymin = min, ymax = max), fill = "#21918c", alpha = 0.5) +
  geom_line(color = "#440154", size = 0.8) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "gray50",size = 0.6) +
  xlab("Day") +
  ylab("Fluconazole total concentration (mg/L)") +
  scale_x_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 1),expand = c(0,0.01)) +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 10),expand = c(0,0.01)) +
  ggtitle("Optimised dosing regimen") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), 
        axis.title = element_text(size = 7, family = "Helvetica"),
        axis.text = element_text(size = 7, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pop_sim_opt

### Save the plot -----------------------------------------------

# Here I use plot_grid function as it works for the Rstudio cluster
# Just note that ggarrange function works just fine too
conc_time_std_opt <- plot_grid(pop_sim_std, 
                               pop_sim_opt,
                               labels = c("a", "b"),
                               align = "h",
                               label_fontfamily = "Helvetica",
                               label_fontface = "bold",
                               label_size = 8)

# Save the plot
setwd('/data/leuven/357/vsc35700/Fluco_revised/R_code/Plots/Pop_sim/')
# SVG
ggsave("conc_time_std_opt.svg", 
       conc_time_std_opt, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")
# JPEG
ggsave("conc_time_std_opt.JPEG", 
       conc_time_std_opt, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")
# EPS
ggsave("conc_time_std_opt.EPS", 
       conc_time_std_opt, 
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

## Box plot & ROC -----------------------------------------

### ROC cuvrve -------------------------------------------

#### Generate datasets --------------------------------------

## Standard dosing

working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Pop_sim/'
pop_box_std_dir <- paste0(working.directory, '/pop_sim_std.1.npctab.dta') 
pop_box_std <- read_nonmem_table(pop_box_std_dir)

# Extract only the trough events
pop_std_dos <- pop_box_std[pop_box_std$TROUGH == 1,]

# Convert ID into numeric
pop_std_dos$ID <- as.integer(pop_std_dos$ID)

# Create a new variable AUC24 for each ID
pop_std_dos <- pop_std_dos %>%
  group_by(ID) %>%
  mutate(AUC24 = AUC2 - lag(AUC2, default = AUC2[1])) %>%
  ungroup()

pop_std_dos <- pop_std_dos %>% mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
pop_std_dos$fAUC <- pop_std_dos$AUC24*0.89

## Optimised dosing

pop_box_opt_dir <- paste0(working.directory, '/pop_sim_opt.1.npctab.dta') 
pop_box_opt <- read_nonmem_table(pop_box_opt_dir)

# Extract only the trough events
pop_opt_dos <- pop_box_opt[pop_box_opt$TROUGH == 1,]

# Convert ID into numeric
pop_opt_dos$ID <- as.integer(pop_opt_dos$ID)

# Create a new variable AUC24 for each ID
pop_opt_dos <- pop_opt_dos %>%
  group_by(ID) %>%
  mutate(AUC24 = AUC2 - lag(AUC2, default = AUC2[1])) %>%
  ungroup()

pop_opt_dos <- pop_opt_dos %>% mutate(AUC24 = ifelse(DAY == "1", AUC2,AUC24))

# Create a new variable fAUC, which equals AUC24*89 (we assume protein binding is 11%)
pop_opt_dos$fAUC <- pop_opt_dos$AUC24*0.89

#### BW's impact --------------------------------------------

TA_std_day1 <- pop_std_dos %>%
  filter(DAY == 1) %>%
  mutate(TA = ifelse(fAUC >= 200, 1, 0))

set.seed(123)
ROC_BW_std_day1 <- pROC::roc(TA_std_day1$TA,
                             TA_std_day1$BW,
                             auc.polygon = FALSE, 
                             max.auc.polygon = FALSE, 
                             grid = FALSE,
                             ci = TRUE, 
                             boot.n = 1000, 
                             ci.alpha = 0.95, 
                             stratified = TRUE,
                             plot = FALSE, 
                             print.auc = TRUE, 
                             print.auc.x = 0.2, 
                             print.auc.y = 0.3,
                             show.thres = TRUE)

# Calculate 95% CI
AUC_std_day1 <- pROC::auc(ROC_BW_std_day1, conf.level = 0.95, ci.alpha = 0.95, method = "bootstrap", boot.n = 1000, stratified = TRUE)

AUC_firstpart_std_day1 <- round(as.numeric(sub(".*(: *)", "//1", AUC_std_day1)), 3)

CI_AUC_std_day1 <- pROC::ci.auc(ROC_BW_std_day1, conf.level = 0.95, ci.alpha = 0.95, method = "bootstrap", boot.n = 1000, stratified = TRUE)

AUC_anotherpart_std_day1 <- sub(".*(: *)", "//1", CI_AUC_std_day1)

AUC_secondpart_std_day1 <- as.numeric(substr(AUC_anotherpart_std_day1[1], 1, 5))

AUC_thirdpart_std_day1 <- as.numeric(substr(AUC_anotherpart_std_day1[3], 1, 5))

AUC_label_std_day1 <- paste0("AUC: ", AUC_firstpart_std_day1, " (", AUC_secondpart_std_day1, "-", AUC_thirdpart_std_day1, ")", sep = "") #needs to figure out how to save this object for next time

# For recommended dosing
TA_opt_day1 <- pop_opt_dos %>%
  filter(DAY == 1) %>% mutate(TA = ifelse(fAUC >= 200, 1, 0))

set.seed(123)
ROC_BW_opt_day1 <- pROC::roc(TA_opt_day1$TA,
                             TA_opt_day1$BW,
                             auc.polygon = FALSE,
                             max.auc.polygon = FALSE,
                             grid = FALSE,
                             ci = TRUE,
                             boot.n = 1000,
                             ci.alpha = 0.95,
                             stratified = TRUE,
                             plot = FALSE,
                             print.auc = TRUE,
                             print.auc.x = 0.2,
                             print.auc.y = 0.3,
                             show.thres = TRUE)

# Calculate 95% CI for the roc curves
AUC_opt_day1 <- pROC::auc(ROC_BW_opt_day1, conf.level = 0.95, ci.alpha = 0.95, method = "bootstrap", boot.n = 1000, stratified = TRUE)

AUC_firstpart_opt_day1 <- round(as.numeric(sub(".*(: *)", "//1", AUC_opt_day1)), 3)

CI_AUC_opt_day1 <- pROC::ci.auc(ROC_BW_opt_day1, conf.level = 0.95, ci.alpha = 0.95, method = "bootstrap", boot.n = 1000, stratified = TRUE)

AUC_anotherpart_opt_day1 <- sub(".*(: *)", "//1", CI_AUC_opt_day1)

AUC_secondpart_opt_day1 <- as.numeric(substr(AUC_anotherpart_opt_day1[1], 1, 5))

AUC_thirdpart_opt_day1 <- as.numeric(substr(AUC_anotherpart_opt_day1[3], 1, 5))

AUC_label_opt_day1 <- paste0("AUC: ", AUC_firstpart_opt_day1, " (", AUC_secondpart_opt_day1, "-", AUC_thirdpart_opt_day1, ")", sep = "")

# Convert roc objects to data frames
roc_std_day1 <- pROC::coords(ROC_BW_std_day1)
roc_opt_day1 <- pROC::coords(ROC_BW_opt_day1)

# Create the plot with DAY 1 title
roc_curve_day1 <- ggplot() +
  geom_line(data = roc_std_day1, aes(x = 1 - specificity, y = sensitivity, color = "Standard dosing regimen", linetype = "Standard dosing regimen"), size = 0.8) +
  geom_line(data = roc_opt_day1, aes(x = 1 - specificity, y = sensitivity, color = "Optimised dosing regimen", linetype = "Optimised dosing regimen"), size = 0.8) +
  scale_color_manual(values = c("Standard dosing regimen" = "#440154",
                                "Optimised dosing regimen" = "#31688e"),
                     breaks = c('Standard dosing regimen', 'Optimised dosing regimen')) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  theme_minimal() +
  ggtitle("Day 1") +
  labs(y = "Sensitivity", x = "1 - Specificity") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), 
        axis.title = element_text(size = 7, family = "Helvetica"),
        axis.text = element_text(size = 7, family = "Helvetica"),
        legend.title = element_text(size = 6, family = "Helvetica"),
        legend.text = element_text(size = 6, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "gray50", fill = NA, size = 0.8),
        axis.line = element_line(colour = "gray50")) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "#35b779", linetype = "solid") +
  annotate(geom = "text", x = 0.60, y = 0.30, size = 7/.pt, 
           label = AUC_label_std_day1, color = "#440154", 
           family = theme_get()$text[["family"]]) +
  annotate(geom = "text", x = 0.60, y = 0.20, size = 7/.pt, 
           label = AUC_label_opt_day1, color = "#31688e",
           family = theme_get()$text[["family"]]) +
  guides(color = guide_legend(override.aes = list(linetype = c("solid", "dashed"))), linetype = "none")

roc_curve_day1





