

# Call libraries ----------------------------------------------------------

library(styler)
library(tidyverse)
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
library(table1) #create descriptive statistic table1

# Explanatory Data Analysis -----------------------------------------------

## Load the data ------------------------------------------------------------

setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Original")
FlucTotIV_clean <- read.csv("FlucTotIV_clean.csv", na.strings = "NA")
View(FlucTotIV_clean)

## Clean the dataset, make IHD from non CRRT to CRRT -------------------------

# When IHD == 1, converts CRRT == 1, otherwise remains

FlucTotIV_clean$CRRT[FlucTotIV_clean$IHD == 1] <- 1

# Convert AMT == 0 back into NA

FlucTotIV_clean$AMT[FlucTotIV_clean$AMT == 0] <- NA

# Convert DV == 0 back into NA

FlucTotIV_clean$DV[FlucTotIV_clean$DV == 0] <- NA

# Convert TAD and RATE == 0 back into NA

FlucTotIV_clean$TAD[FlucTotIV_clean$TAD == 0] <- NA

FlucTotIV_clean$RATE[FlucTotIV_clean$RATE == 0] <- NA

# Export into Fluco_clean_revised

write.csv(FlucTotIV_clean, "Fluco_clean_revised.csv", quote = F, row.names = F)

## Check the CRRT pattern -------------------------

# Filter every patient that has at least one observation with CRRT == 1

CRRT_one <- Fluco_clean_revised |> 
  group_by(ID) |>
  dplyr::filter(any(CRRT == 1)) 
# 34 patients

# Check how many patients have CRRT == 1 all the time
CRRT_one |> 
  group_by(ID) |> 
  summarise(CRRT = all(CRRT == 1)) |> 
  summarise(n = sum(CRRT))
# 18 patients have CRRT == 1 all the time

# Add DAY column, which equals 1 if TIME <=24, 2 if 24 < TIME <= 48, and so on

CRRT_one <- CRRT_one |> 
  mutate(DAY = cut(TIME, breaks = seq(0, 648, 24), labels = 1:27))

# If TIME == 0, DAY == 1

CRRT_one$DAY[CRRT_one$TIME == 0] <- 1

# For those who have both CRRT == 1 and CRRT == 0, calculate the percentage of 
# DAY that CRRT == 1 over all DAY

# First select ID that has CRRT == 1 and CRRT == 0

Miscel_CRRT <- CRRT_one |>
  group_by(ID) |> 
  summarise(CRRT = n_distinct(CRRT)) |> 
  dplyr::filter(CRRT == 2) |> 
  select(ID)

# Percentage of days that CRRT == 1
CRRT_one |> 
  dplyr::filter(ID %in% Miscel_CRRT$ID) |>
  group_by(ID) |> 
  group_by(DAY) |>
  summarise(CRRT = mean(CRRT == 1, na.rm = T)) |> 
  summarise(mean = mean(CRRT, na.rm = T))

# Try another way, calculate the total number of days that CRRT == 1, then divided
# by the total number of days
CRRT_one |> 
  dplyr::filter(ID %in% Miscel_CRRT$ID) |>
  group_by(DAY) |>
  select(DAY, CRRT) |>
  summarise(total_days = n(), CRRT_days = sum(CRRT == 1, na.rm = T)) |>
  summarise(mean = sum(CRRT_days) / sum(total_days))

# Another way, extract the unique combination of ID, DAY and CRRT == 1,
# Sum up the number of days, then divided by the total number of days

# Total CRRT days
CRRT_one |> 
  dplyr::filter(ID %in% Miscel_CRRT$ID) |>
  group_by(ID, DAY) |>
  dplyr::filter(CRRT == 1) |>
  unique() |>
  ungroup() |>
  group_by(ID) |>
  select(ID, DAY) |>
  summarise(CRRT_days = n()) |>
  summarise(Total_CRRT_days = sum(CRRT_days))
# 199 CRRT days

# Total number of days
CRRT_one |> 
  dplyr::filter(ID %in% Miscel_CRRT$ID) |>
  group_by(ID, DAY) |>
  unique() |>
  ungroup() |>
  group_by(ID) |>
  select(ID, DAY) |>
  summarise(total_days = n()) |>
  summarise(Total_days = sum(total_days))
# 340 total days

## Conclude about data for population simulation: I have 18 (10.2%) patients 
## that have CRRT == 1 all the time, and 16 patients (9%) that have CRRT == 1 only 
## for 199 days out of 340 days (58.5% of days). The rest of the patients 
## have CRRT == 0.

## Check the availability of APACHE score ---------------------------------

### The percentage of missing APACHE originally ---------------------------

FlucTotIV_clean |> 
  select(APACHE) |> 
  summary()

FlucTotIV_clean |> 
  summarise(missing_percent = mean(is.na(APACHE)) * 100)
# 65.8%, when using the single imputation, which means 65.8% patients do not 
# have the information of APACHE score.

### After removing duplicates in APACHE ----------------------------------------

FlucTotIV_clean |> 
  distinct(APACHE) |> 
  summarise(n = n())
# Only 34/3170 unique APACHE values, which is slightly more than 1%
# Thus, it's just imfeasible to include APACHE as a covariate

## Report BW at a patient level --------------------------------------------

# Extract data that each ID has only 1 observation

FlucTotIV_clean2 <- FlucTotIV_clean[!duplicated(FlucTotIV_clean$ID),]

# Now summarise BW with table1 at a patient level

table1(~ BW|HOSPITAL, data = FlucTotIV_clean2,render.continuous =
         c(. = "Median [Q1-Q3]"))


# Multiple imputation - 70 models ------------------------------------------

## Create datasets for imputation ----------------------------------------------

### Some data manipulation ---------------------------------------------------

# Converting variables into appropriate types

FlucTotIV_clean$SEX <- as.factor(FlucTotIV_clean$SEX)
FlucTotIV_clean$HOSPITAL <- as.factor(FlucTotIV_clean$HOSPITAL)
FlucTotIV_clean$CRRT <- as.factor(FlucTotIV_clean$CRRT)
FlucTotIV_clean$OCC <- as.factor(FlucTotIV_clean$OCC)
FlucTotIV_clean$EVID <- as.factor(FlucTotIV_clean$EVID)

# Converting all the NAs in TAD, AMT, RATE to 0, and First row of DV to 0

FlucTotIV_clean$TAD[is.na(FlucTotIV_clean$TAD)] <- 0
FlucTotIV_clean$AMT[is.na(FlucTotIV_clean$AMT)] <- 0
FlucTotIV_clean$RATE[is.na(FlucTotIV_clean$RATE)] <- 0
FlucTotIV_clean$DV[FlucTotIV_clean$TIME == 0] <- 0

# Converting BW2 into BW

FlucTotIV_clean$BW <- FlucTotIV_clean$BW2

### Make crrt and non-crrt datasets ------------------------------------------

crrt_patients <- unique(FlucTotIV_clean[FlucTotIV_clean$CRRT == 1, "ID"])

FlucTotIV_clean_crrt <- FlucTotIV_clean |>  
                          dplyr::filter(ID %in% crrt_patients)

FlucTotIV_clean_nocrrt <- FlucTotIV_clean |> 
                            dplyr::filter(!ID %in% crrt_patients)

### Best imputation model for crrt and non-crrt datasets ---------------------

impute_CREAT_DOS_INT_crrt <- FlucTotIV_clean_crrt[, c("ID","TIME","RATE","AMT",
                                                      "DV","AGE","SEX","BW",
                                                      "LENGTH", "CRRT","OCC",
                                                      "HOSPITAL","CREAT_DOS_INT")]

impute_CREAT_DOS_INT_nocrrt <- FlucTotIV_clean_nocrrt[, c("ID","TIME","RATE","AMT",
                                                          "DV","AGE","SEX","BW",
                                                          "LENGTH","OCC","HOSPITAL",
                                                          "CREAT_DOS_INT")]

## Performing MI, 2l.pan ----------------------------------------------

# For crrt dataset

imp0 <- mice(impute_CREAT_DOS_INT_crrt, maxit = 0)
predmat <- imp0$predictorMatrix
predmat["CREAT_DOS_INT","ID"] <- -2
meth <- imp0$method
meth["CREAT_DOS_INT"] <- "2l.pan"
maxit <- 20
nimp <- 70
imputed_CREAT_DOS_INT_crrt <- mice::mice(data = impute_CREAT_DOS_INT_crrt,
                                         method = meth, predictorMatrix = predmat,
                                         maxit = maxit, m = nimp, printFlag = F,
                                         seed = 123)
for (i in 1:70) {
  FlucTotIV_clean_crrt[[paste0("CREAT", sprintf("%02d", i))]] <- (complete(imputed_CREAT_DOS_INT_crrt, action = i))$CREAT_DOS_INT
}

# For non-crrt dataset

imp0 <- mice(impute_CREAT_DOS_INT_nocrrt, maxit = 0)
predmat <- imp0$predictorMatrix
predmat["CREAT_DOS_INT", "ID"] <- -2
meth <- imp0$method
meth["CREAT_DOS_INT"] <- "2l.pan"
maxit <- 20
nimp <- 70
imputed_CREAT_DOS_INT_nocrrt <- mice::mice(data = impute_CREAT_DOS_INT_nocrrt,
                                           method = meth,predictorMatrix = predmat,
                                           maxit = maxit, m = nimp, printFlag = F,
                                           seed = 123)
for (i in 1:70) {
  FlucTotIV_clean_nocrrt[[paste0("CREAT", sprintf("%02d", i))]] <- (complete(imputed_CREAT_DOS_INT_nocrrt, action = i))$CREAT_DOS_INT
}

### BW and LENGTH treatment -------------------------------------------------

## Adding BW and LENGTH to crrt and non-crrt datasets

# Define a function to add variables to a dataset
add_variables <- function(dataset, imputed_data, var_name) {
  for (i in 1:70) {
    new_var_name <- paste0(var_name, sprintf("%02d", i))
    dataset[[new_var_name]] <- complete(imputed_data, action = i)[[var_name]]
  }
  return(dataset)
}

# Apply the function to FlucTotIV_clean_crrt and FlucTotIV_clean_nocrrt
FlucTotIV_clean_crrt <- add_variables(FlucTotIV_clean_crrt, imputed_CREAT_DOS_INT_crrt, "BW")
FlucTotIV_clean_crrt <- add_variables(FlucTotIV_clean_crrt, imputed_CREAT_DOS_INT_crrt, "LENGTH")

FlucTotIV_clean_nocrrt <- add_variables(FlucTotIV_clean_nocrrt, imputed_CREAT_DOS_INT_nocrrt, "BW")
FlucTotIV_clean_nocrrt <- add_variables(FlucTotIV_clean_nocrrt, imputed_CREAT_DOS_INT_nocrrt, "LENGTH")

## Combine 2 datasets into 1 FlucTotIV_clean_imputed dataset

FlucTotIV_clean_imputed <- rbind(FlucTotIV_clean_crrt, FlucTotIV_clean_nocrrt)
FlucTotIV_clean_imputed <- FlucTotIV_clean_imputed[order(FlucTotIV_clean_imputed$ID, 
                                                         FlucTotIV_clean_imputed$HOSPITAL, 
                                                         FlucTotIV_clean_imputed$TIME), ]

## Now average BW and LENGTH per each ID

# Creating subset of hospital 3 and 4 patients
subset <- FlucTotIV_clean_imputed |>  
            dplyr::filter(HOSPITAL %in% c(3,4))

# Creating new variables for mean values
bw_mean_cols <- paste0("BW", sprintf("%02d", 1:70), "_mean")
length_mean_cols <- paste0("LENGTH", sprintf("%02d", 1:70), "_mean")
subset[, bw_mean_cols] <- NA
subset[, length_mean_cols] <- NA

# Filling mean values per ID
for (i in unique(subset$ID)) {
  for (j in 1:70) {
    bw_col <- paste0("BW", sprintf("%02d", 1:70))
    length_col <- paste0("LENGTH", sprintf("%02d", 1:70))
    bw_mean_col <- paste0("BW", sprintf("%02d", 1:70), "_mean")
    length_mean_col <- paste0("LENGTH", sprintf("%02d", 1:70), "_mean")
    subset[subset$ID == i, bw_mean_col[j]] <- mean(subset[subset$ID == i, ][[bw_col[j]]], na.rm = TRUE)
    subset[subset$ID == i, length_mean_col[j]] <- mean(subset[subset$ID == i, ][[length_col[j]]], na.rm = TRUE)
  }
}

# Integrating mean values to the original dataset
FlucTotIV_clean_imputed[, bw_mean_cols] <- NA
FlucTotIV_clean_imputed[, length_mean_cols] <- NA
FlucTotIV_clean_imputed[FlucTotIV_clean_imputed$ID %in% subset$ID, bw_mean_cols] <- subset[, bw_mean_cols]
FlucTotIV_clean_imputed[FlucTotIV_clean_imputed$ID %in% subset$ID, length_mean_cols] <- subset[, length_mean_cols]

# Replace values of 70 columns with mean values in hospital 3 and 4
for (j in 1:70) {
  bw_col <- paste0("BW", sprintf("%02d", j))
  length_col <- paste0("LENGTH", sprintf("%02d", j))
  bw_mean_col <- paste0("BW", sprintf("%02d", j), "_mean")
  length_mean_col <- paste0("LENGTH", sprintf("%02d", j), "_mean")
  FlucTotIV_clean_imputed[FlucTotIV_clean_imputed$HOSPITAL %in% c(3, 4), bw_col] <- FlucTotIV_clean_imputed[FlucTotIV_clean_imputed$HOSPITAL %in% c(3, 4), bw_mean_col]
  FlucTotIV_clean_imputed[FlucTotIV_clean_imputed$HOSPITAL %in% c(3, 4), length_col] <- FlucTotIV_clean_imputed[FlucTotIV_clean_imputed$HOSPITAL %in% c(3, 4), length_mean_col]
}

# Delete the 70 mean columns
bw_mean_cols <- paste0("BW", sprintf("%02d", 1:70), "_mean")
length_mean_cols <- paste0("LENGTH", sprintf("%02d", 1:70), "_mean")
FlucTotIV_clean_imputed <- FlucTotIV_clean_imputed |>  select(-one_of(bw_mean_cols), -one_of(length_mean_cols))

### CREAT treatment ---------------------------------------------------------

# Create a subset of rows with missing CREAT_DOS_INT
subset <- FlucTotIV_clean_imputed |>  
            dplyr::filter(is.na(CREAT_DOS_INT))

# Calculate the mean values per ID and OCC
mean_cols <- paste0("CREAT", sprintf("%02d", 1:70))
subset <- subset |> 
  group_by(ID, OCC) |> 
  mutate(across(mean_cols, mean, na.rm = TRUE)) |> 
  ungroup()

# Replace the mean across ID and OCC in the actual dataset
FlucTotIV_clean_imputed[is.na(FlucTotIV_clean_imputed$CREAT_DOS_INT), mean_cols] <- subset[, mean_cols]

# Fixing all the negative imputed values to the minimum observed value
FlucTotIV_clean_imputed <- FlucTotIV_clean_imputed |> 
  mutate_at(vars(CREAT01:CREAT70), ~ ifelse(. < 0, min(CREAT,na.rm=TRUE), .))

### MI diagnostics ---------------------------------------------------------

## Make a violin plot

# Select the 11 columns by specifying their names explicitly
CREAT_columns <- FlucTotIV_clean_imputed |> 
  select(CREAT_DOS_INT, CREAT01, CREAT02, CREAT03, CREAT04, CREAT05, CREAT06, CREAT07, CREAT08, CREAT09, CREAT10)

# Convert from wide to long
CREAT_long_data <- CREAT_columns |> 
  pivot_longer(cols = everything(), 
               names_to = "Imputation", 
               values_to = "Value")

# Rename the levels of the "Imputation" column
CREAT_long_data$Imputation <- factor(CREAT_long_data$Imputation, 
                                     levels = c("CREAT_DOS_INT", "CREAT01", 
                                                "CREAT02", "CREAT03", "CREAT04", 
                                                "CREAT05", "CREAT06", "CREAT07", 
                                                "CREAT08", "CREAT09", "CREAT10"))

# Remove NA from CREAT_long_data
CREAT_long_data <- na.omit(CREAT_long_data)

# Define a custom color palette
#custom_palette <- c("skyblue", "lightgreen", "pink", "orange", "purple",
                    #"lightblue", "salmon", "gold", "lightpink", "lightgray", "lightcyan")

# Create a violin plot with custom aesthetics
violin_plot<- ggplot(data = CREAT_long_data, aes(x = Imputation, y = Value, fill = Imputation)) +
  geom_violin(trim = TRUE) +
  labs(title = "Distribution of observed and 10 imputed serum creatinine",
       x = NULL,
       y = "Serum Creatinine (mg/dl)") +
  scale_x_discrete(labels = c("Observed", "Imputation 01", "Imputation 02", 
                              "Imputation 03", "Imputation 04", "Imputation 05", 
                              "Imputation 06", "Imputation 07", "Imputation 08", 
                              "Imputation 09", "Imputation 10")) +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, color = "black", size = 16, face = "bold"),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 14, face = "bold"))

# Exporting the plot

setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Plots/MI_diagnostics")
ggsave("obs_vs_imp_violin.png", plot = violin_plot, dpi = 300, width = 13, height = 7)

## Make a density plot

# Define the title and label options
title_options <- list(
  label = "Density plot of 10 imputed SCr vs the original SCr",
  size = 14,
  color = "black",
  hjust = 0.5,
  vjust = 1,
  face = "bold"
)

x_label_options <- list(
  label = "SCr (mg/dl)",
  size = 12,
  color = "black",
  hjust = 0.5,
  vjust = 0
)

y_label_options <- list(
  label = "Density",
  size = 12,
  color = "black",
  hjust = 0,
  vjust = 0.5,
  angle = 0
)

axis_text_options <- list(
  size = 10
)

legend_options <- list(
  title = "Variables",
  size = 8
)

# Density plot first 10 imputations 
density_plot<-ggplot(FlucTotIV_clean_imputed, aes(x = CREAT01)) + 
  geom_density(aes(fill = "CREAT01"), alpha = 0.7) + 
  geom_density(aes(x = CREAT02, fill = "CREAT02"), alpha = 0.7) + 
  geom_density(aes(x = CREAT03, fill = "CREAT03"), alpha = 0.7) + 
  geom_density(aes(x = CREAT04, fill = "CREAT04"), alpha = 0.7) + 
  geom_density(aes(x = CREAT05, fill = "CREAT05"), alpha = 0.7) + 
  geom_density(aes(x = CREAT06, fill = "CREAT06"), alpha = 0.7) + 
  geom_density(aes(x = CREAT07, fill = "CREAT07"), alpha = 0.7) + 
  geom_density(aes(x = CREAT08, fill = "CREAT08"), alpha = 0.7) + 
  geom_density(aes(x = CREAT09, fill = "CREAT09"), alpha = 0.7) + 
  geom_density(aes(x = CREAT10, fill = "CREAT10"), alpha = 0.7) +
  geom_density(aes(x = CREAT_DOS_INT, fill = "CREAT"), alpha = 0.7) +
  theme_bw() +
  scale_fill_viridis_d(
    name = "Variables",
    option = "A"
  ) + 
  labs(title = title_options,
           x = x_label_options,
           y = y_label_options) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(axis_text_options),
    legend.position = "right"
  )

# Exporting the plot

setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Plots/MI_diagnostics")
ggsave("density_plot.png", plot = density_plot, dpi = 300, width = 13, height = 7)

### Calculating the CKDEPI from CREAT -------------------------------

Fluc_NONMEM_mul_impute <- FlucTotIV_clean_imputed

# Define a function to calculate CKDEPI for each CKDEPI column
calculate_CKDEPI <- function(data, creatinine_column, sex_column, age_column, multiplier) {
  CKDEPI <- ifelse(data[[sex_column]] == 1 & data[[creatinine_column]] < 0.9,
                   multiplier * ((data[[creatinine_column]]/0.9)^-0.302) * (0.9938^data[[age_column]]),
                   ifelse(data[[sex_column]] == 1 & data[[creatinine_column]] >= 0.9,
                          multiplier * ((data[[creatinine_column]]/0.9)^-1.200) * (0.9938^data[[age_column]]),
                          ifelse(data[[sex_column]] == 2 & data[[creatinine_column]] < 0.7,
                                 multiplier * ((data[[creatinine_column]]/0.7)^-0.241) * (0.9938^data[[age_column]]) * 1.012,
                                 multiplier * ((data[[creatinine_column]]/0.7)^-1.200) * (0.9938^data[[age_column]]) * 1.012)
                   )
  )
  return(CKDEPI)
}

# Calculate CKDEPI01 to CKDEPI70
for (i in 1:70) {
  CKDEPI_column <- paste0("CKDEPI", sprintf("%02d", i))
  Fluc_NONMEM_mul_impute[[CKDEPI_column]] <- calculate_CKDEPI(Fluc_NONMEM_mul_impute, 
                                                              paste0("CREAT", sprintf("%02d", i)), 
                                                              "SEX", "AGE", 142)
}

# View the resulting dataset
View(Fluc_NONMEM_mul_impute)

## Treatment after MI ----------------------------------------------

# Converting all NAs to -99
Fluc_NONMEM_mul_impute[is.na(Fluc_NONMEM_mul_impute)] <- -99

# Convert DV from 0 & -99 to "."
Fluc_NONMEM_mul_impute$DV <- ifelse(Fluc_NONMEM_mul_impute$DV %in% c(0, -99), ".", Fluc_NONMEM_mul_impute$DV)

# AMT TAD and RATE from 0 to "."
Fluc_NONMEM_mul_impute$TAD <- ifelse(Fluc_NONMEM_mul_impute$TAD == 0, ".", Fluc_NONMEM_mul_impute$TAD)
Fluc_NONMEM_mul_impute$AMT <- ifelse(Fluc_NONMEM_mul_impute$AMT == 0, ".", Fluc_NONMEM_mul_impute$AMT)
Fluc_NONMEM_mul_impute$RATE <- ifelse(Fluc_NONMEM_mul_impute$RATE == 0, ".", Fluc_NONMEM_mul_impute$RATE)

### Export into different datasets for modelling -------------------

# Remove variables that are not needed
FlucTotIV_clean_imputed <- Fluc_NONMEM_mul_impute[, !names(Fluc_NONMEM_mul_impute) 
                                                  %in% c("CMT", "APACHE", "RACE", 
                                                        "CL24", "GGT", "AFT", "ALT", 
                                                        "AST", "BILI", "ALB", "SOFA", 
                                                        "IHD", "UF", "ECMO", "OCC3", 
                                                        "ARCCKD", "ARC24", "ARCAlg", 
                                                        "BMI", "BSA", "dup","FLUID",
                                                        paste0("CREAT", sprintf("%02d", 1:70)),
                                                        paste0("LENGTH", sprintf("%02d", 1:70)),
                                                        "CREAT_DOS_INT","unique_CREAT",
                                                        "CREAT3","CRRT2","dosing_interval","Creat_imputed",
                                                        "OCC", "TAD", "SS","II",
                                                        "ADMIN", "HOSPITAL", "LENGTH",
                                                        "AGE", "SEX", "CREAT",
                                                        "CKDEPI", "BW", "CKDEPI",
                                                        "BW2", "CREAT2")]

# Then, make Fluco_revised_imputed_01 dataset, which is FlucTotIV_clean_imputed
# removing BW19 to BW70 and CKDEPI19 to CKDEPI70

Fluco_revised_imputed_01 <- FlucTotIV_clean_imputed[, !names(FlucTotIV_clean_imputed) %in% 
                                                   c(paste0("BW", sprintf("%02d", 19:70)),
                                                        paste0("CKDEPI", sprintf("%02d", 19:70)))]

# Next, make Fluco_revised_imputed_02 dataset, which is FlucTotIV_clean_imputed
# removing BW01 to BW18 and BW37 to BW70; similarly for CKDEPI

Fluco_revised_imputed_02 <- FlucTotIV_clean_imputed[, !names(FlucTotIV_clean_imputed) %in%
                                                   c(paste0("BW", sprintf("%02d", c(1:18, 37:70))),
                                                        paste0("CKDEPI", sprintf("%02d", c(1:18, 37:70))))]

# Next, make Fluco_revised_imputed_03 dataset, which is FlucTotIV_clean_imputed
# removing BW01 to BW36 and BW55 to BW70; similarly for CKDEPI

Fluco_revised_imputed_03 <- FlucTotIV_clean_imputed[, !names(FlucTotIV_clean_imputed) %in% 
                                                     c(paste0("BW", sprintf("%02d", c(1:36, 55:70))),
                                                          paste0("CKDEPI", sprintf("%02d", c(1:36, 55:70))))]

# Finally, make Fluco_revised_imputed_04 dataset, which is FlucTotIV_clean_imputed
# removing BW01 to BW54 and CKDEPI01 to CKDEPI54

Fluco_revised_imputed_04 <- FlucTotIV_clean_imputed[, !names(FlucTotIV_clean_imputed) %in%
                                                     c(paste0("BW", sprintf("%02d", 01:54)),
                                                          paste0("CKDEPI", sprintf("%02d", 01:54)))]

# Export the datasets

setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Imputed")
write.csv(Fluco_revised_imputed_01, "Fluco_revised_imputed_01.csv", quote = F, row.names = F)
write.csv(Fluco_revised_imputed_02, "Fluco_revised_imputed_02.csv", quote = F, row.names = F)
write.csv(Fluco_revised_imputed_03, "Fluco_revised_imputed_03.csv", quote = F, row.names = F)
write.csv(Fluco_revised_imputed_04, "Fluco_revised_imputed_04.csv", quote = F, row.names = F)

### Generate the overall dataset for making GOF plots and VPC -----------------

# Generate two new columns in Fluc_NONMEM_mul_impute, BW_ave and CKDEPI_ave
# BW_ave is the average of BW01 to BW70, and similarly for CKDEPI_ave

Fluc_NONMEM_mul_impute$BW_ave <- rowMeans(Fluc_NONMEM_mul_impute[,paste0("BW", sprintf("%02d", 1:70))], na.rm = TRUE)
                                                                 
Fluc_NONMEM_mul_impute$CKDEPI_ave <- rowMeans(Fluc_NONMEM_mul_impute[,paste0("CKDEPI", sprintf("%02d", 1:70))], na.rm = TRUE)                                                                 

# Then make Fluco_revised_imputed_overall dataset, which should have TAD & HOSPITAL 
# to make GOF plots
Fluco_revised_imputed_overall <- Fluc_NONMEM_mul_impute[, !names(Fluc_NONMEM_mul_impute) %in% 
                                                        c("CMT", "APACHE", "RACE", 
                                                        "CL24", "GGT", "AFT", "ALT", 
                                                        "AST", "BILI", "ALB", "SOFA", 
                                                        "IHD", "UF", "ECMO", "OCC3", 
                                                        "ARCCKD", "ARC24", "ARCAlg", 
                                                        "BMI", "BSA", "dup","FLUID",
                                                        paste0("CREAT", sprintf("%02d", 1:70)),
                                                        paste0("LENGTH", sprintf("%02d", 1:70)),
                                                        paste0("BW", sprintf("%02d", 1:70)),
                                                        paste0("CKDEPI", sprintf("%02d", 1:70)),
                                                        "CREAT_DOS_INT","unique_CREAT",
                                                        "CREAT3","CRRT2","dosing_interval","Creat_imputed",
                                                        "OCC", "SS","II",
                                                        "ADMIN", "LENGTH",
                                                        "AGE", "SEX", "CREAT",
                                                        "CKDEPI", "BW", "CKDEPI",
                                                        "BW2", "CREAT2")]

# Export the dataset

setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Imputed")
write.csv(Fluco_revised_imputed_overall, "Fluco_revised_imputed_overall.csv", quote = F, row.names = F)


## Pooling model parameters ----------------------------------------------

# Load the dataset
setwd("./Datasets/Pooled_MI")
Pooling_Parameters_2l.pan <- read.csv("Pooling_Parameters_2l.pan_220324.csv",sep = ",")

# Creating Var (Variance) column 
Pooling_Parameters_2l.pan$Var <- (Pooling_Parameters_2l.pan$RSE*Pooling_Parameters_2l.pan$Estimates/100)^2

### First of all, CKDEPI -----------------

# Mean estimate
CKDEPI <- mean(Pooling_Parameters_2l.pan$Estimates[Pooling_Parameters_2l.pan$Parameters == "CKDEPI"]) #0.5316571

W_CKDEPI <- mean(Pooling_Parameters_2l.pan$Var[Pooling_Parameters_2l.pan$Parameters == "CKDEPI"]) #0.03150853 #within imputation variance

B_CKDEPI <- var(Pooling_Parameters_2l.pan$Estimates[Pooling_Parameters_2l.pan$Parameters == "CKDEPI"]) #0.00929304 #between imputation variance

CKDEPI
W_CKDEPI
B_CKDEPI

# Variance estimate
m <- 70
T_CKDEPI <- W_CKDEPI + (1 + 1/m) * B_CKDEPI #0.04093433
T_CKDEPI

# Confidence interval - T distribution
m <- 70 #number of imputations
k <- 12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect parameters where covariance is one of them)
n <- 177

lambda_CKDEPI <- (B_CKDEPI + B_CKDEPI/m)/T_CKDEPI #proportion of variation attributable to the missing data
r_CKDEPI <- (lambda_CKDEPI)/(1 - lambda_CKDEPI) #relative increase in variance due to nonresponse
v_old_CKDEPI <- (m - 1) * (1 + 1/r_CKDEPI^2) #old degree of freedom
v_com = n - k #degree of freedom of parameter estimate (CKDEPI) in the hypothetically complete data
v_obs_CKDEPI = ((v_com + 1)/(v_com + 3)) * v_com * (1 - lambda_CKDEPI) #observed data degrees of freedom that accounts for the missing information
v_CKDEPI <- (v_old_CKDEPI*v_obs_CKDEPI)/(v_old_CKDEPI + v_obs_CKDEPI) 
t_crit_CKDEPI <- qt(0.975, v_CKDEPI) 

lower_bound_CKDEPI <- CKDEPI - t_crit_CKDEPI*sqrt(T_CKDEPI) #0.1306684
upper_bound_CKDEPI <- CKDEPI + t_crit_CKDEPI*sqrt(T_CKDEPI) #0.9326459

lower_bound_CKDEPI
upper_bound_CKDEPI

### Second of all, CLcrrt -----------------

# Mean estimate
CLcrrt <- mean(Pooling_Parameters_2l.pan$Estimates[Pooling_Parameters_2l.pan$Parameters == "CLcrrt"]) #1.561

W_CLcrrt <- mean(Pooling_Parameters_2l.pan$Var[Pooling_Parameters_2l.pan$Parameters == "CLcrrt"]) #0.02459816 #within imputation variance

B_CLcrrt <- var(Pooling_Parameters_2l.pan$Estimates[Pooling_Parameters_2l.pan$Parameters == "CLcrrt"]) #2.362319e-05 #between imputation variance

CLcrrt
W_CLcrrt
B_CLcrrt

# Variance estimate
m <- 70
T_CLcrrt <- W_CLcrrt + (1 + 1/m)*B_CLcrrt #0.02462213
T_CLcrrt

# Confidence interval - T distribution
m <- 70 #number of imputations
k <- 12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect parameters where covariance is one of them)
n <- 177

lambda_CLcrrt <- (B_CLcrrt + B_CLcrrt/m)/T_CLcrrt #proportion of variation attributable to the missing data

r_CLcrrt <- (lambda_CLcrrt)/(1 - lambda_CLcrrt) #relative increase in variance due to nonresponse

v_old_CLcrrt <- (m - 1)*(1 + 1/r_CLcrrt^2) #old degree of freedom

v_com = n - k #degree of freedom of parameter estimate (CLcrrt) in the hypothetically complete data

v_obs_CLcrrt = ((v_com + 1)/(v_com + 3))*v_com*(1 - lambda_CLcrrt) #observed data degrees of freedom that accounts for the missing information

v_CLcrrt <- (v_old_CLcrrt*v_obs_CLcrrt)/(v_old_CLcrrt + v_obs_CLcrrt)

t_crit_CLcrrt <- qt(0.975, v_CLcrrt) 

lower_bound_CLcrrt <- CLcrrt - t_crit_CLcrrt*sqrt(T_CLcrrt) #1.251151
upper_bound_CLcrrt <- CLcrrt + t_crit_CLcrrt*sqrt(T_CLcrrt) #1.870849

lower_bound_CLcrrt
upper_bound_CLcrrt

### Third of all, CLr -----------------

# Mean estimate
CLr <- mean(Pooling_Parameters_2l.pan$Estimates[Pooling_Parameters_2l.pan$Parameters == "CLr"]) #0.6137429

W_CLr <- mean(Pooling_Parameters_2l.pan$Var[Pooling_Parameters_2l.pan$Parameters == "CLr"]) #0.001016242 #within imputation variance

B_CLr <- var(Pooling_Parameters_2l.pan$Estimates[Pooling_Parameters_2l.pan$Parameters == "CLr"]) #8.886046e-05
#between imputation variance

CLr
W_CLr
B_CLr

# Variance estimate
m <- 70
T_CLr <- W_CLr + (1 + 1/m)*B_CLr #0.001106372
T_CLr

# Confidence interval - T distribution
m <- 70 #number of imputations
k <- 12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect parameters where covariance is one of them)
n <- 177

lambda_CLr <- (B_CLr + B_CLr/m)/T_CLr #proportion of variation attributable to the missing data

r_CLr <- (lambda_CLr)/(1 - lambda_CLr) #relative increase in variance due to nonresponse

v_old_CLr <- (m - 1)*(1 + 1/r_CLr^2) #old degree of freedom

v_com = n - k #degree of freedom of parameter estimate (CLr) in the hypothetically complete data

v_obs_CLr = ((v_com + 1)/(v_com + 3))*v_com*(1 - lambda_CLr) #observed data degrees of freedom that accounts for the missing information

v_CLr <- (v_old_CLr*v_obs_CLr)/(v_old_CLr + v_obs_CLr)

t_crit_CLr <- qt(0.975, v_CLr) 

lower_bound_CLr <- CLr - t_crit_CLr*sqrt(T_CLr) #0.54801
upper_bound_CLr <- CLr + t_crit_CLr*sqrt(T_CLr) #0.6794757

lower_bound_CLr
upper_bound_CLr

### Next, V1 -----------------

# Mean estimate
V1 <- mean(Pooling_Parameters_2l.pan$Estimates[Pooling_Parameters_2l.pan$Parameters == "V1"]) #39.39571

W_V1 <- mean(Pooling_Parameters_2l.pan$Var[Pooling_Parameters_2l.pan$Parameters == "V1"]) #10.21739 #within imputation variance

B_V1 <- var(Pooling_Parameters_2l.pan$Estimates[Pooling_Parameters_2l.pan$Parameters == "V1"]) #0.09345963
#between imputation variance

V1
W_V1
B_V1

# Variance estimate
m <- 70
T_V1 <- W_V1 + (1 + 1/m)*B_V1 #10.31219
T_V1

# Confidence interval - T distribution
m <- 70 #number of imputations
k <- 12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect parameters where covariance is one of them)
n <- 177

lambda_V1 <- (B_V1 + B_V1/m)/T_V1 #proportion of variation attributable to the missing data

r_V1 <- (lambda_V1)/(1 - lambda_V1) #relative increase in variance due to nonresponse

v_old_V1 <- (m - 1)*(1 + 1/r_V1^2) #old degree of freedom

v_com = n - k #degree of freedom of parameter estimate (V1) in the hypothetically complete data

v_obs_V1 = ((v_com + 1)/(v_com + 3))*v_com*(1 - lambda_V1) #observed data degrees of freedom that accounts for the missing information

v_V1 <- (v_old_V1*v_obs_V1)/(v_old_V1 + v_obs_V1)

t_crit_V1 <- qt(0.975, v_V1) 

lower_bound_V1 <- V1 - t_crit_V1*sqrt(T_V1) #33.05424
upper_bound_V1 <- V1 + t_crit_V1*sqrt(T_V1) #45.73719

lower_bound_V1
upper_bound_V1

### After that, Q -----------------

# Mean estimate
Q <- mean(Pooling_Parameters_2l.pan$Estimates[Pooling_Parameters_2l.pan$Parameters == "Q"]) #12.11286

W_Q <- mean(Pooling_Parameters_2l.pan$Var[Pooling_Parameters_2l.pan$Parameters == "Q"]) #18.14542 #within imputation variance

B_Q <- var(Pooling_Parameters_2l.pan$Estimates[Pooling_Parameters_2l.pan$Parameters == "Q"]) #0.03620911
#between imputation variance

Q
W_Q
B_Q

# Variance estimate
m <- 70
T_Q <- W_Q + (1 + 1/m)*B_Q #18.18215
T_Q

# Confidence interval - T distribution
m <- 70 #number of imputations
k <- 12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect parameters where covariance is one of them)
n <- 177

lambda_Q <- (B_Q + B_Q/m)/T_Q #proportion of variation attributable to the missing data

r_Q <- (lambda_Q)/(1 - lambda_Q) #relative increase in variance due to nonresponse

v_old_Q <- (m - 1)*(1 + 1/r_Q^2) #old degree of freedom

v_com = n - k #degree of freedom of parameter estimate (Q) in the hypothetically complete data

v_obs_Q = ((v_com + 1)/(v_com + 3))*v_com*(1 - lambda_Q) #observed data degrees of freedom that accounts for the missing information

v_Q <- (v_old_Q * v_obs_Q)/(v_old_Q + v_obs_Q)

t_crit_Q <- qt(0.975, v_Q) 

lower_bound_Q <- Q - t_crit_Q * sqrt(T_Q) #3.692838
upper_bound_Q <- Q + t_crit_Q * sqrt(T_Q) #20.53288

lower_bound_Q
upper_bound_Q

### Next, V2 -----------------

# Mean estimate
V2 <- mean(Pooling_Parameters_2l.pan$Estimates[Pooling_Parameters_2l.pan$Parameters == "V2"]) #8.465714

W_V2 <- mean(Pooling_Parameters_2l.pan$Var[Pooling_Parameters_2l.pan$Parameters == "V2"]) #8.436197 #within imputation variance

B_V2 <- var(Pooling_Parameters_2l.pan$Estimates[Pooling_Parameters_2l.pan$Parameters == "V2"]) #0.006737888
#between imputation variance

V2
W_V2
B_V2

# Variance estimate
m <- 70
T_V2 <- W_V2 + (1 + 1/m)*B_V2 #8.443031
T_V2

# Confidence interval - T distribution
m <- 70 #number of imputations
k <- 12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect parameters where covariance is one of them)
n <- 177

lambda_V2 <- (B_V2 + B_V2/m)/T_V2 #proportion of variation attributable to the missing data

r_V2 <- (lambda_V2)/(1 - lambda_V2) #relative increase in variance due to nonresponse

v_old_V2 <- (m - 1)*(1 + 1/r_V2^2) #old degree of freedom

v_com = n - k #degree of freedom of parameter estimate (V2) in the hypothetically complete data

v_obs_V2 = ((v_com + 1)/(v_com + 3))*v_com*(1 - lambda_V2) #observed data degrees of freedom that accounts for the missing information

v_V2 <- (v_old_V2*v_obs_V2)/(v_old_V2 + v_obs_V2)

t_crit_V2 <- qt(0.975, v_V2) 

lower_bound_V2 <- V2 - t_crit_V2*sqrt(T_V2) #2.728043
upper_bound_V2 <- V2 + t_crit_V2*sqrt(T_V2) #14.20339

lower_bound_V2
upper_bound_V2

### Next, BW -----------------

# Mean estimate
BW <- mean(Pooling_Parameters_2l.pan$Estimates[Pooling_Parameters_2l.pan$Parameters == "BW"]) #0.9076

W_BW <- mean(Pooling_Parameters_2l.pan$Var[Pooling_Parameters_2l.pan$Parameters == "BW"]) #0.04289933 #within imputation variance

B_BW <- var(Pooling_Parameters_2l.pan$Estimates[Pooling_Parameters_2l.pan$Parameters == "BW"]) #0.0001656638
#between imputation variance

BW
W_BW
B_BW

# Variance estimate
m <- 70
T_BW <- W_BW + (1 + 1/m)*B_BW #0.04306736
T_BW

# Confidence interval - T distribution
m <- 70 #number of imputations
k <- 12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect parameters where covariance is one of them)
n <- 177

lambda_BW <- (B_BW + B_BW/m)/T_BW #proportion of variation attributable to the missing data

r_BW <- (lambda_BW)/(1 - lambda_BW) #relative increase in variance due to nonresponse

v_old_BW <- (m - 1)*(1 + 1/r_BW^2) #old degree of freedom

v_com = n - k #degree of freedom of parameter estimate (BW) in the hypothetically complete data

v_obs_BW = ((v_com + 1)/(v_com + 3))*v_com*(1 - lambda_BW) #observed data degrees of freedom that accounts for the missing information

v_BW <- (v_old_BW * v_obs_BW)/(v_old_BW + v_obs_BW)

t_crit_BW <- qt(0.975, v_BW) 

lower_bound_BW <- BW - t_crit_BW*sqrt(T_BW) #0.4978011
upper_bound_BW <- BW + t_crit_BW*sqrt(T_BW) #1.317399

lower_bound_BW
upper_bound_BW

### Next, Cov_CLr_V1 -----------------

# First of all, replacing the Var with the values in RSE column
Pooling_Parameters_2l.pan <- Pooling_Parameters_2l.pan |> 
  mutate(Var = ifelse(Parameters == "Cov_CLr_V1", RSE^2, Var))

# Second of all, replacing the values of the RSE column
Pooling_Parameters_2l.pan <- Pooling_Parameters_2l.pan |> 
  mutate(RSE = ifelse(Parameters == "Cov_CLr_V1", (RSE/Estimates)*100, RSE))

# Mean estimate
Cov_CLr_V1 <- mean(Pooling_Parameters_2l.pan$Estimates[Pooling_Parameters_2l.pan$Parameters == "Cov_CLr_V1"]) #0.06732857

W_Cov_CLr_V1 <- mean(Pooling_Parameters_2l.pan$Var[Pooling_Parameters_2l.pan$Parameters == "Cov_CLr_V1"]) #0.001466198 #within imputation variance

B_Cov_CLr_V1 <- var(Pooling_Parameters_2l.pan$Estimates[Pooling_Parameters_2l.pan$Parameters == "Cov_CLr_V1"]) #1.70554e-05
#between imputation variance

Cov_CLr_V1
W_Cov_CLr_V1
B_Cov_CLr_V1

# Variance estimate
m <- 70
T_Cov_CLr_V1 <- W_Cov_CLr_V1 + (1 + 1/m)*B_Cov_CLr_V1 #0.001483497
T_Cov_CLr_V1

# Confidence interval - T distribution
m <- 70 #number of imputations
k <- 12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect parameters where covariance is one of them)
n <- 177

lambda_Cov_CLr_V1 <- (B_Cov_CLr_V1 + B_Cov_CLr_V1/m)/T_Cov_CLr_V1 #proportion of variation attributable to the missing data

r_Cov_CLr_V1 <- (lambda_Cov_CLr_V1)/(1 - lambda_Cov_CLr_V1) #relative increase in variance due to nonresponse

v_old_Cov_CLr_V1 <- (m - 1)*(1 + 1/r_Cov_CLr_V1^2) #old degree of freedom

v_com = n - k #degree of freedom of parameter estimate (Cov_CLr_V1) in the hypothetically complete data

v_obs_Cov_CLr_V1 = ((v_com + 1)/(v_com + 3))*v_com*(1 - lambda_Cov_CLr_V1) #observed data degrees of freedom that accounts for the missing information

v_Cov_CLr_V1 <- (v_old_Cov_CLr_V1*v_obs_Cov_CLr_V1)/(v_old_Cov_CLr_V1 + v_obs_Cov_CLr_V1)

t_crit_Cov_CLr_V1 <- qt(0.975, v_Cov_CLr_V1) 

lower_bound_Cov_CLr_V1 <- Cov_CLr_V1 - t_crit_Cov_CLr_V1*sqrt(T_Cov_CLr_V1) #-0.008733233
upper_bound_Cov_CLr_V1 <- Cov_CLr_V1 + t_crit_Cov_CLr_V1*sqrt(T_Cov_CLr_V1) #0.1433904

lower_bound_Cov_CLr_V1
upper_bound_Cov_CLr_V1

### Do transformation for log terms -----------------

# create a log-transformed Log column
Pooling_Parameters_2l.pan$Log <- log(Pooling_Parameters_2l.pan$Estimates)

# Variance of the log-transformed variable
Pooling_Parameters_2l.pan$Var_Log <- (1/(Pooling_Parameters_2l.pan$Estimates))^2*Pooling_Parameters_2l.pan$Var 
#Var(log(x))=[1/(mean(x))]^2*var(x) (log(x) in R is actually ln(x))

### IIV_CLcrrt -----------------

# Mean estimate
IIV_CLcrrt_log <- mean(Pooling_Parameters_2l.pan$Log[Pooling_Parameters_2l.pan$Parameters == "IIV_CLcrrt"]) #-1.264295

W_IIV_CLcrrt_log <- mean(Pooling_Parameters_2l.pan$Var_Log[Pooling_Parameters_2l.pan$Parameters == "IIV_CLcrrt"]) #0.1782487 #within imputation variance

B_IIV_CLcrrt_log <- var(Pooling_Parameters_2l.pan$Log[Pooling_Parameters_2l.pan$Parameters == "IIV_CLcrrt"]) #0.0001359124 #between imputation variance

IIV_CLcrrt_log
W_IIV_CLcrrt_log
B_IIV_CLcrrt_log

# Variance estimate
m <- 70
T_IIV_CLcrrt_log <- W_IIV_CLcrrt_log + (1 + 1/m)*B_IIV_CLcrrt_log #0.1783865
T_IIV_CLcrrt_log

# Confidence interval - T distribution
m <- 70 #number of imputations
k <- 12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect parameters where covariance is one of them)
n <- 177

lambda_IIV_CLcrrt_log <- (B_IIV_CLcrrt_log + B_IIV_CLcrrt_log/m)/T_IIV_CLcrrt_log #proportion of variation attributable to the missing data
 
r_IIV_CLcrrt_log <- (lambda_IIV_CLcrrt_log)/(1 - lambda_IIV_CLcrrt_log) #relative increase in variance due to nonresponse

v_old_IIV_CLcrrt_log <- (m - 1)*(1 + 1/r_IIV_CLcrrt_log^2) #old degree of freedom

v_com = n - k #degree of freedom of parameter estimate (IIV_CLcrrt_log) in the hypothetically complete data

v_obs_IIV_CLcrrt_log = ((v_com + 1)/(v_com + 3))*v_com*(1 - lambda_IIV_CLcrrt_log) #observed data degrees of freedom that accounts for the missing information

v_IIV_CLcrrt_log <- (v_old_IIV_CLcrrt_log*v_obs_IIV_CLcrrt_log)/(v_old_IIV_CLcrrt_log + v_obs_IIV_CLcrrt_log)

t_crit_IIV_CLcrrt_log <- qt(0.975, v_IIV_CLcrrt_log) 

lower_bound_IIV_CLcrrt_log <- IIV_CLcrrt_log - t_crit_IIV_CLcrrt_log*sqrt(T_IIV_CLcrrt_log) #-2.098298
upper_bound_IIV_CLcrrt_log <- IIV_CLcrrt_log + t_crit_IIV_CLcrrt_log*sqrt(T_IIV_CLcrrt_log) #-0.4302929

lower_bound_IIV_CLcrrt_log
upper_bound_IIV_CLcrrt_log

# back transform IIV_CLcrrt_log to IIV_CLcrrt
IIV_CLcrrt <- exp(IIV_CLcrrt_log) #0.2824382

lower_bound_IIV_CLcrrt <- exp(lower_bound_IIV_CLcrrt_log) #0.122665

upper_bound_IIV_CLcrrt <- exp(upper_bound_IIV_CLcrrt_log) #0.6503186

IIV_CLcrrt

lower_bound_IIV_CLcrrt

upper_bound_IIV_CLcrrt

### IIV_CLr -----------------

# Mean estimate
IIV_CLr_log <- mean(Pooling_Parameters_2l.pan$Log[Pooling_Parameters_2l.pan$Parameters == "IIV_CLr"]) #-1.526628

W_IIV_CLr_log <- mean(Pooling_Parameters_2l.pan$Var_Log[Pooling_Parameters_2l.pan$Parameters == "IIV_CLr"]) #0.04087116 #within imputation variance

B_IIV_CLr_log <- var(Pooling_Parameters_2l.pan$Log[Pooling_Parameters_2l.pan$Parameters == "IIV_CLr"]) #0.001377065 #between imputation variance

IIV_CLr_log
W_IIV_CLr_log
B_IIV_CLr_log

# Variance estimate
m <- 70
T_IIV_CLr_log <- W_IIV_CLr_log + (1 + 1/m)*B_IIV_CLr_log #0.04226789
T_IIV_CLr_log

# Confidence interval - T distribution
m <- 70 #number of imputations
k <- 12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect parameters where covariance is one of them)
n <- 177

lambda_IIV_CLr_log <- (B_IIV_CLr_log + B_IIV_CLr_log/m)/T_IIV_CLr_log #proportion of variation attributable to the missing data

r_IIV_CLr_log <- (lambda_IIV_CLr_log)/(1 - lambda_IIV_CLr_log) #relative increase in variance due to nonresponse

v_old_IIV_CLr_log <- (m - 1)*(1 + 1/r_IIV_CLr_log^2) #old degree of freedom

v_com = n - k #degree of freedom of parameter estimate (IIV_CLr_log) in the hypothetically complete data

v_obs_IIV_CLr_log = ((v_com + 1)/(v_com + 3))*v_com*(1 - lambda_IIV_CLr_log) #observed data degrees of freedom that accounts for the missing information

v_IIV_CLr_log <- (v_old_IIV_CLr_log*v_obs_IIV_CLr_log)/(v_old_IIV_CLr_log + v_obs_IIV_CLr_log)

t_crit_IIV_CLr_log <- qt(0.975, v_IIV_CLr_log) 

lower_bound_IIV_CLr_log <- IIV_CLr_log - t_crit_IIV_CLr_log*sqrt(T_IIV_CLr_log) #-1.932705
upper_bound_IIV_CLr_log <- IIV_CLr_log + t_crit_IIV_CLr_log*sqrt(T_IIV_CLr_log) #-1.12055

lower_bound_IIV_CLr_log
upper_bound_IIV_CLr_log

# back transform IIV_CLr_log to IIV_CLr
IIV_CLr <- exp(IIV_CLr_log) #0.2172671

lower_bound_IIV_CLr <- exp(lower_bound_IIV_CLr_log) #0.1447561

upper_bound_IIV_CLr <- exp(upper_bound_IIV_CLr_log) #0.3261004

IIV_CLr

lower_bound_IIV_CLr

upper_bound_IIV_CLr

### IIV_V1 -----------------

# Mean estimate
IIV_V1_log <- mean(Pooling_Parameters_2l.pan$Log[Pooling_Parameters_2l.pan$Parameters == "IIV_V1"]) #-1.164145

W_IIV_V1_log <- mean(Pooling_Parameters_2l.pan$Var_Log[Pooling_Parameters_2l.pan$Parameters == "IIV_V1"]) #0.04341324 #within imputation variance

B_IIV_V1_log <- var(Pooling_Parameters_2l.pan$Log[Pooling_Parameters_2l.pan$Parameters == "IIV_V1"]) #0.000160931 #between imputation variance

IIV_V1_log
W_IIV_V1_log
B_IIV_V1_log

# Variance estimate
m <- 70
T_IIV_V1_log <- W_IIV_V1_log + (1 + 1/m)*B_IIV_V1_log #0.04357647
T_IIV_V1_log

# Confidence interval - T distribution
m <- 70 #number of imputations
k <- 12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect parameters where covariance is one of them)
n <- 177

lambda_IIV_V1_log <- (B_IIV_V1_log + B_IIV_V1_log/m)/T_IIV_V1_log #proportion of variation attributable to the missing data

r_IIV_V1_log <- (lambda_IIV_V1_log)/(1 - lambda_IIV_V1_log) #relative increase in variance due to nonresponse

v_old_IIV_V1_log <- (m - 1)*(1 + 1/r_IIV_V1_log^2) #old degree of freedom

v_com = n - k #degree of freedom of parameter estimate (IIV_V1_log) in the hypothetically complete data

v_obs_IIV_V1_log = ((v_com + 1)/(v_com + 3))*v_com*(1 - lambda_IIV_V1_log) #observed data degrees of freedom that accounts for the missing information

v_IIV_V1_log <- (v_old_IIV_V1_log*v_obs_IIV_V1_log)/(v_old_IIV_V1_log + v_obs_IIV_V1_log)

t_crit_IIV_V1_log <- qt(0.975, v_IIV_V1_log) 

lower_bound_IIV_V1_log <- IIV_V1_log - t_crit_IIV_V1_log*sqrt(T_IIV_V1_log) #-1.576358
upper_bound_IIV_V1_log <- IIV_V1_log + t_crit_IIV_V1_log*sqrt(T_IIV_V1_log) #-0.7519313

lower_bound_IIV_V1_log
upper_bound_IIV_V1_log

# back transform IIV_CLr_log to IIV_CLr
IIV_V1 <- exp(IIV_V1_log) #0.3121895

lower_bound_IIV_V1 <- exp(lower_bound_IIV_V1_log) #0.2067266

upper_bound_IIV_V1 <- exp(upper_bound_IIV_V1_log) #0.4714551

IIV_V1

lower_bound_IIV_V1

upper_bound_IIV_V1

### PropErr -----------------

# Mean estimate
PropErr_log <- mean(Pooling_Parameters_2l.pan$Log[Pooling_Parameters_2l.pan$Parameters == "PropErr"]) #-3.588325

W_PropErr_log <- mean(Pooling_Parameters_2l.pan$Var_Log[Pooling_Parameters_2l.pan$Parameters == "PropErr"]) #0.01817211 #within imputation variance

B_PropErr_log <- var(Pooling_Parameters_2l.pan$Log[Pooling_Parameters_2l.pan$Parameters == "PropErr"]) #8.316011e-05
#between imputation variance

PropErr_log
W_PropErr_log
B_PropErr_log

# Variance estimate
m <- 70
T_PropErr_log <- W_PropErr_log + (1 + 1/m)*B_PropErr_log #0.01825646
T_PropErr_log

# Confidence interval - T distribution
m <- 70 #number of imputations
k <- 12 #my model has 12 parameters (including 6 fixed effects parameters and 6 random effect parameters where covariance is one of them)
n <- 177

lambda_PropErr_log <- (B_PropErr_log + B_PropErr_log/m)/T_PropErr_log #proportion of variation attributable to the missing data

r_PropErr_log <- (lambda_PropErr_log)/(1 - lambda_PropErr_log) #relative increase in variance due to nonresponse

v_old_PropErr_log <- (m - 1)*(1 + 1/r_PropErr_log^2) #old degree of freedom

v_com = n - k #degree of freedom of parameter estimate (PropErr_log) in the hypothetically complete data

v_obs_PropErr_log = ((v_com + 1 )/(v_com + 3))*v_com*(1 - lambda_PropErr_log) #observed data degrees of freedom that accounts for the missing information

v_PropErr_log <- (v_old_PropErr_log*v_obs_PropErr_log)/(v_old_PropErr_log + v_obs_PropErr_log)

t_crit_PropErr_log <- qt(0.975, v_PropErr_log) 

lower_bound_PropErr_log <- PropErr_log - t_crit_PropErr_log*sqrt(T_PropErr_log) #-3.855139
upper_bound_PropErr_log <- PropErr_log + t_crit_PropErr_log*sqrt(T_PropErr_log) #-3.321512

lower_bound_PropErr_log
upper_bound_PropErr_log

# back transform IIV_CLr_log to IIV_CLr
PropErr <- exp(PropErr_log) #0.02764459

lower_bound_PropErr <- exp(lower_bound_PropErr_log) #0.02117067

upper_bound_PropErr <- exp(upper_bound_PropErr_log) #0.03609821

PropErr

lower_bound_PropErr

upper_bound_PropErr

# Table 2 - Manucript ------------------------------------------

## Base model ------------------------------------------

# Calculate %RSE of the CV term of Cov_CLr_V1

RSE_CV_Cov_CLr_V1 <- 0.0454/0.111/2*100

## Pooled model ------------------------------------------

# Mean (95% CI) for CV term of IIV_CLcrrt

sqrt(exp(IIV_CLcrrt) - 1)*100 # 57.12791
sqrt(exp(lower_bound_IIV_CLcrrt) - 1)*100 # 36.12557
sqrt(exp(upper_bound_IIV_CLcrrt) - 1)*100 # 95.71579

# Mean (95% CI) for CV term of IIV_CLr

sqrt(exp(IIV_CLr) - 1)*100 # 49.26216
sqrt(exp(lower_bound_IIV_CLr) - 1)*100 # 39.46614
sqrt(exp(upper_bound_IIV_CLr) - 1)*100 # 62.09303

# Mean (95% CI) for CV term of IIV_V1

sqrt(exp(IIV_V1) - 1)*100 # 60.53211
sqrt(exp(lower_bound_IIV_V1) - 1)*100 # 47.92143
sqrt(exp(upper_bound_IIV_V1) - 1)*100 # 77.60954

# Mean (95% CI) for CV term of Cov_CLr_V1

sqrt(exp(Cov_CLr_V1) - 1)*100 # 26.39069
sqrt(exp(-lower_bound_Cov_CLr_V1) - 1)*100 # 9.365618
sqrt(exp(upper_bound_Cov_CLr_V1) - 1)*100 # 39.2658

# Mean (95% CI) for CV term of PropErr

sqrt(exp(PropErr) - 1)*100 # 16.74223
sqrt(exp(lower_bound_PropErr) - 1)*100 # 14.62749
sqrt(exp(upper_bound_PropErr) - 1)*100 # 19.17229

# Dose-finding simulations ------------------------------------------

## Firstly, the effect of CKDEPI on non-CRRT events --------------------------

### Some data manipulations ------------------------------------------

# load the dataset
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Bootstrap and dosing simulation/BW_Simulations/Simulation datasets")
sim_dosing_01 <- read.csv("dosing_08.csv")

# Eliminate all rows with PEAK == 1
sim_dosing_01 <- sim_dosing_01 |>  
                  dplyr::filter(PEAK == 0)

# Eliminate PEAK column
sim_dosing_01$PEAK <- NULL

# Convert CRRT into 0
sim_dosing_01$CRRT <- 0

# summary statistics for CKDEPI in non-CRRT events
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Original")
Fluco_clean_revised <- read.csv("Fluco_clean_revised.csv", na.strings = "NA")
summary(Fluco_clean_revised$CKDEPI[Fluco_clean_revised$CRRT == 0])

# Min 5 - max 215 - so I will try firstly distance = 10, and then = 5.
# Note: Normal kidney function is from 90 to 120 mL/min/1.73 m2

# Count the number of observed CKDEPI in non-CRRT events
Fluco_clean_revised |> 
  dplyr::filter(CRRT == 0) |>
  select(CKDEPI) |> 
  na.omit() |>
  summarize(n = n()) # 485

# Count the number of observed CKDEPI >= 195 in non-CRRT events
Fluco_clean_revised |> 
  dplyr::filter(CRRT == 0 & CKDEPI >= 195) |>
  select(CKDEPI) |> 
  na.omit() |>
  summarize(n = n()) # 8

# The percentage of CKDEPI >= 195 in non-CRRT events - Erwin's question
8/485 # 1.6%

# Count the number of observed CKDEPI <= 25 in non-CRRT events
Fluco_clean_revised |> 
  dplyr::filter(CRRT == 0 & CKDEPI <= 25) |>
  select(CKDEPI) |> 
  na.omit() |>
  summarize(n = n()) # 23 

# The percentage of CKDEPI <= 25 in non-CRRT events
23/485 # 4.7%
# Based on the analysis that follows, these are from 13 patients, only 2 of which
# have CRRT after having CKDEPI <= 25.

# Count the number of observed CKDEPI <= 15 in non-CRRT events
Fluco_clean_revised |> 
  dplyr::filter(CRRT == 0 & CKDEPI <= 15) |>
  select(CKDEPI) |> 
  na.omit() |>
  summarize(n = n()) # 5 - which is only 1.03%

## Exploratory code from chat GPT

# First, identify patient IDs with at least one CKDEPI <= 25 and CRRT == 0
ID_CKDEPI25 <- Fluco_clean_revised |>
  filter(CKDEPI <= 25 & CRRT == 0) |> 
  pull(ID) |>
  unique()

# Next, identify patients having CRRT among those IDs above
CKDEPI_25 <- Fluco_clean_revised |>
  filter(ID %in% ID_CKDEPI25) |>
  filter(CRRT == 1) |>
  pull(ID) |>
  unique()

CKDEPI_25_CRRT <- Fluco_clean_revised |>
  filter(ID %in% CKDEPI_25)
## End of exploratory code from chat GPT

# Creat a new column "CKDEPI" for sim_dosing_01 dataset
sim_dosing_01 <- sim_dosing_01 |> 
  mutate(CKDEPI = ifelse(ID == 1, 5, 5 + (ID - 1) * 10))

# Delete all the rows with CKDEPI > 215
sim_dosing_01 <- sim_dosing_01 |> 
  dplyr::filter(CKDEPI <= 215)

# Change ID into ID_CKDEPI
sim_dosing_01 <- sim_dosing_01 |> 
  rename(ID_CKDEPI = ID)

# Replicate each ID_CKDEPI 100 times
sim_dosing_01 <- sim_dosing_01 |> 
  group_by(ID_CKDEPI) |> 
  slice(rep(1:n(), each = 100)) |> 
  ungroup()

# Assign an ID for each replication
sim_dosing_01 <- sim_dosing_01 |> 
  arrange(TIME) |> 
  group_by(TIME) |> 
  mutate(ID = row_number()) |> 
  ungroup()

##  Assign a unique BW per each ID

# Step 1: Get unique, non-missing values of BW2 from Fluco_clean_revised 
unique_BW2 <- Fluco_clean_revised  |> 
  distinct(BW2, .keep_all = TRUE) |> 
  dplyr::filter(!is.na(BW2)) |> 
  pull(BW2)

# Step 1 - revised 020424: Get unique, non-missing values of BW2 per ID level
unique_BW2 <- Fluco_clean_revised  |>
  group_by(ID) |>
  distinct(BW2, .keep_all = TRUE) |>
  dplyr::filter(!is.na(BW2)) |>
  pull(BW2)

# Create a new dataset "new_dataset" with ID from 1 to 2200
set.seed(123)
new_dataset <- data.frame(ID = 1:2200)

# Step 2: Create an ECDF based on unique BW2 values
ecdf_BW2 <- ecdf(unique_BW2)

# Create a vector to store the sampled BWs
sampled_BWs <- numeric(0)

# Loop 22 times to get 22 sets of 100 random values
for (i in 1:22) {
  # Generate random values from the ECDF
  random_values <- runif(100)
  
  # Map the random values to the empirical distribution using the inverse transform sampling method
  sampled_BW <- round(quantile(unique_BW2, random_values), 1)
  
  # Append the sampled BWs to the vector
  sampled_BWs <- c(sampled_BWs, sampled_BW)
}

# Assuming new_dataset is your data frame with IDs
new_dataset$BW <- sampled_BWs

# Now, you have "new_dataset" with unique BW values for each set of 100 IDs

# Merge "new_dataset" with "sim_dosing_01" dataset based on ID
sim_dosing_01 <- sim_dosing_01 |> 
  left_join(new_dataset, by = "ID")

# Rename BW
sim_dosing_01$BW.x <- NULL
sim_dosing_01 <- sim_dosing_01 |> 
  rename(BW = BW.y)

# Eliminate ID_CKDEPI, change the name of the dataset
sim_CKDEPI_01 <- sim_dosing_01 |> 
  select(ID_CKDEPI, ID, everything())
sim_CKDEPI_01$ID_CKDEPI <- NULL

# Sort sim_CKDEPI_01 according to ID
sim_CKDEPI_01 <- sim_CKDEPI_01 |> 
  arrange(ID)

# Here I have another _min dataset, where all dosing events need to be kept, and eliminate all redundant dummy lines
sim_CKDEPI_01_min <- sim_CKDEPI_01 |> 
  dplyr::filter((DAY %in% c(3, 4, 5, 8, 9, 10, 11, 12) & TAD == 0) | (!DAY %in% c(3, 4,5, 8, 9, 10, 11, 12)))

### Standard dosing regimen ------------------------------------------

# Change the dosing follows standard regimen 800 LD & 400 MD
sim_CKDEPI_01_min$AMT <- ifelse(sim_CKDEPI_01_min$TIME == 0,800,ifelse(sim_CKDEPI_01_min$TAD == 0,400,0))

# Then export into sim_CKDEPI_01.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/CKDEPI_Sim")
write.csv(sim_CKDEPI_01_min, "sim_CKDEPI_01.csv", quote = F, row.names = F)

### 1000 mg LD & 400 mg MD ------------------------------------------

# Change the dosing regimen
sim_CKDEPI_01_min$AMT <- ifelse(sim_CKDEPI_01_min$TIME == 0, 1000, ifelse(sim_CKDEPI_01_min$TAD == 0, 400, 0))

# Then export into sim_CKDEPI_02.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/CKDEPI_Sim")
write.csv(sim_CKDEPI_01_min, "sim_CKDEPI_02.csv", quote = F, row.names = F)

### 1200 mg LD & 400 mg MD ------------------------------------------

# Change the dosing regimen
sim_CKDEPI_01_min$AMT <- ifelse(sim_CKDEPI_01_min$TIME == 0, 1200, ifelse(sim_CKDEPI_01_min$TAD == 0, 400, 0))

# Then export into sim_CKDEPI_03.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/CKDEPI_Sim")
write.csv(sim_CKDEPI_01_min, "sim_CKDEPI_03.csv", quote = F, row.names = F)

### 12 mg/kg LD & 6 mg/kg MD ------------------------------------------

# Change the dosing regimen
sim_CKDEPI_01_min$AMT <- ifelse(sim_CKDEPI_01_min$TIME == 0, 12 * sim_CKDEPI_01_min$BW, ifelse(sim_CKDEPI_01_min$TAD == 0, 6 * sim_CKDEPI_01_min$BW, 0))

# Then export into sim_CKDEPI_04.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/CKDEPI_Sim")
write.csv(sim_CKDEPI_01_min, "sim_CKDEPI_04.csv", quote = F, row.names = F)

### 1200 mg LD & 600 mg MD ------------------------------------------

# Change the dosing regimen
sim_CKDEPI_01_min$AMT <- ifelse(sim_CKDEPI_01_min$TIME == 0, 1200, ifelse(sim_CKDEPI_01_min$TAD == 0, 600, 0))

# Then export into sim_CKDEPI_05.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/CKDEPI_Sim")
write.csv(sim_CKDEPI_01_min, "sim_CKDEPI_05.csv", quote = F, row.names = F)

### 800 mg LD & 200 mg MD ------------------------------------------

# Change the dosing regimen
sim_CKDEPI_01_min$AMT <- ifelse(sim_CKDEPI_01_min$TIME == 0, 800, ifelse(sim_CKDEPI_01_min$TAD == 0, 200, 0))

# Then export into sim_CKDEPI_06.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/CKDEPI_Sim")
write.csv(sim_CKDEPI_01_min, "sim_CKDEPI_06.csv", quote = F, row.names = F)

## Secondly, the effect of BW on non-CRRT events --------------------------

### Some data manipulations ------------------------------------------

# load the dataset
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Bootstrap and dosing simulation/BW_Simulations/Simulation datasets")
sim_dosing_01 <- read.csv("dosing_08.csv")

# Eliminate all rows with PEAK == 1
sim_dosing_01 <- sim_dosing_01 |> dplyr::filter(PEAK == 0)

# Eliminate PEAK column
sim_dosing_01$PEAK <- NULL

# Convert CRRT into 0
sim_dosing_01$CRRT <- 0

# summary statistics for CKDEPI in non-CRRT events
summary(Fluco_clean_revised$BW2) # range: 34 - 142 kg. However Erwin said 30 - 150

# Change ID into ID_BW
sim_dosing_01 <- sim_dosing_01 |>
  rename(ID_BW = ID)

# Replicate each ID_CKDEPI 100 times
sim_dosing_01 <- sim_dosing_01 |>
  group_by(ID_BW) |>
  slice(rep(1:n(), each = 100)) |>
  ungroup()

# Assign an ID for each replication
sim_dosing_01 <- sim_dosing_01 |>
  arrange(TIME) |>
  group_by(TIME) |>
  mutate(ID = row_number()) |>
  ungroup()

##  Assign a unique CKDEPI per each ID

# Step 1: Get unique, non-missing values of CKDEPI from FlucTotIV_clean with CRRT = 0
unique_CKDEPI <- FlucTotIV_clean |> dplyr::filter(CRRT == 0) |>
  dplyr::filter(!is.na(CKDEPI)) |>
  pull(CKDEPI)

# Create a new dataset "new_dataset" with ID from 1 to 2500
set.seed(123)
new_dataset <- data.frame(ID = 1:2500)

# Step 2: Create an ECDF based on unique CKDEPI values
ecdf_CKDEPI <- ecdf(unique_CKDEPI)

# Create a vector to store the sampled CKDEPIs
sampled_CKDEPIs <- numeric(0)

# Loop 25 times to get 25 sets of 100 random values
for (i in 1:25) {
  # Generate random values from the ECDF
  random_values <- runif(100)
  
  # Map the random values to the empirical distribution using the inverse transform sampling method
  sampled_CKDEPI <- round(quantile(unique_CKDEPI, random_values), 1)
  
  # Append the sampled CKDEPIs to the vector
  sampled_CKDEPIs <- c(sampled_CKDEPIs, sampled_CKDEPI)
}

# Assuming new_dataset is your data frame with IDs
new_dataset$CKDEPI <- sampled_CKDEPIs

# Now, you have "new_dataset" with unique BW values for each set of 100 IDs

# Merge "new_dataset" with "sim_dosing_01" dataset based on ID
sim_dosing_01 <- sim_dosing_01 |>
  left_join(new_dataset, by = "ID")

# Eliminate ID_BW, change the name of the dataset
sim_BW_01 <- sim_dosing_01 |>
  select(ID_BW, ID, everything())
sim_BW_01$ID_BW <- NULL

# Sort sim_BW_01 according to ID
sim_BW_01 <- sim_BW_01 |>
  arrange(ID)

# Minimizing the dataset where all dosing events are kept, and all redundant dummy lines can be eliminated
sim_BW_01_min <- sim_BW_01 |>
  dplyr::filter((DAY %in% c(3, 4, 5, 8, 9, 10, 11, 12) & TAD == 0) | (!DAY %in% c(3, 4,5, 8, 9, 10, 11, 12)))

### Standard dosing regimen ------------------------------------------

# Change the dosing follows standard regimen 800 LD & 400 MD
sim_BW_01_min$AMT <- ifelse(sim_BW_01_min$TIME == 0,800,ifelse(sim_BW_01_min$TAD == 0,400,0))

# Then export into sim_BW_noCRRT_01.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/BW_noCRRT")
write.csv(sim_BW_01_min, "sim_BW_noCRRT_01.csv", quote = F, row.names = F)

### 1000 mg LD & 400 mg MD ------------------------------------------

# Change the dosing regimen
sim_BW_01_min$AMT <- ifelse(sim_BW_01_min$TIME == 0, 1000, ifelse(sim_BW_01_min$TAD == 0, 400, 0))

# Then export into sim_BW_noCRRT_02.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/BW_noCRRT")
write.csv(sim_BW_01_min, "sim_BW_noCRRT_02.csv", quote = F, row.names = F)

### 1200 mg LD & 400 mg MD ------------------------------------------

# Change the dosing regimen
sim_BW_01_min$AMT <- ifelse(sim_BW_01_min$TIME == 0, 1200, ifelse(sim_BW_01_min$TAD == 0, 400, 0))

# Then export into sim_BW_noCRRT_03.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/BW_noCRRT")
write.csv(sim_BW_01_min, "sim_BW_noCRRT_03.csv", quote = F, row.names = F)

### 1400 mg LD & 400 mg MD ------------------------------------------

# Change the dosing regimen
sim_BW_01_min$AMT <- ifelse(sim_BW_01_min$TIME == 0, 1400, ifelse(sim_BW_01_min$TAD == 0, 400, 0))

# Then export into sim_BW_noCRRT_04.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/BW_noCRRT")
write.csv(sim_BW_01_min, "sim_BW_noCRRT_04.csv", quote = F, row.names = F)

### 1600 mg LD & 400 mg MD ------------------------------------------

# Change the dosing regimen
sim_BW_01_min$AMT <- ifelse(sim_BW_01_min$TIME == 0, 1600, ifelse(sim_BW_01_min$TAD == 0, 400, 0))

# Then export into sim_BW_noCRRT_05.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/BW_noCRRT")
write.csv(sim_BW_01_min, "sim_BW_noCRRT_05.csv", quote = F, row.names = F)

### 1800 mg LD & 400 mg MD ------------------------------------------

# Change the dosing regimen
sim_BW_01_min$AMT <- ifelse(sim_BW_01_min$TIME == 0, 1800, ifelse(sim_BW_01_min$TAD == 0, 400, 0))

# Then export into sim_BW_noCRRT_06.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/BW_noCRRT")
write.csv(sim_BW_01_min, "sim_BW_noCRRT_06.csv", quote = F, row.names = F)

### 12 mg/kg LD & 6 mg/kg MD ------------------------------------------

# Change the dosing regimen
sim_BW_01_min$AMT <- ifelse(sim_BW_01_min$TIME == 0, 12 * sim_BW_01_min$BW, ifelse(sim_BW_01_min$TAD == 0, 6 * sim_BW_01_min$BW, 0))

# Then export into sim_BW_noCRRT_07.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/BW_noCRRT")
write.csv(sim_BW_01_min, "sim_BW_noCRRT_07.csv", quote = F, row.names = F)

### 800 mg LD & 200 mg MD ------------------------------------------

# Change the dosing regimen
sim_BW_01_min$AMT <- ifelse(sim_BW_01_min$TIME == 0, 800, ifelse(sim_BW_01_min$TAD == 0, 200, 0))

# Then export into sim_BW_noCRRT_08.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/BW_noCRRT")
write.csv(sim_BW_01_min, "sim_BW_noCRRT_08.csv", quote = F, row.names = F)

## Third, the effect of BW on CRRT events --------------------------

### Some data manipulations ------------------------------------------

# load the dataset
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Bootstrap and dosing simulation/BW_Simulations/Simulation datasets")
sim_dosing_01 <- read.csv("dosing_08.csv")

# Eliminate all rows with PEAK == 1
sim_dosing_01 <- sim_dosing_01 |> dplyr::filter(PEAK == 0)

# Eliminate PEAK column
sim_dosing_01$PEAK <- NULL

# BW CRRT sim dataset
sim_BW_CRRT_01 <- sim_dosing_01 |>
  dplyr::filter((DAY %in% c(3, 4, 5, 8, 9, 10, 11, 12) & TAD == 0) | (!DAY %in% c(3, 4,5, 8, 9, 10, 11, 12)))

sim_BW_CRRT_01$CKDEPI <- 0

### Standard dosing regimen ------------------------------------------

# Change the dosing follows standard regimen 800 LD & 400 MD
sim_BW_CRRT_01$AMT <- ifelse(sim_BW_CRRT_01$TIME == 0, 800, ifelse(sim_BW_CRRT_01$TAD == 0, 400, 0))

# Then export into sim_BW_CRRT_01.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/BW_CRRT")
write.csv(sim_BW_CRRT_01, "sim_BW_CRRT_01.csv", quote = F, row.names = F)

### 1000 mg LD & 800 mg MD ------------------------------------------

# Change the dosing regimen
sim_BW_CRRT_01$AMT <- ifelse(sim_BW_CRRT_01$TIME == 0, 1000, ifelse(sim_BW_CRRT_01$TAD == 0, 800, 0))

# Then export into sim_BW_CRRT_02.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/BW_CRRT")
write.csv(sim_BW_CRRT_01, "sim_BW_CRRT_02.csv", quote = F, row.names = F)

### 1200 mg LD & 800 mg MD ------------------------------------------

# Change the dosing regimen
sim_BW_CRRT_01$AMT <- ifelse(sim_BW_CRRT_01$TIME == 0, 1200, ifelse(sim_BW_CRRT_01$TAD == 0, 800, 0))

# Then export into sim_BW_CRRT_03.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/BW_CRRT")
write.csv(sim_BW_CRRT_01, "sim_BW_CRRT_03.csv", quote = F, row.names = F)

### 1400 mg LD & 800 mg MD ------------------------------------------

# Change the dosing regimen
sim_BW_CRRT_01$AMT <- ifelse(sim_BW_CRRT_01$TIME == 0, 1400, ifelse(sim_BW_CRRT_01$TAD == 0, 800, 0))

# Then export into sim_BW_CRRT_04.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/BW_CRRT")
write.csv(sim_BW_CRRT_01, "sim_BW_CRRT_04.csv", quote = F, row.names = F)

### 1600 mg LD & 800 mg MD ------------------------------------------

# Change the dosing regimen
sim_BW_CRRT_01$AMT <- ifelse(sim_BW_CRRT_01$TIME == 0, 1600, ifelse(sim_BW_CRRT_01$TAD == 0, 800, 0))

# Then export into sim_BW_CRRT_03.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/BW_CRRT")
write.csv(sim_BW_CRRT_01, "sim_BW_CRRT_05.csv", quote = F, row.names = F)

### 1800 mg LD & 800 mg MD ------------------------------------------

# Change the dosing regimen
sim_BW_CRRT_01$AMT <- ifelse(sim_BW_CRRT_01$TIME == 0, 1800, ifelse(sim_BW_CRRT_01$TAD == 0, 800, 0))

# Then export into sim_BW_CRRT_03.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/BW_CRRT")
write.csv(sim_BW_CRRT_01, "sim_BW_CRRT_06.csv", quote = F, row.names = F)

### 12 mg/kg LD & 6 mg/kg MD ------------------------------------------

# Change the dosing regimen
sim_BW_CRRT_01$AMT <- ifelse(sim_BW_CRRT_01$TIME == 0, 12 * sim_BW_CRRT_01$BW, ifelse(sim_BW_CRRT_01$TAD == 0, 6 * sim_BW_CRRT_01$BW, 0))

# Then export into sim_BW_CRRT_04.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/BW_CRRT")
write.csv(sim_BW_CRRT_01, "sim_BW_CRRT_07.csv", quote = F, row.names = F)

### 800 mg LD & 600 mg MD ------------------------------------------

# Change the dosing regimen
sim_BW_CRRT_01$AMT <- ifelse(sim_BW_CRRT_01$TIME == 0, 800, ifelse(sim_BW_CRRT_01$TAD == 0, 600, 0))

# Then export into sim_BW_CRRT_03.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/BW_CRRT")
write.csv(sim_BW_CRRT_01, "sim_BW_CRRT_08.csv", quote = F, row.names = F)

## Fourth, verify the BW-based dosing on CKDEPI range --------------------------

### Some data manipulations ------------------------------------------

# load the dataset
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Bootstrap and dosing simulation/BW_Simulations/Simulation datasets")
core_sim <- read.csv("dosing_08.csv") # core_sim stands for core simulation dataset

# Eliminate all rows with PEAK == 1
core_sim <- core_sim |>  
  dplyr::filter(PEAK == 0)

# Eliminate PEAK column
core_sim$PEAK <- NULL

# Convert CRRT into 0
core_sim$CRRT <- 0

# Creat a new column "CKDEPI" for core_sim dataset
core_sim <- core_sim |> 
  mutate(CKDEPI = ifelse(ID == 1, 5, 5 + (ID - 1) * 10))

# Delete all the rows with CKDEPI > 215
core_sim <- core_sim |> 
  dplyr::filter(CKDEPI <= 215)

# Change ID into ID_CKDEPI
core_sim <- core_sim |> 
  rename(ID_CKDEPI = ID)

# Hereon I will create 5 different datasets with different dosing regimens, 
# each will be assigned 20 different BWs. 

### 800 mg LD, 400 MD for BW <= 60kg -------------------------------

# Replicate each ID_CKDEPI 20 times
LD_800 <- core_sim |> 
  group_by(ID_CKDEPI) |> 
  slice(rep(1:n(), each = 20)) |> 
  ungroup()

# Assign an ID for each replication
LD_800 <- LD_800 |> 
  arrange(TIME) |> 
  group_by(TIME) |> 
  mutate(ID = row_number()) |> 
  ungroup()

##  Assign a unique BW per each ID

# Extract only BW values that are <= 60 kg
BW_LD_800 <- Fluco_clean_revised |>
  group_by(ID) |>
  distinct(BW2, .keep_all = TRUE) |>
  dplyr::filter(!is.na(BW2)) |>
  dplyr::filter(BW2 <= 60) |>
  pull(BW2)

# Create a new dataset "LD_800_BW" with ID from 1 to 440
set.seed(123)
LD_800_BW <- data.frame(ID = 1:440)

# Step 2: Create an ECDF based on unique BW2 values
ecdf_LD_800_BW <- ecdf(BW_LD_800)

# Create a vector to store the sampled BWs
sampled_LD_800_BWs <- numeric(0)

# Loop 22 times to get 22 sets of 20 random values
for (i in 1:22) {
  # Generate random values from the ECDF
  random_values <- runif(20)
  
  # Map the random values to the empirical distribution using the inverse transform sampling method
  sampled_LD_800_BW <- round(quantile(BW_LD_800, random_values), 1)
  
  # Append the sampled BWs to the vector
  sampled_LD_800_BWs <- c(sampled_LD_800_BWs, sampled_LD_800_BW)
}

# Assuming new_dataset is your data frame with IDs
LD_800_BW$BW <- sampled_LD_800_BWs

# Now, you have "LD_800_BW" with unique BW values for each set of 20 IDs

# Merge "new_dataset" with "core_sim" dataset based on ID
LD_800 <- LD_800 |> 
  left_join(LD_800_BW, by = "ID")

# Rename BW
LD_800$BW.x <- NULL
LD_800 <- LD_800 |> 
  rename(BW = BW.y)

# Re-arrange the dataset
LD_800_CKDEPI <- LD_800 |> 
  select(ID_CKDEPI, ID, everything())

LD_800_CKDEPI <- LD_800_CKDEPI |> 
  arrange(ID)

# Add regimen column
LD_800_CKDEPI$Regimen <- 1

# Here I keep all dosing events, and eliminate all redundant dummy lines
LD_800_CKDEPI <- LD_800_CKDEPI |> 
  dplyr::filter((DAY %in% c(3, 4, 5, 8, 9, 10, 11, 12) & TAD == 0) | (!DAY %in% c(3, 4,5, 8, 9, 10, 11, 12)))

# Change the dosing follows 800 LD & 400 MD
LD_800_CKDEPI$AMT <- ifelse(LD_800_CKDEPI$TIME == 0, 800, ifelse(LD_800_CKDEPI$TAD == 0, 400, 0))

### 1000 mg LD, 400 MD for 60 < BW <= 80kg -------------------------------

# Replicate each ID_CKDEPI 20 times
LD_1000 <- core_sim |> 
  group_by(ID_CKDEPI) |> 
  slice(rep(1:n(), each = 20)) |> 
  ungroup()

# Assign an ID for each replication
LD_1000 <- LD_1000 |> 
  arrange(TIME) |> 
  group_by(TIME) |> 
  mutate(ID = row_number()) |> 
  ungroup()

##  Assign a unique BW per each ID

# Extract only 60 < BW <= 80 values
BW_LD_1000 <- Fluco_clean_revised |>
  group_by(ID) |>
  distinct(BW2, .keep_all = TRUE) |>
  dplyr::filter(!is.na(BW2)) |>
  dplyr::filter(60 < BW2 & BW2 <= 80) |>
  pull(BW2)

# Create a new dataset "LD_1000_BW" with ID from 1 to 440
set.seed(123)
LD_1000_BW <- data.frame(ID = 1:440)

# Step 2: Create an ECDF based on unique BW2 values
ecdf_LD_1000_BW <- ecdf(BW_LD_1000)

# Create a vector to store the sampled BWs
sampled_LD_1000_BWs <- numeric(0)

# Loop 22 times to get 22 sets of 20 random values
for (i in 1:22) {
  # Generate random values from the ECDF
  random_values <- runif(20)
  
  # Map the random values to the empirical distribution using the inverse transform sampling method
  sampled_LD_1000_BW <- round(quantile(BW_LD_1000, random_values), 1)
  
  # Append the sampled BWs to the vector
  sampled_LD_1000_BWs <- c(sampled_LD_1000_BWs, sampled_LD_1000_BW)
}

# Assuming new_dataset is your data frame with IDs
LD_1000_BW$BW <- sampled_LD_1000_BWs

# Now, you have "LD_1000_BW" with unique BW values for each set of 20 IDs

# Merge "new_dataset" with "core_sim" dataset based on ID
LD_1000 <- LD_1000 |> 
  left_join(LD_1000_BW, by = "ID")

# Rename BW
LD_1000$BW.x <- NULL
LD_1000 <- LD_1000 |> 
  rename(BW = BW.y)

# Re-arrange the dataset
LD_1000_CKDEPI <- LD_1000 |> 
  select(ID_CKDEPI, ID, everything())

LD_1000_CKDEPI <- LD_1000_CKDEPI |> 
  arrange(ID)

# Add regimen column
LD_1000_CKDEPI$Regimen <- 2

# Here I keep all dosing events, and eliminate all redundant dummy lines
LD_1000_CKDEPI <- LD_1000_CKDEPI |> 
  dplyr::filter((DAY %in% c(3, 4, 5, 8, 9, 10, 11, 12) & TAD == 0) | (!DAY %in% c(3, 4,5, 8, 9, 10, 11, 12)))

# Change the dosing follows 1000 LD & 400 MD
LD_1000_CKDEPI$AMT <- ifelse(LD_1000_CKDEPI$TIME == 0, 1000, ifelse(LD_1000_CKDEPI$TAD == 0, 400, 0))

### 1200 mg LD, 400 MD for 80 < BW <= 100kg -------------------------------

# Replicate each ID_CKDEPI 20 times
LD_1200 <- core_sim |> 
  group_by(ID_CKDEPI) |> 
  slice(rep(1:n(), each = 20)) |> 
  ungroup()

# Assign an ID for each replication
LD_1200 <- LD_1200 |> 
  arrange(TIME) |> 
  group_by(TIME) |> 
  mutate(ID = row_number()) |> 
  ungroup()

##  Assign a unique BW per each ID

# Extract only 80 < BW <= 100 values
BW_LD_1200 <- Fluco_clean_revised |>
  group_by(ID) |>
  distinct(BW2, .keep_all = TRUE) |>
  dplyr::filter(!is.na(BW2)) |>
  dplyr::filter(80 < BW2 & BW2 <= 100) |>
  pull(BW2)

# Create a new dataset "LD_1200_BW" with ID from 1 to 440
set.seed(123)
LD_1200_BW <- data.frame(ID = 1:440)

# Step 2: Create an ECDF based on unique BW2 values
ecdf_LD_1200_BW <- ecdf(BW_LD_1200)

# Create a vector to store the sampled BWs
sampled_LD_1200_BWs <- numeric(0)

# Loop 22 times to get 22 sets of 20 random values
for (i in 1:22) {
  # Generate random values from the ECDF
  random_values <- runif(20)
  
  # Map the random values to the empirical distribution using the inverse transform sampling method
  sampled_LD_1200_BW <- round(quantile(BW_LD_1200, random_values), 1)
  
  # Append the sampled BWs to the vector
  sampled_LD_1200_BWs <- c(sampled_LD_1200_BWs, sampled_LD_1200_BW)
}

# Assuming new_dataset is your data frame with IDs
LD_1200_BW$BW <- sampled_LD_1200_BWs

# Now, you have "LD_1200_BW" with unique BW values for each set of 20 IDs

# Merge "new_dataset" with "core_sim" dataset based on ID
LD_1200 <- LD_1200 |> 
  left_join(LD_1200_BW, by = "ID")

# Rename BW
LD_1200$BW.x <- NULL
LD_1200 <- LD_1200 |> 
  rename(BW = BW.y)

# Re-arrange the dataset
LD_1200_CKDEPI <- LD_1200 |> 
  select(ID_CKDEPI, ID, everything())

LD_1200_CKDEPI <- LD_1200_CKDEPI |> 
  arrange(ID)

# Add regimen column
LD_1200_CKDEPI$Regimen <- 3

# Here I keep all dosing events, and eliminate all redundant dummy lines
LD_1200_CKDEPI <- LD_1200_CKDEPI |> 
  dplyr::filter((DAY %in% c(3, 4, 5, 8, 9, 10, 11, 12) & TAD == 0) | (!DAY %in% c(3, 4,5, 8, 9, 10, 11, 12)))

# Change the dosing follows 1200 LD & 400 MD
LD_1200_CKDEPI$AMT <- ifelse(LD_1200_CKDEPI$TIME == 0, 1200, ifelse(LD_1200_CKDEPI$TAD == 0, 400, 0))

### 1400 mg LD, 400 MD for 100 < BW <= 120kg -------------------------------

# Replicate each ID_CKDEPI 20 times
LD_1400 <- core_sim |> 
  group_by(ID_CKDEPI) |> 
  slice(rep(1:n(), each = 20)) |> 
  ungroup()

# Assign an ID for each replication
LD_1400 <- LD_1400 |> 
  arrange(TIME) |> 
  group_by(TIME) |> 
  mutate(ID = row_number()) |> 
  ungroup()

##  Assign a unique BW per each ID

# Extract only 100 < BW <= 120 values
BW_LD_1400 <- Fluco_clean_revised |>
  group_by(ID) |>
  distinct(BW2, .keep_all = TRUE) |>
  dplyr::filter(!is.na(BW2)) |>
  dplyr::filter(100 < BW2 & BW2 <= 120) |>
  pull(BW2)

# Create a new dataset "LD_1400_BW" with ID from 1 to 440
set.seed(123)
LD_1400_BW <- data.frame(ID = 1:440)

# Step 2: Create an ECDF based on unique BW2 values
ecdf_LD_1400_BW <- ecdf(BW_LD_1400)

# Create a vector to store the sampled BWs
sampled_LD_1400_BWs <- numeric(0)

# Loop 22 times to get 22 sets of 20 random values
for (i in 1:22) {
  # Generate random values from the ECDF
  random_values <- runif(20)
  
  # Map the random values to the empirical distribution using the inverse transform sampling method
  sampled_LD_1400_BW <- round(quantile(BW_LD_1400, random_values), 1)
  
  # Append the sampled BWs to the vector
  sampled_LD_1400_BWs <- c(sampled_LD_1400_BWs, sampled_LD_1400_BW)
}

# Assuming new_dataset is your data frame with IDs
LD_1400_BW$BW <- sampled_LD_1400_BWs

# Now, you have "LD_1400_BW" with unique BW values for each set of 20 IDs

# Merge "new_dataset" with "core_sim" dataset based on ID
LD_1400 <- LD_1400 |> 
  left_join(LD_1400_BW, by = "ID")

# Rename BW
LD_1400$BW.x <- NULL
LD_1400 <- LD_1400 |> 
  rename(BW = BW.y)

# Re-arrange the dataset
LD_1400_CKDEPI <- LD_1400 |> 
  select(ID_CKDEPI, ID, everything())

LD_1400_CKDEPI <- LD_1400_CKDEPI |> 
  arrange(ID)

# Add regimen column
LD_1400_CKDEPI$Regimen <- 4

# Here I keep all dosing events, and eliminate all redundant dummy lines
LD_1400_CKDEPI <- LD_1400_CKDEPI |> 
  dplyr::filter((DAY %in% c(3, 4, 5, 8, 9, 10, 11, 12) & TAD == 0) | (!DAY %in% c(3, 4,5, 8, 9, 10, 11, 12)))

# Change the dosing follows 1400 LD & 400 MD
LD_1400_CKDEPI$AMT <- ifelse(LD_1400_CKDEPI$TIME == 0, 1400, ifelse(LD_1400_CKDEPI$TAD == 0, 400, 0))

### 1600 mg LD, 400 MD for 120kg < BW  -------------------------------

# Replicate each ID_CKDEPI 20 times
LD_1600 <- core_sim |> 
  group_by(ID_CKDEPI) |> 
  slice(rep(1:n(), each = 20)) |> 
  ungroup()

# Assign an ID for each replication
LD_1600 <- LD_1600 |> 
  arrange(TIME) |> 
  group_by(TIME) |> 
  mutate(ID = row_number()) |> 
  ungroup()

##  Assign a unique BW per each ID

# Extract only 120 < BW values
BW_LD_1600 <- Fluco_clean_revised |>
  group_by(ID) |>
  distinct(BW2, .keep_all = TRUE) |>
  dplyr::filter(!is.na(BW2)) |>
  dplyr::filter(120 < BW2) |>
  pull(BW2)

# Create a new dataset "LD_1600_BW" with ID from 1 to 440
set.seed(123)
LD_1600_BW <- data.frame(ID = 1:440)

# Step 2: Create an ECDF based on unique BW2 values
ecdf_LD_1600_BW <- ecdf(BW_LD_1600)

# Create a vector to store the sampled BWs
sampled_LD_1600_BWs <- numeric(0)

# Loop 22 times to get 22 sets of 20 random values
for (i in 1:22) {
  # Generate random values from the ECDF
  random_values <- runif(20)
  
  # Map the random values to the empirical distribution using the inverse transform sampling method
  sampled_LD_1600_BW <- round(quantile(BW_LD_1600, random_values), 1)
  
  # Append the sampled BWs to the vector
  sampled_LD_1600_BWs <- c(sampled_LD_1600_BWs, sampled_LD_1600_BW)
}

# Assuming new_dataset is your data frame with IDs
LD_1600_BW$BW <- sampled_LD_1600_BWs

# Now, you have "LD_1600_BW" with unique BW values for each set of 20 IDs

# Merge "new_dataset" with "core_sim" dataset based on ID
LD_1600 <- LD_1600 |> 
  left_join(LD_1600_BW, by = "ID")

# Rename BW
LD_1600$BW.x <- NULL
LD_1600 <- LD_1600 |> 
  rename(BW = BW.y)

# Re-arrange the dataset
LD_1600_CKDEPI <- LD_1600 |> 
  select(ID_CKDEPI, ID, everything())

LD_1600_CKDEPI <- LD_1600_CKDEPI |> 
  arrange(ID)

# Add regimen column
LD_1600_CKDEPI$Regimen <- 5

# Here I keep all dosing events, and eliminate all redundant dummy lines
LD_1600_CKDEPI <- LD_1600_CKDEPI |> 
  dplyr::filter((DAY %in% c(3, 4, 5, 8, 9, 10, 11, 12) & TAD == 0) | (!DAY %in% c(3, 4,5, 8, 9, 10, 11, 12)))

# Change the dosing follows 1600 LD & 400 MD
LD_1600_CKDEPI$AMT <- ifelse(LD_1600_CKDEPI$TIME == 0, 1600, ifelse(LD_1600_CKDEPI$TAD == 0, 400, 0))

### Combining the 5 datasets  -------------------------------

# Remove ID_CKDEPI column from LD_800_CKDEPI, LD_1000_CKDEPI, LD_1200_CKDEPI, 
# LD_1400_CKDEPI, LD_1600_CKDEPI

LD_800_CKDEPI$ID_CKDEPI <- NULL
LD_1000_CKDEPI$ID_CKDEPI <- NULL
LD_1200_CKDEPI$ID_CKDEPI <- NULL
LD_1400_CKDEPI$ID_CKDEPI <- NULL
LD_1600_CKDEPI$ID_CKDEPI <- NULL

# Update ID of LD_1000_CKDEPI, from 1 to 440 to 441 to 880
LD_1000_CKDEPI$ID <- LD_1000_CKDEPI$ID + 440

# Update ID of LD_1200_CKDEPI, from 1 to 440 to 881 to 1320
LD_1200_CKDEPI$ID <- LD_1200_CKDEPI$ID + 880

# Update ID of LD_1400_CKDEPI, from 1 to 440 to 1321 to 1760
LD_1400_CKDEPI$ID <- LD_1400_CKDEPI$ID + 1320

# Update ID of LD_1600_CKDEPI, from 1 to 440 to 1761 to 2200
LD_1600_CKDEPI$ID <- LD_1600_CKDEPI$ID + 1760

# Combine all datasets
BW_stra_LD_CKDEPI <- rbind(LD_800_CKDEPI, LD_1000_CKDEPI, LD_1200_CKDEPI, 
                           LD_1400_CKDEPI, LD_1600_CKDEPI)

# Then export into BW_stra_LD_CKDEPI.csv dataset 
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Dosing_simulations/Dose_finding/CKDEPI_Opt_dose")
write.csv(BW_stra_LD_CKDEPI, "BW_stra_LD_CKDEPI.csv", quote = F, row.names = F)

# VPC ------------------------------------------

## Refine vpc datasets, exclude all observations having either DV>70 or TAD>25.5

# Import the overall dataset
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/Imputed")

Fluco_revised_imputed_overall <- read.csv("Fluco_revised_imputed_overall.csv")

# Convert DV & TAD into numeric
Fluco_revised_imputed_overall$DV <- as.numeric(Fluco_revised_imputed_overall$DV)
Fluco_revised_imputed_overall$TAD <- as.numeric(Fluco_revised_imputed_overall$TAD)

# The refined dataset
Fluco_revised_imputed_overall_refined <- Fluco_revised_imputed_overall  |> 
  dplyr::filter((TAD <= 25.5 | is.na(DV)))

Fluco_revised_imputed_overall_refined <- Fluco_revised_imputed_overall_refined |> 
  dplyr::filter((DV <= 70 | is.na(DV)))
# eliminated 25 concentration - time points

# Replace NAs in DV & TAD with "."
Fluco_revised_imputed_overall_refined$DV[is.na(Fluco_revised_imputed_overall_refined$DV)] <- "."
Fluco_revised_imputed_overall_refined$TAD[is.na(Fluco_revised_imputed_overall_refined$TAD)] <- "."

# Export into refined dataset
setwd("C:/Users/u0164053/OneDrive - KU Leuven/Fluconazole PoPPK/Fluconazol_project/Revision 210324/Datasets/VPC")
write.csv(Fluco_revised_imputed_overall_refined, "Fluco_revised_imputed_overall_refined.csv", quote = F, row.names = F)

# Population level simulation ------------------------------------------

## Create dataset for conc-time plot -----------------------------------

### Generate BW, CKDEPI, and CRRT -------------------------------

# load the dataset
setwd(dirname(getwd())) # set wd to the parent folder
dosing_CRRT <- read.csv("Pop_sim/dosing_03.csv")

## create a dataset of 1000 patients from dosing_CRRT dataset

# Create a new dataset with only patient 13
dosing_CRRT_13 <- dosing_CRRT[dosing_CRRT$ID == 13, ]

# Create 987 copies of the patient 13 dataset and assign IDs 14-1000
new_patients <- lapply(14:1000, function(id){
  new_df <- dosing_CRRT_13
  new_df$ID <- id
  return(new_df)
})

# Combine the new dataset with the original dataset
dosing_CRRT <- rbind(dosing_CRRT, do.call(rbind, new_patients))

# Rename FFM to BW
names(dosing_CRRT)[names(dosing_CRRT) == "FFM"] <- "BW"

## Assign BW to the new dataset

# for reproducibility
set.seed(123)

# Create the empirical cumulative distribution function (ECDF)
pop_BW <- Fluco_clean_revised  |>
  group_by(ID) |>
  distinct(BW2, .keep_all = TRUE) |>
  dplyr::filter(!is.na(BW2)) |>
  ungroup() |>
  select(BW2)
ecdf_BW <- ecdf(pop_BW$BW2)

# Generate random values from the ECDF
random_values <- runif(1000)

# Map the random values to the empirical distribution using the inverse transform sampling method
sampled_BW <- round(quantile(pop_BW$BW2, random_values), 1)

# Assign the sampled BW values to dosing_CRRT dataset
dosing_CRRT$BW <- sampled_BW[match(dosing_CRRT$ID, 1:1000)]

## Assign CKDEPI to the new dataset

# for reproducibility
set.seed(123)

# Create the empirical cumulative distribution function (ECDF)
pop_CKDEPI <- Fluco_clean_revised  |>  
  dplyr::filter(!is.na(CKDEPI)) |>
  select(CKDEPI)
ecdf_CKDEPI <- ecdf(pop_CKDEPI$CKDEPI)

# Generate random values from the ECDF
random_values <- runif(1000)

# Map the random values to the empirical distribution using the inverse transform sampling method
sampled_CKDEPI <- round(quantile(pop_CKDEPI$CKDEPI, random_values), 1)

# Assign the sampled BW values to dosing_CRRT dataset
dosing_CRRT$CKDEPI <- sampled_CKDEPI[match(dosing_CRRT$ID, 1:1000)]

## Create the dosing_CRRT dataset with 10.2% full-time CRRT, and 9% part time 
## CRRT (in which 58.5% DAYS on CRRT)

# Create a dummy dataset with 1000 patient IDs and 14 days:

dummy_data <- data.frame(
  ID = rep(1:1000, each = 14),
  DAY = rep(1:14, times = 1000)
)

# Add the CRRT variable to the dummy dataset and assign it temporarily to 0:
dummy_data$CRRT <- 0

# Randomly select 10.2% of the patient IDs and assign their CRRT to 1:
set.seed(123)  
selected_ids <- sample(unique(dummy_data$ID), size = round(0.102 * 1000))
dummy_data$CRRT[dummy_data$ID %in% selected_ids] <- 1

# Randomly select another 9% of patient IDs that are not in selected_ids:
available_ids <- setdiff(unique(dummy_data$ID), selected_ids)
selected_ids_2 <- sample(available_ids, size = round(0.09 * 1000))

# Assign 58.5% of the selected days to CRRT == 1, and the remaining to CRRT == 0
selected_days <- sample(unique(dummy_data$DAY), size = round(0.585 * 14))
dummy_data$CRRT[dummy_data$ID %in% selected_ids_2 & dummy_data$DAY %in% selected_days] <- 1

# Assign the CRRT information from the dummy dataset to the dosing_CRRT dataset based on ID and DAY:
dosing_CRRT <- merge(dosing_CRRT, dummy_data[, c("ID", "DAY", "CRRT")], by = c("ID", "DAY"), all.x = TRUE)
dosing_CRRT$CRRT.x <- NULL
colnames(dosing_CRRT)[colnames(dosing_CRRT) == "CRRT.y"] <- "CRRT"

### Updating TIME and TAD -------------------------------

## Create a dataset for conc time curve opt dosing

# Repeat the 2nd row of each ID per DAY 14 times
rows = seq(2, 41999, by = 3)
times = 14
dosing_CRRT_repeated <- dosing_CRRT[rep(rows, times),]

# Update the TAD values in dosing_CRRT 
dosing_CRRT_combined <- rbind(dosing_CRRT,dosing_CRRT_repeated)
dosing_CRRT_combined <- arrange(dosing_CRRT_combined, ID, DAY, TAD, TIME)
desired_tad <- c(0, 0.1, 0.25, 0.5, 0.6, 1, 1.5, 2, 3, 4, 6, 8, 12, 18, 21, 22.5, 23.9)
dosing_CRRT_combined$TAD <- rep(desired_tad, 14000)

# Reset row names
rownames(dosing_CRRT_combined) <- 1:238000

# Remove the PEAK and TROUGH columns
dosing_CRRT_combined$PEAK <- NULL
dosing_CRRT_combined$TROUGH <- NULL

# Update TIME
dosing_CRRT_combined$TIME <- (dosing_CRRT_combined$DAY - 1) * 24 + 
  dosing_CRRT_combined$TAD

### Update dosing information -------------------------------

#### Optimised -------------------------------

dosing_CRRT_combined <- dosing_CRRT_combined %>%
  mutate(
    AMT = ifelse(TIME == 0 & BW > 120 & CRRT == 1, 1800,
                 ifelse(TIME == 0 & BW <= 120 & BW > 100 & CRRT == 1, 1600,
                        ifelse(TIME == 0 & BW <= 100 & BW > 80 & CRRT == 1, 1400,
                               ifelse(TIME == 0 & BW <= 80 & BW > 60 & CRRT == 1, 1200, 
                                      ifelse(TIME == 0 & BW <= 60 & CRRT == 1, 1000,
                                             ifelse(TIME == 0 & BW > 120 & CRRT == 0, 1600,
                                                    ifelse(TIME == 0 & BW <= 120 & BW > 100 & CRRT == 0, 1400,
                                                           ifelse(TIME == 0 & BW <= 100 & BW > 80 & CRRT == 0, 1200,
                                                                  ifelse(TIME == 0 & BW <= 80 & BW > 60 & CRRT == 0, 1000,
                                                                         ifelse(TIME == 0 & BW <= 60 & CRRT == 0, 800,
                                                                                ifelse(TIME > 0 & MDV == 1 & CRRT == 1, 800,
                                                                                       ifelse(TIME > 0 & MDV == 1 & CRRT == 0, 400,0                                                                                    ))))
                                                                  )))
                                             )))))
  )

#### Standard -------------------------------

dosing_CRRT_std <- dosing_CRRT_combined %>%
  mutate(
    AMT = ifelse(TIME == 0, 800,
                 ifelse(TIME > 0 & TAD == 0, 400, 0
                 )
    )
  )

# Export pop_opt_conc & pop_std_conc datasets
setwd("./Pop_sim")
write.csv(dosing_CRRT_combined,"pop_opt_conc.csv", quote = F, row.names = F)
write.csv(dosing_CRRT_std,"pop_std_conc.csv", quote = F, row.names = F)

## Create datasets for boxplot & ROC --------------------------------

# Only retain observations with 0 & 23.9 hours TAD
dosing_CRRT_combined <- read.csv("Pop_sim/pop_opt_conc.csv")
dosing_CRRT_std <- read.csv("Pop_sim/pop_std_conc.csv")
pop_box_std <- dosing_CRRT_std[dosing_CRRT_std$TAD %in% c(0, 23.9), ]
pop_box_opt <- dosing_CRRT_combined[dosing_CRRT_combined$TAD %in% c(0, 23.9), ]

# Export pop_box_std & pop_box_opt datasets
setwd("./Pop_sim")
write.csv(pop_box_std,"pop_box_std.csv", quote = F, row.names = F)
write.csv(pop_box_opt,"pop_box_opt.csv", quote = F, row.names = F)

