# Load libraries -------------------------------------------

library(ggplot2)
library(xpose)
library(dplyr)
library(PsNR)
library(magrittr)
library(methods)
library(xpose4)
library(vpc)
library(cowplot)
library(extrafont)

# VPC overall -------------------------------------------

setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/VPC/")

Number_of_simulations_performed <- 1000
PI <- 95
CI <- 95
perc_PI <- c(0 + (1 - PI/100) / 2, 1 - (1 - PI/100) / 2)
perc_CI <- c(0 + (1 - CI/100) / 2, 1 - (1 - CI/100) / 2)
bin_times <- c(0.01549, 1.33, 1.77, 3.665, 5.81, 7.275, 9.6, 10.77, 12.3835, 13.135, 14.425, 15.485, 16.905, 18.035, 19.54, 21.225, 22.225, 23.1585, 24.63, 25.52549)
Number_of_bins <- length(bin_times) - 1

# Load simulation data
simulations_tablefile_mod_pooled <- "vpc_sim_overall_refined.1.npctab.dta" # better use this
dataframe_simulations_mod_pooled <- read_nonmem_table(simulations_tablefile_mod_pooled)
dataframe_simulations_mod_pooled <- dataframe_simulations_mod_pooled[dataframe_simulations_mod_pooled$MDV == 0,]
dataframe_simulations_mod_pooled$replicate <- rep(1:Number_of_simulations_performed, each = nrow(dataframe_simulations_mod_pooled) / Number_of_simulations_performed)

# Set vector with unique replicates
Replicate_vector_mod_pooled <- unique(dataframe_simulations_mod_pooled$replicate)

# Define BIN based on TAD
for (i in Number_of_bins:1) {
  dataframe_simulations_mod_pooled$BIN[dataframe_simulations_mod_pooled$TAD <= bin_times[i + 1]] <- i
}

# Load observations data
observations_tablefile_mod_pooled <- "vpc_original_overall_refined.npctab.dta" #better use the refined
Obs_mod_pooled <- read_nonmem_table(observations_tablefile_mod_pooled)
Obs_mod_pooled <- Obs_mod_pooled[Obs_mod_pooled$MDV == 0,]

# Define BIN based on TAD
for (i in Number_of_bins:1) {
  Obs_mod_pooled$BIN[Obs_mod_pooled$TAD <= bin_times[i + 1]] <- i
}

  # Calculate PRED_BIN and merge
  PRED_BIN_mod_pooled <- dataframe_simulations_mod_pooled[dataframe_simulations_mod_pooled$replicate == 1, ] %>%
    group_by(BIN) %>%
    summarize(PREDBIN = median(PRED))
  
  dataframe_simulations_mod_pooled <- merge(dataframe_simulations_mod_pooled, PRED_BIN_mod_pooled, by = 'BIN')
  
  # Calculate PCDV
  dataframe_simulations_mod_pooled$PCDV <- dataframe_simulations_mod_pooled$DV * (dataframe_simulations_mod_pooled$PREDBIN / dataframe_simulations_mod_pooled$PRED)
  
  # Order data
  dataframe_simulations_mod_pooled <- dataframe_simulations_mod_pooled[order(dataframe_simulations_mod_pooled$replicate, dataframe_simulations_mod_pooled$ID, dataframe_simulations_mod_pooled$TAD),]
  
  # Initialize storage for PI data
  sim_PI_mod_pooled <- NULL
  
  # Loop through replicates
  for (i in Replicate_vector_mod_pooled) {
    sim_vpc_ci_mod_pooled <- dataframe_simulations_mod_pooled %>%
      dplyr::filter(replicate %in% i) %>%
      group_by(BIN) %>%
      summarize(
        C_median = median(PCDV),
        C_lower = quantile(PCDV, perc_PI[1]),
        C_upper = quantile(PCDV, perc_PI[2])
      ) %>%
      mutate(replicate = i)
    
    sim_PI_mod_pooled <- rbind(sim_PI_mod_pooled, sim_vpc_ci_mod_pooled)
  }
  
  # Calculate CI
  sim_CI_mod_pooled <- sim_PI_mod_pooled %>%
    group_by(BIN) %>%
    summarize(
      C_median_CI_lwr = quantile(C_median, perc_CI[1]),
      C_median_CI_upr = quantile(C_median, perc_CI[2]),
      C_low_lwr = quantile(C_lower, perc_CI[1]),
      C_low_upr = quantile(C_lower, perc_CI[2]),
      C_up_lwr = quantile(C_upper, perc_CI[1]),
      C_up_upr = quantile(C_upper, perc_CI[2])
    )
  
  # Assign NA values
  sim_CI_mod_pooled$x1 <- NA
  sim_CI_mod_pooled$x2 <- NA
  
  # Assign BIN values based on bin_times
  for (i in 1:Number_of_bins) {
    sim_CI_mod_pooled$x1[sim_CI_mod_pooled$BIN == i] <- bin_times[i]
    sim_CI_mod_pooled$x2[sim_CI_mod_pooled$BIN == i] <- bin_times[i + 1]
  }
  
    # Merge with PRED_BIN_mod_pooled
  Obs_mod_pooled <- merge(Obs_mod_pooled, PRED_BIN_mod_pooled, by = 'BIN')
  
  # Calculate PCDV
  Obs_mod_pooled$PCDV <- Obs_mod_pooled$DV * (Obs_mod_pooled$PREDBIN / Obs_mod_pooled$PRED)
  
  # Summarize data
  obs_vpc_mod_pooled <- Obs_mod_pooled %>%
    group_by(BIN) %>% 
    summarize(
      C_median = median(PCDV, na.rm = TRUE),
      C_lower = quantile(PCDV, perc_PI[1], na.rm = TRUE),
      C_upper = quantile(PCDV, perc_PI[2], na.rm = TRUE)
    )
  
  bin_middle_mod_pooled <- Obs_mod_pooled %>%
    group_by(BIN) %>%
    summarize(bin_middle = median(TAD, na.rm = TRUE))
  
  bin_middle_mod_pooled <- bin_middle_mod_pooled$bin_middle
  
  obs_vpc_mod_pooled$TAD <- bin_middle_mod_pooled
  
  # Create VPC plot
  CI_VPC_lin_mod_pooled <- ggplot() +
    geom_rect(
      data = sim_CI_mod_pooled,
      mapping = aes(xmin = x1, xmax = x2, ymin = C_median_CI_lwr, ymax = C_median_CI_upr),
      fill = '#35b779',
      alpha = 0.25
    ) +
    geom_rect(
      data = sim_CI_mod_pooled,
      mapping = aes(xmin = x1, xmax = x2, ymin = C_low_lwr, ymax = C_low_upr),
      fill = '#31688e',
      alpha = 0.25
    ) +
    geom_rect(
      data = sim_CI_mod_pooled,
      mapping = aes(xmin = x1, xmax = x2, ymin = C_up_lwr, ymax = C_up_upr),
      fill = '#31688e',
      alpha = 0.25
    ) +
    geom_line(
      data = obs_vpc_mod_pooled,
      aes(TAD, C_median),
      col = '#440154',
      linetype = 'dashed',
      linewidth = 0.8
    ) +
    geom_line(
      data = obs_vpc_mod_pooled,
      aes(TAD, C_lower),
      col = '#440154',
      linetype = 'dashed',
      linewidth = 0.8
    ) +
    geom_line(
      data = obs_vpc_mod_pooled,
      aes(TAD, C_upper),
      col = '#440154',
      linetype = 'dashed',
      linewidth = 0.8
    ) +
    geom_point(
      data = Obs_mod_pooled,
      aes(x = TAD, y = PCDV),
      color = '#440154',
      alpha = 0.3
    ) +
    scale_x_continuous(
      limits = c(0, 26),
      breaks = seq(0, 25, by = 5),
      name = "Time since last dose (hours)",
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 80),
      breaks = seq(0, 80, by = 10),
      name = "Prediction-corrected fluconazole concentration (mg/L)",
      expand = c(0.01, 0)
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
          axis.title = element_text(size = 7, family = "Helvetica"),
          axis.text = element_text(size = 7, family = "Helvetica"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, color = '#31688e')
  
# Save the combine plot
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Plots/VPC/")

ggsave("vpc_overall.svg",
       CI_VPC_lin_mod_pooled,
       dpi = 300, 
       width = 9, 
       height = 9,
       unit = "cm")

ggsave("vpc_overall.JPEG",
       CI_VPC_lin_mod_pooled,
       dpi = 300, 
       width = 9, 
       height = 9,
       unit = "cm")

ggsave("vpc_overall.EPS",
       CI_VPC_lin_mod_pooled,
       dpi = 300, 
       width = 9, 
       height = 9,
       unit = "cm")

# VPC stratified by HOSPITAL -------------------------------------------

setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/VPC/")

Number_of_simulations_performed <- 1000
PI <- 95
CI <- 95
perc_PI <- c(0 + (1 - PI/100) / 2, 1 - (1 - PI/100) / 2)
perc_CI <- c(0 + (1 - CI/100) / 2, 1 - (1 - CI/100) / 2)
bin_times <- c(0.01549, 1.33, 1.77, 3.665, 5.81, 7.275, 9.6, 10.77, 12.3835, 13.135, 14.425, 15.485, 16.905, 18.035, 19.54, 21.225, 22.225, 23.1585, 24.63, 25.52549)
Number_of_bins <- length(bin_times) - 1

# Load simulation data
simulations_tablefile_mod_pooled <- "vpc_sim_overall_refined.1.npctab.dta" # better use this
dataframe_simulations_mod_pooled <- read_nonmem_table(simulations_tablefile_mod_pooled)
dataframe_simulations_mod_pooled <- dataframe_simulations_mod_pooled[dataframe_simulations_mod_pooled$MDV == 0,]
dataframe_simulations_mod_pooled$replicate <- rep(1:Number_of_simulations_performed, each = nrow(dataframe_simulations_mod_pooled) / Number_of_simulations_performed)

# Set vector with unique replicates
Replicate_vector_mod_pooled <- unique(dataframe_simulations_mod_pooled$replicate)

# Define BIN based on TAD
for (i in Number_of_bins:1) {
  dataframe_simulations_mod_pooled$BIN[dataframe_simulations_mod_pooled$TAD <= bin_times[i + 1]] <- i
}

# Create new column called HOSPITAL
dataframe_simulations_mod_pooled$HOSPITAL <- as.numeric(substr(dataframe_simulations_mod_pooled$ID, 1, 1))

# Load observations data
observations_tablefile_mod_pooled <- "vpc_original_overall_refined.npctab.dta" #better use the refined
Obs_mod_pooled <- read_nonmem_table(observations_tablefile_mod_pooled)
Obs_mod_pooled <- Obs_mod_pooled[Obs_mod_pooled$MDV == 0,]

# Define BIN based on TAD
for (i in Number_of_bins:1) {
  Obs_mod_pooled$BIN[Obs_mod_pooled$TAD <= bin_times[i + 1]] <- i
}

# Create new column called HOSPITAL observed dataset
Obs_mod_pooled$HOSPITAL <- as.numeric(substr(Obs_mod_pooled$ID, 1, 1))

## Special case HOSPITAL 1 -------------------------------------------

sim_mod_pooled_flag1 <- dataframe_simulations_mod_pooled %>% 
  dplyr::filter(HOSPITAL == 1)

# Calculate PRED_BIN and merge
PRED_BIN_mod_pooled_flag1 <- sim_mod_pooled_flag1[sim_mod_pooled_flag1$replicate == 1, ] %>%
  group_by(BIN) %>%
  summarize(PREDBIN = median(PRED))

sim_mod_pooled_flag1 <- merge(sim_mod_pooled_flag1, PRED_BIN_mod_pooled_flag1, by = 'BIN')

# Calculate PCDV
sim_mod_pooled_flag1$PCDV <- sim_mod_pooled_flag1$DV * (sim_mod_pooled_flag1$PREDBIN / sim_mod_pooled_flag1$PRED)

# Order data
sim_mod_pooled_flag1 <- sim_mod_pooled_flag1[order(sim_mod_pooled_flag1$replicate, sim_mod_pooled_flag1$ID, sim_mod_pooled_flag1$TAD),]

# Initialize storage for PI data
sim_PI_mod_pooled_flag1 <- NULL

# Loop through replicates
for (i in Replicate_vector_mod_pooled) {
  sim_vpc_ci_mod_pooled_flag1 <- sim_mod_pooled_flag1 %>%
    dplyr::filter(replicate %in% i) %>%
    group_by(BIN) %>%
    summarize(
      C_median = median(PCDV),
      C_lower = quantile(PCDV, perc_PI[1]),
      C_upper = quantile(PCDV, perc_PI[2])
    ) %>%
    mutate(replicate = i)
  
  sim_PI_mod_pooled_flag1 <- rbind(sim_PI_mod_pooled_flag1, sim_vpc_ci_mod_pooled_flag1)
}

# Calculate CI
sim_CI_mod_pooled_flag1 <- sim_PI_mod_pooled_flag1 %>%
  group_by(BIN) %>%
  summarize(
    C_median_CI_lwr = quantile(C_median, perc_CI[1]),
    C_median_CI_upr = quantile(C_median, perc_CI[2]),
    C_low_lwr = quantile(C_lower, perc_CI[1]),
    C_low_upr = quantile(C_lower, perc_CI[2]),
    C_up_lwr = quantile(C_upper, perc_CI[1]),
    C_up_upr = quantile(C_upper, perc_CI[2])
  )

# Assign NA values
sim_CI_mod_pooled_flag1$x1 <- NA
sim_CI_mod_pooled_flag1$x2 <- NA

# Assign BIN values based on bin_times
for (i in 1:Number_of_bins) {
  sim_CI_mod_pooled_flag1$x1[sim_CI_mod_pooled_flag1$BIN == i] <- bin_times[i]
  sim_CI_mod_pooled_flag1$x2[sim_CI_mod_pooled_flag1$BIN == i] <- bin_times[i + 1]
}

# dplyr::filter Obs_mod_pooled for the current hospital
Obs_mod_pooled_flag1 <- Obs_mod_pooled %>% 
  dplyr::filter(HOSPITAL == hosp)

# Merge with PRED_BIN_mod_pooled_flag1
Obs_mod_pooled_flag1 <- merge(Obs_mod_pooled_flag1, PRED_BIN_mod_pooled_flag1, by = 'BIN')

# Calculate PCDV
Obs_mod_pooled_flag1$PCDV <- Obs_mod_pooled_flag1$DV * (Obs_mod_pooled_flag1$PREDBIN / Obs_mod_pooled_flag1$PRED)

# Summarize data
obs_vpc_mod_pooled_flag1 <- Obs_mod_pooled_flag1 %>%
  group_by(BIN) %>% 
  summarize(
    C_median = median(PCDV, na.rm = TRUE),
    C_lower = quantile(PCDV, perc_PI[1], na.rm = TRUE),
    C_upper = quantile(PCDV, perc_PI[2], na.rm = TRUE)
  )

bin_middle_mod_pooled_flag1 <- Obs_mod_pooled_flag1 %>%
  group_by(BIN) %>%
  summarize(bin_middle = median(TAD, na.rm = TRUE))

bin_middle_mod_pooled_flag1 <- bin_middle_mod_pooled_flag1$bin_middle

obs_vpc_mod_pooled_flag1$TAD <- bin_middle_mod_pooled_flag1

# Create VPC plot
CI_VPC_lin_mod_pooled_flag1 <- ggplot() +
  geom_rect(
    data = sim_CI_mod_pooled_flag1,
    mapping = aes(xmin = x1, xmax = x2, ymin = C_median_CI_lwr, ymax = C_median_CI_upr),
    fill = '#35b779',
    alpha = 0.25
  ) +
  geom_rect(
    data = sim_CI_mod_pooled_flag1,
    mapping = aes(xmin = x1, xmax = x2, ymin = C_low_lwr, ymax = C_low_upr),
    fill = '#31688e',
    alpha = 0.25
  ) +
  geom_rect(
    data = sim_CI_mod_pooled_flag1,
    mapping = aes(xmin = x1, xmax = x2, ymin = C_up_lwr, ymax = C_up_upr),
    fill = '#31688e',
    alpha = 0.25
  ) +
  geom_line(
    data = obs_vpc_mod_pooled_flag1,
    aes(TAD, C_median),
    col = '#440154',
    linetype = 'dashed',
    linewidth = 0.8
  ) +
  geom_line(
    data = obs_vpc_mod_pooled_flag1,
    aes(TAD, C_lower),
    col = '#440154',
    linetype = 'dashed',
    linewidth = 0.8
  ) +
  geom_line(
    data = obs_vpc_mod_pooled_flag1,
    aes(TAD, C_upper),
    col = '#440154',
    linetype = 'dashed',
    linewidth = 0.8
  ) +
  geom_point(
    data = Obs_mod_pooled_flag1,
    aes(x = TAD, y = PCDV),
    color = '#440154',
    alpha = 0.3
  ) +
  scale_x_continuous(
    limits = c(0, 26),
    breaks = seq(0, 25, by = 5),
    name = "Time since last dose (hours)",
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 80),
    breaks = seq(0, 80, by = 10),
    name = "Prediction-corrected fluconazole concentration (mg/L)",
    expand = c(0.01, 0)
  ) +
  ggtitle("Van Daele et al.") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
        axis.title = element_text(size = 7, family = "Helvetica"),
        axis.text = element_text(size = 7, family = "Helvetica"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, color = '#31688e')

## HOSPITAL 2 to 8 -------------------------------------------

# List to store plots for HOSPITAL 2 to 8
hospital_plots <- list()

# Title of 7 plots
hospital_titles <- c(
  "Muilwijk et al.", "Bergner et al.", "Buijk et al.", "Sandaradura et al.",
  "Sinnollareddy et al.", "Alobaid et al.", "Boonstra et al."
)

# Loop through HOSPITAL values (2 to 8)
for (hosp in 2:8) { #Hospital 1 is an outlier, therefore I selected hospital 2 to 8 instead
  # Filter simulation data for the current hospital
  sim_mod_pooled_flag <- dataframe_simulations_mod_pooled %>% dplyr::filter(HOSPITAL == hosp)
  
  # Calculate PRED_BIN and merge
  PRED_BIN_mod_pooled_flag <- sim_mod_pooled_flag[sim_mod_pooled_flag$replicate == 1, ] %>%
    group_by(BIN) %>%
    summarize(PREDBIN = median(PRED))
  
  sim_mod_pooled_flag <- merge(sim_mod_pooled_flag, PRED_BIN_mod_pooled_flag, by = 'BIN')
  
  # Calculate PCDV
  sim_mod_pooled_flag$PCDV <- sim_mod_pooled_flag$DV * (sim_mod_pooled_flag$PREDBIN / sim_mod_pooled_flag$PRED)
  
  # Order data
  sim_mod_pooled_flag <- sim_mod_pooled_flag[order(sim_mod_pooled_flag$replicate, sim_mod_pooled_flag$ID, sim_mod_pooled_flag$TAD),]
  
  # Initialize storage for PI data
  sim_PI_mod_pooled_flag <- NULL
  
  # Loop through replicates
  for (i in Replicate_vector_mod_pooled) {
    sim_vpc_ci_mod_pooled_flag <- sim_mod_pooled_flag %>%
      dplyr::filter(replicate %in% i) %>%
      group_by(BIN) %>%
      summarize(
        C_median = median(PCDV),
        C_lower = quantile(PCDV, perc_PI[1]),
        C_upper = quantile(PCDV, perc_PI[2])
      ) %>%
      mutate(replicate = i)
    
    sim_PI_mod_pooled_flag <- rbind(sim_PI_mod_pooled_flag, sim_vpc_ci_mod_pooled_flag)
  }
  
  # Calculate CI
  sim_CI_mod_pooled_flag <- sim_PI_mod_pooled_flag %>%
    group_by(BIN) %>%
    summarize(
      C_median_CI_lwr = quantile(C_median, perc_CI[1]),
      C_median_CI_upr = quantile(C_median, perc_CI[2]),
      C_low_lwr = quantile(C_lower, perc_CI[1]),
      C_low_upr = quantile(C_lower, perc_CI[2]),
      C_up_lwr = quantile(C_upper, perc_CI[1]),
      C_up_upr = quantile(C_upper, perc_CI[2])
    )
  
  # Assign NA values
  sim_CI_mod_pooled_flag$x1 <- NA
  sim_CI_mod_pooled_flag$x2 <- NA
  
  # Assign BIN values based on bin_times
  for (i in 1:Number_of_bins) {
    sim_CI_mod_pooled_flag$x1[sim_CI_mod_pooled_flag$BIN == i] <- bin_times[i]
    sim_CI_mod_pooled_flag$x2[sim_CI_mod_pooled_flag$BIN == i] <- bin_times[i + 1]
  }
  
    # dplyr::filter Obs_mod_pooled for the current hospital
  Obs_mod_pooled_flag <- Obs_mod_pooled %>% 
    dplyr::filter(HOSPITAL == hosp)
  
  # Merge with PRED_BIN_mod_pooled_flag
  Obs_mod_pooled_flag <- merge(Obs_mod_pooled_flag, PRED_BIN_mod_pooled_flag, by = 'BIN')
  
  # Calculate PCDV
  Obs_mod_pooled_flag$PCDV <- Obs_mod_pooled_flag$DV * (Obs_mod_pooled_flag$PREDBIN / Obs_mod_pooled_flag$PRED)
  
  # Summarize data
  obs_vpc_mod_pooled_flag <- Obs_mod_pooled_flag %>%
    group_by(BIN) %>% 
    summarize(
      C_median = median(PCDV, na.rm = TRUE),
      C_lower = quantile(PCDV, perc_PI[1], na.rm = TRUE),
      C_upper = quantile(PCDV, perc_PI[2], na.rm = TRUE)
    )
  
  bin_middle_mod_pooled_flag <- Obs_mod_pooled_flag %>%
    group_by(BIN) %>%
    summarize(bin_middle = median(TAD, na.rm = TRUE))
  
  bin_middle_mod_pooled_flag <- bin_middle_mod_pooled_flag$bin_middle
  
  obs_vpc_mod_pooled_flag$TAD <- bin_middle_mod_pooled_flag
  
  # Create VPC plot
  CI_VPC_lin_mod_pooled_flag <- ggplot() +
    geom_rect(
      data = sim_CI_mod_pooled_flag,
      mapping = aes(xmin = x1, xmax = x2, ymin = C_median_CI_lwr, ymax = C_median_CI_upr),
      fill = '#35b779',
      alpha = 0.25
    ) +
    geom_rect(
      data = sim_CI_mod_pooled_flag,
      mapping = aes(xmin = x1, xmax = x2, ymin = C_low_lwr, ymax = C_low_upr),
      fill = '#31688e',
      alpha = 0.25
    ) +
    geom_rect(
      data = sim_CI_mod_pooled_flag,
      mapping = aes(xmin = x1, xmax = x2, ymin = C_up_lwr, ymax = C_up_upr),
      fill = '#31688e',
      alpha = 0.25
    ) +
    geom_line(
      data = obs_vpc_mod_pooled_flag,
      aes(TAD, C_median),
      col = '#440154',
      linetype = 'dashed',
      linewidth = 0.8
    ) +
    geom_line(
      data = obs_vpc_mod_pooled_flag,
      aes(TAD, C_lower),
      col = '#440154',
      linetype = 'dashed',
      linewidth = 0.8
    ) +
    geom_line(
      data = obs_vpc_mod_pooled_flag,
      aes(TAD, C_upper),
      col = '#440154',
      linetype = 'dashed',
      linewidth = 0.8
    ) +
    geom_point(
      data = Obs_mod_pooled_flag,
      aes(x = TAD, y = PCDV),
      color = '#440154',
      alpha = 0.3
    ) +
    scale_x_continuous(
      limits = c(0, 26),
      breaks = seq(0, 25, by = 5),
      name = "Time since last dose (hours)",
      expand = c(0.01, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 80),
      breaks = seq(0, 80, by = 10),
      name = "Prediction-corrected fluconazole concentration (mg/L)",
      expand = c(0.01, 0)
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
          axis.title = element_text(size = 7, family = "Helvetica"),
          axis.text = element_text(size = 7, family = "Helvetica"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, color = '#31688e')
    
  # Save the plot in the list with the corresponding title
  hospital_plots[[as.character(hosp)]] <- CI_VPC_lin_mod_pooled_flag +
    ggtitle(hospital_titles[hosp - 1])  # Subtract 1 because the index starts from 2
}
  
# Add this after making plot for hospital 1
hospital_plots <- c(list("1" = CI_VPC_lin_mod_pooled_flag1), hospital_plots)

# Combine all the plots
combined_plot <- cowplot::plot_grid(
  plotlist = lapply(names(hospital_plots), function(hosp) {
    return(hospital_plots[[hosp]])
  }),
  nrow = 2,
  labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'),
  label_size = 8
)

# Save the combine plot
setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Plots/VPC/")

ggsave("vpc_per_study.svg",
       combined_plot,
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

ggsave("vpc_per_study.JPEG",
       combined_plot,
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

ggsave("vpc_per_study.EPS",
       combined_plot,
       dpi = 300, 
       width = 19, 
       height = 19,
       unit = "cm")

# VPC stratified by CRRT status -------------------------------------------
 
  setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/VPC/")
  
  Number_of_simulations_performed <- 1000
  PI <- 95
  CI <- 95
  perc_PI <- c(0 + (1 - PI/100) / 2, 1 - (1 - PI/100) / 2)
  perc_CI <- c(0 + (1 - CI/100) / 2, 1 - (1 - CI/100) / 2)
  bin_times <- c(0.01549, 1.33, 1.77, 3.665, 5.81, 7.275, 9.6, 10.77, 12.3835, 13.135, 14.425, 15.485, 16.905, 18.035, 19.54, 21.225, 22.225, 23.1585, 24.63, 25.52549)
  Number_of_bins <- length(bin_times) - 1
  
  # Load simulation data
  simulations_tablefile_mod_pooled <- "vpc_sim_overall_refined.1.npctab.dta" # better use this
  dataframe_simulations_mod_pooled <- read_nonmem_table(simulations_tablefile_mod_pooled)
  dataframe_simulations_mod_pooled <- dataframe_simulations_mod_pooled[dataframe_simulations_mod_pooled$MDV == 0,]
  dataframe_simulations_mod_pooled$replicate <- rep(1:Number_of_simulations_performed, each = nrow(dataframe_simulations_mod_pooled) / Number_of_simulations_performed)
  
  # Set vector with unique replicates
  Replicate_vector_mod_pooled <- unique(dataframe_simulations_mod_pooled$replicate)
  
  # Define BIN based on TAD
  for (i in Number_of_bins:1) {
    dataframe_simulations_mod_pooled$BIN[dataframe_simulations_mod_pooled$TAD <= bin_times[i + 1]] <- i
  }
  
  # Load observations data
  observations_tablefile_mod_pooled <- "vpc_original_overall_refined.npctab.dta" #better use the refined
  Obs_mod_pooled <- read_nonmem_table(observations_tablefile_mod_pooled)
  Obs_mod_pooled <- Obs_mod_pooled[Obs_mod_pooled$MDV == 0,]
  
  # Define BIN based on TAD
  for (i in Number_of_bins:1) {
    Obs_mod_pooled$BIN[Obs_mod_pooled$TAD <= bin_times[i + 1]] <- i
  }
  
  # List to store plots for HOSPITAL 2 to 8
  crrt_plots <- list()
  
  # Title of 7 plots
  crrt_titles <- c(
    "Off-CRRT","On-CRRT"
  )
  
  # Loop through CRRT values (0 to 1)
  for (crrt in 0:1) { 
    sim_mod_pooled_flag <- dataframe_simulations_mod_pooled %>% 
      dplyr::filter(CRRT == crrt)
    
    # Calculate PRED_BIN and merge
    PRED_BIN_mod_pooled_flag <- sim_mod_pooled_flag[sim_mod_pooled_flag$replicate == 1, ] %>%
      group_by(BIN) %>%
      summarize(PREDBIN = median(PRED))
    
    sim_mod_pooled_flag <- merge(sim_mod_pooled_flag, PRED_BIN_mod_pooled_flag, by = 'BIN')
    
    # Calculate PCDV
    sim_mod_pooled_flag$PCDV <- sim_mod_pooled_flag$DV * (sim_mod_pooled_flag$PREDBIN / sim_mod_pooled_flag$PRED)
    
    # Order data
    sim_mod_pooled_flag <- sim_mod_pooled_flag[order(sim_mod_pooled_flag$replicate, sim_mod_pooled_flag$ID, sim_mod_pooled_flag$TAD),]
    
    # Initialize storage for PI data
    sim_PI_mod_pooled_flag <- NULL
    
    # Loop through replicates
    for (i in Replicate_vector_mod_pooled) {
      sim_vpc_ci_mod_pooled_flag <- sim_mod_pooled_flag %>%
        dplyr::filter(replicate %in% i) %>%
        group_by(BIN) %>%
        summarize(
          C_median = median(PCDV),
          C_lower = quantile(PCDV, perc_PI[1]),
          C_upper = quantile(PCDV, perc_PI[2])
        ) %>%
        mutate(replicate = i)
      
      sim_PI_mod_pooled_flag <- rbind(sim_PI_mod_pooled_flag, sim_vpc_ci_mod_pooled_flag)
    }
    
    # Calculate CI
    sim_CI_mod_pooled_flag <- sim_PI_mod_pooled_flag %>%
      group_by(BIN) %>%
      summarize(
        C_median_CI_lwr = quantile(C_median, perc_CI[1]),
        C_median_CI_upr = quantile(C_median, perc_CI[2]),
        C_low_lwr = quantile(C_lower, perc_CI[1]),
        C_low_upr = quantile(C_lower, perc_CI[2]),
        C_up_lwr = quantile(C_upper, perc_CI[1]),
        C_up_upr = quantile(C_upper, perc_CI[2])
      )
    
    # Assign NA values
    sim_CI_mod_pooled_flag$x1 <- NA
    sim_CI_mod_pooled_flag$x2 <- NA
    
    # Assign BIN values based on bin_times
    for (i in 1:Number_of_bins) {
      sim_CI_mod_pooled_flag$x1[sim_CI_mod_pooled_flag$BIN == i] <- bin_times[i]
      sim_CI_mod_pooled_flag$x2[sim_CI_mod_pooled_flag$BIN == i] <- bin_times[i + 1]
    }
    
    # dplyr::filter Obs_mod_pooled for the current hospital
    Obs_mod_pooled_flag <- Obs_mod_pooled %>% 
      dplyr::filter(CRRT == crrt)
    
    # Merge with PRED_BIN_mod_pooled_flag
    Obs_mod_pooled_flag <- merge(Obs_mod_pooled_flag, PRED_BIN_mod_pooled_flag, by = 'BIN')
    
    # Calculate PCDV
    Obs_mod_pooled_flag$PCDV <- Obs_mod_pooled_flag$DV * (Obs_mod_pooled_flag$PREDBIN / Obs_mod_pooled_flag$PRED)
    
    # Summarize data
    obs_vpc_mod_pooled_flag <- Obs_mod_pooled_flag %>%
      group_by(BIN) %>% 
      summarize(
        C_median = median(PCDV, na.rm = TRUE),
        C_lower = quantile(PCDV, perc_PI[1], na.rm = TRUE),
        C_upper = quantile(PCDV, perc_PI[2], na.rm = TRUE)
      )
    
    bin_middle_mod_pooled_flag <- Obs_mod_pooled_flag %>%
      group_by(BIN) %>%
      summarize(bin_middle = median(TAD, na.rm = TRUE))
    
    bin_middle_mod_pooled_flag <- bin_middle_mod_pooled_flag$bin_middle
    
    obs_vpc_mod_pooled_flag$TAD <- bin_middle_mod_pooled_flag
    
    # Create VPC plot
    CI_VPC_lin_mod_pooled_flag <- ggplot() +
      geom_rect(
        data = sim_CI_mod_pooled_flag,
        mapping = aes(xmin = x1, xmax = x2, ymin = C_median_CI_lwr, ymax = C_median_CI_upr),
        fill = '#35b779',
        alpha = 0.25
      ) +
      geom_rect(
        data = sim_CI_mod_pooled_flag,
        mapping = aes(xmin = x1, xmax = x2, ymin = C_low_lwr, ymax = C_low_upr),
        fill = '#31688e',
        alpha = 0.25
      ) +
      geom_rect(
        data = sim_CI_mod_pooled_flag,
        mapping = aes(xmin = x1, xmax = x2, ymin = C_up_lwr, ymax = C_up_upr),
        fill = '#31688e',
        alpha = 0.25
      ) +
      geom_line(
        data = obs_vpc_mod_pooled_flag,
        aes(TAD, C_median),
        col = '#440154',
        linetype = 'dashed',
        linewidth = 0.8
      ) +
      geom_line(
        data = obs_vpc_mod_pooled_flag,
        aes(TAD, C_lower),
        col = '#440154',
        linetype = 'dashed',
        linewidth = 0.8
      ) +
      geom_line(
        data = obs_vpc_mod_pooled_flag,
        aes(TAD, C_upper),
        col = '#440154',
        linetype = 'dashed',
        linewidth = 0.8
      ) +
      geom_point(
        data = Obs_mod_pooled_flag,
        aes(x = TAD, y = PCDV),
        color = '#440154',
        alpha = 0.3
      ) +
      scale_x_continuous(
        limits = c(0, 26),
        breaks = seq(0, 25, by = 5),
        name = "Time since last dose (hours)",
        expand = c(0.01, 0)
      ) +
      scale_y_continuous(
        limits = c(0, 80),
        breaks = seq(0, 80, by = 10),
        name = "Prediction-corrected fluconazole concentration (mg/L)",
        expand = c(0.01, 0)
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "Helvetica"), # following Nature: label: 8 pt, other text: 7 pt, min: 5 pt
            axis.title = element_text(size = 7, family = "Helvetica"),
            axis.text = element_text(size = 7, family = "Helvetica"),
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", size = 0.6, color = '#31688e')
    
    # Save the plot in the list with the corresponding title
    crrt_plots[[as.character(crrt)]] <- CI_VPC_lin_mod_pooled_flag +
      ggtitle(crrt_titles[crrt + 1])  
  }
  
  
  # Combine all the plots
  combined_plot <- cowplot::plot_grid(
    plotlist = lapply(names(crrt_plots), function(crrt) {
      return(crrt_plots[[crrt]])
    }),
    nrow = 1,
    labels = c('a', 'b'),
    label_size = 8
  )
  
  # Save the combine plot
  setwd("/lustre1/scratch/357/vsc35700/Fluco_revised/Plots/VPC/")
  
  ggsave("vpc_per_crrt_status.svg",
         combined_plot,
         dpi = 300, 
         width = 19, 
         height = 19,
         unit = "cm")
  
  ggsave("vpc_per_crrt_status.JPEG",
         combined_plot,
         dpi = 300, 
         width = 19, 
         height = 19,
         unit = "cm")
  
  ggsave("vpc_per_crrt_status.EPS",
         combined_plot,
         dpi = 300, 
         width = 19, 
         height = 19,
         unit = "cm")




  