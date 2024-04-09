## Conc Time plot -----------------------------------------

### Standard dosing -------------------------------------------

# Load the dataset
working.directory <- '/lustre1/scratch/357/vsc35700/Fluco_revised/Pop_sim/'
std_tablefile <- paste0(working.directory, '/pop_std_conc.1.npctab.dta') 
pop_dosing_std <- read_nonmem_table(std_tablefile)

# Extract only TIME and DV columns
DV_TIME_std <- pop_dosing_std %>% select(DV, TIME)
setwd('/lustre1/scratch/357/vsc35700/Fluconazole_project/Fluc_Sim_201023/')
write.csv(DV_TIME_std, "DV_TIME_std.csv",quote = F,row.names = FALSE)

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
setwd('/lustre1/scratch/357/vsc35700/Fluco_revised/Dose_finding_simulations/Pop_sim/')
write.csv(DV_TIME_opt, "DV_TIME_opt.csv",quote = F,row.names = FALSE)

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