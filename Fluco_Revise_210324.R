

# Call libraries ----------------------------------------------------------

library(tidyverse)
library(conflicted)

# Explanatory Data Analysis -----------------------------------------------

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



