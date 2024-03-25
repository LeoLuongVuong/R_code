# In this script, I do simulations with my heavy datasets, utilising the HPC
# resources :D

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