# ---------------------------
# Purpose of script: Look at correlation of TPC parameters
#
# What this script does:
# 1. Extracts the variance covariance matrix of TPC parameters
# 2. Looks at these across different models
# 3. Are there any patterns?
#
# Author: Dr. Daniel Padfield
#
# Date Created: 2024-06-11
#
# Copyright (c) Daniel Padfield, 2024
#
# ---------------------------
#
# Notes:
#
# ---------------------------

# if librarian is not installed, install it
if (!requireNamespace("librarian", quietly = TRUE)){
  install.packages("librarian")
}
# load packages
librarian::shelf(rTPC, nls.multstart, broom, tidyverse)

## ---------------------------

# extra functions ####

# write function to get correlation matrix between parameters of an nls fit
get_cor_df <- function(model){
  temp <- vcov(model) %>%
    cov2cor()
  temp[!lower.tri(temp)] <- NA
  temp <- data.frame(temp) %>%
    rownames_to_column(var = 'param1') %>%
    pivot_longer(cols = -param1, names_to = 'param2', values_to = 'correlation') %>%
    filter(!is.na(correlation))
  return(temp)
}