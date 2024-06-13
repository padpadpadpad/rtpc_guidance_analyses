# ---------------------------
# Purpose of script: Look at whether TPC fits with limits are up against their lower and upper limits
#
# What this script does:
# 1. Extracts upper and lower limits for the TPC model fit
# 2. Extracts the coefficients for the TPC model fit
# 3. Compares coefficient values for the same model with and without limits
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
librarian::shelf(rTPC, broom, nls.multstart, tidyverse)

## ---------------------------

# extra functions ####

# function to grab lower limits and upper limits from a nls fit if present
# if not present return NA
get_fit_limits <- function(mod){
  
  #grab lower and upper limits if present
  temp_lower <- data.frame(lower_lim = mod$call$lower) %>%
    rownames_to_column(var = 'param')
  # if not present create an empty data frame with NA
  if(is.null(mod$call$lower)){
    temp_lower <- data.frame(param = names(coef(mod)), upper_lim = NA)
  }
  
  # grab upper limits if present
  temp_upper <- data.frame(upper_lim = mod$call$upper) %>%
    rownames_to_column(var = 'param')
  
  # if not present create an empty data frame with NA
  if(is.null(mod$call$upper)){
    temp_upper <- data.frame(param = names(coef(mod)), upper_lim = NA)
  }
  
  # join and return
  temp <- full_join(temp_lower, temp_upper, by = 'param')
  
  return(temp)
}

# read in models with limits ####
fits_with_limits <- readRDS('data/chlorella_allTPCfits_withlimits.rds') %>%
  pivot_longer(c(beta: weibull), names_to = 'model_name', values_to = 'fit') %>%
  mutate(., params = map(fit, calc_params, .progress = TRUE))

# count the number of NULL fits
fits_with_limits %>%
  summarise(n = sum(map_lgl(fit, is.null)))

# get limits and coefficients for all fits
# put in a single dataframe and unnest
d_check <- mutate(fits_with_limits, limits = map(fit, get_fit_limits),
                  coef = map(fit, broom::tidy),
                  combined = map2(coef, limits, ~ cbind(.x, .y))) %>%
  select(curve_id, model_name, combined) %>%
  unnest(combined)

# add column if any of the coefficients estimates are the same as lower or upper limits
d_check <- mutate(d_check, check = case_when(lower_lim == estimate ~ 'lower',
                                            upper_lim == estimate ~ 'upper',
                                            TRUE ~ 'none'))

# create dataframe with all params for each model with a check column
d_params <- d_check %>%
  select(model_name, param) %>%
  distinct() %>%
  group_by_all() %>%
  do(data.frame(check = c('lower', 'upper', 'none')))

# for each model count the number of each check
d_check2 <- mutate(d_check, check = as.factor(check)) %>%
  group_by(model_name, param, check) %>%
  count(.drop = FALSE)

d_check2 <- left_join(d_params, d_check2) %>%
  replace_na(list(n = 0))

# read in models with no limits ####
fits_no_limits <- readRDS('data/chlorella_allTPCfits_nolimits.rds') %>%
  pivot_longer(c(beta: weibull), names_to = 'model_name', values_to = 'fit') %>%
  mutate(., params = map(fit, calc_params, .progress = TRUE))

# count the number of NULL fits
fits_no_limits %>%
  summarise(n = sum(map_lgl(fit, is.null)))

# filter data for when fit is null
filter(fits_no_limits, map_lgl(fit, is.null))
# all delong

# compare the AIC scores for fits with and without limits
d_no_limits_aic <- filter(fits_no_limits, !map_lgl(fit, is.null)) %>%
  mutate(., aic = map(fit, AIC)) %>%
  select(curve_id, model_name, aic) %>%
  unnest(aic) %>%
  select(curve_id, model_name, aic_nolimit = aic)

d_limit_aic <- mutate(fits_with_limits, aic = map(fit, AIC)) %>%
  select(curve_id, model_name, aic) %>%
  unnest(aic) %>%
  select(curve_id, model_name, aic_limit = aic)

d_aic <- full_join(d_limit_aic, d_no_limits_aic, by = c('curve_id', 'model_name'))

# compare the coefficients for the fits with and without limits
d_no_limits_coefs <- mutate(fits_no_limits, coef = map(fit, broom::tidy)) %>%
  select(curve_id, model_name, coef) %>%
  unnest(coef) %>%
  select(curve_id, model_name, term, estimate_nolimit = estimate)

d_limit_coefs <- mutate(fits_with_limits, coef = map(fit, broom::tidy)) %>%
  select(curve_id, model_name, coef) %>%
  unnest(coef) %>%
  select(curve_id, model_name, term, estimate_limit = estimate)

# join the two dataframes
d_compare <- full_join(d_limit_coefs, d_no_limits_coefs, by = c('curve_id', 'model_name', 'term')) %>%
  left_join(., select(d_check, curve_id, model_name, term, check))

d_compare2 <- filter(d_compare, round(estimate_limit, 1) != round(estimate_nolimit, 1))
