# ---------------------------
# Purpose of script: Look at correlation of parameters of different TPC models
#
# What this script does:
# 1. Fits all TPC models to a all the curves from Chlorella TPC
# 2. Fits them with and without upper and lower limits
# 3. Saves them out
#
# Author: Dr. Daniel Padfield
#
# Date Created: 2024-05-29
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
librarian::shelf(tidyverse, rTPC, nls.multstart, furrr, progressr, microbenchmark)

## ---------------------------

# read in Chlorella TPC data
data("chlorella_tpc")

# fit a few models to the data

# load in data
data("chlorella_tpc")

d <- chlorella_tpc

# compare future_map with map ####

# compare using nls.multstart with 2 models and 10 curves
check_purrr <- microbenchmark(
  purrr = filter(d, curve_id <= 10) %>%
    nest(., data = c(temp, rate)) %>%
    mutate(beta = map(data, possibly(~nls_multstart(rate~beta_2012(temp = temp, a, b, c, d, e),
                                                    data = .x,
                                                    iter = c(6,6,6,6,6),
                                                    start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') * 0.5,
                                                    start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') * 1.5,
                                                    lower = get_lower_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                                    upper = get_upper_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                                    supp_errors = 'Y',
                                                    convergence_count = FALSE)), NA),
           boatman = map(data, possibly(~nls_multstart(rate~boatman_2017(temp = temp, rmax, tmin, tmax, a,b),
                                                       data = .x,
                                                       iter = c(5,5,5,5,5),
                                                       start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') * 0.5,
                                                       start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') * 1.5,
                                                       lower = get_lower_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                                       upper = get_upper_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                                       supp_errors = 'Y',
                                                       convergence_count = FALSE)), NA)),
  times = 10
)

check_purrr

check_furrr <- microbenchmark(
  furrr = {
  plan(multisession, workers = 3)
  
  filter(d, curve_id <= 10) %>%
    nest(., data = c(temp, rate)) %>%
    mutate(beta = future_map(data, possibly(~nls_multstart(rate~beta_2012(temp = temp, a, b, c, d, e),
                                                                     data = .x,
                                                                     iter = c(6,6,6,6,6),
                                                                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') * 0.5,
                                                                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') * 1.5,
                                                                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                                                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                                                     supp_errors = 'Y',
                                                                     convergence_count = FALSE, p = p)), NA),
           boatman = future_map(data, possibly(~nls_multstart(rate~boatman_2017(temp = temp, rmax, tmin, tmax, a,b),
                                                                        data = .x,
                                                                        iter = c(5,5,5,5,5),
                                                                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') * 0.5,
                                                                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') * 1.5,
                                                                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                                                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                                                        supp_errors = 'Y',
                                                                        convergence_count = FALSE, p = p)), NA))},
  times = 10
)

check_furrr
check_purrr