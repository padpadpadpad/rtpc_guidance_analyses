# ---------------------------
# Purpose of script: Benchmark time to fit each model in rTPC
#
# What this script does:
# 1. runs microbenchmark on each model in rTPC for the first 10 TPCs
# 2. aggregates the data together
# 3. saves it out
#
# Author: Dr. Daniel Padfield
#
# Date Created: 2024-06-12
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
librarian::shelf(tidyverse, nls.multstart, rTPC, microbenchmark, furrr)

## ---------------------------

# read in Chlorella TPC data
data("chlorella_tpc")

# fit a few models to the data

# load in data
data("chlorella_tpc")

d <- chlorella_tpc

# create dataframe of the first 10 models
d10 <- d %>% 
  filter(curve_id <= 10) %>%
  nest(data = c(temp, rate))

# beta
check_beta <- microbenchmark(
  beta_furrr = {
    plan(multisession, workers = 3)
    
    d10 = mutate(d10, beta = future_map(data, possibly(~nls_multstart(rate~beta_2012(temp = temp, a, b, c, d, e),
                                                                data = .x,
                                                                iter = c(5,5,5,5,5),
                                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') * 0.25,
                                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') * 1.75,
                                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                                                supp_errors = 'Y',
                                                                convergence_count = FALSE)), NA))},
  beta_purrr = {d10 = mutate(d10, beta = map(data, possibly(~nls_multstart(rate~beta_2012(temp = temp, a, b, c, d, e),
                                                                                  data = .x,
                                                                                  iter = c(5,5,5,5,5),
                                                                                  start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') * 0.25,
                                                                                  start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') * 1.75,
                                                                                  lower = get_lower_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                                                                  upper = get_upper_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                                                                  supp_errors = 'Y',
                                                                                  convergence_count = FALSE)), NA))},
  times = 1
)

plan(sequential)

check_beta

# check boatman with both furrr and purrr
check_boatman <- microbenchmark(
  boatman_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, boatman = future_map(data, possibly(~nls_multstart(rate~boatman_2017(temp = temp, rmax, tmin, tmax, a,b),
                                                                            data = .x,
                                                                            iter = c(4,4,4,4,4),
                                                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') *.25,
                                                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') *1.75 ,
                                                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                                                            supp_errors = 'Y',
                                                                            convergence_count = FALSE), NA)))},
  boatman_purrr = {
    d10 = mutate(d10, boatman = map(data, possibly(~nls_multstart(rate~boatman_2017(temp = temp, rmax, tmin, tmax, a,b),
                                                                            data = .x,
                                                                            iter = c(4,4,4,4,4),
                                                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') *.25,
                                                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') *1.75 ,
                                                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                                                            supp_errors = 'Y',
                                                                            convergence_count = FALSE), NA)))},
  times = 1
)

check_boatman

plan(sequential)

# check briere2 with both furrr and purrr
check_briere2 <- microbenchmark(
  briere2_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, briere2 = future_map(data, possibly(~nls_multstart(rate~briere2_1999(temp = temp, tmin, tmax, a,b),
                                                                            data = .x,
                                                                            iter = c(4,4,4,4),
                                                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') * 0.25,
                                                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') * 1.75,
                                                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                                                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                                                                            supp_errors = 'Y',
                                                                            convergence_count = FALSE), NA)))},
  briere2_purrr = {
    d10 = mutate(d10, briere2 = map(data, possibly(~nls_multstart(rate~briere2_1999(temp = temp, tmin, tmax, a,b),
                                                                            data = .x,
                                                                            iter = c(4,4,4,4),
                                                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') * 0.25,
                                                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') * 1.75,
                                                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                                                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                                                                            supp_errors = 'Y',
                                                                            convergence_count = FALSE), NA)))},
  times = 1
)

plan(sequential)

check_briere2

# check delong with both purrr and furrr on d_60
check_delong <- microbenchmark(
  delong_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, delong = future_map(data, ~nls_multstart(rate~delong_2017(temp = temp, c, eb, ef, tm, ehc),
                                                        data = .x,
                                                        iter = c(5,5,5,5,5),
                                                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'delong_2017') * 0.25,
                                                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'delong_2017') * 1.75,
                                                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'delong_2017'),
                                                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'delong_2017'),
                                                        supp_errors = 'Y',
                                                        convergence_count = FALSE)))},
  delong_purrr = {
    d10 = mutate(d10, delong = map(data, ~nls_multstart(rate~delong_2017(temp = temp, c, eb, ef, tm, ehc),
                                                        data = .x,
                                                        iter = c(5,5,5,5,5),
                                                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'delong_2017') * 0.25,
                                                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'delong_2017') * 1.75,
                                                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'delong_2017'),
                                                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'delong_2017'),
                                                        supp_errors = 'Y',
                                                        convergence_count = FALSE)))},
  times = 1
)

plan(sequential)

# check deutsch_2008 on d_60 using both purrr and furrr
check_deutsch <- microbenchmark(
  deutsch_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, deutsch = future_map(data, possibly(~nls_multstart(rate~deutsch_2008(temp = temp, rmax, topt, ctmax, a),
                                                        data = .x,
                                                        iter = c(4,4,4,4),
                                                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'deutsch_2008') * 0.25,
                                                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'deutsch_2008') * 1.75,
                                                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'deutsch_2008'),
                                                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'deutsch_2008'),
                                                        supp_errors = 'Y',
                                                        convergence_count = FALSE), NA)))},
  deutsch_purrr = {
    d10 = mutate(d10, deutsch = map(data, possibly(~nls_multstart(rate~deutsch_2008(temp = temp, rmax, topt, ctmax, a),
                                                        data = .x,
                                                        iter = c(4,4,4,4),
                                                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'deutsch_2008') * 0.25,
                                                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'deutsch_2008') * 1.75,
                                                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'deutsch_2008'),
                                                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'deutsch_2008'),
                                                        supp_errors = 'Y',
                                                        convergence_count = FALSE), NA)))},
  times = 1
)

plan(sequential)

# check flinn_1991, it has parameters temp, a, b, c, using both purrr and furrr on d_10
check_flinn <- microbenchmark(
  flinn_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, flinn = future_map(data, ~nls_multstart(rate~flinn_1991(temp = temp, a, b, c),
                                                        data = .x,
                                                        iter = c(4,4,4),
                                                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') * 0.25,
                                                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') * 1.75,
                                                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                                                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                                                        supp_errors = 'Y',
                                                        convergence_count = FALSE)))},
  flinn_purrr = {
    d10 = mutate(d10, flinn = map(data, ~nls_multstart(rate~flinn_1991(temp = temp, a, b, c),
                                                        data = .x,
                                                        iter = c(4,4,4),
                                                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') * 0.25,
                                                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') * 1.75,
                                                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                                                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                                                        supp_errors = 'Y',
                                                        convergence_count = FALSE)))},
  times = 1
)

plan(sequential)

# check gaussian_1987, it has parameters temp, rmax, topt, a using both purrr and furrr on d_60
check_gaussian <- microbenchmark(
  gaussian_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, gaussian = future_map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                                            data = .x,
                                                            iter = c(4,4,4),
                                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') * 0.25,
                                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') * 1.75,
                                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                                            supp_errors = 'Y',
                                                            convergence_count = FALSE)))},
  gaussian_purrr = {
    d10 = mutate(d10, gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                                            data = .x,
                                                            iter = c(4,4,4),
                                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') * 0.25,
                                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') * 1.75,
                                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                                            supp_errors = 'Y',
                                                            convergence_count = FALSE)))},
  times = 1
)

plan(sequential)

# check hinshelwood_1947, it has parameters a, e, b, eh using furrr and purrr on d_60
check_hinshelwood <- microbenchmark(
  hinshelwood_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, hinshelwood = future_map(data, ~nls_multstart(rate~hinshelwood_1947(temp = temp, a, e, b, eh),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'hinshelwood_1947') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'hinshelwood_1947') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'hinshelwood_1947'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'hinshelwood_1947'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  hinshelwood_purrr = {
    d10 = mutate(d10, hinshelwood = map(data, ~nls_multstart(rate~hinshelwood_1947(temp = temp, a, e, b, eh),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'hinshelwood_1947') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'hinshelwood_1947') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'hinshelwood_1947'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'hinshelwood_1947'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  times = 1
)

plan(sequential)

# check joehnk_2008, parameters rmax, topt, a, b, c using furrr and purrr on d_60
check_joehnk <- microbenchmark(
  joehnk_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, joehnk = future_map(data, ~nls_multstart(rate~joehnk_2008(temp = temp, rmax, topt, a, b, c),
                                                        data = .x,
                                                        iter = c(4,4,4,4,4),
                                                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') * 0.25,
                                                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') * 1.75,
                                                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'joehnk_2008'),
                                                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'joehnk_2008'),
                                                        supp_errors = 'Y',
                                                        convergence_count = FALSE)))},
  joehnk_purrr = {
    d10 = mutate(d10, joehnk = map(data, ~nls_multstart(rate~joehnk_2008(temp = temp, rmax, topt, a, b, c),
                                                        data = .x,
                                                        iter = c(4,4,4,4,4),
                                                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') * 0.25,
                                                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') * 1.75,
                                                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'joehnk_2008'),
                                                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'joehnk_2008'),
                                                        supp_errors = 'Y',
                                                        convergence_count = FALSE)))},
  times = 1
)

plan(sequential)

# check johnsonlewin_1946, params r0, e, eh, topt using furrr and purrr on d_60
check_johsonlewin <- microbenchmark(
  johnsonlewin_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, johnsonlewin = future_map(data, ~nls_multstart(rate~johnsonlewin_1946(temp = temp, r0, e, eh, topt),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'johnsonlewin_1946'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'johnsonlewin_1946'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  johnsonlewin_purrr = {
    d10 = mutate(d10, johnsonlewin = map(data, ~nls_multstart(rate~johnsonlewin_1946(temp = temp, r0, e, eh, topt),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'johnsonlewin_1946'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'johnsonlewin_1946'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  times = 1
)

plan(sequential)

# check kamykowski_1985, params tmin, tmax, a, b, c using furrr and purrr on d_60
check_kamykoswki <- microbenchmark(
  kamykowski_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, kamykowski = future_map(data, ~nls_multstart(rate~kamykowski_1985(temp = temp, tmin, tmax, a, b, c),
                                                               data = .x,
                                                               iter = c(4,4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'kamykowski_1985'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'kamykowski_1985'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  kamykowski_purrr = {
    d10 = mutate(d10, kamykowski = map(data, ~nls_multstart(rate~kamykowski_1985(temp = temp, tmin, tmax, a, b, c),
                                                               data = .x,
                                                               iter = c(4,4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'kamykowski_1985'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'kamykowski_1985'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  times = 1
)

plan(sequential)

# check lactin2_1995, params a, b, tmax, delta_t using furrr and purrr on d_60
check_lactin2 <- microbenchmark(
  lactin2_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, lactin2 = future_map(data, ~nls_multstart(rate~lactin2_1995(temp = temp, a, b, tmax, delta_t),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  lactin2_purrr = {
    d10 = mutate(d10, lactin2 = map(data, ~nls_multstart(rate~lactin2_1995(temp = temp, a, b, tmax, delta_t),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  times = 1
)
  
plan(sequential)

# check lrf_1991, params rmax, topt, tmin, tmax using furrr and purrr on d_60
check_lrf <- microbenchmark(
  lrf_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, lrf = future_map(data, ~nls_multstart(rate~lrf_1991(temp = temp, rmax, topt, tmin, tmax),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lrf_1991') * 0.25,
                                                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lrf_1991') * 1.75,
                                                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'lrf_1991'),
                                                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'lrf_1991'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)))},
  lrf_purrr = {
    d10 = mutate(d10, lrf = map(data, ~nls_multstart(rate~lrf_1991(temp = temp, rmax, topt, tmin, tmax),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lrf_1991') * 0.25,
                                                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lrf_1991') * 1.75,
                                                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'lrf_1991'),
                                                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'lrf_1991'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)))},
  times = 1
)

plan(sequential)

# check modifiedgauassian_2006, params rmax, topt, a, b using furrr and purrr on d_60
check_modifiedgaussian <- microbenchmark(
  modifiedgaussian_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, modifiedgaussian = future_map(data, ~nls_multstart(rate~modifiedgaussian_2006(temp = temp, rmax, topt, a, b),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  modifiedgaussian_purrr = {
    d10 = mutate(d10, modifiedgaussian = map(data, ~nls_multstart(rate~modifiedgaussian_2006(temp = temp, rmax, topt, a, b),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  times = 1
)

plan(sequential)

# check oneill_1972, params rmax, ctmax, topt, q10 using furrr and purrr on d_60
check_oneill <- microbenchmark(
  oneill_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, oneill = future_map(data, ~nls_multstart(rate~oneill_1972(temp = temp, rmax, ctmax, topt, q10),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  oneill_purrr = {
    d10 = mutate(d10, oneill = map(data, ~nls_multstart(rate~oneill_1972(temp = temp, rmax, ctmax, topt, q10),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  times = 1
)

plan(sequential)

# check pawar_2018, params r_tref, e, eh, topt, tref = 15 using both furrr and purr on d_60
check_pawar <- microbenchmark(
  pawar_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, pawar = future_map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref, e, eh, topt, tref = 15),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'pawar_2018') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'pawar_2018') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'pawar_2018'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'pawar_2018'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  pawar_purrr = {
    d10 = mutate(d10, pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref, e, eh, topt, tref = 15),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'pawar_2018') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'pawar_2018') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'pawar_2018'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'pawar_2018'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  times = 1
)

plan(sequential)

# check quadratic_2008, params rmax, a, b c using both furrr and purr on d_60
check_quadratic <- microbenchmark(
  quadratic_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, quadratic = future_map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, rmax, a, b, c),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  quadratic_purrr = {
    d10 = mutate(d10, quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, rmax, a, b, c),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  times = 1
)

plan(sequential)

# check ratkowsky_1983 params tmin tmax a b
check_ratkowsky <- microbenchmark(
  ratkowsky_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, ratkowsky = future_map(data, ~nls_multstart(rate~ratkowsky_1983(temp = temp, tmin, tmax, a, b),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  ratkowsky_purrr = {
    d10 = mutate(d10, ratkowsky = map(data, ~nls_multstart(rate~ratkowsky_1983(temp = temp, tmin, tmax, a, b),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  times = 1
)

plan(sequential)

# check rezende_2019 params q10,a,b,c
check_rezende <- microbenchmark(
  rezende_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, rezende = future_map(data, ~nls_multstart(rate~rezende_2019(temp = temp, q10, a, b, c),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  rezende_purrr = {
    d10 = mutate(d10, rezende = map(data, ~nls_multstart(rate~rezende_2019(temp = temp, q10, a, b, c),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  times = 1
)

plan(sequential)

# check sharpeschoolfull_1981 params r_tref, e, el, tl, eh, th, tref =15
check_ssfull <- microbenchmark(
  sharpeschoolfull_furrr = {
    plan(multisession, workers = 3)
    d10 = mutate(d10, sharpeschoolfull = future_map(data, ~nls_multstart(rate~sharpeschoolfull_1981(temp = temp, r_tref, e, el, tl, eh, th, tref = 15),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  sharpeschoolfull_purrr = {
    d10 = mutate(d10, sharpeschoolfull = map(data, ~nls_multstart(rate~sharpeschoolfull_1981(temp = temp, r_tref, e, el, tl, eh, th, tref = 15),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') * 0.25,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') * 1.75,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981'),
                                                               supp_errors = 'Y',
                                                               convergence_count = FALSE)))},
  times = 1
)

plan(sequential)

rbind(check_beta, check_boatman, check_briere2, check_delong, check_deutsch, check_flinn, check_gaussian, check_hinshelwood, check_joehnk)

, check_johnsonlewin, check_kamykoswki, check_lactin2, check_logan, check_lrf, check_quadratic, check_ratkowsky, check_rezende, check_ssfull)





