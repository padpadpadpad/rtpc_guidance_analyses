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
librarian::shelf(tidyverse, rTPC, nls.multstart, furrr, progressr)

## ---------------------------

# read in Chlorella TPC data
data("chlorella_tpc")

# fit a few models to the data

# load in data
data("chlorella_tpc")

d <- chlorella_tpc

#-------------------------------------------#
# fit all models with and without limits ####
#-------------------------------------------#

# when scaling up our code to fit hundreds of models, its nice to have a progress bar
# using future_map will speed this up too
# https://furrr.futureverse.org/articles/progress.html
# specifically when using future_map with dplyr

# make nls_multstart progress function compatible with furrr
nls_multstart_progress <- function(formula, data = parent.frame(), iter, start_lower, 
                                   start_upper, supp_errors = c("Y", "N"), convergence_count = 100, 
                                   control, modelweights, p, ...){
  p()
  
  nls_multstart(formula = formula, data = data, iter = iter, start_lower = start_lower, 
                start_upper = start_upper, supp_errors = supp_errors, convergence_count = convergence_count, 
                control = control, modelweights = modelweights, ...)
}

# set number of workers
plan(multisession, workers = 4)

# set up custom progress bar
handlers(handler_progress(clear = FALSE,
                          format ="[:bar] :percent :elapsedfull"))
# setup progress bar
with_progress({
  
  p <- progressor(steps = 1*60)

  d_fits <- nest(d, data = c(temp, rate)) %>%
    mutate(delong = future_map(data, possibly(~nls_multstart_progress(rate~delong_2017(temp = temp, c, eb, ef, tm, ehc),
                                                           data = .x,
                                                           iter = 10,
                                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'delong_2017') -10,
                                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'delong_2017') * 1.5,
                                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'delong_2017'),
                                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'delong_2017'),
                                                           supp_errors = 'Y',
                                                           convergence_count = 100, p=p), NA)))
})

# find number of models being fitted in total
number_of_models <- 24
number_of_curves <- length(unique(d$curve_id))

# omit hinshelwood and delong as they appear to be slowing everything down

# set number of workers
plan(sequential)
plan(multisession, workers = 3)

# set up custom progress bar
handlers(handler_progress(clear = FALSE,
                          format ="[:bar] :percent :elapsedfull"))
# setup progress bar
with_progress({
  
  p <- progressor(steps = number_of_models*number_of_curves)
  
  d_fits <- nest(d, data = c(temp, rate)) %>%
    mutate(beta = future_map(data, possibly(~nls_multstart_progress(rate~beta_2012(temp = temp, a, b, c, d, e),
                                                             data = .x,
                                                             iter = 10,
                                                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') * 0.5,
                                                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') * 1.5,
                                                             lower = get_lower_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                                             upper = get_upper_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                                             supp_errors = 'Y',
                                                             convergence_count = 100, p=p), NA)),
           boatman = future_map(data, possibly(~nls_multstart_progress(rate~boatman_2017(temp = temp, rmax, tmin, tmax, a,b),
                                                                data = .x,
                                                                iter = 10,
                                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') * 0.5,
                                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') * 1.5,
                                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                                                supp_errors = 'Y',
                                                                convergence_count = 100, p=p), NA)),
           briere2 = future_map(data, possibly(~nls_multstart_progress(rate~briere2_1999(temp = temp, tmin, tmax, a,b),
                                                                data = .x,
                                                                iter = 10,
                                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') * 0.5,
                                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') * 1.5,
                                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                                                                supp_errors = 'Y',
                                                                convergence_count = 100, p=p), NA)),
           deutsch = future_map(data, possibly(~nls_multstart_progress(rate~deutsch_2008(temp = temp, rmax, topt, ctmax, a),
                                                                data = .x,
                                                                iter = 10,
                                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'deutsch_2008') * 0.5,
                                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'deutsch_2008') * 1.5,
                                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'deutsch_2008'),
                                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'deutsch_2008'),
                                                                supp_errors = 'Y',
                                                                convergence_count = 100, p=p), NA)),
           flinn = future_map(data, possibly(~nls_multstart_progress(rate~flinn_1991(temp = temp, a, b, c),
                                                              data = .x,
                                                              iter = 10,
                                                              start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') * 0.5,
                                                              start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') * 1.5,
                                                              lower = get_lower_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                                                              upper = get_upper_lims(.x$temp, .x$rate, model_name = 'flinn_1991'),
                                                              supp_errors = 'Y',
                                                              convergence_count = 100, p=p), NA)),
           gaussian = future_map(data, possibly(~nls_multstart_progress(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                                                 data = .x,
                                                                 iter = 10,
                                                                 start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') * 0.5,
                                                                 start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') * 1.5,
                                                                 lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                                                 upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                                                 supp_errors = 'Y',
                                                                 convergence_count = 100, p=p), NA)),
           joehnk = future_map(data, possibly(~nls_multstart_progress(rate~joehnk_2008(temp = temp, rmax, topt, a, b, c),
                                                               data = .x,
                                                               iter = 10,
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') * 0.5,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') * 1.5,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'joehnk_2008'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'joehnk_2008'),
                                                               supp_errors = 'Y',
                                                               convergence_count = 100, p=p), NA)),
           johnson_lewin = future_map(data, possibly(~nls_multstart_progress(rate~ johnsonlewin_1946(temp = temp, r0, e, eh, topt),
                                                                      data = .x,
                                                                      iter = 10,
                                                                      start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') * 0.5,
                                                                      start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') * 1.5,
                                                                      lower = get_lower_lims(.x$temp, .x$rate, model_name = 'johnsonlewin_1946'),
                                                                      upper = get_upper_lims(.x$temp, .x$rate, model_name = 'johnsonlewin_1946'),
                                                                      supp_errors = 'Y',
                                                                      convergence_count = 100, p=p), NA)),
           kamykowski = future_map(data, possibly(~nls_multstart_progress(rate~kamykowski_1985(temp = temp, tmin, tmax, a, b, c),
                                                                   data = .x,
                                                                   iter = 10,
                                                                   start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') * 0.5,
                                                                   start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') * 1.5,
                                                                   lower = get_lower_lims(.x$temp, .x$rate, model_name = 'kamykowski_1985'),
                                                                   upper = get_upper_lims(.x$temp, .x$rate, model_name = 'kamykowski_1985'),
                                                                   supp_errors = 'Y',
                                                                   convergence_count = 100, p=p), NA)),
           lactin2 = future_map(data, possibly(~nls_multstart_progress(rate~lactin2_1995(temp = temp, a, b, tmax, delta_t),
                                                                data = .x,
                                                                iter = 10,
                                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') * 0.5,
                                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') * 1.5,
                                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
                                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'lactin2_1995'),
                                                                supp_errors = 'Y',
                                                                convergence_count = 100, p=p), NA)),
           lrf = future_map(data, possibly(~nls_multstart_progress(rate~lrf_1991(temp = temp, rmax, topt, tmin, tmax),
                                                            data = .x,
                                                            iter = 10,
                                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lrf_1991') * 0.5,
                                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lrf_1991') * 1.5,
                                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'lrf_1991'),
                                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'lrf_1991'),
                                                            supp_errors = 'Y',
                                                            convergence_count = 100, p=p), NA)),
           modifiedgaussian = future_map(data, possibly(~nls_multstart_progress(rate~modifiedgaussian_2006(temp = temp, rmax, topt, a, b),
                                                                         data = .x,
                                                                         iter = 10,
                                                                         start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') * 0.5,
                                                                         start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') * 1.5,
                                                                         lower = get_lower_lims(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006'),
                                                                         upper = get_upper_lims(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006'),
                                                                         supp_errors = 'Y',
                                                                         convergence_count = 100, p=p), NA)),
           oneill = future_map(data, possibly(~nls_multstart_progress(rate~oneill_1972(temp = temp, rmax, ctmax, topt, q10),
                                                               data = .x,
                                                               iter = 10,
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') * 0.5,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') * 1.5,
                                                               lower = get_lower_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                                                               upper = get_upper_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                                                               supp_errors = 'Y',
                                                               convergence_count = 100, p=p), NA)),
           pawar = future_map(data, possibly(~nls_multstart_progress(rate~pawar_2018(temp = temp, r_tref, e, eh, topt, tref = 15),
                                                              data = .x,
                                                              iter = 10,
                                                              start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'pawar_2018') * 0.5,
                                                              start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'pawar_2018') * 1.5,
                                                              lower = get_lower_lims(.x$temp, .x$rate, model_name = 'pawar_2018'),
                                                              upper = get_upper_lims(.x$temp, .x$rate, model_name = 'pawar_2018'),
                                                              supp_errors = 'Y',
                                                              convergence_count = 100, p=p), NA)),
           quadratic = future_map(data, possibly(~nls_multstart_progress(rate~quadratic_2008(temp = temp, a, b, c),
                                                                  data = .x,
                                                                  iter = 10,
                                                                  start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') * 0.5,
                                                                  start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') * 1.5,
                                                                  lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                                                  upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                                                  supp_errors = 'Y',
                                                                  convergence_count = 100, p=p), NA)),
           ratkowsky = future_map(data, possibly(~nls_multstart_progress(rate~ratkowsky_1983(temp = temp, tmin, tmax, a, b),
                                                                  data = .x,
                                                                  iter = 10,
                                                                  start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') * 0.5,
                                                                  start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') * 1.5,
                                                                  lower = get_lower_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
                                                                  upper = get_upper_lims(.x$temp, .x$rate, model_name = 'ratkowsky_1983'),
                                                                  supp_errors = 'Y',
                                                                  convergence_count = 100, p=p), NA)),
           rezende = future_map(data, possibly(~nls_multstart_progress(rate~rezende_2019(temp = temp, q10, a,b,c),
                                                                data = .x,
                                                                iter = 10,
                                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') * 0.5,
                                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') * 1.5,
                                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
                                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
                                                                supp_errors = 'Y',
                                                                convergence_count = 100, p=p), NA)),
           sharpeschoolfull = future_map(data, possibly(~nls_multstart_progress(rate~sharpeschoolfull_1981(temp = temp, r_tref,e,el,tl,eh,th, tref = 15),
                                                                         data = .x,
                                                                         iter = 10,
                                                                         start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') * 0.5,
                                                                         start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') * 1.5,
                                                                         lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981'),
                                                                         upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981'),
                                                                         supp_errors = 'Y',
                                                                         convergence_count = 100, p=p), NA)),
           sharpeschoolhigh = future_map(data, possibly(~nls_multstart_progress(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                                                         data = .x,
                                                                         iter = 10,
                                                                         start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') * 0.5,
                                                                         start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') * 1.5,
                                                                         lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                                         upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                                         supp_errors = 'Y',
                                                                         convergence_count = 100, p=p), NA)),
           sharpeschoollow = future_map(data, possibly(~nls_multstart_progress(rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15),
                                                                        data = .x,
                                                                        iter = 10,
                                                                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') * 0.5,
                                                                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') * 1.5,
                                                                        lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
                                                                        upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981'),
                                                                        supp_errors = 'Y',
                                                                        convergence_count = 100, p=p), NA)),
           spain = future_map(data, possibly(~nls_multstart_progress(rate~spain_1982(temp = temp, a,b,c,r0),
                                                              data = .x,
                                                              iter = 10,
                                                              start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') * 0.5,
                                                              start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') * 1.5,
                                                              lower = get_lower_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                                                              upper = get_upper_lims(.x$temp, .x$rate, model_name = 'spain_1982'),
                                                              supp_errors = 'Y',
                                                              convergence_count = 100, p=p), NA)),
           thomas1 = future_map(data, possibly(~nls_multstart_progress(rate~thomas_2012(temp = temp, a,b,c,tref),
                                                                data = .x,
                                                                iter = 10,
                                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2012') * 0.5,
                                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2012') * 1.5,
                                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'thomas_2012'),
                                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'thomas_2012'),
                                                                supp_errors = 'Y',
                                                                convergence_count = 100, p=p), NA)),
           thomas2 = future_map(data, possibly(~nls_multstart_progress(rate~thomas_2017(temp = temp, a,b,c,d,e),
                                                                data = .x,
                                                                iter = 10,
                                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') * 0.5,
                                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') * 1.5,
                                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'thomas_2017'),
                                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'thomas_2017'),
                                                                supp_errors = 'Y',
                                                                convergence_count = 100, p=p), NA)),
           weibull = future_map(data, possibly(~nls_multstart_progress(rate~weibull_1995(temp = temp, a,topt,b,c),
                                                                data = .x,
                                                                iter = 10,
                                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') * 0.5,
                                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') * 1.5,
                                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'weibull_1995'),
                                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'weibull_1995'),
                                                                supp_errors = 'Y',
                                                                convergence_count = 100, p=p), NA)))
  
  
})

# save out all model fits
saveRDS(d_fits, 'data/chlorella_allTPCfits_withlimits.rds')






# fit with no limits

# set number of workers
plan(multisession, workers = 3)

# set up custom progress bar
handlers(handler_progress(clear = FALSE,
                          format ="[:bar] :percent :elapsedfull"))
# setup progress bar
with_progress({
  
  p <- progressor(steps = number_of_models*number_of_curves)

  d_fits2 <- nest(d, data = c(temp, rate)) %>%
    mutate(beta = future_map(data, possibly(~nls_multstart_progress(rate~beta_2012(temp = temp, a, b, c, d, e),
                                                             data = .x,
                                                             iter = 10,
                                                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') * 0.5,
                                                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') * 1.5,
                                                             supp_errors = 'Y',
                                                             convergence_count = 100, p=p)), NA),
           boatman = future_map(data, possibly(~nls_multstart_progress(rate~boatman_2017(temp = temp, rmax, tmin, tmax, a,b),
                                                                data = .x,
                                                                iter = c(5,5,5,5,5),
                                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') * 0.5,
                                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') * 1.5,
                                                                supp_errors = 'Y',
                                                                convergence_count = 100, p=p)), NA),
           briere2 = future_map(data, possibly(~nls_multstart_progress(rate~briere2_1999(temp = temp, tmin, tmax, a,b),
                                                                data = .x,
                                                                iter = c(4,4,4,4),
                                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') * 0.5,
                                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') * 1.5,
                                                                supp_errors = 'Y',
                                                                convergence_count = 100, p=p)), NA),
           delong = future_map(data, possibly(~nls_multstart_progress(rate~delong_2017(temp = temp, c, eb, ef, tm, ehc),
                                                               data = .x,
                                                               iter = c(6,6,6,6,6),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'delong_2017') -10,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'delong_2017') * 1.5,
                                                               supp_errors = 'Y',
                                                               convergence_count = 100, p=p)), NA),
           deutsch = future_map(data, possibly(~nls_multstart_progress(rate~deutsch_2008(temp = temp, rmax, topt, ctmax, a),
                                                                data = .x,
                                                                iter = c(4,4,4,4),
                                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'deutsch_2008') * 0.5,
                                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'deutsch_2008') * 1.5,
                                                                supp_errors = 'Y',
                                                                convergence_count = 100, p=p)), NA),
           flinn = future_map(data, possibly(~nls_multstart_progress(rate~flinn_1991(temp = temp, a, b, c),
                                                              data = .x,
                                                              iter = c(5,5,5),
                                                              start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') * 0.5,
                                                              start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'flinn_1991') * 1.5,
                                                              supp_errors = 'Y',
                                                              convergence_count = 100, p=p)), NA),
           gaussian = future_map(data, possibly(~nls_multstart_progress(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                                                 data = .x,
                                                                 iter = c(4,4,4),
                                                                 start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') * 0.5,
                                                                 start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') * 1.5,
                                                                 supp_errors = 'Y',
                                                                 convergence_count = 100, p=p)), NA),
           hinshelwood = future_map(data, possibly(~nls_multstart_progress(rate~hinshelwood_1947(temp = temp, a, e, b, eh),
                                                                    data = .x,
                                                                    iter = c(5,5,5,5),
                                                                    start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'hinshelwood_1947') * 0.5,
                                                                    start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'hinshelwood_1947') * 1.5,
                                                                    supp_errors = 'Y',
                                                                    convergence_count = 100, p=p)), NA),
           joehnk = future_map(data, possibly(~nls_multstart_progress(rate~joehnk_2008(temp = temp, rmax, topt, a, b, c),
                                                               data = .x,
                                                               iter = c(4,4,4,4, 4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') * 0.5,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'joehnk_2008') * 1.5,
                                                               supp_errors = 'Y',
                                                               convergence_count = 100, p=p)), NA),
           johnson_lewin = future_map(data, possibly(~nls_multstart_progress(rate~ johnsonlewin_1946(temp = temp, r0, e, eh, topt),
                                                                      data = .x,
                                                                      iter = c(4,4,4,4),
                                                                      start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') * 0.5,
                                                                      start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'johnsonlewin_1946') * 1.5,
                                                                      supp_errors = 'Y',
                                                                      convergence_count = 100, p=p)), NA),
           kamykowski = future_map(data, possibly(~nls_multstart_progress(rate~kamykowski_1985(temp = temp, tmin, tmax, a, b, c),
                                                                   data = .x,
                                                                   iter = c(4,4,4,4,4),
                                                                   start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') * 0.5,
                                                                   start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'kamykowski_1985') * 1.5,
                                                                   supp_errors = 'Y',
                                                                   convergence_count = 100, p=p)), NA),
           lactin2 = future_map(data, possibly(~nls_multstart_progress(rate~lactin2_1995(temp = temp, a, b, tmax, delta_t),
                                                                data = .x,
                                                                iter = c(5,5,5,5),
                                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') * 0.5,
                                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lactin2_1995') * 1.5,
                                                                supp_errors = 'Y',
                                                                convergence_count = 100, p=p)), NA),
           lrf = future_map(data, possibly(~nls_multstart_progress(rate~lrf_1991(temp = temp, rmax, topt, tmin, tmax),
                                                            data = .x,
                                                            iter = c(3,3,3,3),
                                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'lrf_1991') * 0.5,
                                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'lrf_1991') * 1.5,
                                                            supp_errors = 'Y',
                                                            convergence_count = 100, p=p)), NA),
           modifiedgaussian = future_map(data, possibly(~nls_multstart_progress(rate~modifiedgaussian_2006(temp = temp, rmax, topt, a, b),
                                                                         data = .x,
                                                                         iter = c(4,4,4,4),
                                                                         start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') * 0.5,
                                                                         start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'modifiedgaussian_2006') * 1.5,
                                                                         supp_errors = 'Y',
                                                                         convergence_count = 100, p=p)), NA),
           oneill = future_map(data, possibly(~nls_multstart_progress(rate~oneill_1972(temp = temp, rmax, ctmax, topt, q10),
                                                               data = .x,
                                                               iter = c(4,4,4,4),
                                                               start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') * 0.5,
                                                               start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') * 1.5,
                                                               supp_errors = 'Y',
                                                               convergence_count = 100, p=p), NA)),
           pawar = future_map(data, possibly(~nls_multstart_progress(rate~pawar_2018(temp = temp, r_tref, e, eh, topt, tref = 15),
                                                              data = .x,
                                                              iter = c(4,4,4,4),
                                                              start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'pawar_2018') * 0.5,
                                                              start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'pawar_2018') * 1.5,
                                                              supp_errors = 'Y',
                                                              convergence_count = 100, p=p)), NA),
           quadratic = future_map(data, possibly(~nls_multstart_progress(rate~quadratic_2008(temp = temp, a, b, c),
                                                                  data = .x,
                                                                  iter = c(4,4,4),
                                                                  start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') * 0.5,
                                                                  start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') * 1.5,
                                                                  supp_errors = 'Y',
                                                                  convergence_count = 100, p=p)), NA),
           ratkowsky = future_map(data, possibly(~nls_multstart_progress(rate~ratkowsky_1983(temp = temp, tmin, tmax, a, b),
                                                                  data = .x,
                                                                  iter = c(4,4,4,4),
                                                                  start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') * 0.5,
                                                                  start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'ratkowsky_1983') * 1.5,
                                                                  supp_errors = 'Y',
                                                                  convergence_count = 100, p=p)), NA),
           rezende = future_map(data, possibly(~nls_multstart_progress(rate~rezende_2019(temp = temp, q10, a,b,c),
                                                                data = .x,
                                                                iter = c(4,4,4,4),
                                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') * 0.5,
                                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') * 1.5,
                                                                supp_errors = 'Y',
                                                                convergence_count = 100, p=p)), NA),
           sharpeschoolfull = future_map(data, possibly(~nls_multstart_progress(rate~sharpeschoolfull_1981(temp = temp, r_tref,e,el,tl,eh,th, tref = 15),
                                                                         data = .x,
                                                                         iter = c(4,4,4,4,4,4),
                                                                         start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') * 0.5,
                                                                         start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolfull_1981') * 1.5,
                                                                         supp_errors = 'Y',
                                                                         convergence_count = 100, p=p)), NA),
           sharpeschoolhigh = future_map(data, possibly(~nls_multstart_progress(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                                                         data = .x,
                                                                         iter = c(4,4,4,4),
                                                                         start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') * 0.5,
                                                                         start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') * 1.5,
                                                                         supp_errors = 'Y',
                                                                         convergence_count = 100, p=p)), NA),
           sharpeschoollow = future_map(data, possibly(~nls_multstart_progress(rate~sharpeschoollow_1981(temp = temp, r_tref,e,el,tl, tref = 15),
                                                                        data = .x,
                                                                        iter = c(4,4,4,4),
                                                                        start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') * 0.5,
                                                                        start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoollow_1981') * 1.5,
                                                                        supp_errors = 'Y',
                                                                        convergence_count = 100, p=p)), NA),
           spain = future_map(data, possibly(~nls_multstart_progress(rate~spain_1982(temp = temp, a,b,c,r0),
                                                              data = .x,
                                                              iter = c(4,4,4,4),
                                                              start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') * 0.5,
                                                              start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'spain_1982') * 1.5,
                                                              supp_errors = 'Y',
                                                              convergence_count = 100, p=p)), NA),
           thomas1 = future_map(data, possibly(~nls_multstart_progress(rate~thomas_2012(temp = temp, a,b,c,tref),
                                                                data = .x,
                                                                iter = c(4,4,4,4),
                                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2012') * 0.5,
                                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2012') * 1.5,
                                                                supp_errors = 'Y',
                                                                convergence_count = 100, p=p)), NA),
           thomas2 = future_map(data, possibly(~nls_multstart_progress(rate~thomas_2017(temp = temp, a,b,c,d,e),
                                                                data = .x,
                                                                iter = c(3,3,3,3,3),
                                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') * 0.5,
                                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'thomas_2017') * 1.5,
                                                                supp_errors = 'Y',
                                                                convergence_count = 100, p=p)), NA),
           weibull = future_map(data, possibly(~nls_multstart_progress(rate~weibull_1995(temp = temp, a,topt,b,c),
                                                                data = .x,
                                                                iter = c(4,4,4,4),
                                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') * 0.5,
                                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') * 1.5,
                                                                supp_errors = 'Y',
                                                                convergence_count = 100, p=p), NA)))
  
})


# save out all model fits
saveRDS(d_fits2, 'data/chlorella_allTPCfits_nolimits.rds')

