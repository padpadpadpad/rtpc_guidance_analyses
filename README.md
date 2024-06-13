# Thermal performance curves: A practical guide for choosing which models to fit and doing model selection

Thermal performance curves (TPCs) are routinely used to understand how individuals, populations, and communities respond to climate change. Previously, the fitting of models to these data was difficult and bespoke solutions were created across research groups. There are now multiple workflows that exist to make it easier and more reproducible to fit TPCs, such as [rTPC](https://github.com/padpadpadpad/rTPC) and [bayestpc](https://github.com/johnwilliamsmithjr/bayesTPC). However, these workflows are not a panacea for poor experimental design, and careful consideration still needs to be taken to make sure enough data is collected, and during the analysis in terms of what models to fit and how to do model selection. This repository contains scripts and analyses that underpin a sensible set of criteria against which researchers use to choose which models to fit to their data, allowing them to use the models to best answer the questions they are interested in, while also returning the best model fits and parameter estimates. There is no one size fits all method for analysing TPCs, and we want to highlight problems and pitfalls so researchers can use current tools with greater awareness and confidence.

**This repository is in development and the code and analyses are being actively worked on.**

## Scripts

1. **furrr_purrr_comparison.R** - checks whether parallelising model fitting using **furrr::future_map()** speeds up model fitting. 

2. **check_rtpc_model_performance.R** - checks each model individually to see if **furrr::future_map()** speeds up model fitting. Also gives an inclination about which models take a much longer time to fit than others. Generally the more parameters, the longer it takes, especially when using the grid start approach in **nls.multstart**. 

3. **all_tpcs_chlorella_parallel.R** - trying to fit all the TPC models to all the curves in chlorella_tpc using **furrr::future_map()**. This does not work and gets stuck but cannot work out why.

4. **all_tpcs_chlorella.R** - fits all the TPC models to all the curves in chlorella_tpc using **map()**. Fits models with and without limits on parameters.

5. **compare_fits_limits.R** - first post fit analysis. Checks whether any fits with limits are up against either the upper or lower limit. This is a red flag that the model may not be optimal, but also that the data may not be able to give sensible parameter estimates.

6. **compare_fits.R** - how often do fits with limits return a better fit than fits without limits?

6. **variance_covariance_analysis.R** - look at the variance-covariance matrix of the parameters of each model. High covariance between parameters may indicate that the model is not identifiable.