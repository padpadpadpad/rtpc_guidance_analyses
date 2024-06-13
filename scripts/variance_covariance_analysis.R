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


d_fits <- readRDS('data/chlorella_allTPCfits_nolimits.rds') %>%
  pivot_longer(c(beta:weibull), names_to = 'model', values_to = 'fit')

# get_correlation matrix for each model
d_cor <- d_fits %>%
  filter(!is.na(fit)) %>%
  mutate(cor_df = map(fit, get_cor_df)) %>%
  select(curve_id, model, cor_df) %>%
  unnest(cor_df)

head(d_cor)

# plot correlation matrices
d_cor <- mutate(d_cor, param_comp = paste(param1, param2, sep = ' - '))
d_cor_label <- group_by(d_cor, model, param_comp) %>%
  summarise(correlation = max(correlation), .groups = 'drop')

d_cor %>% ggplot(aes(model, correlation, group = param_comp)) +
  ylim(c(-1, 2)) +
  geom_hline(aes(yintercept = 1)) +
  geom_hline(aes(yintercept = -1)) +
  geom_point(aes(col = param_comp), position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.05), shape = 21, fill = 'white', alpha = 0.5, show.legend = FALSE) +
  geom_text(aes(x = model, y = correlation+(0.1), label = param_comp), position = position_dodge(width = 0.5), vjust = 0.4, hjust = 0, angle = 90, d_cor_label, size = MicrobioUoE::pts(6)) +
  theme_bw() +
  scale_color_manual(values = rep('black', times = length(unique(d_cor$param_comp)))) +
  facet_wrap(~model, scales = 'free_x')

ggsave('plots/correlation_matrix.png', width = 14, height = 10)

# grab an example fit
filter(d_fits, model == 'flinn') %>%
  pull(fit) %>%
  .[[1]]
