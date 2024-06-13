# ---------------------------
# Purpose of script: Create discovery curve of TPC models
#
# What this script does:
# 1. reads in bib file from Dimitrios' paper
# 2. make discovery curve to look at how models accumulated through time
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
librarian::shelf(tidyverse, bib2df)

## ---------------------------

# read in bib file
d <- bib2df('data/bibliography.bib') %>%
  janitor::clean_names()

# select year and create cumulative sum
d_year <- select(d, year) %>%
  mutate(year = as.numeric(year)) %>%
  arrange(year) %>%
  mutate(n = 1) %>%
  group_by(year) %>%
  summarise(n = sum(n), .groups = 'drop') %>%
  mutate(cumulative = cumsum(n))

# plot
ggplot(d_year, aes(x=year, y=cumulative)) +
  geom_line(size = 2, colour = 'dodgerblue4')+
  theme_bw(base_size = 14)+
  theme(legend.position = 'none') +
  labs(x = 'Year',
       y = 'Number of TPC models')
