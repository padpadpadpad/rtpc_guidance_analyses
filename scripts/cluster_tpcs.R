# ---------------------------
# Purpose of script: Cluster analysis of TPC models
#
# What this script does:
# 1. Creates predictions of TPC models
# 2. Creates distance matrix of all model predictions
# 3. Clusters models based on distance matrix
#
# Author: Dr. Daniel Padfield
#
# Date Created: 2024-06-14
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
librarian::shelf(tidyverse, rTPC, nls.multstart, widyr, cluster, NbClust, fpc, dendextend, factoextra)

## ---------------------------

# load in models
d_fits <- readRDS('data/chlorella_allTPCfits_nolimits.rds') %>%
  pivot_longer(beta:weibull, names_to = 'model', values_to = 'fit')

# make new data set from data that makes gives 100 points for temp between the max and min temp inside the data column
d_fits <- d_fits %>%
  filter(!is.na(fit)) %>%
  mutate(new_data = map(data, \(x) data.frame(temp = seq(min(x$temp), max(x$temp), length.out = 50))),
         preds = map2(fit, new_data, \(x, y) broom::augment(x, newdata = y)))

d_fits$preds[[1]]

scale_vec <- function(x){scale(x) %>% as.vector()}

# unnest the predictions for each model
d_preds <- d_fits %>%
  unnest(preds) %>%
  select(model, curve_id, temp, .fitted) %>%
  group_by(curve_id) %>%
  # create scaled fitted values for each curve
  mutate(scaled_fitted = scale_vec(.fitted))

# create distance matrix
# scale predictions before calculating distance for each curve
d_dist <- d_preds %>%
  pairwise_dist(., model, temp, scaled_fitted, method = 'euclidean') %>%
  ungroup()
  
# calculate average distance between models
d_avg <- d_dist %>%
  group_by(item1, item2) %>%
  summarise(avg_dist = median(distance), .groups = 'drop')

d_dist2 <- ungroup(d_preds) %>% 
  mutate(id = paste(model, curve_id, sep = '_')) %>%
  pairwise_dist(., id, temp, scaled_fitted, method = 'euclidean') %>%
  ungroup()

# change to matrix
d_matrix <- d_avg %>%
  pivot_wider(names_from = item2, values_from = avg_dist) %>%
  column_to_rownames('item1')
  # arrange columns alphabetically

d_matrix <- d_matrix[,order(colnames(d_matrix))] %>%
  as.matrix()

# check column names and row names are the same
all(colnames(d_matrix) == rownames(d_matrix))

# replace NAs with 0
d_matrix[is.na(d_matrix)] <- 0

# try and cluster models ####

# do K medoids clustering
pamfun = function(x, k){list(cluster = pam(x, k, cluster.only = TRUE))}
pam_clusters = clusGap(d_matrix, FUN = pamfun, K.max = 25, B = 300, verbose = TRUE)

d_gap <- data.frame(pam_clusters$Tab, k=1:nrow(pam_clusters$Tab)) %>%
  data.frame()

# it thinks that 25 clusters is best! This seems silly - something must be going wrong.
ggplot(d_gap, aes(k, gap)) +
  geom_line() +
  geom_point(size = 5) +
  geom_linerange(aes(ymin = gap - SE.sim, ymax = gap + SE.sim)) +
  theme_bw(base_size = 14) +
  labs(y = 'Gap statistic',
       x = 'Number of clusters')

# do hierarchical clustering ####

hier_clust1 <- agnes(d_matrix,  method = 'ward')
hier_clust2 <- agnes(d_matrix,  method = 'average')
plot(hier_clust1)
plot(hier_clust2)

hier_clusgap <- clusGap(d_matrix, FUN = factoextra::hcut, K.max = 25, B = 500, hc_func = 'agnes', hc_method = 'complete')
d_gap <- data.frame(hier_clusgap$Tab, k=1:nrow(hier_clusgap$Tab)) %>%
  data.frame()

# calculate the optimal number of clusters
maxSE(d_gap$gap, d_gap$SE.sim, method = 'firstSEmax')
maxSE(d_gap$gap, d_gap$SE.sim, method = 'Tibs2001SEmax')

# calculate the optimal number of clusters
maxSE(d_gap$gap, d_gap$SE.sim, method = 'firstSEmax')

ggplot(d_gap, aes(k, gap)) +
  geom_line() +
  geom_point(size = 5) +
  geom_linerange(aes(ymin = gap - SE.sim, ymax = gap + SE.sim)) +
  theme_bw(base_size = 14) +
  labs(y = 'Gap statistic',
       x = 'Number of clusters')

# do NbClustering on Ward method
hier_nbclust <- NbClust::NbClust(d_matrix, method = 'single', max.nc = 4)
summary(hier_nbclust)
# The TSS matrix is indefinite. There must be too many missing values. The index cannot be calculated.

x <- hclust(as.dist(d_matrix), method = 'complete')
plot(x)

# the NbClust method fails
# there must be something going wrong currently in the clustering but not really sure what