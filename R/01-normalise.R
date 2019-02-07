# Prepare DGEList and normalise counts using TMM
# Author: Stuart Lee
library(tidyverse)
library(edgeR)

# the counts object was prepared earlier by Matt
# this also loads an object that has three elements
# counts: a matrix of gene counts by samples
# annotation: a data.frame of gene annotations
# targets: a data.frame of sample metadata
load(here::here("data", "counts.rda"))

# correct sample PERI1_10Gy counts, this sample was split accross
# three lanes, so we merge its counts
peri1_lane_split <- stringr::str_detect(colnames(counts$counts), "PERI1_10GY")
total_peri1_counts <- rowSums(counts$counts[, peri1_lane_split])
# update counts matrix, with merged samples
counts_mat <- matrix(0L, 
                     nrow = nrow(counts$counts),
                     ncol = sum(!peri1_lane_split) + 1)
counts_mat[, seq_len(sum(!peri1_lane_split))] <- counts$counts[, !peri1_lane_split]
counts_mat[, ncol(counts_mat)] <- total_peri1_counts
# update column names
colnames(counts_mat) <- c(colnames(counts$counts)[!peri1_lane_split],
                          colnames(counts$counts)[peri1_lane_split][1L])
# filter genes based on cpm
gene_filter <- rowSums(edgeR::cpm(counts_mat) > 0.5) >= 3
counts_mat <- counts_mat[gene_filter, ]

# create colData
samples <- dplyr::tibble(sample_names = colnames(counts_mat)) %>%
  dplyr::mutate(sample_meta = stringr::str_split(sample_names, "_"),
                celltype =  purrr::map_chr(sample_meta, ~ .x[1]) %>% 
                  stringr::str_replace(., "1|2", ""),
                replicate = purrr::map_chr(sample_meta, ~ .x[1]) %>% 
                  stringr::str_extract(., "1|2") %>% 
                  as.integer(.),
                treatment = purrr::map_chr(sample_meta, ~.x[2]),
                lane = purrr::map_chr(sample_meta, ~.x[length(.x)]) %>%
                  stringr::str_replace_all(., "L00|\\.sam", "") %>% 
                  as.integer(.),
                library.size = colSums(counts_mat),
                norm.factors = 1,
                group = as.factor(paste0(celltype, "_", treatment))
  ) %>% 
  dplyr::select(-sample_meta)

# create rowData
geneanno_update <- geneanno[gene_filter, ] %>% 
  dplyr::select(-EntrezID) 

# now have all the required information for DGEList
elist <- edgeR::DGEList(counts = counts_mat,
                        samples = as.data.frame(samples),
                        group = samples$group,
                        genes =  as.data.frame(geneanno_update))

# perform TMM normalisation
elist <- calcNormFactors(elist, method = "TMM")

readr::write_rds(elist,
                 here::here("data", "processed_counts.rds"))
