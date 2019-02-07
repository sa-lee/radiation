# process_geo.R
# Prepare summary tables for GEO submission
load("data/counts.rda")
library(tidyverse)
assay <- counts
counts_mat <- assay$counts
# raw samples
samples <- dplyr::tiblble(sample_names = colnames(counts_mat)) %>%
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
                sample_title = paste(celltype, 
                                     dplyr::if_else(treatment == "0Gy", 
                                                    "untreated", "treated"),
                                     "technical replicate", replicate),
                organism = "Homo sapiens",
                source =  paste(celltype, "RNAseq radiotherapy timecourse"),
                processed_file = paste0(celltype, "_", treatment, "_", replicate, ".txt"),
                raw_file = stringr::str_replace(sample_names, ".sam$", "_R1.fastq.gz"),
                
                
  ) %>% 
  dplyr::select(-sample_meta)

readr::write_excel_csv(samples, here::here("geo_submission", "samples.txt"))

# proccessed counts
processed_counts <- readr::read_rds("geo_submission/processed_counts.rds")

for (i in seq_len(ncol(processed_counts))) {
  
  sample_name <- processed_counts$samples$sample_names[[i]]
  group <- processed_counts$samples$group[[i]]
  replicate <- processed_counts$samples$replicate[[i]]
  sample_counts <- processed_counts$counts[, sample_name]
  gene_counts <- tibble(entrez_gene_id = processed_counts$genes$GeneID,
                            gene_length = processed_counts$genes$Length,
                            read_count = sample_counts)
  readr::write_tsv(gene_counts, 
            path = here::here("geo_submission", 
                              "summarised_counts",
                               paste0(group, "_", replicate, ".txt"))
  )
}
