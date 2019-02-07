# Exploratory data analysis with voom
# Author: Stuart Lee
library(tidyverse)
library(patchwork)
library(limma)
theme_set(theme_bw())
elist <- readr::read_rds(here::here("data", "processed_counts.rds"))

# -- MDS (with ggplot2)
mds_obj <- plotMDS(elist, ndim = 10, plot = FALSE) 

mds_data <- mds_obj$cmdscale.out %>% 
  as.data.frame(.) %>% 
  rownames_to_column(var = "sample_names") %>% 
  left_join(elist$samples)

# compute proportion of variance explained
mds_redo <- as.dist(mds_obj$distance.matrix) %>% 
  cmdscale(., k = 10, eig = TRUE) 

prop_variance <- tibble(dim = seq_len(10),
                        variance_explained = mds_redo$eig[seq_len(10)] / sum(mds_redo$eig))


# DE within contrasts
mds_lec <- plotMDS(elist[, grepl("LEC|HMEC", elist$samples$group)], 
                   ndim = 2, 
                   plot = FALSE)$cmdscale.out %>%
  as.data.frame(.) %>% 
  rownames_to_column(., var = "sample_names") %>% 
  left_join(elist$samples)

# plots 
bar_plot <- ggplot(prop_variance, aes(x = dim, y = variance_explained)) +
  geom_bar(stat = "identity") + 
  scale_x_continuous(breaks = seq_len(10), labels = seq_len(10)) +
  xlab("Dimension") + ylab("Proportion variance explained") +
  theme(axis.title.y = element_text(size = 8))

mds_by_cellline <- ggplot(mds_data, 
                          aes(x = V1, y = V2, colour = celltype, label = group)) +
  geom_point(stroke = 0) + 
  xlab("Leading logFC dim 1") + ylab("Leading logFC dim 2") +
  scale_color_brewer(palette = "Set2") +
  guides(colour = guide_legend("Cell line")) + 
  theme(legend.position="right", 
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))

mds_by_lec <- ggplot(mds_lec, 
                     aes(x = V1, y = V2, label = group, colour = as.factor(replicate))) +
  geom_point(stroke = 0) +
  ggrepel::geom_label_repel(size = 2) +
  guides(colour = guide_legend("replicate number")) +
  xlab("Leading logFC dim 1") + ylab("Leading logFC dim 2") +
  theme(legend.position = "bottom", 
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))

# voom
treatment <- as.factor(elist$samples$treatment)
celltype <- as.factor(elist$samples$celltype)
replicate <- as.factor(elist$samples$replicate)
design <- model.matrix(~0 + treatment + celltype + replicate)
#colnames(design) <- c(as.character(unique(group)), "replicate")

contr <-  makeContrasts(treatment2GYx5-treatment0GY,
                        treatment10GY-treatment2GYx5,
                        treatment10GY-treatment0GY,
                        levels = colnames(design))

readr::write_rds(
  contr,
  here::here("data", "contrasts.rds")
)

colnames(contr) <- gsub(" ", "", colnames(contr))

vfit <- voomWithQualityWeights(elist, 
                               plot = FALSE, 
                               design = design,
                               save.plot = TRUE)

readr::write_rds(vfit, 
                 here::here("data", "voom.rds"))

mean_variance_df <- tibble(mean_signal = vfit$voom.xy$x,
                               variance_signal = vfit$voom.xy$y,
                               loess_line_x = vfit$voom.line$x,
                               loess_line_y = vfit$voom.line$y)



voom_plot <- ggplot(mean_variance_df, aes(x = mean_signal, y = variance_signal)) +
  geom_point(stroke = 0) +
  geom_line(aes(x = loess_line_x, y = loess_line_y), colour = "blue") +
  xlab("log2(count size + 0.5)") + ylab("sqrt(standard deviation)") +
  theme(axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8))

# combine all figures
pl <- wrap_plots(
  bar_plot,
  mds_by_cellline,
  mds_by_lec,
  voom_plot,
  nrow = 2,
  ncol = 2,
  tag_level = "new") +
  plot_annotation(tag_levels = "A") 



ggsave(here::here("figures", "figure-01.png"), pl)
