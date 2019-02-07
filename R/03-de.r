library(tidyverse)
library(edgeR)

# DE via lmFit
vfit <- read_rds(here::here("data", "voom.rds"))

contr <- read_rds(here::here("data", "contrasts.rds"))

lm_fit <- lmFit(vfit)
lm_fit <- contrasts.fit(lm_fit, contrasts = contr)
lm_fit <- eBayes(lm_fit)

corrected_variance_plot <- data_frame(mean_signal = lm_fit$Amean,
                                      variance_signal = sqrt(lm_fit$sigma),
                                      trend = sqrt(lm_fit$s2.prior)) %>% 
  ggplot(., aes(x = mean_signal, y = variance_signal)) +
  geom_point() +
  geom_line(aes(x = mean_signal, y = trend), colour = "blue") +
  xlab("Average log-expression") +
  ylab("sqrt(sigma)")


# -- top Tables
dt_all <- decideTests(lm_fit, method = "global", p.value = 0.05 )
dt_by_contrast <- decideTests(lm_fit, method = "separate", p.value = 0.05)


top_tables <- lapply(seq_len(ncol(lm_fit)), 
                     function(i) 
                       as_tibble(topTable(lm_fit, 
                                              coef = i, 
                                              n = nrow(lm_fit))) %>%
                       mutate(contrast = colnames(lm_fit)[i])) %>% 
  bind_rows() %>%
  group_by(contrast) 

top_tables %>% 
  select(Chr, GeneID, Symbols, GeneName, contrast, logFC, adj.P.Val) %>% 
  arrange(adj.P.Val, desc(abs(logFC)), .by_group = TRUE) %>% 
  write_csv(here::here("results", "top_tables_by_contrast.csv"))


## - volcano plot
rvc_summary <- top_tables %>% 
  filter(contrast == "treatment10GY-treatment0GY") %>%
  mutate(de = as.factor(adj.P.Val < 0.05)) 

large_FC <- rvc_summary %>% filter(abs(logFC) > 1.0, adj.P.Val < 0.05)
volcano_rvc <- ggplot(data = rvc_summary, 
                      aes(x = logFC, y = -log10(P.Value), colour = de)) +
  geom_point(stroke = 0) +
  ggrepel::geom_text_repel(data = large_FC, 
                           aes(x = logFC, y = -log10(P.Value), 
                               colour = de, label = Symbols),
                           size = 3) +
  guides(colour = FALSE) + xlab("log Fold Change") + ylab("-log10(P Value)") +
  theme(axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8))

ggsave(here::here("figures", "figure-02.png"), volcano_rvc)


# -- de genes count
top_list_by_contrast <- top_tables %>% 
  mutate(updown = sign(logFC)) %>% 
  filter(adj.P.Val < 0.05) 

top_de_all <-top_list_by_contrast %>% 
  count()
top_de_updown <- top_list_by_contrast %>% 
  group_by(contrast, updown) %>% 
  count()