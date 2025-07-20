library(dplyr)
source("utils/eval.R")
source("utils/mini_import.R")

# ---- epianeufinder ----
invisible(list2env(readRDS("intmd/Combined_outputs_2025-02-21.rds"),environment()))
wgsCall_gr <- map(wgs_call,\(w) w[,c("seqnames","start","end")])
CNV2eval <- map(combinedout$numbat_call[["DL3267"]],\(d) d %>% 
                    select(chr=CHROM,start=seg_start,
                           end=seg_end,
                           eventType=cnv_state_post))
CNV2eval$epianeufinder <- fread("benchmark/epianeufinder_DL3267_cnv_calls.tsv") %>% 
  select(-n_cells) %>% 
  mutate(seqnames=gsub("chr","",seqnames))

map(CNV2eval,\(n)eval_call(wgsCall_gr[["DL3267"]],n)) %>%
  bind_rows() %>% 
  mutate(mode = names(CNV2eval))

library(ggplot2)
library(tidyr)

# Assuming the previous code produces a data frame like this:
pr_df <- map(CNV2eval, \(n) eval_call(wgsCall_gr[["DL3267"]], n)) %>%
  bind_rows() %>%
  mutate(mode = names(CNV2eval))

# Reshape for plotting
pr_long <- pr_df %>%
  select(mode, precision, recall, f1) %>%
  pivot_longer(cols = c("precision", "recall", "f1"), names_to = "metric", values_to = "value")

# Plot
ggplot(pr_long, aes(x = mode, y = value, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Precision, Recall, F1 for DL3267", y = "Score", x = "Mode") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")

library(patchwork) # already loaded

p_f1 <- ggplot(filter(pr_long, metric == "f1"), aes(x = mode, y = value, fill = mode)) +
  geom_bar(stat = "identity") +
  labs(title = "F1 Score", y = "F1", x = NULL) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  guides(fill = guide_legend(title = "Mode"))

p_precision <- ggplot(filter(pr_long, metric == "precision"), aes(x = mode, y = value, fill = mode)) +
  geom_bar(stat = "identity") +
  labs(title = "Precision", y = "Precision", x = NULL) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  guides(fill = guide_legend(title = "Mode"))

p_recall <- ggplot(filter(pr_long, metric == "recall"), aes(x = mode, y = value, fill = mode)) +
  geom_bar(stat = "identity") +
  labs(title = "Recall", y = "Recall", x = NULL) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  guides(fill = guide_legend(title = "Mode"))

# Combine the three plots with a shared legend
(p_f1 / p_precision / p_recall) + plot_layout(guides = "collect")

# ---- CONGASP ----
library(purrr)
paired_results <- map(intersect(names(wgsCall_gr), names(numbat_call)), function(sample) {
  if (!"comb_bincnt" %in% names(numbat_call[[sample]])) return(NULL)
  congasp_file <- file.path('congasp/cnv_congasp', paste0(sample, '_cnv_segments_formatted.tsv'))
  if (!file.exists(congasp_file)) return(NULL)
  congasp_dat <- read.delim(congasp_file, stringsAsFactors = FALSE) %>%
    select(seqnames, start, end, cnv, cluster) %>%
    mutate(seqnames = gsub('chr', '', seqnames))
  numbat_dat <- numbat_call[[sample]][["comb_bincnt"]] %>%
    select(seqnames = CHROM, start = seg_start, end = seg_end, cnv = cnv_state_post) %>%
    mutate(seqnames = as.character(seqnames))
  # Evaluate byCNV=TRUE for numbat
  pr_numbat <- eval_call(wgsCall_gr[[sample]], numbat_dat, byCNV = TRUE) %>%
    filter(cnv %in% c('amp','del'))
  pr_numbat$method <- 'Numbat-multiome\nCombined Bin'; pr_numbat$sample <- sample
  # For CONGASP: evaluate per cluster, take max for each metric and cnv
  pr_congasp_all <- map(split(congasp_dat, congasp_dat$cluster), function(df) {
    eval_call(wgsCall_gr[[sample]], df %>% select(seqnames, start, end, cnv), byCNV = TRUE)
  }) %>% bind_rows() %>% filter(cnv %in% c('amp','del'))
  pr_congasp_max <- pr_congasp_all %>%
    group_by(cnv) %>%
    summarise(precision = max(precision, na.rm=TRUE),
              recall = max(recall, na.rm=TRUE),
              f1 = max(f1, na.rm=TRUE)) %>%
    ungroup()
  pr_congasp_max$method <- 'CONGAS+'; pr_congasp_max$sample <- sample
  list(congasp = pr_congasp_max, numbat = pr_numbat)
})
paired_all_df <- map_dfr(paired_results, ~bind_rows(.x), .id = NULL)
paired_all_long <- paired_all_df %>%
  select(sample, method, cnv, precision, recall, f1) %>%
  pivot_longer(cols = c('precision', 'recall', 'f1'), 
    names_to = 'metric', values_to = 'value')

# Plot: stacked precision, recall, F1 vertically using patchwork
library(ggplot2)
library(patchwork)
plot_metric <- function(met) {
  paired_all_long %>% filter(metric == met) %>%
    ggplot(aes(x = method, y = value, group = sample, color = sample)) +
    geom_boxplot(aes(group = method), width = 0.5, fill = 'gray90', color = 'black', alpha = 0.4, outlier.shape = NA) +
    geom_point(size = 2) +
    geom_line() +
    facet_wrap(~cnv) +
    labs(title = paste('Paired', toupper(met), 'for all samples: CONGAS+ vs Numbat-multiome Combined Bin'),
         y = toupper(met), x = NULL) +
    theme_minimal() +
    theme(legend.position = 'none')
}
# Plot: stacked precision, recall, F1 vertically using ggarrange
library(ggpubr)
p_f1 <- plot_metric('f1')
p_precision <- plot_metric('precision')
p_recall <- plot_metric('recall')
ggarrange(p_f1, p_precision, p_recall, ncol = 1, nrow = 3)
