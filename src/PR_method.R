library(dplyr)
library(data.table)
library(purrr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(ggpubr)
source("utils/eval.R")
source("utils/mini_import.R")

# Load combined outputs
invisible(list2env(readRDS("intmd/Combined_outputs_2025-02-21.rds"),environment()))
wgsCall_gr <- map(wgs_call,\(w) w[,c("seqnames","start","end","eventType")])

# Get available samples that have both epianeufinder and numbat results
epianeufinder_samples <- list.files("epianeufinder/cnv_epianeufinder", pattern = "epianeufinder_.*_cnv_calls.tsv") %>%
  gsub("epianeufinder_", "", .) %>%
  gsub("_cnv_calls.tsv", "", .)

congasp_samples <- list.files("congasp/cnv_congasp", pattern = ".*_cnv_segments_formatted.tsv") %>%
  gsub("_cnv_segments_formatted.tsv", "", .)

# Find samples with both epianeufinder and numbat ATAC results
common_samples_epi <- intersect(epianeufinder_samples, names(numbat_call))
common_samples_congasp <- intersect(congasp_samples, names(numbat_call))

print(paste("Samples with epianeufinder results:", paste(common_samples_epi, collapse = ", ")))
print(paste("Samples with congasp results:", paste(common_samples_congasp, collapse = ", ")))

# ---- epianeufinder vs ATAC_bin evaluation ----
epianeufinder_results <- map(common_samples_epi, function(sample) {
  # Load epianeufinder results
  epi_file <- file.path('epianeufinder/cnv_epianeufinder', paste0('epianeufinder_', sample, '_cnv_calls.tsv'))
  epi_dat <- fread(epi_file) %>% 
    select(seqnames, start, end, cnv) %>%
    mutate(seqnames = gsub('chr', '', seqnames))
  
  # Load numbat ATAC results
  if (!"ATAC_bincnt" %in% names(numbat_call[[sample]])) return(NULL)
  numbat_atac_dat <- numbat_call[[sample]][["ATAC_bincnt"]] %>%
    select(seqnames = CHROM, start = seg_start, end = seg_end, cnv = cnv_state_post) %>%
    mutate(seqnames = as.character(seqnames))
  
  # Evaluate epianeufinder
  pr_epi <- eval_call(wgsCall_gr[[sample]], epi_dat, byCNV = TRUE) %>%
    filter(cnv %in% c('amp','del'))
  pr_epi$method <- 'epianeufinder'; pr_epi$sample <- sample
  
  # Evaluate numbat ATAC
  pr_numbat_atac <- eval_call(wgsCall_gr[[sample]], numbat_atac_dat, byCNV = TRUE) %>%
    filter(cnv %in% c('amp','del'))
  pr_numbat_atac$method <- 'Numbat-multiome\nATAC Bin'; pr_numbat_atac$sample <- sample
  
  list(epianeufinder = pr_epi, numbat_atac = pr_numbat_atac)
}) %>% compact()

# ---- CONGASP vs comb_bin evaluation ----
congasp_results <- map(common_samples_congasp, function(sample) {
  # Load congasp results
  congasp_file <- file.path('congasp/cnv_congasp', paste0(sample, '_cnv_segments_formatted.tsv'))
  congasp_dat <- fread(congasp_file, stringsAsFactors = FALSE) %>%
    select(seqnames, start, end, cnv, cluster) %>%
    mutate(seqnames = gsub('chr', '', seqnames))
  
  # Load numbat combined results
  if (!"comb_bincnt" %in% names(numbat_call[[sample]])) return(NULL)
  numbat_comb_dat <- numbat_call[[sample]][["comb_bincnt"]] %>%
    select(seqnames = CHROM, start = seg_start, end = seg_end, cnv = cnv_state_post) %>%
    mutate(seqnames = as.character(seqnames))
  
  # Evaluate numbat combined
  pr_numbat_comb <- eval_call(wgsCall_gr[[sample]], numbat_comb_dat, byCNV = TRUE) %>%
    filter(cnv %in% c('amp','del'))
  pr_numbat_comb$method <- 'Numbat-multiome\nCombined Bin'; pr_numbat_comb$sample <- sample
  
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
  
  list(congasp = pr_congasp_max, numbat_comb = pr_numbat_comb)
}) %>% compact()

# Combine all results
epi_all_df <- map_dfr(epianeufinder_results, ~bind_rows(.x), .id = NULL)
congasp_all_df <- map_dfr(congasp_results, ~bind_rows(.x), .id = NULL)

# Combine into one dataframe for plotting
all_results <- bind_rows(epi_all_df, congasp_all_df)

# Reshape for plotting
all_results_long <- all_results %>%
  select(sample, method, cnv, precision, recall, f1) %>%
  pivot_longer(cols = c('precision', 'recall', 'f1'), 
               names_to = 'metric', values_to = 'value')

# Create paired boxplots
plot_metric_paired <- function(met, comparison_type) {
  if (comparison_type == "epianeufinder") {
    data_subset <- all_results_long %>% 
      filter(metric == met, 
             method %in% c('epianeufinder', 'Numbat-multiome\nATAC Bin'))
  } else {
    data_subset <- all_results_long %>% 
      filter(metric == met, 
             method %in% c('CONGAS+', 'Numbat-multiome\nCombined Bin'))
  }
  
  # Capitalize only first letter for y-axis label
  y_label <- paste0(toupper(substr(met, 1, 1)), tolower(substr(met, 2, nchar(met))))
  
  ggplot(data_subset, aes(x = method, y = value, group = sample)) +
    geom_boxplot(aes(group = method), width = 0.5, fill = 'gray90', color = 'black', alpha = 0.4, outlier.shape = NA) +
    geom_point(size = 2, color = 'darkred') +
    geom_line(color = 'darkred') +
    facet_wrap(~cnv) +
    labs(y = y_label, x = NULL) +
    theme_minimal() +
    theme(legend.position = 'none')
}

# Create plots for epianeufinder comparison (reordered: precision, recall, f1)
p_epi_precision <- plot_metric_paired('precision', 'epianeufinder')
p_epi_recall <- plot_metric_paired('recall', 'epianeufinder')
p_epi_f1 <- plot_metric_paired('f1', 'epianeufinder')

# Create plots for congasp comparison (reordered: precision, recall, f1)
p_congasp_precision <- plot_metric_paired('precision', 'congasp')
p_congasp_recall <- plot_metric_paired('recall', 'congasp')
p_congasp_f1 <- plot_metric_paired('f1', 'congasp')

# Combine plots using ggarrange and save as high-resolution PNG
print("epianeufinder vs Numbat-multiome ATAC Bin comparison:")
epi_combined <- ggarrange(p_epi_precision, p_epi_recall, p_epi_f1, ncol = 1, nrow = 3)
print(epi_combined)

print("CONGAS+ vs Numbat-multiome Combined Bin comparison:")
congasp_combined <- ggarrange(p_congasp_precision, p_congasp_recall, p_congasp_f1, ncol = 1, nrow = 3)
print(congasp_combined)

# Save plots as high-resolution PNG files with 6:4 ratio
ggsave("epianeufinder_vs_numbat_atac_comparison.png", 
       plot = epi_combined, 
       width = 6, height = 8, 
       dpi = 300, 
       units = "in")

ggsave("congasp_vs_numbat_combined_comparison.png", 
       plot = congasp_combined, 
       width = 6, height = 8, 
       dpi = 300, 
       units = "in")

# Summary statistics
summary_stats <- all_results %>%
  group_by(method, cnv) %>%
  summarise(
    mean_precision = mean(precision, na.rm = TRUE),
    mean_recall = mean(recall, na.rm = TRUE),
    mean_f1 = mean(f1, na.rm = TRUE),
    sd_precision = sd(precision, na.rm = TRUE),
    sd_recall = sd(recall, na.rm = TRUE),
    sd_f1 = sd(f1, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  )

print("Summary statistics:")
print(summary_stats)


