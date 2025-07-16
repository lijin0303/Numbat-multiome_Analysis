library(dplyr)
source("utils/eval.R")

# Load WGS (patA) calls
wgs <- read.delim("pat_comb/segs_consensus_2.tsv", stringsAsFactors = FALSE)
wgs <- wgs %>%
  filter(!is.na(cnv_state), cnv_state != "neu") %>%
  transmute(seqnames = CHROM, start = seg_start, 
            end = seg_end, cnv = cnv_state) %>% 
  filter(seqnames==7)

# Load CONGASp calls
congasp <- read.delim("congasp/cnv_segments_formatted.tsv", stringsAsFactors = FALSE)

# Get unique clusters
clusters <- sort(unique(congasp$cluster))

# Loop through clusters and evaluate
for (cl in clusters) {
  cat("\nCluster", cl, ":\n")
  congasp_cl <- congasp %>% filter(cluster == cl) %>% 
    select(seqnames, start, end, cnv) %>% 
    unique() %>% 
    mutate(seqnames=gsub("chr","",seqnames))
  print(eval_call(wgs, congasp_cl, byCNV=TRUE))
}

numbat_seg <- map(combinedout$numbat_call[["patA"]],\(d) d %>% 
                    select(chr=CHROM,start=seg_start,
                           end=seg_end,
                           eventType=cnv_state_post) %>% 
                    filter(chr==7))

map(numbat_seg,\(x) eval_call(wgs, x, byCNV=TRUE))


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
