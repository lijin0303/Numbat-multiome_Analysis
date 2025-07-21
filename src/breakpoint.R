source("utils/mini_import.R")
source("utils/vis.R")
pacman::p_load(karyoploteR,
               ggplotify,ggplot2,cowplot,
               ggpubr)
wgs_df <- fread("intmd/patA_wgs_seg.tsv")
dirs <- fs::dir_ls(path = "binsize", 
                   regexp = "patA_bin\\d+kb_outputs$", 
                   type = "directory")
numbat_segL <- dirs %>% 
  set_names(~stringr::str_extract(.x, "bin\\d+kb")) %>%
  map(~{
    file_path <- file.path(.x, "segs_consensus_2.tsv")
    if (file.exists(file_path)) {
      readr::read_tsv(file_path, show_col_types = FALSE) %>%
        select(chr = CHROM, start = seg_start, end = seg_end, 
               eventType = cnv_state_post) %>%
        mutate(chr = paste0("chr", chr))
    } else {
      NULL
    }
  }) %>% 
  compact()
ords <- paste0("bin",c(50,80,100,200,300,500),"kb")
segD <- c(list("WGS"=wgs_df),numbat_segL[ords]) 
chrLen <- numbat::chrom_sizes_hg38$size %>% 
  setNames(paste0("chr",numbat::chrom_sizes_hg38$CHROM))
pp <- getDefaultPlotParams(1)
pp$data2height <- 50
karyo_cnv <- as.ggplot(expression(
  kp <- plotKaryotype(plot.type=2, chromosomes=c("chr7"),
                      cex=0.9,
                      plot.params = pp),
  nparts <- names(segD),
  for(i in seq_along(nparts)) {
    CNV_D <- segD[[i]]
    at <- autotrack(i,length(nparts),margin = 0.05)
    kpRect(kp, chr=CNV_D$chr,
           x0=CNV_D$start, 
           x1=CNV_D$end, 
           y0=0.1, y1=0.85, 
           border="black", 
           lty=1, lwd=0.5, 
           r0=at$r0, r1=at$r1,
           col=numbat:::cnv_colors[c(CNV_D$eventType)])
    kpRect(kp, chr=CNV_D$chr, x0=0, x1=chrLen[CNV_D$chr], y0=0, y1=1)
  }))+
  theme(plot.margin=unit(c(0,0.03,0,0), "null"))
ypos=0.42
karyo_cnv <- karyo_cnv+
  annotate(geom="text", 
           x=0.04, y=seq(ypos,ypos+0.31,length.out=length(segD)), 
           label=c("WGS",ords),
           size=rel(2.3),
           color="black",hjust = 0,fontface="bold") 
karyo_cnv

# Save karyotype plot as high-resolution PNG
ggsave("_panel/karyotype_chr7_binsize_comparison.png", 
       plot = karyo_cnv, 
       width = 12, height = 8, 
       dpi = 300, 
       units = "in")

cat("Karyotype plot saved as: karyotype_chr7_binsize_comparison.png\n")

# ============================================================================
# Bin Size vs Runtime Analysis
# ============================================================================

# Load runtime data
runtime_data <- fread("binsize/runtime.tsv") %>%
  mutate(bin_size_kb = as.numeric(size_kb),
         runtime_minutes = as.numeric(minutes))

# Create scatterplot with connected lines
runtime_plot <- ggplot(runtime_data, aes(x = bin_size_kb, y = runtime_minutes)) +
  geom_point(size = 4, color = "darkred", alpha = 0.8) +
  geom_line(color = "darkred", linewidth = 1.2, alpha = 0.7) +
  geom_text(aes(label = paste0(bin_size_kb, "kb")), 
            vjust = -0.8, hjust = 0.5, size = 3, fontface = "bold") +
  scale_x_log10(breaks = c(50, 80, 100, 200, 300, 500),
                labels = c("50kb", "80kb", "100kb", "200kb", "300kb", "500kb")) +
  scale_y_continuous(limits = c(0, max(runtime_data$runtime_minutes) * 1.1)) +
  labs(title = "Bin Size vs Runtime Trade-off",
       subtitle = "Computational cost decreases with larger bin sizes",
       x = "Bin Size (log scale)",
       y = "Runtime (minutes)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12, color = "gray40"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        panel.grid.minor = element_blank())

# Calculate and display key statistics
cat("=== Bin Size vs Runtime Analysis ===\n")
cat("Runtime scaling factor (500kb vs 50kb):", 
    round(runtime_data$runtime_minutes[runtime_data$bin_size_kb == 500] / 
          runtime_data$runtime_minutes[runtime_data$bin_size_kb == 50], 2), "x faster\n")
cat("Recommended stable range: 100-200kb (runtime: 11-15 minutes)\n")
cat("High-resolution option: 50kb (runtime: 21 minutes)\n")
cat("Fast screening option: 300-500kb (runtime: 6-9 minutes)\n\n")

# Display the runtime plot
print(runtime_plot)