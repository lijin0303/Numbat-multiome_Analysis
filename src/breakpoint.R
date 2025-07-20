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