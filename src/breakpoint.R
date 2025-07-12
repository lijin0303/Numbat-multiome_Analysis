sample <- "patA"
wgs_seg <- combinedout$wgs_call[[sample]] %>% 
  select(chr=seqnames,start,end,eventType) %>% 
  mutate(chr=paste0("chr",chr))
numbat_seg <- map(combinedout$numbat_call[[sample]],\(d) d %>% 
                    select(chr=CHROM,start=seg_start,
                           end=seg_end,
                           eventType=cnv_state_post) %>% 
                    mutate(chr=paste0("chr",chr)))
names(numbat_seg) <- mode_label
segD <- c(list("WGS"=wgs_seg),numbat_seg) 
chrLen <- numbat::chrom_sizes_hg38$size %>% 
  setNames(paste0("chr",numbat::chrom_sizes_hg38$CHROM))
pp$data2height <- 50
karyo_cnv <- as.ggplot(expression(
  kp <- plotKaryotype(plot.type=2, chromosomes=c("chr7"),cex=0.9,
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
ypos=0.901
karyo_cnv <- karyo_cnv+
  annotate(geom="text", 
           x=0.12, y=seq(ypos,ypos+0.046,length.out=5), 
           label=c("WGS",mode_label),
           size=rel(2.3),
           color=c("black",mode_cols),hjust = 0,fontface="bold") 