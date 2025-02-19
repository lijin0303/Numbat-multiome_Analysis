##### set up #####
setwd("~/numbat_EpiMultiome/Numbat-multiome_Analysis/")
source("utils/mini_import.R")
source("utils/vis.R")
pacman::p_load(karyoploteR,
               ggplotify,ggplot2,cowplot,
               ggpubr)
combinedout <- readRDS("intmd/Combined_outputs.rds")
pp <- getDefaultPlotParams(1)
pp$topmargin <- 500
meta_cohort <- readRDS("intmd/Analysis_meta.rds")
invisible(list2env(meta_cohort,environment()))
##### Karyotype to show CNVs #####
sample <- "pM10114"
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
chrs <- Reduce("intersect",map(segD,\(x) x$chr))
chrLevs <- c(paste0("chr",1:22))
relev_chrs <- chrLevs[chrLevs %in% chrs]
karyo_cnv <- as.ggplot(expression(
  kp <- plotKaryotype(plot.type=2, chromosomes=relev_chrs,cex=0.9,
                      plot.params = pp),
  kpDataBackground(kp, data.panel = 2, color = "#BBBBFF"),
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
    kpRect(kp, chr=CNV_D$chr, x0=0, x1=25*10e6, y0=0, y1=1)
    }))+
  theme(plot.margin=unit(c(-0.14,0.03,0,0), "null"))
ypos=0.783
karyo_cnv <- karyo_cnv+
  annotate(geom="text", 
           x=0.12, y=seq(ypos,ypos+0.046,length.out=5), 
           label=c("WGS",mode_label),
           size=rel(2.4),
           color=c("black",mode_cols),hjust = 0,fontface="bold") 
##### Sample-level Evals #####
samplerun_eval <- readRDS("intmd/samplerun_eval.rds")
samplerun_eval$mode <- factor(samplerun_eval$mode,levels=names(mode_label))
sample_eval <- samplerun_eval%>%
  gather(metric,value,-mode,-sample)%>%
  mutate(metric=stringr::str_to_title(metric)) %>% 
  ggVBS(mode,value,mode,legendpos="top",guider=4)+
  facet_wrap(~metric)+
  scale_fill_manual(labels=mode_label,values=mode_cols)+
  labs(title = "",x = "",y = "")+ 
  scale_x_discrete(labels = gsub(" ","\n ",mode_label))
##### CNV type stratified Evals #####
samplerun_cnv_eval <- readRDS("intmd/samplerun_CNV_eval.rds")
cnvCnt <- table(samplerun_cnv_eval$cnv)
cnv_color <- numbat:::cnv_colors[names(cnvCnt)]
cnv_map <- setNames(c("Amp","bAmp","Del","CNLoH"),names(cnv_color))
names(cnv_color)[1:4] <- as.character(cnv_map)

cnv_xLab <- paste0(names(cnv_color),"\n(N=",cnvCnt,")")
state_eval <- samplerun_cnv_eval%>%
  mutate(cnv = factor(cnv_map[cnv],levels=names(cnv_color))) %>%
  gather(metric,value,-mode,-sample,-cnv)%>%
  mutate(metric=stringr::str_to_title(metric)) %>% 
  ggVBS(cnv,value,cnv,legendpos="top",guider=4)+
  # ggVBS2(cnv,value,cnv,mode,legendpos="top",guider=4)+
  facet_wrap(~metric)+
  scale_fill_manual(values=cnv_color)+
  labs(title = "",x = "",y = "")+
  scale_x_discrete(labels = cnv_xLab)
##### upset plots #####
samplerun_event_eval <- readRDS("intmd/samplerun_event_eval.rds")
modeHit <- samplerun_event_eval %>% 
  filter(sample!="pM9916") %>%
  select(sample,cnvID,recall,mode) %>%
  unite("sample_cnv", c("sample","cnvID"),sep=":") %>% 
  filter(recall>0.8) %>% 
  mutate(mode = mode_label[mode]) %>% 
  group_nest(sample_cnv) %>% 
  mutate(mode = map(data,\(x) x$mode))
pacman::p_load(ggpubr,ggupset)
upset_detect <- ggplot(modeHit, aes(x=mode)) + 
  geom_bar(fill=c("gray30","#f2cc8f","#f2cc8f","gray30","gray30"),color="black") + 
  theme_pubr() + 
  scale_x_upset()+
  scale_y_continuous(expand = c(0, 0),limits = c(0,85))+
  theme(axis.ticks = element_blank())+
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  labs(x="",title="# of detected CNV events \n by mode combinations",y="")
##### Event length <> Recall #####
eventLen_scatter <- samplerun_event_eval %>% 
  filter(sample!="pM9916") %>%
  mutate(mode_label = mode_label[mode]) %>% 
  mutate(cnvID = gsub("_"," ",cnvID)) %>% 
  unite("sample_cnv", c("sample","cnvID_label"),sep=":",remove=F) %>%
  mutate(eventLen = eventLen/1000) %>% 
  mutate(mode = factor(mode,levels=modes)) %>% 
  ggscatter(x = "eventLen", y = "recall",
            xscale="log10",
            shape="mode",
            color = "mode",
            size = rel(4))+
  #geom_vline(xintercept = 200,linetype="dashed")+
  scale_shape_manual(labels=mode_label,values=mode_shapes)+
  scale_color_manual(labels=mode_label,values=mode_cols)+
  theme_bw()+
  ggrepel::geom_label_repel(
    show_guide  = F,
    alpha=0.7,
    aes(label = if_else(recall<0.6, gsub(":","\n",sample_cnv), ""),
        fill=mode),
    max.overlaps = Inf,
    alpha = 0.6,
    seed = 1234
  ) +
    ggrepel::geom_label_repel(
     show_guide  = F,
     aes(label = if_else(recall<0.6, 
                         gsub(":","\n",sample_cnv), "")),
    max.overlaps = Inf,
     alpha = 1,
     fill = NA,
     seed = 1234)+
  scale_fill_manual(labels=mode_label,values=mode_cols)+
  rremove("legend.title")+
  theme(
    axis.line = element_line(colour = "black",size=1),
    axis.ticks = element_line(size=1,color="black"),
    axis.text = element_text(color="black"),
    axis.ticks.length=unit(0.2,"cm"),
    legend.position = "top")+
  font("xylab",size=rel(1.3))+  
  font("xy",size=rel(1.3))+ 
  font("xy.text", size = rel(1.5)) +  
  font("legend.text",size = rel(1.2))+
  guides(color = guide_legend(
    override.aes = list(alpha = 1,size=rel(4)),ncol = 2))+
  labs(x="CNV event length in kb (axis in log10 scale)")

list("sample_eval"=sample_eval,
     "state_eval"=state_eval,
     "karyotype"=karyo_cnv,
     "upset_plot"=upset_detect,
     "eventLen_scatter"=eventLen_scatter) %>% 
  saveRDS("Figures/figure_obj.rds")
##### Put everything together #####
vbs2 <- plot_grid(sample_eval, state_eval,ncol=1, labels=LETTERS[1:2])
combined_eval <- plot_grid(vbs2, karyo_cnv,ncol=2, labels=c("","C"),
                           rel_widths=c(1.8,1))
event_eval <- plot_grid(upset_detect, eventLen_scatter,
                        ncol=2, labels=c("D","E"),
                        rel_widths=c(0.7,1))
extrafont::loadfonts()
p <- plot_grid(combined_eval, event_eval,ncol=1,rel_heights =c(2.1,1))
ggsave(p,filename = "Figures/Figure2.pdf",width=12,height=13,units = "in")

