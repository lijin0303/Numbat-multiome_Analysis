##### set up #####
# setwd("~/numbatm/Numbat-multiome_Analysis/")
source("utils/mini_import.R")
source("utils/vis.R")
pacman::p_load(karyoploteR,
               ggplotify,ggplot2,cowplot,
               ggpubr)
combinedout <- readRDS("intmd/Combined_outputs_2025-02-21.rds")
pp <- getDefaultPlotParams(1)
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
chrLen <- numbat::chrom_sizes_hg38$size %>% 
  setNames(paste0("chr",numbat::chrom_sizes_hg38$CHROM))
pp$data2height <- 50
karyo_cnv <- as.ggplot(expression(
  kp <- plotKaryotype(plot.type=2, chromosomes=relev_chrs,cex=0.9,
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
##### Sample-level Evals #####
samplerun_eval <- readRDS("intmd/samplerun_eval_2025-02-21.rds")
samplerun_eval$mode <- factor(samplerun_eval$mode,levels=names(mode_label))
library(lme4)
lmer_anova <- map_dbl(c("f1","precision","recall"),\(x){
  anova_lmer <- car::Anova(lmer(as.formula(glue("{x}~mode+ (1 | sample)")),samplerun_eval))
  return(anova_lmer$`Pr(>Chisq)`)
})
sigTest <- data.frame(metric=stringr::str_to_title(c("f1","precision","recall")),
                      pval = paste0("p=",round(lmer_anova,2)))
sample_eval <- samplerun_eval%>%
  gather(metric,value,-mode,-sample)%>%
  mutate(metric=stringr::str_to_title(metric)) %>% 
  ggVBS(mode,value,mode,legendpos="top",guider=4)+
  facet_wrap(~metric)+
  scale_fill_manual(labels=mode_label,values=mode_cols)+
  labs(title = "",x = "",y = "score")+ 
  scale_x_discrete(labels = gsub(" ","\n ",mode_label))+ 
  geom_text(
    inherit.aes = F,
    data    = sigTest,
    mapping = aes(x = 1.4, y =0.88, label = pval))
summary_df <- samplerun_eval %>% 
  group_by(mode) %>% 
  summarise(across(where(is.numeric),
                list(mean = ~mean(.),
                     median = ~median(.),
                     q25 = ~quantile(., 0.25),
                     q75 = ~quantile(., 0.75)),
                .names = "{.col}_{.fn}"))
##### CNV type stratified Evals #####
samplerun_cnv_eval <- readRDS("intmd/samplerun_CNV_eval_2025-02-21.rds")
cnvCnt <- table(samplerun_cnv_eval$cnv)
cnv_color <- numbat:::cnv_colors[names(cnvCnt)]
cnv_map <- setNames(c("Amp","bAmp","Del","CNLoH"),names(cnv_color))
names(cnv_color)[1:4] <- as.character(cnv_map)
cnv_xLab <- paste0(names(cnv_color),"\n(N=",cnvCnt,")")
samplerun_cnv_eval%<>%
  mutate(cnv = factor(cnv_map[cnv],levels=names(cnv_color)))
library(lme4)
lmer_anova <- map_dbl(c("f1","precision","recall"),\(x){
  anova_lmer <- car::Anova(lmer(as.formula(glue("{x}~cnv+ (1 | mode)")),samplerun_cnv_eval))
  return(anova_lmer$`Pr(>Chisq)`)
})
sigTest <- data.frame(metric=stringr::str_to_title(c("f1","precision","recall")),
                      pval = paste0("p=",round(lmer_anova,4)))
my_comparisons <- list( c("Amp", "bAmp"), c("Amp", "Del"), c("Amp", "CNLoH") )
state_eval <- samplerun_cnv_eval %>%
  gather(metric,value,-mode,-sample,-cnv)%>%
  mutate(metric=stringr::str_to_title(metric)) %>% 
  ggVBS(cnv,value,cnv,legendpos="top",guider=4)+
  stat_compare_means(comparisons=my_comparisons,label.y = seq(1,1.07,length.out=3),
                     method="wilcox.test")+
  facet_wrap(~metric)+
  scale_fill_manual(values=cnv_color)+
  labs(title = "",x = "",y = "score")+
  scale_x_discrete(labels = cnv_xLab)+ 
  scale_y_continuous(limits=c(0.7,1.12),breaks = c(0.8,0.9,1))+
  geom_text(
    inherit.aes = F,
    data    = sigTest,
    mapping = aes(x = 1.4, y =0.75, label = pval))
##### upset plots #####
samplerun_event_eval <- readRDS("intmd/samplerun_event_eval_2025-02-21.rds")
modeHit <- samplerun_event_eval %>% 
  filter(sample!="pM9916") %>%
  select(sample,cnvID,recall,mode) %>%
  unite("sample_cnv", c("sample","cnvID"),sep=":") %>% 
  filter(recall>0.8) %>% 
  mutate(mode = mode_label[mode]) %>% 
  group_nest(sample_cnv) %>% 
  mutate(mode = map(data,\(x) x$mode))
pacman::p_load(ggpubr,ggupset)
levCol <- c("gray30","#bc6c25","#778da9","#a3b18a")
upset_detect <- ggplot(modeHit, aes(x=mode)) + 
  geom_bar(fill=levCol,color="black") + 
  theme_pubr() + 
  scale_x_upset()+
  scale_y_log10()+
  scale_y_continuous(expand = c(0, 0),limits = c(0,90))+
  theme(axis.ticks = element_blank())+
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  labs(x="",title="# of detected CNV events (recall > 0.8)",y="")+
  theme(plot.title = element_text(hjust = 0.5),
    legend.position = "bottom")

upset_cats <- modeHit %>% 
  mutate(detect_cat =  map_chr(mode,\(x) paste0(sort(x),collapse = ";"))) %>% 
  distinct(sample_cnv,detect_cat)
upset_levs <- table(upset_cats$detect_cat) %>% sort(decreasing = T)
levCol <- setNames(levCol,names(upset_levs))
upset_cats$detect_cat <- factor(upset_cats$detect_cat,levels=names(upset_levs))
##### Event length <> Recall #####
ylpos = c(0.1,0.66,0.47,0.53)
xlpos = c(rep(10000,2),10^5*0.24,10^5*0.43)
ids <- qs::qread("_data/maskids.qs")
ids <- c(ids,"patA"="patA")
eventLen_scatterD <- samplerun_event_eval %>% 
  filter(sample!="pM9916") %>% 
  mutate(mode_label = mode_label[mode]) %>% 
  unite("sample_cnv", c("sample","cnvID"),sep=":",remove=F) %>% 
  inner_join(upset_cats,by="sample_cnv") %>% 
  mutate(sampleM = ids[sample]) %>% 
  unite("sample_cnv", c("sampleM","cnvID_label"),sep=":",remove=F) %>%
  mutate(eventLen = eventLen/1000) %>% 
  mutate(mode = factor(mode,levels=modes)) 
labelD <- eventLen_scatterD%>% 
  filter(recall<0.6) %>% 
  mutate(label = gsub(":","\n",sample_cnv)) %>% 
  arrange(eventLen,recall) %>%
  filter(mode=="comb_bincnt") %>% # it turns out those event fail for both mode very similarly
  mutate(ypos =ylpos) %>% 
  mutate(xpos = xlpos)
SegD <- labelD %>%
  filter(cnvID=="21q_arm_amp") %>% 
  select(label,
         x_end = xpos,
         x_start = eventLen,
         y_end = ypos,
         y_start = recall) %>% 
  mutate(x_end=x_end*1.2,y_end = y_end+0.06)
mode_shapes[1:4] <- c(0,1,2,4)
eventLen_scatterD$revSize = 1/(eventLen_scatterD$recall+0.54)
eventLen_scatter <-   ggscatter(eventLen_scatterD,
            x = "eventLen", y = "recall",
            shape="mode",
            color = "detect_cat",
            size =  "revSize")+
  geom_hline(yintercept = 0.8,color="gray50")+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  scale_x_log10(
    breaks=c(10000,30000,100000),
    labels=paste0(c(10000,30000,100000),"kb"))+
  scale_size_continuous(
    range = c(3,5),
    guide="none"
  )+
  scale_shape_manual(labels=mode_label,
                     values=mode_shapes)+
  scale_color_manual(values=levCol,guide="none")+
  theme_bw()+
  geom_segment(inherit.aes = FALSE,
               data = SegD,
               aes(x = x_start, xend = x_end, 
                   y = y_start, yend = y_end),
               arrow = arrow(ends = "last",length = unit(0.01, "npc")),
               color = 'black',
               size = 0.6)+
  geom_text(
    data = labelD,
    aes(x = xpos,y=ypos,label = label),
    inherit.aes = FALSE,
    show_guide  = F,
    color="#bc6c25",face="bold",
    hjust = 0,
  )+
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
  guides(shape = guide_legend(
    override.aes = list(alpha = 1,size=rel(4)),ncol = 2))+
  labs(x="CNV event length in kb (axis in log10 scale)")
eventLen_scatter
scatterD <- eventLen_scatterD %>% 
  select(sample_cnv,sample,cnvID_label,
         recall,mode,detect_cat,revSize)
saveRDS(scatterD,"intmd/scatterD_cnv.rds")
##### Put everything together #####
vbs2 <- plot_grid(sample_eval, state_eval,ncol=1, labels=LETTERS[2:3])
combined_eval <- plot_grid(karyo_cnv,vbs2,ncol=2, labels=c("A",""),
                           rel_widths=c(1,1.8))
event_eval <- plot_grid(upset_detect, eventLen_scatter,
                        ncol=2, labels=c("D","E"),
                        rel_widths=c(0.75,1))
extrafont::loadfonts()
p <- plot_grid(combined_eval, event_eval,ncol=1,rel_heights =c(1.9,1.2))
ggsave(p,filename = "Figures/fig2.pdf",width=12,height=13.5,units = "in")

