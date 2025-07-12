setwd("~/numbatm/Numbat-multiome_Analysis/")
source("utils/mini_import.R")
source("utils/vis.R")
source("utils/eval.R")
i <- "DL3267"
c(tumor_clone,cnvs)%<-%readRDS(glue("intmd/{i}_CNVsignals.rds"))
zoomin <- 1:22
c(D,baf_segs,lfc_segs)%<-% CNV_plotD(tumor_clone,cnvs,zoomin,bamp_aware=T)
pop_CNV1 <- CNVp_base(D)+
  CNVp_exclude(zoomin)+
  CNVp_seg(lfc_segs,"logFC")+
  CNVp_seg(baf_segs,"pHF")+
  CNVp_theme()+
  theme(plot.background = element_blank(),
        plot.margin = margin(0.5,0.5,0,0.5, "cm"),
        axis.title.y = element_text(angle = 0,vjust = 0.5,face="bold",size=rel(1.8),
                                    margin = margin(0,20,0,0)))+
  labs(x="    ",y="Time point 1")

i <- "pM9916"
c(tumor_clone,cnvs)%<-%readRDS(glue("intmd/{i}_CNVsignals.rds"))
c(D,baf_segs,lfc_segs)%<-% CNV_plotD(tumor_clone,cnvs,zoomin,bamp_aware=T)
pop_CNV2 <- CNVp_base(D)+
  CNVp_exclude(zoomin)+
  CNVp_seg(lfc_segs,"logFC")+
  CNVp_seg(baf_segs,"pHF")+
  CNVp_theme()+
  theme(plot.background = element_blank(),
        plot.margin = margin(0,0.5,0,0.5, "cm"),
        strip.text.x = element_blank(),
        axis.title.y = element_text(angle = 0,hjust = 1,vjust = 0.5,face="bold",size=rel(1.8),
                                    margin = margin(0,20,0,0)))+
  labs(y="Time point 2 \n (+1yr)")

library(ggpubr)
comb_cnv <- ggpubr::ggarrange(pop_CNV1,pop_CNV2,ncol=1,
                             title="MM3",align = "v",heights = c(1.3,1))+
  theme(plot.margin = margin(0,0.5,-6,0.5, "cm"))
ggsave("Figures/fig3A.pdf",comb_cnv,width=11.5,height=4.5)