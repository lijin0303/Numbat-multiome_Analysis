##### Serial sample in MM #####
i <- "DL3267"
i <- "pM9916"
c(tumor_clone,cnvs)%<-%readRDS(glue("{i}_CNVsignals.rds"))
c(D,baf_segs,lfc_segs)%<-% CNV_plotD(tumor_clone,cnvs,zoomin,bamp_aware=T)
pop_CNV <- CNVp_base(D)+
  CNVp_exclude(zoomin)+
  CNVp_seg(lfc_segs,"logFC")+
  CNVp_seg(baf_segs,"pHF")+
  CNVp_theme()+
  theme(plot.background = element_blank(),
        plot.margin = margin(0.5,0.5,0.7,0.5, "cm"))
##### CLL: clone evolution #####
i <- "patA"
numbatRuns <- paste0(i,"_",modes)
x <- numbatRuns[4]
bulkL <- readRDS(glue("intmd/{x}_bulk_clones_final.rds"))
sample_segs <- readRDS(glue('intmd/{i}_cloneSegs.rds'))[[x]]
numbatClones <- readRDS(glue('intmd/{i}_cloneAssignment.rds'))[[x]]
CloneD <- numbatClones %>% 
  mutate(mode=ifelse(grepl("146p",cell),"ATAC","RNA")) %>% 
  group_by(clone,mode) %>% 
  summarise(ncell=n()) %>% 
  group_by(clone)%>%
  mutate(ntotal = sum(ncell),
         percent = round(ncell*100/ntotal)) %>% 
  as.data.frame() %>% 
  group_nest(clone,ntotal) %>% 
  mutate(clonelabel = paste0("clone ",as.numeric(clone)-1," (N = ",ntotal,")"))
pieD <- setNames(CloneD$data,CloneD$clonelabel)
ratios <- c(1,1,1,2.2)
walk(2:5,\(pop){
  zoomin <- c(1,5,7,9,10,11,12,16,17)
  c(D,baf_segs,lfc_segs)%<-% CNV_plotD(
    bulkL[[pop]],sample_segs[[pop]],zoomin)
  pop_CNV <- CNVp_base(D)+
    CNVp_exclude(zoomin)+
    CNVp_seg(lfc_segs,"logFC")+
    CNVp_seg(baf_segs,"pHF")+
    CNVp_theme()+
    theme(plot.background = element_blank(),
          plot.margin = margin(0.5,0.5,0.7,0.5, "cm"))
  
  clone <- pop -1
  circleClone <- ggClone(clonecols[pop],clone,c(30,6)*ratios[clone])
  p = ggarrange(circleClone,pop_CNV,widths = c(1,5))
  p <-  annotate_figure(p, top = text_grob(names(pieD)[pop], 
                                           color = "black", face = "bold", size = 20,
                                           vjust = 0.8))
  p <- p+ theme(plot.background = element_rect(colour = clonecols[pop], fill=NA, size=5))
  ggsave(p,
         filename = glue("Figures/patA_clone{clone}.pdf"),
         height=3.2,width=1.4*length(zoomin))
})
##### CLL:lost sensitivity #####
numbatClonesAll <- readRDS(glue('intmd/{i}_cloneAssignment.rds')) %>% bind_rows()
heatmapL <- map(3:1,\(m) mode_check(numbatClonesAll,c(m,4)))
heatcomb <- ggpubr::ggarrange(plotlist = heatmapL,nrow=1,legend = "top",common.legend = T)
ggsave(heatcomb,filename = "Figures/Clone_sensitivity_decreased.pdf",width=7,height=3)
##### CLL:clone composition #####
samplesOrder <- c("NIH-A-CLL-BM-01",
          "NIH-A-CLL-PB-01",
          "NIH-A-CLL-PB-02",
          "NIH-A-RT-BM-01"
)
x <- "patA_comb_bincnt"
min_cells <- 50
if(!file.exists("_data/cellAnnot_patA.rds")){
  cellAnnot_atac <- fread("_data/ATAC_cellannot.tsv") %>% 
    unite(cell,CB:pool,sep="#") %>% 
    select(cell,sampleID = `sample.ID`) %>% 
    mutate(source="ATAC")
  cellAnnot_rna <- fread("_data/RNA_cellannot_patA.csv",header = T) %>%
    select(cell=V1) %>% 
    mutate(sampleID=gsub(".*\\#","",cell),.after="cell") %>% 
    mutate(source="RNA")
  rbind(cellAnnot_atac,cellAnnot_rna) %>% saveRDS("_data/cellAnnot_patA.rds")
}else{
  cellAnnot <- readRDS("_data/cellAnnot_patA.rds")
}
i <- "patA"
x <- numbatRuns[4]
numbatClones <- readRDS(glue('intmd/{i}_cloneAssignment.rds'))[[x]]
numbatClones$clone <- paste0("Clone ",as.numeric(numbatClones$clone)-1)
numbatClones%<>% 
  inner_join(cellAnnot,by="cell") %>% 
  mutate(stage = case_when(sampleID=="NIH-A-CLL-BM-01"~"CLL",TRUE~"RT"))
fwrite(numbatClones[,1:2],"intmd/CNVclonecall_numbatmultiome.tsv",sep="\t")
names(clonecols) <- paste0("Clone ",(1:5)-1)  
stage_clone <- ggplot(numbatClones %>% 
                   group_by(stage,clone) %>% 
                   summarise(n=n()) %>% 
                    as.data.frame() %>% 
                   mutate(stage = as.numeric(factor(stage))) %>% 
                    filter(clone!="Clone 0"), 
                 aes(fill=clone, y=n, x=stage))+
  geom_area(position="fill", stat="identity")+
  scale_fill_manual(values= clonecols) +
  ggtitle(expression(Delta  * " in clone compositions during RT")) +
  scale_x_continuous(breaks=c(1,2),labels = c("Before","After"))+
  labs(y="",x="Richters' transformation")
clone_source <- ggplot(numbatClones%>% 
                    group_by(clone,source) %>% 
                    summarise(n=n()) %>% 
                    as.data.frame() %>% 
                      filter(clone!="Clone 0") %>% 
                    mutate(clone = as.numeric(factor(clone))), 
                  aes(fill=source, y=n, x=clone))+
  geom_area(position="fill", stat="identity")+
  scale_fill_manual(values= sourcecol) +
  ggtitle("Source of cells per clone") +
  labs(x="Inferred clone Number",y="")

p1 <- boxplot_theme(stage_clone,"top",gr=2)
p2 <- boxplot_theme(clone_source,"top")
clone_dynamics = ggarrange(p1,p2,nrow=1,widths = c(1.2,1))

ggsave("Figures/clone_composition.pdf",clone_dynamics,height=5,width=8)

