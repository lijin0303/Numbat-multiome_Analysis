##### set up #####
setwd("~/numbatm/Numbat-multiome_Analysis/")
source("utils/mini_import.R")
source("utils/vis.R")
pacman::p_load(karyoploteR,
               ggplotify,ggplot2,cowplot,
               ggpubr)
meta_cohort <- readRDS("intmd/Analysis_meta.rds")
invisible(list2env(meta_cohort,environment()))
df2gr <- function(df){
  return(df %>% makeGRangesFromDataFrame(keep.extra.columns = T))
}
zoomin <- c(1,5,7,9:12,16,17)
i <- "patA"
numbatRuns <- paste0(i,"_",modes)
x <- numbatRuns[4]
bulkL <- readRDS(glue("intmd/{x}_bulk_clones_final.rds"))
sample_segs <- readRDS(glue('intmd/{i}_cloneSegs.rds'))[[x]]
numbatClones <- readRDS(glue('intmd/{i}_cloneAssignment.rds'))[[x]] %>% 
  mutate(source=ifelse(grepl("146p",cell),"ATAC","RNA")) %>% 
  mutate(pop = clone) %>% 
  filter(pop != "1")%>% 
  mutate(clone = as.numeric(factor(clone)))
  
##### CLL: clone evolution #####
CloneD <- numbatClones %>% 
  group_by(pop,source) %>% 
  summarise(ncell=n()) %>% 
  spread(source,ncell) %>% 
  mutate(total=RNA+ATAC) %>% 
  group_split(pop)
ratios <- c(1,1,1,1.7)*0.7
library(ggtext) 
clonep <- function(pop){
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
  pop_CNV <- pop_CNV +
    labs(
      title = glue("**Clone {clone} 
    <span style='font-size:11pt'>(N =
    <span style='color:{sourcecol['ATAC']};'>{CloneD[[clone]]$ATAC}</span> + 
    <span style='color:{sourcecol['RNA']};'>{CloneD[[clone]]$RNA}</span> = 
    <span style='color:black;'>{CloneD[[clone]]$total})</span>
    </span>**")
    ) +
    theme(
      plot.title = element_markdown(lineheight = 1.1,hjust = 0.3),
      legend.text = element_markdown(size = 11)
    )
  
  circleClone <- ggClone(clonecols[pop],clone,c(30,6)*ratios[clone])
  p = ggarrange(circleClone,pop_CNV,widths = c(1,5))
  return(p)
}
cloneL <- map(2:5,clonep)
ggarrange(plotlist = cloneL,ncol=1) %>% 
  ggsave(filename="Figures/fig3B.pdf",width=7,height=10)
##### CLL: clone-modality #####
rnap <- table(numbatClones$source)[2]/nrow(numbatClones)
errorbarD <- numbatClones%>% 
  group_nest(clone)  %>%
  mutate(
    fit = map(data,\(x) binom.test(sum(x$source=="RNA"), nrow(x), p = rnap)),
    tidied = map(fit, broom::glance)
  ) %>% 
  tidyr::unnest(tidied) %>% 
  select(clone,cil=conf.low,cih=conf.high)

clone_source <- ggplot(numbatClones%>% 
                         group_by(clone,source) %>% 
                         summarise(n=n()) %>% 
                         as.data.frame(), 
                       aes(fill=source, y=n, x=clone))+
  geom_bar(position="fill", stat="identity")+
  # geom_area(position="fill", stat="identity")+
  scale_fill_manual(values= sourcecol) +
  ggtitle("Source of cells per clone") +
  labs(x="Inferred clone Number",y="")+
  geom_errorbar(inherit.aes = FALSE,
                aes(x=clone,ymin=cil, ymax=cih), 
                data=errorbarD,width=0.2)
clone_source <- boxplot_theme(clone_source,"top",2)
ggsave("Figures/fig3D.pdf",clone_source,width = 5,height = 5)
##### CLL:lost sensitivity (grayscale) #####
numbatClonesAll <- readRDS(glue('intmd/{i}_cloneAssignment.rds')) %>% bind_rows()
numbatClonesAll%<>%mutate(
  clone = case_when(mode=="patA_comb_bincnt" & clone=="1"~"N",
                    mode=="patA_comb_bincnt" & clone!="1"~paste0("Cl",as.numeric(clone)-1),
                    TRUE~clone)
)
heatmapL <- map(3:1,\(m) mode_check(numbatClonesAll,c(m,4)))
heatcomb <- ggpubr::ggarrange(plotlist = heatmapL,ncol=1,legend = "top",common.legend = T)
ggsave("Figures/fig3C.pdf",heatcomb,width = 4,height = 10)
##### CLL:clone composition #####
samplesOrder <- c("NIH-A-CLL-BM-01",
          "NIH-A-CLL-PB-01",
          "NIH-A-CLL-PB-02",
          "NIH-A-RT-BM-01"
)
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
numbatClones$cloneL <- paste0("Clone ",numbatClones$clone)
numbatClones%<>% 
  inner_join(cellAnnot,by="cell") %>% 
  mutate(stage = case_when(sampleID=="NIH-A-CLL-BM-01"~"CLL",TRUE~"RT"))
fwrite(numbatClones[,1:2],"intmd/CNVclonecall_numbatmultiome.tsv",sep="\t")
names(clonecols) <- paste0("Clone ",(1:5)-1)  
stage_clone <- ggplot(numbatClones %>% 
                   group_by(stage,cloneL) %>% 
                   summarise(n=n()) %>% 
                   as.data.frame() %>% 
                   mutate(stage = as.numeric(factor(stage))), 
                 aes(fill=cloneL, y=n, x=stage))+
  geom_area(position="fill", stat="identity")+
  scale_fill_manual(values= clonecols) +
  ggtitle(expression(Delta  * " in clone compositions")) +
  scale_x_continuous(breaks=c(1,2),labels = c("CLL","RS"))+
  labs(y="",x="transformation")

ggsave("Figures/fig3E.pdf",
       boxplot_theme(stage_clone,"right",gr=1),
       height=6,width=5.5)



