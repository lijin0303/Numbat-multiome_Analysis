##### CNV event overview ####
setwd("~/numbatm/numbat-multiome_Analysis/")
source("utils/mini_import.R")
source("utils/vis.R")
samplerun_event_eval <- readRDS("intmd/2025-02-01/samplerun_event_eval.rds")
levs <- paste0(rep(c(1:22),each=2),rep(c("p","q"),22))
levs <- levs[levs %in% samplerun_event_eval$arm]

labelOrd <- samplerun_event_eval %>% 
  filter(sample!="pM9916") %>%
  distinct(arm,annot3,eventType,cnvID_label) %>% 
  mutate(arm = factor(arm,levels=levs)) %>% 
  mutate(eventType = factor(eventType)) %>% 
  arrange(arm,annot3,eventType) %>% 
  pull(cnvID_label)
labels <- samplerun_event_eval %>% 
  filter(sample!="pM9916") %>%
  distinct(sample,cnvID,cnvID_label) %>% 
  mutate(cnvID2 = gsub("bamp","amp",cnvID))
CNV_binary <- map(list.files("intmd","_CNV_coverage",full.names = T),readRDS) %>% 
  bind_rows() %>% 
  group_by(event,sample) %>% 
  summarise(avg_recall=mean(CNV_recall)) %>%
  filter(sample!="pM9916") %>% 
  mutate(CNV=1) %>%
  mutate(sample2 = gsub("-.*","",sample)) %>% 
  mutate(event2=gsub("bamp","amp",event)) %>% 
  inner_join(labels,by=c("sample2"="sample","event2"="cnvID2")) %>% 
  filter(!(cnvID_label=="11q23.1-q21(del)" & sample=="patA-RT")) %>% 
  as.data.frame() %>% 
  select(sample,cnvID_label,CNV) %>% 
  spread(cnvID_label,CNV) %>% 
  column_to_rownames("sample") %>% 
  replace(is.na(.), 0) %>% 
  t() 
colors <- c("white", "black")  # Blue, White, Red
colord <- c("pM10975","pM10109","pM10114",
            "pM11004" ,"pM1835",
            "DL3267",
            "patA-CLL","patA-RS")
pdf("Figures/supp1.pdf",width=5.5,height=7)
p <- pheatmap::pheatmap(CNV_binary[labelOrd,colord],
                   color = c("white", "black"),
                   cluster_rows = F,cluster_cols = F,
                   border_color = "grey60",
                   fontsize_col = rel(8),
                   fontsize_row = rel(9),
                   angle_col = 0,legend=F)
print(p)
dev.off()


CNV_binary <- readRDS("intmd/cohort_CNVevents.rds")
ids <- qs::qread("_data/maskids.qs")
colnames(CNV_binary)[1:6] <- ids[colnames(CNV_binary)[1:6]]
colnames(CNV_binary)[8] <- "patA-RS" 
pdf("Figures/suppFig1.pdf",width=5.5,height=7)
p <- pheatmap::pheatmap(CNV_binary,
                        color = c("white", "black"),
                        cluster_rows = F,cluster_cols = F,
                        border_color = "grey60",
                        fontsize_col = rel(8),
                        fontsize_row = rel(9),
                        angle_col = 0,legend=F)
print(p)
dev.off()


##### arm signal info (not pursued) #####
# gene_cor <- readRDS("_data/gtf_gcnt.rds")%>% 
#   mutate(genecor = paste0("chr",CHROM,":",gene_start,"-",gene_end))
# rownames(ref_comb) <- setNames(gene_cor$genecor,gene_cor$gene)[rownames(ref_comb)]
setwd("~/numbatm/numbat-multiome_Analysis/")
source("utils/mini_import.R")
armF <- "_data/chrom_arm.rds"
armdf <- readRDS(armF) %>%
  as.data.frame() %>%
  mutate(seqnames=paste0("chr",seqnames)) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)
sum_meta <- function(grA, grB) {
  # Find overlaps
  hits <- findOverlaps(grB, grA)
  
  # Extract metadata columns from grA
  meta_cols <- colnames(mcols(grA))
  
  # Initialize output matrix
  result <- matrix(0, nrow = length(grB), ncol = length(meta_cols))
  colnames(result) <- meta_cols
  
  # For each metadata column, sum values from grA into grB
  for (col in meta_cols) {
    vals <- mcols(grA)[[col]]
    if (!is.numeric(vals)) next  # Only process numeric columns
    
    sums <- tapply(vals[subjectHits(hits)],
                   INDEX = queryHits(hits),
                   FUN = sum)
    
    result[as.integer(names(sums)), col] <- sums
  }
  
  # Add the results to grB metadata
  mcols(grB) <- cbind(mcols(grB), as.data.frame(result))
  return(grB)
}

arm_info <- sum_meta(ref_gr,armdf)
arm_signalD <- as.data.frame(arm_info)
arm_signalD$chr_arm <- factor(arm_signalD$chr_arm,levels=arm_info$chr_arm)
##### CNV event baseline signal calculation #####
ref_comb <- readRDS("_data/lambdas_comb_bincnt.rds") %>% as.matrix()
sum_overlap <- function(grA, grB){
  hits <- findOverlaps(grB, grA)
  return(colSums(as.data.frame(mcols(grA)[subjectHits(hits),])))
}
binmat2gr <- function(x){
  as.data.frame(x)%>%
    rownames_to_column("id") %>% 
    separate(id,into = c("chr","cors"),sep=":",remove=F)%>% 
    separate(cors,into = c("start","end"),sep="-") %>% 
    select(-id) %>% 
    mutate(strand = "*") %>% 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
}
ref_gr <- binmat2gr(ref_comb)
eventSegD <- readRDS("intmd/eventSegD.rds") %>% 
  mutate(sample=gsub("-.*","",sample))%>%
  filter(!(sample=="pM10975" & event=="15q_arm_amp")) %>% 
  mutate(event = gsub("15q_arm_amp","15q_arm_bamp",event))
samplerun_event_eval <- readRDS("intmd/2025-02-01/samplerun_event_eval.rds") 
eventSegD%<>%
  left_join(samplerun_event_eval %>% filter(sample!="pM9916") %>% 
               select(event=cnvID,sample,cnvID_label) %>% 
               distinct()) %>% 
  select(-event)
eventSegD$grs <- map(eventSegD$data,\(x) x%>% 
                       mutate(seqnames=paste0("chr",seqnames)) %>% 
                       makeGRangesFromDataFrame())
eventSignalD <- eventSegD %>% 
  select(sample,cnvID_label) %>% 
  cbind(map(eventSegD$grs,\(d) sum_overlap(ref_gr,d)) %>% bind_rows()) %>% 
  unite("combID",sample:cnvID_label,sep=":",remove = F)
refsig <- eventSignalD%>%
  mutate(ATAC_refSig = rowMeans(eventSignalD[,4:8]),
         RNA_refSig = rowMeans(eventSignalD[,9:13])) %>% 
  select(sample,cnvID_label,ATAC_refSig,RNA_refSig)
saveRDS(refsig,"intmd/Baseline_SumTotalCov_Event.rds")
##### CNV event baseline signal plot #####
scatterD <- readRDS("intmd/scatterD_cnv.rds")
refsig <- readRDS("intmd/Baseline_SumTotalCov_Event.rds")
atac_scatter <- scatterD %>% 
  inner_join(refsig,by=c("sample","cnvID_label")) %>% 
  ggscatter(x = "ATAC_refSig", y = "recall",
          shape="mode",
          color = "detect_cat",
          size =  "revSize")+
  geom_hline(yintercept = 0.8,color="gray50")+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  scale_size_continuous(
    range = c(3,5),
    guide="none"
  )+
  scale_shape_manual(labels=mode_label,
                     values=mode_shapes)+
  scale_color_manual(values=levCol,guide="none")+
  theme_bw()+
  rremove("legend.title")+
  theme(
    axis.line = element_line(colour = "black",size=1),
    axis.ticks = element_line(size=1,color="black"),
    axis.text = element_text(color="black"),
    axis.ticks.length=unit(0.2,"cm"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right")+
  font("xylab",size=rel(1))+  
  font("xy",size=rel(1))+ 
  font("xy.text", size = rel(1.5)) +  
  font("legend.text",size = rel(1.2))+
  guides(shape = guide_legend(
    override.aes = list(alpha = 1,size=rel(4)),ncol = 1))+
  labs(x="Normalized total molecular abundance",
       title = "ATAC")
rna_scatter <- scatterD %>% 
  inner_join(refsig,by=c("sample","cnvID_label")) %>% 
  ggscatter(x = "RNA_refSig", y = "recall",
            shape="mode",
            color = "detect_cat",
            size =  "revSize")+
  geom_hline(yintercept = 0.8,color="gray50")+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  scale_size_continuous(
    range = c(3,5),
    guide="none"
  )+
  scale_shape_manual(labels=mode_label,
                     values=mode_shapes)+
  scale_color_manual(values=levCol,guide="none")+
  theme_bw()+
  rremove("legend.title")+
  theme(
    axis.line = element_line(colour = "black",size=1),
    axis.ticks = element_line(size=1,color="black"),
    axis.text = element_text(color="black"),
    axis.ticks.length=unit(0.2,"cm"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none")+
  font("xylab",size=rel(1))+  
  font("xy",size=rel(1))+ 
  font("xy.text", size = rel(1.5)) +  
  font("legend.text",size = rel(1.2))+
  labs(x="Normalized total molecular abundance",
       title = "RNA")
pscatter <- ggarrange(rna_scatter,atac_scatter,nrow=1,widths = c(1,1.7))
saveRDS(pscatter,"intmd/scatter_sumref.rds")
pdf("Figures/SuppFig2.pdf",width=7,height=4)
annotate_figure(pscatter, 
                top = text_grob("Summed over binned reference profiles \n overlapping with each CNV event", 
                color = "gray35", face = "bold", size = 14))
dev.off()
##### verify the weird bin coverage (not pursued) ######
hits <- findOverlaps(armdf, ref_gr)
as.data.frame(hits) %>% 
  group_by(queryHits) %>% 
  summarise(n=n()) %>% View()

varbins <- readRDS("_data/gtf_gcnt.rds")
hits <- findOverlaps(armdf, varbins)
as.data.frame(hits) %>% 
  group_by(queryHits) %>% 
  summarise(n=n()) %>% View()

load("~/Desktop/hg38_grangeslist.rda")
hg38_grangeslist$hg38_200kb
hits <- findOverlaps(armdf, hg38_grangeslist$hg38_100kb)
as.data.frame(hits) %>% 
  group_by(queryHits) %>% 
  summarise(n=n()) %>% View()

load("~/Desktop/sysdata.rda")

clonecols <- c("gray89",c("#fab81b","#06d6a0","#2176ff","#ef476f"))
names(clonecols) <- paste0("Cl",(1:5)-1)  
ggClone_plotL <- purrr::imap(clonecols[-1],\(x,i) ggplot() +
                              geom_point(data = data.frame(x = 0, y = 0), 
                                         aes(x = x, y = y), 
                                         size = rel(14), shape = 21, 
                                         fill = x, color = "black") +
                              geom_text(aes(x = 0, y = 0, label = i),olor = "black") +
                              coord_fixed()+
                              theme_void())

pdf("Figures/ggclone.pdf",width=2,height=7)
ggpubr::ggarrange(plotlist = ggClone_plotL,ncol=1)
dev.off()

