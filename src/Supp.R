##### CNV event overview ####
setwd("~/numbat_EpiMultiome/numbat-multiome_Analysis/")
source("utils/mini_import.R")
source("utils/vis.R")
samplerun_event_eval <- readRDS("intmd/samplerun_event_eval.rds")
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
            "patA-CLL","patA-RT")
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
