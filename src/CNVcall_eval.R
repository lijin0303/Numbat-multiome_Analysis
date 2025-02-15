setwd("~/numbat_EpiMultiome/numbat-multiome_Analysis")
# setwd("~/numbat_EpiMultiome/manuscript/")
source("mini_import.R")
##### Evaluate per numbat Run #####
source("utils/eval.R")
invisible(list2env(readRDS("Combined_outputs.rds"),environment()))
wgsCall_gr <- map(wgs_call,\(w) w[,c("seqnames","start","end")])
samplerun_eval <- map(names(wgsCall_gr),\(s){
  w <- wgsCall_gr[[s]]
  ns <- numbat_call[[s]]
  map(ns,\(n)eval_call(w,n)) %>%
    bind_rows() %>% 
    mutate(mode = names(ns)) %>% 
    mutate(sample=s)
})%>%
  bind_rows()%>%
  mutate_if(is.numeric,round,digits=2)%>% 
  replace(is.na(.), 0)
saveRDS(samplerun_eval,"intmd/samplerun_eval.rds")
##### Split the analysis by CNV types #####
wgsCall_gr <- map(wgs_call,\(w) w[,c("seqnames","start","end","eventType")])
samplerun_byCNV_eval <- map(names(wgsCall_gr),\(s){
  w <- wgsCall_gr[[s]]
  ns <- numbat_call[[s]]
  imap(ns,\(n,m) eval_call(w,n,byCNV = T) %>% mutate(mode=m)) %>%
    bind_rows() %>% 
    mutate(sample=s)})%>%
  bind_rows()%>%
  mutate_if(is.numeric,round,digits=2)%>% 
  replace(is.na(.), 0)
saveRDS(samplerun_byCNV_eval,"intmd/samplerun_CNV_eval.rds")
##### event by event evaluation #####
wgsCall_gr <- map(wgs_call,\(w) w %>% 
                    unite("CNV_ID",c("arm","annot3","eventType"),remove=F) %>% 
                    select(seqnames,start,end,CNV_ID) %>% 
                    split(as.factor(.$CNV_ID)))

cytoband_gr <- getCytobands("hg38")
cytoband_gr <- cytoband_gr[seqnames(cytoband_gr) %in% paste0("chr",1:22)]
levs <- unique(cytoband_gr$name)
gr <- imap(wgs_call,\(w,s) w %>% 
             unite("cnvID",c("arm","annot3","eventType"),remove=F) %>% 
             select(seqnames,start,end,cnvID,width,arm_perc) %>% mutate(sample=s)) %>% 
  bind_rows() %>% 
  group_by(seqnames,cnvID,sample) %>% 
  summarise(start=min(start),end=max(end)) %>% 
  mutate(seqnames=paste0("chr",seqnames)) %>%
  filter(grepl("focal",cnvID)) %>% 
  df2gr()
overlaps <- findOverlaps(gr, cytoband_gr)
annotated_focal <- as.data.frame(gr) %>% 
  mutate(queryHits = 1:30) %>% 
  inner_join(as.data.frame(overlaps),by="queryHits") %>% 
  mutate(cytobands = factor(cytoband_gr$name[subjectHits],levels=levs,ordered=TRUE)) %>% 
  group_by(sample,cnvID) %>% 
  summarise(cytoband_all = paste0(min(cytobands),"-",max(cytobands)))

event_len <- imap(wgs_call,\(w,s) w %>%     
      unite("cnvID",c("arm","annot3","eventType"),remove=F) %>% 
      select(seqnames,start,end,cnvID,width,arm_perc) %>% mutate(sample=s)) %>% 
  bind_rows() %>% 
  group_by(sample,cnvID,arm_perc) %>% 
  summarise(eventLen=sum(width))%>%
  left_join(annotated_focal,by=c("sample","cnvID")) %>% 
  separate("cnvID",c("arm","annot3","eventType"),remove=F) %>% 
  mutate(cytobandA = case_when(is.na(cytoband_all)~"",
                                  TRUE~cytoband_all)) %>% 
  mutate(arm2 = case_when(cytobandA!=""~ gsub("p|q","",arm),
                          TRUE~arm)) %>% 
  mutate(cnvID_label = paste0(arm2,cytobandA,"(",eventType,")"))

samplerun_byEvent_eval <- map(names(wgsCall_gr),\(s){
  w_cnv <- wgsCall_gr[[s]]
  ns <- numbat_call[[s]]
  sample_event <- imap(w_cnv,\(w,c){
    imap(ns,\(n,m) eval_call(w,n) %>% 
           select(recall) %>% 
           mutate(mode=m)) %>%
      bind_rows()%>% 
      mutate(cnvID=c)})%>%
    bind_rows()%>% 
    mutate(sample=s)})%>% 
  bind_rows()%>%
  mutate_if(is.numeric,round,digits=5) 
samplerun_byEvent_eval%>%
  inner_join(event_len,by=c("sample","cnvID")) %>% 
  saveRDS("intmd/samplerun_event_eval.rds")