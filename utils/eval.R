##### GR reference Load in ####
require(GenomicRanges)
armF <- "_data/chrom_arm.rds"
gapF <- "_data/gaps.rds"
teloF <- "_data/telo.rds"
igF <- "_data/ig.rds"
blacklistF <- "_data/blacklist.rds"

if(!file.exists(armF)){
  arm_gr <- numbat::acen_hg38 %>% 
    gather(ind,pos,-CHROM) %>% 
    mutate(arm = case_when(ind=="end"~"p",TRUE~"q"))
  arm_gr%<>% 
    rbind(arm_gr%>% filter(arm=="p") %>% 
            select(CHROM,arm) %>%
            mutate(pos=10^4,ind="start"))%>% 
    rbind(numbat::chrom_sizes_hg38 %>% 
            select(CHROM,pos=size)%>%
            mutate(pos=pos-10^4) %>% 
            mutate(arm="q",ind="end")) %>% 
    spread(ind,pos) %>% 
    mutate(chr_arm = paste0(CHROM,arm)) %>% 
    makeGRangesFromDataFrame(keep.extra.columns = T)
  saveRDS(arm_gr,"_data/chrom_arm.rds")
}else{
  arm_gr <- readRDS(armF)
}

if(!file.exists(gapF)){
  gap_gr = numbat::gaps_hg38 %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  saveRDS(gap_gr,gapF)
}else{
  gap_gr <- readRDS(gapF)
}

if(!file.exists(teloF)){
  telo_gr <- Alleloscope::telomere.GRCh38 |> 
    filter(V8=="telomere") |> 
    select(V2,V3,V4) |> 
    mutate(V2 = gsub("chr","",V2)) %>% 
    set_colnames(c("seqnames","start","end")) |> 
    makeGRangesFromDataFrame()
  saveRDS(telo_gr,teloF)
}else{
  telo_gr <- readRDS(teloF)
}

if(!file.exists(igF)){
 IG_gr <- Reduce(rbind,list(c("IGK",2,88857361,90235368),
                  c("IGL",22,22026076,22922913),
                  c("IGH",14,105586437,106879844))) %>% 
  as.data.frame() %>% 
  set_rownames(NULL) %>% 
  set_colnames(c("gene","chr","start","end")) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)
 saveRDS(IG_gr,igF)
}else{
  IG_gr <- readRDS(igF)  
}

if(!file.exists(blacklistF)){
  blacklist_gr <- fread("_data/hg38-blacklist.v2.bed") %>% 
    mutate(V1 = gsub("chr","",V1)) %>% 
    set_colnames(c("seqnames","start","end","reason")) %>%  
    makeGRangesFromDataFrame(keep.extra.columns = T)
  saveRDS(blacklist_gr,blacklistF)
}else{
  blacklist_gr <- readRDS(blacklistF)  
}
  
##### gr utils #####
df2gr <- function(df){
  return(df %>% makeGRangesFromDataFrame(keep.extra.columns = T))
}
nearest_gr <- function(gr1,gr2){
  nearest_indices <- nearest(gr1,gr2)
  # Calculate distances
  distances <- distance(gr1, gr2[nearest_indices])
  return(distances)
}
perc_GR <- function(query, target){
  # Find overlaps
  overlaps <- findOverlaps(query, target)
  
  # Calculate the overlap width
  overlap_width <- pmin(end(query[queryHits(overlaps)]), end(target[subjectHits(overlaps)])) -
    pmax(start(query[queryHits(overlaps)]), start(target[subjectHits(overlaps)])) + 1
  
  # Add percentages for both query and target
  query_length <- width(query[queryHits(overlaps)])
  target_length <- width(target[subjectHits(overlaps)])
  
  percent_query <- (overlap_width / query_length) * 100
  percent_target <- (overlap_width / target_length) * 100
  
  # Combine results into a data frame
  result <- data.frame(
    query_id = queryHits(overlaps),
    target_id = subjectHits(overlaps),
    overlap_width = overlap_width,
    percent_query = round(percent_query,2),
    percent_target = round(percent_target,2)
  )
  return(result)
}
GR_annot <- function(gr1, gr2,...) {
  hits <- findOverlaps(gr1, gr2,select="all",ignore.strand=T,...) 
  shared <- gr1[queryHits(hits),]
  mcols(shared) <- cbind(mcols(shared),mcols(gr2)[subjectHits(hits),,drop=F])
  names(shared) <-  NULL
  return(unique(shared))
} 
##### segment comparison #####
liftover2hg38 <- function(gr){
  chFile <- "~/Data/epiConsortium/hg19ToHg38.over.chain"
  if(!file.exists(chFile)){
    hg19ToHg38 <- "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz"
    system(glue("wget -O - {hg19ToHg38} | gunzip -c > ~/Data/epiConsortium/hg19ToHg38.over.chain"))
  }
  pacman::p_load(rtracklayer)
  ch = import.chain(chFile)
  lifted_gr <- liftOver(gr, ch)
  lifted_gr <-unique(unlist(lifted_gr))
  lifted_gr$id <- 1:length(lifted_gr)
  return(lifted_gr)
}
PR_calc <- function(cnv_gr,call_gr){
  cnv_gr%<>%numbat:::subtract_ranges(gap_gr)
  call_gr%<>%numbat:::subtract_ranges(gap_gr)
  overlap_total = GenomicRanges::intersect(cnv_gr, call_gr) %>% 
    as.data.frame() %>% pull(width) %>% sum
  dna_total = cnv_gr %>% as.data.frame() %>% pull(width) %>% sum
  call_total = call_gr %>% as.data.frame() %>% pull(width) %>% sum
  pre = overlap_total/call_total
  rec = overlap_total/dna_total
  F1 <- 2*(pre*rec)/(pre+rec)
  return(data.frame(precision = pre, recall = rec,f1=F1))
}
# eval_call <- function(cnvs_dna, cnvs_call){
#   colnames(cnvs_dna)[1:3] <- c("seqnames","start","end")
#   colnames(cnvs_call)[1:3] <- c("seqnames","start","end")
#   cnv_gr = cnvs_dna %>% makeGRangesFromDataFrame(keep.extra.columns = T)
#   call_gr = cnvs_call %>% makeGRangesFromDataFrame(keep.extra.columns = T)
#   return(PR_calc(cnv_gr, call_gr))
# }

eval_call <- function(cnvs_dna, cnvs_call,byCNV=F){
  if(byCNV){
    colnames(cnvs_dna)[1:4] <- c("seqnames","start","end","cnv")
    colnames(cnvs_call)[1:4] <- c("seqnames","start","end","cnv")
    cnv_gr = map(split(cnvs_dna,cnvs_dna$cnv),\(x) df2gr(x))
    call_gr = map(split(cnvs_call,cnvs_call$cnv),\(x) df2gr(x))
    cn_states <- intersect(names(cnv_gr),names(call_gr))
    return(map(cn_states,\(c) PR_calc(cnv_gr[[c]],call_gr[[c]]) %>% mutate(cnv=c)) %>%
      bind_rows())
  }else{
    colnames(cnvs_dna)[1:3] <- c("seqnames","start","end")
    colnames(cnvs_call)[1:3] <- c("seqnames","start","end")
    cnv_gr = cnvs_dna %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    call_gr = cnvs_call %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    return(PR_calc(cnv_gr, call_gr))
  }
}


CNV_annot <- function(CNV_events,min_p=10,arm_p=90,focalSeg=F){
  arm_hg38 <- arm_gr
  CNV_armOverlap <- perc_GR(CNV_events,arm_hg38)
  CNV_armOverlap%<>% 
    mutate(arm = arm_hg38$chr_arm[target_id]) %>% 
    select(query_id,arm,seg_perc=percent_target)
  annotated_cnv <- cbind(CNV_armOverlap,as.data.frame(CNV_events)[CNV_armOverlap$query_id,])
  if(!"eventType" %in% colnames(annotated_cnv)){
    annotated_cnv %<>% 
      mutate(eventType=ifelse(cn>2,"amp","del"))
  }
  
  
  # want to implement the way that unique only when not occur different values
  arm_annot <- annotated_cnv%>% 
    group_by(arm,eventType) %>% 
    summarise(arm_perc = round(sum(seg_perc),1)) %>% 
    filter(arm_perc>min_p) %>% 
    mutate(annot = case_when(arm_perc>arm_p~"arm",
                             TRUE~"focal"))
  annotated_cnv%<>%
    inner_join(arm_annot,by=c("arm","eventType"))%>% 
    group_by(arm,eventType) %>% 
    arrange(start) %>% 
    mutate(annot2 = paste0("seg",1:n())) 
  if(focalSeg){
    annotated_cnv%<>%mutate(annot3 = case_when(annot=="arm"~annot,TRUE~annot2))
    }else{
    annotated_cnv%<>%mutate(annot3 = annot)
      } 
    
  annotated_cnv%<>% 
      unite("event",c("arm","annot3","eventType"),remove=F) %>% 
      select(-query_id) %>% 
      arrange(seqnames,event) %>% 
      as.data.frame()
  
  return(annotated_cnv)
}
cnv_coverage <- function(seg_dna,cnvs_call,low_t = 1.3,high_t = 2.5,...){
  cnvs_call <- bind_rows(cnvs_call)
  colnames(seg_dna)[1:4] <- c("seqnames","start","end","cn")
  colnames(cnvs_call)[1:3] <- c("seqnames","start","end")
  col2keep <- intersect(c("seqnames","start","end","cn","eventType"),colnames(seg_dna))
  if("eventType" %in% col2keep){
    cnvs_dna <- as.data.frame(seg_dna)[,col2keep] %>% filter(eventType!="neu")
  }else{
    cnvs_dna <- as.data.frame(seg_dna)[,col2keep] %>% filter(cn<low_t | cn > high_t)
  }
  call_gr = cnvs_call%>% makeGRangesFromDataFrame(keep.extra.columns = T)
  cnv_eventD = cnvs_dna %>% 
    makeGRangesFromDataFrame(keep.extra.columns = T) %>% 
    CNV_annot(...) %>% 
    group_nest(event) %>% 
    mutate(CNV_recall =map_dbl(data,\(cnv){
      PR_calc(df2gr(cnv),call_gr)[2]}))
  #%>%unnest_wider(CNV_recall)
  return(cnv_eventD)
}
##### clone assignment #####
shannon_M <- function(t){
  entropy <- apply(t,2,\(x) x/sum(x)) %>% apply(2,\(p) -1*sum(p*log2(p+1e-5)))
  return(mean(entropy))
}
check_cloneAssign <- function(x){
  overlapped <- inner_join(x[[1]],x[[2]],by="cell")
  ari <- mclust ::adjustedRandIndex(overlapped[,2],overlapped[,3])
  cohen_kappa <- irr::kappa2(overlapped %>% select(-cell))$value
  confusion <- table(overlapped[,2],overlapped[,3])
  avg_shannon <- shannon_M(confusion)
  return(list(metric = data.frame(ari,cohen_kappa,avg_shannon),
              compareTab = overlapped))
}
reAssign <- function(x,mapL){
  x <- factor(x)
  recode_arg <- c(mapL,list(.f=x))
  x_recoded <- do.call(forcats::fct_collapse,recode_arg)
  return(x_recoded)
}