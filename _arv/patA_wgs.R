combinedout <- readRDS("intmd/Combined_outputs_2025-02-21.rds")
wgs_seg <- combinedout$wgs_call[["patA"]] %>% 
  select(chr=seqnames,start,end,eventType) %>% 
  mutate(chr=paste0("chr",chr))
wgs_seg1 <- wgs_seg %>% filter(start<10^8)
gr <- makeGRangesFromDataFrame(wgs_seg1 %>% filter(chr=="chr7"), keep.extra.columns = TRUE)
max_gap <- 100000
merged <- reduce(gr, with.revmap=TRUE, ignore.strand=TRUE, min.gapwidth = max_gap + 1)
merged$eventType <- sapply(merged$revmap, function(idxs) {
  subs  <- gr[idxs]
  # widths of each sub-segment
  w     <- width(subs)
  # total bp per eventType
  sums  <- tapply(w, subs$eventType, sum)
  # choose the name of the eventType with max coverage
  names(sums)[which.max(sums)]
})

# 5) back to a data.frame
final_df <- as.data.frame(merged)[, c("seqnames","start","end","eventType")]
colnames(final_df)[1] <- "chr"
wgs_df = rbind(final_df,wgs_seg %>% filter(chr=="chr7")%>% filter(start>10^8))
data.table::fwrite(wgs_df, "intmd/patA_wgs_seg.tsv")