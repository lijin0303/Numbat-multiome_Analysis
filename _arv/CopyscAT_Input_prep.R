sample_CBs <- split(proj$cellNames,proj$Sample)
saveRDS(sample_CBs,"Intermediate/SampleCB.rds")
purrr:::iwalk(sample_CBs,\(x,i) write.table(c("Cell_id",gsub(".*#","",x)),
                                            glue("~/CopyscAT/CB/{i}_plasma.txt"),
                                            col.names = F,row.names = F,
                                            quote=F))
for (sample in names(sample_CBs)){
  system(glue("~/CopyscAT/Refine_binfrag.sh {sample}"))}

##### Use immune cells from NBM10 #####
#~300 cells from T cell; monocyte annotation
annotR <- readRDS("~/QCr2_47/Intermediate/clusterID_annot.rds")
immune_CBs <- annotR |> 
  filter(ct1 %in% c("T","monocyte")) |> 
  filter(grepl("NBM",Sample)) |> 
  rownames_to_column("cb") |> 
  mutate(cellname = gsub(".*#","",cb)) |> 
  filter(!cellname %in% gsub(".*#","",proj$cellNames)) |> 
  filter(Sample=="NBM10")
write.table(c("Cell_id",immune_CBs$cellname),
            glue("~/CopyscAT/CB/Ctrl_immune.txt"),
            col.names = F,row.names = F,
            quote=F)
system("awk 'FNR==NR{ arr[$1]; next }$1 in arr' ~/CopyscAT/CB/Ctrl_immune.txt ~/CopyscAT/ProcInput/NBM10_copyscat.tsv \
> ~/CopyscAT/ProcInput2/Ctrl_binfrag.tsv")

##### Use subsampled plasma cells from NBM15 #####
# NBM15 is the only one NBM donors with good sample quality
# NBM 5, 9 , 10 are pretty bad
# NBM14 is fairly acceptable
# instead of random sampling; run on NBM15 with NMF cluster identification
# choose the normal cluster & remove cycling cell
nmf_normal <- read.csv("~/CopyscAT/normal/NBM15_nmf_clusters.csv")
nmf_clusters <- split(nmf_normal$Barcode,nmf_normal$nmf_results.cellAssigns)
cycling <- read.delim("~/CopyscAT/normal/NBM15_cycling_cells.tsv",header = F)
pseudo_candidates <- setdiff(nmf_clusters[["2"]],cycling[cycling$V2,"V1"])
annotR <- readRDS("~/QCr3_47/Intermediate/SampleCB.rds")
set.seed(0)
NBM15_PC <- gsub(".*#","",annotR[["NBM15"]])
otherSamples <- gsub(".*#","",unlist(annotR[setdiff(names(annotR),"NBM15")]))
uniqueCB <- setdiff(NBM15_PC,otherSamples) # to avoid CB confusion
cleanCB <- intersect(uniqueCB,pseudo_candidates)
write.table(c("Cell_id",sample(cleanCB,400)),
            glue("~/CopyscAT/CB/Ctrl_PC.txt"),
            col.names = F,row.names = F,
            quote=F)
system(paste0("awk 'FNR==NR{ arr[$1]; next }$1 in arr' ~/CopyscAT/CB/Ctrl_PC.txt ~/CopyscAT/ProcInput/NBM15_copyscat.tsv ",
" > ~/CopyscAT/ProcInput2/Ctrl_binfrag2.tsv"))

##### Use subsampled plasma cells from NBM15 #####
