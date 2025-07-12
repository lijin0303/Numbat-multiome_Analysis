annot <- data.table::fread("RNA_cellannot.csv",header = T)
CBs <- readLines("normal_RNA.txt")
group_CB <- annot[annot$V1 %in% CBs,]
colnames(group_CB) <- c("cell","group")
write.table(group_CB, file = "normal_RNA_annot.tsv",
            sep = "\t", row.names = F, quote = F)

annot <- data.table::fread("ATAC_cellannot.csv",header = T)
CBs <- readLines("normal_ATAC.txt")
group_CB <- annot[annot$barcode %in% CBs,c("barcode","celltype")]
colnames(group_CB) <- c("cell","group")
write.table(group_CB, file = "normal_ATAC_annot.tsv",
            sep = "\t", row.names = F, quote = F)