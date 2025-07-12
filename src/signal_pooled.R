set.seed(1)
setwd("~/numbatm/Numbat-multiome_Analysis/")
source("utils/mini_import.R")
joint_inference <- fread("intmd/joint_post_2.tsv")
sig1 <- joint_inference %>% 
  distinct(cell,p_cnv_x,seg_label) %>%
  spread(seg_label,p_cnv_x) %>% 
  column_to_rownames("cell")
colnames(sig1) <- paste0(colnames(sig1),"_x")
sig2 <- joint_inference %>% 
  distinct(cell,p_cnv_y,seg_label) %>%
  spread(seg_label,p_cnv_y) %>% 
  column_to_rownames("cell")
colnames(sig2) <- paste0(colnames(sig2),"_y")
cnv_signals <- cbind(sig1,sig2[rownames(sig1),])

library(Matrix)
library(irlba)
library(uwot)
library(igraph)
library(leidenAlg)
pca_result <- irlba(as.matrix(cnv_signals), nv = 27, center = TRUE, scale = FALSE)
pca_scores <- pca_result$u %*% diag(pca_result$d)
umap_result <- umap(pca_scores, n_neighbors = 20, min_dist = 0.1, metric = "euclidean")

plot(umap_result[,1], umap_result[,2], col = clusters,
     pch = 19, cex = 0.7, xlab = "UMAP 1", ylab = "UMAP 2",
     main = "UMAP embedding with Leiden clusters")

metaCell <- data.table::fread("intmd/CNVclonecall_numbatmultiome.tsv") %>% 
  mutate(source=ifelse(grepl("146p",cell),"ATAC","RNA")) %>% 
  as.data.frame()
rownames(metaCell) <- metaCell$cell
metaCell <- metaCell[rownames(cnv_signals),]
umap_df <- as.data.frame(umap_result) %>% 
  set_colnames(c("umap1","umap2")) %>% 
  mutate(source=factor(metaCell$source),
         clone = factor(metaCell$clone)) %>% 
  na.omit()
require(ggpubr)
g <- ggscatter(umap_df, x = "umap1", y = "umap2", 
          color = "clone",
          palette = c("#fab81b","#06d6a0","#2176ff","#ef476f"),
          alpha=0.7,facet.by = "source")+
  theme_void()+
  theme(legend.position = "none",
        strip.text.x = element_text(size=rel(1.8)))

ggsave(filename="Figures/suppfig3.pdf",g,width=8,height=4)

