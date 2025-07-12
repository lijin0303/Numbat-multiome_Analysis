combinedBin <- readRDS("../intmd/patA.rds")
library(Matrix)
library(irlba)
library(uwot)
library(igraph)
library(leidenAlg)
n_pcs <- 50  
pca_result <- irlba(t(combinedBin), nv = n_pcs, center = TRUE, scale = FALSE)
pca_scores <- pca_result$u %*% diag(pca_result$d)
dim(pca_scores)
umap_result <- umap(pca_scores, n_neighbors = 20, min_dist = 0.1, metric = "euclidean")
dim(umap_result)
# Find nearest neighbors (kNN graph)
nn <- uwot::umap(pca_scores, n_neighbors = 20, ret_nn = TRUE, n_components = 2)
knn_indices <- nn$nn$euclidean$idx
edges <- cbind(
  rep(1:nrow(knn_indices), ncol(knn_indices)),
  as.vector(knn_indices)
)
graph <- graph_from_edgelist(edges, directed = FALSE)
graph <- simplify(graph)  # remove duplicate edges
# Run Leiden clustering:
partition <- leidenAlg::leiden.community(graph, resolution = 1)
clusters <- membership(partition)
table(clusters)

plot(umap_result[,1], umap_result[,2], col = clusters,
     pch = 19, cex = 0.7, xlab = "UMAP 1", ylab = "UMAP 2",
     main = "UMAP embedding with Leiden clusters")

source("utils/mini_import.R")
metaCell <- data.table::fread("intmd/CNVclonecall_numbatmultiome.tsv") %>% 
  mutate(source=ifelse(grepl("146p",cell),"ATAC","RNA")) %>% 
  as.data.frame()
rownames(metaCell) <- metaCell$cell
metaCell <- metaCell[colnames(combinedBin),]
plot(umap_result[,1], umap_result[,2], col = factor(metaCell$source),
     pch = 19, cex = 0.7, xlab = "UMAP 1", ylab = "UMAP 2",
     main = "UMAP embedding with Leiden clusters")

  