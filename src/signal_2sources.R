##### Libraries and fxns #####
set.seed(1)
setwd("~/numbatm/Numbat-multiome_Analysis/")
source("utils/mini_import.R")
source("utils/vis.R")
library(ggpubr)
pseudobulk_force <- function(cnt,geno,ref,gtf,pars,consensus,diploidRegion,subtrees){
  invisible(list2env(pars,environment()))
  count_mat = numbat:::check_matrix(cnt)
  df_allele = numbat:::annotate_genes(geno, gtf)
  df_allele = numbat:::check_allele_df(df_allele)
  lambdas_ref = numbat:::check_exp_ref(ref)
  # filter for annotated genes
  genes_annotated = unique(gtf$gene) %>% 
    intersect(rownames(count_mat)) %>%
    intersect(rownames(lambdas_ref))
  count_mat = count_mat[genes_annotated,,drop=FALSE]
  lambdas_ref = lambdas_ref[genes_annotated,,drop=FALSE]
  bulk_subtrees = numbat:::make_group_bulks(
    groups = subtrees,
    count_mat = count_mat,
    df_allele = df_allele, 
    lambdas_ref = lambdas_ref,
    gtf = gtf,
    ncores = 6)
  tesbulk_subtrees <- numbat:::run_group_hmms(
    bulk_subtrees,
    t = t,
    gamma = gamma,
    alpha = 1e-4,
    nu = 1,
    min_genes = 10,
    common_diploid = FALSE,
    diploid_chroms = c(2,3,4,6,8,13,14,15,18,19))
  consensus$CHROM <- factor(consensus$CHROM,1:22)
  retested_subtrees = numbat:::retest_bulks(
    tesbulk_subtrees,
    consensus, 
    diploid_chroms = c(2,3,4,6,8,13,14,15,18,19), 
    gamma = gamma,
    min_LLR = min_LLR,
    ncores = 6
  )
  sample_bulk <- retested_subtrees %>% group_split(sample)
  return(sample_bulk)
}
clonep <- function(pop,bulkL,sample_segs,zoomin){
  c(D,baf_segs,lfc_segs)%<-% CNV_plotD(
    bulkL[[pop]],sample_segs[[pop]],zoomin)
  pop_CNV <- CNVp_base(D)+
    CNVp_exclude(zoomin)+
    CNVp_seg(lfc_segs,"logFC")+
    CNVp_seg(baf_segs,"pHF")+
    CNVp_theme()+
    theme(plot.background = element_blank(),
          plot.margin = margin(0.2,0.5,0.7,0.2, "cm"))
  return(pop_CNV)
}
##### Format fixed clone structure #####
numbatClones <- readRDS(glue('intmd/patA_cloneAssignment.rds'))[[4]] %>% 
  mutate(source=ifelse(grepl("146p",cell),"ATAC","RNA")) %>% 
  mutate(pop = clone)
cellSplit <- split(numbatClones$cell,list(numbatClones$source,numbatClones$pop))
segs=c("",
       "11b,12a",
       "1a,1b,7a,7c,7d,9a,10b,11b,12a,16b,17b",
       "1a,1b,7a,7c,7d,9a,10b,11a,11b,11c,12a,16b,17b",
       "1a,1b,5a,7a,7b,7c,7d,9a,10b,11a,11b,11c,16b,17b")
RNA_subtree <- map(1:5,\(x) list(sample=x,
                  members=segs[x],
                  cells=cellSplit[[glue("RNA.{x}")]],
                  size = length(cellSplit[[glue("RNA.{x}")]]))) %>% 
  setNames(1:5)
ATAC_subtree <- map(1:5,\(x) list(sample=x,
                                 members=segs[x],
                                 cells=cellSplit[[glue("ATAC.{x}")]],
                                 size = length(cellSplit[[glue("ATAC.{x}")]]))) %>% 
  setNames(1:5)
saveRDS(list("RNA"=RNA_subtree,"ATAC" = ATAC_subtree),"intmd/patA_bulk_subtree.rds")
##### Load in raw inputs for pseudobulk #####
count_mat <- readRDS("pat_comb/patA.rds")
df_allele <- fread("pat_comb/patA_allele_counts.tsv.gz", header=TRUE, sep="\t")
lambdas_ref <- readRDS("pat_comb/lambdas_comb_bincnt.rds") %>% as.matrix()
gtf <- readRDS("pat_comb/gtf_bincnt.rds")
pars <- readRDS("pat_comb/par_comb_bincnt.rds")
segs_consensus <- fread("pat_comb/segs_consensus_2.tsv")
diploidRegion <- c(2,3,4,6,8,13,14,15,18,19)
##### Pseudobulk RNA, ATAC separately #####
RNA_pseudo <- pseudobulk_force(count_mat,df_allele,lambdas_ref,gtf,pars,
                               segs_consensus,diploidRegion,RNA_subtree)

ATAC_pseudo <- pseudobulk_force(count_mat,df_allele,lambdas_ref,gtf,pars,
                               segs_consensus,diploidRegion,ATAC_subtree)
saveRDS(list("RNA"=RNA_pseudo,"ATAC" = ATAC_pseudo),"intmd/patA_pseudobulks.rds")
##### Load in inputs for visualizing pseudobulks #####
zoomin <- c(1,5,7,9:12,16,17)
sample_segs <- readRDS(glue('intmd/patA_cloneSegs.rds'))[[4]]
##### Visualizing pseudobulks #####
cloneL <- map(2:5,\(x) clonep(x,RNA_pseudo,sample_segs,zoomin))
RNA_clonesp <- ggarrange(plotlist = cloneL,ncol=1) 

cloneL <- map(2:5,\(x) clonep(x,ATAC_pseudo,sample_segs,zoomin))
ATAC_clonesp <- ggarrange(plotlist = cloneL,ncol=1) 

ggarrange(RNA_clonesp,ATAC_clonesp,nrow=1,labels = c("RNA","ATAC")) %>% 
  ggsave(filename="Figures/Pseudobulk_sources.pdf",width=12,height=10)