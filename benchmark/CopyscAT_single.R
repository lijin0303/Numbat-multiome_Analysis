#load in packages required
source("~/InHouse/copyscat_utils.R")
## set command line arguments
option_list = list(make_option("--sample_id", default = NA, type = 'character'),
                   make_option("--ctrl", default = T, type = 'logical'))
opt = parse_args(OptionParser(option_list = option_list))
sample_id = opt$sample_id
pseudo  = opt$ctrl
##### REGULAR WORKFLOW #####
initialiseEnvironment(genomeFile="Reference/hg38_chrom_sizes.tsv",
                      cytobandFile="Reference/hg38_1e+06_cytoband_densities_granges.tsv",
                      cpgFile="Reference/hg38_1e+06_cpg_densities.tsv",
                      binSize=1e6,
                      minFrags=1e4,
                      cellSuffix=c("-1","-2"),
                      lowerTrim=0.5,
                      upperTrim=0.8)
setOutputFile("Output4/",glue("{sample_id}"))
#output2 for just randomly sampled 300
zoomed_barcodes = readr::read_lines(glue("CB/{sample_id}_plasma.txt"))
scData <-readInputTable(glue("ProcInput2/{sample_id}_binfrag.tsv"))
scData <- scData[intersect(rownames(scData),zoomed_barcodes),]

if(pseudo){
  #normal_barcodes = readr::read_lines(glue("CB/Ctrl_immune.txt"))
  normal_barcodes = readr::read_lines(glue("CB/Ctrl_PC.txt"))
  scData_Ctrl <-readInputTable(glue("ProcInput2/Ctrl_binfrag2.tsv"))
  if(length(normal_barcodes)>=length(zoomed_barcodes)){
    set.seed(0)
    normal_barcodes <- sample(normal_barcodes,size=as.integer(length(zoomed_barcodes)*0.85))
    scData_Ctrl <- scData_Ctrl[intersect(rownames(scData_Ctrl),normal_barcodes),]
  }
  if(all.equal(colnames(scData),colnames(scData_Ctrl))){scData <- rbind(scData_Ctrl,scData)
  }else{message("Please check the column name of inputs!")}
}

scData_k_norm <- normalizeMatrixN(scData,
                                  logNorm = FALSE,maxZero=2000,
                                  imputeZeros = FALSE,blacklistProp = 0.8,
                                  blacklistCutoff=125,dividingFactor=1,
                                  upperFilterQuantile = 0.95)
#collapse into chromosome arm level
scData_collapse<-collapseChrom3N(scData_k_norm,summaryFunction=cutAverage,
                                 binExpand = 1,minimumChromValue = 100,logTrans = FALSE,
                                 tssEnrich = 1,
                                 logBase=2,minCPG=300,powVal=0.73) 
#apply additional filters
scData_collapse<-filterCells(scData_collapse,minimumSegments = 40,minDensity = 0.1)
#select just the normal barcodes that passed the copyscAT filtering, or itll throw it error
normal_barcodes_manual_qc = normal_barcodes[normal_barcodes %in% colnames(scData_collapse)]
median_iqr <- computeCenters(scData_collapse %>% select(chrom,normal_barcodes_manual_qc),
                             summaryFunction=cutAverage)
# setting medianQuantileCutoff to -1 and feeding non-neoplastic barcodes 
# in as normalCells can improve accuracy of CNV calls
# candidate_cnvs<-identifyCNVClusters(scData_collapse,median_iqr,useDummyCells = TRUE,
#                                     propDummy=0.25,minMix=0.01,
#                                     deltaMean = 0.03,deltaBIC2 = 0.25,
#                                     bicMinimum = 0.1, 
#                                     subsetSize=min(800,as.integer(ncol(scData_collapse)*0.9)),
#                                     fakeCellSD = 0.09, 
#                                     uncertaintyCutoff = 0.65,
#                                     summaryFunction=cutAverage,
#                                     maxClust = 4,mergeCutoff = 3,
#                                     IQRCutoff = 0.25,medianQuantileCutoff = -1,
#                                     normalCells=normal_barcodes_manual_qc)
candidate_cnvs<-identifyCNVClusters2(scData_collapse,median_iqr,useDummyCells = TRUE,
                                    propDummy=0.25,minMix=0.01,
                                    deltaMean = 0.03,deltaBIC2 = 0.25,
                                    bicMinimum = 0.1,
                                    subsetSize=min(800,as.integer(ncol(scData_collapse)*0.9)),
                                    fakeCellSD = 0.09,
                                    uncertaintyCutoff = 0.65,
                                    summaryFunction=cutAverage,
                                    maxClust = 4,mergeCutoff = 3,
                                    IQRCutoff = 0.25,medianQuantileCutoff = -1,
                                    normalCells=normal_barcodes_manual_qc)
candidate_cnvs_clean<-clusterCNV(initialResultList = candidate_cnvs,
                                 medianIQR = candidate_cnvs[[3]],minDiff=1.0) 
#to save this data you can use annotateCNV4 as per usual
final_cnv_list<-annotateCNV4(candidate_cnvs_clean, saveOutput=TRUE,
                             outputSuffix = "",sdCNV = 0.6,
                             filterResults=TRUE,filterRange=0.4)

