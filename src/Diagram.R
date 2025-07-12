library(Gviz)
library(GenomicRanges)
gen <- "hg38"
chr <- "chr7"
geneBin <- readRDS("~/Desktop/geneBinmap.rds")
geneGTF <- readRDS("~/Desktop/gtf_gcnt.rds")
itrack <- IdeogramTrack(genome = gen, chromosome = chr)
data(geneModels)
selected_tx <- c("ENST00000242109",
                 "ENST00000345317",
                 "ENST00000457000")
geneModels%<>%filter(transcript %in% selected_tx)
# geneModels%<>%select(-gene,-exon,-transcript)
# grtrack <- GeneRegionTrack(geneModels, genome = gen,
#                            chromosome = chr, name = "Gene Model", 
#                            transcriptAnnotation = "symbol",
#                            background.panel = "#FFFEDB",
#                            background.title = "white")
grtrack <- GeneRegionTrack(geneModels, genome = gen,
                           chromosome = chr)
gtrack <- GenomceAxisTrack()
pdf("~/Desktop/grtrack.pdf",width=5,height=2)
plotTracks(list(itrack, gtrack),
           cex.main = 100, fontface.main = 100, 
           from = 26472740, to = 27000000)
dev.off()
gencode_ver <- 37
my_url <- paste0(
  "ftp://ftp.ebi.ac.uk/pub/databases/gencode/",
  "Gencode_human/release_", 
  gencode_ver,"/gencode.v", gencode_ver,
  ".annotation.gtf.gz")
my_txdb <- makeTxDbFromGFF(my_url)
my_start <- 7661779 - 5000
my_end <- 7687538 + 5000
geneTrack <- GeneRegionTrack(my_txdb, 
                             chromosome="chr17", 
                             from=my_start, to=my_end)

