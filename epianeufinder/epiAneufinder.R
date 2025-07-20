require(optparse)
require(epiAneufinder)
options(stringsAsFactors = F)

# Get number of available cores and set default to 2/3
available_cores <- parallel::detectCores()
default_cores <- max(1, floor(available_cores * 2/3))

option_list = list(make_option("--frag", default = NULL),
                   make_option("--blacklist", default = NULL),
                   make_option("--core", default = default_cores),
                   make_option("--outdir", default = "."),
                   make_option("--selected_cells", default = NULL))
opt = parse_args(OptionParser(option_list = option_list))
frag = opt$frag
blacklist = opt$blacklist
Ncore = opt$core
selected_cells = opt$selected_cells

cat("Running epiAneufinder with:\n")
cat("  Fragment file:", frag, "\n")
cat("  Blacklist:", blacklist, "\n")
cat("  Cores:", Ncore, "\n")
cat("  Selected cells file:", ifelse(is.null(selected_cells), "NULL", selected_cells), "\n")

#prefix <- gsub(".*\\/|\\_chrPCfragments.tsv.gz","",frag)
epiAneufinder(input=frag, 
              outdir=opt$outdir, 
              blacklist=blacklist,
              windowSize=1e5, 
              genome="BSgenome.Hsapiens.UCSC.hg38", 
              exclude=c('chrX','chrY','chrM'), 
              reuse.existing=TRUE,
              title_karyo="Karyogram", 
              ncores=Ncore,
              minFrags=20000,
              plotKaryo=FALSE,
              selected_cells=selected_cells)
