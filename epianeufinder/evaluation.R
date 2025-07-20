# use the min cell as low threshold to use th detected regions as aneu and comapre
library(dplyr)
eval_MM <- readRDS("~/numbat_EpiMultiome/Numbat-multiome_Analysis/intmd/Combined_outputs.rds")
samples2check <- names(eval_MM$wgs_call)[1:7]
inputVec <-  readRDS("~/earlyChromMM/Analysis/Annotation/sampleTab/latest_file2profile.rds") |> 
  filter((SampleID %in% samples2check)) |> 
  pull(cnt_input)
writeLines(inputVec,"epiAneufinder_inputs.txt")
