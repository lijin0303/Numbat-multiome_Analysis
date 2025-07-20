library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# ============================================================================
# SECTION 1: CNV Inference Workflow
# ============================================================================

cat("Generating CNV Inference Workflow Diagram...\n")

dot_cnv <- "
digraph {
  rankdir = TB                   // top-to-bottom flow
  nodesep=0.05; ranksep=0.15;    // ultra tight
  node [shape=box, style=solid, fontsize=18, width=2.2, height=0.5, margin=\"0.03,0.02\", penwidth=1.4, fontname=\"Arial\"]
  edge [penwidth=1.2, arrowsize=1]

  A [label = \"Phased SNP Genotype Data\"]
  B [shape=diamond, label = \"Extract Paternal\\nAllele Counts (Yj)\"]
  C [label = \"Gene Expression\\nCount Matrix\"]
  D [shape=diamond, label = \"Extract Gene Expression\\nCounts (Xi)\"]
  E [label = \"Beta-Binomial Model\"]
  F [label = \"Poisson Log-Normal Model\"]
  G [shape=diamond, label = \"Joint Hidden Markov\\nModel (HMM)\"]
  H [shape=diamond, label = \"Viterbi Algorithm\"]
  I [label = \"Allele-Specific Copy\\nNumber Inference\"]
  J [label = \"Iterative Subclonal\\nPhylogeny Reconstruction\"]
  K [label = \"Final Inferred CNVs and\\nClonal Architecture\"]

  A -> B
  C -> D
  B -> E [label = \"Allelic Imbalance\\nSignal\"]
  D -> F [label = \"Molecule Abundance\\nSignal\"]
  E -> G [label = \"Likelihood\\np(Yj|CNV state)\"]
  F -> G [label = \"Likelihood\\np(Xi|CNV state)\"]
  G -> H [label = \"15 Hidden States\\n(CNV Configurations)\"]
  G -> H [label = \"Transition Probabilities\\n(t, ps)\"]
  H -> I [label = \"Most Probable\\nCNV State\"]
  I -> J
  J -> K
}
"

g_cnv       <- grViz(dot_cnv)
svg_txt_cnv <- export_svg(g_cnv)

# Save CNV workflow as tall, crisp PNG and PDF
rsvg_png(charToRaw(svg_txt_cnv),
         file   = "CNV_workflow_diagram.png",
         width  = 2000,
         height = 5000)

rsvg_pdf(charToRaw(svg_txt_cnv),
         file   = "CNV_workflow_diagram.pdf",
         width  = 2000,
         height = 5000)

# Display the CNV diagram
print(g_cnv)

# ============================================================================
# SECTION 2: Metacell Analysis Workflow
# ============================================================================

cat("\nGenerating Metacell Analysis Workflow Diagram...\n")

dot_metacell <- "
digraph {
  rankdir = TB                   // top-to-bottom flow
  nodesep=0.05; ranksep=0.15;    // ultra tight
  node [shape=box, style=solid, fontsize=18, width=2.2, height=0.5, margin=\"0.03,0.02\", penwidth=1.4, fontname=\"Arial\"]
  edge [penwidth=1.2, arrowsize=1]

  A [label = \"ScRNA-seq Data\\n(Cell-by-Gene)\"]
  B [label = \"ScATAC-seq Data\\n(Cell-by-Fragment/Peak)\"]
  C [shape=diamond, label = \"Find Common Annotations\\n(key_to_group_by)\"]
  D [shape=diamond, label = \"Random Sample N Cells\\nper Group per Modality\"]
  E [label = \"Generate RNA Pseudocells\\n(Average Counts)\"]
  F [label = \"Generate ATAC Pseudocells\\n(Average Counts)\"]
  G [shape=diamond, label = \"Match Metacell Names\\n(meta_cell_split='_')\"]
  H [label = \"SCENICPLUS Object\\n(X_EXP + X_ACC)\"]
  I [label = \"Separate PCA\\n(RNA + ATAC)\"]
  J [label = \"Separate UMAP\\n(RNA + ATAC)\"]
  K [label = \"Integrated Analysis\\n(Matched Metacells)\"]

  A -> C [label = \"Cell Annotations\"]
  B -> C [label = \"Cell Annotations\"]
  C -> D [label = \"Common Groups\"]
  D -> E [label = \"Sampled RNA Cells\\n(nr_cells_per_metacells)\"]
  D -> F [label = \"Sampled ATAC Cells\\n(nr_cells_per_metacells)\"]
  E -> G [label = \"RNA Metacells\"]
  F -> G [label = \"ATAC Metacells\"]
  G -> H [label = \"Matched Metacells\"]
  H -> I [label = \"Separate Matrices\"]
  I -> J [label = \"Separate Embeddings\"]
  J -> K
}
"

g_metacell       <- grViz(dot_metacell)
svg_txt_metacell <- export_svg(g_metacell)

# Save metacell workflow as tall, crisp PNG and PDF
rsvg_png(charToRaw(svg_txt_metacell),
         file   = "metacell_analysis_diagram.png",
         width  = 2000,
         height = 5000)

rsvg_pdf(charToRaw(svg_txt_metacell),
         file   = "metacell_analysis_diagram.pdf",
         width  = 2000,
         height = 5000)

# Display the metacell diagram
print(g_metacell)

cat("\nBoth diagrams have been generated and saved!\n")
cat("- CNV workflow: CNV_workflow_diagram.png/pdf\n")
cat("- Metacell analysis: metacell_analysis_diagram.png/pdf\n") 