#CENSUS API
install.packages(
  "cellxgene.census",
  repos=c('https://chanzuckerberg.r-universe.dev', 'https://cloud.r-project.org')
)
?cellxgene.census::get_seurat

if (!requireNamespace("devtools", quietly = TRUE))
   install.packages("devtools")
devtools::install_github("montilab/OmicSignature") # Install OmicSignature from GitHub

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("SummarizedExperiment", "SingleCellExperiment", "MAST"), update = FALSE, ask = FALSE)
