#CENSUS API
install.packages(
  "cellxgene.census",
  repos=c('https://chanzuckerberg.r-universe.dev', 'https://cloud.r-project.org')
)
?cellxgene.census::get_seurat

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TabulaMurisSenisData")
