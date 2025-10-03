library(Biobase)
library(SummarizedExperiment)
library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(Matrix)
library(DelayedArray)


library(TabulaMurisSenisData)


# Bulk RNA-seq:
bulk <- TabulaMurisSenisBulk()
# Single-cell droplet data (10x):
droplet <- TabulaMurisSenisDroplet(tissues = "All")
# Single-cell FACS data (smart-seq2):
facs <- TabulaMurisSenisFACS(tissues = "All")
#see all tissues
TabulaMurisSenisFACS(tissues = NULL, infoOnly = TRUE)


# Extract the data in the bulk
bulk_metadata <- colData(bulk) %>% data.frame() #get metadata
#table(bulk_metadata$organ)
bulk_annot <- rowData(bulk) %>% data.frame()#get annotation
bulk_gene_matrix <- assay(bulk)
#bulk_gene_matrix <- assays(bulk)$counts


# Extract the data in the droplet
droplet_metadata <- colData(droplet$All) %>% data.frame() #get metadata
#table(droplet_metadata$organ)
droplet_annot <- rowData(droplet$All) %>% data.frame()#get annotation
droplet_gene_matrix <- assay(droplet$All)
#droplet_gene_matrix <- assays(droplet$All)$counts
#droplet_gene_matrix@seed


# Extract the data in the facs
facs_metadata <- colData(facs$All) %>% data.frame() #get metadata
#table(facs_metadata$organ)
facs_annot <- rowData(facs$All) %>% data.frame()#get annotation
facs_gene_matrix <- assay(facs$All)
#facs_gene_matrix <- assays(facs$All)$counts
#facs_gene_matrix@seed


#Convert the single cell experiment data to Seurat

facs_liver <- TabulaMurisSenisFACS(tissues = "Liver")
class(facs_liver)
names(facs_liver)

facs_liver <- facs_liver[["Liver"]]
x <- assay(facs_liver) %>% data.frame()
class(facs_liver)
seurat_liver <- as.Seurat(facs_liver)
seurat_liver <- as.Seurat(facs_liver, counts = "counts")
assayNames(facs_liver)
seurat_liver <- as.Seurat(facs_liver, counts = "counts", data = NULL)




# Convert entire object
seurat_facs <- as.Seurat(facs, counts = "counts")

# Or a single tissue
seurat_liver <- as.Seurat(facs_liver, counts = "counts")





sce_to_seurat <- function(sce, counts_assay = "counts", add_data = FALSE) {
  # Check assays
  if (!(counts_assay %in% assayNames(sce))) {
    stop(paste("Assay", counts_assay, "not found in SCE. Available:", 
               paste(assayNames(sce), collapse = ", ")))
  }
  
  # Extract counts
  counts <- assay(sce, counts_assay)
  
  # If DelayedArray, convert to sparse dgCMatrix (memory efficient)
  if (inherits(counts, "DelayedArray") || inherits(counts, "DelayedMatrix")) {
    message("Converting DelayedArray to sparse dgCMatrix...")
    counts <- as(counts, "dgCMatrix")
  }
  
  # Put counts back into SCE
  assay(sce, counts_assay) <- counts
  
  # Convert to Seurat
  seu <- as.Seurat(sce, counts = counts_assay, 
                   data = if (add_data) counts_assay else NULL)
  
  return(seu)
}

# Convert to Seurat
seurat_liver <- sce_to_seurat(facs_liver)

#check
seurat_liver

DimPlot(seurat_liver, group.by = "age")
DimPlot(seurat_liver, group.by = "sex")



