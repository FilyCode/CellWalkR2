# R Script: Make OmicSignatures and OmicSignatureCollection from Replog Dataset Phenotypes

# This script loads regression results from perturbational scRNA-seq experiments,
# specifically from the Replogle et al. 2022 dataset for K562 and RPE1 cell lines.
# For each perturbed gene, it filters the differential expression data,
# creates a signature of top responding genes, and encapsulates these
# into OmicSignature objects. Finally, it organizes these into
# OmicSignatureCollection objects for each cell line and a combined one.

# --- 1. Setup and Load Libraries ---
library(tidyverse)    # For data manipulation
library(OmicSignature) # To work with OmicSignature objects and collections
library(biomaRt)      # For gene ID mapping
# For faster gene ID mapping, use a local annotation database.
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db) # Local Human Gene annotation database


# --- 2. Define Paths and Parameters ---
data_input_path <- file.path("/restricted/projectnb/agedisease/projects/challenge2025/data/perturbational_sigs/replogle_2022")
k562_data_file <- file.path(data_input_path, "k562_ps_sig_all.rds")
rpe1_data_file <- file.path(data_input_path, "rpe1_ps_sig_all.rds")

output_base_path <- file.path(data_input_path)

k562_collection_output_file <- file.path(output_base_path, "Replogle_K562_Perturb_OmicSignatureCollection.rds")
rpe1_collection_output_file <- file.path(output_base_path, "Replogle_RPE1_Perturb_OmicSignatureCollection.rds")
combined_collection_output_file <- file.path(output_base_path, "Replogle_Perturb_Combined_OmicSignatureCollection.rds")
gene_mapping_cache_file <- file.path(output_base_path, "gene_symbol_to_ensembl_map.rds") # Cache file for gene mapping

# Define filter parameters
adj_p_cutoff <- 0.05                # Adjusted p-value cutoff for significant genes in signature
log2fc_abs_cutoff <- 0.25           # Absolute log2FC cutoff for significant genes in signature
max_genes_in_signature <- 500       # Max. number of significant genes saved in the signature part of the OmicSignature object


# --- 3. Load Data Files ---
message("--- Loading Perturb-seq Data ---")
k562_ps_sig_all <- readRDS(k562_data_file)
rpe1_ps_sig_all <- readRDS(rpe1_data_file)

message(paste0("Loaded K562 data: ", length(k562_ps_sig_all), " perturbation experiments."))
message(paste0("Loaded RPE1 data: ", length(rpe1_ps_sig_all), " perturbation experiments."))


# --- 4. Gene Symbol to Ensembl ID Mapping (Optimized with org.Hs.eg.db + biomaRt fallback) ---
message("\n--- Mapping Gene Symbols to Ensembl IDs ---")

# Extract all unique gene symbols
all_gene_symbols <- unique(c(
  unlist(lapply(k562_ps_sig_all, rownames)),
  unlist(lapply(rpe1_ps_sig_all, rownames))
))
all_gene_symbols <- all_gene_symbols[all_gene_symbols != "" & !is.na(all_gene_symbols)]
message(paste0("Found ", length(all_gene_symbols), " unique gene symbols across all datasets."))

# Check if a cached mapping exists and load it
gene_symbol_to_ensembl <- NULL
if (file.exists(gene_mapping_cache_file)) {
  message("  Loading gene mapping from cache...")
  gene_symbol_to_ensembl <- readRDS(gene_mapping_cache_file)
  # Filter to only symbols needed for this run
  gene_symbol_to_ensembl <- gene_symbol_to_ensembl[names(gene_symbol_to_ensembl) %in% all_gene_symbols]
  # Identify any new symbols not in cache
  symbols_to_map <- setdiff(all_gene_symbols, names(gene_symbol_to_ensembl))
} else {
  symbols_to_map <- all_gene_symbols
}


if (length(symbols_to_map) > 0) {
  # --- Step 1: Map using org.Hs.eg.db (local and fast) ---
  message(paste0("  Attempting to map ", length(symbols_to_map), " symbols using org.Hs.eg.db..."))
  ensembl_ids_from_db <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = symbols_to_map,
    column = "ENSEMBL",
    keytype = "SYMBOL",
    multiVals = "first" # Take the first Ensembl ID if a symbol maps to multiple
  )
  
  # Remove NAs from the result (symbols not found in org.Hs.eg.db)
  ensembl_ids_from_db <- ensembl_ids_from_db[!is.na(ensembl_ids_from_db)]
  
  # Update the master mapping
  if (is.null(gene_symbol_to_ensembl)) {
    gene_symbol_to_ensembl <- ensembl_ids_from_db
  } else {
    gene_symbol_to_ensembl <- c(gene_symbol_to_ensembl, ensembl_ids_from_db)
  }
  
  mapped_by_db_count <- length(ensembl_ids_from_db)
  symbols_remaining_for_biomart <- setdiff(symbols_to_map, names(ensembl_ids_from_db))
  
  message(paste0("  Mapped ", mapped_by_db_count, " symbols using org.Hs.eg.db."))
  
  # --- Step 2: Fallback to biomaRt for any remaining unmapped symbols ---
  if (length(symbols_remaining_for_biomart) > 0) {
    message(paste0("  Attempting to map ", length(symbols_remaining_for_biomart), " remaining symbols using biomaRt (EnsemblIDs v114)..."))
    
    # Connect to Ensembl BioMart
    ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 114)
    
    gene_id_map_df_biomart <- getBM(
      attributes = c("hgnc_symbol", "ensembl_gene_id"),
      filters = "hgnc_symbol",
      values = symbols_remaining_for_biomart,
      mart = ensembl
    )
    
    ensembl_ids_from_biomart <- setNames(gene_id_map_df_biomart$ensembl_gene_id, gene_id_map_df_biomart$hgnc_symbol)
    
    # Update the master mapping
    gene_symbol_to_ensembl <- c(gene_symbol_to_ensembl, ensembl_ids_from_biomart)
    message(paste0("  Mapped ", length(ensembl_ids_from_biomart), " symbols using biomaRt."))
  }
  
  # Save the updated complete mapping for future runs
  saveRDS(gene_symbol_to_ensembl, file = gene_mapping_cache_file)
  message(paste0("  Saved updated gene mapping cache to '", gene_mapping_cache_file, "'."))
  
} else {
  message("  All gene symbols already mapped and present in cache.")
}

final_mapped_count <- length(unique(names(gene_symbol_to_ensembl)))
unmapped_overall_count <- length(all_gene_symbols) - final_mapped_count

message(paste0("Total successfully mapped symbols: ", final_mapped_count, " out of ", length(all_gene_symbols), "."))
if (unmapped_overall_count > 0) {
  message(paste0("  (Note: ", unmapped_overall_count, " gene symbols could not be mapped to Ensembl IDs and will use the gene symbol as ID)."))
}


# --- 5. Helper Function to Create OmicSignature for a Single Perturbation ---
create_perturb_omic_signature <- function(
    perturbation_df, perturbation_gene_symbol, cell_line,
    gene_symbol_to_ensembl_map,
    adj_p_cutoff, log2fc_abs_cutoff, max_genes_in_signature
) {
  message(paste0("  Processing signature for ", perturbation_gene_symbol, " in ", cell_line, "..."))
  
  # Add 'gene_symbol' column from rownames and convert to tibble for easier manipulation
  current_results_df <- perturbation_df %>%
    rownames_to_column("gene_symbol") %>%
    as_tibble()
  
  # Apply gene ID mapping
  current_results_df <- current_results_df %>%
    dplyr::mutate(
      ensembl_id = gene_symbol_to_ensembl_map[gene_symbol],
      # probe_id: unique identifier, use Ensembl ID if available, otherwise gene_symbol
      probe_id = ifelse(is.na(ensembl_id) | ensembl_id == "", gene_symbol, ensembl_id),
      # feature_name: human-readable feature name, Ensembl ID if available, otherwise gene_symbol
      feature_name = ifelse(is.na(ensembl_id) | ensembl_id == "", gene_symbol, ensembl_id)
    )
  
  # Remove rows that might have originated from unmapped or empty gene symbols if any.
  current_results_df <- current_results_df %>%
    filter(!is.na(gene_symbol) & gene_symbol != "" & !is.na(probe_id) & probe_id != "")
  
  # If no genes remain after mapping, skip this signature
  if (nrow(current_results_df) == 0) {
    message(paste0("    Skipped: No valid genes after ID mapping for ", perturbation_gene_symbol, " (", cell_line, ")."))
    return(NULL)
  }
  
  # --- Prepare 'difexp' dataframe for OmicSignature object ---
  # Use 'avg_log2FC' as the primary 'score' for OmicSignature, as it represents magnitude and direction.
  difexp_data <- current_results_df %>%
    dplyr::mutate(
      score = avg_log2FC * p_val_adj, # This will be the main score in OmicSignature
      logfc = avg_log2FC,
      p_value = p_val,
      adj_p = p_val_adj,
      group_label = as.factor(
        ifelse(avg_log2FC > 0, "Increased_by_Perturbation", "Decreased_by_Perturbation")
      )
    ) %>%
    # Select only the columns required for the OmicSignature difexp slot
    dplyr::select(probe_id, feature_name, score, p_value, adj_p, logfc, group_label, gene_symbol)
  
  # --- Prepare 'signature' dataframe for OmicSignature object ---
  # 1. Filter genes based on defined cutoffs
  significant_genes <- difexp_data %>%
    dplyr::filter(adj_p <= adj_p_cutoff & abs(logfc) >= log2fc_abs_cutoff)
  
  if (nrow(significant_genes) == 0) {
    message(paste0("    Skipped: No significant genes found for ", perturbation_gene_symbol, " (", cell_line, ") with current cutoffs (adj.p <= ", adj_p_cutoff, ", |log2FC| >= ", log2fc_abs_cutoff, ")."))
    return(NULL)
  }
  
  # 3. Rank by absolute score and take the top N genes
  significant_genes_ranked <- significant_genes %>%
    dplyr::arrange(desc(abs(score))) %>%
    dplyr::slice_head(n = max_genes_in_signature)
  
  if (nrow(significant_genes_ranked) == 0) {
    message(paste0("    Skipped: No genes remain for signature for ", perturbation_gene_symbol, " (", cell_line, ") after ranking and limiting to ", max_genes_in_signature, " genes."))
    return(NULL)
  }
  
  # Select and format columns for the 'signature' slot
  signature_data <- significant_genes_ranked %>%
    dplyr::select(probe_id, feature_name, score, group_label) %>%
    # Ensure probe_id uniqueness, picking the first entry if duplicates exist
    distinct(probe_id, .keep_all = TRUE)
  
  if (nrow(signature_data) == 0) {
    message(paste0("    Skipped: No unique probe IDs remain for signature for ", perturbation_gene_symbol, " (", cell_line, ")."))
    return(NULL)
  }
  
  # --- Create Metadata for the OmicSignature object ---
  # Determine sample_type for metadata description
  phenotype_desc <- switch(cell_line,
                           "K562" = "chronic myeloid leukemia (CML) K562 cells",
                           "RPE1" = "retinal pigment epithelial RPE1 cells",
                           paste0(cell_line, " cells")
  )
  
  metadata_object <- OmicSignature::createMetadata(
    signature_name = paste0(perturbation_gene_symbol, " Knockdown Signature - ", cell_line),
    organism = "Homo sapiens",
    direction_type = "bi-directional",
    assay_type = "transcriptomics",
    phenotype = paste0(perturbation_gene_symbol, " knockdown"),
    author = "BU_Bioinformatics_ChallengeProject2025",
    year = as.numeric(format(Sys.Date(), "%Y")),
    platform = "single-cell RNA-seq (CRISPRi Perturb-seq)",
    sample_type = phenotype_desc,
    description = paste0(
      "Transcriptional signature showing gene expression changes upon CRISPRi-mediated knockdown of ",
      perturbation_gene_symbol, " in ", phenotype_desc, ". ",
      "Data from Replogle et al. 2022 perturbational scRNA-seq regression analysis. ",
      "Signature genes filtered by adjusted p-value <= ", adj_p_cutoff, " and absolute log2FC >= ", log2fc_abs_cutoff, ". ",
      "The top ", max_genes_in_signature, " genes were selected by ranking based on abs(avg_log2FC * adj_p)."
    ),
    adj_p_cutoff = adj_p_cutoff,
    logfc_cutoff = log2fc_abs_cutoff,
    score_cutoff = NULL, # No specific score cutoff, as ranking is by custom score, not 'score' directly
    keywords = c("Perturb-seq", "CRISPRi", "knockdown", perturbation_gene_symbol, cell_line, "scRNA-seq", "gene expression")
  )
  
  # Create the OmicSignature object
  omic_sig_obj <- tryCatch({
    OmicSignature$new(
      metadata = metadata_object,
      signature = signature_data,
      difexp = difexp_data
    )
  }, error = function(e) {
    message(paste0("    ERROR creating OmicSignature for ", perturbation_gene_symbol, " (", cell_line, "): ", e$message))
    return(NULL)
  })
  
  return(omic_sig_obj)
}


# --- 6. Process Each Perturbation to Create OmicSignatures ---

message("\n--- Generating K562 OmicSignatures ---")
k562_omicsigs_list <- list()
for (gene_sym in names(k562_ps_sig_all)) {
  omic_sig <- create_perturb_omic_signature(
    perturbation_df = k562_ps_sig_all[[gene_sym]],
    perturbation_gene_symbol = gene_sym,
    cell_line = "K562",
    gene_symbol_to_ensembl_map = gene_symbol_to_ensembl,
    adj_p_cutoff = adj_p_cutoff,
    log2fc_abs_cutoff = log2fc_abs_cutoff,
    max_genes_in_signature = max_genes_in_signature
  )
  if (!is.null(omic_sig)) {
    k562_omicsigs_list[[omic_sig$metadata$signature_name]] <- omic_sig
  }
}
message(paste0("Finished K562. Created ", length(k562_omicsigs_list), " OmicSignatures."))


message("\n--- Generating RPE1 OmicSignatures ---")
rpe1_omicsigs_list <- list()
for (gene_sym in names(rpe1_ps_sig_all)) {
  omic_sig <- create_perturb_omic_signature(
    perturbation_df = rpe1_ps_sig_all[[gene_sym]],
    perturbation_gene_symbol = gene_sym,
    cell_line = "RPE1",
    gene_symbol_to_ensembl_map = gene_symbol_to_ensembl,
    adj_p_cutoff = adj_p_cutoff,
    log2fc_abs_cutoff = log2fc_abs_cutoff,
    max_genes_in_signature = max_genes_in_signature
  )
  if (!is.null(omic_sig)) {
    rpe1_omicsigs_list[[omic_sig$metadata$signature_name]] <- omic_sig
  }
}
message(paste0("Finished RPE1. Created ", length(rpe1_omicsigs_list), " OmicSignatures."))


# --- 7. Create and Save OmicSignatureCollections ---

# Metadata for collections
collection_base_description <- paste0(
  "Collection of transcriptional signatures from CRISPRi Perturb-seq experiments. ",
  "Each signature captures gene expression changes upon knockdown of a specific gene. ",
  "Signatures were filtered by adjusted p-value <= ", adj_p_cutoff, " and absolute log2FC >= ", log2fc_abs_cutoff, ". ",
  "The top ", max_genes_in_signature, " genes were selected by ranking based on abs(avg_log2FC * adj_p)."
)
collection_base_keywords <- c("Perturb-seq", "CRISPRi", "knockdown", "scRNA-seq", "gene expression", "Replogle 2022")


# K562 OmicSignatureCollection
if (length(k562_omicsigs_list) > 0) {
  message("\n--- Creating K562 OmicSignatureCollection ---")
  k562_collection_metadata <- OmicSignature::createMetadata(
    collection_name = "Replogle K562 Perturb-seq Signatures",
    description = paste0("Signatures derived from K562 cells. ", collection_base_description),
    organism = "Homo sapiens",
    direction_type = "bi-directional",
    phenotype = "CRISPRi Perturbation in K562",
    assay_type = "transcriptomics",
    platform = "single-cell RNA-seq (CRISPRi Perturb-seq)",
    author = "BU_Bioinformatics_ChallengeProject2025",
    year = as.numeric(format(Sys.Date(), "%Y")),
    keywords = c("K562", collection_base_keywords)
  )
  
  k562_omic_collection <- OmicSignatureCollection$new(
    metadata = k562_collection_metadata,
    OmicSigList = k562_omicsigs_list
  )
  saveRDS(k562_omic_collection, file = k562_collection_output_file)
  message(paste0("Saved K562 OmicSignatureCollection with ", length(k562_omic_collection$OmicSigList), " signatures to '", k562_collection_output_file, "'."))
} else {
  message("\nNo K562 OmicSignatures were successfully created. K562 OmicSignatureCollection not generated.")
}

# RPE1 OmicSignatureCollection
if (length(rpe1_omicsigs_list) > 0) {
  message("\n--- Creating RPE1 OmicSignatureCollection ---")
  rpe1_collection_metadata <- OmicSignature::createMetadata(
    collection_name = "Replogle RPE1 Perturb-seq Signatures",
    description = paste0("Signatures derived from RPE1 cells. ", collection_base_description),
    organism = "Homo sapiens",
    direction_type = "bi-directional",
    phenotype = "CRISPRi Perturbation in RPE1",
    assay_type = "transcriptomics",
    platform = "single-cell RNA-seq (CRISPRi Perturb-seq)",
    author = "BU_Bioinformatics_ChallengeProject2025",
    year = as.numeric(format(Sys.Date(), "%Y")),
    keywords = c("RPE1", collection_base_keywords)
  )
  
  rpe1_omic_collection <- OmicSignatureCollection$new(
    metadata = rpe1_collection_metadata,
    OmicSigList = rpe1_omicsigs_list
  )
  saveRDS(rpe1_omic_collection, file = rpe1_collection_output_file)
  message(paste0("Saved RPE1 OmicSignatureCollection with ", length(rpe1_omic_collection$OmicSigList), " signatures to '", rpe1_collection_output_file, "'."))
} else {
  message("\nNo RPE1 OmicSignatures were successfully created. RPE1 OmicSignatureCollection not generated.")
}

# Combined OmicSignatureCollection
all_omicsigs_combined <- c(k562_omicsigs_list, rpe1_omicsigs_list)
if (length(all_omicsigs_combined) > 0) {
  message("\n--- Creating Combined OmicSignatureCollection ---")
  combined_collection_metadata <- OmicSignature::createMetadata(
    collection_name = "Replogle Perturb-seq Signatures (K562 & RPE1)",
    description = paste0("Combined collection of signatures from K562 and RPE1 cells. ",
                         "Each signature name specifies its cell line origin. ", collection_base_description),
    organism = "Homo sapiens",
    direction_type = "bi-directional",
    phenotype = "CRISPRi Perturbation",
    assay_type = "transcriptomics",
    platform = "single-cell RNA-seq (CRISPRi Perturb-seq)",
    author = "BU_Bioinformatics_ChallengeProject2025",
    year = as.numeric(format(Sys.Date(), "%Y")),
    keywords = c("K562", "RPE1", collection_base_keywords)
  )
  
  combined_omic_collection <- OmicSignatureCollection$new(
    metadata = combined_collection_metadata,
    OmicSigList = all_omicsigs_combined
  )
  saveRDS(combined_omic_collection, file = combined_collection_output_file)
  message(paste0("Saved Combined OmicSignatureCollection with ", length(combined_omic_collection$OmicSigList), " signatures to '", combined_collection_output_file, "'."))
} else {
  message("\nNo OmicSignatures available for K562 or RPE1. Combined OmicSignatureCollection not generated.")
}

message("\nScript finished successfully!")