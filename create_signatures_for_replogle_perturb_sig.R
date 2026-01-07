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
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")

# Create local gene annotation mapping
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("ensembldb")
library(org.Hs.eg.db) # Local Human Gene annotation database
library(ensembldb)

gtf <- "/restricted/projectnb/agedisease/projects/challenge2025/data/Homo_sapiens.GRCh38.114.gtf"
outdb <- "/restricted/projectnb/agedisease/projects/challenge2025/data/EnsDb.Hsapiens.GRCh38.114.sqlite"
DBfile <- ensembldb::ensDbFromGtf(
  gtf = gtf,
  outfile = outdb,
  path = dirname(outdb),
  organism = "Homo sapiens",
  genomeVersion = "GRCh38",
  version = 114
)

# Load the local SQLite EnsDb file to annotate the gene symbols
edb <- EnsDb(outdb)
metadata = metadata(edb)
message("Created EnsDb with organism: ", organism(edb),
        "; genomeVersion: ", metadata[metadata$name == "genome_build","value"],
        "; ensemblVersion: ", ensemblVersion(edb))

# --- 2. Define Paths and Parameters ---
data_input_path <- file.path("/restricted/projectnb/agedisease/projects/challenge2025/data/perturbational_sigs/replogle_2022")
k562_data_file <- file.path(data_input_path, "k562_ps_sig_all.rds")
rpe1_data_file <- file.path(data_input_path, "rpe1_ps_sig_all.rds")

output_base_path <- file.path("/restricted/projectnb/agedisease/projects/challenge2025/results/perturbational_omic_sigs/replogle_2022")

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
  message(paste0("  Loaded ", length(gene_symbol_to_ensembl), " genes mapped from cache..."))
} else {
  symbols_to_map <- all_gene_symbols
}


if (length(symbols_to_map) > 0) {
  # Step 1: Map using org.Hs.eg.db (local and fast)
  message(paste0("  Attempting to map ", length(symbols_to_map), " symbols using org.Hs.eg.db..."))
  ensembl_ids_from_db <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = symbols_to_map,
    column = "ENSEMBL",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  ensembl_ids_from_db <- ensembl_ids_from_db[!is.na(ensembl_ids_from_db)]
  if (is.null(gene_symbol_to_ensembl)) {
    gene_symbol_to_ensembl <- ensembl_ids_from_db
  } else {
    gene_symbol_to_ensembl <- c(gene_symbol_to_ensembl, ensembl_ids_from_db)
  }
  mapped_by_db_count <- length(ensembl_ids_from_db)
  
  # Resolve aliases/deprecated symbols to boost mapping
  symbols_remaining <- setdiff(symbols_to_map, names(ensembl_ids_from_db))
  message(paste0("  Mapped ", mapped_by_db_count, " symbols using org.Hs.eg.db. ",
                 "Resolving aliases for ", length(symbols_remaining), " remaining..."))
  alias2symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = symbols_remaining,
    keytype = "ALIAS",
    column = "SYMBOL",
    multiVals = "first"
  )
  symbols_resolved <- ifelse(!is.na(alias2symbol), alias2symbol, symbols_remaining)
  
  # Map remaining using EnsDb v114 offline
  if (length(symbols_resolved) > 0) {
    message(paste0("  Attempting to map ", length(symbols_resolved), " symbols using EnsDb.Hsapiens.v114..."))
    res <- ensembldb::select(
      edb,
      keys = symbols_resolved,
      keytype = "SYMBOL",
      columns = c("GENEID", "SYMBOL")
    )
    # Build named vector SYMBOL -> GENEID
    ensdb_map <- setNames(res$GENEID, res$SYMBOL)
    ensdb_map <- ensdb_map[!is.na(ensdb_map) & ensdb_map != ""]
    
    # Prefer mappings for symbols that were originally requested
    ensdb_map <- ensdb_map[names(ensdb_map) %in% symbols_remaining]
    
    # Merge and deduplicate (prefer existing org.Hs.eg.db mappings)
    gene_symbol_to_ensembl <- c(gene_symbol_to_ensembl, ensdb_map)
    gene_symbol_to_ensembl <- tapply(gene_symbol_to_ensembl, names(gene_symbol_to_ensembl), `[`, 1)
  }
  
  # Strip version suffix from Ensembl IDs (e.g., ENSG... .xx)
  strip_ver <- function(x) sub("\\.\\d+$","", x)
  gene_symbol_to_ensembl <- strip_ver(gene_symbol_to_ensembl)
  
  # Save cache
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
  
  # Add 'gene_symbol' column from rownames and ensure it's a plain data.frame early
  current_results_df <- as.data.frame(perturbation_df) %>%
    rownames_to_column("gene_symbol") %>%
    as_tibble() # Keeping as_tibble here for dplyr operations, but convert back later
  
  # Apply gene ID mapping
  current_results_df <- current_results_df %>%
    dplyr::mutate(
      ensembl_id = gene_symbol_to_ensembl_map[gene_symbol],
      probe_id = ifelse(is.na(ensembl_id) | ensembl_id == "", gene_symbol, ensembl_id),
      feature_name = ifelse(is.na(ensembl_id) | ensembl_id == "", gene_symbol, ensembl_id)
    ) %>% # Ensure character type for probe_id and feature_name
    dplyr::mutate(
      probe_id = as.character(probe_id),
      feature_name = as.character(feature_name),
      gene_symbol = as.character(gene_symbol)
    )
  
  # Remove rows that might have originated from unmapped or empty gene symbols if any.
  current_results_df <- current_results_df %>%
    dplyr::filter(!is.na(gene_symbol) & gene_symbol != "" & !is.na(probe_id) & probe_id != "")
  
  # If no genes remain after mapping, skip this signature
  if (nrow(current_results_df) == 0) {
    message(paste0("    Skipped: No valid genes after ID mapping for ", perturbation_gene_symbol, " (", cell_line, ")."))
    return(NULL)
  }
  
  # --- Prepare 'difexp' dataframe for OmicSignature object ---
  # Create the dataframe, perform all calculations and explicit type coercions,
  # and ensure non-finite values are handled, all within a pipeline that results in a plain data.frame.
  difexp_data_clean <- current_results_df %>%
    dplyr::mutate(
      # Ensure initial avg_log2FC and p_val_adj are explicitly numeric and handle non-finite
      temp_avg_log2FC = as.numeric(avg_log2FC),
      temp_p_val_adj = as.numeric(p_val_adj),
      temp_p_val = as.numeric(p_val),
      
      # Handle non-finite values for logFC, adj_p, and p_value as specified:
      logFC = ifelse(!is.finite(temp_avg_log2FC), 0, temp_avg_log2FC),
      adj_p = ifelse(!is.finite(temp_p_val_adj), 1, temp_p_val_adj),
      p_value = ifelse(!is.finite(temp_p_val), 1, temp_p_val)
    ) %>%
    dplyr::mutate(
      # Calculate raw score using cleaned logFC and adj_p
      raw_score_calc = logFC * (-log(adj_p)),
      # Final 'score': replace any non-finite raw_score_calc with 0
      score = ifelse(!is.finite(raw_score_calc), 0, raw_score_calc),
      
      # Group label based on the cleaned logFC
      group_label = as.factor(
        ifelse(logFC > 0, "Increased_by_Perturbation", "Decreased_by_Perturbation")
      )
    ) %>%
    # Select final columns and convert to a plain data.frame
    dplyr::select(probe_id, feature_name, score, p_value, adj_p, logFC, group_label, gene_symbol) %>%
    as.data.frame(stringsAsFactors = FALSE) # CRITICAL: Convert to plain data.frame here
  
  # Optional: Final check for non-finite values in 'score' after all operations
  if (any(!is.finite(difexp_data_clean$score))) {
    message(paste0("    RE-WARNING: Non-finite values in 'score' after full cleaning pipeline for ", perturbation_gene_symbol, " (", cell_line, "). This is highly unexpected. Setting remaining to 0."))
    difexp_data_clean$score[!is.finite(difexp_data_clean$score)] <- 0
  }
  
  # --- Prepare 'signature' dataframe for OmicSignature object ---
  # 1. Filter genes based on defined cutoffs
  significant_genes <- difexp_data_clean %>% # Use difexp_data_clean for filtering too
    dplyr::filter(adj_p <= adj_p_cutoff & abs(logFC) >= log2fc_abs_cutoff)
  
  if (nrow(significant_genes) == 0) {
    message(paste0("    Skipped: No significant genes found for ", perturbation_gene_symbol, " (", cell_line, ") with current cutoffs (adj.p <= ", adj_p_cutoff, ", |logFC| >= ", log2fc_abs_cutoff, ")."))
    return(NULL)
  }
  
  # 2. Rank by absolute score and take the top N genes
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
    distinct(probe_id, .keep_all = TRUE) %>%
    as.data.frame(stringsAsFactors = FALSE) # Also convert signature_data to plain data.frame
  
  if (nrow(signature_data) == 0) {
    message(paste0("    Skipped: No unique probe IDs remain for signature for ", perturbation_gene_symbol, " (", cell_line, ")."))
    return(NULL)
  }
  
  # --- Create Metadata for the OmicSignature object ---
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
      "The top ", max_genes_in_signature, " genes were selected by ranking based on abs(logFC * -log(adj_p))."
    ),
    adj_p_cutoff = adj_p_cutoff,
    logfc_cutoff = log2fc_abs_cutoff,
    score_cutoff = NULL,
    keywords = c("Perturb-seq", "CRISPRi", "knockdown", perturbation_gene_symbol, cell_line, "scRNA-seq", "gene expression")
  )
  
  # Create the OmicSignature object
  omic_sig_obj <- tryCatch({
    OmicSignature$new(
      metadata = metadata_object,
      signature = signature_data,
      difexp = difexp_data_clean # Pass the aggressively cleaned data frame
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
  k562_collection_metadata <- list(
    collection_name = "Replogle_K562_PerturbSeq_Signatures",
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
  rpe1_collection_metadata <- list(
    collection_name = "Replogle_RPE1_PerturbSeq_Signatures",
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
  combined_collection_metadata <- list(
    collection_name = "Replogle_PerturbSeq_Signatures_(K562_and_RPE1)",
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
