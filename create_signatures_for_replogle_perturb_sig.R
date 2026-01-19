# R Script: Make OmicSignatures and OmicSignatureCollection from Replog Dataset Phenotypes

# This script loads regression results from perturbational scRNA-seq experiments,
# specifically from the Replogle et al. 2022 dataset for K562 and RPE1 cell lines.
# For each perturbed gene, it filters the differential expression data,
# creates a signature of top responding genes, and encapsulates these
# into OmicSignature objects. Finally, it organizes these into
# OmicSignatureCollection objects for each cell line and a combined one.

# --- 0. Helper Function for Named Vector Deduplication ---
# Merges named vectors, keeping the first occurrence for duplicate names.
unique_names_merge <- function(vec1, vec2) {
  combined_vec <- c(vec1, vec2)
  # Deduplicates by name, keeping the first value encountered
  tapply(combined_vec, names(combined_vec), `[`, 1)
}

# --- 1. Setup Environment and Libraries ---
#' @title Setup R environment and load necessary libraries
#' @description Loads required R packages and initializes the EnsDb object from a GTF file.
#' @param gtf_path Path to the GTF file.
#' @param output_db_path Desired path for the EnsDb SQLite database file.
#' @return An EnsDb object for gene annotation.
setup_environment <- function(gtf_path, output_db_path) {
  # Load core libraries, suppressing startup messages for cleaner output
  suppressPackageStartupMessages({
    library(tidyverse)     # For data manipulation
    library(OmicSignature) # To work with OmicSignature objects
    library(biomaRt)       # For gene ID mapping (fallback, if needed)
    library(org.Hs.eg.db)  # Local Human Gene annotation database
    library(ensembldb)     # For EnsDb functionalities
    library(gprofiler2)    # For functional enrichment (loaded, not used in this script)
    library(rentrez)       # For NCBI E-utilities (loaded, not used in this script)
    library(anndata)       # For reading H5AD files
  })
  
  message("--- Setting up R Environment and EnsDb ---")
  
  # Create local EnsDb from GTF if it doesn't already exist
  if (!file.exists(output_db_path)) {
    message(paste0("  Creating EnsDb from GTF: ", gtf_path))
    ensembldb::ensDbFromGtf(
      gtf = gtf_path,
      outfile = output_db_path,
      path = dirname(output_db_path),
      organism = "Homo sapiens",
      genomeVersion = "GRCh38",
      version = 114
    )
  } else {
    message(paste0("  EnsDb already exists at: ", output_db_path))
  }
  
  # Load the local SQLite EnsDb file
  edb <- EnsDb(output_db_path)
  metadata_edb <- metadata(edb)
  message(paste0("  EnsDb loaded: Organism=", organism(edb),
                 "; GenomeVersion=", metadata_edb[metadata_edb$name == "genome_build","value"],
                 "; EnsemblVersion=", ensemblVersion(edb)))
  return(edb)
}

# --- 2. Configuration Parameters ---
# Centralized list for all file paths and analysis parameters.
config <- list(
  # GTF and EnsDb paths
  gtf_file = "/restricted/projectnb/agedisease/projects/challenge2025/data/Homo_sapiens.GRCh38.114.gtf",
  ensdb_output_db = "/restricted/projectnb/agedisease/projects/challenge2025/data/EnsDb.Hsapiens.GRCh38.114.sqlite",
  
  # Input data paths for Perturb-seq experiments
  data_input_dir = file.path("/restricted/projectnb/agedisease/projects/challenge2025/data/perturbational_sigs/replogle_2022"),
  k562_data_file = "k562_ps_sig_all.rds",
  rpe1_data_file = "rpe1_ps_sig_all.rds",
  
  # Output paths for processed data and gene mapping cache
  output_base_dir = file.path("/restricted/projectnb/agedisease/projects/challenge2025/results/perturbational_omic_sigs/replogle_2022"),
  k562_collection_output_file = "Replogle_K562_Perturb_OmicSignatureCollection.rds",
  rpe1_collection_output_file = "Replogle_RPE1_Perturb_OmicSignatureCollection.rds",
  combined_collection_output_file = "Replogle_Perturb_Combined_OmicSignatureCollection.rds",
  gene_mapping_cache_file = "gene_symbol_to_ensembl_map.rds",
  
  # External H5AD data paths for additional gene mapping
  k562_h5ad_file = "/restricted/projectnb/agedisease/CBMrepositoryData/replogle_2022/K562_essential_raw_singlecell_01.h5ad",
  rpe1_h5ad_file = "/restricted/projectnb/agedisease/CBMrepositoryData/replogle_2022/rpe1_raw_singlecell_01.h5ad",
  
  # Filter parameters for OmicSignature objects (not used in this mapping script)
  adj_p_cutoff = 0.05,
  log2fc_abs_cutoff = 0.25,
  max_genes_in_signature = 500
)

# --- 3. Data Loading ---
#' @title Load Perturb-seq datasets
#' @description Loads K562 and RPE1 Perturb-seq data from RDS files.
#' @param config A list containing data file paths.
#' @return A list with loaded K562 and RPE1 data.
load_perturb_data <- function(config) {
  message("\n--- Loading Perturb-seq Data ---")
  k562_data <- readRDS(file.path(config$data_input_dir, config$k562_data_file))
  rpe1_data <- readRDS(file.path(config$data_input_dir, config$rpe1_data_file))
  
  message(paste0("  Loaded K562 data: ", length(k562_data), " perturbation experiments."))
  message(paste0("  Loaded RPE1 data: ", length(rpe1_data), " perturbation experiments."))
  return(list(k562 = k562_data, rpe1 = rpe1_data))
}

# --- 4. Gene Symbol to Ensembl ID Mapping ---
#' @title Map gene symbols to Ensembl IDs using a layered approach
#' @description Implements a robust strategy for mapping gene symbols to Ensembl IDs,
#'   using cached mappings, org.Hs.eg.db, EnsDb, and external H5AD files.
#' @param all_gene_symbols A character vector of unique gene symbols to map.
#' @param edb An EnsDb object for gene annotation.
#' @param config A list containing configuration parameters, especially paths for cache and H5AD files.
#' @return A named character vector where names are gene symbols and values are Ensembl IDs.
map_gene_symbols_to_ensembl <- function(all_gene_symbols, edb, config) {
  message("\n--- Mapping Gene Symbols to Ensembl IDs ---")
  # Clean up gene symbols: remove empty strings and NAs
  all_gene_symbols <- all_gene_symbols[all_gene_symbols != "" & !is.na(all_gene_symbols)]
  message(paste0("  Found ", length(all_gene_symbols), " unique gene symbols across all datasets."))
  
  gene_mapping_cache_file <- file.path(config$output_base_dir, config$gene_mapping_cache_file)
  current_gene_map <- character(0) # Initialize an empty named vector for accumulating mappings
  
  # Load existing mapping from cache to avoid re-mapping known symbols
  if (file.exists(gene_mapping_cache_file)) {
    message("  Loading gene mapping from cache...")
    cached_map <- readRDS(gene_mapping_cache_file)
    # Filter cached map to only include symbols relevant for the current run
    current_gene_map <- cached_map[names(cached_map) %in% all_gene_symbols]
    message(paste0("  ", length(current_gene_map), " genes mapped from cache."))
  }
  
  # Identify symbols that still need mapping
  symbols_to_map <- setdiff(all_gene_symbols, names(current_gene_map))
  
  if (length(symbols_to_map) > 0) {
    message(paste0("  ", length(symbols_to_map), " symbols still need mapping. Starting layered approach."))
    
    # Step 1: Map using org.Hs.eg.db (direct SYMBOL -> ENSEMBL)
    mapped_by_orghs <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = symbols_to_map,
                                             column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
    mapped_by_orghs <- mapped_by_orghs[!is.na(mapped_by_orghs)]
    current_gene_map <- unique_names_merge(current_gene_map, mapped_by_orghs)
    message(paste0("  Mapped ", length(mapped_by_orghs), " symbols using org.Hs.eg.db (direct)."))
    
    symbols_to_map <- setdiff(all_gene_symbols, names(current_gene_map)) # Update remaining symbols
    
    # Step 2: Resolve aliases using org.Hs.eg.db (ALIAS -> SYMBOL -> ENSEMBL)
    if (length(symbols_to_map) > 0) {
      message(paste0("  Attempting to resolve aliases for ", length(symbols_to_map), " remaining symbols."))
      alias_to_symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = symbols_to_map,
                                               keytype = "ALIAS", column = "SYMBOL", multiVals = "first")
      alias_to_symbol <- alias_to_symbol[!is.na(alias_to_symbol)]
      
      if (length(alias_to_symbol) > 0) {
        # Map these canonical symbols to Ensembl IDs
        canonical_symbols <- unique(alias_to_symbol)
        mapped_canonical_to_ensembl <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = canonical_symbols,
                                                             column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
        mapped_canonical_to_ensembl <- mapped_canonical_to_ensembl[!is.na(mapped_canonical_to_ensembl)]
        
        # Reconstruct ALIAS -> ENSEMBL map from successful resolutions
        resolved_alias_map <- character(0)
        for (alias_key in names(alias_to_symbol)) {
          canonical_sym <- alias_to_symbol[[alias_key]]
          if (canonical_sym %in% names(mapped_canonical_to_ensembl)) {
            resolved_alias_map[alias_key] <- mapped_canonical_to_ensembl[[canonical_sym]]
          }
        }
        resolved_alias_map <- resolved_alias_map[!is.na(resolved_alias_map)]
        current_gene_map <- unique_names_merge(current_gene_map, resolved_alias_map)
        message(paste0("  Resolved and mapped ", length(resolved_alias_map), " symbols via aliases using org.Hs.eg.db."))
      }
    }
    
    symbols_to_map <- setdiff(all_gene_symbols, names(current_gene_map)) # Update remaining symbols
    
    # Step 3: Map using local EnsDb v114 (SYMBOL -> GENEID)
    if (length(symbols_to_map) > 0) {
      message(paste0("  Attempting to map ", length(symbols_to_map), " symbols using local EnsDb (v114)."))
      res_ensdb <- ensembldb::select(edb, keys = symbols_to_map, keytype = "SYMBOL", columns = c("GENEID", "SYMBOL"))
      ensdb_map <- setNames(res_ensdb$GENEID, res_ensdb$SYMBOL)
      ensdb_map <- ensdb_map[!is.na(ensdb_map) & ensdb_map != ""]
      current_gene_map <- unique_names_merge(current_gene_map, ensdb_map)
      message(paste0("  Mapped ", length(ensdb_map), " symbols using EnsDb."))
    }
    
    symbols_to_map <- setdiff(all_gene_symbols, names(current_gene_map)) # Update remaining symbols
    
    # Step 4: Map using external H5AD files (SYMBOL -> ENSG)
    if (length(symbols_to_map) > 0) {
      message(paste0("  Attempting to map ", length(symbols_to_map), " symbols using external h5ad data."))
      k562_h5ad <- anndata::read_h5ad(config$k562_h5ad_file)
      rpe1_h5ad <- anndata::read_h5ad(config$rpe1_h5ad_file)
      
      # Extract and combine mappings from H5AD files
      map_k562 <- data.frame(symbol = as.character(k562_h5ad$obs$gene), ensg = as.character(k562_h5ad$obs$gene_id), stringsAsFactors = FALSE)
      map_rpe1 <- data.frame(symbol = as.character(rpe1_h5ad$obs$gene), ensg = as.character(rpe1_h5ad$obs$gene_id), stringsAsFactors = FALSE)
      
      combined_h5ad_map_df <- unique(rbind(map_k562, map_rpe1))
      # Remove rows with NA/empty symbol or Ensembl ID
      combined_h5ad_map_df <- combined_h5ad_map_df[!(is.na(combined_h5ad_map_df$symbol) | combined_h5ad_map_df$symbol == "" |
                                                       is.na(combined_h5ad_map_df$ensg)   | combined_h5ad_map_df$ensg == ""), ]
      # Create one-to-one mapping, preferring first ENSG per symbol
      h5ad_symbol_to_ensg <- tapply(combined_h5ad_map_df$ensg, combined_h5ad_map_df$symbol, function(x) unique(x)[1])
      
      # Apply this map to the currently remaining symbols
      mapped_by_h5ad <- h5ad_symbol_to_ensg[symbols_to_map]
      mapped_by_h5ad <- mapped_by_h5ad[!is.na(mapped_by_h5ad)]
      current_gene_map <- unique_names_merge(current_gene_map, mapped_by_h5ad)
      message(paste0("  Mapped ", length(mapped_by_h5ad), " symbols using h5ad data."))
      
      # Report symbols still unmapped after H5AD attempt
      unmapped_h5ad_post <- symbols_to_map[!symbols_to_map %in% names(mapped_by_h5ad)]
      if (length(unmapped_h5ad_post) > 0) {
        warning(paste0("  ", length(unmapped_h5ad_post), " symbols were not found in the combined h5ad mapping (e.g., ",
                       paste(head(unmapped_h5ad_post, 5), collapse = ", "),
                       if (length(unmapped_h5ad_post) > 5) paste0(" ... (+", length(unmapped_h5ad_post) - 5, " more)")))
      }
    }
    
    # Step 5: Strip version suffix from Ensembl IDs (e.g., ENSG00000123456.10 -> ENSG00000123456)
    current_gene_map <- sub("\\.\\d+$","", current_gene_map)
    
    # Save updated cache to disk
    saveRDS(current_gene_map, file = gene_mapping_cache_file)
    message(paste0("  Saved updated gene mapping cache to '", gene_mapping_cache_file, "'."))
    
  } else {
    message("  All gene symbols already mapped and present in cache. No new mapping performed.")
  }
  return(current_gene_map)
}


# --- Main Script Execution ---
# This section executes the primary workflow of the script.

# 1. Setup Environment: Load libraries and initialize EnsDb object
ens_db_obj <- setup_environment(config$gtf_file, config$ensdb_output_db)

# 2. Load Perturb-seq Data for K562 and RPE1 cell lines
perturb_data <- load_perturb_data(config)

# 3. Extract all unique gene symbols from the loaded datasets
all_gene_symbols_k562 <- unlist(lapply(perturb_data$k562, rownames))
all_gene_symbols_rpe1 <- unlist(lapply(perturb_data$rpe1, rownames))
all_unique_symbols <- unique(c(all_gene_symbols_k562, all_gene_symbols_rpe1))

# 4. Perform Gene Symbol to Ensembl ID Mapping using the layered approach
gene_symbol_to_ensembl <- map_gene_symbols_to_ensembl(all_unique_symbols, ens_db_obj, config)

# 5. Final Mapping Summary: Report statistics on mapped and unmapped symbols
final_mapped_count <- length(gene_symbol_to_ensembl)
overall_symbols_count <- length(all_unique_symbols)
unmapped_overall_count <- overall_symbols_count - final_mapped_count

message("\n--- Gene Mapping Summary ---")
message(paste0("  Total unique gene symbols requested: ", overall_symbols_count))
message(paste0("  Successfully mapped to Ensembl IDs: ", final_mapped_count))
message(paste0("  Symbols remaining unmapped: ", unmapped_overall_count))

if (unmapped_overall_count > 0) {
  unmapped_symbols <- setdiff(all_unique_symbols, names(gene_symbol_to_ensembl))
  message(paste0("  First 10 unmapped symbols: ", paste(head(unmapped_symbols, 10), collapse = ", ")))
}




# Create signatures for perturbation data 
k562_ps_sig_all <- perturb_data$k562
rpe1_ps_sig_all <- perturb_data$rpe1

adj_p_cutoff <- 0.05                # Adjusted p-value cutoff for significant genes in signature
log2fc_abs_cutoff <- 0.25           # Absolute log2FC cutoff for significant genes in signature
max_genes_in_signature <- 500       # Max. number of significant genes saved in the signature part of the OmicSignature object

output_base_path <- file.path("/restricted/projectnb/agedisease/projects/challenge2025/results/perturbational_omic_sigs/replogle_2022")

k562_collection_output_file <- file.path(output_base_path, "Replogle_K562_Perturb_OmicSignatureCollection.rds")
rpe1_collection_output_file <- file.path(output_base_path, "Replogle_RPE1_Perturb_OmicSignatureCollection.rds")
combined_collection_output_file <- file.path(output_base_path, "Replogle_Perturb_Combined_OmicSignatureCollection.rds")


# --- Helper Function to Create OmicSignature for a Single Perturbation ---
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

