# --- Load Libraries ---
# List of packages to check and install
packages <- c("BiocManager", "NMF", "SummarizedExperiment", "matrixStats")

# Function to check if a package is installed
load_packages <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% rownames(installed.packages())) {
      print(paste0("Package already installed: ", pkg))
    } else if (pkg == "BiocManager") {
      install.packages("BiocManager")
      print("Installed BiocManager from CRAN")
    } else if (pkg %in% c("NMF", "SummarizedExperiment")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
      print(paste0("Installed from Bioconductor: ", pkg))
    } else {
      install.packages(pkg)
      print(paste0("Installed from CRAN: ", pkg))
    }
  }
  if (pkg == "BiocManager" && !requireNamespace("BiocManager", quietly = TRUE)) {
    library(BiocManager)
  }
  library(pkg, character.only = TRUE)
  print(paste0("Loaded: ", pkg))
}

# Apply function to each package
invisible(lapply(packages, load_packages))

# --- Configuration ---
expr_rds_file <- "BRCA_mrna_df.rds"
clinical_file <- "NIHMS958212-supplement-2.csv"

# Define parameters
study_filter <- "BRCA"
subtypes_to_process <- c("BRCA.LumA", "BRCA.Basal", "BRCA.Normal", "BRCA.LumB", "BRCA.Her2")
num_var_genes <- 2000        # Number of top variable genes to select
nmf_ranks <- 2:5             # Range of ranks (k) for NMF
nmf_nrun <- 10               # Number of NMF runs for each rank
nmf_seed <- 123456           # Seed for NMF reproducibility
global_seed <- 123           # Seed for overall script reproducibility
num_basis_genes <- 50        # Number of top basis genes to save per factor

# Define output directory
output_dir <- "nmf_subtype_analysis_output"

# --- Setup ---

set.seed(global_seed)

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created output directory: ", output_dir)
}

# --- Load and Prepare Data ---
# Load gene expression data
message("Loading gene expression data from: ", expr_rds_file)
stopifnot(file.exists(expr_rds_file)) # Check if file exists
gene_expression_se <- readRDS(expr_rds_file)

# --- Extract expression matrix ---
expr_mat <- NULL

if (inherits(gene_expression_se, "SummarizedExperiment")) {
  message("Input is SummarizedExperiment. Checking assays...")
  assay_names <- assayNames(gene_expression_se)
  if (length(assay_names) > 0) {
    assay_name_to_use <- assay_names[1]
    message("Attempting to extract assay: '", assay_name_to_use, "'")
    expr_mat <- tryCatch({
      assay(gene_expression_se, assay_name_to_use)
    }, error = function(e){
      warning("Failed to extract assay '", assay_name_to_use, "': ", e$message)
      return(NULL)
    })
  } else {
    warning("SummarizedExperiment contains no assays!")
  }
  
} else if (is.matrix(gene_expression_se) || is.data.frame(gene_expression_se)) {
  message("Input is matrix or data.frame. Converting to matrix.")
  expr_mat <- as.matrix(gene_expression_se)
  if(!is.matrix(expr_mat)){
    warning("Conversion to matrix failed for non-SE object.")
    expr_mat <- NULL
  }
  
} else {
  warning("Input RDS file does not contain a recognized expression data structure (e.g., SummarizedExperiment, matrix, data.frame). Class: ", class(gene_expression_se)[1])
}

# --- VALIDATION STEP ---
# Check if expr_mat was successfully created and is a matrix with dimensions
if (is.null(expr_mat)) {
  stop("Failed to extract or create a valid expression matrix ('expr_mat' is NULL). Check input file and assay names.")
} else if (!is.matrix(expr_mat)) {
  stop("Extracted 'expr_mat' is not a matrix. Class: ", class(expr_mat)[1], ". Cannot proceed.")
} else if (length(dim(expr_mat)) != 2) {
  # Handle edge case where dim might be NULL for a 0x0 matrix from as.matrix(data.frame())
  if (nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
    stop("'expr_mat' has zero rows or zero columns. Dimensions: ", nrow(expr_mat), " x ", ncol(expr_mat), ". Cannot proceed.")
  } else {
    stop("'expr_mat' does not have two dimensions. Dimensions: ", paste(dim(expr_mat), collapse=" x "), ". Cannot proceed.")
  }
} else if (nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
  stop("'expr_mat' has zero rows or zero columns. Dimensions: ", paste(dim(expr_mat), collapse=" x "), ". Cannot proceed.")
} else {
  message("Successfully obtained expression matrix with dimensions: ", paste(dim(expr_mat), collapse=" x "))
}

storage.mode(expr_mat) <- "numeric"

# Standardize sample names
if (!is.null(colnames(expr_mat))) {
  colnames(expr_mat) <- gsub("\\.", "-", colnames(expr_mat))
} else {
  warning("Column names not found in expression matrix. Skipping standardization.")
}
message("Expression data dimensions: ", nrow(expr_mat), " genes, ", ncol(expr_mat), " samples.")

# Load clinical data
message("Loading clinical data from: ", clinical_file)
stopifnot(file.exists(clinical_file))
clinical_data <- read.csv(clinical_file, header = TRUE, stringsAsFactors = FALSE)

# Filter clinical data for the relevant study
clinical_data <- clinical_data[clinical_data$TCGA.Study == study_filter, ]
message("Filtered clinical data to ", nrow(clinical_data), " samples for study: ", study_filter)

# Check required columns exist
required_clin_cols <- c("TCGA.Participant.Barcode", "TCGA.Subtype")
if (!all(required_clin_cols %in% colnames(clinical_data))) {
  stop("Clinical data is missing required columns: ", paste(setdiff(required_clin_cols, colnames(clinical_data)), collapse=", "))
}

# --- Perform NMF for each subtype ---
message("\nStarting NMF analysis for defined subtypes...")

for (i in 1:length(subtypes_to_process)) {
  
  subtype <- subtypes_to_process[i]
  
  message("\n----------------------------------------")
  message("Processing subtype: ", subtype)
  message("----------------------------------------")
  
  # Subset clinical data for the current subtype
  clinical_sub <- clinical_data[clinical_data$TCGA.Subtype == subtype, ]
  subtype_samples <- clinical_sub$TCGA.Participant.Barcode
  
  if (length(subtype_samples) == 0) {
    warning("No samples found for subtype: ", subtype, ". Skipping.")
    next
  }
  
  # Subset expression matrix for the samples in this subtype
  samples_in_expr <- intersect(subtype_samples, colnames(expr_mat))
  num_found_samples <- length(samples_in_expr)
  
  if (num_found_samples < length(subtype_samples)) {
    warning("Not all clinical samples for subtype ", subtype, " found in expression data. Using ", num_found_samples, " overlapping samples.")
  }
  
  if (num_found_samples <= max(nmf_ranks)) {
    warning("Insufficient samples found (", num_found_samples, ") for subtype: ", subtype, " after matching. Minimum required > max rank (", max(nmf_ranks),"). Skipping NMF.")
    next
  }
  
  expr_mat_sub <- expr_mat[, samples_in_expr, drop = FALSE]
  message("Subsetting expression data: ", nrow(expr_mat_sub), " genes, ", ncol(expr_mat_sub), " samples.")
  
  # --- Feature Selection: Select top n most variable genes ---
  message("Calculating gene variances...")

  gene_variances <- matrixStats::rowVars(expr_mat_sub, na.rm = TRUE)
  
  # Check for valid variances
  valid_variances <- gene_variances[!is.na(gene_variances) & gene_variances > 1e-6]
  if(length(valid_variances) == 0) {
    warning("No genes with sufficient variance found for subtype: ", subtype, ". Skipping NMF.")
    next
  }
  
  # Select top n genes
  num_genes_available <- length(valid_variances)
  num_genes_to_select <- min(num_var_genes, num_genes_available)
  
  if(num_genes_to_select < num_var_genes){
    warning("Fewer than ", num_var_genes, " genes with sufficient variance available (", num_genes_to_select,"). Using available genes.")
  }

  if(num_genes_to_select <= max(nmf_ranks)){
    warning("Number of selected genes (", num_genes_to_select, ") is less than or equal to the maximum rank (", max(nmf_ranks),"). NMF may be unstable or fail. Skipping subtype ", subtype, ".")
    next
  }
  
  high_var_genes <- names(sort(valid_variances, decreasing = TRUE))[1:num_genes_to_select]
  expr_mat_filtered <- expr_mat_sub[high_var_genes, , drop = FALSE]
  message("Filtered expression matrix to top ", nrow(expr_mat_filtered), " variable genes.")
  
  # --- Handle Negative Values ---
  min_expr_val <- min(expr_mat_filtered, na.rm = TRUE)
  if (min_expr_val < 0) {
    warning("Negative values found in filtered data for subtype: ", subtype, ". Applying offset (adding ", abs(min_expr_val), ") to make data non-negative.")
    expr_mat_filtered <- expr_mat_filtered - min_expr_val
  }
  
  # --- Handle NAs ---
  if (anyNA(expr_mat_filtered)) {
    warning("NAs detected in the final filtered matrix for subtype ", subtype, ". NMF may fail or produce unexpected results. Consider imputation if errors occur.")
  }
  
  # Final check on dimensions before running NMF
  if(nrow(expr_mat_filtered) <= max(nmf_ranks) || ncol(expr_mat_filtered) <= max(nmf_ranks)){
    warning("Final matrix dimensions too small after filtering/preprocessing for subtype: ", subtype, " (", nrow(expr_mat_filtered),"x",ncol(expr_mat_filtered), "). Skipping NMF.")
    next
  }
  
  # --- Run NMF ---
  message("Running NMF for ranks: ", paste(nmf_ranks, collapse=", "), " (nrun=", nmf_nrun, ", seed=", nmf_seed, ")")
  nmf_result <- tryCatch({
    nmf(expr_mat_filtered, rank = nmf_ranks, nrun = nmf_nrun, seed = nmf_seed, .options="tv")
  }, error = function(e) {
    warning("NMF run failed for subtype: ", subtype, ". Error: ", e$message)
    return(NULL)
  })
  
  # --- Process results ---
  if (!is.null(nmf_result)) {
    
    # Determine the class of the result
    nmf_result_class <- class(nmf_result)[1]
    message("NMF run finished. Result class: ", nmf_result_class)
    
    # Check if the result structure is usable
    fit_list <- NULL
    if (nmf_result_class == "NMF.rank" && is.list(nmf_result) && !is.null(nmf_result$fit)) {
      fit_list <- nmf_result$fit
    } else if (inherits(nmf_result, "NMFfitX")) {
      possible_ranks <- try(rank(nmf_result), silent=TRUE)
      if(!inherits(possible_ranks, "try-error")) {
        message("NMF result object seems to be NMFfitX compatible.")
      } else {
        warning("Could not determine ranks from NMF result object of class ", nmf_result_class)
      }
    } else {
      warning("NMF result object (Class: ", nmf_result_class, ") is not of a recognized type containing rank fits. Skipping output generation.")
      next
    }

    # --- Save Rank Survey Plot ---
    rank_survey_file <- file.path(output_dir, paste0(subtype, "_NMF_rank_survey.pdf"))
    pdf(rank_survey_file, width = 10, height = 7)
    tryCatch({
      if(hasMethod("plot", nmf_result_class)){
        print(plot(nmf_result))
      } else if (nmf_result_class == "NMF.rank" && !is.null(nmf_result$measures)) {
        plot(nmf_result$measures$rank, nmf_result$measures$cophenetic, type='b',
             xlab="Rank (k)", ylab="Cophenetic Coefficient", main=paste("NMF Rank Survey (Metrics) -", subtype))
        warning("Standard plot method not found for class ", nmf_result_class, ". Plotting cophenetic metric.")
      } else {
        plot.new(); title(main="Rank Survey Plot Unavailable")
        warning("Cannot generate rank survey plot for class ", nmf_result_class)
      }
    }, error = function(e) {
      warning("Failed rank survey plot generation: ", e$message)
      plot.new(); title(main="Rank Survey Plot Failed")
    })
    dev.off()
    message("Saved rank survey plot attempt: ", rank_survey_file)
    
    # --- Save Combined Consensus Map Plot ---
    consensus_map_file <- file.path(output_dir, paste0(subtype, "_consensusmap.pdf"))
    pdf(consensus_map_file, width = 15, height = 10)
    tryCatch({
      consensusmap(nmf_result, labCol = NA, labRow = NA,
                   main = paste("Consensus Maps -", subtype))
    }, error = function(e){
      warning("Consensus map generation failed: ", e$message)
      plot.new(); title(main="Consensus Map Generation Failed")
    })
    dev.off()
    message("Saved combined consensus map plot: ", consensus_map_file)
    
    
    # --- Loop through ranks to save basis genes ---
    processed_ranks <- NULL
    if (!is.null(fit_list)) {
      processed_ranks <- as.integer(names(fit_list))
    } else if (exists("possible_ranks") && !inherits(possible_ranks, "try-error")) {
      processed_ranks <- possible_ranks
    }
    
    processed_ranks <- intersect(processed_ranks, nmf_ranks)
    processed_ranks <- processed_ranks[!is.na(processed_ranks)]
    
    if(length(processed_ranks) == 0){
      warning("No valid/processed ranks found to extract basis genes for subtype ", subtype, ".")
    } else {
      message("Extracting basis genes for ranks: ", paste(processed_ranks, collapse=", "))
      
      for (k in processed_ranks) {
        message("  Extracting for rank k = ", k)
        k_char <- as.character(k)
        
        # --- Get the specific NMFfit object for rank k ---
        actual_nmf_fit <- NULL
        extract_error <- NULL
        tryCatch({
          if (!is.null(fit_list)) {
            if (k_char %in% names(fit_list)) {
              actual_nmf_fit <- fit_list[[k_char]]
            }
          } else if (inherits(nmf_result, "NMFfitX")) {
            fit_container <- nmf_result[k_char]
            if (inherits(fit_container, "NMFfit")) {
              actual_nmf_fit <- fit_container
            } else if (is.list(fit_container) && length(fit_container) == 1 && inherits(fit_container[[1]], "NMFfit")) {
              actual_nmf_fit <- fit_container[[1]]
            }
          }
        }, error = function(e) { extract_error <<- e$message })
        
        # Check if extraction was successful and object is correct type
        if (!is.null(extract_error)) {
          warning("    Error extracting fit object for rank ", k, ": ", extract_error)
          next
        } else if (is.null(actual_nmf_fit)) {
          warning("    Could not extract fit object for rank ", k, ".")
          next
        } else if (!inherits(actual_nmf_fit, "NMFfit")) {
          warning("    Extracted object for rank ", k, " is not class NMFfit (Actual: ", class(actual_nmf_fit)[1], ").")
          next
        }
        
        # --- Extract and Save Basis Genes for rank k ---
        basis_genes_file <- file.path(output_dir, paste0(subtype, "_k", k, "_basis_genes.tsv"))
        tryCatch({
          # Use extractFeatures on the validated NMFfit object
          top_features_indices <- extractFeatures(actual_nmf_fit, n = num_basis_genes)

          if (!is.list(top_features_indices)) {
            if(is.vector(top_features_indices)){
              top_features_indices <- list(Factor_1 = top_features_indices)
              warning("    extractFeatures returned a vector for rank ", k, "; wrapped in a list.")
            } else {
              stop("Unexpected format returned by extractFeatures: ", class(top_features_indices)[1])
            }
          }
          
          # Get feature names
          all_feature_names <- featureNames(actual_nmf_fit)
          if (is.null(all_feature_names) || length(all_feature_names) != nrow(basis(actual_nmf_fit))) {
            all_feature_names <- rownames(basis(actual_nmf_fit)) # Fallback
          }
          
          # Check if feature names are available
          feature_name_mode <- !is.null(all_feature_names)
          if(!feature_name_mode) { warning("    Could not get feature names for rank ", k, ". Using indices.") }
          
          # Map indices to names
          top_features_mapped <- lapply(top_features_indices, function(indices) {
            if(feature_name_mode && is.numeric(indices)){
              valid_indices <- indices[indices > 0 & indices <= length(all_feature_names)]
              if(length(valid_indices) < length(indices)) warning("    Some feature indices were out of bounds for rank ", k, ".")
              return(all_feature_names[valid_indices])
            } else {
              return(paste("Index:", indices))
            }
          })
          
          names(top_features_mapped) <- paste0("Factor_", seq_along(top_features_mapped))
          
          # Save to TSV
          con <- file(basis_genes_file, "w")
          writeLines(paste0("# Top ", num_basis_genes, " Basis Genes per Factor for Subtype: ", subtype, ", Rank: ", k), con)
          for (factor_name in names(top_features_mapped)) {
            genes_string <- paste(top_features_mapped[[factor_name]], collapse = "\t")
            cat(factor_name, "\t", genes_string, "\n", file = con, sep = "")
          }
          close(con)
          message("    Saved top basis genes to: ", basis_genes_file)
          
        }, error = function(e) {
          warning("    Failed to extract or save basis genes for rank ", k, ": ", e$message)
          if(file.exists(basis_genes_file)) {
            try(file.remove(basis_genes_file), silent = TRUE)
          }
        })
      }
    }
    
  } else {
    message("Skipping plot generation and basis gene extraction for subtype ", subtype, " due to NMF run failure or NULL result.")
  }
  
}

message("\n----------------------------------------")
message("NMF analysis complete. Results saved in: ", output_dir)
message("----------------------------------------")