# ==============================================================================
# PXD069424 proteomics analysis
# ==============================================================================
# Standalone script: edit DATA_DIR and OUTPUT_DIR below, then source this file.
# Required inputs: original protein/peptide reports, S1 workbook and M3 gene table.
# Cohort: HC=36, NFHP=31, FHP=28, IPF=30; corrected labels are retained.
# Candidate selection: HC/IPF detection >=60%, |log2FC| >= log2(1.20), raw P<0.05.
# Protein tests use minimum-imputed values; boxplots/peptide tests use observed values.

# 1. Paths ----------------------------------------------------------------------
DATA_DIR <- "~/PXD069424"
OUTPUT_DIR <- "~"

PROTEIN_FILE <- file.path(DATA_DIR, "A190_DIA_20250104_Report_protein.tsv")
PEPTIDE_FILE <- file.path(DATA_DIR, "A190_DIA_20250104_Report_peptide.tsv")
S1_FILE <- file.path(DATA_DIR, "12967_2026_8643_MOESM5_ESM.xlsx")
M3_FILE <- file.path(DATA_DIR, "brown_module_high_correlation_genes.csv")

FIGURE_DIR <- OUTPUT_DIR
DIR_GS_MM <- file.path(FIGURE_DIR, "00_GS_MM")
DIR_QC <- file.path(FIGURE_DIR, "00_QC")
DIR_IPF <- file.path(FIGURE_DIR, "01_IPF_vs_HC")
DIR_FHP <- file.path(FIGURE_DIR, "02_FHP_vs_HC")
DIR_NFHP <- file.path(FIGURE_DIR, "03_NFHP_vs_HC")
DIR_FOUR_GROUP <- file.path(FIGURE_DIR, "04_FourGroup")
DIR_COMBINED_FOREST <- file.path(FIGURE_DIR, "05_CombinedDiseaseForest")

# 2. Dependencies and reproducibility -------------------------------------------
# If needed, run once:
# install.packages(c("data.table", "readxl", "ggplot2", "ggrepel", "scales"))
required_packages <- c("data.table", "readxl", "ggplot2", "ggrepel", "scales")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Install the missing R packages before running: ",
       paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})
options(stringsAsFactors = FALSE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(8848)

for (directory in c(FIGURE_DIR, DIR_GS_MM, DIR_QC, DIR_IPF, DIR_FHP,
                    DIR_NFHP, DIR_FOUR_GROUP, DIR_COMBINED_FOREST)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}


# 3. Study settings -------------------------------------------------------------
GROUP_LEVELS <- c("HC", "NFHP", "FHP", "IPF")

FULL_EXPECTED_COUNTS <- c(HC = 36, NFHP = 31, FHP = 28, IPF = 30)

SAMPLE_NAME_CORRECTIONS <- c(FHP_2958 = "NFHP_2958", FHP_3389 = "NFHP_3389")

# Legacy upstream annotation only; candidate selection is computed from data.
FINAL_11 <- c("COMP", "CXCL14", "THY1", "CTHRC1", "COL1A1", "POSTN", "COL1A2", "IGF1", "CXCL12",
  "TGFB3", "ADAM12")

CANDIDATE_DETECTION_CUTOFF <- 0.6

FOCUS_LINEAR_FC_CUTOFF <- 1.2

FOCUS_WILCOXON_P_CUTOFF <- 0.05

SELECTED_GENES <- character(0)

MISSING_RATE_CUTOFF <- 0.5

P_CUTOFF <- 0.05

FC_CUTOFF <- log2(FOCUS_LINEAR_FC_CUTOFF)

PEPTIDE_MIN_N_PER_GROUP <- 5

PEPTIDE_TOP_PER_GENE <- 4

PEPTIDE_GENE_PRIORITY <- c("FBLN1", "POSTN", "JCHAIN")

FOREST_GENE_PRIORITY <- c("FBLN1", "POSTN", "JCHAIN")

BOOT_REPS <- 2000

GROUP_COLORS <- c(HC = "#4C78A8", NFHP = "#54A24B", FHP = "#ECA82C", IPF = "#E45756")

# 4. Data processing and statistical functions ----------------------------------
check_file <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

safe_bool <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

parse_quantity_label <- function(x, level = c("PG", "PEP")) {
  level <- match.arg(level)
  y <- sub("^\\[[0-9]+\\]\\s*", "", x)
  y <- sub("^HFX_DIA_A190_", "", y)
  y <- sub(paste0("\\.raw\\.", level, "\\.Quantity$"), "", y)
  sub("^HC_", "HC", y)
}

parse_group <- function(x) {
  ifelse(grepl("^HC[0-9]+$", x), "HC", ifelse(grepl("^NFHP_[0-9]+$", x), "NFHP", ifelse(grepl("^FHP_[0-9]+$",
    x), "FHP", ifelse(grepl("^IPF_[0-9]+$", x), "IPF", NA_character_))))
}

parse_subject_id <- function(x) {
  ifelse(grepl("^HC[0-9]+$", x), sub("^HC", "", x), ifelse(grepl("^(NFHP|FHP|IPF)_[0-9]+$",
    x), sub("^.*_", "", x), NA_character_))
}

split_tokens <- function(x) {
  if (is.na(x) || !nzchar(trimws(as.character(x)))) {
    return(character(0))
  }
  z <- trimws(unlist(strsplit(as.character(x), "[;,]")))
  unique(z[nzchar(z) & !z %in% c("NA", "nan")])
}

first_token <- function(x, fallback = "") {
  z <- split_tokens(x)
  if (length(z) == 0) {
    return(fallback)
  }
  z[1]
}

contains_gene <- function(annotation_string, target_gene) {
  target_gene %in% split_tokens(annotation_string)
}

make_raw_matrix <- function(data_table, columns, sample_names) {
  values <- data_table[, columns, with = FALSE]
  values <- values[, lapply(.SD, safe_numeric)]
  matrix_object <- as.matrix(values)
  matrix_object[!is.finite(matrix_object) | matrix_object <= 0] <- NA_real_
  rownames(matrix_object) <- seq_len(nrow(matrix_object))
  colnames(matrix_object) <- sample_names
  matrix_object
}

preprocess_protein_cohort <- function(raw_matrix, group, cohort_name) {
  existing_groups <- levels(droplevels(group))
  detection_by_group <- sapply(existing_groups, function(g) {
    rowMeans(!is.na(raw_matrix[, group == g, drop = FALSE]))
  })
  if (is.null(dim(detection_by_group))) {
    detection_by_group <- matrix(detection_by_group, ncol = 1, dimnames = list(rownames(raw_matrix),
      existing_groups))
  }
  keep <- apply(detection_by_group, 1, function(x) {
    any((1 - x) < MISSING_RATE_CUTOFF)
  })
  filtered_raw <- raw_matrix[keep, , drop = FALSE]
  complete_rows <- rowSums(is.na(filtered_raw)) == 0
  complete_n <- sum(complete_rows)
  if (complete_n < 10) {
    stop(cohort_name, ": fewer than 10 complete proteins were available for normalization.")
  }
  sample_medians <- apply(filtered_raw[complete_rows, , drop = FALSE], 2, median, na.rm = TRUE)
  if (any(!is.finite(sample_medians)) || any(sample_medians <= 0)) {
    stop(cohort_name, ": invalid protein sample medians.")
  }
  global_median <- median(sample_medians, na.rm = TRUE)
  normalization_factors <- (global_median/sample_medians)
  normalized_linear <- sweep(filtered_raw, 2, normalization_factors, FUN = "*")
  global_minimum <- min(normalized_linear[is.finite(normalized_linear) & normalized_linear >
    0], na.rm = TRUE)
  imputed_linear <- normalized_linear
  imputed_linear[is.na(imputed_linear)] <- global_minimum
  list(keep = keep, filtered_raw = filtered_raw, normalized_linear = normalized_linear, observed_log2 = log2(normalized_linear),
    imputed_log2 = log2(imputed_linear), detection_by_group = detection_by_group, complete_n = complete_n,
    sample_medians = sample_medians, global_median = global_median, normalization_factors = normalization_factors,
    global_minimum = global_minimum)
}

preprocess_peptide_cohort <- function(raw_matrix, cohort_name) {
  sample_medians <- apply(raw_matrix, 2, median, na.rm = TRUE)
  if (any(!is.finite(sample_medians)) || any(sample_medians <= 0)) {
    stop(cohort_name, ": invalid peptide sample medians.")
  }
  global_median <- median(sample_medians, na.rm = TRUE)
  normalization_factors <- (global_median/sample_medians)
  normalized_linear <- sweep(raw_matrix, 2, normalization_factors, FUN = "*")
  list(normalized_linear = normalized_linear, observed_log2 = log2(normalized_linear), sample_medians = sample_medians,
    global_median = global_median, normalization_factors = normalization_factors)
}

bootstrap_se_ci <- function(control_values, case_values, estimate, reps = BOOT_REPS) {
  bootstrap_values <- replicate(reps, {
    mean(sample(case_values, length(case_values), replace = TRUE)) - mean(sample(control_values,
      length(control_values), replace = TRUE))
  })
  se_value <- stats::sd(bootstrap_values, na.rm = TRUE)
  if (!is.finite(se_value)) {
    se_value <- 0
  }
  c(estimate - 1.96 * se_value, estimate + 1.96 * se_value)
}

build_full_cohort_meta <- function(raw_meta, article_sample_names, data_level) {
  output <- raw_meta[!is.na(raw_meta$raw_group) & !is.na(raw_meta$subject_id), , drop = FALSE]
  output$analysis_sample <- output$raw_sample
  output$analysis_group <- output$raw_group
  output$group_source <- "Protein/peptide report label"
  correction_index <- match(output$raw_sample, names(SAMPLE_NAME_CORRECTIONS))
  correction_hit <- !is.na(correction_index)
  if (any(correction_hit)) {
    corrected_names <- unname(SAMPLE_NAME_CORRECTIONS[correction_index[correction_hit]])
    output$analysis_sample[correction_hit] <- corrected_names
    output$analysis_group[correction_hit] <- parse_group(corrected_names)
    output$group_source[correction_hit] <- "Article-corrected classification"
  }
  required_corrected_names <- unname(SAMPLE_NAME_CORRECTIONS)
  if (!all(required_corrected_names %in% article_sample_names)) {
    stop(data_level, ": corrected NFHP labels were not both found in the S1 article sample sheet.")
  }
  if (anyDuplicated(output$analysis_sample)) {
    duplicated_names <- unique(output$analysis_sample[duplicated(output$analysis_sample) |
      duplicated(output$analysis_sample, fromLast = TRUE)])
    stop(data_level, ": duplicated analysis sample names after correction: ", paste(duplicated_names,
      collapse = ", "))
  }
  observed_counts <- table(factor(output$analysis_group, levels = GROUP_LEVELS))
  if (nrow(output) != sum(FULL_EXPECTED_COUNTS) || !all(as.integer(observed_counts) == as.integer(FULL_EXPECTED_COUNTS))) {
    stop(data_level, ": full-cohort counts do not match HC=36, NFHP=31, FHP=28, IPF=30. ",
      "Observed: ", paste(paste0(names(FULL_EXPECTED_COUNTS), "=", as.integer(observed_counts)),
        collapse = ", "))
  }
  output
}

screen_m3_candidates <- function(m3_table, annotation_filtered, observed_log2, imputed_log2,
  group) {
  hc_idx <- which(group == "HC")
  ipf_idx <- which(group == "IPF")
  output_list <- list()
  for (gene_name in unique(m3_table$Gene)) {
    hit_rows <- which(vapply(annotation_filtered$Genes, contains_gene, logical(1), target_gene = gene_name))
    if (length(hit_rows) == 0) {
      next
    }
    detected_n <- rowSums(is.finite(observed_log2[hit_rows, , drop = FALSE]))
    selected_idx <- hit_rows[which.max(detected_n)]
    protein_row <- annotation_filtered$protein_row[selected_idx]
    matrix_row <- match(protein_row, rownames(observed_log2))
    if (is.na(matrix_row)) {
      next
    }
    hc_observed <- as.numeric(observed_log2[matrix_row, hc_idx, drop = TRUE])
    ipf_observed <- as.numeric(observed_log2[matrix_row, ipf_idx, drop = TRUE])
    hc_detection_rate <- mean(is.finite(hc_observed))
    ipf_detection_rate <- mean(is.finite(ipf_observed))
    hc_imputed <- as.numeric(imputed_log2[matrix_row, hc_idx, drop = TRUE])
    ipf_imputed <- as.numeric(imputed_log2[matrix_row, ipf_idx, drop = TRUE])
    test_result <- tryCatch(suppressWarnings(stats::wilcox.test(ipf_imputed, hc_imputed,
      alternative = "two.sided", exact = FALSE, correct = TRUE)), error = function(e) {
      NULL
    })
    p_value <- if (is.null(test_result)) {
      NA_real_
    } else {
      unname(test_result$p.value)
    }
    mean_log2fc <- mean(ipf_imputed) - mean(hc_imputed)
    m3_row <- match(gene_name, m3_table$Gene)
    mm_value <- m3_table$ModuleMembership[m3_row]
    gs_value <- m3_table$GeneSignificance[m3_row]
    pass_detection <- (hc_detection_rate >= CANDIDATE_DETECTION_CUTOFF && ipf_detection_rate >=
      CANDIDATE_DETECTION_CUTOFF)
    linear_fc <- 2^mean_log2fc
    pass_fc <- (is.finite(mean_log2fc) && abs(mean_log2fc) >= FC_CUTOFF)
    pass_p <- (is.finite(p_value) && p_value < FOCUS_WILCOXON_P_CUTOFF)
    detected_overlap <- (pass_detection)
    selected_candidate <- (detected_overlap && pass_fc && pass_p)
    selection_source <- if (selected_candidate) {
      "Detected M3 protein + IPF-vs-HC differential"
    } else if (detected_overlap) {
      "Detected M3 protein"
    } else {
      "Not selected"
    }
    output_list[[length(output_list) + 1]] <- data.frame(Gene = gene_name, protein_row = protein_row,
      ModuleMembership = mm_value, GeneSignificance = gs_value, DetectedOverlap = detected_overlap,
      HC_Detected = sum(is.finite(hc_observed)), HC_Total = length(hc_observed), HC_DetectionRate = hc_detection_rate,
      IPF_Detected = sum(is.finite(ipf_observed)), IPF_Total = length(ipf_observed),
      IPF_DetectionRate = ipf_detection_rate, mean_log2FC = mean_log2fc, linear_FC = linear_fc,
      P.Value = p_value, PassDetection = pass_detection, PassFC = pass_fc, PassP = pass_p,
      Selected = selected_candidate, SelectionSource = selection_source, Genes = annotation_filtered$Genes[selected_idx],
      ProteinAccessions = annotation_filtered$ProteinAccessions[selected_idx], ProteinNames = annotation_filtered$ProteinNames[selected_idx],
      FullCohortObservedCount = detected_n[which.max(detected_n)], stringsAsFactors = FALSE)
  }
  if (length(output_list) == 0) {
    stop("No M3 genes could be mapped to filtered protein groups.")
  }
  output <- do.call(rbind, output_list)
  output$adj.P.Val <- p.adjust(output$P.Value, method = "BH")
  output <- output[order(!output$DetectedOverlap, output$P.Value, -output$mean_log2FC), ,
    drop = FALSE]
  rownames(output) <- NULL
  output
}

format_gene_stats <- function(x) {
  if (is.null(x) || nrow(x) == 0) {
    return("none")
  }
  paste(paste0(x$Gene, " [FC=", sprintf("%.2f", x$linear_FC), ", log2FC=", sprintf("%.3f",
    x$mean_log2FC), ", P=", format.pval(x$P.Value, digits = 3, eps = 1e-04), "]"), collapse = "; ")
}

build_unified_selected_gene_order <- function(selected_table) {
  if (nrow(selected_table) == 0) {
    return(character(0))
  }
  up_table <- selected_table[selected_table$mean_log2FC > 0, , drop = FALSE]
  down_table <- selected_table[selected_table$mean_log2FC < 0, , drop = FALSE]
  up_table <- up_table[order(!up_table$Gene %in% FOREST_GENE_PRIORITY, match(up_table$Gene,
    FOREST_GENE_PRIORITY), up_table$P.Value, -abs(up_table$mean_log2FC)), , drop = FALSE]
  down_table <- down_table[order(down_table$P.Value, -abs(down_table$mean_log2FC)), , drop = FALSE]
  unique(c(up_table$Gene, down_table$Gene))
}

assert_selected_gene_set <- function(observed_genes, object_name) {
  observed_genes <- unique(as.character(observed_genes))
  observed_genes <- observed_genes[!is.na(observed_genes) & nzchar(observed_genes)]
  if (!setequal(observed_genes, SELECTED_GENES)) {
    stop(object_name, " does not contain exactly the automatically selected gene set. ",
      "Expected: ", paste(SELECTED_GENES, collapse = ", "), "; observed: ", paste(observed_genes,
        collapse = ", "))
  }
  invisible(TRUE)
}

run_protein_wilcoxon <- function(imputed_log2, group, annotation, control_group, case_group,
  comparison_name) {
  control_idx <- which(group == control_group)
  case_idx <- which(group == case_group)
  result_list <- vector("list", nrow(imputed_log2))
  for (i in seq_len(nrow(imputed_log2))) {
    control_values <- as.numeric(imputed_log2[i, control_idx, drop = TRUE])
    case_values <- as.numeric(imputed_log2[i, case_idx, drop = TRUE])
    test_result <- tryCatch(suppressWarnings(stats::wilcox.test(case_values, control_values,
      alternative = "two.sided", exact = FALSE, correct = TRUE)), error = function(e) {
      NULL
    })
    p_value <- if (is.null(test_result)) {
      NA_real_
    } else {
      unname(test_result$p.value)
    }
    mean_log2fc <- (mean(case_values) - mean(control_values))
    result_list[[i]] <- data.frame(protein_row = rownames(imputed_log2)[i], mean_log2FC = mean_log2fc,
      linear_FC = 2^mean_log2fc, median_difference = (median(case_values) - median(control_values)),
      P.Value = p_value, n_control = length(control_values), n_case = length(case_values),
      stringsAsFactors = FALSE)
  }
  result <- do.call(rbind, result_list)
  result$adj.P.Val <- p.adjust(result$P.Value, method = "BH")
  result <- merge(annotation, result, by = "protein_row", all.y = TRUE, sort = FALSE)
  result$Comparison <- comparison_name
  result$ControlGroup <- control_group
  result$CaseGroup <- case_group
  result$VolcanoClass <- "Not significant"
  result$VolcanoClass[result$P.Value < P_CUTOFF & result$mean_log2FC >= FC_CUTOFF] <- "Up"
  result$VolcanoClass[result$P.Value < P_CUTOFF & result$mean_log2FC <= -FC_CUTOFF] <- "Down"
  result
}

make_forest_data <- function(result, imputed_log2, group, control_group, case_group) {
  output_list <- list()
  for (i in seq_len(nrow(candidate_selection))) {
    protein_row <- candidate_selection$protein_row[i]
    result_row <- match(protein_row, result$protein_row)
    matrix_row <- match(protein_row, rownames(imputed_log2))
    if (is.na(result_row) || is.na(matrix_row)) {
      next
    }
    control_values <- as.numeric(imputed_log2[matrix_row, group == control_group, drop = TRUE])
    case_values <- as.numeric(imputed_log2[matrix_row, group == case_group, drop = TRUE])
    estimate <- (mean(case_values) - mean(control_values))
    ci_values <- bootstrap_se_ci(control_values, case_values, estimate, reps = BOOT_REPS)
    output_list[[length(output_list) + 1]] <- data.frame(Gene = candidate_selection$Gene[i],
      protein_row = protein_row, mean_log2FC = estimate, linear_FC = 2^estimate, CI_low = ci_values[1],
      CI_high = ci_values[2], P.Value = result$P.Value[result_row], adj.P.Val = result$adj.P.Val[result_row],
      stringsAsFactors = FALSE)
  }
  output <- do.call(rbind, output_list)
  output$Evidence <- ifelse(output$P.Value < P_CUTOFF, "Wilcoxon P < 0.05", "Wilcoxon P >= 0.05")
  output
}

run_peptide_statistics <- function(peptide_log2, peptide_group, control_group, case_group) {
  control_idx <- which(peptide_group == control_group)
  case_idx <- which(peptide_group == case_group)
  result_list <- list()
  for (gene_name in SELECTED_GENES) {
    accessions <- candidate_accession_map[[gene_name]]
    if (length(accessions) == 0) {
      next
    }
    hit_rows <- which(vapply(peptide_accession_tokens, function(z) {
      length(intersect(z, accessions)) > 0
    }, logical(1)))
    if (length(hit_rows) == 0) {
      next
    }
    for (row_idx in hit_rows) {
      control_values <- as.numeric(peptide_log2[row_idx, control_idx, drop = TRUE])
      case_values <- as.numeric(peptide_log2[row_idx, case_idx, drop = TRUE])
      control_values <- control_values[is.finite(control_values)]
      case_values <- case_values[is.finite(case_values)]
      if (length(control_values) < PEPTIDE_MIN_N_PER_GROUP || length(case_values) < PEPTIDE_MIN_N_PER_GROUP) {
        next
      }
      test_result <- tryCatch(suppressWarnings(stats::wilcox.test(case_values, control_values,
        alternative = "two.sided", exact = FALSE, correct = TRUE)), error = function(e) {
        NULL
      })
      if (is.null(test_result)) {
        next
      }
      estimate <- (mean(case_values) - mean(control_values))
      ci_values <- bootstrap_se_ci(control_values, case_values, estimate, reps = BOOT_REPS)
      result_list[[length(result_list) + 1]] <- data.frame(Gene = gene_name, Peptide = peptide_annotation$Peptide[row_idx],
        IsProteotypic = peptide_annotation$IsProteotypic[row_idx], IsProteinGroupSpecific = peptide_annotation$IsProteinGroupSpecific[row_idx],
        mean_log2FC = estimate, linear_FC = 2^estimate, CI_low = ci_values[1], CI_high = ci_values[2],
        P.Value = unname(test_result$p.value), n_control = length(control_values),
        n_case = length(case_values), stringsAsFactors = FALSE)
    }
  }
  if (length(result_list) == 0) {
    return(data.frame())
  }
  result <- do.call(rbind, result_list)
  result$adj.P.Val <- p.adjust(result$P.Value, method = "BH")
  result
}

make_box_data <- function(observed_log2, group, selected_groups) {
  box_list <- list()
  for (i in seq_len(nrow(candidate_selection))) {
    protein_row <- candidate_selection$protein_row[i]
    matrix_row <- match(protein_row, rownames(observed_log2))
    if (is.na(matrix_row)) {
      next
    }
    selected_index <- which(group %in% selected_groups)
    values <- as.numeric(observed_log2[matrix_row, selected_index, drop = TRUE])
    box_list[[length(box_list) + 1]] <- data.frame(Gene = candidate_selection$Gene[i],
      protein_row = protein_row, sample = colnames(observed_log2)[selected_index], group = as.character(group[selected_index]),
      abundance = values, detected = is.finite(values), stringsAsFactors = FALSE)
  }
  do.call(rbind, box_list)
}

make_detection_summary <- function(box_data) {
  detection_summary <- aggregate(detected ~ Gene + group, data = box_data, FUN = function(x) {
    c(detected_n = sum(x), total_n = length(x), detection_rate = mean(x))
  })
  data.frame(Gene = detection_summary$Gene, group = detection_summary$group, detected_n = detection_summary$detected[,
    "detected_n"], total_n = detection_summary$detected[, "total_n"], detection_rate = detection_summary$detected[,
    "detection_rate"], stringsAsFactors = FALSE)
}

make_qc_data <- function(observed_log2, imputed_log2, group, selected_groups) {
  selected_index <- which(group %in% selected_groups)
  selected_group <- factor(group[selected_index], levels = selected_groups)
  observed_subset <- observed_log2[, selected_index, drop = FALSE]
  imputed_subset <- imputed_log2[, selected_index, drop = FALSE]
  qc_table <- data.frame(sample = colnames(observed_subset), group = selected_group, detected_proteins = colSums(is.finite(observed_subset)),
    missing_rate = colMeans(!is.finite(observed_subset)), median_observed_log2 = apply(observed_subset,
      2, median, na.rm = TRUE), stringsAsFactors = FALSE)
  row_variance <- apply(imputed_subset, 1, stats::var, na.rm = TRUE)
  pca_input <- imputed_subset[is.finite(row_variance) & row_variance > 0, , drop = FALSE]
  pca_table <- data.frame()
  if (nrow(pca_input) >= 2 && ncol(pca_input) >= 3) {
    pca_fit <- stats::prcomp(t(pca_input), center = TRUE, scale. = TRUE)
    pca_variance <- summary(pca_fit)$importance[2, 1:2] * 100
    pca_table <- data.frame(sample = colnames(pca_input), group = selected_group, PC1 = pca_fit$x[,
      1], PC2 = pca_fit$x[, 2], PC1_variance = pca_variance[1], PC2_variance = pca_variance[2],
      stringsAsFactors = FALSE)
  }
  list(qc_table = qc_table, pca_table = pca_table)
}

run_comparison_data <- function(protein_pre, protein_group, protein_annotation_filtered, peptide_pre,
  peptide_group, control_group, case_group, comparison_name) {
  protein_result <- run_protein_wilcoxon(protein_pre$imputed_log2, protein_group, protein_annotation_filtered,
    control_group, case_group, comparison_name)
  forest_data <- make_forest_data(protein_result, protein_pre$imputed_log2, protein_group,
    control_group, case_group)
  peptide_result <- run_peptide_statistics(peptide_pre$observed_log2, peptide_group, control_group,
    case_group)
  box_data <- make_box_data(protein_pre$observed_log2, protein_group, c(control_group, case_group))
  detection_summary <- make_detection_summary(box_data)
  qc_data <- make_qc_data(protein_pre$observed_log2, protein_pre$imputed_log2, protein_group,
    c(control_group, case_group))
  list(comparison_name = comparison_name, control_group = control_group, case_group = case_group,
    protein_result = protein_result, forest_data = forest_data, peptide_result = peptide_result,
    box_data = box_data, detection_summary = detection_summary, qc_table = qc_data$qc_table,
    pca_table = qc_data$pca_table)
}

# 5. Read inputs and run all comparisons ----------------------------------------
for (path in c(PROTEIN_FILE, PEPTIDE_FILE, S1_FILE, M3_FILE)) {
  check_file(path)
}

s1_header <- readxl::read_excel(S1_FILE, sheet = "493proteins", skip = 1, n_max = 1, .name_repair = "minimal")

article_samples <- names(s1_header)[grepl("^(HC[0-9]+|NFHP_[0-9]+|FHP_[0-9]+|IPF_[0-9]+)$",
  names(s1_header))]

# Load the protein report, correct sample groups, then normalize the full cohort.
protein_dt <- data.table::fread(PROTEIN_FILE, sep = "\t", header = TRUE, check.names = FALSE,
  na.strings = c("", "NA", "NaN", "nan"))

required_protein_columns <- c("PG.ProteinAccessions", "PG.Genes", "PG.ProteinDescriptions",
  "PG.ProteinNames")

missing_protein_columns <- setdiff(required_protein_columns, names(protein_dt))

if (length(missing_protein_columns) > 0) {
  stop("Missing protein columns: ", paste(missing_protein_columns, collapse = ", "))
}

protein_quantity_columns <- grep("\\.PG\\.Quantity$", names(protein_dt), value = TRUE)

protein_raw_meta <- data.frame(raw_column = protein_quantity_columns, raw_sample = vapply(protein_quantity_columns,
  parse_quantity_label, character(1), level = "PG"), stringsAsFactors = FALSE)

protein_raw_meta$raw_group <- parse_group(protein_raw_meta$raw_sample)

protein_raw_meta$subject_id <- parse_subject_id(protein_raw_meta$raw_sample)

protein_full_meta <- build_full_cohort_meta(protein_raw_meta, article_samples, "Protein report")

protein_full_raw <- make_raw_matrix(protein_dt, protein_full_meta$raw_column, protein_full_meta$analysis_sample)

protein_full_group <- factor(protein_full_meta$analysis_group, levels = GROUP_LEVELS)

protein_annotation <- data.frame(protein_row = paste0("PG_", seq_len(nrow(protein_full_raw))),
  original_row = seq_len(nrow(protein_full_raw)), Genes = as.character(protein_dt$PG.Genes),
  ProteinAccessions = as.character(protein_dt$PG.ProteinAccessions), ProteinDescriptions = as.character(protein_dt$PG.ProteinDescriptions),
  ProteinNames = as.character(protein_dt$PG.ProteinNames), stringsAsFactors = FALSE)

protein_annotation$GeneLabel <- vapply(protein_annotation$Genes, first_token, character(1),
  fallback = "NoGene")

protein_annotation$AccessionLabel <- vapply(protein_annotation$ProteinAccessions, first_token,
  character(1), fallback = "NoAccession")

rownames(protein_full_raw) <- protein_annotation$protein_row

protein_full_pre <- preprocess_protein_cohort(protein_full_raw, protein_full_group, "Full-cohort protein")

protein_full_anno_filtered <- protein_annotation[protein_full_pre$keep, , drop = FALSE]

# Map the upstream M3 genes to protein groups and select candidates once.
m3 <- data.table::fread(M3_FILE)

required_m3_columns <- c("Gene", "ModuleMembership", "GeneSignificance", "Module")

missing_m3_columns <- setdiff(required_m3_columns, names(m3))

if (length(missing_m3_columns) > 0) {
  stop("Missing M3 columns: ", paste(missing_m3_columns, collapse = ", "))
}

m3[, `:=`(Gene, trimws(as.character(Gene)))]

m3[, `:=`(ModuleMembership, safe_numeric(ModuleMembership))]

m3[, `:=`(GeneSignificance, safe_numeric(GeneSignificance))]

m3[, `:=`(Final11, (Gene %in% FINAL_11))]

candidate_screening <- screen_m3_candidates(m3, protein_full_anno_filtered, protein_full_pre$observed_log2,
  protein_full_pre$imputed_log2, protein_full_group)

overlap_screening <- candidate_screening[candidate_screening$DetectedOverlap, , drop = FALSE]

if (nrow(overlap_screening) == 0) {
  stop("No M3 gene satisfied the detected-overlap detection criterion.")
}

selected_screening <- candidate_screening[candidate_screening$Selected, , drop = FALSE]

if (nrow(selected_screening) == 0) {
  stop("No detected M3-overlap gene met the focused protein-differential criterion. ",
    "Check input provenance and the study thresholds near the top of this script.")
}

OVERLAP_GENES <- overlap_screening$Gene

effect_size_screening <- overlap_screening[is.finite(overlap_screening$mean_log2FC) & abs(overlap_screening$mean_log2FC) >=
  FC_CUTOFF, , drop = FALSE]

effect_size_screening <- effect_size_screening[order(-abs(effect_size_screening$mean_log2FC),
  effect_size_screening$P.Value), , drop = FALSE]

up_effect_screening <- effect_size_screening[effect_size_screening$mean_log2FC >= FC_CUTOFF,
  , drop = FALSE]

down_effect_screening <- effect_size_screening[effect_size_screening$mean_log2FC <= -FC_CUTOFF,
  , drop = FALSE]

small_effect_significant_screening <- overlap_screening[is.finite(overlap_screening$P.Value) &
  overlap_screening$P.Value < FOCUS_WILCOXON_P_CUTOFF & is.finite(overlap_screening$mean_log2FC) &
  abs(overlap_screening$mean_log2FC) < FC_CUTOFF, , drop = FALSE]

selected_up_screening <- selected_screening[selected_screening$mean_log2FC > 0, , drop = FALSE]

selected_down_screening <- selected_screening[selected_screening$mean_log2FC < 0, , drop = FALSE]

SELECTED_GENES <- build_unified_selected_gene_order(selected_screening)

selected_screening$PlotOrder <- match(selected_screening$Gene, SELECTED_GENES)

selected_screening <- selected_screening[order(selected_screening$PlotOrder), , drop = FALSE]

overlap_selection <- overlap_screening[, c("Gene", "protein_row", "Genes", "ProteinAccessions",
  "ProteinNames", "FullCohortObservedCount"), drop = FALSE]

candidate_selection <- selected_screening[, c("Gene", "protein_row", "Genes", "ProteinAccessions",
  "ProteinNames", "FullCohortObservedCount"), drop = FALSE]

candidate_accession_map <- lapply(seq_len(nrow(candidate_selection)), function(i) {
  split_tokens(candidate_selection$ProteinAccessions[i])
})

names(candidate_accession_map) <- candidate_selection$Gene

message("M3 proteins detected in BOTH HC and IPF at >= ", CANDIDATE_DETECTION_CUTOFF * 100,
  "%: ", nrow(overlap_selection), " genes")

message("FC >= ", FOCUS_LINEAR_FC_CUTOFF, " (up; regardless of P): ", format_gene_stats(up_effect_screening))

message("FC <= ", round(1/FOCUS_LINEAR_FC_CUTOFF, 3), " (down; regardless of P): ", format_gene_stats(down_effect_screening))

message("Primary UP candidates (effect cutoff + Wilcoxon P < ", FOCUS_WILCOXON_P_CUTOFF, "): ",
  format_gene_stats(selected_up_screening))

message("Primary DOWN candidates (effect cutoff + Wilcoxon P < ", FOCUS_WILCOXON_P_CUTOFF,
  "): ", format_gene_stats(selected_down_screening))

message("P < ", FOCUS_WILCOXON_P_CUTOFF, " but |FC| < 1.20-fold (exploratory small-effect set): ",
  format_gene_stats(small_effect_significant_screening))

detected_gene_set <- unique(unlist(lapply(protein_full_anno_filtered$Genes, split_tokens)))

m3_plot_data <- as.data.frame(m3)

m3_plot_data$Category <- "Other M3 gene"

m3_plot_data$Category[m3_plot_data$Gene %in% OVERLAP_GENES] <- "M3 protein detected in both groups (>=60%)"

m3_plot_data$Category[m3_plot_data$Gene %in% SELECTED_GENES] <- "IPF-HC differential candidate"

# Load peptide quantities and specificity annotations using the same cohort.
peptide_dt <- data.table::fread(PEPTIDE_FILE, sep = "\t", header = TRUE, check.names = FALSE,
  na.strings = c("", "NA", "NaN", "nan"))

required_peptide_columns <- c("PG.ProteinAccessions", "PEP.StrippedSequence", "PEP.IsProteotypic",
  "PEP.IsProteinGroupSpecific")

missing_peptide_columns <- setdiff(required_peptide_columns, names(peptide_dt))

if (length(missing_peptide_columns) > 0) {
  stop("Missing peptide columns: ", paste(missing_peptide_columns, collapse = ", "))
}

peptide_quantity_columns <- grep("\\.PEP\\.Quantity$", names(peptide_dt), value = TRUE)

peptide_raw_meta <- data.frame(raw_column = peptide_quantity_columns, raw_sample = vapply(peptide_quantity_columns,
  parse_quantity_label, character(1), level = "PEP"), stringsAsFactors = FALSE)

peptide_raw_meta$raw_group <- parse_group(peptide_raw_meta$raw_sample)

peptide_raw_meta$subject_id <- parse_subject_id(peptide_raw_meta$raw_sample)

peptide_full_meta <- build_full_cohort_meta(peptide_raw_meta, article_samples, "Peptide report")

peptide_full_raw <- make_raw_matrix(peptide_dt, peptide_full_meta$raw_column, peptide_full_meta$analysis_sample)

peptide_full_group <- factor(peptide_full_meta$analysis_group, levels = GROUP_LEVELS)

peptide_full_pre <- preprocess_peptide_cohort(peptide_full_raw, "Full-cohort peptide")

peptide_annotation <- data.frame(peptide_row = seq_len(nrow(peptide_dt)), ProteinAccessions = as.character(peptide_dt$PG.ProteinAccessions),
  Peptide = as.character(peptide_dt$PEP.StrippedSequence), IsProteotypic = safe_bool(peptide_dt$PEP.IsProteotypic),
  IsProteinGroupSpecific = safe_bool(peptide_dt$PEP.IsProteinGroupSpecific), stringsAsFactors = FALSE)

peptide_accession_tokens <- lapply(peptide_annotation$ProteinAccessions, split_tokens)

comparisons <- list(IPF_vs_HC = run_comparison_data(protein_full_pre, protein_full_group, protein_full_anno_filtered,
  peptide_full_pre, peptide_full_group, "HC", "IPF", "IPF vs HC"), FHP_vs_HC = run_comparison_data(protein_full_pre,
  protein_full_group, protein_full_anno_filtered, peptide_full_pre, peptide_full_group, "HC",
  "FHP", "FHP vs HC"), NFHP_vs_HC = run_comparison_data(protein_full_pre, protein_full_group,
  protein_full_anno_filtered, peptide_full_pre, peptide_full_group, "HC", "NFHP", "NFHP vs HC"))

four_group_box_data <- make_box_data(protein_full_pre$observed_log2, protein_full_group, GROUP_LEVELS)

four_group_detection_summary <- make_detection_summary(four_group_box_data)

four_group_qc <- make_qc_data(protein_full_pre$observed_log2, protein_full_pre$imputed_log2,
  protein_full_group, GROUP_LEVELS)

assert_selected_gene_set(candidate_selection$Gene, "candidate_selection")

assert_selected_gene_set(comparisons$IPF_vs_HC$forest_data$Gene, "IPF_vs_HC forest data")

assert_selected_gene_set(comparisons$FHP_vs_HC$forest_data$Gene, "FHP_vs_HC forest data")

assert_selected_gene_set(comparisons$NFHP_vs_HC$forest_data$Gene, "NFHP_vs_HC forest data")

assert_selected_gene_set(four_group_box_data$Gene, "four-group boxplot data")

assert_selected_gene_set(four_group_detection_summary$Gene, "four-group detection data")

message("Selected genes successfully propagated to all candidate-level downstream data objects.")

# 6. Plot functions -------------------------------------------------------------
save_pdf <- function(plot_object, filename, width, height) {
  grDevices::pdf(filename, width = width, height = height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot_object)
  invisible(filename)
}

theme_nice <- function(base_size = 11) {
  theme_bw(base_size = base_size) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    plot.title = element_blank(), plot.subtitle = element_blank(), axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"), legend.title = element_blank(), strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = "black", margin = margin(b = 5)),
    legend.position = "bottom", plot.margin = margin(10, 18, 10, 10))
}

format_p <- function(p) {
  if (!is.finite(p)) {
    return("P=NA")
  }
  if (p < 1e-04) {
    return("P<1e-4")
  }
  if (p < 0.001) {
    return(paste0("P=", formatC(p, format = "e", digits = 2)))
  }
  paste0("P=", sprintf("%.3f", p))
}

get_public_label <- function(comparison_name) {
  x <- comparison_name
  x <- gsub(" - AllReport", "", x, fixed = TRUE)
  x <- gsub(" - Article", "", x, fixed = TRUE)
  x
}

make_gs_mm_plot <- function(m3_plot_data, output_directory) {
  m3_plot_data$Category <- factor(m3_plot_data$Category, levels = c("Other M3 gene", "M3 protein detected in both groups (>=60%)",
    "IPF-HC differential candidate"))
  label_data <- m3_plot_data[m3_plot_data$Gene %in% SELECTED_GENES, , drop = FALSE]
  plot_object <- ggplot(m3_plot_data, aes(x = ModuleMembership, y = GeneSignificance, color = Category,
    shape = Category)) + geom_point(alpha = 0.78, size = 1.9, stroke = 0.9) + ggrepel::geom_text_repel(data = label_data,
    aes(label = Gene), color = "black", fontface = "bold", size = 3.35, max.overlaps = Inf,
    box.padding = 0.42, point.padding = 0.25, min.segment.length = 0, seed = 8848, show.legend = FALSE) +
    scale_color_manual(values = c(`Other M3 gene` = "#BDBDBD", `M3 protein detected in both groups (>=60%)` = "#2A8C82",
      `IPF-HC differential candidate` = "black"), name = NULL) + scale_shape_manual(values = c(`Other M3 gene` = 16,
    `M3 protein detected in both groups (>=60%)` = 17, `IPF-HC differential candidate` = 16),
    name = NULL) + guides(color = guide_legend(order = 1, override.aes = list(alpha = 1,
    size = 2.8)), shape = guide_legend(order = 1, override.aes = list(alpha = 1, size = 2.8))) +
    labs(x = "Module Membership (MM)", y = "Gene Significance (GS)") + theme_nice()
  save_pdf(plot_object, file.path(output_directory, "M3_GS_MM_candidate_scatter.pdf"), 7.8,
    6)
}

make_volcano_plot <- function(result, output_directory) {
  volcano_data <- result[result$protein_row %in% overlap_selection$protein_row & is.finite(result$mean_log2FC) &
    is.finite(result$P.Value), , drop = FALSE]
  volcano_data$minus_log10_P <- -log10(pmax(volcano_data$P.Value, 1e-300))
  volcano_data$VolcanoClass <- factor(volcano_data$VolcanoClass, levels = c("Down", "Not significant",
    "Up"))
  label_rows <- volcano_data[volcano_data$protein_row %in% candidate_selection$protein_row,
    , drop = FALSE]
  label_rows$DisplayLabel <- candidate_selection$Gene[match(label_rows$protein_row, candidate_selection$protein_row)]
  label_rows <- label_rows[order(label_rows$mean_log2FC, label_rows$P.Value), , drop = FALSE]
  plot_object <- ggplot(volcano_data, aes(x = mean_log2FC, y = minus_log10_P, color = VolcanoClass)) +
    geom_point(alpha = 0.68, size = 1.6) + geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF),
    linetype = "dashed", linewidth = 0.45) + geom_hline(yintercept = -log10(P_CUTOFF),
    linetype = "dashed", linewidth = 0.45) + ggrepel::geom_text_repel(data = label_rows,
    aes(label = DisplayLabel), color = "black", fontface = "bold", size = 3.4, max.overlaps = Inf,
    box.padding = 0.45, point.padding = 0.25, min.segment.length = 0, seed = 8848, show.legend = FALSE) +
    scale_color_manual(values = c(Down = "#4C78A8", `Not significant` = "#CFCFCF", Up = "#E45756")) +
    labs(x = "Mean log2 fold change", y = "-log10(Wilcoxon P value)") + theme_nice()
  save_pdf(plot_object, file.path(output_directory, "Volcano_Wilcoxon_MinImpute.pdf"), 7.6,
    6.1)
}

make_forest_plot <- function(forest_data, output_directory) {
  forest_data$Gene <- factor(forest_data$Gene, levels = rev(SELECTED_GENES))
  forest_data$ResultLabel <- paste0("FC=", sprintf("%.2f", forest_data$linear_FC), "; ",
    vapply(forest_data$P.Value, format_p, character(1)))
  effect_range <- diff(range(c(forest_data$CI_low, forest_data$CI_high), na.rm = TRUE))
  if (!is.finite(effect_range) || effect_range <= 0) {
    effect_range <- 1
  }
  label_gap <- max(0.06, 0.025 * effect_range)
  forest_data$Label_x <- (forest_data$CI_high + label_gap)
  forest_x_min <- (min(forest_data$CI_low, na.rm = TRUE) - 0.08 * effect_range)
  forest_x_max <- (max(forest_data$Label_x, na.rm = TRUE) + max(0.55, 0.28 * effect_range))
  plot_object <- ggplot(forest_data, aes(x = mean_log2FC, y = Gene, color = Evidence)) +
    geom_vline(xintercept = 0, linewidth = 0.55, color = "black") + geom_vline(xintercept = c(-FC_CUTOFF,
    FC_CUTOFF), linetype = "dotted", linewidth = 0.45, color = "#555555") + geom_errorbar(aes(xmin = CI_low,
    xmax = CI_high), orientation = "y", width = 0.16, linewidth = 0.68) + geom_point(size = 2.9) +
    geom_text(aes(x = Label_x, label = ResultLabel), color = "black", hjust = 0, size = 3,
      fontface = "bold", show.legend = FALSE) + scale_color_manual(values = c(`Wilcoxon P < 0.05` = "#B2182B",
    `Wilcoxon P >= 0.05` = "#7F7F7F")) + coord_cartesian(xlim = c(forest_x_min, forest_x_max),
    clip = "off") + labs(x = "Mean log2 fold change with 95% CI", y = NULL) + theme_nice() +
    theme(plot.margin = margin(8, 34, 8, 8))
  save_pdf(plot_object, file.path(output_directory, "Forest_Wilcoxon_MinImpute.pdf"), 9.1,
    5.6)
}

make_boxplot <- function(box_data, selected_groups, output_directory, four_group = FALSE) {
  plot_data <- box_data[box_data$detected & is.finite(box_data$abundance), , drop = FALSE]
  plot_data$Gene <- factor(plot_data$Gene, levels = SELECTED_GENES)
  plot_data$group <- factor(plot_data$group, levels = selected_groups)
  if (four_group) {
    plot_object <- ggplot(plot_data, aes(x = group, y = abundance, fill = group)) + geom_boxplot(width = 0.6,
      outlier.shape = NA, alpha = 0.72, color = "black") + geom_jitter(width = 0.1, size = 1.05,
      alpha = 0.74) + facet_wrap(~Gene, scales = "free_y", ncol = 3) + scale_fill_manual(values = GROUP_COLORS) +
      labs(x = NULL, y = "Median-normalized log2 abundance") + theme_nice() + theme(legend.position = "none",
      axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"))
    save_pdf(plot_object, file.path(output_directory, "Boxplot_FourGroup_ObservedOnly_NoPvalue.pdf"),
      11.6, 7.6)
  } else {
    plot_object <- ggplot(plot_data, aes(x = group, y = abundance, fill = group)) + geom_boxplot(width = 0.58,
      outlier.shape = NA, alpha = 0.72, color = "black") + geom_jitter(width = 0.1, size = 1.25,
      alpha = 0.78) + facet_wrap(~Gene, scales = "free_y", nrow = 1) + scale_fill_manual(values = GROUP_COLORS) +
      labs(x = NULL, y = "Median-normalized log2 abundance") + theme_nice() + theme(legend.position = "none",
      axis.text.x = element_text(face = "bold"))
    save_pdf(plot_object, file.path(output_directory, "Boxplot_ObservedOnly_NoPvalue.pdf"),
      12.8, 4.9)
  }
}

make_detection_plot <- function(detection_summary, selected_groups, output_directory, four_group = FALSE) {
  detection_summary$Gene <- factor(detection_summary$Gene, levels = SELECTED_GENES)
  detection_summary$group <- factor(detection_summary$group, levels = selected_groups)
  detection_summary$DetectionLabel <- paste0(detection_summary$detected_n, "/", detection_summary$total_n)
  detection_summary$Label_y <- pmax(detection_summary$detection_rate - 0.015, 0.08)
  dodge_width <- ifelse(four_group, 0.84, 0.76)
  bar_width <- ifelse(four_group, 0.7, 0.66)
  label_size <- ifelse(four_group, 2.15, 2.45)
  plot_object <- ggplot(detection_summary, aes(x = Gene, y = detection_rate, fill = group)) +
    geom_col(position = position_dodge(width = dodge_width), width = bar_width, color = "black",
      linewidth = 0.22) + geom_text(aes(y = Label_y, label = DetectionLabel), position = position_dodge(width = dodge_width),
    angle = 90, hjust = 1, vjust = 0.5, size = label_size, color = "black", show.legend = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.04),
      breaks = seq(0, 1, by = 0.2), expand = expansion(mult = c(0, 0.01))) + scale_fill_manual(values = GROUP_COLORS) +
    labs(x = NULL, y = "Protein detection rate") + theme_nice() + theme(axis.text.x = element_text(angle = 25,
    hjust = 1, face = "bold"))
  save_pdf(plot_object, file.path(output_directory, ifelse(four_group, "ProteinDetectionRate_FourGroup.pdf",
    "ProteinDetectionRate.pdf")), ifelse(four_group, 12.8, 9.2), ifelse(four_group, 5.8,
    5.2))
}

make_peptide_plot <- function(peptide_result, output_directory) {
  if (nrow(peptide_result) == 0) {
    return(invisible(NULL))
  }
  plot_data <- peptide_result[peptide_result$IsProteotypic & peptide_result$IsProteinGroupSpecific &
    is.finite(peptide_result$mean_log2FC) & is.finite(peptide_result$CI_low) & is.finite(peptide_result$CI_high),
    , drop = FALSE]
  if (nrow(plot_data) == 0) {
    return(invisible(NULL))
  }
  plot_data <- do.call(rbind, lapply(split(plot_data, plot_data$Gene), function(x) {
    x <- x[order(x$P.Value), , drop = FALSE]
    head(x, PEPTIDE_TOP_PER_GENE)
  }))
  peptide_gene_order_top_to_bottom <- SELECTED_GENES[SELECTED_GENES %in% unique(as.character(plot_data$Gene))]
  plot_data$GeneOrder <- match(as.character(plot_data$Gene), peptide_gene_order_top_to_bottom)
  plot_data <- plot_data[order(plot_data$GeneOrder, plot_data$P.Value, -abs(plot_data$mean_log2FC)),
    , drop = FALSE]
  plot_data$Gene <- factor(plot_data$Gene, levels = peptide_gene_order_top_to_bottom)
  plot_data$Evidence <- ifelse(plot_data$P.Value < P_CUTOFF, "Wilcoxon P < 0.05", "Wilcoxon P >= 0.05")
  plot_data$Peptide_label <- paste0(plot_data$Gene, " | ", plot_data$Peptide)
  peptide_label_levels_bottom_to_top <- rev(unique(plot_data$Peptide_label))
  plot_data$Peptide_label <- factor(plot_data$Peptide_label, levels = peptide_label_levels_bottom_to_top)
  peptide_range <- diff(range(c(plot_data$CI_low, plot_data$CI_high), na.rm = TRUE))
  if (!is.finite(peptide_range) || peptide_range <= 0) {
    peptide_range <- 1
  }
  peptide_x_min <- (min(plot_data$CI_low, na.rm = TRUE) - 0.06 * peptide_range)
  peptide_x_max <- (max(plot_data$CI_high, na.rm = TRUE) + 0.06 * peptide_range)
  plot_object <- ggplot(plot_data, aes(x = mean_log2FC, y = Peptide_label, color = Evidence)) +
    geom_vline(xintercept = 0, linewidth = 0.5) + geom_errorbar(aes(xmin = CI_low, xmax = CI_high),
    orientation = "y", width = 0.18, linewidth = 0.65) + geom_point(size = 2.4) + scale_color_manual(values = c(`Wilcoxon P < 0.05` = "#B2182B",
    `Wilcoxon P >= 0.05` = "#7F7F7F")) + coord_cartesian(xlim = c(peptide_x_min, peptide_x_max),
    clip = "off") + labs(x = "Peptide mean log2 fold change with 95% CI", y = NULL) + theme_nice(base_size = 10)
  save_pdf(plot_object, file.path(output_directory, "PeptideDirection_Observed_Wilcoxon.pdf"),
    9.2, max(5.6, 0.3 * nrow(plot_data) + 2.6))
}

make_qc_plots <- function(qc_table, pca_table, output_directory) {
  qc_table$group <- factor(qc_table$group, levels = unique(as.character(qc_table$group)))
  detected_plot <- ggplot(qc_table, aes(x = group, y = detected_proteins, fill = group)) +
    geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.7, color = "black") + geom_jitter(width = 0.1,
    size = 1.45, alpha = 0.76) + scale_fill_manual(values = GROUP_COLORS) + labs(x = NULL,
    y = "Detected protein groups") + theme_nice() + theme(legend.position = "none", axis.text.x = element_text(face = "bold"))
  save_pdf(detected_plot, file.path(output_directory, "QC_DetectedProteinGroups.pdf"), 5.4,
    4.8)
  missing_plot <- ggplot(qc_table, aes(x = group, y = missing_rate, fill = group)) + geom_boxplot(width = 0.58,
    outlier.shape = NA, alpha = 0.7, color = "black") + geom_jitter(width = 0.1, size = 1.45,
    alpha = 0.76) + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = GROUP_COLORS) + labs(x = NULL, y = "Protein missing rate") +
    theme_nice() + theme(legend.position = "none", axis.text.x = element_text(face = "bold"))
  save_pdf(missing_plot, file.path(output_directory, "QC_MissingRate.pdf"), 5.4, 4.8)
  intensity_plot <- ggplot(qc_table, aes(x = group, y = median_observed_log2, fill = group)) +
    geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.7, color = "black") + geom_jitter(width = 0.1,
    size = 1.45, alpha = 0.76) + scale_fill_manual(values = GROUP_COLORS) + labs(x = NULL,
    y = "Median observed log2 abundance") + theme_nice() + theme(legend.position = "none",
    axis.text.x = element_text(face = "bold"))
  save_pdf(intensity_plot, file.path(output_directory, "QC_MedianObservedIntensity.pdf"),
    5.4, 4.8)
  if (nrow(pca_table) > 0) {
    pca_table$group <- factor(pca_table$group, levels = unique(as.character(pca_table$group)))
    pc1_variance <- unique(pca_table$PC1_variance)[1]
    pc2_variance <- unique(pca_table$PC2_variance)[1]
    pca_plot <- ggplot(pca_table, aes(x = PC1, y = PC2, color = group)) + geom_point(size = 2.25,
      alpha = 0.82) + stat_ellipse(type = "norm", linewidth = 0.6, show.legend = FALSE) +
      scale_color_manual(values = GROUP_COLORS) + labs(x = paste0("PC1 (", round(pc1_variance,
      1), "%)"), y = paste0("PC2 (", round(pc2_variance, 1), "%)")) + theme_nice()
    save_pdf(pca_plot, file.path(output_directory, "QC_PCA.pdf"), 6.4, 5.2)
  }
}

make_combined_disease_forest <- function(comparisons, output_directory) {
  combined_result <- rbind(transform(comparisons$IPF_vs_HC$forest_data, Comparison = "IPF vs HC"),
    transform(comparisons$FHP_vs_HC$forest_data, Comparison = "FHP vs HC"), transform(comparisons$NFHP_vs_HC$forest_data,
      Comparison = "NFHP vs HC"))
  display_order_top_to_bottom <- SELECTED_GENES
  combined_result$Gene <- factor(combined_result$Gene, levels = display_order_top_to_bottom)
  combined_result$Comparison <- factor(combined_result$Comparison, levels = c("IPF vs HC",
    "FHP vs HC", "NFHP vs HC"))
  combined_result$ResultLabel <- paste0("FC=", sprintf("%.2f", combined_result$linear_FC),
    "; ", vapply(combined_result$P.Value, format_p, character(1)))
  base_y <- rev(seq_along(display_order_top_to_bottom))
  names(base_y) <- display_order_top_to_bottom
  comparison_offsets <- c(`IPF vs HC` = 0.24, `FHP vs HC` = 0, `NFHP vs HC` = -0.24)
  combined_result$y <- (base_y[as.character(combined_result$Gene)] + comparison_offsets[as.character(combined_result$Comparison)])
  effect_range <- diff(range(c(combined_result$CI_low, combined_result$CI_high), na.rm = TRUE))
  if (!is.finite(effect_range) || effect_range <= 0) {
    effect_range <- 1
  }
  label_gap <- max(0.06, 0.02 * effect_range)
  combined_result$Label_x <- (combined_result$CI_high + label_gap)
  forest_x_min <- (min(combined_result$CI_low, na.rm = TRUE) - 0.08 * effect_range)
  forest_x_max <- (max(combined_result$Label_x, na.rm = TRUE) + max(0.6, 0.3 * effect_range))
  plot_object <- ggplot(combined_result, aes(x = mean_log2FC, y = y, color = Comparison,
    shape = Comparison)) + geom_vline(xintercept = 0, linewidth = 0.55, color = "black") +
    geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF), linetype = "dotted", linewidth = 0.45,
      color = "#555555") + geom_segment(aes(x = CI_low, xend = CI_high, yend = y), linewidth = 0.78) +
    geom_point(size = 2.8) + geom_text(aes(x = Label_x, label = ResultLabel), color = "black",
    hjust = 0, size = 2.55, fontface = "bold", show.legend = FALSE) + scale_y_continuous(breaks = base_y,
    labels = names(base_y), limits = c(0.55, length(display_order_top_to_bottom) + 0.45),
    expand = expansion(mult = c(0, 0))) + scale_color_manual(values = c(`IPF vs HC` = "#D95F59",
    `FHP vs HC` = "#C97B1D", `NFHP vs HC` = "#3A924A")) + scale_shape_manual(values = c(`IPF vs HC` = 15,
    `FHP vs HC` = 16, `NFHP vs HC` = 17)) + coord_cartesian(xlim = c(forest_x_min, forest_x_max),
    clip = "off") + labs(x = "Mean log2 fold change with 95% CI", y = NULL) + theme_nice() +
    theme(plot.margin = margin(8, 34, 8, 8), axis.text.y = element_text(face = "bold"))
  save_pdf(plot_object, file.path(output_directory, "Forest_IPF_FHP_NFHP_vs_HC_Combined_FullCohort_Wilcoxon_MinImpute.pdf"),
    10.2, 6.6)
  invisible(combined_result)
}

# 7. Export figures -------------------------------------------------------------
# Preserve the original plotting order and random jitter sequence.
make_gs_mm_plot(m3_plot_data, DIR_GS_MM)

comparison_output_map <- list(IPF_vs_HC = DIR_IPF, FHP_vs_HC = DIR_FHP, NFHP_vs_HC = DIR_NFHP)

for (comparison_id in names(comparison_output_map)) {
  comparison_object <- comparisons[[comparison_id]]
  output_directory <- comparison_output_map[[comparison_id]]
  selected_groups <- c(comparison_object$control_group, comparison_object$case_group)
  make_volcano_plot(comparison_object$protein_result, output_directory)
  make_forest_plot(comparison_object$forest_data, output_directory)
  make_boxplot(comparison_object$box_data, selected_groups, output_directory, four_group = FALSE)
  make_detection_plot(comparison_object$detection_summary, selected_groups, output_directory,
    four_group = FALSE)
  make_peptide_plot(comparison_object$peptide_result, output_directory)
  make_qc_plots(comparison_object$qc_table, comparison_object$pca_table, output_directory)
}

make_boxplot(four_group_box_data, GROUP_LEVELS, DIR_FOUR_GROUP, four_group = TRUE)

make_detection_plot(four_group_detection_summary, GROUP_LEVELS, DIR_FOUR_GROUP, four_group = TRUE)

make_qc_plots(four_group_qc$qc_table, four_group_qc$pca_table, DIR_QC)

combined_forest_data <- make_combined_disease_forest(comparisons, DIR_COMBINED_FOREST)

assert_selected_gene_set(combined_forest_data$Gene, "combined three-disease forest data")

