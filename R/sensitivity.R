# Reported component, unit, cell-signal, and donor-influence analyses in R.
# Shared validators and bulk helpers are loaded by run.R without rerunning them.

sensitivity_development_components <- function(development) {
  panel <- development_gene_panel()
  panel <- panel[seq_len(7L), , drop = FALSE]
  pooled <- development$pooled
  rows <- vapply(panel$gene_id, function(identifier) {
    hits <- which(pooled$stable_gene_id == identifier)
    if (length(hits) != 1L) stop("A development score component lacks a unique source row.")
    hits[[1L]]
  }, integer(1))
  metadata <- development$analysis
  ids <- as.character(metadata$participant_id)
  columns <- match(ids, colnames(pooled$counts))
  normalization <- match(ids, as.character(pooled$tmm$participant_id))
  if (anyNA(columns) || anyNA(normalization) || anyDuplicated(ids)) stop("Misaligned development component identities.")
  block <- pooled$counts[rows, columns, drop = FALSE]
  dimnames(block) <- list(panel$gene_symbol, ids)
  tmm <- pooled$tmm[normalization, , drop = FALSE]
  if (any(!is.finite(block)) || any(block < 0) ||
      any(tmm$effective_library_size <= 0) || any(!is.finite(tmm$effective_library_size))) {
    stop("Invalid development components or effective libraries.")
  }
  numerator <- colSums(block[seq_len(4L), , drop = FALSE])
  denominator <- colSums(block[5L:7L, , drop = FALSE])
  score <- log2((numerator + 0.5) / (denominator + 0.5))
  if (max(abs(score - metadata$IGG_IGAM_LOG_RATIO)) > 1e-10) stop("Development components do not reconstruct the validated score.")
  if (!"diagnosis_harmonized" %in% names(metadata)) stop("Missing validated development diagnosis labels.")
  out <- data.frame(participant_id = metadata$participant_id, diagnosis = as.character(metadata$diagnosis_harmonized),
    n_biopsies = metadata$n_rectum_biopsies, raw_library_size = tmm$raw_library_size,
    tmm_norm_factor = tmm$tmm_norm_factor, effective_library_size = tmm$effective_library_size,
    stringsAsFactors = FALSE)
  if (nrow(out) != 665L || anyNA(out$diagnosis) || any(!out$diagnosis %in% c("CD", "UC"))) {
    stop("Invalid development component population.")
  }
  count_columns <- paste0(panel$gene_symbol, "__", panel$gene_id, "__count")
  for (i in seq_along(count_columns)) out[[count_columns[[i]]]] <- as.numeric(block[i, ])
  out$score_rebuilt <- as.numeric(score)
  out$score_stage29 <- metadata$IGG_IGAM_LOG_RATIO
  out$difference <- out$score_rebuilt - out$score_stage29
  list(table = out, panel = block, G = as.numeric(numerator), D = as.numeric(denominator),
    diagnosis = out$diagnosis, is_uc = metadata$is_uc, y_cd = metadata$y_cd,
    raw_library_size = tmm$raw_library_size, effective_library_size = tmm$effective_library_size,
    score = out$score_stage29)
}

sensitivity_component_rows <- function(cohort, role, ids, diagnosis, group, native_unit, normalized_unit, G, D, k, saved) {
  if (length(k) == 1L) k <- rep(k, length(ids))
  if (any(!is.finite(c(G, D, k, saved))) || any(G <= 0) || any(D <= 0) || any(k <= 0) ||
      anyDuplicated(ids) || length(G) != length(ids) || length(D) != length(ids) || length(k) != length(ids)) {
    stop("Component descriptions require aligned positive native sums and scale factors.")
  }
  diagnosis <- as.character(diagnosis)
  group <- as.character(group)
  diagnosis[diagnosis == "CONTROL"] <- "Control"
  group[group == "CONTROL"] <- "Control"
  out <- data.frame(cohort = cohort, role = role, participant_id = as.character(ids),
    diagnosis = diagnosis, descriptive_group = group, native_unit = native_unit, normalized_unit = normalized_unit,
    G = G, D = D, scale_k = k, G_normalized = k * G, D_normalized = k * D, scaled_pseudocount = k * 0.5,
    logG_native = log2(G + 0.5), logD_native = log2(D + 0.5),
    logG_normalized_exact = log2(k * (G + 0.5)), logD_normalized_exact = log2(k * (D + 0.5)),
    logG_normalized_no_pseudo = log2(k * G), logD_normalized_no_pseudo = log2(k * D),
    score = log2((G + 0.5) / (D + 0.5)), score_saved = saved, stringsAsFactors = FALSE)
  out$reconstruction_error <- out$score - out$score_saved
  out$identity_error <- out$score - (out$logG_normalized_exact - out$logD_normalized_exact)
  if (max(abs(out$reconstruction_error), abs(out$identity_error)) > 1e-10) stop("A native or normalized component identity failed.")
  out
}

sensitivity_component_analysis <- function(development, external) {
  ext <- external$external$all
  context <- external$context$all
  group <- ext$diagnosis
  group[ext$group == "ILEAL_CD"] <- "CD L3 (source: ileal CD)"
  group[ext$group == "COLONIC_CD"] <- "CD L2 (source: colonic CD)"
  rows <- rbind(
    sensitivity_component_rows("GSE193677", "development", development$table$participant_id,
      development$diagnosis, development$diagnosis, "estimated counts", "TMM-adjusted CPM",
      development$G, development$D, 1e6 / development$effective_library_size, development$score),
    sensitivity_component_rows("GSE117993", "external", ext$sample_id, ext$diagnosis, group, "TPM", "TPM",
      ext$igg_sum, ext$igam_sum, 1, ext$IGG_IGAM_LOG_RATIO),
    sensitivity_component_rows("GSE111889", "context only", context$participant_id, context$diagnosis,
      context$diagnosis, "pooled counts", "pooled CPM", context$igg_sum, context$igam_sum,
      1e6 / context$raw_library_size, context$IGG_IGAM_LOG_RATIO))
  if (nrow(rows) != 945L || anyDuplicated(paste(rows$cohort, rows$participant_id, sep = "||"))) {
    stop("Unexpected combined component population.")
  }
  metrics <- c("G", "D", "G_normalized", "D_normalized", "logG_native", "logD_native",
    "logG_normalized_exact", "logD_normalized_exact", "logG_normalized_no_pseudo", "logD_normalized_no_pseudo", "score")
  descriptions <- list()
  for (grouping in c("diagnosis", "descriptive_group")) {
    for (cohort in unique(rows$cohort)) {
      data <- rows[rows$cohort == cohort, , drop = FALSE]
      for (label in unique(data[[grouping]])) {
        selected <- data[data[[grouping]] == label, , drop = FALSE]
        for (metric in metrics) {
          values <- selected[[metric]]
          quartiles <- as.numeric(stats::quantile(values, c(0.25, 0.75), type = 7L))
          descriptions[[length(descriptions) + 1L]] <- data.frame(cohort = cohort, grouping = grouping,
            group = label, metric = metric, n = length(values), mean = mean(values), sd = stats::sd(values),
            median = stats::median(values), q25 = quartiles[[1L]], q75 = quartiles[[2L]],
            minimum = min(values), maximum = max(values), stringsAsFactors = FALSE)
        }
      }
    }
  }
  scales <- data.frame(scale = c("native_exact", "normalized_exact", "normalized_no_pseudocount_sensitivity"),
    gcol = c("logG_native", "logG_normalized_exact", "logG_normalized_no_pseudo"),
    dcol = c("logD_native", "logD_normalized_exact", "logD_normalized_no_pseudo"), stringsAsFactors = FALSE)
  decomposition <- list()
  for (cohort in unique(rows$cohort)) {
    data <- rows[rows$cohort == cohort, , drop = FALSE]
    uc <- data[data$diagnosis == "UC", , drop = FALSE]
    cd <- data[data$diagnosis == "CD", , drop = FALSE]
    delta <- function(column) mean(uc[[column]]) - mean(cd[[column]])
    for (i in seq_len(nrow(scales))) {
      dg <- delta(scales$gcol[[i]])
      dd <- delta(scales$dcol[[i]])
      ds <- delta("score")
      error <- dg - dd - ds
      if (i <= 2L && abs(error) > 1e-10) stop("Exact mean-score component decomposition failed.")
      decomposition[[length(decomposition) + 1L]] <- data.frame(cohort = cohort, role = data$role[[1L]],
        scale = scales$scale[[i]], normalized_unit = data$normalized_unit[[1L]], n_UC = nrow(uc), n_CD = nrow(cd),
        delta_mean_logG = dg, delta_mean_logD = dd, numerator_contribution = dg, denominator_contribution = -dd,
        sum_contributions = dg - dd, delta_mean_original_score = ds, difference_from_original_score = error,
        G_geometric_mean_ratio_UC_over_CD = 2^dg, D_geometric_mean_ratio_UC_over_CD = 2^dd,
        stringsAsFactors = FALSE)
    }
  }
  decomposition <- do.call(rbind, decomposition)
  fields <- c("cohort", "role", "scale", "n_UC", "n_CD", "delta_mean_logG", "delta_mean_logD",
    "delta_mean_original_score", "sum_contributions", "G_geometric_mean_ratio_UC_over_CD", "D_geometric_mean_ratio_UC_over_CD")
  list(participants = rows, descriptions = do.call(rbind, descriptions), decomposition = decomposition,
    table_s5 = decomposition[, fields, drop = FALSE])
}

sensitivity_rank_auc <- function(score, diagnosis) {
  if (any(!is.finite(score)) || anyNA(diagnosis) || any(!diagnosis %in% c("UC", "CD"))) stop("Invalid sensitivity AUC inputs.")
  positive <- diagnosis == "UC"
  n1 <- sum(positive)
  n0 <- sum(!positive)
  if (!n1 || !n0) stop("Both diseases are required for a sensitivity AUC.")
  rank_auc <- (sum(rank(score, ties.method = "average")[positive]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  reference <- as.numeric(pROC::auc(uc_roc(score, as.integer(positive))))
  if (abs(rank_auc - reference) > 1e-12) stop("Rank-based sensitivity AUC disagrees with pROC.")
  rank_auc
}

sensitivity_unit_analysis <- function(development, external) {
  score <- function(a, b, prior) log2((a + prior) / (b + prior))
  one <- function(cohort, unit, a, b, diagnosis) {
    if (any(a <= 0) || any(b <= 0)) stop("Positive sums are required for the zero-pseudocount sensitivity.")
    frozen <- score(a, b, 0.5)
    baseline <- sensitivity_rank_auc(frozen, diagnosis)
    do.call(rbind, lapply(c(0, 0.1, 0.5, 1), function(prior) {
      values <- score(a, b, prior)
      observed <- sensitivity_rank_auc(values, diagnosis)
      data.frame(cohort = cohort, unit = unit, n = length(values), n_uc = sum(diagnosis == "UC"),
        n_cd = sum(diagnosis == "CD"), pseudocount = prior, auc_uc_positive = observed,
        auc_change_from_frozen = observed - baseline, spearman_with_frozen = stats::cor(values, frozen, method = "spearman"),
        max_abs_score_change = max(abs(values - frozen)), stringsAsFactors = FALSE)
    }))
  }
  ext <- external$external$all
  ext <- ext[ext$diagnosis %in% c("UC", "CD"), , drop = FALSE]
  context <- external$context$all
  context <- context[context$diagnosis %in% c("UC", "CD"), , drop = FALSE]
  pseudocount <- rbind(one("GSE193677_DEVELOPMENT", "AGGREGATED_COUNTS", development$G, development$D, development$diagnosis),
    one("GSE117993_PRIMARY_EXTERNAL", "TPM", ext$igg_sum, ext$igam_sum, ext$diagnosis),
    one("GSE111889_CONTEXT", "AGGREGATED_COUNTS", context$igg_sum, context$igam_sum, context$diagnosis))
  k <- 1e6 / development$effective_library_size
  primary <- score(development$G, development$D, 0.5)
  variants <- list(counts_prior_0.5 = primary,
    tmm_cpm_prior_0.5 = score(development$G * k, development$D * k, 0.5),
    tmm_cpm_prior_scaled = score(development$G * k, development$D * k, 0.5 * k))
  scaling <- do.call(rbind, lapply(names(variants), function(label) {
    value <- variants[[label]]
    observed <- sensitivity_rank_auc(value, development$diagnosis)
    data.frame(variant = label, auc_uc_positive = observed,
      auc_change_from_frozen = observed - sensitivity_rank_auc(primary, development$diagnosis),
      spearman_with_frozen = stats::cor(value, primary, method = "spearman"),
      max_abs_score_change = max(abs(value - primary)), stringsAsFactors = FALSE)
  }))
  if (max(abs(variants$tmm_cpm_prior_scaled - primary)) > 1e-12) stop("Scaled-pseudocount identity failed.")
  # Historical illustrative lengths; these are not estimated effective lengths.
  panel <- development_gene_panel()[seq_len(7L), , drop = FALSE]
  lengths_bp <- c(1590, 1560, 1707, 1554, 1587, 1584, 1876)
  lengths <- data.frame(gene_symbol = panel$gene_symbol, gene_id = panel$gene_id,
    transcript_length_bp = lengths_bp, role = c(rep("IGG_NUMERATOR", 4L), rep("IGAM_DENOMINATOR", 3L)),
    stringsAsFactors = FALSE)
  rate <- development$panel / lengths_bp
  tpm_like <- sweep(rate, 2L, development$raw_library_size / 1000, "/") * 1e6
  cpm_like <- sweep(development$panel, 2L, development$effective_library_size / 1e6, "/")
  from_matrix <- function(block) score(colSums(block[seq_len(4L), , drop = FALSE]), colSums(block[5L:7L, , drop = FALSE]), 0.5)
  emulated <- data.frame(participant_id = development$table$participant_id, score_counts = from_matrix(development$panel),
    score_cpm_like = from_matrix(cpm_like), score_tpm_like = from_matrix(tpm_like),
    is_uc = development$is_uc, y_cd = development$y_cd, IGG_IGAM_LOG_RATIO = development$score)
  robustness <- do.call(rbind, lapply(c("score_counts", "score_cpm_like", "score_tpm_like"), function(label) {
    result <- uc_auc_summary(emulated[[label]], development$is_uc)
    data.frame(score_variant = label, auc_uc_positive = result$auc, auc_ci_low = result$ci95_low,
      auc_ci_high = result$ci95_high, auc_difference_from_counts = result$auc - sensitivity_rank_auc(primary, development$diagnosis),
      spearman_with_counts = stats::cor(emulated[[label]], primary, method = "spearman"),
      pearson_with_counts = stats::cor(emulated[[label]], primary), stringsAsFactors = FALSE)
  }))
  list(pseudocount = pseudocount, scaling = scaling, illustrative_lengths = lengths,
    illustrative_scores = emulated, illustrative_robustness = robustness,
    interpretation = "Retrospective diagnostics; the historical length-weighting illustration is not actual TPM and does not establish measurement equivalence.")
}

sensitivity_cell_panel <- function(path, discovery, emit) {
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  obs <- h5[["obs"]]
  native <- lapply(c("Patient", "sample_id", "Disease", "Site", "Treatment", "major", "minor", "Inflammation", "total_counts"),
    function(field) read_h5ad_column(obs, field))
  names(native) <- c("patient", "sample_id", "disease_raw", "site_raw", "treatment_raw", "major", "minor", "inflammation", "total_counts")
  n_cells <- length(native$patient)
  if (any(lengths(native) != n_cells)) stop("H5AD cell metadata lengths differ.")
  cell_id <- as.character(read_h5ad_column(obs, if ("_index" %in% names(obs)) "_index" else "cellbarcode"))
  if (length(cell_id) != n_cells) stop("H5AD cell identifiers are misaligned.")
  external_unique_key(cell_id, "native myeloid cell identifiers")
  metadata <- data.frame(cell_index = seq_len(n_cells), patient = as.character(native$patient),
    sample_id = as.character(native$sample_id), disease_harmonized = normalize_disease(native$disease_raw),
    site_harmonized = normalize_site(native$site_raw), treatment_harmonized = normalize_treatment(native$treatment_raw),
    major = as.character(native$major), minor = as.character(native$minor), stringsAsFactors = FALSE)
  candidate <- metadata$disease_harmonized %in% c("UC", "CD") & !is.na(metadata$site_harmonized) &
    metadata$site_harmonized == "LeftColon" & metadata$treatment_harmonized != "Post" &
    !is.na(metadata$patient) & nzchar(metadata$patient)
  selected <- which(candidate)
  cells <- metadata[selected, , drop = FALSE]
  previous <- discovery$candidate_cell_metadata$myeloid
  if (nrow(cells) != 3601L || !identical(cells$cell_index, previous$cell_index) ||
      !identical(cells$patient, previous$patient) || !identical(cells$sample_id, previous$sample_id) ||
      !identical(cells$disease_harmonized, previous$disease) || !identical(cells$minor, previous$minor)) {
    stop("New native cell selection disagrees with the validated discovery selection.")
  }
  cells$lineage <- "myeloid"
  cells$cell_id <- cell_id[selected]
  cells$inflammation_native <- as.character(native$inflammation[selected])
  cells$total_counts_native <- as.numeric(native$total_counts[selected])
  if (any(!is.finite(cells$total_counts_native)) || any(cells$total_counts_native <= 0)) stop("Invalid native cell libraries.")
  identifiers <- as.character(read_h5ad_column(h5[["var"]], "gene_id"))
  symbols <- as.character(read_h5ad_column(h5[["var"]], "gene_symbol"))
  if (length(identifiers) != 33075L || length(symbols) != length(identifiers) || anyNA(identifiers) || anyNA(symbols) ||
      anyDuplicated(identifiers) || anyDuplicated(symbols)) stop("Invalid native myeloid feature annotation.")
  panel <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM", "IGKC", "IGLC2",
    "JCHAIN", "MZB1", "DERL3", "SDC1", "MS4A1", "CD79A", "CD79B", "C1QA", "C1QB", "C1QC",
    "CD68", "LST1", "TYROBP", "FCER1G", "LYZ", "CD3D", "EPCAM", "COL1A1")
  panel_rows <- match(panel, symbols)
  present <- !is.na(panel_rows)
  panel_ids <- rep("MISSING", length(panel))
  panel_ids[present] <- identifiers[panel_rows[present]]
  if (!all(c("IGHG3", "C1QA", "CD68") %in% panel[present]) ||
      sub("\\..*$", "", panel_ids[panel == "IGHG3"]) != "ENSG00000211897") {
    stop("A required cell-signal marker is missing or mismapped.")
  }
  expression <- h5[["X"]]
  if (!inherits(expression, "H5Group") || !all(c("data", "indices", "indptr") %in% names(expression))) stop("Expected native CSR counts.")
  pointers <- as.numeric(expression[["indptr"]]$read())
  if (length(pointers) != n_cells + 1L || pointers[[1L]] != 0 || any(!is.finite(pointers)) ||
      any(diff(pointers) < 0) || any(pointers != floor(pointers))) stop("Invalid cell-signal CSR pointers.")
  values_ds <- expression[["data"]]
  indices_ds <- expression[["indices"]]
  extracted <- matrix(NA_real_, nrow(cells), length(panel))
  extracted[, present] <- 0
  totals <- numeric(nrow(cells))
  detected <- integer(nrow(cells))
  for (i in seq_len(nrow(cells))) {
    row <- selected[[i]]
    start <- pointers[[row]] + 1
    end <- pointers[[row + 1L]]
    value <- numeric()
    gene <- integer()
    if (end >= start) {
      if (end > .Machine$integer.max) stop("CSR index exceeds supported integer range.")
      at <- seq.int(as.integer(start), as.integer(end))
      value <- as.numeric(values_ds[at])
      gene <- as.integer(indices_ds[at]) + 1L
      if (length(value) != length(gene) || any(!is.finite(value)) || any(value < 0) ||
          any(abs(value - round(value)) > 1e-8) || anyNA(gene) || any(gene < 1L | gene > length(identifiers)) ||
          anyDuplicated(gene)) stop("Invalid native raw-count cell row.")
      hit <- match(panel_rows[present], gene)
      block <- numeric(sum(present))
      block[!is.na(hit)] <- value[hit[!is.na(hit)]]
      extracted[i, present] <- block
    }
    totals[[i]] <- sum(value)
    detected[[i]] <- sum(value > 0)
    if (i %% 500L == 0L || i == nrow(cells)) emit(paste("INFO | cell_panel | Native cells read:", i, "/", nrow(cells)))
  }
  if (any(totals != cells$total_counts_native)) stop("Native cell row sums do not conserve total_counts.")
  columns <- paste0(panel, "__", panel_ids, "__count")
  colnames(extracted) <- columns
  macro <- !is.na(cells$minor) & cells$minor == "Macro"
  macro_n <- table(cells$patient[macro])
  cells$n_candidate_macro_for_donor <- as.integer(macro_n[match(cells$patient, names(macro_n))])
  cells$primary_macro_eligible <- macro & !is.na(cells$n_candidate_macro_for_donor) & cells$n_candidate_macro_for_donor >= 10L
  primary <- discovery$models$MACROPHAGE_PRIMARY$donor_metadata
  if (sum(macro) != 880L || sum(cells$primary_macro_eligible) != 848L ||
      !setequal(unique(cells$patient[cells$primary_macro_eligible]), primary$patient)) stop("Cell-level primary macrophage eligibility differs.")
  pseudobulk <- discovery$pseudobulk$MACROPHAGE_PRIMARY
  donors <- colnames(pseudobulk$counts)
  rebuilt <- vapply(donors, function(donor) colSums(extracted[macro & cells$patient == donor, present, drop = FALSE]), numeric(sum(present)))
  pb_rows <- match(panel_ids[present], pseudobulk$gene_id)
  if (anyNA(pb_rows) || any(rebuilt != pseudobulk$counts[pb_rows, donors, drop = FALSE])) {
    stop("Cell panel sums disagree with the newly computed native macrophage pseudobulk.")
  }
  cells$total_counts_rebuilt <- totals
  cells$n_detected_genes_rebuilt <- detected
  cells <- cbind(cells, as.data.frame(extracted, check.names = FALSE))
  use <- cells$primary_macro_eligible
  gene <- extracted[, match("IGHG3", panel)]
  macrophage_marker <- extracted[, match("C1QA", panel)] > 0 | extracted[, match("CD68", panel)] > 0
  detection <- do.call(rbind, lapply(c("UC", "CD"), function(disease) {
    selected <- use & cells$disease_harmonized == disease
    data.frame(disease = disease, n_cells = sum(selected), n_donors = length(unique(cells$patient[selected])),
      IGHG3_positive_cells = sum(gene[selected] > 0), IGHG3_counts = sum(gene[selected]),
      positive_with_C1QA_or_CD68 = sum(gene[selected] > 0 & macrophage_marker[selected]))
  }))
  donor_summary <- do.call(rbind, lapply(as.character(primary$patient), function(donor) {
    selected <- use & cells$patient == donor
    data.frame(patient = donor, disease = unique(cells$disease_harmonized[selected]), n_cells = sum(selected),
      IGHG3_positive_cells = sum(gene[selected] > 0), IGHG3_counts = sum(gene[selected]), stringsAsFactors = FALSE)
  }))
  uc_counts <- donor_summary$IGHG3_counts[donor_summary$disease == "UC"]
  top_three_fraction <- sum(utils::head(sort(uc_counts, decreasing = TRUE), 3L)) / sum(uc_counts)
  list(cells = cells, panel_mapping = data.frame(gene_symbol = panel, gene_id = panel_ids, present = present),
    detection = detection, donors = donor_summary, uc_top_three_count_fraction = top_three_fraction,
    positive_co_detection = sum(use & gene > 0 & macrophage_marker), maximum_library_discrepancy = max(abs(totals - cells$total_counts_native)))
}

sensitivity_donor_influence <- function(discovery, record, emit) {
  model <- discovery$models$MACROPHAGE_PRIMARY
  dge <- model$dge
  design <- model$design
  donor <- colnames(dge$counts)
  metadata <- model$donor_metadata
  metadata <- metadata[match(donor, metadata$patient), , drop = FALSE]
  if (nrow(dge$counts) != 2101L || length(donor) != 24L || !identical(donor, rownames(design)) ||
      qr(design)$rank != ncol(design) || anyNA(metadata$patient) ||
      sum(metadata$disease == "UC") != 12L || sum(metadata$disease == "CD") != 12L) {
    stop("The validated macrophage model lacks the required gene universe or donor geometry.")
  }
  coefficient <- match("diseaseUC", colnames(design))
  gene <- which(sub("\\..*$", "", rownames(dge$counts)) == "ENSG00000211897")
  if (is.na(coefficient) || length(gene) != 1L) stop("Missing donor-influence coefficient or IGHG3 row.")
  # Recreate the pre-dispersion object so cached fitted components cannot persist.
  clean <- edgeR::DGEList(counts = dge$counts, genes = dge$genes)
  clean$samples <- dge$samples
  fit_one <- function(omit = NA_character_) {
    keep <- if (is.na(omit)) rep(TRUE, length(donor)) else donor != omit
    row <- data.frame(omitted_donor = if (is.na(omit)) "NONE_BASELINE" else omit,
      omitted_disease = if (is.na(omit)) NA_character_ else as.character(metadata$disease[donor == omit]),
      n_donors = sum(keep), n_uc = sum(metadata$disease[keep] == "UC"), n_cd = sum(metadata$disease[keep] == "CD"),
      n_genes = nrow(clean$counts), logFC = NA_real_, PValue = NA_real_, FDR_BH_within_fit = NA_real_,
      status = "NOT_RUN", warning_text = "", error_text = "", stringsAsFactors = FALSE)
    warnings <- character()
    table <- NULL
    tryCatch(withCallingHandlers({
      reduced_design <- design[keep, , drop = FALSE]
      if (qr(reduced_design)$rank != ncol(reduced_design)) stop("Donor-deletion design is rank deficient.")
      y <- if (is.na(omit)) clean else edgeR::calcNormFactors(clean[, keep, keep.lib.sizes = FALSE], method = "TMM")
      y <- edgeR::estimateDisp(y, design = reduced_design)
      fitted <- edgeR::glmQLFit(y, design = reduced_design, robust = TRUE)
      test <- edgeR::glmQLFTest(fitted, coef = coefficient)
      table <- result_table(test, "MACROPHAGE_PRIMARY", "PRIMARY", "QLF")
      target <- table[table$gene_id == rownames(clean$counts)[gene], , drop = FALSE]
      if (nrow(table) != 2101L || nrow(target) != 1L || any(!is.finite(c(target$logFC, target$PValue, target$FDR)))) {
        stop("Unexpected donor-influence test results.")
      }
      row$logFC <- target$logFC
      row$PValue <- target$PValue
      row$FDR_BH_within_fit <- target$FDR
      row$status <- "OK"
    }, warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }), error = function(error) {
      row$status <<- "FAILED"
      row$error_text <<- conditionMessage(error)
    })
    row$warning_text <- paste(unique(warnings), collapse = " | ")
    list(row = row, table = table)
  }
  emit("INFO | donor_influence | Recomputing the 2101-gene baseline with the validated frozen normalization and design.")
  baseline <- fit_one()
  if (baseline$row$status != "OK") stop("Donor-influence baseline failed: ", baseline$row$error_text)
  if (!isTRUE(compare_table(baseline$table, model$table, "gene_id", "donor_influence_full_baseline", record))) {
    stop("The donor-influence baseline does not reproduce the validated discovery model.")
  }
  results <- list(baseline$row)
  for (i in seq_along(donor)) {
    emit(paste("INFO | donor_influence | Omitting donor", donor[[i]], ":", i, "/", length(donor)))
    result <- fit_one(donor[[i]])$row
    results[[i + 1L]] <- result
    if (result$status != "OK") record("FAIL", paste0("donor_fit_", donor[[i]]), result$error_text)
    if (nzchar(result$warning_text)) record("WARN", paste0("donor_warning_", donor[[i]]), result$warning_text)
  }
  table <- do.call(rbind, results)
  if (nzchar(baseline$row$warning_text)) record("WARN", "donor_baseline_warning", baseline$row$warning_text)
  deleted <- table[table$omitted_donor != "NONE_BASELINE", , drop = FALSE]
  passed <- nrow(deleted) == 24L && all(deleted$status == "OK") && all(deleted$n_donors == 23L) &&
    all(deleted$n_genes == 2101L) && all(deleted$logFC > 0)
  record(if (passed) "PASS" else "FAIL", "donor_influence_geometry", "All 24 donor deletions; original 2101 genes and design columns; fresh library sizes, TMM, and dispersion in each deletion.")
  ranges <- data.frame(metric = c("minimum_logFC", "maximum_logFC", "minimum_FDR", "maximum_FDR"),
    observed = c(min(deleted$logFC), max(deleted$logFC), min(deleted$FDR_BH_within_fit), max(deleted$FDR_BH_within_fit)),
    reported = c(3.640, 5.502, 1.67e-7, 0.001891),
    comparison = c("round to 3 decimal places", "round to 3 decimal places", "3 significant digits", "round to 6 decimal places"),
    stringsAsFactors = FALSE)
  rounded <- c(round(ranges$observed[1L:2L], 3L), signif(ranges$observed[[3L]], 3L), round(ranges$observed[[4L]], 6L))
  ranges$matches_reported <- is.finite(rounded) & rounded == ranges$reported
  record(if (all(ranges$matches_reported)) "PASS" else "FAIL", "donor_influence_reported_ranges",
    "Compared at the precision printed in the current manuscript; the full historical 25-row result export is unavailable in the supplied archives.")
  list(results = table, reported_ranges = ranges,
    comparison_scope = "Full newly computed baseline compared across 2101 genes; deletion ranges compared at manuscript precision; no claim of full historical deletion-table equality.")
}

check_sensitivity_helpers <- function() {
  G <- c(4, 1)
  D <- c(2, 8)
  k <- c(0.01, 5)
  saved <- log2((G + 0.5) / (D + 0.5))
  rows <- sensitivity_component_rows("fixture", "fixture", c("a", "b"), c("UC", "CD"), c("UC", "CD"),
    "counts", "scaled", G, D, k, saved)
  scaled <- log2((k * G + k * 0.5) / (k * D + k * 0.5))
  fixed <- log2((k * G + 0.5) / (k * D + 0.5))
  if (max(abs(scaled - saved), abs(rows$identity_error)) > 1e-12 || max(abs(fixed - saved)) < 0.1) {
    stop("Component identity and scaled-pseudocount self-check failed.")
  }
  for (fixture in list(c(1, 2, 3, 4), c(4, 3, 2, 1), c(1, 1, 1, 1))) {
    observed <- sensitivity_rank_auc(fixture, c("CD", "CD", "UC", "UC"))
    expected <- if (length(unique(fixture)) == 1L) 0.5 else if (fixture[[1L]] < fixture[[4L]]) 1 else 0
    if (observed != expected) stop("Sensitivity rank AUC self-check failed.")
  }
  invisible(TRUE)
}

validate_sensitivity <- function(bundle, project, record) {
  archive <- file.path(project, "reference", "sensitivity_expected.zip")
  files <- c("component.tsv", "cells.tsv.gz", "33c13_participant_components.tsv", "33c13_group_descriptions.tsv",
    "33c13_UC_CD_decomposition.tsv", "30d_component_transcript_lengths.tsv", "30d_emulated_scores.tsv", "30d_unit_robustness.tsv")
  if (!identical(unname(tools::md5sum(archive)), "c3f405f14cfa76a12a171c63a5ce28a7")) stop("The sensitivity comparison archive failed its checksum.")
  listing <- utils::unzip(archive, list = TRUE)$Name
  if (anyDuplicated(listing) || !setequal(listing, files)) stop("Unexpected sensitivity comparison archive members.")
  directory <- tempfile(pattern = "ibd_sensitivity_expected_")
  if (!dir.create(directory)) stop("Cannot create a temporary sensitivity comparison directory.")
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(archive, files = files, exdir = directory, unzip = "internal")
  read <- function(name) {
    connection <- if (endsWith(name, ".gz")) gzfile(file.path(directory, name), "rt") else file(file.path(directory, name), "rt")
    on.exit(close(connection), add = TRUE)
    utils::read.delim(connection, check.names = FALSE, stringsAsFactors = FALSE)
  }
  outcomes <- logical()
  compare <- function(observed, name, keys, fields = NULL, exact = character()) {
    expected <- read(name)
    if (is.null(fields)) fields <- names(expected)
    outcomes <<- c(outcomes, compare_table(observed, expected, keys, name, record, fields = fields, exact_fields = exact))
  }
  component <- bundle$development_components
  native_fields <- c("n_biopsies", "raw_library_size", names(component)[endsWith(names(component), "__count")])
  compare(component, "component.tsv", "participant_id", exact = native_fields)
  cells <- bundle$cell_signal$cells
  fields <- names(cells)
  exact <- c("cell_index", "n_candidate_macro_for_donor", "total_counts_native", "total_counts_rebuilt", "n_detected_genes_rebuilt",
    names(cells)[endsWith(names(cells), "__count")])
  compare(cells, "cells.tsv.gz", "cell_id", fields = fields, exact = exact)
  compare(bundle$components$participants, "33c13_participant_components.tsv", c("cohort", "participant_id"))
  compare(bundle$components$descriptions, "33c13_group_descriptions.tsv", c("cohort", "grouping", "group", "metric"), exact = "n")
  compare(bundle$components$decomposition, "33c13_UC_CD_decomposition.tsv", c("cohort", "scale"), exact = c("n_UC", "n_CD"))
  compare(bundle$unit_sensitivity$illustrative_lengths, "30d_component_transcript_lengths.tsv", "gene_id", exact = "transcript_length_bp")
  compare(bundle$unit_sensitivity$illustrative_scores, "30d_emulated_scores.tsv", "participant_id", exact = c("is_uc", "y_cd"))
  robustness <- bundle$unit_sensitivity$illustrative_robustness
  numeric_fields <- setdiff(names(robustness), "score_variant")
  robustness <- development_archive_view(robustness, decimal = stats::setNames(rep(6L, length(numeric_fields)), numeric_fields))
  compare(robustness, "30d_unit_robustness.tsv", "score_variant")
  final_zip <- file.path(project, "reference", "tables_reference.zip")
  mapping <- utils::read.delim(file.path(project, "config", "tables.tsv"), colClasses = "character", check.names = FALSE)
  expected_md5 <- unique(mapping$reference_zip_md5)
  if (length(expected_md5) != 1L || !identical(unname(tools::md5sum(final_zip)), expected_md5)) stop("The supplementary reference archive failed its checksum.")
  final <- "Table_S5_component_decomposition.tsv"
  utils::unzip(final_zip, files = final, exdir = directory, unzip = "internal")
  outcomes <- c(outcomes, compare_table(bundle$components$table_s5, read(final), c("cohort", "scale"),
    "current_Table_S5_complete", record, exact_fields = c("n_UC", "n_CD")))
  # These are comparison targets from the current article, never analysis inputs.
  detection <- bundle$cell_signal$detection
  expected_detection <- data.frame(disease = c("UC", "CD"), n_cells = c(339L, 509L),
    n_donors = c(12L, 12L), IGHG3_positive_cells = c(105L, 27L))
  outcomes <- c(outcomes, compare_table(detection, expected_detection, "disease", "current_article_cell_detection", record,
    exact_fields = c("n_cells", "n_donors", "IGHG3_positive_cells")))
  co_detection_passed <- bundle$cell_signal$positive_co_detection == 130L
  concentration_passed <- round(100 * bundle$cell_signal$uc_top_three_count_fraction, 1L) == 78.5
  record(if (co_detection_passed) "PASS" else "FAIL", "current_article_marker_codetection", "IGHG3-positive primary macrophages co-expressing C1QA or CD68; exact target 130.")
  record(if (concentration_passed) "PASS" else "FAIL", "current_article_donor_concentration", "Top three UC donors' IGHG3 count fraction compared at the reported one-decimal percent precision.")
  all(c(outcomes, co_detection_passed, concentration_passed))
}

run_sensitivity <- function(context) {
  project <- context$project
  work <- file.path(project, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  report_path <- file.path(work, "sensitivity_report.txt")
  connection <- file(report_path, open = "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  counters <- c(PASS = 0L, WARN = 0L, FAIL = 0L)
  emit <- function(lines) {
    lines <- as.character(lines)
    writeLines(lines, connection)
    flush(connection)
    cat(paste(lines, collapse = "\n"), "\n", sep = "")
  }
  record <- function(status, label, detail) {
    counters[[status]] <<- counters[[status]] + 1L
    emit(paste(status, label, detail, sep = " | "))
  }
  status <- 1L
  tryCatch(withCallingHandlers({
    emit(c("IBD component and sensitivity reproduction report", "Code revision: stage06.0",
      paste("Started UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      "Analysis inputs: validated newly computed discovery/development/external caches and the native myeloid H5AD.",
      "All component, unit-sensitivity, cell-signal and donor-influence computations are in R.",
      "Post hoc status and the original score, gene universe, disease direction, and donor rules are retained.",
      "Archived tables and manuscript numbers are comparison targets only."))
    check_validation_logic()
    record("PASS", "validator_self_checks", "Shared exact-count, identifier, numerical, missingness, and probability checks.")
    required <- c(edgeR = "4.4.2", limma = "3.62.2", statmod = "1.5.2", hdf5r = "1.3.12", Matrix = "1.7-1", pROC = "1.19.0.1")
    for (package in names(required)) {
      if (!requireNamespace(package, quietly = TRUE) || utils::packageVersion(package) != package_version(required[[package]])) {
        stop("Missing or incompatible package: ", package, ".")
      }
      record("PASS", paste0("package_", package), as.character(utils::packageVersion(package)))
    }
    check_sensitivity_helpers()
    record("PASS", "sensitivity_helper_self_checks", "Exact component identities, scaled versus fixed pseudocounts, and perfect/reversed/tied rank AUCs.")
    if (getRversion() != package_version("4.4.2")) record("WARN", "R_version", paste("Observed", as.character(getRversion()), "; archived analyses used 4.4.2."))
    if (!isTRUE(context$settings$verify_md5)) stop("Sensitivity reconstruction requires checksum verification.")
    manifest <- context$read_manifest(file.path(project, "config", "inputs.tsv"))
    cache_md5 <- character()
    read_cache <- function(stage, revision, identifiers) {
      path <- file.path(work, paste0(stage, ".rds"))
      if (!file.exists(path)) stop("Missing prerequisite cache: ", stage, ".")
      emit(paste("INFO | cache | Reading the computed", stage, "bundle."))
      bundle <- readRDS(path)
      if (!isTRUE(bundle$validation_passed) || !identical(bundle$code_revision, revision)) stop("Unvalidated or incompatible cache: ", stage, ".")
      for (identifier in identifiers) {
        specification <- manifest[manifest$input_id == identifier, , drop = FALSE]
        if (nrow(specification) != 1L || !identical(unname(bundle$source_md5[identifier]), specification$expected_md5[[1L]])) {
          stop("Cache source identities disagree with the input definitions: ", stage, ".")
        }
      }
      cache_md5[stage] <<- unname(tools::md5sum(path))
      record("PASS", paste0(stage, "_cache"), paste("Validated", revision, "; cache MD5:", cache_md5[[stage]]))
      bundle
    }
    discovery <- read_cache("discovery", "stage02.0", c("discovery_myeloid", "discovery_fibperi"))
    development <- read_cache("development", "stage04.2", c("development_counts", "gse193677_sample_metadata"))
    external <- read_cache("external", "stage05.0", c("external_archive", "gse117993_sample_metadata", "context_counts",
      "context_metadata", "nonestimable_archive", "nonestimable_rpkm"))
    if (!identical(external$development_cache_md5, unname(cache_md5["development"]))) stop("The external bundle was computed from a different development cache.")
    development_components <- sensitivity_development_components(development)
    rm(development)
    invisible(gc())
    components <- sensitivity_component_analysis(development_components, external)
    record("PASS", "component_identities", paste("Participants:", nrow(components$participants),
      "; maximum score reconstruction error:", format(max(abs(components$participants$reconstruction_error)), scientific = TRUE),
      "; maximum normalized identity error:", format(max(abs(components$participants$identity_error)), scientific = TRUE)))
    units <- sensitivity_unit_analysis(development_components, external)
    record("PASS", "unit_sensitivity_calculations", "Four native pseudocounts in each cohort and three development scaling variants; the original score is unchanged.")
    emit("INFO | historical_length_illustration | Reproduced as a historical numerical illustration only; its denominator is not actual all-gene length-normalized TPM.")
    specification <- manifest[manifest$input_id == "discovery_myeloid", , drop = FALSE]
    if (nrow(specification) != 1L || !context$safe_relative(specification$relative_path[[1L]])) stop("Missing myeloid input definition.")
    path <- file.path(context$settings$legacy_root, specification$relative_path[[1L]])
    if ("discovery_myeloid" %in% names(context$settings$file_overrides)) path <- context$settings$file_overrides[["discovery_myeloid"]]
    path <- normalizePath(path.expand(path), winslash = "/", mustWork = TRUE)
    source_md5 <- unname(tools::md5sum(path))
    if (file.info(path)$size != as.numeric(specification$expected_bytes[[1L]]) || !identical(source_md5, specification$expected_md5[[1L]])) {
      stop("Native myeloid source bytes or checksum differ.")
    }
    record("PASS", "discovery_myeloid", paste("Archived source MD5:", source_md5))
    cells <- sensitivity_cell_panel(path, discovery, emit)
    record("PASS", "cell_panel_conservation", "3601 eligible cells; 27-gene panel; exact native library sums and panel-to-donor pseudobulk count conservation; 848 primary macrophages.")
    influence <- sensitivity_donor_influence(discovery, record, emit)
    bundle <- list(code_revision = "stage06.0", source_md5 = c(discovery_myeloid = source_md5), cache_md5 = cache_md5,
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), development_components = development_components$table,
      components = components, unit_sensitivity = units, cell_signal = cells, donor_influence = influence,
      supplementary = list(table_s5 = components$table_s5),
      settings = list(component_inference = "descriptive; no component significance tests", mean_difference = "UC minus CD",
        ratio_pseudocount = 0.5, auc_positive = "UC", auc_direction = "<", unit_variants = "retrospective diagnostics; not score selection",
        cell_inference = "descriptive detection and co-detection; donor remains the biological replicate",
        donor_influence = "fixed 2101-gene universe and original design columns; deletion-specific normalization and dispersion",
        numerical_atol = 1e-10, probability_atol = 1e-14, numerical_rtol = 1e-6), session = utils::sessionInfo())
    emit(c("", "Comparison with archived component/cell values and current supplementary Table S5"))
    passed <- validate_sensitivity(bundle, project, record)
    bundle$validation_passed <- isTRUE(passed) && counters[["FAIL"]] == 0L
    bundle$validation_counts <- counters
    cache <- file.path(work, "sensitivity.rds")
    saveRDS(bundle, cache, compress = "xz")
    emit(c("", "Newly computed component decomposition", utils::capture.output(print(components$table_s5, row.names = FALSE, digits = 10)),
      "", "Native-unit pseudocount sensitivity", utils::capture.output(print(units$pseudocount, row.names = FALSE, digits = 10)),
      "", "Primary macrophage detection", utils::capture.output(print(cells$detection, row.names = FALSE)),
      paste("Top-three UC donor count fraction:", format(cells$uc_top_three_count_fraction, digits = 12)),
      "", "Donor-deletion ranges", utils::capture.output(print(influence$reported_ranges, row.names = FALSE, digits = 10)),
      paste("Comparison scope:", influence$comparison_scope), paste("Cache:", normalizePath(cache, winslash = "/", mustWork = TRUE)),
      "The cache contains newly computed values at full precision. Later stages must check validation_passed.",
      "A fatal run does not refresh a previous cache; review the current report before advancing.",
      "", "R session information", utils::capture.output(utils::sessionInfo())))
    status <- if (isTRUE(bundle$validation_passed)) 0L else 1L
  }, warning = function(warning) {
    record("WARN", "R_warning", conditionMessage(warning))
    invokeRestart("muffleWarning")
  }), error = function(error) {
    record("FAIL", "sensitivity_error", conditionMessage(error))
    status <<- 1L
  })
  emit(c("", paste("Checks: PASS =", counters[["PASS"]], "; WARN =", counters[["WARN"]], "; FAIL =", counters[["FAIL"]]),
    paste("Sensitivity/component reproduction:", if (status == 0L) "PASS" else "BLOCKED"),
    "Final workbook, figures, and isolated publication-release verification: NOT_RUN",
    paste("Finished UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Report:", normalizePath(report_path, winslash = "/", mustWork = TRUE))))
  status
}
