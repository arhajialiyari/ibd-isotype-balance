# Reconstruct the GSE193677 development analyses from estimated source counts.
# Retains the manuscript analyses from stages 18, 19, 27, 27b, 28, 29, and 33c15.
# Shared input-label and comparison helpers are loaded from discovery.R by run.R.

development_gene_panel <- function() {
  panel <- data.frame(
    gene_symbol = c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHM",
      "IGHD", "IGKC", "IGLC1", "IGLC2", "IGLC3", "JCHAIN", "MZB1", "DERL3", "XBP1",
      "PRDM1", "SDC1", "TNFRSF17", "PTPRC", "EPCAM", "COL1A1", "CD68", "C1QA", "LYZ"),
    gene_id = c("ENSG00000211896", "ENSG00000211893", "ENSG00000211897", "ENSG00000211892",
      "ENSG00000211895", "ENSG00000211890", "ENSG00000211899", "ENSG00000211898",
      "ENSG00000211592", "ENSG00000211675", "ENSG00000211677", "ENSG00000211679",
      "ENSG00000132465", "ENSG00000170476", "ENSG00000099958", "ENSG00000100219",
      "ENSG00000057657", "ENSG00000115884", "ENSG00000048462", "ENSG00000081237",
      "ENSG00000119888", "ENSG00000108821", "ENSG00000129226", "ENSG00000173372", "ENSG00000090382"),
    gene_class = c(rep("IG_CONSTANT", 12L), rep("PLASMA_CELL_PROGRAM", 7L),
      rep("LINEAGE_REFERENCE", 3L), rep("MYELOID_REFERENCE", 3L)),
    role = c(rep("IgG numerator", 4L), rep("IgA plus IgM denominator", 3L),
      rep("Additional total Ig", 5L), rep("Plasma cell comparator", 7L),
      rep("Lineage reference", 6L)), stringsAsFactors = FALSE)
  if (anyDuplicated(panel$gene_id) || anyDuplicated(panel$gene_symbol)) stop("Duplicate panel definitions.")
  panel
}

development_characteristics <- function(metadata) {
  required <- c("sample_id", "title", "characteristics_raw")
  if (!all(required %in% names(metadata)) || anyNA(metadata$sample_id) ||
      anyDuplicated(metadata$sample_id) || anyNA(metadata$characteristics_raw)) {
    stop("Source development metadata lack unique GSMs or characteristics strings.")
  }
  missing_to_na <- function(value) {
    value <- trimws(as.character(value))
    value[value %in% c("", "NA", "na", "N/A", "NULL", ".")] <- NA_character_
    value
  }
  parts <- strsplit(metadata$characteristics_raw, "\\s*\\|\\s*")
  keys <- sort(unique(unlist(lapply(parts, function(value) {
    trimws(sub(":.*$", "", value[grepl(":", value, fixed = TRUE)]))
  }), use.names = FALSE)))
  wide <- data.frame(sample_id = metadata$sample_id, stringsAsFactors = FALSE)
  for (key in keys) {
    wide[[key]] <- vapply(parts, function(value) {
      value <- value[grepl(":", value, fixed = TRUE)]
      hit <- which(trimws(sub(":.*$", "", value)) == key)
      if (!length(hit)) return(NA_character_)
      values <- missing_to_na(sub("^[^:]*:", "", value[hit]))
      if (length(unique(values)) > 1L) stop("Conflicting repeated characteristic: ", key)
      values[[1L]]
    }, character(1))
  }
  wide
}

select_development_samples <- function(path, read_manifest) {
  metadata <- read_manifest(path)
  if (nrow(metadata) != 2490L || !"study_id" %in% names(metadata) ||
      anyNA(metadata$study_id) || any(metadata$study_id != "GSE193677")) {
    stop("The development source metadata do not describe the expected 2,490 samples.")
  }
  wide <- development_characteristics(metadata)
  needed <- c("ibd_disease", "regionre", "typere", "ibd_clinicianmeasure_inactive_active",
    "demographics_gender", "study_eligibility_age_at_endo", "ibd_endoseverity_4levels",
    "nancyindex", "crp_jjmgl_log2", "log2_fecalcalpro_mgperg")
  if (!all(needed %in% names(wide))) stop("Missing source development characteristics: ",
    paste(setdiff(needed, names(wide)), collapse = ", "))
  disease <- normalize_disease(wide$ibd_disease)
  anatomy <- normalize_site(wide$regionre)
  selected <- !is.na(anatomy) & anatomy == "Rectum" & disease %in% c("CD", "UC")
  primary <- metadata[selected, , drop = FALSE]
  wide <- wide[selected, , drop = FALSE]
  if (nrow(primary) != 679L) stop("The source selection does not contain 679 rectal CD/UC biopsies.")
  primary$diagnosis_harmonized <- disease[selected]
  primary$anatomy_harmonized <- anatomy[selected]
  title_matches <- regmatches(primary$title, regexec("MSCCR_reGRID_([0-9]+)", primary$title))
  primary$participant_id_candidate <- vapply(title_matches, function(value) {
    if (length(value) == 2L) value[[2L]] else NA_character_
  }, character(1))
  if (anyNA(primary$participant_id_candidate) ||
      length(unique(primary$participant_id_candidate)) != 665L) {
    stop("Source titles do not map the selected biopsies to 665 participants.")
  }
  # Reconstruct the historical clinical-activity label only for its audit.
  clinical <- wide$ibd_clinicianmeasure_inactive_active
  frozen <- rep(NA_character_, length(clinical))
  frozen[clinical %in% "Active"] <- "Active_or_Inflamed"
  frozen[clinical %in% "Inactive"] <- "Inactive_or_NonInflamed"
  primary$inflammation_harmonized <- frozen
  selection <- primary[, c("study_id", "sample_id", "title", "diagnosis_harmonized",
    "anatomy_harmonized", "inflammation_harmonized", "participant_id_candidate"), drop = FALSE]
  geometry <- do.call(rbind, lapply(c("CD", "UC"), function(diagnosis) {
    data <- selection[selection$diagnosis_harmonized == diagnosis, , drop = FALSE]
    data.frame(study_id = "GSE193677", anatomy = "Rectum", diagnosis = diagnosis,
      n_biopsies = nrow(data), n_participants = length(unique(data$participant_id_candidate)),
      repeated_biopsy_rows = sum(duplicated(data$participant_id_candidate)),
      inflammation_values = paste(sort(unique(data$inflammation_harmonized)), collapse = ";"),
      stringsAsFactors = FALSE)
  }))
  list(selection = selection, characteristics = wide, geometry = geometry)
}

import_development_counts <- function(path, selection, forward, emit) {
  read_header <- function(path) {
    connection <- gzfile(path, open = "rt", encoding = "UTF-8")
    on.exit(close(connection), add = TRUE)
    readLines(connection, n = 2L, warn = FALSE)
  }
  lines <- read_header(path)
  if (length(lines) != 2L) stop("The development source lacks a header and first gene row.")
  parse <- function(value) scan(text = value, what = character(), sep = "", quote = "\"", quiet = TRUE)
  header <- parse(lines[[1L]])
  if (length(header) != 2490L || anyDuplicated(header) || length(parse(lines[[2L]])) != 2491L) {
    stop("Unexpected development source header or data-field geometry.")
  }
  tokens <- trimws(sub(",.*$", "", selection$title))
  position <- match(tokens, header) + 1L
  if (anyNA(position) || anyDuplicated(position) || anyDuplicated(tokens)) {
    stop("Source-title tokens do not map one-to-one onto native count fields.")
  }
  mapping <- data.frame(sample_id = selection$sample_id, source_sample_token = tokens,
    participant_id = selection$participant_id_candidate, raw_data_field_position = position,
    header_token_verified = TRUE, stringsAsFactors = FALSE)
  mapping <- mapping[order(mapping$raw_data_field_position), , drop = FALSE]
  if (!identical(unname(header[mapping$raw_data_field_position - 1L]), mapping$source_sample_token)) {
    stop("Native field positions do not reproduce the selected source tokens.")
  }
  # Base-R decompression avoids a hidden R.utils or external gzip dependency.
  temporary <- tempfile(pattern = "ibd_development_counts_", fileext = ".txt")
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  decompress <- function() {
    input <- gzfile(path, open = "rb")
    on.exit(close(input), add = TRUE)
    output <- file(temporary, open = "wb")
    on.exit(close(output), add = TRUE)
    repeat {
      block <- readBin(input, what = "raw", n = 8L * 1024L * 1024L)
      if (!length(block)) break
      writeBin(block, output)
    }
  }
  emit("INFO | import | Decompressing the source into temporary storage; importing 679 declared biopsy columns.")
  decompress()
  input <- data.table::fread(temporary, header = FALSE, skip = 1L,
    select = c(1L, mapping$raw_data_field_position), sep = " ", quote = "\"",
    check.names = FALSE, showProgress = FALSE, data.table = FALSE)
  if (nrow(input) != 56632L || ncol(input) != 680L) {
    stop("Imported development counts do not have 56,632 genes and 679 biopsies.")
  }
  names(input) <- c("raw_gene_id", mapping$source_sample_token)
  raw_gene_id <- as.character(input[[1L]])
  if (anyNA(raw_gene_id) || any(!nzchar(raw_gene_id)) || anyDuplicated(raw_gene_id)) {
    stop("The source gene identifiers are incomplete or duplicated.")
  }
  selected_data <- input[, -1L, drop = FALSE]
  if (!all(vapply(selected_data, is.numeric, logical(1)))) stop("A source count column is nonnumeric.")
  biopsies <- as.matrix(selected_data)
  storage.mode(biopsies) <- "double"
  dimnames(biopsies) <- list(raw_gene_id, mapping$source_sample_token)
  rm(input, selected_data)
  unlink(temporary, force = TRUE)
  stable_gene_id <- sub("\\..*$", "", raw_gene_id)
  forward_rows <- vapply(forward$gene_id, function(id) {
    hits <- which(stable_gene_id == id)
    if (length(hits) != 1L) stop("A forward gene lacks unique native stable-Ensembl support: ", id)
    hits[[1L]]
  }, integer(1))
  fractional_entries <- 0L
  fractional_biopsies <- integer(ncol(biopsies))
  for (start in seq.int(1L, nrow(biopsies), by = 1000L)) {
    rows <- seq.int(start, min(start + 999L, nrow(biopsies)))
    block <- biopsies[rows, , drop = FALSE]
    if (any(!is.finite(block)) || any(block < 0)) stop("Nonfinite or negative estimated source counts.")
    fractional <- abs(block - round(block)) > 1e-8
    fractional_entries <- fractional_entries + sum(fractional)
    fractional_biopsies <- fractional_biopsies + colSums(fractional)
  }
  forward_fractional <- sum(abs(biopsies[forward_rows, , drop = FALSE] -
    round(biopsies[forward_rows, , drop = FALSE])) > 1e-8)
  participant_ids <- sort(unique(mapping$participant_id))
  counts <- matrix(0, nrow = nrow(biopsies), ncol = length(participant_ids),
    dimnames = list(raw_gene_id, participant_ids))
  aggregation <- vector("list", length(participant_ids))
  for (i in seq_along(participant_ids)) {
    participant <- participant_ids[[i]]
    index <- which(mapping$participant_id == participant)
    counts[, i] <- if (length(index) == 1L) biopsies[, index] else rowSums(biopsies[, index, drop = FALSE])
    aggregation[[i]] <- data.frame(participant_id = participant, n_biopsies_aggregated = length(index),
      source_sample_tokens = paste(mapping$source_sample_token[index], collapse = ";"),
      aggregation_rule = "GENE_WISE_RAW_SUM_BEFORE_NORMALIZATION",
      disease_label_used_in_aggregation = FALSE, inflammation_used_in_aggregation = FALSE,
      stringsAsFactors = FALSE)
  }
  source_total <- sum(biopsies)
  pooled_total <- sum(counts)
  if (!isTRUE(all.equal(pooled_total, source_total, tolerance = 0))) {
    stop("Estimated-count conservation failed during participant aggregation.")
  }
  rm(biopsies, block, fractional)
  invisible(gc())
  emit(paste("INFO | aggregation | Participants:", length(participant_ids),
    "; fractional source entries preserved:", fractional_entries, "; count conservation: exact."))
  library_size <- colSums(counts)
  if (any(!is.finite(library_size)) || any(library_size <= 0)) stop("A pooled library is not positive and finite.")
  dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts), method = "TMM")
  tmm <- data.frame(participant_id = participant_ids, raw_library_size = dge$samples$lib.size,
    tmm_norm_factor = dge$samples$norm.factors,
    effective_library_size = dge$samples$lib.size * dge$samples$norm.factors,
    stringsAsFactors = FALSE)
  if (any(!is.finite(tmm$effective_library_size)) || any(tmm$effective_library_size <= 0)) {
    stop("TMM normalization produced an invalid effective library size.")
  }
  forward_logcpm <- edgeR::cpm(dge[forward_rows, , keep.lib.sizes = TRUE],
    log = TRUE, prior.count = 0.5, normalized.lib.sizes = TRUE)
  if (nrow(forward) != 1L || forward$gene_symbol[[1L]] != "IGHG3") {
    stop("The validated forward set differs from the single-gene development endpoint in the manuscript.")
  }
  feature <- data.frame(participant_id = participant_ids, gene_id = forward$gene_id[[1L]],
    gene_symbol = forward$gene_symbol[[1L]], raw_count_sum = as.numeric(counts[forward_rows, ]),
    library_size_full_native = tmm$raw_library_size, TMM_norm_factor = tmm$tmm_norm_factor,
    log2CPM_prior_0.5 = as.numeric(forward_logcpm), orientation_rule = forward$orientation_rule[[1L]],
    stringsAsFactors = FALSE)
  numeric_audit <- data.frame(metric = c("raw_gene_rows", "frozen_biopsy_columns", "participant_columns",
    "fractional_biopsy_entries", "fractional_biopsies", "frozen_feature_fractional_entries",
    "all_counts_finite_nonnegative", "counts_rounded", "disease_used_in_aggregation"),
    value = c(nrow(counts), nrow(mapping), ncol(counts), fractional_entries,
      sum(fractional_biopsies > 0), forward_fractional, "TRUE", "FALSE", "FALSE"),
    stringsAsFactors = FALSE)
  list(counts = counts, raw_gene_id = raw_gene_id, stable_gene_id = stable_gene_id,
    mapping = mapping, aggregation = do.call(rbind, aggregation), tmm = tmm,
    forward_feature = feature, numeric_audit = numeric_audit)
}

development_covariates <- function(selected) {
  selection <- selected$selection
  wide <- selected$characteristics
  index <- order(selection$sample_id)
  selection <- selection[index, , drop = FALSE]
  wide <- wide[index, , drop = FALSE]
  if (!identical(selection$sample_id, wide$sample_id)) stop("Characteristics lost their GSM alignment.")
  site <- ifelse(wide$typere == "NonI", "SITE_NONINFLAMED",
                 ifelse(wide$typere == "I", "SITE_INFLAMED", NA_character_))
  title <- ifelse(grepl("Rectum\\s+NonI\\s+tissue", selection$title), "SITE_NONINFLAMED",
                  ifelse(grepl("Rectum\\s+I\\s+tissue", selection$title), "SITE_INFLAMED", NA_character_))
  if (anyNA(site) || anyNA(title) || any(site != title)) {
    stop("The site-inflammation field does not reproduce every selected source title.")
  }
  number <- function(value) suppressWarnings(as.numeric(value))
  biopsy <- data.frame(sample_id = selection$sample_id,
    participant_id = as.integer(selection$participant_id_candidate), diagnosis = selection$diagnosis_harmonized,
    site_status = site, clinical_activity = wide$ibd_clinicianmeasure_inactive_active,
    sex = wide$demographics_gender, age = number(wide$study_eligibility_age_at_endo),
    endoscopic_severity = wide$ibd_endoseverity_4levels, nancy_index = number(wide$nancyindex),
    crp_log2 = number(wide$crp_jjmgl_log2), calprotectin_log2 = number(wide$log2_fecalcalpro_mgperg),
    stringsAsFactors = FALSE)
  participant_ids <- sort(unique(biopsy$participant_id))
  groups <- lapply(participant_ids, function(id) which(biopsy$participant_id == id))
  first <- function(field) vapply(groups, function(index) {
    values <- biopsy[[field]][index]
    values <- values[!is.na(values)]
    if (length(values)) as.character(values[[1L]]) else NA_character_
  }, character(1))
  aggregate_number <- function(field, fun) vapply(groups, function(index) {
    values <- biopsy[[field]][index]
    values <- values[is.finite(values)]
    if (length(values)) as.numeric(fun(values)) else NA_real_
  }, numeric(1))
  for (index in groups) {
    if (length(unique(biopsy$diagnosis[index])) != 1L) stop("Diagnosis differs between a participant's biopsies.")
  }
  severity_order <- c("Inactive", "Mild", "Moderate", "Severe")
  biopsy$endoscopic_rank <- match(biopsy$endoscopic_severity, severity_order)
  participant <- data.frame(participant_id = participant_ids,
    diagnosis_harmonized = first("diagnosis"), n_rectum_biopsies = lengths(groups),
    site_inflammation = vapply(groups, function(index) {
      if (any(biopsy$site_status[index] == "SITE_INFLAMED")) "ANY_SITE_INFLAMED" else "ALL_SITE_NONINFLAMED"
    }, character(1)),
    clinical_activity = first("clinical_activity"), sex = first("sex"),
    endoscopic_severity = first("endoscopic_severity"), age = aggregate_number("age", mean),
    nancy_index = aggregate_number("nancy_index", max), crp_log2 = aggregate_number("crp_log2", mean),
    calprotectin_log2 = aggregate_number("calprotectin_log2", mean),
    endo_rank_max = aggregate_number("endoscopic_rank", max), stringsAsFactors = FALSE)
  participant$endoscopic_severity_4 <- factor(severity_order[participant$endo_rank_max], levels = severity_order)
  participant$site_inflammation <- factor(participant$site_inflammation,
    levels = c("ALL_SITE_NONINFLAMED", "ANY_SITE_INFLAMED"))
  participant$sex <- factor(participant$sex)
  participant$y_cd <- as.integer(participant$diagnosis_harmonized == "CD")
  participant$is_uc <- as.integer(participant$diagnosis_harmonized == "UC")
  participant$crp_z <- as.numeric(scale(participant$crp_log2))
  if (anyNA(participant[, c("sex", "age", "endoscopic_severity_4")])) {
    stop("The manuscript's complete-case development covariates are unexpectedly incomplete.")
  }
  biopsy$endoscopic_rank <- NULL
  list(biopsy = biopsy, participant = participant)
}

construct_development_features <- function(pooled, covariates, panel) {
  rows <- vapply(panel$gene_id, function(id) {
    index <- which(pooled$stable_gene_id == id)
    if (length(index) != 1L) stop("Panel gene lacks unique stable-Ensembl support: ", id)
    index[[1L]]
  }, integer(1))
  dge <- edgeR::DGEList(counts = pooled$counts[rows, , drop = FALSE],
    lib.size = pooled$tmm$raw_library_size, norm.factors = pooled$tmm$tmm_norm_factor)
  logcpm <- t(edgeR::cpm(dge, log = TRUE, prior.count = 0.5, normalized.lib.sizes = TRUE))
  colnames(logcpm) <- panel$gene_symbol
  native <- t(pooled$counts[rows, , drop = FALSE])
  colnames(native) <- panel$gene_symbol
  igg <- rowSums(native[, c("IGHG1", "IGHG2", "IGHG3", "IGHG4"), drop = FALSE])
  iga <- rowSums(native[, c("IGHA1", "IGHA2"), drop = FALSE])
  igm <- native[, "IGHM"]
  total <- rowSums(native[, panel$gene_class == "IG_CONSTANT", drop = FALSE])
  effective <- pooled$tmm$effective_library_size
  features <- data.frame(participant_id = as.integer(rownames(logcpm)), logcpm,
    TOTAL_IG = log2((total + 0.5) / (effective + 1) * 1e6),
    IGHG3_FRACTION_OF_IG = log2((native[, "IGHG3"] + 0.5) / (total - native[, "IGHG3"] + 0.5)),
    PLASMA_CELL_SCORE = rowMeans(scale(logcpm[, panel$gene_class == "PLASMA_CELL_PROGRAM", drop = FALSE])),
    IGG_IGA_LOG_RATIO = log2((igg + 0.5) / (iga + 0.5)),
    IGG_IGAM_LOG_RATIO = log2((igg + 0.5) / (iga + igm + 0.5)),
    IGG_FRACTION_OF_IG = log2((igg + 0.5) / (total - igg + 0.5)),
    IGG_ABUNDANCE = log2((igg + 0.5) / (effective + 1) * 1e6),
    G_native = igg, A_native = iga, M_native = igm, D_native = iga + igm, T_native = total,
    effective_library_size = effective, normalization_multiplier = 1e6 / effective,
    check.names = FALSE, stringsAsFactors = FALSE)
  work <- covariates$participant
  features <- features[match(work$participant_id, features$participant_id), , drop = FALSE]
  work <- cbind(work, features[, setdiff(names(features), "participant_id"), drop = FALSE])
  work$log2CPM_prior_0.5 <- pooled$forward_feature$log2CPM_prior_0.5[
    match(work$participant_id, pooled$forward_feature$participant_id)]
  work$ighg3_z <- as.numeric(scale(work$log2CPM_prior_0.5))
  work$total_ig_z <- as.numeric(scale(work$TOTAL_IG))
  work$plasma_z <- as.numeric(scale(work$PLASMA_CELL_SCORE))
  work$ig_fraction_z <- as.numeric(scale(work$IGHG3_FRACTION_OF_IG))
  if (nrow(work) != 665L || anyNA(work$IGHG3) ||
      max(abs(work$IGHG3 - work$log2CPM_prior_0.5)) > 1e-12) {
    stop("Panel normalization disagrees with the independently constructed forward-gene feature.")
  }
  numeric_fields <- c(panel$gene_symbol, "TOTAL_IG", "IGHG3_FRACTION_OF_IG", "PLASMA_CELL_SCORE",
    "IGG_IGA_LOG_RATIO", "IGG_IGAM_LOG_RATIO", "IGG_FRACTION_OF_IG", "IGG_ABUNDANCE")
  if (any(!is.finite(as.matrix(work[numeric_fields])))) stop("A development score is nonfinite.")
  work
}

uc_roc <- function(score, is_uc) {
  if (length(score) != length(is_uc) || anyNA(score) || any(!is.finite(score)) ||
      anyNA(is_uc) || !all(is_uc %in% c(0L, 1L)) ||
      sum(is_uc == 1L) < 2L || sum(is_uc == 0L) < 2L) {
    stop("ROC requires finite aligned scores and at least two participants of each disease.")
  }
  pROC::roc(response = is_uc, predictor = score, levels = c(0L, 1L), direction = "<", quiet = TRUE)
}

uc_auc_summary <- function(score, is_uc) {
  roc <- uc_roc(score, is_uc)
  ci <- as.numeric(pROC::ci.auc(roc, method = "delong", conf.level = 0.95))
  # Independent placement values retain half-credit for ties and the unclipped SE.
  placement <- outer(score[is_uc == 1L], score[is_uc == 0L],
    function(a, b) as.numeric(a > b) + 0.5 * as.numeric(a == b))
  standard_error <- sqrt(stats::var(rowMeans(placement)) / nrow(placement) +
    stats::var(colMeans(placement)) / ncol(placement))
  auc <- as.numeric(pROC::auc(roc))
  normal_ci <- pmin(1, pmax(0, auc + c(-1, 1) * stats::qnorm(0.975) * standard_error))
  if (abs(auc - mean(placement)) > 1e-12 || max(abs(normal_ci - ci[c(1L, 3L)])) > 1e-10) {
    stop("pROC and independent placement-based DeLong calculations disagree.")
  }
  data.frame(auc = auc, ci95_low = ci[[1L]], ci95_high = ci[[3L]],
    n_uc = sum(is_uc == 1L), n_cd = sum(is_uc == 0L), standard_error = standard_error)
}

cohort_score_comparisons <- function(scores, is_uc, cohort, role) {
  features <- c("primary_ratio", "IGHG3_abundance", "total_Ig_abundance", "IGHG3_fraction", "summed_IgG_abundance")
  if (!all(features %in% names(scores))) stop("Missing a declared cohort score or abundance comparator.")
  aucs <- do.call(rbind, lapply(features, function(feature) {
    result <- uc_auc_summary(scores[[feature]], is_uc)
    data.frame(cohort = cohort, role = role, feature = feature, result, stringsAsFactors = FALSE)
  }))
  primary_roc <- uc_roc(scores$primary_ratio, is_uc)
  paired <- do.call(rbind, lapply(features[-1L], function(comparator) {
    other_roc <- uc_roc(scores[[comparator]], is_uc)
    test <- pROC::roc.test(primary_roc, other_roc, paired = TRUE, method = "delong",
      alternative = "two.sided", conf.level = 0.95)
    interval <- as.numeric(test$conf.int)
    if (length(interval) != 2L || any(!is.finite(interval))) stop("Paired DeLong interval is unavailable.")
    data.frame(cohort = cohort, role = role, comparator = comparator,
      n_uc = sum(is_uc == 1L), n_cd = sum(is_uc == 0L), auc_primary = as.numeric(pROC::auc(primary_roc)),
      auc_comparator = as.numeric(pROC::auc(other_roc)),
      delta_primary_minus_comparator = as.numeric(pROC::auc(primary_roc) - pROC::auc(other_roc)),
      standard_error = diff(interval) / (2 * stats::qnorm(0.975)),
      ci95_low = interval[[1L]], ci95_high = interval[[2L]], p_two_sided = test$p.value,
      stringsAsFactors = FALSE)
  }))
  paired$p_holm_4_within_cohort <- stats::p.adjust(paired$p_two_sided, method = "holm", n = 4L)
  paired$p_holm_8_external_context_sensitivity <- NA_real_
  list(aucs = aucs, paired = paired)
}

development_background <- function(pooled, work, emit) {
  counts <- pooled$counts
  effective <- pooled$tmm$effective_library_size
  keep <- logical(nrow(counts))
  minimum <- ceiling(0.10 * ncol(counts))
  for (start in seq.int(1L, nrow(counts), by = 1000L)) {
    rows <- seq.int(start, min(start + 999L, nrow(counts)))
    cpm <- sweep(counts[rows, , drop = FALSE], 2L, effective / 1e6, "/")
    keep[rows] <- rowSums(cpm >= 1) >= minimum
  }
  indices <- which(keep)
  if (!length(indices)) stop("The condition-blind transcriptome filter retained no genes.")
  emit(paste("INFO | transcriptome | Retained genes:", length(indices),
    "; rule: TMM-adjusted CPM >= 1 in at least", minimum, "participants."))
  aligned <- match(as.integer(colnames(counts)), work$participant_id)
  if (anyNA(aligned)) stop("Transcriptome columns do not align to participant outcomes.")
  is_uc <- work$is_uc[aligned] == 1L
  n_uc <- sum(is_uc)
  n_cd <- sum(!is_uc)
  auc <- numeric(length(indices))
  for (start in seq.int(1L, length(indices), by = 1000L)) {
    positions <- seq.int(start, min(start + 999L, length(indices)))
    dge <- edgeR::DGEList(counts = counts[indices[positions], , drop = FALSE],
      lib.size = pooled$tmm$raw_library_size, norm.factors = pooled$tmm$tmm_norm_factor)
    expression <- edgeR::cpm(dge, log = TRUE, prior.count = 0.5)
    auc[positions] <- apply(expression, 1L, function(value) {
      (sum(rank(value, ties.method = "average")[is_uc]) - n_uc * (n_uc + 1) / 2) / (n_uc * n_cd)
    })
    emit(paste("INFO | transcriptome | Gene AUCs computed:", max(positions), "/", length(indices)))
  }
  deviation <- abs(auc - 0.5)
  probability <- c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 0.999, 1)
  distribution <- data.frame(quantile = c("min", "p05", "p25", "median", "p75", "p95", "p99", "p99.9", "max"),
    auc_uc_positive = as.numeric(stats::quantile(auc, probability)),
    absolute_deviation_from_0.5 = as.numeric(stats::quantile(deviation, probability)),
    stringsAsFactors = FALSE)
  all_genes <- data.frame(raw_gene_id = pooled$raw_gene_id[indices], gene_id = pooled$stable_gene_id[indices],
    auc_uc_positive = auc, absolute_deviation = deviation, stringsAsFactors = FALSE)
  # The archived top-100 ordering used deviations already rounded to six decimals.
  top <- all_genes[utils::head(order(-round(deviation, 6)), 100L),
    c("gene_id", "auc_uc_positive", "absolute_deviation"), drop = FALSE]
  top$status <- "CONTEXT_ONLY_NEVER_PROMOTED"
  universe <- data.frame(item = c("filter_rule", "genes_before_filter", "genes_after_filter", "uses_disease_label"),
    value = c("CPM >= 1 in at least 10 percent of participants", nrow(counts), sum(keep), "FALSE"),
    stringsAsFactors = FALSE)
  list(keep = keep, all_genes = all_genes, distribution = distribution, top = top,
       universe = universe, minimum_participants = minimum)
}

check_development_background <- function() {
  # Exercise the actual isolated helper before any large source file is opened.
  counts <- rbind(fixture_up = c(1, 2, 3, 4), fixture_down = c(4, 3, 2, 1),
                  fixture_tied = c(1, 1, 1, 1))
  colnames(counts) <- as.character(seq_len(ncol(counts)))
  library_size <- colSums(counts)
  fixture <- list(counts = counts, raw_gene_id = rownames(counts), stable_gene_id = rownames(counts),
    tmm = data.frame(raw_library_size = library_size, tmm_norm_factor = 1,
                     effective_library_size = library_size))
  outcomes <- data.frame(participant_id = seq_len(ncol(counts)), is_uc = c(0L, 0L, 1L, 1L))
  observed <- development_background(fixture, outcomes, function(...) NULL)
  if (!identical(as.numeric(observed$all_genes$auc_uc_positive), c(1, 0, 0.5)) ||
      !all(observed$keep) || nrow(observed$top) != 3L ||
      !setequal(observed$top$gene_id, rownames(counts))) {
    stop("Development background self-check failed for perfect, reversed, or tied scores.")
  }
  invisible(TRUE)
}

development_glm_table <- function(model, label) {
  if (!isTRUE(model$converged) || anyNA(stats::coef(model))) stop("A development association model did not converge or is rank deficient.")
  coefficients <- summary(model)$coefficients
  ci <- stats::confint.default(model, level = 0.95)
  data.frame(model_label = label, term = rownames(coefficients),
    estimate_log_odds_cd = coefficients[, 1L], standard_error = coefficients[, 2L],
    z_value = coefficients[, 3L], p_value = coefficients[, 4L],
    odds_ratio_cd = exp(coefficients[, 1L]), or_ci_low = exp(ci[, 1L]), or_ci_high = exp(ci[, 2L]),
    model_n = stats::nobs(model), stringsAsFactors = FALSE)
}

development_final_wald <- function(model, label) {
  if (!isTRUE(model$converged) || model$family$family != "binomial" || model$family$link != "logit" ||
      any(model$prior.weights != 1) || any(model$offset != 0)) {
    stop("Final-information intervals require an unweighted, converged binary logistic model without offsets.")
  }
  design <- stats::model.matrix(model)
  response <- model$y
  if (model$rank != ncol(design) || length(response) != nrow(design) ||
      anyNA(response) || any(!response %in% c(0, 1))) stop("Invalid final-information model design or response.")
  # Reuse the same complete cases and design; refine numerical convergence only.
  refined <- stats::glm.fit(x = design, y = response, start = stats::coef(model),
    family = stats::binomial(), control = stats::glm.control(epsilon = 1e-12, maxit = 100L))
  if (!isTRUE(refined$converged) || any(!is.finite(refined$coefficients))) {
    stop("The final-information logistic fit did not converge.")
  }
  coefficient <- as.numeric(refined$coefficients)
  probability <- stats::plogis(drop(design %*% coefficient))
  # Evaluate curvature at the final coefficients, not the preceding IWLS weights.
  information <- crossprod(design, design * as.numeric(probability * (1 - probability)))
  covariance <- chol2inv(chol(information))
  standard_error <- sqrt(diag(covariance))
  correction <- drop(covariance %*% crossprod(design, response - probability))
  if (any(!is.finite(standard_error)) || any(standard_error <= 0) ||
      any(!is.finite(correction)) || max(abs(correction)) > 1e-9) {
    stop("The final-information fit failed its coefficient convergence check.")
  }
  z <- coefficient / standard_error
  width <- stats::qnorm(0.975) * standard_error
  data.frame(model_label = label, term = colnames(design), estimate_log_odds_cd = coefficient,
    standard_error = standard_error, z_value = z, p_value = 2 * stats::pnorm(-abs(z)),
    odds_ratio_cd = exp(coefficient), or_ci_low = exp(coefficient - width),
    or_ci_high = exp(coefficient + width), model_n = nrow(design),
    maximum_final_newton_step = max(abs(correction)), stringsAsFactors = FALSE)
}

check_development_final_wald <- function() {
  fixture <- data.frame(x = rep(c(-1, 0, 1), each = 10L),
    y = c(rep(0, 7L), rep(1, 3L), rep(0, 5L), rep(1, 5L), rep(0, 3L), rep(1, 7L)))
  model <- stats::glm(y ~ x, data = fixture, family = stats::binomial())
  observed <- development_final_wald(model, "analytic_fixture")
  expected_coefficient <- c(0, log(7 / 3))
  expected_se <- sqrt(c(1 / 6.7, 1 / 4.2))
  expected_low <- exp(expected_coefficient - stats::qnorm(0.975) * expected_se)
  expected_high <- exp(expected_coefficient + stats::qnorm(0.975) * expected_se)
  if (!identical(observed$term, c("(Intercept)", "x")) ||
      max(abs(observed$estimate_log_odds_cd - expected_coefficient)) > 1e-10 ||
      max(abs(observed$standard_error - expected_se)) > 1e-10 ||
      max(abs(observed$or_ci_low - expected_low), abs(observed$or_ci_high - expected_high)) > 1e-10 ||
      any(observed$model_n != 30L)) stop("Final-information Wald self-check failed against an analytic logistic fixture.")
  invisible(TRUE)
}

fit_development_associations <- function(work, background) {
  specifications <- c(
    A2U_UNADJUSTED_REFERENCE = "y_cd ~ ighg3_z",
    A2_CONTRACTED_ADJUSTED = "y_cd ~ ighg3_z + site_inflammation + n_rectum_biopsies",
    A2B_DECLARED_EXTENDED_ADJUSTED = "y_cd ~ ighg3_z + site_inflammation + n_rectum_biopsies + sex + age",
    M3A_ADJUSTED_FOR_TOTAL_IG = "y_cd ~ ighg3_z + total_ig_z + site_inflammation",
    M3B_ADJUSTED_FOR_PLASMA_CELL_SCORE = "y_cd ~ ighg3_z + plasma_z + site_inflammation",
    M4_IGHG3_PROPORTION_OF_IG = "y_cd ~ ig_fraction_z + site_inflammation",
    M5A_GRADED_ENDOSCOPIC_SEVERITY = "y_cd ~ ighg3_z + endoscopic_severity_4 + sex + age",
    M5B_GRADED_SEVERITY_PLUS_CRP = "y_cd ~ ighg3_z + endoscopic_severity_4 + crp_z + sex + age")
  models <- do.call(rbind, lapply(names(specifications), function(label) {
    model <- stats::glm(stats::as.formula(specifications[[label]]), data = work,
      family = stats::binomial(), na.action = stats::na.omit)
    development_glm_table(model, label)
  }))
  features <- c("IGHG3", "IGHG3_FRACTION_OF_IG", "TOTAL_IG", "PLASMA_CELL_SCORE", "LYZ",
                "COL1A1", "IGG_IGA_LOG_RATIO", "IGG_IGAM_LOG_RATIO", "IGG_FRACTION_OF_IG")
  feature_models <- feature_models_final <- list()
  context <- do.call(rbind, lapply(features, function(feature) {
    data <- work
    data$feature_z <- as.numeric(scale(data[[feature]]))
    crude <- stats::glm(y_cd ~ feature_z, data = data, family = stats::binomial(), na.action = stats::na.omit)
    adjusted <- stats::glm(y_cd ~ feature_z + endoscopic_severity_4 + sex + age,
      data = data, family = stats::binomial(), na.action = stats::na.omit)
    crude_table <- development_glm_table(crude, feature)
    adjusted_table <- development_glm_table(adjusted, feature)
    feature_models[[feature]] <<- adjusted_table
    final_table <- development_final_wald(adjusted, feature)
    feature_models_final[[feature]] <<- final_table
    crude_term <- crude_table[crude_table$term == "feature_z", , drop = FALSE]
    adjusted_term <- adjusted_table[adjusted_table$term == "feature_z", , drop = FALSE]
    final_term <- final_table[final_table$term == "feature_z", , drop = FALSE]
    auc <- uc_auc_summary(data[[feature]], data$is_uc)
    deviation <- abs(background$all_genes$auc_uc_positive - 0.5)
    data.frame(feature = feature, auc_uc_positive = auc$auc, auc_ci_low = auc$ci95_low,
      auc_ci_high = auc$ci95_high, transcriptome_percentile = 100 * mean(deviation <= abs(auc$auc - 0.5)),
      genes_more_extreme = sum(deviation > abs(auc$auc - 0.5)),
      crude_log_odds_cd = crude_term$estimate_log_odds_cd, crude_p = crude_term$p_value,
      adjusted_log_odds_cd = adjusted_term$estimate_log_odds_cd, adjusted_p = adjusted_term$p_value,
      adjusted_odds_ratio_cd = adjusted_term$odds_ratio_cd, model_n = adjusted_term$model_n,
      adjusted_ci_low = adjusted_term$or_ci_low, adjusted_ci_high = adjusted_term$or_ci_high,
      final_odds_ratio_cd = final_term$odds_ratio_cd,
      final_ci_low = final_term$or_ci_low, final_ci_high = final_term$or_ci_high,
      stringsAsFactors = FALSE)
  }))
  context <- context[order(-abs(context$auc_uc_positive - 0.5)), , drop = FALSE]
  strata <- do.call(rbind, lapply(levels(work$site_inflammation), function(stratum) {
    data <- work[work$site_inflammation == stratum, , drop = FALSE]
    auc <- uc_auc_summary(data$IGHG3, data$is_uc)
    data.frame(auc = auc$auc, ci_low = auc$ci95_low, ci_high = auc$ci95_high,
      n_uc = auc$n_uc, n_cd = auc$n_cd, ci_method = "DELONG", stratum = stratum,
      endpoint_id = "A3", stringsAsFactors = FALSE)
  }))
  list(models = models, model_specifications = specifications, feature_context = context,
       feature_models = do.call(rbind, feature_models),
       feature_models_final = do.call(rbind, feature_models_final), strata = strata)
}

format_interval <- function(low, high) {
  paste(formatC(low, format = "f", digits = 3L, decimal.mark = "."), "to",
        formatC(high, format = "f", digits = 3L, decimal.mark = "."))
}

development_supplementary_tables <- function(panel, associations, comparisons, n_participants) {
  table_s1 <- data.frame(gene_symbol = panel$gene_symbol, ensembl_id = panel$gene_id, role = panel$role,
    stringsAsFactors = FALSE)
  context <- associations$feature_context
  table_s3a <- data.frame(feature = context$feature, auc = context$auc_uc_positive,
    interval = format_interval(context$auc_ci_low, context$auc_ci_high),
    transcriptome_percentile = context$transcriptome_percentile, genes_more_extreme = context$genes_more_extreme,
    adjusted_odds_ratio = context$final_odds_ratio_cd,
    adjusted_OR_ci_low = context$final_ci_low, adjusted_OR_ci_high = context$final_ci_high,
    n_analyzed = context$model_n, adjusted_OR_outcome = "CD=1; UC=0; per one SD of feature",
    adjusted_OR_model = "y_cd ~ feature_z + endoscopic_severity_4 + sex + age",
    adjusted_OR_interval_type = "Nominal 95% Wald; no multiplicity adjustment", stringsAsFactors = FALSE)
  auc_all <- comparisons$aucs[comparisons$aucs$feature == "IGHG3_abundance", , drop = FALSE]
  auc_row <- function(endpoint, value, low, high, n, note) {
    data.frame(endpoint = endpoint, estimate = value, interval = format_interval(low, high),
      scale = "ROC-AUC, ulcerative colitis positive", model_formula = "Not applicable; direct score AUC",
      interval_type = "Pointwise 95% DeLong", n_analyzed = n, population_note = note,
      stringsAsFactors = FALSE)
  }
  rows <- list(auc_row("A1 discrimination, corrected orientation", auc_all$auc,
    auc_all$ci95_low, auc_all$ci95_high, auc_all$n_uc + auc_all$n_cd,
    paste("All", n_participants, "development participants")))
  labels <- c(
    A2U_UNADJUSTED_REFERENCE = "A2U association, unadjusted",
    A2_CONTRACTED_ADJUSTED = "A2 association, adjusted for site inflammation",
    A2B_DECLARED_EXTENDED_ADJUSTED = "A2B association, plus sex and age",
    M5A_GRADED_ENDOSCOPIC_SEVERITY = "M5A association, graded endoscopic severity",
    M5B_GRADED_SEVERITY_PLUS_CRP = "M5B association, graded severity plus CRP",
    M3A_ADJUSTED_FOR_TOTAL_IG = "M3A association, conditional on total immunoglobulin")
  for (label in names(labels)) {
    term <- associations$models[associations$models$model_label == label &
      associations$models$term == "ighg3_z", , drop = FALSE]
    if (nrow(term) != 1L) stop("A supplementary IGHG3 endpoint has no unique model term.")
    note <- if (label == "M5B_GRADED_SEVERITY_PLUS_CRP") {
      paste(n_participants - term$model_n, "excluded for missing CRP; complete cases")
    } else paste("All", n_participants, "development participants")
    rows[[length(rows) + 1L]] <- data.frame(endpoint = labels[[label]], estimate = term$odds_ratio_cd,
      interval = format_interval(term$or_ci_low, term$or_ci_high),
      scale = "odds ratio for Crohn disease per standard deviation",
      model_formula = associations$model_specifications[[label]],
      interval_type = "Nominal 95% Wald; no multiplicity adjustment", n_analyzed = term$model_n,
      population_note = note, stringsAsFactors = FALSE)
  }
  for (stratum in c("ALL_SITE_NONINFLAMED", "ANY_SITE_INFLAMED")) {
    auc <- associations$strata[associations$strata$stratum == stratum, , drop = FALSE]
    label <- if (stratum == "ALL_SITE_NONINFLAMED") "A3 discrimination, non inflamed stratum" else
      "A3 discrimination, inflamed stratum"
    rows[[length(rows) + 1L]] <- auc_row(label, auc$auc, auc$ci_low, auc$ci_high,
      auc$n_uc + auc$n_cd, "Inflammation stratum")
  }
  list(table_s1 = table_s1, table_s3a = table_s3a, table_s3b = do.call(rbind, rows))
}

development_archive_view <- function(table, decimal = numeric(), significant = character()) {
  for (field in names(decimal)) table[[field]] <- round(table[[field]], decimal[[field]])
  for (field in significant) table[[field]] <- signif(table[[field]], 6L)
  table
}

validate_development <- function(bundle, project, record) {
  archive <- file.path(project, "reference", "development_expected.zip")
  if (!identical(unname(tools::md5sum(archive)), "1a4e2e88f4139bea9fbd0844329fc095")) {
    stop("The archived development comparison tables failed their checksum.")
  }
  expected_files <- c("18_gse193677_sample_selection.tsv", "18_gse193677_participant_geometry.tsv",
    "19_source_token_mapping_audit.tsv", "19_participant_aggregation_audit.tsv",
    "19_tmm_normalization_qc.tsv", "19_gse193677_IGHG3_participant_feature.tsv",
    "19_gse193677_numeric_audit.tsv", "27_biopsy_level_covariates.tsv",
    "27_A2_adjusted_association.tsv", "27_A3_stratified_discrimination.tsv",
    "28_panel_gene_availability.tsv", "28_M3_M5_conditional_models.tsv",
    "29_analysis_table_with_composites.tsv", "29_feature_percentile_context.tsv",
    "29_transcriptome_wide_auc_distribution.tsv", "29_top_genes_by_auc_context_only.tsv",
    "Table_S15_expanded_paired_auc.tsv", "Table_S16_comparator_auc.tsv", "29_gene_universe_definition.tsv",
    "independent_association_checks.tsv")
  listing <- utils::unzip(archive, list = TRUE)$Name
  if (anyDuplicated(listing) || !setequal(listing, expected_files)) {
    stop("The development comparison archive has unexpected members.")
  }
  directory <- tempfile(pattern = "ibd_development_expected_")
  if (!dir.create(directory)) stop("Cannot create the temporary comparison directory.")
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(archive, files = expected_files, exdir = directory, unzip = "internal")
  read <- function(name) utils::read.delim(file.path(directory, name),
    check.names = FALSE, stringsAsFactors = FALSE)
  outcomes <- logical()
  compare <- function(observed, filename, keys, fields = NULL, exact = character(), expected_filter = NULL) {
    expected <- read(filename)
    if (!is.null(expected_filter)) expected <- expected_filter(expected)
    if (is.null(fields)) fields <- names(expected)
    outcomes <<- c(outcomes, compare_table(observed, expected, keys, filename, record,
      fields = fields, exact_fields = exact))
  }
  compare(bundle$selection, expected_files[[1L]], "sample_id", fields = names(bundle$selection),
    exact = "participant_id_candidate")
  compare(bundle$geometry, expected_files[[2L]], "diagnosis",
    exact = c("n_biopsies", "n_participants", "repeated_biopsy_rows"))
  compare(bundle$pooled$mapping, expected_files[[3L]], "sample_id", exact = "raw_data_field_position")
  compare(bundle$pooled$aggregation, expected_files[[4L]], "participant_id", exact = "n_biopsies_aggregated")
  compare(bundle$pooled$tmm, expected_files[[5L]], "participant_id", exact = "raw_library_size")
  feature_fields <- setdiff(names(read(expected_files[[6L]])), "expression_values_read_in_stage19")
  compare(bundle$pooled$forward_feature, expected_files[[6L]], "participant_id", fields = feature_fields,
    exact = c("raw_count_sum", "library_size_full_native"))
  compare(bundle$pooled$numeric_audit, expected_files[[7L]], "metric", expected_filter = function(table) {
    table[table$metric != "model_fitted_in_stage19", , drop = FALSE]
  })
  compare(bundle$biopsy_covariates, expected_files[[8L]], "sample_id", exact = "participant_id")
  models <- bundle$associations$models
  a2 <- models[startsWith(models$model_label, "A2"), , drop = FALSE]
  rounded_model <- function(table) development_archive_view(table,
    decimal = c(estimate_log_odds_cd = 6L, standard_error = 6L, z_value = 6L,
      odds_ratio_cd = 6L, or_ci_low = 6L, or_ci_high = 6L), significant = "p_value")
  compare(rounded_model(a2), expected_files[[9L]], c("model_label", "term"), exact = "model_n")
  strata <- development_archive_view(bundle$associations$strata, decimal = c(auc = 6L, ci_low = 6L, ci_high = 6L))
  compare(strata, expected_files[[10L]], "stratum", exact = c("n_uc", "n_cd"))
  compare(bundle$panel_availability, expected_files[[11L]], "gene_id")
  conditional <- rounded_model(models[!startsWith(models$model_label, "A2") & models$term != "(Intercept)", , drop = FALSE])
  names(conditional)[names(conditional) == "estimate_log_odds_cd"] <- "log_odds_cd"
  compare(conditional, expected_files[[12L]], c("model_label", "term"), exact = "model_n")
  compare(bundle$analysis, expected_files[[13L]], "participant_id",
    exact = c("n_rectum_biopsies", "y_cd", "is_uc", "endo_rank_max"))
  context <- development_archive_view(bundle$associations$feature_context,
    decimal = c(auc_uc_positive = 6L, auc_ci_low = 6L, auc_ci_high = 6L, transcriptome_percentile = 4L,
      crude_log_odds_cd = 6L, adjusted_log_odds_cd = 6L, adjusted_odds_ratio_cd = 6L),
    significant = c("crude_p", "adjusted_p"))
  compare(context, expected_files[[14L]], "feature", exact = c("genes_more_extreme", "model_n"))
  compare(development_archive_view(bundle$background$distribution,
    decimal = c(auc_uc_positive = 6L, absolute_deviation_from_0.5 = 6L)), expected_files[[15L]], "quantile")
  compare(development_archive_view(bundle$background$top,
    decimal = c(auc_uc_positive = 6L, absolute_deviation = 6L)), expected_files[[16L]], "gene_id")
  development_only <- function(table) table[table$cohort == "GSE193677", , drop = FALSE]
  compare(bundle$comparisons$paired, expected_files[[17L]], c("cohort", "comparator"),
    exact = c("n_uc", "n_cd"), expected_filter = development_only)
  compare(bundle$comparisons$aucs, expected_files[[18L]], c("cohort", "feature"),
    exact = c("n_uc", "n_cd"), expected_filter = development_only)
  compare(bundle$background$universe, expected_files[[19L]], "item")
  final_models <- bundle$associations$feature_models_final
  final_models <- final_models[final_models$term == "feature_z", , drop = FALSE]
  final_comparison <- data.frame(feature = final_models$model_label, n = final_models$model_n,
    OR = final_models$odds_ratio_cd, ci_low = final_models$or_ci_low, ci_high = final_models$or_ci_high,
    coefficient = final_models$estimate_log_odds_cd, se = final_models$standard_error)
  compare(final_comparison, expected_files[[20L]], "feature", fields = names(final_comparison), exact = "n",
    expected_filter = function(table) table[table$model == "Stage29_adjusted", , drop = FALSE])
  final_zip <- file.path(project, "reference", "tables_reference.zip")
  map <- utils::read.delim(file.path(project, "config", "tables.tsv"),
    colClasses = "character", check.names = FALSE)
  reference_md5 <- unique(map$reference_zip_md5)
  if (length(reference_md5) != 1L || !identical(unname(tools::md5sum(final_zip)), reference_md5)) {
    stop("The supplementary comparison archive failed its checksum.")
  }
  final_names <- c("Table_S1_gene_panel.tsv", "Table_S3A_development_comparisons.tsv",
                   "Table_S3B_IGHG3_development_endpoints.tsv")
  utils::unzip(final_zip, files = final_names, exdir = directory, unzip = "internal")
  outcomes <- c(outcomes, compare_table(bundle$supplementary$table_s1, read(final_names[[1L]]),
    "ensembl_id", "current_Table_S1_complete", record))
  # Published numeric rounding is applied only to the comparison view.
  s3a <- development_archive_view(bundle$supplementary$table_s3a,
    decimal = c(auc = 3L, transcriptome_percentile = 1L, adjusted_odds_ratio = 3L))
  outcomes <- c(outcomes, compare_table(s3a, read(final_names[[2L]]), "feature", "current_Table_S3A_complete",
    record, exact_fields = c("genes_more_extreme", "n_analyzed")))
  s3b <- development_archive_view(bundle$supplementary$table_s3b, decimal = c(estimate = 3L))
  outcomes <- c(outcomes, compare_table(s3b, read(final_names[[3L]]), "endpoint", "current_Table_S3B_complete",
    record, exact_fields = "n_analyzed"))
  all(outcomes)
}

run_development <- function(context) {
  project <- context$project
  work <- file.path(project, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  report_path <- file.path(work, "development_report.txt")
  connection <- file(report_path, open = "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  counters <- c(PASS = 0L, WARN = 0L, FAIL = 0L)
  emit <- function(lines) {
    lines <- as.character(lines)
    writeLines(lines, connection)
    flush(connection)
    cat(paste(lines, collapse = "\n"), "\n", sep = "")
  }
  record <- function(status, name, detail) {
    counters[[status]] <<- counters[[status]] + 1L
    emit(paste(status, name, detail, sep = " | "))
  }
  status <- 1L
  tryCatch(withCallingHandlers({
    emit(c("IBD bulk development source-to-result reproduction report", "Code revision: stage04.2",
      paste("Started UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      "Analysis inputs: the native estimated-count matrix and declared source sample metadata for GSE193677.",
      "The validated newly computed replication cache supplies the forward feature and source gene annotation.",
      "Archived output tables are opened only for comparison after computation.",
      "Development score exploration and later comparator additions retain their retrospective status."))
    check_validation_logic()
    record("PASS", "validator_self_checks", "Exact identifiers/counts and numerical, probability, missingness, and duplicate-key checks.")
    required <- c(edgeR = "4.4.2", limma = "3.62.2", Matrix = "1.7-1", data.table = "1.18.4", pROC = "1.19.0.1")
    for (package in names(required)) {
      if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package)
      observed <- utils::packageVersion(package)
      if (observed != package_version(required[[package]])) stop("Package version differs for ", package,
        ": observed ", as.character(observed), "; expected ", required[[package]])
      record("PASS", paste0("package_", package), as.character(observed))
    }
    check_development_background()
    record("PASS", "background_helper_self_checks",
      "Isolated helper executes with perfect, reversed, and tied synthetic scores; AUCs are 1, 0, and 0.5.")
    check_development_final_wald()
    record("PASS", "final_information_wald_self_checks",
      "Both coefficients, standard errors, and interval limits agree with an analytic binary logistic fixture.")
    if (getRversion() != package_version("4.4.2")) record("WARN", "R_version",
      paste("Observed", as.character(getRversion()), "; archived analyses used 4.4.2."))
    if (!isTRUE(context$settings$verify_md5)) stop("Development requires checksum verification.")
    manifest <- context$read_manifest(file.path(project, "config", "inputs.tsv"))
    cache_path <- file.path(work, "replication.rds")
    if (!file.exists(cache_path)) stop("Run and validate replication first; work/replication.rds is required.")
    replication <- readRDS(cache_path)
    if (!isTRUE(replication$validation_passed) || !identical(replication$code_revision, "stage03.0") ||
        !nrow(replication$primary_freeze) || any(replication$primary_freeze$feature_status != "FROZEN_PRIMARY_BULK_FEATURE")) {
      stop("The replication cache is incomplete, incompatible, or not validated.")
    }
    for (identifier in c("replication_archive", "replication_annotations", "gse214695_sample_metadata")) {
      spec <- manifest[manifest$input_id == identifier, , drop = FALSE]
      if (nrow(spec) != 1L || !identical(unname(replication$source_md5[identifier]), spec$expected_md5[[1L]])) {
        stop("The replication cache source identities disagree with the input definitions.")
      }
    }
    forward <- replication$primary_freeze
    if (any(forward$orientation_rule != "positive_logFC_means_UC_gt_CD")) stop("The forward-gene orientation changed.")
    replication_md5 <- unname(tools::md5sum(cache_path))
    record("PASS", "replication_cache", paste("Validated stage03.0; forward genes:",
      paste(forward$gene_symbol, collapse = ", "), "; cache MD5:", replication_md5))
    panel <- development_gene_panel()
    source_features <- replication$pseudobulk$features
    source_index <- match(panel$gene_id, source_features$gene_id)
    if (anyNA(source_index) || anyDuplicated(source_features$gene_id)) {
      stop("The fixed panel lacks unique stable-Ensembl support in the verified replication annotation.")
    }
    source_symbols <- source_features$gene_symbol[source_index]
    if (identical(source_symbols, panel$gene_symbol)) {
      record("PASS", "panel_source_annotation", "All 25 fixed panel gene identifiers and symbols match the verified source feature annotation.")
    } else {
      difference <- which(source_symbols != panel$gene_symbol)
      record("WARN", "panel_source_symbols", paste("Stable identifiers match; source symbol differences:",
        paste(paste0(panel$gene_symbol[difference], "=", source_symbols[difference]), collapse = "; ")))
    }
    rm(replication, source_features)
    paths <- source_md5 <- character()
    for (identifier in c("development_counts", "gse193677_sample_metadata")) {
      spec <- manifest[manifest$input_id == identifier, , drop = FALSE]
      if (nrow(spec) != 1L || !context$safe_relative(spec$relative_path[[1L]]) || is.na(spec$expected_md5[[1L]])) {
        stop("Missing or invalid input definition: ", identifier)
      }
      path <- file.path(context$settings$legacy_root, spec$relative_path[[1L]])
      if (identifier %in% names(context$settings$file_overrides)) path <- context$settings$file_overrides[[identifier]]
      path <- normalizePath(path.expand(path), winslash = "/", mustWork = TRUE)
      observed <- unname(tools::md5sum(path))
      if (file.info(path)$size != as.numeric(spec$expected_bytes[[1L]]) || !identical(observed, spec$expected_md5[[1L]])) {
        stop("Source byte count or checksum differs for ", identifier, ".")
      }
      paths[identifier] <- path
      source_md5[identifier] <- observed
      record("PASS", identifier, paste("Archived source MD5:", observed))
    }
    selected <- select_development_samples(paths[["gse193677_sample_metadata"]], context$read_manifest)
    record("PASS", "source_selection", paste("Selected biopsies:", nrow(selected$selection),
      "; participants:", length(unique(selected$selection$participant_id_candidate))))
    pooled <- import_development_counts(paths[["development_counts"]], selected$selection, forward, emit)
    record("PASS", "estimated_count_aggregation", paste("Native gene rows:", nrow(pooled$counts),
      "; pooled participants:", ncol(pooled$counts), "; repeated-biopsy participants:",
      sum(pooled$aggregation$n_biopsies_aggregated > 1L), "; full-universe TMM; source fractions preserved."))
    covariates <- development_covariates(selected)
    record("PASS", "site_inflammation_provenance", "Every selected typere value agrees with its source-title site token; clinical activity is retained separately.")
    analysis <- construct_development_features(pooled, covariates, panel)
    record("PASS", "panel_and_scores", "25 uniquely mapped genes; full-library abundance transforms and native-count isotype ratios constructed.")
    emit(paste("INFO | participants | UC:", sum(analysis$is_uc), "; CD:", sum(analysis$y_cd),
      "; noninflamed site:", sum(analysis$site_inflammation == "ALL_SITE_NONINFLAMED"),
      "; inflamed site:", sum(analysis$site_inflammation == "ANY_SITE_INFLAMED"),
      "; missing CRP:", sum(is.na(analysis$crp_log2))))
    scores <- data.frame(primary_ratio = analysis$IGG_IGAM_LOG_RATIO, IGHG3_abundance = analysis$IGHG3,
      total_Ig_abundance = analysis$TOTAL_IG, IGHG3_fraction = analysis$IGHG3_FRACTION_OF_IG,
      summed_IgG_abundance = analysis$IGG_ABUNDANCE)
    comparisons <- cohort_score_comparisons(scores, analysis$is_uc, "GSE193677", "exploratory development")
    record("PASS", "delong_implementation", "All five AUC intervals agree with independent tie-aware placement calculations; four explicitly paired comparisons completed.")
    background <- development_background(pooled, analysis, emit)
    associations <- fit_development_associations(analysis, background)
    supplementary <- development_supplementary_tables(panel, associations, comparisons, nrow(analysis))
    availability <- panel[, c("gene_symbol", "gene_id", "gene_class"), drop = FALSE]
    availability$present <- panel$gene_id %in% pooled$stable_gene_id
    bundle <- list(code_revision = "stage04.2", source_md5 = source_md5,
      replication_cache_md5 = replication_md5, created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      selection = selected$selection, geometry = selected$geometry, pooled = pooled,
      biopsy_covariates = covariates$biopsy, analysis = analysis, panel_availability = availability,
      background = background, associations = associations, comparisons = comparisons, supplementary = supplementary,
      settings = list(unit = "participant", cohort = "GSE193677", anatomy = "Rectum",
        count_unit = "native estimated counts; fractions retained", normalization = "TMM on all 56632 native rows",
        log2cpm_prior_count = 0.5, ratio_pseudocount = 0.5,
        abundance_sum_formula = "log2(1e6 * (X + 0.5) / (L_eff + 1))",
        ratio_formula = "log2((G + 0.5) / (A + M + 0.5))",
        auc_positive = "UC", auc_direction = "<", association_positive = "CD",
        association_intervals = "nominal 95 percent Wald via confint.default",
        table_s3a_intervals = "same nine logistic designs; refined convergence; final-coefficient Fisher information; nominal 95 percent Wald",
        background = "TMM-adjusted CPM >= 1 in at least ceiling(0.10*n) participants; observed-label ranks",
        numerical_atol = 1e-10, probability_atol = 1e-14, numerical_rtol = 1e-6), session = utils::sessionInfo())
    emit(c("", "Comparison with archived participant values, models, and current supplementary tables"))
    passed <- validate_development(bundle, project, record)
    bundle$validation_passed <- isTRUE(passed) && counters[["FAIL"]] == 0L
    bundle$validation_counts <- counters
    cache <- file.path(work, "development.rds")
    emit("INFO | cache | Compressing one computed bundle, including the full native participant-count matrix.")
    saveRDS(bundle, cache, compress = "xz")
    emit(c("", "Newly computed development discrimination", utils::capture.output(print(comparisons$aucs,
      row.names = FALSE, digits = 10)), "", "Newly computed development feature context",
      utils::capture.output(print(associations$feature_context[, c("feature", "auc_uc_positive",
        "transcriptome_percentile", "genes_more_extreme", "adjusted_odds_ratio_cd")], row.names = FALSE, digits = 10)),
      paste("Background genes:", nrow(background$all_genes)), paste("Cache:", normalizePath(cache, winslash = "/", mustWork = TRUE)),
      "The cache retains full calculated precision. Rounding is limited to reference views and display intervals.",
      "Later stages must check validation_passed; a fatal run does not refresh a previous cache.",
      "", "R session information", utils::capture.output(utils::sessionInfo())))
    status <- if (isTRUE(bundle$validation_passed)) 0L else 1L
  }, warning = function(warning) {
    record("WARN", "R_warning", conditionMessage(warning))
    invokeRestart("muffleWarning")
  }), error = function(error) {
    record("FAIL", "development_error", conditionMessage(error))
    status <<- 1L
  })
  emit(c("", paste("Checks: PASS =", counters[["PASS"]], "; WARN =", counters[["WARN"]], "; FAIL =", counters[["FAIL"]]),
    paste("Development reproduction:", if (status == 0L) "PASS" else "BLOCKED"),
    "External/context source reconstruction, remaining sensitivity analyses, final workbook, and figures: NOT_RUN",
    paste("Finished UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Report:", normalizePath(report_path, winslash = "/", mustWork = TRUE))))
  status
}
