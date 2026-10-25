# Reconstruct external TPM and HMP2 context results from declared source files.
# Shared validators and bulk statistical helpers are loaded by run.R.

external_clean <- function(value) gsub('^"|"$', "", trimws(as.character(value)))

external_unique_key <- function(value, label) {
  if (anyNA(value) || any(!nzchar(value)) || anyDuplicated(value)) {
    stop("Missing or duplicate identifiers: ", label, ".")
  }
  invisible(TRUE)
}

external_ig_panel <- function() {
  panel <- development_gene_panel()
  panel$gene_symbol[panel$gene_class == "IG_CONSTANT"]
}

external_components <- function(block) {
  panel <- external_ig_panel()
  if (!is.matrix(block) || any(!is.finite(block)) || any(block < 0) ||
      anyDuplicated(rownames(block)) || !setequal(rownames(block), panel)) {
    stop("The external Ig block must contain finite, nonnegative values for all 12 unique panel genes.")
  }
  data.frame(igg_sum = colSums(block[c("IGHG1", "IGHG2", "IGHG3", "IGHG4"), , drop = FALSE]),
    iga_sum = colSums(block[c("IGHA1", "IGHA2"), , drop = FALSE]),
    igm_sum = as.numeric(block["IGHM", ]),
    igam_sum = colSums(block[c("IGHA1", "IGHA2", "IGHM"), , drop = FALSE]),
    total_ig_sum = colSums(block), ighg3_count = as.numeric(block["IGHG3", ]),
    row.names = NULL, stringsAsFactors = FALSE)
}

external_score_table <- function(components, unit, library_size = NULL) {
  if (!unit %in% c("TPM", "COUNTS")) stop("Unknown external expression unit.")
  if (unit == "COUNTS") {
    if (length(library_size) != nrow(components) || any(!is.finite(library_size)) ||
        any(library_size <= 0) || any(components$total_ig_sum > library_size)) {
      stop("Pooled context libraries must be finite, positive, and contain the Ig counts.")
    }
    multiplier <- 1e6 / library_size
    ighg3_abundance <- 1e6 * components$ighg3_count / library_size
    total_ig_abundance <- 1e6 * components$total_ig_sum / library_size
    igg_abundance <- 1e6 * components$igg_sum / library_size
  } else {
    multiplier <- rep(1, nrow(components))
    ighg3_abundance <- components$ighg3_count
    total_ig_abundance <- components$total_ig_sum
    igg_abundance <- components$igg_sum
  }
  result <- components
  result$IGHG3_value <- components$ighg3_count
  result$TOTAL_IG_value <- components$total_ig_sum
  result$IGG_IGAM_LOG_RATIO <- log2((components$igg_sum + 0.5) / (components$igam_sum + 0.5))
  result$IGHG3_FRACTION_OF_IG <- log2((components$ighg3_count + 0.5) /
    (components$total_ig_sum - components$ighg3_count + 0.5))
  result$IGHG3_LOG <- log2(components$ighg3_count + 0.5)
  result$TOTAL_IG_LOG <- log2(components$total_ig_sum + 0.5)
  result$IGHG3_ABUNDANCE <- log2(ighg3_abundance + 0.5)
  result$TOTAL_IG_ABUNDANCE <- log2(total_ig_abundance + 0.5)
  result$IGG_ABUNDANCE <- log2(igg_abundance + 0.5)
  result$normalization_multiplier <- multiplier
  if (any(!is.finite(as.matrix(result)))) stop("An external score or component is nonfinite.")
  result
}

external_unit_audit <- function(genes, values, study, member) {
  required <- external_ig_panel()[seq_len(7L)]
  total <- sum(values)
  integer_share <- mean(abs(values - round(values)) < 1e-8)
  unit <- if (integer_share > 0.999) "RAW_COUNTS" else if (abs(total - 1e6) / 1e6 < 0.02) {
    "PER_MILLION_TPM_OR_CPM"
  } else if (abs(total - 100) / 100 < 0.02) "PERCENTAGE" else "LENGTH_NORMALISED_OR_OTHER"
  positive <- values[values > 0]
  data.frame(study_id = study, member_name = member, n_rows = length(genes), n_finite = length(values),
    column_sum = total, integer_share = integer_share,
    min_positive = if (length(positive)) min(positive) else NA_real_, max_value = max(values),
    inferred_unit = unit, panel_genes_found = sum(required %in% genes),
    duplicated_panel_rows = sum(duplicated(genes[genes %in% required])),
    igg_sum = sum(values[genes %in% required[seq_len(4L)]]),
    igam_sum = sum(values[genes %in% required[5L:7L]]), stringsAsFactors = FALSE)
}

external_archive_import <- function(path, study, expected_members, safe_relative, emit) {
  members <- utils::untar(path, list = TRUE, tar = "internal")
  if (length(members) != expected_members || anyDuplicated(members) ||
      !all(vapply(members, safe_relative, logical(1))) ||
      any(!grepl("^GSM[0-9]+_[^/\\\\]+\\.txt\\.gz$", members))) {
    stop("Unexpected or unsafe expression archive members for ", study, ".")
  }
  samples <- sub("_.*$", "", members)
  external_unique_key(samples, paste(study, "archive GSMs"))
  directory <- tempfile(pattern = paste0("ibd_", study, "_"))
  if (!dir.create(directory)) stop("Cannot create temporary expression storage.")
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  emit(paste("INFO | archive_extraction |", study, ":", length(members), "compressed members; temporary storage only."))
  status <- utils::untar(path, files = members, exdir = directory, tar = "internal")
  if (!identical(as.integer(status), 0L) || !all(file.exists(file.path(directory, members)))) {
    stop("Expression archive extraction failed for ", study, ".")
  }
  panel <- external_ig_panel()
  block <- matrix(NA_real_, nrow = length(panel), ncol = length(members), dimnames = list(panel, samples))
  audits <- vector("list", length(members))
  feature_order <- NULL
  for (i in seq_along(members)) {
    member_path <- file.path(directory, members[[i]])
    connection <- gzfile(member_path, open = "rt")
    table <- tryCatch(utils::read.delim(connection, sep = "\t", header = TRUE,
      colClasses = "character", quote = "", comment.char = "", check.names = FALSE, na.strings = character(),
      stringsAsFactors = FALSE), finally = close(connection))
    if (ncol(table) != 2L) stop("Expected exactly two columns in ", members[[i]], ".")
    genes <- external_clean(table[[1L]])
    values <- as.numeric(external_clean(table[[2L]]))
    if (any(!is.finite(values)) || any(values < 0)) {
      stop("Invalid source expression values in ", members[[i]], ".")
    }
    expected_rows <- if (study == "GSE117993") 58524L else 23228L
    if (length(genes) != expected_rows) stop("Unexpected native feature count in ", members[[i]], ".")
    if (is.null(feature_order)) feature_order <- genes
    audits[[i]] <- external_unit_audit(genes, values, study, members[[i]])
    if (study == "GSE117993") {
      if (!all(panel %in% genes) || abs(sum(values) - 1e6) / 1e6 >= 0.02) {
        stop("The external file lacks the complete panel or its expected per-million scale.")
      }
      # Preserve the original sum rule if a panel symbol has multiple native rows.
      block[, i] <- vapply(panel, function(symbol) sum(values[genes %in% symbol]), numeric(1))
    } else if (any(panel[seq_len(7L)] %in% genes)) {
      stop("The archived missing-component disposition of GSE83687 requires review.")
    }
    unlink(member_path)
    if (i %% 25L == 0L || i == length(members)) {
      emit(paste("INFO | expression_import |", study, ":", i, "/", length(members)))
    }
  }
  list(panel = if (study == "GSE117993") block else NULL,
    member_index = data.frame(study_id = study, member_name = members, extension = "txt",
      looks_like_gsm = TRUE, stringsAsFactors = FALSE),
    member_audit = do.call(rbind, audits),
    unit_probes = do.call(rbind, audits[round(seq(1, length(members), length.out = 5L))]),
    source_features = feature_order,
    disposition = if (study == "GSE117993") "ESTIMABLE_TPM" else "NON_ESTIMABLE_MISSING_CONSTANT_GENES")
}

external_characteristics <- function(metadata) {
  parts <- strsplit(external_clean(metadata$characteristics_raw), "\\s*\\|\\s*")
  get <- function(key) vapply(parts, function(row) {
    hit <- row[startsWith(row, paste0(key, ":"))]
    values <- unique(trimws(sub("^[^:]*:", "", hit)))
    if (length(values) != 1L || is.na(values) || !nzchar(values)) {
      stop("Missing or conflicting external source characteristic: ", key, ".")
    }
    values[[1L]]
  }, character(1))
  data.frame(sample_id = metadata$sample_id, tissue = get("tissue"), diagnosis_characteristic = get("diagnosis"),
    montreal_location = get("location"), disease_subtype = get("disease subtype"),
    population = get("population"), biopsy_location = get("biopsy location"), stringsAsFactors = FALSE)
}

external_metadata_join <- function(imported, path, read_manifest) {
  metadata <- read_manifest(path)
  required <- c("sample_id", "study_id", "title", "diagnosis_raw", "anatomy_raw", "characteristics_raw",
    "inflammation_raw", "treatment_raw")
  if (!all(required %in% names(metadata)) || nrow(metadata) != 190L ||
      anyNA(metadata$study_id) || any(metadata$study_id != "GSE117993")) stop("Unexpected external source metadata.")
  external_unique_key(metadata$sample_id, "external source GSMs")
  if (!setequal(metadata$sample_id, colnames(imported$panel))) stop("External expression and metadata GSMs disagree.")
  metadata <- metadata[match(sort(metadata$sample_id), metadata$sample_id), , drop = FALSE]
  characteristic <- external_characteristics(metadata)
  diagnosis <- toupper(external_clean(metadata$diagnosis_raw))
  if (anyNA(diagnosis) || any(!diagnosis %in% c("UC", "CD", "CONTROL")) ||
      any(toupper(characteristic$diagnosis_characteristic) != diagnosis) ||
      any(toupper(metadata$anatomy_raw) != "RECTUM") || any(toupper(characteristic$tissue) != "RECTUM") ||
      any(toupper(characteristic$biopsy_location) != "RECTAL")) stop("External diagnosis or rectal-site annotations disagree.")
  source_token <- sub("_RISK_.*$", "", metadata$title)
  if (any(!grepl("^[0-9]+$", source_token))) stop("Cannot recover the external source title tokens.")
  external_unique_key(source_token, "external source title tokens")
  block <- imported$panel[, metadata$sample_id, drop = FALSE]
  scores <- external_score_table(external_components(block), "TPM")
  out <- data.frame(sample_id = metadata$sample_id, scores,
    diagnosis_source = metadata$diagnosis_raw, inflammation_source = metadata$inflammation_raw,
    diagnosis = diagnosis, inflamed = NA_character_,
    characteristic[, c("montreal_location", "disease_subtype", "population"), drop = FALSE],
    stringsAsFactors = FALSE)
  out$cd_subtype <- ifelse(out$diagnosis != "CD", NA_character_,
    ifelse(grepl("ILEAL", toupper(out$disease_subtype)), "ILEAL_CD",
      ifelse(grepl("COLONIC", toupper(out$disease_subtype)), "COLONIC_CD", "OTHER_CD")))
  out$group <- ifelse(out$diagnosis == "CD", out$cd_subtype, out$diagnosis)
  if (any(out$diagnosis == "CD" & ((out$cd_subtype == "ILEAL_CD" & out$montreal_location != "L3") |
      (out$cd_subtype == "COLONIC_CD" & out$montreal_location != "L2") | out$cd_subtype == "OTHER_CD"))) {
    stop("The external CD phenotype labels and original L2/L3 codes disagree.")
  }
  missingness <- do.call(rbind, lapply(c("inflammation_raw", "treatment_raw"), function(variable) {
    do.call(rbind, lapply(sort(unique(metadata$diagnosis_raw)), function(label) {
      value <- metadata[[variable]][metadata$diagnosis_raw == label]
      missing <- is.na(value) | tolower(trimws(value)) %in% c("", "na", "n/a", "null", ".")
      data.frame(cohort = "GSE117993", variable = variable, diagnosis = label, n = length(value),
        n_missing = sum(missing), percent_missing = 100 * mean(missing), scope = "Used metadata only")
    }))
  }))
  if (any(missingness$n_missing != missingness$n)) {
    stop("External clinical metadata availability changed; review the adjustment specification.")
  }
  list(all = out, panel = block, source_metadata = metadata, characteristics = characteristic,
    missingness = missingness, participant_tokens = data.frame(sample_id = metadata$sample_id,
      source_title_token = source_token, stringsAsFactors = FALSE))
}

external_pool_context <- function(rectal, block) {
  external_unique_key(rectal$external_id, "rectal external IDs")
  if (anyNA(rectal$participant_id) || any(!nzchar(rectal$participant_id)) ||
      !identical(rectal$external_id, colnames(block))) stop("Context count columns and participant mapping are misaligned.")
  ids <- sort(unique(rectal$participant_id))
  groups <- lapply(ids, function(id) which(rectal$participant_id == id))
  if (any(vapply(groups, function(i) length(unique(rectal$diagnosis[i])) != 1L, logical(1)))) {
    stop("Conflicting diagnoses within a context participant.")
  }
  pooled <- vapply(groups, function(i) rowSums(block[, i, drop = FALSE]), numeric(nrow(block)))
  dimnames(pooled) <- list(rownames(block), ids)
  library_size <- vapply(groups, function(i) sum(rectal$raw_library_size[i]), numeric(1))
  if (sum(pooled) != sum(block) || sum(library_size) != sum(rectal$raw_library_size)) {
    stop("Context pooling does not conserve source counts and library totals.")
  }
  scores <- external_score_table(external_components(pooled), "COUNTS", library_size)
  out <- data.frame(participant_id = ids,
    diagnosis = vapply(groups, function(i) rectal$diagnosis[i[[1L]]], character(1)),
    n_rectal_samples = lengths(groups), scores, raw_library_size = library_size,
    n_zero_library_rectal_samples = vapply(groups, function(i) sum(rectal$raw_library_size[i] == 0), integer(1)),
    stringsAsFactors = FALSE)
  out$n_positive_library_rectal_samples <- out$n_rectal_samples - out$n_zero_library_rectal_samples
  out$IGHG3_CPM <- 1e6 * out$ighg3_count / out$raw_library_size
  out$TOTAL_IG_CPM <- 1e6 * out$total_ig_sum / out$raw_library_size
  out$IGHG3_LOG2_CPM_PLUS_0_5 <- log2(out$IGHG3_CPM + 0.5)
  out$TOTAL_IG_LOG2_CPM_PLUS_0_5 <- log2(out$TOTAL_IG_CPM + 0.5)
  out$cohort_role <- "CONTEXT_ONLY_PRIOR_PROJECT_EXPOSURE"
  list(participants = out, panel = pooled)
}

external_context_import <- function(count_path, metadata_path, emit) {
  emit("INFO | context_import | Reading all native HMP2 gene-count rows and 251 sample columns.")
  connection <- gzfile(count_path, open = "rt")
  table <- tryCatch(utils::read.delim(connection, sep = "\t", header = TRUE, quote = "", comment.char = "",
    check.names = FALSE, colClasses = c("character", rep("numeric", 251L)), stringsAsFactors = FALSE),
    finally = close(connection))
  if (nrow(table) != 55765L || ncol(table) != 252L) stop("Unexpected native HMP2 matrix dimensions.")
  genes <- external_clean(table[[1L]])
  sample_ids <- names(table)[-1L]
  external_unique_key(sample_ids, "context matrix samples")
  if (anyNA(genes) || any(!nzchar(genes)) ||
      any(grepl("^__|^N_(unmapped|multimapping|noFeature|ambiguous)$", genes)) ||
      !all(vapply(table[-1L], function(value) all(is.finite(value) & value >= 0 & value == round(value)), logical(1)))) {
    stop("Invalid context count values, symbols, or counting-summary rows.")
  }
  panel <- external_ig_panel()
  if (!all(vapply(panel, function(symbol) sum(genes == symbol) == 1L, logical(1)))) {
    stop("The context matrix must have exactly one native row for each of the 12 Ig genes.")
  }
  library_size <- vapply(table[-1L], sum, numeric(1))
  unit <- external_unit_audit(genes, table[[2L]], "GSE111889", basename(count_path))
  unit$igg_sum <- unit$igam_sum <- NA_real_
  block <- as.matrix(table[match(panel, genes), -1L, drop = FALSE])
  dimnames(block) <- list(panel, sample_ids)
  rm(table)
  sample <- data.frame(external_id = sample_ids, raw_library_size = unname(library_size),
    external_components(block), zero_library = unname(library_size) == 0, stringsAsFactors = FALSE)
  if (any(sample$total_ig_sum > sample$raw_library_size)) stop("Context sample Ig exceeds its full library.")
  metadata <- utils::read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("External ID", "Participant ID", "data_type", "diagnosis", "biopsy_location")
  if (!all(required %in% names(metadata))) stop("Missing required HMP2 source metadata columns.")
  slim <- as.data.frame(lapply(metadata[required], external_clean), stringsAsFactors = FALSE)
  names(slim) <- c("external_id", "participant_id", "data_type", "diagnosis_source", "biopsy_location")
  host <- slim[!is.na(slim$data_type) & slim$data_type == "host_transcriptomics", , drop = FALSE]
  mappings <- unique(host[!is.na(host$external_id) & host$external_id %in% sample_ids, , drop = FALSE])
  external_unique_key(mappings$external_id, "context metadata mappings; conflicting duplicates are not allowed")
  if (!setequal(sample_ids, mappings$external_id)) stop("A context count sample lacks a unique source metadata mapping.")
  joined <- merge(sample, mappings, by = "external_id", sort = FALSE)
  rectal <- joined[!is.na(joined$biopsy_location) & grepl("RECTUM", toupper(joined$biopsy_location)), , drop = FALSE]
  diagnosis <- toupper(rectal$diagnosis_source)
  rectal$diagnosis <- ifelse(grepl("^UC$|ULCERATIVE", diagnosis), "UC",
    ifelse(grepl("^CD$|CROHN", diagnosis), "CD",
      ifelse(grepl("NONIBD|CONTROL", diagnosis), "CONTROL", NA_character_)))
  if (nrow(rectal) != 100L || anyNA(rectal$diagnosis)) stop("Unexpected rectal context sample population.")
  pooled <- external_pool_context(rectal, block[, rectal$external_id, drop = FALSE])
  if (nrow(pooled$participants) != 90L ||
      sum(pooled$participants$diagnosis == "UC") != 25L || sum(pooled$participants$diagnosis == "CD") != 43L ||
      sum(pooled$participants$diagnosis == "CONTROL") != 22L) stop("Unexpected context participant geometry.")
  source_metadata <- metadata[!is.na(slim$data_type) & slim$data_type == "host_transcriptomics" &
    !is.na(slim$external_id) & slim$external_id %in% sample_ids, , drop = FALSE]
  list(all = pooled$participants, panel = pooled$panel, sample = sample, rectal = rectal,
    source_metadata = source_metadata, source_features = genes, unit_audit = unit,
    join_report = data.frame(item = c("matrix_columns", "host_transcriptomics_rows", "joined_rows", "unmatched_columns"),
      value = c(length(sample_ids), length(unique(host$external_id)), nrow(joined),
        length(setdiff(sample_ids, mappings$external_id)))),
    count_conservation = list(panel_counts = sum(block[, rectal$external_id, drop = FALSE]),
      pooled_panel_counts = sum(pooled$panel), rectal_library_total = sum(rectal$raw_library_size),
      pooled_library_total = sum(pooled$participants$raw_library_size)))
}

external_zero_component_audit <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  header <- strsplit(readLines(connection, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
  if (length(header) != 256L || anyDuplicated(header[-c(1L, 2L)])) stop("Unexpected GSE57945 sample header.")
  rows <- list()
  genes <- first_sample <- vector("list", 0L)
  offset <- 0L
  required <- external_ig_panel()[seq_len(7L)]
  repeat {
    lines <- readLines(connection, n = 2000L, warn = FALSE)
    if (!length(lines)) break
    fields <- strsplit(lines, "\t", fixed = TRUE)
    if (any(lengths(fields) != length(header))) stop("An RPKM source row has an inconsistent field count.")
    symbols <- external_clean(vapply(fields, function(value) value[[2L]], character(1)))
    first <- as.numeric(external_clean(vapply(fields, function(value) value[[3L]], character(1))))
    if (any(!is.finite(first)) || any(first < 0)) stop("Invalid source RPKM values.")
    genes[[length(genes) + 1L]] <- symbols
    first_sample[[length(first_sample) + 1L]] <- first
    for (i in which(symbols %in% required)) {
      value_text <- external_clean(fields[[i]][-c(1L, 2L)])
      values <- as.numeric(value_text)
      if (any(!is.finite(values)) || any(values < 0)) stop("Invalid RPKM score-component values.")
      rows[[length(rows) + 1L]] <- data.frame(row_index = offset + i,
        entrez_id = external_clean(fields[[i]][[1L]]), gene_symbol = symbols[[i]],
        n_fields = length(fields[[i]]), n_finite = length(values), n_positive = sum(values > 0),
        value_median = stats::median(values), value_max = max(values),
        first_five_values = paste(utils::head(value_text, 5L), collapse = " | "), stringsAsFactors = FALSE)
    }
    offset <- offset + length(lines)
  }
  genes <- unlist(genes, use.names = FALSE)
  first_sample <- unlist(first_sample, use.names = FALSE)
  audit <- do.call(rbind, rows)
  if (offset != 36372L || is.null(audit) || nrow(audit) != 7L ||
      anyDuplicated(audit$gene_symbol) || !setequal(audit$gene_symbol, required) || any(audit$n_positive != 0L)) {
    stop("The archived GSE57945 zero-component disposition requires review.")
  }
  unit <- external_unit_audit(genes, first_sample, "GSE57945", basename(path))
  unit$igg_sum <- unit$igam_sum <- NA_real_
  list(component_audit = audit, unit_audit = unit, n_samples = length(header) - 2L,
    disposition = "NON_ESTIMABLE_ZERO_COMPONENTS", score = NA_real_, auc = NA_real_)
}

external_auc <- function(score, positive) {
  result <- uc_auc_summary(score, positive)
  first <- score[positive == 1L]
  second <- score[positive == 0L]
  exact <- length(first) < 50L && length(second) < 50L && !anyDuplicated(c(first, second))
  test <- stats::wilcox.test(first, second, alternative = "two.sided", exact = exact, correct = TRUE)
  data.frame(auc = result$auc, ci_low = result$ci95_low, ci_high = result$ci95_high,
    n_positive = result$n_uc, n_negative = result$n_cd, p_value = test$p.value)
}

external_descriptions <- function(external, context) {
  old_features <- c("IGG_IGAM_LOG_RATIO", "IGHG3_FRACTION_OF_IG", "IGHG3_LOG", "TOTAL_IG_LOG")
  comparator_results <- function(table, cohort) {
    use <- table$diagnosis %in% c("UC", "CD")
    rows <- do.call(rbind, lapply(old_features, function(feature) {
      data.frame(external_auc(table[[feature]][use], as.integer(table$diagnosis[use] == "UC")),
        feature = feature, cohort = cohort, stringsAsFactors = FALSE)
    }))
    if (cohort == "GSE111889") rows$role <- "CONTEXT_ONLY_PRIOR_PROJECT_EXPOSURE"
    rows
  }
  group_summary <- function(table) do.call(rbind, lapply(split(table, table$diagnosis), function(data) {
    data.frame(diagnosis = data$diagnosis[[1L]], n = nrow(data), score_mean = mean(data$IGG_IGAM_LOG_RATIO),
      score_median = stats::median(data$IGG_IGAM_LOG_RATIO), score_sd = stats::sd(data$IGG_IGAM_LOG_RATIO))
  }))
  subtype_summary <- do.call(rbind, lapply(split(external, external$group), function(data) {
    data.frame(group = data$group[[1L]], n = nrow(data), score_mean = mean(data$IGG_IGAM_LOG_RATIO),
      score_median = stats::median(data$IGG_IGAM_LOG_RATIO), score_sd = stats::sd(data$IGG_IGAM_LOG_RATIO),
      total_ig_mean = mean(data$TOTAL_IG_LOG), ighg3_fraction_mean = mean(data$IGHG3_FRACTION_OF_IG))
  }))
  geometry <- as.data.frame(table(external$group, external$montreal_location, useNA = "ifany"))
  names(geometry) <- c("group", "montreal_location", "n")
  geometry <- geometry[geometry$n > 0L, , drop = FALSE]
  specifications <- list(c("UC", "ILEAL_CD"), c("UC", "COLONIC_CD"))
  subtype_contrasts <- do.call(rbind, lapply(specifications, function(pair) {
    data <- external[external$group %in% pair, , drop = FALSE]
    data.frame(external_auc(data$IGG_IGAM_LOG_RATIO, as.integer(data$group == pair[[1L]])),
      positive_group = pair[[1L]], negative_group = pair[[2L]], feature = "IGG_IGAM_LOG_RATIO",
      analysis_class = "POST_HOC_CONFOUNDER_PROBE", stringsAsFactors = FALSE)
  }))
  list(external_comparators_legacy = comparator_results(external, "GSE117993"),
    context_comparators_legacy = comparator_results(context, "GSE111889"),
    external_group_summary = group_summary(external), context_group_summary = group_summary(context),
    subtype_group_summary = subtype_summary, subtype_geometry = geometry, subtype_contrasts = subtype_contrasts)
}

external_comparisons <- function(external, context, development_comparisons) {
  calculate <- function(data, cohort, role) {
    data <- data[data$diagnosis %in% c("UC", "CD"), , drop = FALSE]
    scores <- data.frame(primary_ratio = data$IGG_IGAM_LOG_RATIO, IGHG3_abundance = data$IGHG3_ABUNDANCE,
      total_Ig_abundance = data$TOTAL_IG_ABUNDANCE, IGHG3_fraction = data$IGHG3_FRACTION_OF_IG,
      summed_IgG_abundance = data$IGG_ABUNDANCE)
    cohort_score_comparisons(scores, as.integer(data$diagnosis == "UC"), cohort, role)
  }
  first <- calculate(external, "GSE117993", "external score cohort; this added analysis is post hoc")
  second <- calculate(context, "GSE111889", "context only due to prior exposure")
  paired <- rbind(first$paired, second$paired)
  paired$p_holm_8_external_context_sensitivity <- stats::p.adjust(paired$p_two_sided, method = "holm", n = 8L)
  aucs <- rbind(first$aucs, second$aucs)
  historical <- paired[paired$comparator != "summed_IgG_abundance", , drop = FALSE]
  historical$p_holm_3_within_cohort <- NA_real_
  for (cohort in unique(historical$cohort)) {
    use <- historical$cohort == cohort
    historical$p_holm_3_within_cohort[use] <- stats::p.adjust(historical$p_two_sided[use], method = "holm", n = 3L)
  }
  historical$p_holm_6_sensitivity <- stats::p.adjust(historical$p_two_sided, method = "holm", n = 6L)
  historical$role <- ifelse(historical$cohort == "GSE117993", "PRIMARY_EXTERNAL_VALIDATION",
    "CONTEXT_ONLY_PRIOR_PROJECT_EXPOSURE")
  input_columns <- function(feature, cohort) {
    value <- c(primary_ratio = "IGG_IGAM_LOG_RATIO", IGHG3_abundance = "IGHG3_LOG",
      total_Ig_abundance = "TOTAL_IG_LOG", IGHG3_fraction = "IGHG3_FRACTION_OF_IG")[feature]
    value[cohort == "GSE111889" & feature == "IGHG3_abundance"] <- "IGHG3_LOG2_CPM_PLUS_0_5"
    value[cohort == "GSE111889" & feature == "total_Ig_abundance"] <- "TOTAL_IG_LOG2_CPM_PLUS_0_5"
    unname(value)
  }
  historical$input_column <- input_columns(historical$comparator, historical$cohort)
  historical_aucs <- aucs[aucs$feature != "summed_IgG_abundance", , drop = FALSE]
  historical_aucs$role <- ifelse(historical_aucs$cohort == "GSE117993", "PRIMARY_EXTERNAL_VALIDATION",
    "CONTEXT_ONLY_PRIOR_PROJECT_EXPOSURE")
  historical_aucs$input_column <- input_columns(historical_aucs$feature, historical_aucs$cohort)
  list(paired = paired, aucs = aucs, historical_paired = historical, historical_aucs = historical_aucs,
    all_paired = rbind(development_comparisons$paired, paired),
    all_aucs = rbind(development_comparisons$aucs, aucs))
}

external_supplementary <- function(external, comparisons) {
  table_s4 <- comparisons$paired[, c("cohort", "role", "comparator", "n_uc", "n_cd", "p_two_sided",
    "p_holm_4_within_cohort", "p_holm_8_external_context_sensitivity"), drop = FALSE]
  groups <- c("CONTROL", "ILEAL_CD", "COLONIC_CD", "UC")
  labels <- c("Control", "CD L3 (source: ileal CD)", "CD L2 (source: colonic CD)", "UC")
  table_s6 <- do.call(rbind, lapply(seq_along(groups), function(i) {
    score <- external$IGG_IGAM_LOG_RATIO[external$group == groups[[i]]]
    quartiles <- as.numeric(stats::quantile(score, probs = c(0.25, 0.75), type = 7L))
    data.frame(group = labels[[i]], n = length(score), mean = mean(score), sd = stats::sd(score),
      median = stats::median(score), q25 = quartiles[[1L]], q75 = quartiles[[2L]], stringsAsFactors = FALSE)
  }))
  list(table_s4 = table_s4, table_s6 = table_s6)
}

check_external_helpers <- function() {
  panel <- external_ig_panel()
  block <- matrix(c(rep(1, length(panel)), rep(2, length(panel)), rep(0, length(panel))),
    nrow = length(panel), dimnames = list(panel, c("s1", "s2", "zero")))
  components <- external_components(block[, 1L:2L, drop = FALSE])
  counts <- external_score_table(components, "COUNTS", c(100, 200))
  tpm <- external_score_table(components, "TPM")
  expected <- log2(c(4.5 / 3.5, 8.5 / 6.5))
  if (max(abs(counts$IGG_IGAM_LOG_RATIO - expected), abs(tpm$IGG_IGAM_LOG_RATIO - expected)) > 1e-12 ||
      max(abs(counts$IGHG3_ABUNDANCE - log2(10000 + 0.5))) > 1e-12 ||
      max(abs(tpm$IGHG3_ABUNDANCE - log2(c(1, 2) + 0.5))) > 1e-12) {
    stop("External score and expression-unit self-check failed.")
  }
  rectal <- data.frame(external_id = colnames(block), participant_id = c("p1", "p2", "p1"),
    diagnosis = c("CD", "UC", "CD"), raw_library_size = c(100, 200, 0))
  pooled <- external_pool_context(rectal, block)$participants
  if (!identical(pooled$participant_id, c("p1", "p2")) ||
      any(pooled$n_rectal_samples != c(2L, 1L)) || any(pooled$n_zero_library_rectal_samples != c(1L, 0L)) ||
      any(pooled$raw_library_size != c(100, 200)) ||
      max(abs(pooled$IGG_IGAM_LOG_RATIO - expected)) > 1e-12) stop("Context zero-sample pooling self-check failed.")
  invisible(TRUE)
}

validate_external <- function(bundle, project, record) {
  filenames <- c("31_GSE117993_scores.tsv", "31_GSE117993_analysis_table.tsv", "31b_GSE117993_panel_scores.tsv",
    "31c_GSE117993_subtype_analysis_table.tsv", "31c_GSE117993_characteristics.tsv", "31b_E5_external_comparators.tsv",
    "31b_E5_external_paired_auc.tsv",
    "31b_GSE117993_group_summary.tsv", "31c_subtype_geometry.tsv", "31c_subtype_group_summary.tsv",
    "31c_subtype_pairwise_contrasts.tsv", "31c_GSE57945_component_row_diagnosis.tsv",
    "32_GSE111889_sample_component_counts.tsv", "32_GSE111889_participant_table.tsv", "32_E4_context_cohort_results.tsv",
    "32_GSE111889_group_summary.tsv", "32_GSE111889_join_report.tsv", "30b_tar_member_index.tsv",
    "30c_external_unit_audit.tsv", "33c6_auc_intervals.tsv", "33c6_paired_auc.tsv", "Table_S15_expanded_paired_auc.tsv",
    "Table_S16_comparator_auc.tsv", "hmp2.tsv", "Table_S14_external_missingness.tsv")
  archive <- file.path(project, "reference", "external_expected.zip")
  if (!identical(unname(tools::md5sum(archive)), "a6f58ed22e71c5d1c6925930e6006880")) {
    stop("The external comparison archive failed its checksum.")
  }
  listing <- utils::unzip(archive, list = TRUE)$Name
  if (anyDuplicated(listing) || !setequal(listing, filenames)) stop("Unexpected external comparison archive members.")
  directory <- tempfile(pattern = "ibd_external_expected_")
  if (!dir.create(directory)) stop("Cannot create the temporary external comparison directory.")
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(archive, files = filenames, exdir = directory, unzip = "internal")
  read <- function(name) utils::read.delim(file.path(directory, name), check.names = FALSE, stringsAsFactors = FALSE)
  outcomes <- logical()
  compare <- function(observed, filename, keys, exact = character(), fields = NULL, expected_filter = NULL) {
    expected <- read(filename)
    if (!is.null(expected_filter)) expected <- expected_filter(expected)
    if (is.null(fields)) fields <- names(expected)
    outcomes <<- c(outcomes, compare_table(observed, expected, keys, filename, record,
      fields = fields, exact_fields = exact))
  }
  rounded_auc <- function(table) development_archive_view(table,
    decimal = c(auc = 6L, ci_low = 6L, ci_high = 6L), significant = "p_value")
  rounded_summary <- function(table) {
    fields <- setdiff(names(table)[vapply(table, is.numeric, logical(1))], "n")
    development_archive_view(table, decimal = stats::setNames(rep(4L, length(fields)), fields))
  }
  external <- bundle$external$all
  scores <- data.frame(sample_id = external$sample_id, components_found = 7L, igg_sum = external$igg_sum,
    igam_sum = external$igam_sum, score = external$IGG_IGAM_LOG_RATIO)
  compare(scores, "31_GSE117993_scores.tsv", "sample_id", exact = "components_found")
  analysis <- cbind(scores, external[, c("diagnosis_source", "diagnosis"), drop = FALSE])
  analysis <- analysis[analysis$diagnosis %in% c("UC", "CD"), , drop = FALSE]
  analysis$is_uc <- as.integer(analysis$diagnosis == "UC")
  compare(analysis, "31_GSE117993_analysis_table.tsv", "sample_id", exact = c("components_found", "is_uc"))
  compare(external, "31b_GSE117993_panel_scores.tsv", "sample_id")
  compare(external, "31c_GSE117993_subtype_analysis_table.tsv", "sample_id")
  compare(bundle$external$characteristics, "31c_GSE117993_characteristics.tsv", "sample_id")
  descriptions <- bundle$descriptions
  compare(rounded_auc(descriptions$external_comparators_legacy), "31b_E5_external_comparators.tsv", "feature",
    exact = c("n_positive", "n_negative"))
  historical <- bundle$comparisons$historical_paired
  use <- historical$cohort == "GSE117993"
  old_paired <- data.frame(cohort = historical$cohort[use], comparator = historical$input_column[use],
    auc_frozen_composite = round(historical$auc_primary[use], 6L), auc_comparator = round(historical$auc_comparator[use], 6L),
    delong_p_value = signif(historical$p_two_sided[use], 6L))
  compare(old_paired, "31b_E5_external_paired_auc.tsv", c("cohort", "comparator"))
  compare(rounded_summary(descriptions$external_group_summary), "31b_GSE117993_group_summary.tsv", "diagnosis", exact = "n")
  compare(descriptions$subtype_geometry, "31c_subtype_geometry.tsv", c("group", "montreal_location"), exact = "n")
  compare(rounded_summary(descriptions$subtype_group_summary), "31c_subtype_group_summary.tsv", "group", exact = "n")
  compare(rounded_auc(descriptions$subtype_contrasts), "31c_subtype_pairwise_contrasts.tsv",
    c("positive_group", "negative_group", "feature"), exact = c("n_positive", "n_negative"),
    expected_filter = function(table) table[table$positive_group == "UC" &
      table$negative_group %in% c("ILEAL_CD", "COLONIC_CD") & table$feature == "IGG_IGAM_LOG_RATIO", , drop = FALSE])
  compare(bundle$nonestimable$gse57945$component_audit, "31c_GSE57945_component_row_diagnosis.tsv", "gene_symbol",
    exact = c("row_index", "entrez_id", "n_fields", "n_finite", "n_positive", "value_median", "value_max"))
  component_fields <- c("igg_sum", "igam_sum", "total_ig_sum", "ighg3_count")
  compare(bundle$context$sample, "32_GSE111889_sample_component_counts.tsv", "external_id", exact = component_fields)
  compare(bundle$context$all, "32_GSE111889_participant_table.tsv", "participant_id",
    exact = c("n_rectal_samples", component_fields))
  compare(rounded_auc(descriptions$context_comparators_legacy), "32_E4_context_cohort_results.tsv", "feature",
    exact = c("n_positive", "n_negative"))
  compare(rounded_summary(descriptions$context_group_summary), "32_GSE111889_group_summary.tsv", "diagnosis", exact = "n")
  compare(bundle$context$join_report, "32_GSE111889_join_report.tsv", "item", exact = "value")
  compare(bundle$member_index, "30b_tar_member_index.tsv", c("study_id", "member_name"))
  unit <- development_archive_view(bundle$unit_audit, decimal = c(integer_share = 6L))
  for (field in c("column_sum", "max_value", "igg_sum", "igam_sum")) unit[[field]] <- signif(unit[[field]], 8L)
  unit$min_positive <- signif(unit$min_positive, 6L)
  compare(unit, "30c_external_unit_audit.tsv", c("study_id", "member_name"),
    exact = c("n_rows", "n_finite", "panel_genes_found", "duplicated_panel_rows"))
  compare(bundle$comparisons$historical_aucs, "33c6_auc_intervals.tsv", c("cohort", "feature"), exact = c("n_uc", "n_cd"))
  compare(historical, "33c6_paired_auc.tsv", c("cohort", "comparator"), exact = c("n_uc", "n_cd"))
  compare(bundle$comparisons$all_paired, "Table_S15_expanded_paired_auc.tsv", c("cohort", "comparator"), exact = c("n_uc", "n_cd"))
  compare(bundle$comparisons$all_aucs, "Table_S16_comparator_auc.tsv", c("cohort", "feature"), exact = c("n_uc", "n_cd"))
  compare(bundle$context$all, "hmp2.tsv", "participant_id", exact = c("n_rectal_samples", component_fields,
    "raw_library_size", "n_zero_library_rectal_samples", "n_positive_library_rectal_samples"))
  compare(bundle$external$missingness, "Table_S14_external_missingness.tsv", c("cohort", "variable", "diagnosis"),
    exact = c("n", "n_missing", "percent_missing"))
  final_zip <- file.path(project, "reference", "tables_reference.zip")
  mapping <- utils::read.delim(file.path(project, "config", "tables.tsv"), colClasses = "character", check.names = FALSE)
  expected_md5 <- unique(mapping$reference_zip_md5)
  if (length(expected_md5) != 1L || !identical(unname(tools::md5sum(final_zip)), expected_md5)) {
    stop("The supplementary comparison archive failed its checksum.")
  }
  current <- c("Table_S4_multiplicity_sensitivity.tsv", "Table_S6_external_group_summaries.tsv")
  utils::unzip(final_zip, files = current, exdir = directory, unzip = "internal")
  outcomes <- c(outcomes, compare_table(bundle$supplementary$table_s4, read(current[[1L]]),
    c("cohort", "comparator"), "current_Table_S4_complete", record, exact_fields = c("n_uc", "n_cd")))
  outcomes <- c(outcomes, compare_table(bundle$supplementary$table_s6, read(current[[2L]]),
    "group", "current_Table_S6_complete", record, exact_fields = "n"))
  all(outcomes)
}

run_external <- function(context) {
  project <- context$project
  work <- file.path(project, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  report_path <- file.path(work, "external_report.txt")
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
    emit(c("IBD external and HMP2 context source-to-result reproduction report", "Code revision: stage05.0",
      paste("Started UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      "Analysis inputs: deposited external TPM members, HMP2 native counts, and declared source metadata.",
      "GSE111889 remains context only because of prior project exposure.",
      "GSE83687 and GSE57945 source component availability is audited; no substitute score is constructed.",
      "Archived outputs are read only after computation; all calculations are in R."))
    check_validation_logic()
    record("PASS", "validator_self_checks", "Shared exact-count, identifier, numerical, missingness, and probability checks.")
    if (!requireNamespace("pROC", quietly = TRUE) || utils::packageVersion("pROC") != package_version("1.19.0.1")) {
      stop("This stage requires pROC 1.19.0.1.")
    }
    record("PASS", "package_pROC", as.character(utils::packageVersion("pROC")))
    check_external_helpers()
    record("PASS", "external_helper_self_checks", "Native ratios, TPM/CPM abundance transforms, pooling, and retained zero-library samples.")
    if (getRversion() != package_version("4.4.2")) record("WARN", "R_version",
      paste("Observed", as.character(getRversion()), "; archived analyses used 4.4.2."))
    if (!isTRUE(context$settings$verify_md5)) stop("External reconstruction requires checksum verification.")
    manifest <- context$read_manifest(file.path(project, "config", "inputs.tsv"))
    development_path <- file.path(work, "development.rds")
    if (!file.exists(development_path)) stop("Run and validate development first; work/development.rds is required.")
    development <- readRDS(development_path)
    if (!isTRUE(development$validation_passed) || !identical(development$code_revision, "stage04.2")) {
      stop("The development bundle is incomplete, incompatible, or not validated.")
    }
    for (identifier in c("development_counts", "gse193677_sample_metadata")) {
      specification <- manifest[manifest$input_id == identifier, , drop = FALSE]
      if (nrow(specification) != 1L ||
          !identical(unname(development$source_md5[identifier]), specification$expected_md5[[1L]])) {
        stop("Development cache source identities disagree with the input definitions.")
      }
    }
    development_md5 <- unname(tools::md5sum(development_path))
    development_comparisons <- development$comparisons
    if (!identical(development$settings$auc_positive, "UC") ||
        !identical(development$settings$auc_direction, "<") || development$settings$ratio_pseudocount != 0.5) {
      stop("The computed development score direction or pseudocount is incompatible.")
    }
    rm(development)
    record("PASS", "development_cache", paste("Validated stage04.2; cache MD5:", development_md5))
    paths <- source_md5 <- character()
    identifiers <- c("external_archive", "gse117993_sample_metadata", "context_counts", "context_metadata",
      "nonestimable_archive", "nonestimable_rpkm")
    for (identifier in identifiers) {
      specification <- manifest[manifest$input_id == identifier, , drop = FALSE]
      if (nrow(specification) != 1L || !context$safe_relative(specification$relative_path[[1L]]) ||
          is.na(specification$expected_md5[[1L]])) stop("Missing or invalid input definition: ", identifier)
      path <- file.path(context$settings$legacy_root, specification$relative_path[[1L]])
      if (identifier %in% names(context$settings$file_overrides)) path <- context$settings$file_overrides[[identifier]]
      path <- normalizePath(path.expand(path), winslash = "/", mustWork = TRUE)
      observed <- unname(tools::md5sum(path))
      if (file.info(path)$size != as.numeric(specification$expected_bytes[[1L]]) ||
          !identical(observed, specification$expected_md5[[1L]])) stop("Source bytes or checksum differ for ", identifier, ".")
      paths[identifier] <- path
      source_md5[identifier] <- observed
      record("PASS", identifier, paste("Archived source MD5:", observed))
    }
    imported <- external_archive_import(paths[["external_archive"]], "GSE117993", 190L, context$safe_relative, emit)
    external <- external_metadata_join(imported, paths[["gse117993_sample_metadata"]], context$read_manifest)
    if (sum(external$all$diagnosis == "UC") != 43L || sum(external$all$diagnosis == "CD") != 92L ||
        sum(external$all$diagnosis == "CONTROL") != 55L) stop("Unexpected external diagnosis counts.")
    record("PASS", "external_source_reconstruction", "190 GSMs; unique source title tokens; all 12 Ig genes; 43 UC, 92 CD, 55 controls; all biopsies rectal.")
    record("PASS", "external_metadata_missingness", "Inflammation and treatment are absent for all 190 records; no adjusted external association is fitted.")
    context_cohort <- external_context_import(paths[["context_counts"]], paths[["context_metadata"]], emit)
    record("PASS", "context_source_reconstruction", paste("Native rows:", length(context_cohort$source_features),
      "; matrix samples:", nrow(context_cohort$sample), "; rectal samples:", nrow(context_cohort$rectal),
      "; participants:", nrow(context_cohort$all), "; zero rectal libraries retained:", sum(context_cohort$rectal$zero_library),
      "; all pooled libraries positive; native counts conserved exactly."))
    unavailable <- external_archive_import(paths[["nonestimable_archive"]], "GSE83687", 134L, context$safe_relative, emit)
    zero <- external_zero_component_audit(paths[["nonestimable_rpkm"]])
    record("PASS", "nonestimable_source_audits", "All 134 GSE83687 members lack the seven required constant genes; all seven GSE57945 component rows are zero across 254 samples. No score AUC is assigned.")
    descriptions <- external_descriptions(external$all, context_cohort$all)
    comparisons <- external_comparisons(external$all, context_cohort$all, development_comparisons)
    record("PASS", "auc_and_multiplicity_calculations", "Fixed UC-positive DeLong AUCs and paired differences; legacy 3/6 and expanded 4/8 Holm families; development excluded from joint families.")
    supplementary <- external_supplementary(external$all, comparisons)
    bundle <- list(code_revision = "stage05.0", source_md5 = source_md5, development_cache_md5 = development_md5,
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), external = external, context = context_cohort,
      member_index = rbind(imported$member_index, unavailable$member_index),
      member_audit = rbind(imported$member_audit, unavailable$member_audit),
      unit_audit = rbind(imported$unit_probes, unavailable$unit_probes, context_cohort$unit_audit, zero$unit_audit),
      nonestimable = list(gse83687 = list(disposition = unavailable$disposition, source_features = unavailable$source_features),
        gse57945 = zero), descriptions = descriptions, comparisons = comparisons, supplementary = supplementary,
      settings = list(external_unit = "TPM; per-sample scale retained", context_unit = "native counts pooled within participant",
        context_abundance = "log2(1e6 * pooled component / pooled full-gene library + 0.5)",
        ratio_formula = "log2((G + 0.5) / (A + M + 0.5))", ratio_pseudocount = 0.5,
        auc_positive = "UC", auc_direction = "<", context_role = "CONTEXT_ONLY_PRIOR_PROJECT_EXPOSURE",
        descriptive_controls_excluded_from_primary_auc = TRUE,
        numerical_atol = 1e-10, probability_atol = 1e-14, numerical_rtol = 1e-6), session = utils::sessionInfo())
    emit(c("", "Comparison with archived external/context results and current supplementary tables"))
    passed <- validate_external(bundle, project, record)
    bundle$validation_passed <- isTRUE(passed) && counters[["FAIL"]] == 0L
    bundle$validation_counts <- counters
    cache <- file.path(work, "external.rds")
    saveRDS(bundle, cache, compress = "xz")
    emit(c("", "Newly computed external and context discrimination",
      utils::capture.output(print(comparisons$aucs, row.names = FALSE, digits = 10)),
      "", "Newly computed external phenotype and control summaries",
      utils::capture.output(print(supplementary$table_s6, row.names = FALSE, digits = 10)),
      paste("Cache:", normalizePath(cache, winslash = "/", mustWork = TRUE)),
      "Computed values retain full precision. Later stages must check validation_passed.",
      "A fatal run does not refresh a previous cache; review the current report before advancing.",
      "", "R session information", utils::capture.output(utils::sessionInfo())))
    status <- if (isTRUE(bundle$validation_passed)) 0L else 1L
  }, warning = function(warning) {
    record("WARN", "R_warning", conditionMessage(warning))
    invokeRestart("muffleWarning")
  }), error = function(error) {
    record("FAIL", "external_error", conditionMessage(error))
    status <<- 1L
  })
  emit(c("", paste("Checks: PASS =", counters[["PASS"]], "; WARN =", counters[["WARN"]], "; FAIL =", counters[["FAIL"]]),
    paste("External/context reproduction:", if (status == 0L) "PASS" else "BLOCKED"),
    "Remaining sensitivity analyses, component decomposition, final workbook, and figures: NOT_RUN",
    paste("Finished UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Report:", normalizePath(report_path, winslash = "/", mustWork = TRUE))))
  status
}
