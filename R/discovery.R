# Reconstruct TAURUS discovery from the native H5AD matrices.
# Consolidates the computations previously split across stages 09 through 14.
# Archived tables are opened only after the new results have been computed.

read_h5ad_column <- function(group, name) {
  object <- group[[name]]
  read_vector <- function(dataset) as.vector(dataset$read())
  if (inherits(object, "H5D")) return(read_vector(object))
  if (!inherits(object, "H5Group")) stop("Unsupported H5AD column: ", name)
  children <- names(object)
  if (all(c("categories", "codes") %in% children)) {
    categories <- as.character(read_vector(object[["categories"]]))
    codes <- as.integer(read_vector(object[["codes"]]))
    if (anyNA(codes) || any(codes < -1L | codes >= length(categories))) {
      stop("Invalid categorical codes in H5AD column: ", name)
    }
    result <- rep(NA_character_, length(codes))
    present <- codes >= 0L
    result[present] <- categories[codes[present] + 1L]
    return(result)
  }
  if (all(c("values", "mask") %in% children)) {
    values <- read_vector(object[["values"]])
    mask <- as.logical(read_vector(object[["mask"]]))
    if (length(values) != length(mask) || anyNA(mask)) {
      stop("Invalid nullable H5AD column: ", name)
    }
    values[mask] <- NA
    return(values)
  }
  stop("Unsupported H5AD column encoding: ", name)
}

normalize_token <- function(value) {
  gsub("[^a-z0-9]", "", tolower(trimws(as.character(value))), perl = TRUE)
}

normalize_disease <- function(value) {
  token <- normalize_token(value)
  result <- rep(NA_character_, length(token))
  result[token %in% c("uc", "ulcerativecolitis") | grepl("ulcerative", token)] <- "UC"
  result[token %in% "cd" | grepl("crohn", token)] <- "CD"
  result[token %in% c("hc", "healthy", "control") |
           grepl("healthy|control", token)] <- "HC"
  result
}

normalize_site <- function(value) {
  raw <- as.character(value)
  token <- normalize_token(value)
  result <- raw
  result[grepl("descending|leftcolon", token) | token == "left"] <- "LeftColon"
  result[grepl("ascending|rightcolon", token) | token == "right"] <- "RightColon"
  result[grepl("transverse", token)] <- "Transverse"
  result[grepl("sigmoid", token)] <- "Sigmoid"
  result[grepl("rect", token)] <- "Rectum"
  result[grepl("ileum", token)] <- "Ileum"
  result[is.na(raw) | !nzchar(trimws(raw))] <- NA_character_
  result
}

normalize_treatment <- function(value) {
  token <- normalize_token(value)
  result <- rep("Unknown", length(token))
  # The supplied sources contain Pre, Post, and missing labels.
  # Keep pre-treatment aliases from matching the substring in "treated".
  before <- grepl("^pre$|pretreatment|baseline|beforetreatment|untreated|naive", token)
  after <- grepl("^post$|posttreatment|aftertreatment|treated", token) & !before
  result[before] <- "Pre"
  result[after] <- "Post"
  result[is.na(value) | !nzchar(trimws(as.character(value)))] <- "Unknown"
  result
}

normalize_inflammation <- function(value) {
  token <- gsub("[^a-z]", "", tolower(as.character(value)), perl = TRUE)
  result <- rep(NA_character_, length(token))
  result[token %in% c("inflamed", "active", "involved")] <- "Inflamed"
  result[token %in% c("noninflamed", "noninvolved", "uninflamed")] <- "NonInflamed"
  result
}

aggregate_lineage <- function(path, lineage, emit) {
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  if (!all(c("X", "obs", "var") %in% names(h5))) {
    stop(lineage, ": missing /X, /obs, or /var.")
  }
  expression <- h5[["X"]]
  if (!inherits(expression, "H5Group") ||
      !all(c("data", "indices", "indptr") %in% names(expression))) {
    stop(lineage, ": the native /X matrix is not CSR encoded.")
  }
  pointers <- as.numeric(expression[["indptr"]]$read())
  if (length(pointers) < 2L || anyNA(pointers) || pointers[[1L]] != 0 ||
      any(diff(pointers) < 0) || any(pointers != floor(pointers))) {
    stop(lineage, ": invalid CSR row pointers.")
  }
  n_cells <- length(pointers) - 1L
  gene_id <- as.character(read_h5ad_column(h5[["var"]], "gene_id"))
  gene_symbol <- as.character(read_h5ad_column(h5[["var"]], "gene_symbol"))
  n_genes <- length(gene_id)
  if (n_genes != 33075L || length(gene_symbol) != n_genes ||
      anyNA(gene_id) || anyNA(gene_symbol) ||
      any(!nzchar(gene_id)) || any(!nzchar(gene_symbol)) ||
      anyDuplicated(gene_id) || anyDuplicated(gene_symbol)) {
    stop(lineage, ": incomplete, duplicated, or unexpected feature metadata.")
  }
  columns <- c("Disease", "Site", "Treatment", "Patient", "sample_id",
               "major", "minor", "Inflammation", "total_counts")
  if (!all(columns %in% names(h5[["obs"]]))) {
    stop(lineage, ": required cell metadata columns are missing.")
  }
  values <- lapply(columns, function(name) read_h5ad_column(h5[["obs"]], name))
  names(values) <- columns
  if (any(vapply(values, length, integer(1)) != n_cells)) {
    stop(lineage, ": expression rows and cell metadata lengths disagree.")
  }
  metadata <- data.frame(
    cell_index = seq_len(n_cells), patient = as.character(values$Patient),
    sample_id = as.character(values$sample_id),
    disease = normalize_disease(values$Disease),
    site = normalize_site(values$Site),
    treatment = normalize_treatment(values$Treatment),
    major = as.character(values$major), minor = as.character(values$minor),
    inflammation = as.character(values$Inflammation),
    total_counts = as.numeric(values$total_counts),
    stringsAsFactors = FALSE
  )
  candidate <- metadata$disease %in% c("CD", "UC") &
    !is.na(metadata$site) & metadata$site == "LeftColon" &
    metadata$treatment != "Post" &
    !is.na(metadata$patient) & nzchar(metadata$patient)
  if (!any(candidate)) stop(lineage, ": no eligible discovery cells.")
  target_names <- if (lineage == "myeloid") {
    c("MACROPHAGE_PRIMARY", "MONO_MACRO_BROAD_SENSITIVITY")
  } else c("FIBROBLAST_PRIMARY", "FIBPERI_BROAD_STROMAL_SENSITIVITY")
  masks <- if (lineage == "myeloid") {
    list(candidate & !is.na(metadata$minor) & metadata$minor == "Macro",
         candidate & !is.na(metadata$major) & metadata$major == "Mono_macro")
  } else {
    list(candidate & !is.na(metadata$minor) & metadata$minor != "Pericyte",
         candidate)
  }
  names(masks) <- target_names
  donor_key <- paste(metadata$disease, metadata$patient, sep = "::")
  targets <- list()
  donor_maps <- list()
  for (target in target_names) {
    rows <- which(masks[[target]])
    if (!length(rows)) stop("No cells for target: ", target)
    keys <- sort(unique(donor_key[rows]))
    donors <- data.frame(
      donor_key = keys, patient = sub("^[^:]+::", "", keys),
      disease = sub("::.*$", "", keys), stringsAsFactors = FALSE
    )
    unique_label <- function(key, field) {
      labels <- metadata[[field]][rows][donor_key[rows] == key]
      labels <- unique(labels[!is.na(labels) & nzchar(labels)])
      if (length(labels) > 1L) stop("Conflicting ", field, " labels for ", key)
      if (length(labels)) labels[[1L]] else NA_character_
    }
    donors$inflammation <- vapply(keys, unique_label, character(1), field = "inflammation")
    donors$sample_id <- vapply(keys, unique_label, character(1), field = "sample_id")
    donors <- donors[order(donors$disease, donors$patient), , drop = FALSE]
    if (anyDuplicated(donors$patient)) stop("Nonunique donor identifiers: ", target)
    accumulator <- new.env(parent = emptyenv())
    accumulator$counts <- matrix(0, nrow = nrow(donors), ncol = n_genes)
    accumulator$donors <- donors
    accumulator$n_cells <- integer(nrow(donors))
    accumulator$expected_total <- numeric(nrow(donors))
    targets[[target]] <- accumulator
    donor_maps[[target]] <- match(donor_key, donors$donor_key)
  }
  data_dataset <- expression[["data"]]
  index_dataset <- expression[["indices"]]
  candidate_rows <- which(candidate)
  max_difference <- 0
  emit(paste("INFO |", lineage, "| eligible cells:", length(candidate_rows),
             "; donors:", length(unique(donor_key[candidate_rows]))))
  for (position in seq_along(candidate_rows)) {
    row <- candidate_rows[[position]]
    start <- pointers[[row]] + 1
    end <- pointers[[row + 1L]]
    if (end >= start) {
      if (end > .Machine$integer.max) stop("CSR index exceeds supported integer range.")
      selection <- seq.int(as.integer(start), as.integer(end))
      counts <- as.numeric(data_dataset[selection])
      indices <- as.integer(index_dataset[selection]) + 1L
      if (length(counts) != length(indices) || anyNA(counts) || anyNA(indices) ||
          any(!is.finite(counts)) || any(counts < 0) ||
          any(abs(counts - round(counts)) > 1e-8) ||
          any(indices < 1L | indices > n_genes) || anyDuplicated(indices)) {
        stop(lineage, ": invalid raw CSR counts or indices at row ", row)
      }
    } else {
      counts <- numeric()
      indices <- integer()
    }
    difference <- abs(sum(counts) - metadata$total_counts[[row]])
    if (!is.finite(difference) || difference >= 1e-6) {
      stop(lineage, ": raw row sum does not match total_counts at row ", row)
    }
    max_difference <- max(max_difference, difference)
    for (target in target_names) {
      if (!masks[[target]][[row]]) next
      donor <- donor_maps[[target]][[row]]
      if (is.na(donor)) stop("Missing donor mapping: ", target)
      if (length(indices)) {
        targets[[target]]$counts[donor, indices] <-
          targets[[target]]$counts[donor, indices] + counts
      }
      targets[[target]]$n_cells[[donor]] <- targets[[target]]$n_cells[[donor]] + 1L
      targets[[target]]$expected_total[[donor]] <-
        targets[[target]]$expected_total[[donor]] + metadata$total_counts[[row]]
    }
    if (position %% 1000L == 0L || position == length(candidate_rows)) {
      emit(paste("INFO |", lineage, "| raw cells processed:",
                 position, "/", length(candidate_rows)))
    }
  }
  bundles <- list()
  for (target in target_names) {
    item <- targets[[target]]
    counts <- t(item$counts)
    if (any(!is.finite(counts)) || any(counts < 0) ||
        any(abs(counts - round(counts)) > 1e-8) ||
        any(counts > .Machine$integer.max)) {
      stop("Invalid pseudobulk matrix: ", target)
    }
    storage.mode(counts) <- "integer"
    rownames(counts) <- gene_id
    colnames(counts) <- item$donors$patient
    donors <- item$donors
    donors$n_cells <- item$n_cells
    donors$expected_total_counts <- item$expected_total
    donors$pseudobulk_library_size <- colSums(counts)
    donors$library_sum_difference <- donors$pseudobulk_library_size - donors$expected_total_counts
    donors$library_sum_exact <- abs(donors$library_sum_difference) < 1e-6
    donors$n_detected_features <- colSums(counts > 0)
    donors$target <- target
    donors$lineage <- lineage
    if (!all(donors$library_sum_exact)) stop("Pseudobulk conservation failure: ", target)
    bundles[[target]] <- list(
      target = target, lineage = lineage, counts = counts,
      donor_metadata = donors, gene_id = gene_id, gene_symbol = gene_symbol
    )
  }
  list(bundles = bundles, n_candidate_cells = length(candidate_rows),
       n_candidate_donors = length(unique(donor_key[candidate_rows])),
       max_cell_sum_difference = max_difference,
       cell_metadata = metadata[candidate_rows, setdiff(names(metadata), "total_counts"), drop = FALSE])
}

result_table <- function(test, analysis_id, hierarchy, test_type) {
  table <- as.data.frame(edgeR::topTags(
    test, n = Inf, adjust.method = "BH", sort.by = "PValue"
  )$table, stringsAsFactors = FALSE)
  feature_names <- rownames(table)
  rownames(table) <- NULL
  data.frame(
    analysis_id = analysis_id, hierarchy = hierarchy, test_type = test_type,
    feature_rowname = feature_names, table,
    direction = ifelse(table$logFC > 0, "UC_gt_CD",
                       ifelse(table$logFC < 0, "CD_gt_UC", "zero")),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

fit_target <- function(bundle) {
  target <- bundle$target
  donors <- bundle$donor_metadata
  donors$disease <- factor(donors$disease, levels = c("CD", "UC"))
  donors$inflammation <- factor(normalize_inflammation(donors$inflammation),
                                levels = c("NonInflamed", "Inflamed"))
  if (anyNA(donors$disease) || anyNA(donors$inflammation)) {
    stop("Missing disease or inflammation model labels: ", target)
  }
  included <- donors[donors$n_cells >= 10L, , drop = FALSE]
  included <- included[order(included$disease, included$patient), , drop = FALSE]
  counts <- bundle$counts[, included$patient, drop = FALSE]
  design <- stats::model.matrix(~ inflammation + disease, data = included)
  rownames(design) <- included$patient
  if (qr(design)$rank != ncol(design) || !"diseaseUC" %in% colnames(design)) {
    stop("The disease model is not estimable: ", target)
  }
  dge <- edgeR::DGEList(
    counts = counts,
    samples = data.frame(row.names = included$patient,
                         disease = included$disease, inflammation = included$inflammation),
    genes = data.frame(gene_id = bundle$gene_id, gene_symbol = bundle$gene_symbol,
                       row.names = bundle$gene_id, stringsAsFactors = FALSE)
  )
  keep <- edgeR::filterByExpr(dge, design = design, min.count = 10,
                              min.total.count = 15, large.n = 10, min.prop = 0.70)
  filtered <- edgeR::calcNormFactors(dge[keep, , keep.lib.sizes = FALSE], method = "TMM")
  if (!nrow(filtered$counts) || any(filtered$samples$lib.size <= 0)) {
    stop("Empty filtered universe or nonpositive library: ", target)
  }
  estimated <- edgeR::estimateDisp(filtered, design = design)
  fit <- edgeR::glmQLFit(estimated, design = design, robust = TRUE)
  coefficient <- match("diseaseUC", colnames(design))
  test <- edgeR::glmQLFTest(fit, coef = coefficient)
  primary <- target %in% c("MACROPHAGE_PRIMARY", "FIBROBLAST_PRIMARY")
  hierarchy <- if (primary) "PRIMARY" else "BROAD_SENSITIVITY"
  main_table <- result_table(test, target, hierarchy, "QLF")
  strict_metadata <- included[included$inflammation == "Inflamed", , drop = FALSE]
  strict_design <- stats::model.matrix(~ disease, data = strict_metadata)
  rownames(strict_design) <- strict_metadata$patient
  if (qr(strict_design)$rank != ncol(strict_design) ||
      !"diseaseUC" %in% colnames(strict_design)) {
    stop("The inflamed-only model is not estimable: ", target)
  }
  treat <- strict <- NULL
  if (primary) {
    treat_test <- edgeR::glmTreat(fit, coef = coefficient, lfc = 0.5)
    treat <- result_table(treat_test, target, "PRIMARY_EFFECT_SIZE", "glmTreat_lfc_0.5")
    strict_dge <- filtered[, strict_metadata$patient, keep.lib.sizes = FALSE]
    strict_dge <- edgeR::calcNormFactors(strict_dge, method = "TMM")
    strict_dge <- edgeR::estimateDisp(strict_dge, design = strict_design)
    strict_fit <- edgeR::glmQLFit(strict_dge, design = strict_design, robust = TRUE)
    strict_test <- edgeR::glmQLFTest(
      strict_fit, coef = match("diseaseUC", colnames(strict_design)))
    strict <- list(
      dge = strict_dge, fit = strict_fit, design = strict_design,
      table = result_table(strict_test, paste0(target, "_INFLAMED_ONLY"),
                           "STRICT_INFLAMED_SENSITIVITY", "QLF")
    )
  }
  effective <- filtered$samples$lib.size * filtered$samples$norm.factors
  filter_qc <- data.frame(
    target = target, n_native_features = length(keep),
    n_features_retained = sum(keep), n_features_removed = sum(!keep),
    fraction_retained = mean(keep), min_effective_lib_size = min(effective),
    median_effective_lib_size = stats::median(effective),
    max_effective_lib_size = max(effective),
    min_TMM_factor = min(filtered$samples$norm.factors),
    max_TMM_factor = max(filtered$samples$norm.factors)
  )
  list(dge = estimated, fit = fit, design = design, table = main_table,
       treat = treat, strict = strict, donor_metadata = included,
       strict_metadata = strict_metadata, gene_keep = keep, filter_qc = filter_qc)
}

select_candidates <- function(primary, treat, strict, broad) {
  candidates <- primary[!is.na(primary$FDR) & primary$FDR < 0.05 &
                          !is.na(primary$logFC) & abs(primary$logFC) >= 0.5, , drop = FALSE]
  key <- paste(candidates$analysis_id, candidates$gene_id, sep = "||")
  treat_index <- match(key, paste(treat$analysis_id, treat$gene_id, sep = "||"))
  strict_index <- match(
    paste0(candidates$analysis_id, "_INFLAMED_ONLY||", candidates$gene_id),
    paste(strict$analysis_id, strict$gene_id, sep = "||"))
  candidates$glmTreat_FDR <- treat$FDR[treat_index]
  candidates$glmTreat_logFC <- treat$logFC[treat_index]
  candidates$glmTreat_supported <- !is.na(candidates$glmTreat_FDR) & candidates$glmTreat_FDR < 0.05
  candidates$strict_inflamed_FDR <- strict$FDR[strict_index]
  candidates$strict_inflamed_logFC <- strict$logFC[strict_index]
  candidates$strict_direction_consistent <-
    !is.na(candidates$strict_inflamed_logFC) &
    sign(candidates$strict_inflamed_logFC) == sign(candidates$logFC)
  broad_id <- ifelse(candidates$analysis_id == "MACROPHAGE_PRIMARY",
                     "MONO_MACRO_BROAD_SENSITIVITY", "FIBPERI_BROAD_STROMAL_SENSITIVITY")
  broad_index <- match(paste(broad_id, candidates$gene_id, sep = "||"),
                       paste(broad$analysis_id, broad$gene_id, sep = "||"))
  candidates$broad_analysis_id <- broad_id
  candidates$broad_FDR <- broad$FDR[broad_index]
  candidates$broad_logFC <- broad$logFC[broad_index]
  candidates$broad_direction_consistent <- !is.na(candidates$broad_logFC) &
    sign(candidates$broad_logFC) == sign(candidates$logFC)
  candidates$evidence_class <- ifelse(candidates$glmTreat_supported,
    "CORE_EFFECT_SIZE_CONFIRMED", "PRIMARY_QL_HIGH_CONFIDENCE_EFFECT_SIZE_UNCONFIRMED")
  candidates <- candidates[order(candidates$analysis_id,
    candidates$evidence_class != "CORE_EFFECT_SIZE_CONFIRMED",
    candidates$FDR, candidates$gene_symbol), , drop = FALSE]
  candidates$signature_rank <- integer(nrow(candidates))
  for (target in unique(candidates$analysis_id)) {
    indices <- which(candidates$analysis_id == target)
    candidates$signature_rank[indices] <- seq_along(indices)
  }
  rownames(candidates) <- NULL
  candidates
}

compare_table <- function(observed, expected, keys, label, record,
                          fields = names(expected), exact_fields = character()) {
  missing <- setdiff(c(keys, fields), names(observed))
  if (length(missing)) {
    record("FAIL", label, paste("Missing columns:", paste(missing, collapse = ", ")))
    return(FALSE)
  }
  make_key <- function(table) {
    do.call(paste, c(lapply(table[keys], as.character), sep = "||"))
  }
  actual_keys <- make_key(observed)
  expected_keys <- make_key(expected)
  if (anyDuplicated(actual_keys) || anyDuplicated(expected_keys) ||
      !setequal(actual_keys, expected_keys) || nrow(observed) != nrow(expected)) {
    record("FAIL", label, "Row identifiers, row counts, or uniqueness differ.")
    return(FALSE)
  }
  observed <- observed[match(expected_keys, actual_keys), , drop = FALSE]
  issues <- character()
  largest_difference <- 0
  for (field in fields) {
    actual <- observed[[field]]
    reference <- expected[[field]]
    if (!identical(is.na(actual), is.na(reference))) {
      issues <- c(issues, paste0(field, ": missing-value pattern differs"))
      next
    }
    present <- !is.na(reference)
    if (is.numeric(reference)) {
      actual <- as.numeric(actual[present])
      reference <- as.numeric(reference[present])
      if (any(!is.finite(actual)) || any(!is.finite(reference))) {
        issues <- c(issues, paste0(field, ": a nonfinite numerical value was found"))
        next
      }
      difference <- abs(actual - reference)
      if (length(difference)) largest_difference <- max(largest_difference, difference)
      if (field %in% exact_fields) {
        if (!identical(actual, reference)) {
          issues <- c(issues, paste0(field, ": exact numerical values differ"))
        }
      } else {
        absolute_tolerance <- if (grepl("(PValue|FDR|_p|_q|p_value|p_two_sided)$|^BH_|^p_holm_",
                                        field, ignore.case = TRUE)) 1e-14 else 1e-10
        if (any(difference > absolute_tolerance + 1e-6 * abs(reference))) {
          issues <- c(issues, paste0(field, ": numerical tolerance exceeded (maximum absolute difference ",
            format(max(difference), digits = 8L, scientific = TRUE), ")"))
        }
      }
    } else if (!identical(as.character(actual[present]), as.character(reference[present]))) {
      issues <- c(issues, paste0(field, ": exact values differ"))
    }
  }
  passed <- !length(issues)
  detail <- if (passed) {
    paste("Rows:", nrow(expected), "; maximum numerical absolute difference:",
          format(largest_difference, digits = 8, scientific = TRUE))
  } else paste(issues, collapse = "; ")
  record(if (passed) "PASS" else "FAIL", label, detail)
  passed
}

check_validation_logic <- function() {
  expected <- data.frame(id = c("a", "b"), count = c(1000000L, 2L),
                          estimate = c(1.25, NA_real_), stringsAsFactors = FALSE)
  observed <- data.frame(id = c("b", "a"), count = c(2, 1000000),
                          estimate = c(NA_real_, 1.25), stringsAsFactors = FALSE)
  compare <- function(actual) compare_table(
    actual, expected, "id", "validation_self_check", function(...) NULL,
    exact_fields = "count")
  if (!isTRUE(compare(observed))) stop("Validation self-check failed: aligned exact values.")
  changed_count <- observed
  changed_count$count[[1L]] <- 3
  if (isTRUE(compare(changed_count))) stop("Validation failed to detect a count difference.")
  changed_estimate <- observed
  changed_estimate$estimate[[2L]] <- 1.251
  if (isTRUE(compare(changed_estimate))) stop("Validation failed to detect an estimate difference.")
  changed_missingness <- observed
  changed_missingness$estimate[[1L]] <- 0
  if (isTRUE(compare(changed_missingness))) stop("Validation failed to detect missingness.")
  duplicate_keys <- observed
  duplicate_keys$id[] <- "a"
  if (isTRUE(compare(duplicate_keys))) stop("Validation failed to detect duplicate keys.")
  for (field in c("replication_PValue", "replication_P", "family_BH_q", "BH_compartment", "BH_all_13",
                  "p_two_sided", "p_holm_4_within_cohort", "crude_p", "adjusted_p")) {
    expected_probability <- data.frame(id = "p", value = 1e-7, stringsAsFactors = FALSE)
    names(expected_probability)[[2L]] <- field
    changed_probability <- expected_probability
    changed_probability[[field]] <- 1.0005e-7
    if (isTRUE(compare_table(changed_probability, expected_probability, "id",
      "probability_tolerance_self_check", function(...) NULL))) {
      stop("Validation did not apply the probability tolerance to ", field, ".")
    }
  }
  if (!isTRUE(package_version("1.7-1") == package_version("1.7.1"))) {
    stop("Package-version normalization self-check failed.")
  }
  invisible(TRUE)
}

validate_discovery <- function(bundle, project, record) {
  archive <- file.path(project, "reference", "discovery_expected.zip")
  expected_md5 <- "227811ef48fc04d8d84b8e17200f46c5"
  if (!identical(unname(tools::md5sum(archive)), expected_md5)) {
    stop("The archived discovery comparison tables failed their checksum.")
  }
  directory <- tempfile(pattern = "ibd_discovery_expected_")
  if (!dir.create(directory)) stop("Cannot create the temporary comparison directory.")
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  expected_files <- c(
    "11_taurus_pseudobulk_library_qc.tsv",
    "12_primary_design_donor_geometry.tsv", "12_strict_inflamed_donor_geometry.tsv",
    "12_abundance_filter_tmm_qc.tsv", "13_primary_ql_all_genes.tsv",
    "13_primary_glmTreat_lfc0p5_all_genes.tsv",
    "13_strict_inflamed_ql_all_genes.tsv", "13_broad_sensitivity_ql_all_genes.tsv",
    "14_discovery_signature_freeze.tsv")
  listing <- utils::unzip(archive, list = TRUE)$Name
  if (!setequal(listing, expected_files) || anyDuplicated(listing)) {
    stop("The discovery comparison archive has unexpected members.")
  }
  utils::unzip(archive, files = expected_files, exdir = directory, unzip = "internal")
  read <- function(name) utils::read.delim(
    file.path(directory, name), check.names = FALSE, stringsAsFactors = FALSE)
  outcomes <- logical()
  compare <- function(observed, filename, keys, fields = NULL, exact_fields = character()) {
    expected <- read(filename)
    if (is.null(fields)) fields <- names(expected)
    passed <- compare_table(observed, expected, keys, filename, record,
                             fields, exact_fields)
    outcomes <<- c(outcomes, passed)
  }
  compare(bundle$library_qc, expected_files[[1L]], c("target", "patient"),
          exact_fields = names(bundle$library_qc))
  compare(bundle$donor_geometry, expected_files[[2L]], c("target", "disease"),
          exact_fields = names(bundle$donor_geometry))
  compare(bundle$strict_geometry, expected_files[[3L]], c("target", "disease"),
          exact_fields = names(bundle$strict_geometry))
  compare(bundle$filter_qc, expected_files[[4L]], "target",
          exact_fields = c("n_native_features", "n_features_retained", "n_features_removed"))
  compare(bundle$tables$primary, expected_files[[5L]], c("analysis_id", "gene_id"))
  compare(bundle$tables$treat, expected_files[[6L]], c("analysis_id", "gene_id"))
  compare(bundle$tables$strict, expected_files[[7L]], c("analysis_id", "gene_id"))
  compare(bundle$tables$broad, expected_files[[8L]], c("analysis_id", "gene_id"))
  signature_fields <- setdiff(names(read(expected_files[[9L]])),
    c("replication_status", "labels_used_for_freeze", "expression_values_read_in_stage14"))
  compare(bundle$candidates, expected_files[[9L]], c("analysis_id", "gene_id"),
          fields = signature_fields, exact_fields = "signature_rank")
  final_zip <- file.path(project, "reference", "tables_reference.zip")
  map <- utils::read.delim(file.path(project, "config", "tables.tsv"),
                          colClasses = "character", check.names = FALSE)
  reference_md5 <- unique(map$reference_zip_md5)
  if (length(reference_md5) != 1L ||
      !identical(unname(tools::md5sum(final_zip)), reference_md5)) {
    stop("The final supplementary reference archive failed its checksum.")
  }
  final_name <- "Table_S2_candidate_replication.tsv"
  utils::unzip(final_zip, files = final_name, exdir = directory, unzip = "internal")
  final <- read(final_name)
  observed <- data.frame(
    compartment = ifelse(bundle$candidates$analysis_id == "MACROPHAGE_PRIMARY",
                          "Macrophage", "Fibroblast"),
    gene_symbol = bundle$candidates$gene_symbol,
    ensembl_id = bundle$candidates$gene_id,
    discovery_log2FC = bundle$candidates$logFC, stringsAsFactors = FALSE)
  outcomes <- c(outcomes, compare_table(
    observed, final, c("compartment", "ensembl_id"), "current_Table_S2_discovery",
    record, fields = names(observed)))
  all(outcomes)
}

run_discovery <- function(context) {
  project <- context$project
  work <- file.path(project, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  report_path <- file.path(work, "discovery_report.txt")
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
    emit(c("IBD single-cell discovery reproduction report",
           "Code revision: stage02.0",
           paste("Started UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
           "Analysis inputs: the two declared native TAURUS H5AD files.",
           "Archived outputs are used only for comparison after computation."))
    check_validation_logic()
    record("PASS", "validator_self_checks",
           "Row alignment, exact counts, numerical discrepancies, missingness, duplicate keys, and version notation.")
    required <- c(edgeR = "4.4.2", limma = "3.62.2", hdf5r = "1.3.12",
                  statmod = "1.5.2", Matrix = "1.7-1")
    for (package in names(required)) {
      if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package)
      observed_version <- utils::packageVersion(package)
      if (observed_version != package_version(required[[package]])) {
        stop("Package version differs from the reproduction target: ", package,
             " observed ", as.character(observed_version), " expected ", required[[package]])
      }
      record("PASS", paste0("package_", package), as.character(observed_version))
    }
    if (getRversion() != package_version("4.4.2")) {
      record("WARN", "R_version", paste("Observed", as.character(getRversion()),
                                       "; archived analyses used 4.4.2."))
    }
    settings <- context$settings
    if (!isTRUE(settings$verify_md5)) stop("Discovery requires checksum verification.")
    manifest <- context$read_manifest(file.path(project, "config", "inputs.tsv"))
    paths <- character()
    source_md5 <- character()
    for (identifier in c("discovery_myeloid", "discovery_fibperi")) {
      spec <- manifest[manifest$input_id == identifier, , drop = FALSE]
      if (nrow(spec) != 1L || !context$safe_relative(spec$relative_path[[1L]])) {
        stop("Missing or invalid input definition: ", identifier)
      }
      path <- file.path(settings$legacy_root, spec$relative_path[[1L]])
      if (identifier %in% names(settings$file_overrides)) path <- settings$file_overrides[[identifier]]
      path <- normalizePath(path.expand(path), winslash = "/", mustWork = TRUE)
      if (file.info(path)$size != as.numeric(spec$expected_bytes[[1L]])) {
        stop("Input byte count differs: ", identifier)
      }
      digest <- unname(tools::md5sum(path))
      if (is.na(digest) || !identical(digest, spec$expected_md5[[1L]])) {
        stop("Input checksum differs: ", identifier)
      }
      paths[[identifier]] <- path
      source_md5[[identifier]] <- digest
      record("PASS", identifier, paste("Source MD5:", digest))
    }
    myeloid <- aggregate_lineage(paths[["discovery_myeloid"]], "myeloid", emit)
    fibroblast <- aggregate_lineage(paths[["discovery_fibperi"]], "fibperi", emit)
    pseudobulk <- c(myeloid$bundles, fibroblast$bundles)
    if (!identical(myeloid$bundles[[1L]]$gene_id, fibroblast$bundles[[1L]]$gene_id) ||
        !identical(myeloid$bundles[[1L]]$gene_symbol, fibroblast$bundles[[1L]]$gene_symbol)) {
      stop("Native feature identities or ordering differ between lineages.")
    }
    record("PASS", "native_features", "Both lineages share 33075 unique ordered features.")
    library_qc <- do.call(rbind, lapply(pseudobulk, function(item) item$donor_metadata))
    record("PASS", "count_conservation",
           paste("Donor-target libraries checked:", nrow(library_qc),
                 "; maximum row-sum discrepancy:",
                 max(myeloid$max_cell_sum_difference, fibroblast$max_cell_sum_difference)))
    models <- list()
    for (target in names(pseudobulk)) {
      emit(paste("INFO | Fitting", target))
      models[[target]] <- fit_target(pseudobulk[[target]])
      model <- models[[target]]
      emit(paste("INFO |", target, "| donors:", nrow(model$donor_metadata),
                 "; CD:", sum(model$donor_metadata$disease == "CD"),
                 "; UC:", sum(model$donor_metadata$disease == "UC"),
                 "; retained genes:", nrow(model$dge$counts)))
    }
    primary_targets <- c("MACROPHAGE_PRIMARY", "FIBROBLAST_PRIMARY")
    broad_targets <- c("MONO_MACRO_BROAD_SENSITIVITY", "FIBPERI_BROAD_STROMAL_SENSITIVITY")
    tables <- list(
      primary = do.call(rbind, lapply(models[primary_targets], function(item) item$table)),
      treat = do.call(rbind, lapply(models[primary_targets], function(item) item$treat)),
      strict = do.call(rbind, lapply(models[primary_targets], function(item) item$strict$table)),
      broad = do.call(rbind, lapply(models[broad_targets], function(item) item$table))
    )
    candidates <- select_candidates(tables$primary, tables$treat, tables$strict, tables$broad)
    geometry <- strict_geometry <- list()
    for (target in names(models)) {
      model <- models[[target]]
      all_donors <- pseudobulk[[target]]$donor_metadata
      geometry[[target]] <- data.frame(
        target = target, disease = c("CD", "UC"),
        n_donors_included = c(sum(model$donor_metadata$disease == "CD"),
                              sum(model$donor_metadata$disease == "UC")),
        n_donors_excluded = c(sum(all_donors$disease == "CD" & all_donors$n_cells < 10L),
                              sum(all_donors$disease == "UC" & all_donors$n_cells < 10L)))
      strict_geometry[[target]] <- data.frame(
        target = target, disease = c("CD", "UC"),
        n_strict_inflamed_donors = c(sum(model$strict_metadata$disease == "CD"),
                                    sum(model$strict_metadata$disease == "UC")))
    }
    bundle <- list(
      code_revision = "stage02.0", source_md5 = source_md5,
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      pseudobulk = pseudobulk, models = models, tables = tables,
      candidates = candidates, library_qc = library_qc,
      donor_geometry = do.call(rbind, geometry),
      strict_geometry = do.call(rbind, strict_geometry),
      filter_qc = do.call(rbind, lapply(models, function(item) item$filter_qc)),
      candidate_cell_metadata = list(myeloid = myeloid$cell_metadata,
                                      fibperi = fibroblast$cell_metadata),
      settings = list(min_cells = 10L, primary_design = "~ inflammation + disease",
        strict_design = "~ disease", reference_disease = "CD",
        filter = list(min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.70),
        fdr_cutoff = 0.05, abs_log2fc_cutoff = 0.5, glmTreat_lfc = 0.5,
        ql_robust = TRUE, numerical_atol = 1e-10,
        probability_atol = 1e-14, numerical_rtol = 1e-6),
      session = utils::sessionInfo()
    )
    emit(c("", "Comparison with archived results and the current supplementary table"))
    passed <- validate_discovery(bundle, project, record)
    bundle$validation_passed <- isTRUE(passed) && counters[["FAIL"]] == 0L
    bundle$validation_counts <- counters
    cache <- file.path(work, "discovery.rds")
    saveRDS(bundle, cache, compress = "xz")
    emit(c("", "Newly computed primary discovery candidates",
           utils::capture.output(print(candidates[, c("analysis_id", "gene_symbol",
             "logFC", "FDR", "glmTreat_FDR", "strict_inflamed_FDR", "broad_FDR")],
             row.names = FALSE, digits = 10)),
           paste("Computed candidate count:", nrow(candidates)),
           paste("Cache:", normalizePath(cache, winslash = "/", mustWork = TRUE)),
           "The cache contains computed values even if comparison fails; later stages must check validation_passed.",
           "", "R session information", utils::capture.output(utils::sessionInfo())))
    status <- if (isTRUE(bundle$validation_passed)) 0L else 1L
  }, warning = function(warning) {
    record("WARN", "R_warning", conditionMessage(warning))
    invokeRestart("muffleWarning")
  }), error = function(error) {
    record("FAIL", "discovery_error", conditionMessage(error))
    status <<- 1L
  })
  emit(c("", paste("Checks: PASS =", counters[["PASS"]], "; WARN =", counters[["WARN"]],
                    "; FAIL =", counters[["FAIL"]]),
         paste("Discovery reproduction:", if (status == 0L) "PASS" else "BLOCKED"),
         "Replication, bulk analyses, final workbook, and figures: NOT_RUN",
         paste("Finished UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
         paste("Report:", normalizePath(report_path, winslash = "/", mustWork = TRUE))))
  status
}
