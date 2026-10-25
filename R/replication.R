# Reconstruct GSE214695 candidate replication from the deposited 10x files.
# Consolidates original stages 15--17 and the later pooled-family BH sensitivity.
# Shared comparison and label helpers are loaded from discovery.R by run.R.

replication_sample_map <- function() {
  data.frame(
    sample = c(paste0("HC", 1:6), paste0("UC", 1:6), paste0("CD", 1:6)),
    sample_id = c(paste0("HC-", 1:6), paste0("UC-", 1:6), paste0("CD-", 1:6)),
    gsm = paste0("GSM", 6614348:6614365), stringsAsFactors = FALSE)
}

read_replication_metadata <- function(path, sample_map, read_manifest) {
  metadata <- read_manifest(path)
  needed <- c("study_id", "sample_id", "title", "diagnosis_raw", "anatomy_raw")
  if (!all(needed %in% names(metadata)) || nrow(metadata) != nrow(sample_map) ||
      anyNA(metadata$sample_id) || anyDuplicated(metadata$sample_id) ||
      !setequal(metadata$sample_id, sample_map$gsm) ||
      anyNA(metadata$study_id) || any(metadata$study_id != "GSE214695")) {
    stop("The source GSE214695 metadata do not define 18 unique expected GSMs.")
  }
  metadata <- metadata[match(sample_map$gsm, metadata$sample_id), , drop = FALSE]
  title_token <- gsub("[[:space:]]", "", sub(",.*$", "", metadata$title))
  if (anyNA(title_token) || !identical(title_token, sample_map$sample)) {
    stop("Source titles disagree with the annotation sample-to-GSM map.")
  }
  disease <- normalize_disease(metadata$diagnosis_raw)
  # Preserve the original stage-02 anatomy harmonization, including sigma.
  raw <- metadata$anatomy_raw
  anatomy <- rep(NA_character_, length(raw))
  anatomy[grepl("rectum|rectal", raw, ignore.case = TRUE)] <- "Rectum"
  anatomy[grepl("sigmoid|(^|[^a-z])sigma([^a-z]|$)", raw,
                ignore.case = TRUE)] <- "Sigmoid"
  anatomy[grepl("descending|left.?colon", raw, ignore.case = TRUE)] <- "LeftColon"
  anatomy[grepl("ascending|right.?colon", raw, ignore.case = TRUE)] <- "RightColon"
  anatomy[grepl("transverse", raw, ignore.case = TRUE)] <- "Transverse"
  region <- ifelse(anatomy %in% c("LeftColon", "Sigmoid", "Rectum"), "Distal",
                   ifelse(anatomy %in% c("RightColon", "Transverse"), "Proximal", NA_character_))
  if (anyNA(disease) || anyNA(anatomy) || anyNA(region) ||
      any(disease != sub("-[0-9]+$", "", sample_map$sample_id))) {
    stop("Source diagnosis or anatomy cannot be mapped consistently for every donor.")
  }
  data.frame(sample_id = sample_map$sample_id, study_id = metadata$study_id,
    diagnosis_harmonized = disease, anatomy_harmonized = anatomy,
    gsm = sample_map$gsm, broad_region = region, stringsAsFactors = FALSE)
}

prepare_replication_annotation <- function(path, metadata, sample_map) {
  connection <- gzfile(path, open = "rt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  annotation <- utils::read.table(connection, header = FALSE, sep = ",", quote = "\"",
    comment.char = "", check.names = FALSE, stringsAsFactors = FALSE,
    colClasses = "character", na.strings = c("", "NA"), strip.white = TRUE, fill = FALSE)
  if (ncol(annotation) != 5L || nrow(annotation) < 2L) {
    stop("Unexpected source cell-annotation schema.")
  }
  names(annotation) <- c("annotation_row_id", "sample", "cell_id",
                         "annotation", "nanostring_reference")
  # The deposited first record is the header, also skipped by the original reader.
  if (isTRUE(grepl("^[ACGT]{16}-1$", annotation$cell_id[[1L]])) ||
      annotation$sample[[1L]] %in% sample_map$sample) {
    stop("The first annotation record unexpectedly contains a cell rather than a header.")
  }
  annotation <- annotation[-1L, , drop = FALSE]
  rownames(annotation) <- NULL
  annotation$annotation_row <- seq_len(nrow(annotation))
  if (nrow(annotation) != 46702L ||
      anyNA(annotation[, c("sample", "cell_id", "annotation")]) ||
      any(!nzchar(annotation$annotation)) ||
      !all(grepl("^[ACGT]{16}-1$", annotation$cell_id))) {
    stop("Source annotation dimensions, labels, or 10x barcodes are invalid.")
  }
  index <- match(annotation$sample, sample_map$sample)
  if (anyNA(index)) stop("An annotation sample has no source GSM mapping.")
  annotation$sample_id <- sample_map$sample_id[index]
  annotation$gsm <- sample_map$gsm[index]
  annotation$composite_cell_id <- paste(annotation$gsm, annotation$cell_id, sep = "|")
  duplicate_keys <- unique(annotation$composite_cell_id[
    duplicated(annotation$composite_cell_id) | duplicated(annotation$composite_cell_id, fromLast = TRUE)])
  conflicting <- vapply(duplicate_keys, function(key) {
    length(unique(annotation$annotation[annotation$composite_cell_id == key])) > 1L
  }, logical(1))
  # No label precedence: every conflicting identity is excluded from both masks.
  conflict_keys <- duplicate_keys[conflicting]
  if (length(duplicate_keys) != length(conflict_keys)) {
    stop("Unexpected identical-label duplicates need explicit source review.")
  }
  annotation$conflict_excluded <- annotation$composite_cell_id %in% conflict_keys
  conflicts <- annotation[annotation$conflict_excluded, c("annotation_row", "gsm",
    "sample_id", "cell_id", "annotation", "composite_cell_id"), drop = FALSE]
  label <- tolower(gsub("[[:space:]]+", " ", trimws(annotation$annotation)))
  compact <- gsub("[^a-z0-9]", "", label)
  macrophage <- (
    grepl("macroph", label) |
    compact %in% c("m0", "m0ribhi", "m0ribhigh", "m2", "m22", "m1", "m1acod1",
                   "m1cxcl5", "ida", "idamacrophage", "idamacrophages") |
    grepl("^m0rib", compact) | grepl("^m1(acod1|cxcl5)$", compact) |
    grepl("^ida(macrophage|macrophages)?$", compact)
  ) & !grepl("monocyte|\\bmono\\b|dendritic|\\bdc\\b|neutroph|eosinoph|mast|plasma|lymph|epithel|fibro|pericyte|endothel|glia", label)
  fibroblast <- (
    grepl("fibroblast|myofibro", label) |
    compact %in% c("s1", "s2", "s2a", "s2b", "s3", "frc", "frcs",
                   "inflammatoryfibroblast", "inflammatoryfibroblasts") |
    grepl("^s[123](a|b)?$", compact)
  ) & !grepl("pericyte|smooth[ _-]*muscle|\\bsmc\\b|endothel|glia|macroph|monocyte|neutroph|eosinoph|mast|epithel", label)
  masks <- list(MACROPHAGE_PRIMARY = macrophage & !annotation$conflict_excluded,
                FIBROBLAST_PRIMARY = fibroblast & !annotation$conflict_excluded)
  if (any(masks[[1L]] & masks[[2L]])) stop("The two source compartment masks overlap.")
  membership <- do.call(rbind, lapply(names(masks), function(compartment) {
    data.frame(compartment = compartment, annotation[masks[[compartment]],
      c("annotation_row", "gsm", "sample_id", "cell_id", "annotation",
        "composite_cell_id", "conflict_excluded"), drop = FALSE], stringsAsFactors = FALSE)
  }))
  if (anyDuplicated(paste(membership$compartment, membership$composite_cell_id, sep = "||"))) {
    stop("Target membership contains a duplicated compartment-cell identity.")
  }
  support <- do.call(rbind, lapply(names(masks), function(compartment) {
    result <- data.frame(compartment = compartment, metadata, stringsAsFactors = FALSE)
    result$n_target_cells <- vapply(result$sample_id, function(sample_id) {
      sum(membership$sample_id == sample_id & membership$compartment == compartment)
    }, integer(1))
    result$support_ge_10 <- result$n_target_cells >= 10L
    result
  }))
  support <- support[order(support$compartment, support$sample_id), , drop = FALSE]
  eligible <- support[support$diagnosis_harmonized %in% c("CD", "UC") &
                        support$support_ge_10, , drop = FALSE]
  readiness <- do.call(rbind, lapply(split(eligible, eligible$compartment), function(data) {
    data$diagnosis_harmonized <- factor(data$diagnosis_harmonized, levels = c("CD", "UC"))
    data$broad_region <- factor(data$broad_region, levels = c("Distal", "Proximal"))
    design <- stats::model.matrix(~ broad_region + diagnosis_harmonized, data = data)
    n_cd <- sum(data$diagnosis_harmonized == "CD")
    n_uc <- sum(data$diagnosis_harmonized == "UC")
    data.frame(compartment = unique(data$compartment), supported_CD = n_cd,
      supported_UC = n_uc, supported_total = nrow(data), design_rank = qr(design)$rank,
      design_columns = paste(colnames(design), collapse = "|"),
      design_full_rank = qr(design)$rank == ncol(design),
      primary_ready = n_cd >= 4L && n_uc >= 4L && nrow(data) >= 8L &&
        qr(design)$rank == ncol(design), stringsAsFactors = FALSE)
  }))
  if (nrow(readiness) != 2L || !all(readiness$primary_ready) || nrow(eligible) != 24L) {
    stop("Source annotation does not support the original two 12-donor compartment models.")
  }
  list(conflicts = conflicts, membership = membership, support = support,
       eligible = eligible, readiness = readiness, n_annotation_rows = nrow(annotation),
       n_conflicting_identities = length(conflict_keys))
}

aggregate_replication <- function(archive, annotation, sample_map, emit) {
  suffixes <- c(barcodes = "_barcodes.tsv.gz", features = "_features.tsv.gz",
                matrix = "_matrix.mtx.gz")
  all_members <- unlist(lapply(seq_len(nrow(sample_map)), function(i) {
    paste0(sample_map$gsm[[i]], "_", sample_map$sample_id[[i]], suffixes)
  }), use.names = FALSE)
  listed <- utils::untar(archive, list = TRUE, tar = "internal")
  if (anyDuplicated(listed) || length(listed) != length(all_members) ||
      !setequal(listed, all_members)) {
    stop("The verified GSE214695 archive does not contain the 54 expected 10x members.")
  }
  eligible <- annotation$eligible
  donor_ids <- sort(unique(eligible$sample_id))
  targets <- c("MACROPHAGE_PRIMARY", "FIBROBLAST_PRIMARY")
  if (length(donor_ids) != 12L || !all(vapply(targets, function(target) {
    setequal(eligible$sample_id[eligible$compartment == target], donor_ids)
  }, logical(1)))) stop("The supported donor set differs between compartments.")
  donor_map <- sample_map[match(donor_ids, sample_map$sample_id), , drop = FALSE]
  required_members <- unlist(lapply(seq_len(nrow(donor_map)), function(i) {
    paste0(donor_map$gsm[[i]], "_", donor_map$sample_id[[i]], suffixes)
  }), use.names = FALSE)
  directory <- tempfile(pattern = "ibd_replication_10x_")
  if (!dir.create(directory)) stop("Cannot create the temporary 10x extraction directory.")
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  emit("INFO | archive_extraction | Extracting only the 36 CD/UC members into temporary storage.")
  utils::untar(archive, files = required_members, exdir = directory, tar = "internal")
  if (!all(file.exists(file.path(directory, required_members)))) {
    stop("Temporary extraction did not produce every required 10x member.")
  }
  read_features <- function(path) {
    connection <- gzfile(path, open = "rt", encoding = "UTF-8")
    on.exit(close(connection), add = TRUE)
    data <- utils::read.delim(connection, header = FALSE, sep = "\t", quote = "",
      comment.char = "", colClasses = "character", check.names = FALSE,
      stringsAsFactors = FALSE)
    if (ncol(data) < 2L || !nrow(data)) stop("Incomplete 10x feature annotation.")
    features <- data.frame(gene_id = data[[1L]], gene_symbol = data[[2L]],
                           stringsAsFactors = FALSE)
    if (anyNA(features) || any(!nzchar(features$gene_id)) ||
        anyDuplicated(features$gene_id)) stop("Invalid or duplicated 10x feature identifiers.")
    features
  }
  read_barcodes <- function(path) {
    connection <- gzfile(path, open = "rt", encoding = "UTF-8")
    on.exit(close(connection), add = TRUE)
    result <- trimws(readLines(connection, warn = FALSE))
    if (!length(result) || anyNA(result) || any(!nzchar(result)) || anyDuplicated(result)) {
      stop("Invalid or duplicated 10x barcodes.")
    }
    result
  }
  read_matrix <- function(path) {
    connection <- gzfile(path, open = "rt", encoding = "UTF-8")
    on.exit(close(connection), add = TRUE)
    result <- methods::as(Matrix::readMM(connection), "CsparseMatrix")
    values <- methods::slot(result, "x")
    if (nrow(result) == 0L || ncol(result) == 0L || any(!is.finite(values)) ||
        any(values < 0) || any(abs(values - round(values)) > 1e-8)) {
      stop("The deposited Matrix Market values are not nonnegative integer raw counts.")
    }
    result
  }
  counts <- list()
  library_qc <- list()
  feature_reference <- NULL
  for (i in seq_along(donor_ids)) {
    donor <- donor_ids[[i]]
    emit(paste("INFO | raw_pseudobulk | Donor", i, "of", length(donor_ids), ":", donor))
    paths <- file.path(directory, paste0(donor_map$gsm[[i]], "_", donor, suffixes))
    names(paths) <- names(suffixes)
    features <- read_features(paths[["features"]])
    barcodes <- read_barcodes(paths[["barcodes"]])
    matrix <- read_matrix(paths[["matrix"]])
    if (nrow(matrix) != nrow(features) || ncol(matrix) != length(barcodes)) {
      stop("10x dimensions disagree for donor ", donor, ".")
    }
    if (is.null(feature_reference)) feature_reference <- features
    if (!identical(feature_reference, features)) {
      stop("10x feature identity or ordering differs between donors.")
    }
    for (target in targets) {
      cells <- annotation$membership[annotation$membership$sample_id == donor &
        annotation$membership$compartment == target, , drop = FALSE]
      index <- match(cells$cell_id, barcodes)
      if (nrow(cells) < 10L || anyNA(index) || anyDuplicated(index)) {
        stop("Target support or one-to-one barcode matching failed for ", target, "/", donor, ".")
      }
      selected <- matrix[, index, drop = FALSE]
      summed <- as.numeric(Matrix::rowSums(selected))
      expected_total <- sum(Matrix::colSums(selected))
      if (sum(summed) != expected_total || sum(summed) <= 0 ||
          any(abs(summed - round(summed)) > 1e-8)) {
        stop("Pseudobulk raw-count conservation failed for ", target, "/", donor, ".")
      }
      if (is.null(counts[[target]])) {
        counts[[target]] <- base::matrix(0, nrow = nrow(features), ncol = length(donor_ids),
          dimnames = list(features$gene_id, donor_ids))
      }
      counts[[target]][, i] <- summed
      library_qc[[length(library_qc) + 1L]] <- data.frame(compartment = target,
        sample_id = donor, n_target_cells = nrow(cells), pseudobulk_library_size = sum(summed),
        n_genes_nonzero = sum(summed > 0), integer_raw_counts = TRUE,
        stringsAsFactors = FALSE)
    }
    # Close readers before deleting these temporary compressed members on Windows.
    unlink(paths, force = TRUE)
    rm(matrix, selected, barcodes)
  }
  metadata <- unique(eligible[, c("sample_id", "gsm", "diagnosis_harmonized",
                                  "anatomy_harmonized", "broad_region"), drop = FALSE])
  if (nrow(metadata) != length(donor_ids)) stop("A donor has conflicting source metadata.")
  metadata <- metadata[match(donor_ids, metadata$sample_id), , drop = FALSE]
  metadata$disease <- factor(metadata$diagnosis_harmonized, levels = c("CD", "UC"))
  metadata$broad_region <- factor(metadata$broad_region, levels = c("Distal", "Proximal"))
  design <- stats::model.matrix(~ broad_region + disease, data = metadata)
  if (anyNA(metadata$disease) || anyNA(metadata$broad_region) ||
      qr(design)$rank != ncol(design) || !"diseaseUC" %in% colnames(design)) {
    stop("The original replication design is incomplete or not full rank.")
  }
  list(counts = counts, features = feature_reference, metadata = metadata, design = design,
       library_qc = do.call(rbind, library_qc), n_archive_members = length(listed),
       n_expression_files_opened = length(required_members))
}

fit_replication_target <- function(target, pseudobulk, signature) {
  family_size <- c(MACROPHAGE_PRIMARY = 3L, FIBROBLAST_PRIMARY = 10L)
  design <- pseudobulk$design
  counts <- pseudobulk$counts[[target]]
  dge <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge, design = design, min.count = 10,
                             min.total.count = 15, large.n = 10, min.prop = 0.70)
  if (!any(keep)) stop("No genes pass the original replication filter in ", target, ".")
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  dge <- edgeR::estimateDisp(dge, design = design)
  fit <- edgeR::glmQLFit(dge, design = design, robust = TRUE)
  test <- edgeR::glmQLFTest(fit, coef = match("diseaseUC", colnames(design)))
  table <- as.data.frame(edgeR::topTags(test, n = Inf, sort.by = "none")$table,
                          stringsAsFactors = FALSE)
  table$gene_id <- rownames(table)
  table$gene_symbol <- pseudobulk$features$gene_symbol[match(table$gene_id, pseudobulk$features$gene_id)]
  candidates <- signature[signature$analysis_id == target,
    c("analysis_id", "gene_id", "gene_symbol", "direction", "logFC", "FDR"), drop = FALSE]
  if (nrow(candidates) != family_size[[target]]) stop("The discovery candidate family changed.")
  names(candidates)[4:6] <- c("TAURUS_direction", "TAURUS_logFC", "TAURUS_FDR")
  index <- match(candidates$gene_id, table$gene_id)
  raw_index <- match(candidates$gene_id, pseudobulk$features$gene_id)
  candidates$mapped_to_raw_feature <- !is.na(raw_index)
  candidates$passed_frozen_filter <- !is.na(index) & !is.na(table$F[index])
  candidates$estimable <- candidates$mapped_to_raw_feature & candidates$passed_frozen_filter
  if (any(candidates$mapped_to_raw_feature &
      candidates$gene_symbol != pseudobulk$features$gene_symbol[raw_index], na.rm = TRUE)) {
    stop("A discovery gene symbol disagrees with the replication source gene identifier.")
  }
  candidates$replication_logFC <- table$logFC[index]
  candidates$replication_logCPM <- table$logCPM[index]
  candidates$replication_F <- table$F[index]
  candidates$replication_PValue <- table$PValue[index]
  candidates$TAURUS_sign <- sign(candidates$TAURUS_logFC)
  candidates$replication_sign <- sign(candidates$replication_logFC)
  candidates$direction_concordant <- candidates$estimable &
    candidates$TAURUS_sign == candidates$replication_sign
  candidates$family_BH_q <- NA_real_
  finite <- is.finite(candidates$replication_PValue)
  if (any(finite)) candidates$family_BH_q[finite] <- stats::p.adjust(
    candidates$replication_PValue[finite], method = "BH", n = family_size[[target]])
  candidates$primary_replication_status <- ifelse(!candidates$estimable, "NON_ESTIMABLE",
    ifelse(!candidates$direction_concordant, "DIRECTION_DISCORDANT",
      ifelse(!is.na(candidates$family_BH_q) & candidates$family_BH_q <= 0.05,
        "MULTIPLICITY_SUPPORTED_REPLICATION",
        ifelse(!is.na(candidates$replication_PValue) & candidates$replication_PValue < 0.05,
          "NOMINAL_DIRECTIONAL_SUPPORT", "DIRECTION_ONLY_DESCRIPTIVE"))))
  candidates$primary_bulk_feature_eligible <- candidates$primary_replication_status ==
    "MULTIPLICITY_SUPPORTED_REPLICATION"
  candidates$nominal_support_cannot_forward <- candidates$primary_replication_status ==
    "NOMINAL_DIRECTIONAL_SUPPORT"
  candidates$sensitivity_cannot_promote <- TRUE
  candidates$frozen_family_size <- family_size[[target]]
  candidates$multiplicity_method <- "BH_with_original_frozen_family_size"
  list(dge = dge, fit = fit, keep = keep, table = table, candidates = candidates,
       n_native_features = nrow(counts), n_retained_features = nrow(dge$counts))
}

synthesize_replication <- function(candidates) {
  # This retrospective sensitivity does not change the historical forward rule.
  candidates$pooled_13_BH_q <- NA_real_
  finite <- is.finite(candidates$replication_PValue)
  if (any(finite)) candidates$pooled_13_BH_q[finite] <- stats::p.adjust(
    candidates$replication_PValue[finite], method = "BH", n = nrow(candidates))
  candidates$pooled_13_same_direction_pass <- candidates$direction_concordant &
    !is.na(candidates$pooled_13_BH_q) & candidates$pooled_13_BH_q <= 0.05
  family_summary <- do.call(rbind, lapply(split(candidates, candidates$analysis_id), function(data) {
    data.frame(compartment = unique(data$analysis_id), frozen_family_size = unique(data$frozen_family_size),
      n_estimable = sum(data$estimable), n_same_direction = sum(data$direction_concordant),
      n_multiplicity_supported = sum(data$primary_bulk_feature_eligible),
      n_nominal_directional = sum(data$primary_replication_status == "NOMINAL_DIRECTIONAL_SUPPORT"),
      n_direction_only = sum(data$primary_replication_status == "DIRECTION_ONLY_DESCRIPTIVE"),
      n_direction_discordant = sum(data$primary_replication_status == "DIRECTION_DISCORDANT"),
      n_nonestimable = sum(!data$estimable), stringsAsFactors = FALSE)
  }))
  endpoint <- candidates[, c("analysis_id", "gene_id", "gene_symbol", "TAURUS_direction",
    "TAURUS_logFC", "replication_logFC", "replication_PValue", "family_BH_q",
    "direction_concordant"), drop = FALSE]
  endpoint$replication_status <- candidates$primary_replication_status
  endpoint$evidence_tier <- ifelse(candidates$primary_bulk_feature_eligible,
    "REPLICATION_MULTIPLICITY_SUPPORTED", ifelse(
      candidates$primary_replication_status == "NOMINAL_DIRECTIONAL_SUPPORT", "NOMINAL_DIRECTION_ONLY",
      ifelse(candidates$primary_replication_status == "DIRECTION_ONLY_DESCRIPTIVE",
        "DIRECTION_ONLY_DESCRIPTIVE", "DISCORDANT_OR_NON_FORWARD")))
  endpoint$primary_bulk_feature_eligible <- candidates$primary_bulk_feature_eligible
  endpoint$sensitivity_cannot_promote <- TRUE
  endpoint$final_feature_status <- ifelse(candidates$primary_bulk_feature_eligible,
    "FORWARD_TO_GSE193677_BULK_DEVELOPMENT", paste0("FROZEN_NON_FORWARD_", candidates$primary_replication_status))
  forward <- candidates[candidates$primary_bulk_feature_eligible,
    c("analysis_id", "gene_id", "gene_symbol", "TAURUS_direction", "TAURUS_logFC",
      "replication_logFC", "replication_PValue", "family_BH_q", "direction_concordant",
      "primary_replication_status", "multiplicity_method"), drop = FALSE]
  forward$feature_rank <- seq_len(nrow(forward))
  forward$orientation_rule <- rep("positive_logFC_means_UC_gt_CD", nrow(forward))
  forward$feature_status <- rep("FROZEN_PRIMARY_BULK_FEATURE", nrow(forward))
  supplementary <- data.frame(
    compartment = ifelse(candidates$analysis_id == "MACROPHAGE_PRIMARY", "Macrophage", "Fibroblast"),
    gene_symbol = candidates$gene_symbol, ensembl_id = candidates$gene_id,
    discovery_log2FC = candidates$TAURUS_logFC, replication_log2FC = candidates$replication_logFC,
    replication_P = candidates$replication_PValue, BH_compartment = candidates$family_BH_q,
    BH_all_13 = candidates$pooled_13_BH_q, direction_concordant = candidates$direction_concordant,
    stringsAsFactors = FALSE)
  list(candidates = candidates, family_summary = family_summary, endpoint = endpoint,
       primary_freeze = forward, supplementary_table_s2 = supplementary)
}

validate_replication <- function(bundle, project, record) {
  archive <- file.path(project, "reference", "replication_expected.zip")
  if (!identical(unname(tools::md5sum(archive)), "c79e6091da9b344bfcbd37049dd6bcea")) {
    stop("The archived replication comparison tables failed their checksum.")
  }
  expected_files <- c("15_annotation_conflict_audit.tsv", "15_target_cell_support_by_donor.tsv",
    "15_replication_compartment_readiness.tsv", "15_target_cell_membership.tsv",
    "16_primary_candidate_replication.tsv", "16_pseudobulk_library_qc.tsv",
    "16_primary_design_qc.tsv", "16_replication_family_summary.tsv",
    "17_replication_synthesis_endpoint.tsv", "17_primary_bulk_feature_freeze.tsv")
  listing <- utils::unzip(archive, list = TRUE)$Name
  if (anyDuplicated(listing) || !setequal(listing, expected_files)) {
    stop("The replication comparison archive has unexpected members.")
  }
  directory <- tempfile(pattern = "ibd_replication_expected_")
  if (!dir.create(directory)) stop("Cannot create the temporary comparison directory.")
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(archive, files = expected_files, exdir = directory, unzip = "internal")
  read <- function(name) utils::read.delim(file.path(directory, name),
    check.names = FALSE, stringsAsFactors = FALSE)
  outcomes <- logical()
  compare <- function(observed, filename, keys, fields = NULL, exact = character()) {
    expected <- read(filename)
    if (is.null(fields)) fields <- names(expected)
    outcomes <<- c(outcomes, compare_table(observed, expected, keys, filename, record,
                                          fields = fields, exact_fields = exact))
  }
  compare(bundle$annotation$conflicts, expected_files[[1L]], "annotation_row",
          exact = "annotation_row")
  compare(bundle$annotation$support, expected_files[[2L]], c("compartment", "sample_id"),
          exact = "n_target_cells")
  compare(bundle$annotation$readiness, expected_files[[3L]], "compartment",
          exact = c("supported_CD", "supported_UC", "supported_total", "design_rank"))
  compare(bundle$annotation$membership, expected_files[[4L]], c("compartment", "composite_cell_id"),
          exact = "annotation_row")
  compare(bundle$candidates, expected_files[[5L]], c("analysis_id", "gene_id"),
          exact = c("TAURUS_sign", "replication_sign", "frozen_family_size"))
  compare(bundle$library_qc, expected_files[[6L]], c("compartment", "sample_id"),
          exact = names(bundle$library_qc))
  compare(bundle$design_qc, expected_files[[7L]], "design_formula", exact = "design_rank")
  compare(bundle$family_summary, expected_files[[8L]], "compartment",
          exact = names(bundle$family_summary))
  # The old stage-17 process flag is not a measurement in the consolidated run.
  endpoint_fields <- setdiff(names(read(expected_files[[9L]])), "expression_values_read_in_stage17")
  compare(bundle$endpoint, expected_files[[9L]], c("analysis_id", "gene_id"), fields = endpoint_fields)
  freeze_fields <- setdiff(names(read(expected_files[[10L]])), "expression_values_read_in_stage17")
  compare(bundle$primary_freeze, expected_files[[10L]], c("analysis_id", "gene_id"),
          fields = freeze_fields, exact = "feature_rank")
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
  logical_text <- tolower(as.character(final$direction_concordant))
  if (anyNA(logical_text) || any(!logical_text %in% c("true", "false"))) {
    stop("The supplementary reference contains an invalid concordance flag.")
  }
  final$direction_concordant <- logical_text == "true"
  outcomes <- c(outcomes, compare_table(bundle$supplementary_table_s2, final,
    c("compartment", "ensembl_id"), "current_Table_S2_complete", record))
  all(outcomes)
}

run_replication <- function(context) {
  project <- context$project
  work <- file.path(project, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  report_path <- file.path(work, "replication_report.txt")
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
    emit(c("IBD independent single-cell candidate replication report", "Code revision: stage03.0",
      paste("Started UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      "Analysis inputs: deposited GSE214695 raw archive, annotation, and source sample metadata.",
      "Candidate identities and discovery estimates come from the validated, newly computed discovery cache.",
      "Archived output tables are used only for comparison after computation."))
    check_validation_logic()
    record("PASS", "validator_self_checks", "Shared numerical and identifier comparison checks.")
    required <- c(edgeR = "4.4.2", limma = "3.62.2", Matrix = "1.7-1", statmod = "1.5.2")
    for (package in names(required)) {
      if (!requireNamespace(package, quietly = TRUE)) stop("Missing package: ", package)
      observed <- utils::packageVersion(package)
      if (observed != package_version(required[[package]])) {
        stop("Package version differs from the reproduction target: ", package,
             " observed ", as.character(observed), " expected ", required[[package]])
      }
      record("PASS", paste0("package_", package), as.character(observed))
    }
    if (getRversion() != package_version("4.4.2")) {
      record("WARN", "R_version", paste("Observed", as.character(getRversion()),
                                         "; archived analyses used 4.4.2."))
    }
    if (!isTRUE(context$settings$verify_md5)) stop("Replication requires checksum verification.")
    manifest <- context$read_manifest(file.path(project, "config", "inputs.tsv"))
    cache_path <- file.path(work, "discovery.rds")
    if (!file.exists(cache_path)) stop("Run the discovery stage first; work/discovery.rds is required.")
    discovery <- readRDS(cache_path)
    signature_columns <- c("analysis_id", "gene_id", "gene_symbol", "direction", "logFC", "FDR")
    if (!isTRUE(discovery$validation_passed) || !identical(discovery$code_revision, "stage02.0") ||
        !all(signature_columns %in% names(discovery$candidates)) || nrow(discovery$candidates) != 13L ||
        anyDuplicated(paste(discovery$candidates$analysis_id, discovery$candidates$gene_id, sep = "||"))) {
      stop("The discovery cache is incomplete, incompatible, or did not pass validation.")
    }
    for (identifier in c("discovery_myeloid", "discovery_fibperi")) {
      spec <- manifest[manifest$input_id == identifier, , drop = FALSE]
      if (nrow(spec) != 1L || !identical(unname(discovery$source_md5[identifier]), spec$expected_md5[[1L]])) {
        stop("The discovery cache source identities disagree with the input definitions.")
      }
    }
    signature <- discovery$candidates[, signature_columns, drop = FALSE]
    discovery_md5 <- unname(tools::md5sum(cache_path))
    record("PASS", "discovery_cache", paste("Validated stage02.0; newly computed candidates:", nrow(signature),
                                             "; cache MD5:", discovery_md5))
    rm(discovery)
    paths <- source_md5 <- character()
    for (identifier in c("replication_archive", "replication_annotations", "gse214695_sample_metadata")) {
      spec <- manifest[manifest$input_id == identifier, , drop = FALSE]
      if (nrow(spec) != 1L || !context$safe_relative(spec$relative_path[[1L]]) ||
          is.na(spec$expected_md5[[1L]])) stop("Missing or incomplete input definition: ", identifier)
      path <- file.path(context$settings$legacy_root, spec$relative_path[[1L]])
      if (identifier %in% names(context$settings$file_overrides)) {
        path <- context$settings$file_overrides[[identifier]]
      }
      path <- normalizePath(path.expand(path), winslash = "/", mustWork = TRUE)
      observed <- unname(tools::md5sum(path))
      if (file.info(path)$size != as.numeric(spec$expected_bytes[[1L]]) ||
          !identical(observed, spec$expected_md5[[1L]])) {
        stop("Source byte count or checksum differs for ", identifier, ".")
      }
      paths[identifier] <- path
      source_md5[identifier] <- observed
      record("PASS", identifier, paste("Archived source MD5:", observed))
    }
    sample_map <- replication_sample_map()
    metadata <- read_replication_metadata(paths[["gse214695_sample_metadata"]], sample_map, context$read_manifest)
    record("PASS", "source_sample_mapping", "18 unique GSMs; source titles, diagnosis, anatomy, and annotation sample IDs agree.")
    annotation <- prepare_replication_annotation(paths[["replication_annotations"]], metadata, sample_map)
    emit(paste("INFO | source_annotations | Records:", annotation$n_annotation_rows,
      "; conflicting cell identities excluded:", annotation$n_conflicting_identities,
      "; target membership records including controls:", nrow(annotation$membership)))
    pseudobulk <- aggregate_replication(paths[["replication_archive"]], annotation, sample_map, emit)
    record("PASS", "raw_pseudobulk", paste("12 donors; two compartments; native features:",
      nrow(pseudobulk$features), "; integer counts, feature order, barcode identity, and count conservation checked."))
    models <- list()
    for (target in c("MACROPHAGE_PRIMARY", "FIBROBLAST_PRIMARY")) {
      emit(paste("INFO | Fitting", target, "| design: ~ broad_region + disease"))
      models[[target]] <- fit_replication_target(target, pseudobulk, signature)
      emit(paste("INFO |", target, "| retained genes:", models[[target]]$n_retained_features))
    }
    synthesis <- synthesize_replication(do.call(rbind, lapply(models, function(model) model$candidates)))
    design_qc <- data.frame(design_formula = "~ broad_region + disease",
      design_columns = paste(colnames(pseudobulk$design), collapse = "|"),
      design_rank = qr(pseudobulk$design)$rank,
      design_full_rank = qr(pseudobulk$design)$rank == ncol(pseudobulk$design), stringsAsFactors = FALSE)
    bundle <- c(list(code_revision = "stage03.0", source_md5 = source_md5,
      discovery_cache_md5 = discovery_md5, created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      annotation = annotation, pseudobulk = pseudobulk, models = models,
      library_qc = pseudobulk$library_qc, design_qc = design_qc,
      settings = list(min_cells = 10L, design = "~ broad_region + disease",
        reference_disease = "CD", reference_region = "Distal",
        filter = list(min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.70),
        ql_robust = TRUE, primary_families = c(macrophage = 3L, fibroblast = 10L),
        pooled_sensitivity_family = 13L, primary_q_cutoff = 0.05,
        numerical_atol = 1e-10, probability_atol = 1e-14, numerical_rtol = 1e-6),
      session = utils::sessionInfo()), synthesis)
    emit(c("", "Comparison with archived results and the complete current supplementary Table S2"))
    passed <- validate_replication(bundle, project, record)
    bundle$validation_passed <- isTRUE(passed) && counters[["FAIL"]] == 0L
    bundle$validation_counts <- counters
    cache <- file.path(work, "replication.rds")
    saveRDS(bundle, cache, compress = "xz")
    emit(c("", "Newly computed candidate replication and pooled-family sensitivity",
      utils::capture.output(print(bundle$supplementary_table_s2, row.names = FALSE, digits = 10)),
      paste("Estimable candidates:", sum(bundle$candidates$estimable)),
      paste("Direction-concordant candidates:", sum(bundle$candidates$direction_concordant)),
      paste("Compartment-family rule passing candidates:", sum(bundle$candidates$primary_bulk_feature_eligible)),
      paste("Primary forward genes:", if (nrow(bundle$primary_freeze)) {
        paste(bundle$primary_freeze$gene_symbol, collapse = ", ")
      } else "NONE"),
      paste("Pooled 13-candidate same-direction passing candidates:",
            sum(bundle$candidates$pooled_13_same_direction_pass)),
      "The pooled sensitivity does not change the historical compartment-family forward rule.",
      "The available historical files do not independently establish prospective specification of the families.",
      paste("Cache:", normalizePath(cache, winslash = "/", mustWork = TRUE)),
      "Later stages must check validation_passed; a fatal run does not refresh the previous cache.",
      "", "R session information", utils::capture.output(utils::sessionInfo())))
    status <- if (isTRUE(bundle$validation_passed)) 0L else 1L
  }, warning = function(warning) {
    record("WARN", "R_warning", conditionMessage(warning))
    invokeRestart("muffleWarning")
  }), error = function(error) {
    record("FAIL", "replication_error", conditionMessage(error))
    status <<- 1L
  })
  emit(c("", paste("Checks: PASS =", counters[["PASS"]], "; WARN =", counters[["WARN"]],
    "; FAIL =", counters[["FAIL"]]),
    paste("Replication reproduction:", if (status == 0L) "PASS" else "BLOCKED"),
    "Bulk analyses, final workbook, and figures: NOT_RUN",
    paste("Finished UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Report:", normalizePath(report_path, winslash = "/", mustWork = TRUE))))
  status
}
