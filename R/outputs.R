# Assemble the current publication tables and figures from validated R results.
# Statistical calculations and computational figures remain in R.
# Figure 1 is the existing AI-generated graphical abstract, outside this stage.

output_one <- function(data, column, value) {
  rows <- which(!is.na(data[[column]]) & data[[column]] == value)
  if (length(rows) != 1L) stop("Expected one publication row for ", column, " = ", value, ".")
  data[rows, , drop = FALSE]
}

output_tables <- function(development, replication, external, sensitivity, mapping) {
  tables <- list(development$supplementary$table_s1, replication$supplementary_table_s2,
    development$supplementary$table_s3a, development$supplementary$table_s3b,
    external$supplementary$table_s4, sensitivity$supplementary$table_s5,
    external$supplementary$table_s6)
  names(tables) <- mapping$source_file
  for (i in seq_along(tables)) {
    table <- tables[[i]]
    fields <- strsplit(mapping$source_columns[[i]], "|", fixed = TRUE)[[1L]]
    if (!is.data.frame(table) || !identical(names(table), fields) ||
        nrow(table) != as.integer(mapping$expected_rows[[i]]) ||
        ncol(table) != as.integer(mapping$expected_columns[[i]])) {
      stop("Computed supplementary schema differs: ", mapping$source_file[[i]], ".")
    }
    for (field in names(table)) {
      if (is.factor(table[[field]])) table[[field]] <- as.character(table[[field]])
      if (!is.atomic(table[[field]]) ||
          (is.numeric(table[[field]]) && any(!is.finite(table[[field]][!is.na(table[[field]])])))) {
        stop("Unsupported or nonfinite workbook column: ", field, ".")
      }
    }
    tables[[i]] <- table
  }
  tables
}

output_reference_flags <- function(value) {
  # Apply the same strict True/False interpretation used by validated replication.
  text <- tolower(as.character(value))
  if (!length(text) || anyNA(text) || any(!text %in% c("true", "false"))) {
    stop("The supplementary reference contains an invalid concordance flag.")
  }
  text == "true"
}

output_reference_flag_self_check <- function() {
  canonical <- c(TRUE, FALSE, TRUE)
  spellings <- list(c("True", "False", "True"), c("TRUE", "FALSE", "TRUE"),
    c("true", "false", "true"), canonical)
  for (value in spellings) {
    if (!identical(output_reference_flags(value), canonical)) stop("Reference flag interpretation failed.")
  }
  invalid <- list(c("True", NA_character_), c("True", "yes"), c("True", "1"),
    c("True", ""), character())
  for (value in invalid) {
    rejected <- tryCatch({ output_reference_flags(value); FALSE }, error = function(error) TRUE)
    if (!isTRUE(rejected)) stop("An invalid reference flag was accepted.")
  }
  observed <- data.frame(key = c("A", "B", "C"), direction_concordant = canonical,
    stringsAsFactors = FALSE)
  expected <- data.frame(key = observed$key,
    direction_concordant = output_reference_flags(spellings[[1L]]), stringsAsFactors = FALSE)
  quiet <- function(...) NULL
  if (!compare_table(observed, expected, "key", "reference_flag_fixture", quiet)) {
    stop("Canonical reference flags did not pass the actual comparison helper.")
  }
  expected$direction_concordant[[1L]] <- FALSE
  if (compare_table(observed, expected, "key", "changed_flag_fixture", quiet)) {
    stop("The actual comparison helper accepted a changed concordance value.")
  }
  invisible(TRUE)
}

output_reference_check <- function(tables, mapping, project, record) {
  archive <- file.path(project, "reference", "tables_reference.zip")
  if (!identical(unname(tools::md5sum(archive)), unique(mapping$reference_zip_md5))) {
    stop("The author-provided supplementary reference archive failed its checksum.")
  }
  members <- utils::unzip(archive, list = TRUE)$Name
  members <- members[!endsWith(members, "/")]
  if (anyDuplicated(members) || !setequal(members, mapping$source_file)) {
    stop("The supplementary reference members differ from the mapping.")
  }
  directory <- tempfile("ibd_output_reference_")
  if (!dir.create(directory)) stop("Cannot create reference-check storage.")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  utils::unzip(archive, files = mapping$source_file, exdir = directory, unzip = "internal")
  keys <- list("ensembl_id", c("compartment", "ensembl_id"), "feature", "endpoint",
    c("cohort", "comparator"), c("cohort", "scale"), "group")
  exact <- c("n", "n_uc", "n_cd", "n_UC", "n_CD", "n_analyzed", "genes_more_extreme")
  passed <- logical(length(tables))
  for (i in seq_along(tables)) {
    observed <- tables[[i]]
    # Match published display precision only in the comparison copy.
    if (i == 3L) observed <- development_archive_view(observed,
      decimal = c(auc = 3L, transcriptome_percentile = 1L, adjusted_odds_ratio = 3L))
    if (i == 4L) observed <- development_archive_view(observed, decimal = c(estimate = 3L))
    expected <- utils::read.delim(file.path(directory, mapping$source_file[[i]]),
      sep = "\t", header = TRUE, quote = "\"", comment.char = "", check.names = FALSE,
      stringsAsFactors = FALSE, na.strings = "NA")
    if (mapping$source_file[[i]] == "Table_S2_candidate_replication.tsv") {
      if (!is.logical(observed$direction_concordant) || anyNA(observed$direction_concordant)) {
        stop("Computed concordance flags must remain complete logical values.")
      }
      expected$direction_concordant <- output_reference_flags(expected$direction_concordant)
    }
    passed[[i]] <- compare_table(observed, expected, keys[[i]],
      paste0("publication_reference_", mapping$table[[i]], "_", mapping$part[[i]]),
      record, exact_fields = intersect(exact, names(expected)))
  }
  if (!all(passed)) stop("Final supplementary reference validation failed; outputs were not promoted.")
}

output_xml_nodes <- function(text, tag) {
  regmatches(text, gregexpr(paste0("<", tag, "\\b[^>]*/>"), text, perl = TRUE))[[1L]]
}

output_xml_attribute <- function(node, name) {
  value <- regmatches(node, regexec(paste0("\\b", name, '="([^"]*)"'), node, perl = TRUE))[[1L]]
  if (length(value) != 2L) stop("Missing generated workbook XML attribute: ", name, ".")
  value[[2L]]
}

output_package_target <- function(base, target) {
  path <- if (startsWith(target, "/")) substring(target, 2L) else paste0(base, target)
  pieces <- strsplit(path, "/", fixed = TRUE)[[1L]]
  result <- character()
  for (piece in pieces) {
    if (piece %in% c("", ".")) next
    if (piece == "..") {
      if (!length(result)) stop("A workbook relationship escapes its package root.")
      result <- result[-length(result)]
    } else result <- c(result, piece)
  }
  paste(result, collapse = "/")
}

output_workbook_package <- function(path, dimensions) {
  # This table-only writer uses pinned openxlsx 4.2.8.1. Correct its stale A1
  # dimensions and unused drawing relationships without changing cell records.
  members <- utils::unzip(path, list = TRUE)$Name
  members <- members[!endsWith(members, "/")]
  safe <- function(value) !grepl("^(/|[A-Za-z]:)", value) &&
    !any(strsplit(value, "/", fixed = TRUE)[[1L]] %in% c("", ".."))
  if (anyDuplicated(members) || !all(vapply(members, safe, logical(1))) ||
      any(grepl("^xl/(drawings|charts|media)/", members))) stop("Unexpected content in the generated table-only workbook.")
  directory <- tempfile("ibd_workbook_package_")
  repacked <- tempfile(fileext = ".xlsx")
  if (!dir.create(directory)) stop("Cannot create workbook package verification storage.")
  on.exit(unlink(c(directory, repacked), recursive = TRUE), add = TRUE)
  utils::unzip(path, exdir = directory, unzip = "internal")
  read <- function(member) paste(readLines(file.path(directory, member), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  write <- function(text, member) writeLines(text, file.path(directory, member), useBytes = TRUE)
  for (member in names(dimensions)) {
    text <- read(member)
    if (grepl("<(drawing|legacyDrawing|legacyDrawingHF)\\b", text, perl = TRUE)) stop("Unexpected drawing element in a table-only worksheet.")
    nodes <- output_xml_nodes(text, "dimension")
    if (length(nodes) != 1L) stop("Generated worksheet lacks a unique dimension declaration.")
    corrected <- paste0('<dimension ref="', dimensions[[member]], '"/>')
    updated <- sub(nodes[[1L]], corrected, text, fixed = TRUE)
    # Replacing this one metadata node must leave all worksheet cell XML intact.
    if (!identical(sub(nodes[[1L]], "", text, fixed = TRUE), sub(corrected, "", updated, fixed = TRUE))) {
      stop("Worksheet content changed outside its dimension declaration.")
    }
    write(updated, member)
  }
  relationships <- members[endsWith(members, ".rels")]
  for (member in relationships) {
    text <- read(member)
    nodes <- output_xml_nodes(text, "Relationship")
    base <- if (member == "_rels/.rels") "" else sub("_rels/[^/]+[.]rels$", "", member)
    for (node in nodes) {
      type <- output_xml_attribute(node, "Type")
      if (grepl("/(drawing|vmlDrawing)$", type)) {
        target <- output_xml_attribute(node, "Target")
        resolved <- output_package_target(base, target)
        if (!grepl("^xl/drawings/(drawing[0-9]+[.]xml|vmlDrawing[0-9]+[.]vml)$", resolved) || resolved %in% members) {
          stop("Unexpected drawing relationship in the generated workbook.")
        }
        text <- sub(node, "", text, fixed = TRUE)
      }
    }
    write(text, member)
  }
  types <- read("[Content_Types].xml")
  for (node in output_xml_nodes(types, "Override")) {
    target <- output_xml_attribute(node, "PartName")
    if (grepl("^/xl/drawings/drawing[0-9]+[.]xml$", target) && !substring(target, 2L) %in% members) {
      types <- sub(node, "", types, fixed = TRUE)
    }
  }
  write(types, "[Content_Types].xml")
  # Every retained internal relationship and content-type override must resolve.
  for (member in relationships) {
    nodes <- output_xml_nodes(read(member), "Relationship")
    ids <- vapply(nodes, output_xml_attribute, character(1), name = "Id")
    if (anyDuplicated(ids)) stop("Duplicate workbook relationship identifiers.")
    base <- if (member == "_rels/.rels") "" else sub("_rels/[^/]+[.]rels$", "", member)
    for (node in nodes) {
      if (grepl('TargetMode="External"', node, fixed = TRUE)) next
      target <- output_package_target(base, output_xml_attribute(node, "Target"))
      if (!target %in% members) stop("Unresolved workbook relationship: ", target, ".")
    }
  }
  for (node in output_xml_nodes(types, "Override")) {
    target <- substring(output_xml_attribute(node, "PartName"), 2L)
    if (!target %in% members) stop("Unresolved workbook content-type declaration: ", target, ".")
  }
  zip::zipr(repacked, files = members, root = directory, include_directories = FALSE, mode = "mirror")
  restored <- utils::unzip(repacked, list = TRUE)$Name
  if (anyDuplicated(restored) || !setequal(restored, members)) stop("Workbook repacking changed package membership.")
  if (!file.copy(repacked, path, overwrite = TRUE)) stop("Cannot retain the structurally checked workbook.")
  invisible(TRUE)
}

output_workbook <- function(path, tables, mapping) {
  wb <- openxlsx::createWorkbook(creator = "IBD publication workflow",
    title = "IBD supplementary tables", subject = "Computed R results; six worksheets")
  openxlsx::modifyBaseFont(wb, fontName = "Calibri", fontSize = 11)
  sheets <- paste0("Table_S", seq_len(6L))
  titles <- c("S1. Gene panel and analytical roles", "S2. Discovery candidates and independent replication",
    "S3A. Development comparisons and adjusted associations", "S3B. IGHG3 development endpoints",
    "S4. Paired AUC multiplicity sensitivity", "S5. Descriptive component decomposition",
    "S6. External phenotype and control summaries")
  for (sheet in sheets) {
    openxlsx::addWorksheet(wb, sheet, gridLines = FALSE)
    openxlsx::freezePane(wb, sheet, firstActiveRow = 5L)
    openxlsx::pageSetup(wb, sheet, orientation = "landscape", fitToWidth = 1L, fitToHeight = 0L)
  }
  heading <- openxlsx::createStyle(fontSize = 14, fontColour = "#17384B", textDecoration = "bold", wrapText = TRUE)
  note <- openxlsx::createStyle(fontSize = 10, fontColour = "#536572", wrapText = TRUE)
  header <- openxlsx::createStyle(fgFill = "#17384B", fontColour = "#FFFFFF",
    textDecoration = "bold", wrapText = TRUE, valign = "center", border = "bottom")
  body <- openxlsx::createStyle(valign = "top", wrapText = TRUE)
  next_row <- stats::setNames(rep(1L, length(sheets)), sheets)
  widths_by_sheet <- lapply(sheets, function(x) numeric())
  names(widths_by_sheet) <- sheets
  for (i in seq_along(tables)) {
    table <- tables[[i]]
    sheet <- mapping$worksheet[[i]]
    row <- next_row[[sheet]]
    last_column <- ncol(table)
    openxlsx::writeData(wb, sheet, titles[[i]], startRow = row, colNames = FALSE)
    openxlsx::mergeCells(wb, sheet, cols = seq_len(last_column), rows = row)
    openxlsx::addStyle(wb, sheet, heading, rows = row, cols = 1L)
    openxlsx::setRowHeights(wb, sheet, row, 28)
    openxlsx::writeData(wb, sheet, paste("Source comparison:", mapping$source_file[[i]],
      "| Numerical cells retain computed precision."), startRow = row + 1L, colNames = FALSE)
    openxlsx::mergeCells(wb, sheet, cols = seq_len(last_column), rows = row + 1L)
    openxlsx::addStyle(wb, sheet, note, rows = row + 1L, cols = 1L)
    openxlsx::setRowHeights(wb, sheet, row + 1L, 28)
    header_row <- row + 3L
    data_rows <- seq.int(header_row + 1L, header_row + nrow(table))
    region <- tools::file_path_sans_ext(mapping$source_file[[i]])
    openxlsx::writeData(wb, sheet, table, startRow = header_row, startCol = 1L,
      colNames = TRUE, rowNames = FALSE, keepNA = TRUE, na.string = "NA", name = region,
      headerStyle = header, withFilter = FALSE)
    openxlsx::setRowHeights(wb, sheet, header_row, 65)
    openxlsx::addStyle(wb, sheet, body, rows = data_rows, cols = seq_len(last_column), gridExpand = TRUE)
    widths <- vapply(seq_along(table), function(j) {
      field <- names(table)[[j]]
      if (is.numeric(table[[j]])) return(if (nchar(field) > 25L) 24 else 19)
      if (is.logical(table[[j]])) return(21)
      longest <- max(nchar(c(field, as.character(table[[j]]))), na.rm = TRUE)
      min(40, max(23, longest * 0.85))
    }, numeric(1))
    previous <- widths_by_sheet[[sheet]]
    length(previous) <- max(length(previous), length(widths))
    previous[is.na(previous)] <- 0
    previous[seq_along(widths)] <- pmax(previous[seq_along(widths)], widths)
    widths_by_sheet[[sheet]] <- previous
    for (j in seq_along(table)) {
      if (is.numeric(table[[j]])) {
        field <- names(table)[[j]]
        count <- field %in% c("n", "n_uc", "n_cd", "n_UC", "n_CD", "n_analyzed", "genes_more_extreme")
        probability <- grepl("(^BH_|^p_|_P$)", field)
        format <- if (count) "0" else if (probability) "0.000E+00" else
          if (field == "transcriptome_percentile") "0.0" else "0.000"
        style <- openxlsx::createStyle(numFmt = format, valign = "top")
        openxlsx::addStyle(wb, sheet, style, rows = data_rows, cols = j, stack = TRUE, gridExpand = TRUE)
      }
    }
    next_row[[sheet]] <- max(data_rows) + 4L
  }
  for (sheet in sheets) {
    widths <- widths_by_sheet[[sheet]]
    openxlsx::setColWidths(wb, sheet, cols = seq_along(widths), widths = widths)
    # Shared S3 column widths are fixed before calculating wrapped row heights.
    for (i in which(mapping$worksheet == sheet)) {
      table <- tables[[i]]
      start <- if (mapping$part[[i]] == "B") nrow(tables[[i - 1L]]) + 12L else 5L
      for (r in seq_len(nrow(table))) {
        lines <- vapply(seq_along(table), function(j) {
          value <- as.character(table[[j]][[r]])
          if (is.na(value) || is.numeric(table[[j]]) || is.logical(table[[j]])) return(1)
          max(1, ceiling(nchar(value) / max(8, widths[[j]] - 4)))
        }, numeric(1))
        openxlsx::setRowHeights(wb, sheet, start + r - 1L, max(22, 16 * max(lines) + 8))
      }
    }
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  dimensions <- stats::setNames(vapply(sheets, function(sheet) {
    paste0("A1:", openxlsx::int2col(length(widths_by_sheet[[sheet]])), next_row[[sheet]] - 4L)
  }, character(1)), paste0("xl/worksheets/sheet", seq_along(sheets), ".xml"))
  output_workbook_package(path, dimensions)
  invisible(path)
}

output_roundtrip <- function(path, tables, mapping) {
  if (!identical(openxlsx::getSheetNames(path), paste0("Table_S", seq_len(6L)))) {
    stop("The saved workbook does not have exactly the six ordered supplementary sheets.")
  }
  regions <- tools::file_path_sans_ext(mapping$source_file)
  actual_regions <- as.character(openxlsx::getNamedRegions(path))
  if (length(actual_regions) != length(regions) || anyDuplicated(actual_regions) ||
      !setequal(actual_regions, regions)) stop("Workbook blocks are missing or duplicated.")
  largest <- 0
  exact <- c("n", "n_uc", "n_cd", "n_UC", "n_CD", "n_analyzed", "genes_more_extreme")
  for (i in seq_along(tables)) {
    actual <- openxlsx::read.xlsx(path, namedRegion = regions[[i]], colNames = TRUE,
      check.names = FALSE, skipEmptyRows = FALSE, skipEmptyCols = FALSE,
      detectDates = FALSE, na.strings = "NA")
    expected <- tables[[i]]
    if (!identical(dim(actual), dim(expected)) || !identical(names(actual), names(expected))) {
      stop("Saved workbook dimensions or columns differ: ", regions[[i]], ".")
    }
    for (field in names(expected)) {
      a <- actual[[field]]
      e <- expected[[field]]
      if (!identical(is.na(a), is.na(e))) stop("Workbook missingness differs: ", field, ".")
      present <- !is.na(e)
      if (!any(present)) next
      if (is.numeric(e)) {
        if (!is.numeric(a)) stop("Workbook numeric type lost: ", field, ".")
        difference <- abs(a[present] - e[present])
        largest <- max(largest, difference)
        tolerance <- if (grepl("(^BH_|^p_|_P$)", field)) 1e-14 else 1e-12
        if (any(!is.finite(a[present])) || any(difference > tolerance + 1e-12 * abs(e[present])) ||
            (field %in% exact && any(difference != 0))) stop("Workbook numerical round-trip failed: ", field, ".")
      } else {
        if (is.character(e) && !is.character(a)) stop("Workbook text type lost: ", field, ".")
        if (is.logical(e) && !is.logical(a)) stop("Workbook logical type lost: ", field, ".")
        if (!identical(a[present], e[present])) stop("Workbook exact values differ: ", field, ".")
      }
    }
  }
  largest
}

output_workbook_self_check <- function() {
  mapping <- data.frame(source_file = paste0("Fixture_", seq_len(7L), ".tsv"),
    worksheet = c("Table_S1", "Table_S2", "Table_S3", "Table_S3", "Table_S4", "Table_S5", "Table_S6"),
    part = c("main", "main", "A", "B", "main", "main", "main"), stringsAsFactors = FALSE)
  fixture <- data.frame(identifier = c("001", "002", "003"), value = c(1 / 3, 1e-12, NA_real_),
    flag = c(TRUE, FALSE, NA), stringsAsFactors = FALSE)
  tables <- rep(list(fixture), 7L)
  tables[[4L]] <- data.frame(endpoint = c("first", "second"), estimate = c(-1.25, 2.5),
    n_analyzed = c(10L, 12L), note = c("Independent S3B columns", NA_character_), stringsAsFactors = FALSE)
  path <- tempfile(fileext = ".xlsx")
  on.exit(unlink(path), add = TRUE)
  output_workbook(path, tables, mapping)
  output_roundtrip(path, tables, mapping)
  invisible(TRUE)
}

output_plot_data <- function(discovery, replication, development, external, sensitivity, record) {
  components <- sensitivity$components$participants
  cohorts <- c("GSE193677", "GSE117993", "GSE111889")
  if (nrow(components) != 945L || anyDuplicated(paste(components$cohort, components$participant_id)) ||
      any(!is.finite(components$score)) || max(abs(components$score - components$score_saved)) > 1e-10) {
    stop("Invalid component-derived publication populations or scores.")
  }
  aucs <- external$comparisons$all_aucs
  paired <- external$comparisons$all_paired
  if (nrow(aucs) != 15L || nrow(paired) != 12L) stop("Incomplete comparator publication tables.")
  primary <- aucs[aucs$feature == "primary_ratio", , drop = FALSE]
  primary <- primary[match(cohorts, primary$cohort), , drop = FALSE]
  for (cohort in cohorts) {
    rows <- components[components$cohort == cohort & components$diagnosis %in% c("UC", "CD"), , drop = FALSE]
    expected <- output_one(primary, "cohort", cohort)
    actual <- uc_auc_summary(rows$score, as.integer(rows$diagnosis == "UC"))
    if (!compare_table(actual, expected, c("n_uc", "n_cd"), paste0("figure_primary_auc_", cohort), record,
      fields = names(actual), exact_fields = c("n_uc", "n_cd"))) stop("Figure score AUC disagrees with validated results.")
  }
  contrasts <- external$descriptions$subtype_contrasts
  subtypes <- c("ILEAL_CD", "COLONIC_CD")
  subtype_rows <- do.call(rbind, lapply(subtypes, function(group) {
    contrast <- output_one(contrasts, "negative_group", group)
    data <- external$external$all
    data <- data[data$group %in% c("UC", group), , drop = FALSE]
    summary <- uc_auc_summary(data$IGG_IGAM_LOG_RATIO, as.integer(data$group == "UC"))
    if (max(abs(summary$auc - contrast$auc), abs(summary$ci95_low - contrast$ci_low),
      abs(summary$ci95_high - contrast$ci_high)) > 1e-10 || summary$n_uc != contrast$n_positive ||
      summary$n_cd != contrast$n_negative) stop("A phenotype forest row differs from its computed source.")
    data.frame(cohort = "GSE117993", role = "post hoc CD phenotype comparison", feature = "primary_ratio", summary)
  }))
  forest <- rbind(primary[1:2, ], subtype_rows, primary[3L, ])
  forest$label <- c("Development: UC vs CD", "External: UC vs CD", "External: UC vs CD L3",
    "External: UC vs CD L2", "HMP2: UC vs CD (context)")
  background <- development$background$all_genes
  if (nrow(background) != 18940L || anyDuplicated(background$raw_gene_id) ||
      anyNA(background$auc_uc_positive) || any(background$auc_uc_positive < 0 | background$auc_uc_positive > 1)) {
    stop("Invalid observed-label transcriptome background for Figure S1.")
  }
  context <- development$associations$feature_context
  feature_order <- c("IGG_IGAM_LOG_RATIO", "IGHG3_FRACTION_OF_IG", "IGHG3", "TOTAL_IG", "PLASMA_CELL_SCORE")
  context <- context[match(feature_order, context$feature), , drop = FALSE]
  if (anyNA(context$feature)) stop("A required Figure S1 feature is absent.")
  macro <- discovery$models$MACROPHAGE_PRIMARY$donor_metadata
  fibro <- discovery$models$FIBROBLAST_PRIMARY$donor_metadata
  list(primary = primary, forest = forest, paired = paired, aucs = aucs,
    external = external$external$all[, c("sample_id", "diagnosis", "group", "IGG_IGAM_LOG_RATIO")],
    components = components, background = background, feature_context = context,
    replication_geometry = table(replication$pseudobulk$metadata$disease),
    development_biopsies = nrow(development$selection), context_rectal_samples = nrow(external$context$rectal),
    discovery_geometry = rbind(data.frame(compartment = "Macrophage", n_uc = sum(macro$disease == "UC"), n_cd = sum(macro$disease == "CD")),
      data.frame(compartment = "Fibroblast", n_uc = sum(fibro$disease == "UC"), n_cd = sum(fibro$disease == "CD"))))
}

output_theme <- function(mar = c(4.5, 5, 3.5, 1.5), oma = c(0, 0, 0, 0)) {
  graphics::par(mar = mar, oma = oma, family = "sans", las = 1, bty = "l",
    col.axis = "#34434B", col.lab = "#233746", col.main = "#17384B",
    fg = "#34434B", cex.axis = 0.95, cex.lab = 1.05, font.main = 2,
    mgp = c(2.8, 0.7, 0), tcl = -0.25, cex = 1)
}

output_jitter <- function(ids, amplitude) {
  # Deterministic offsets preserve all points without changing R's RNG state.
  if (anyNA(ids) || anyDuplicated(ids)) stop("Jitter identities must be unique and nonmissing.")
  if (length(amplitude) != 1L || !is.finite(amplitude) || amplitude < 0) stop("Invalid jitter amplitude.")
  rank <- match(as.character(ids), sort(as.character(ids)))
  # R gives %% higher precedence than multiplication; enclose the product first.
  fraction <- (rank * ((sqrt(5) - 1) / 2)) %% 1
  offsets <- (2 * fraction - 1) * amplitude
  if (any(!is.finite(offsets)) || any(abs(offsets) > amplitude + 1e-12)) stop("Jitter escaped its declared group width.")
  offsets
}

output_jitter_self_check <- function() {
  ids <- sprintf("fixture_%04d", seq_len(1000L))
  forward <- output_jitter(ids, 0.2)
  backward <- output_jitter(rev(ids), 0.2)
  if (length(forward) != length(ids) || any(abs(forward) > 0.2) ||
      !any(forward < 0) || !any(forward > 0) || !identical(forward, rev(backward))) {
    stop("Deterministic jitter bounds or identifier alignment failed.")
  }
  invisible(TRUE)
}

output_box <- function(values, ids, x, colour, width = 0.5, point_cex = 0.55) {
  if (!length(values) || length(values) != length(ids) || any(!is.finite(values))) stop("Invalid boxplot values.")
  quartile <- as.numeric(stats::quantile(values, c(0.25, 0.5, 0.75), type = 7L))
  iqr <- quartile[[3L]] - quartile[[1L]]
  inside <- values >= quartile[[1L]] - 1.5 * iqr & values <= quartile[[3L]] + 1.5 * iqr
  whiskers <- range(values[inside])
  graphics::rect(x - width / 2, quartile[[1L]], x + width / 2, quartile[[3L]],
    col = grDevices::adjustcolor(colour, alpha.f = 0.30), border = colour, lwd = 1.5)
  graphics::segments(x, whiskers[[1L]], x, quartile[[1L]], col = colour)
  graphics::segments(x, quartile[[3L]], x, whiskers[[2L]], col = colour)
  graphics::segments(x - width / 4, whiskers, x + width / 4, whiskers, col = colour)
  graphics::points(x + output_jitter(ids, width * 0.36), values, pch = 16,
    col = grDevices::adjustcolor(colour, alpha.f = 0.42), cex = point_cex)
  graphics::segments(x - width / 2, quartile[[2L]], x + width / 2, quartile[[2L]], col = "#152C3B", lwd = 2)
}

output_forest <- function(data) {
  rows <- data$forest
  output_theme(mar = c(5.2, 14.5, 5.0, 1.0))
  graphics::plot(NA, xlim = c(0.48, 1.53), ylim = c(0.35, 6.1), xlab = "UC-positive AUC (95% CI)",
    ylab = "", yaxt = "n", xaxt = "n", bty = "n", xaxs = "i", yaxs = "i")
  graphics::axis(1, at = seq(0.5, 1, 0.1), labels = sprintf("%.1f", seq(0.5, 1, 0.1)))
  y <- rev(seq_len(nrow(rows)))
  graphics::axis(2, at = y, labels = rows$label, tick = FALSE, las = 1, cex.axis = 1)
  graphics::segments(0.5, 0.6, 0.5, 5.5, lty = 2, col = "#94A6AD")
  graphics::segments(0.5, y - 0.45, 1.5, y - 0.45, col = "#EDF1F3")
  graphics::segments(rows$ci95_low, y, rows$ci95_high, y, lwd = 2, col = "#34434B")
  graphics::segments(rows$ci95_low, y - 0.09, rows$ci95_low, y + 0.09, col = "#34434B")
  graphics::segments(rows$ci95_high, y - 0.09, rows$ci95_high, y + 0.09, col = "#34434B")
  graphics::points(rows$auc, y, pch = c(15, 18, 16, 16, 17), cex = 1.35, col = "#17384B")
  graphics::text(1.115, 5.72, "UC / CD", font = 2, cex = 0.94)
  graphics::text(1.355, 5.72, "AUC (95% CI)", font = 2, cex = 0.94)
  graphics::text(1.115, y, paste(rows$n_uc, rows$n_cd, sep = " / "), cex = 0.97)
  graphics::text(1.355, y, sprintf("%.3f (%.3f-%.3f)", rows$auc, rows$ci95_low, rows$ci95_high), cex = 0.92)
  graphics::title("Discrimination across cohorts and CD phenotypes", line = 2.4, cex.main = 1.25)
  graphics::mtext("All samples are rectal; CD L3 and L2 are source phenotype labels. HMP2 is context only.",
    side = 1, line = 3.7, cex = 0.8)
}

output_external_distribution <- function(data) {
  table <- data$external
  groups <- c("CONTROL", "ILEAL_CD", "COLONIC_CD", "UC")
  colours <- c("#87949C", "#5C8CCC", "#3F948C", "#D36B2D")
  if (nrow(table) != 190L || any(!table$group %in% groups) || anyNA(table$group)) stop("External figure groups are incomplete.")
  output_theme(mar = c(6.5, 7.2, 4.3, 1.5))
  limits <- range(table$IGG_IGAM_LOG_RATIO)
  span <- diff(limits)
  graphics::plot(NA, xlim = c(0.45, 4.55), ylim = limits + c(-0.1, 0.14) * span,
    xaxt = "n", xlab = "", ylab = expression(log[2] * frac(IgG + 0.5, IgA + IgM + 0.5)))
  graphics::abline(h = graphics::axTicks(2), col = "#EEF2F4", lwd = 1)
  labels <- c("Control", "CD L3", "CD L2", "UC")
  for (i in seq_along(groups)) {
    rows <- table[table$group == groups[[i]], , drop = FALSE]
    output_box(rows$IGG_IGAM_LOG_RATIO, rows$sample_id, i, colours[[i]], width = 0.58, point_cex = 0.7)
    graphics::text(i, max(rows$IGG_IGAM_LOG_RATIO) + 0.055 * span, paste0("n = ", nrow(rows)), cex = 0.95)
  }
  graphics::axis(1, at = seq_along(groups), labels = labels, tick = FALSE, cex.axis = 1.1)
  graphics::title("Rectal score distributions in GSE117993", line = 2, cex.main = 1.3)
  graphics::mtext("CD L3: source ileal CD; CD L2: source colonic CD. All biopsies are rectal.", side = 1, line = 3, cex = 0.82)
  graphics::mtext("Boxes: type-7 quartiles and median; whiskers: 1.5 IQR; points: all participants. Descriptive only.",
    side = 1, line = 4.4, cex = 0.79)
}

output_paired <- function(data) {
  cohorts <- c("GSE193677", "GSE117993", "GSE111889")
  features <- c("IGHG3_abundance", "total_Ig_abundance", "IGHG3_fraction", "summed_IgG_abundance")
  labels <- c("IGHG3 abundance", "Total Ig abundance", "IGHG3 / other Ig", "Summed IgG abundance")
  titles <- c("A. Development: exploratory", "B. External: GSE117993", "C. HMP2: context only")
  limits <- range(c(data$paired$ci95_low, data$paired$ci95_high, 0))
  span <- diff(limits)
  limits <- limits + c(-0.1, 0.27) * span
  graphics::par(mfrow = c(3, 1))
  output_theme(mar = c(3.4, 12.7, 3, 1), oma = c(3, 0, 2, 0))
  for (i in seq_along(cohorts)) {
    rows <- data$paired[data$paired$cohort == cohorts[[i]], , drop = FALSE]
    rows <- rows[match(features, rows$comparator), , drop = FALSE]
    if (anyNA(rows$comparator)) stop("A paired-comparison figure row is absent.")
    y <- 4:1
    graphics::plot(NA, xlim = limits, ylim = c(0.5, 4.6), xlab = "", ylab = "", yaxt = "n", bty = "n")
    graphics::abline(v = 0, col = "#7E929D", lty = 2)
    graphics::axis(2, at = y, labels = labels, tick = FALSE, cex.axis = 0.98)
    colour <- c(rep("#34434B", 3L), "#C3672B")
    graphics::segments(rows$ci95_low, y, rows$ci95_high, y, col = colour, lwd = 2)
    graphics::segments(rows$ci95_low, y - 0.10, rows$ci95_low, y + 0.10, col = colour)
    graphics::segments(rows$ci95_high, y - 0.10, rows$ci95_high, y + 0.10, col = colour)
    graphics::points(rows$delta_primary_minus_comparator, y, pch = 16, cex = 1.15, col = colour)
    graphics::text(limits[[2L]], y + 0.25, sprintf("%+.3f", rows$delta_primary_minus_comparator),
      adj = 1, cex = 0.88, col = colour)
    graphics::title(titles[[i]], adj = 0, cex.main = 1.1)
  }
  graphics::mtext("Paired AUC difference: seven-gene ratio minus comparator (95% CI)", side = 1, outer = TRUE, line = 0.7, cex = 1)
  graphics::mtext("Pointwise paired DeLong intervals; multiplicity-adjusted P values are in Table S4.",
    side = 1, outer = TRUE, line = 2.0, cex = 0.79)
  graphics::mtext("Expanded paired discrimination comparisons", side = 3, outer = TRUE, line = 0.5, cex = 1.25, font = 2)
}

output_component_distribution <- function(data) {
  cohorts <- c("GSE193677", "GSE117993", "GSE111889")
  roles <- c("Development", "External", "HMP2 context only")
  units <- c("log2 TMM-adjusted CPM", "log2 TPM", "log2 pooled CPM")
  colours <- c(CD = "#63829B", UC = "#CE772F")
  graphics::par(mfrow = c(3, 1))
  output_theme(mar = c(3.1, 5.7, 3.3, 1.5), oma = c(3, 0, 2, 0))
  for (i in seq_along(cohorts)) {
    table <- data$components[data$components$cohort == cohorts[[i]] & data$components$diagnosis %in% c("CD", "UC"), , drop = FALSE]
    limits <- range(table$logG_normalized_no_pseudo, table$logD_normalized_no_pseudo)
    span <- diff(limits)
    graphics::plot(NA, xlim = c(0.45, 2.55), ylim = limits + c(-0.08, 0.20) * span,
      xaxt = "n", xlab = "", ylab = units[[i]])
    graphics::abline(h = graphics::axTicks(2), col = "#EEF2F4")
    for (j in 1:2) {
      field <- c("logG_normalized_no_pseudo", "logD_normalized_no_pseudo")[[j]]
      for (disease in c("CD", "UC")) {
        rows <- table[table$diagnosis == disease, , drop = FALSE]
        position <- j + if (disease == "CD") -0.18 else 0.18
        output_box(rows[[field]], rows$participant_id, position, colours[[disease]], width = 0.30, point_cex = 0.43)
      }
    }
    graphics::axis(1, at = 1:2, labels = c("IgG (G)", "IgA + IgM (D)"), tick = FALSE)
    graphics::legend("topright", legend = c("CD", "UC"), pch = 15, col = colours, bty = "n", horiz = TRUE, cex = 0.9)
    graphics::title(sprintf("%s. %s: %s | CD n = %d; UC n = %d", LETTERS[[i]], roles[[i]], cohorts[[i]],
      sum(table$diagnosis == "CD"), sum(table$diagnosis == "UC")), adj = 0, cex.main = 1.04)
  }
  graphics::mtext("Normalized component distributions", side = 3, outer = TRUE, line = 0.4, cex = 1.3, font = 2)
  graphics::mtext("No pseudocounts; all component sums are positive. Units and axes are cohort-specific.", side = 1, outer = TRUE, line = 0.6, cex = 0.85)
  graphics::mtext("Points show all UC/CD participants. Controls are excluded. These are descriptive abundance distributions.",
    side = 1, outer = TRUE, line = 1.9, cex = 0.78)
}

output_background <- function(data) {
  values <- data$background$auc_uc_positive
  features <- data$feature_context
  breaks <- seq(floor(min(values, 0.25) * 200) / 200, ceiling(max(values, 0.75) * 200) / 200, by = 0.005)
  histogram <- graphics::hist(values, breaks = breaks, plot = FALSE, include.lowest = TRUE)
  if (sum(histogram$counts) != length(values)) stop("The empirical histogram omitted genes.")
  limits <- range(breaks)
  graphics::layout(matrix(1:2, ncol = 1L), heights = c(1.15, 1))
  output_theme(mar = c(2.5, 12.8, 3.6, 1.5), oma = c(3.2, 0, 1, 0))
  graphics::plot(histogram, freq = TRUE, col = "#CBD3D8", border = "#FFFFFF", main = "",
    xlab = "", ylab = "Number of genes", xlim = limits)
  graphics::abline(v = 0.5, lty = 2, col = "#7F929B")
  graphics::title("Empirical transcriptome background in development", line = 2, cex.main = 1.18)
  graphics::mtext(sprintf("Observed labels; %s genes; expression filter applied without disease labels", format(length(values), big.mark = ",")),
    side = 3, line = 0.4, cex = 0.86)
  graphics::par(mar = c(4.5, 12.8, 1.5, 1.5))
  graphics::plot(NA, xlim = limits, ylim = c(0.5, 5.5), xlab = "UC-positive AUC", ylab = "", yaxt = "n", bty = "n")
  graphics::abline(v = 0.5, lty = 2, col = "#7F929B")
  labels <- c("Seven-gene ratio", "IGHG3 / other Ig", "IGHG3 abundance", "Total Ig abundance", "Plasma-cell score")
  graphics::axis(2, at = 5:1, labels = labels, tick = FALSE, cex.axis = 0.98)
  colour <- c("#267C78", rep("#647D8A", 4L))
  graphics::points(features$auc_uc_positive, 5:1, pch = 16, col = colour, cex = 1.2)
  graphics::text(features$auc_uc_positive, 5:1, sprintf("%.3f", features$auc_uc_positive),
    pos = ifelse(features$auc_uc_positive < 0.5, 2, 4), offset = 0.6, col = colour, cex = 0.98)
  graphics::mtext("Descriptive observed-label context; no permutation null or correction for score exploration/selection.",
    side = 1, outer = TRUE, line = 0.7, cex = 0.82)
  primary <- output_one(features, "feature", "IGG_IGAM_LOG_RATIO")
  graphics::mtext(sprintf("Seven-gene ratio: %.1fth percentile by absolute distance from 0.5; %d genes more extreme.",
    primary$transcriptome_percentile, primary$genes_more_extreme), side = 1, outer = TRUE, line = 2, cex = 0.82)
}

output_device <- function(path, format, width, height, draw, data) {
  if (format == "png") {
    if (capabilities("cairo")) grDevices::png(path, width = width, height = height, units = "in", res = 300, type = "cairo", bg = "white") else
      grDevices::png(path, width = width, height = height, units = "in", res = 300, bg = "white")
  } else {
    if (capabilities("cairo")) grDevices::cairo_pdf(path, width = width, height = height, family = "sans", bg = "white") else
      grDevices::pdf(path, width = width, height = height, family = "Helvetica", useDingbats = FALSE, version = "1.4", bg = "white")
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  draw(data)
  invisible(path)
}

output_file_check <- function(path) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size < 100L) stop("Missing or empty output: ", basename(path), ".")
  extension <- tolower(tools::file_ext(path))
  if (extension %in% c("pdf", "png")) {
    connection <- file(path, "rb")
    on.exit(close(connection), add = TRUE)
    header <- readBin(connection, "raw", n = 8L)
    valid <- if (extension == "png") identical(header, as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))) else
      identical(rawToChar(header[1:5]), "%PDF-")
    if (!valid) stop("Invalid figure signature: ", basename(path), ".")
  }
  invisible(TRUE)
}

output_manuscript_tables <- function(data, path) {
  cohorts <- c("GSE193677", "GSE117993", "GSE111889")
  geometry <- data$discovery_geometry
  populations <- c(paste(sprintf("%s: %d / %d", geometry$compartment, geometry$n_uc, geometry$n_cd), collapse = "; "),
    sprintf("%d / %d", data$replication_geometry[["UC"]], data$replication_geometry[["CD"]]),
    sprintf("%d / %d", data$primary$n_uc, data$primary$n_cd), "Not estimable", "Not estimable")
  table1 <- data.frame(Dataset = c("GSE282122 (TAURUS)", "GSE214695", "GSE193677", "GSE117993",
    "GSE111889 / HMP2", "GSE83687", "GSE57945"),
    `Study role` = c("Discovery", "Secondary candidate replication", "Bulk development",
      "Principal external score evaluation", "Context-only evaluation", "Nonestimable score endpoint", "Nonestimable score endpoint"),
    `Material and anatomical site` = c("scRNA-seq; left colon; before study adalimumab initiation", "scRNA-seq; colon",
      "Bulk RNA-seq; rectum", "Bulk RNA-seq; rectum", "Bulk RNA-seq; rectum", "Bulk RNA-seq; surgical intestinal specimens", "Bulk RNA-seq; ileum"),
    `UC / CD, n` = populations,
    `Expression input` = c("Raw counts", "Raw counts", "Estimated counts", "TPM",
      "Counts; pooled CPM for abundance comparators", "FPKM", "RPKM"), check.names = FALSE, stringsAsFactors = FALSE)
  features <- c("primary_ratio", "summed_IgG_abundance", "IGHG3_abundance", "total_Ig_abundance", "IGHG3_fraction")
  labels <- c("Seven-gene ratio", "Summed IgG abundance", "IGHG3 abundance", "Total Ig abundance", "IGHG3-to-other-Ig ratio")
  table2 <- do.call(rbind, lapply(cohorts, function(cohort) {
    aucs <- data$aucs[data$aucs$cohort == cohort, , drop = FALSE]
    aucs <- aucs[match(features, aucs$feature), , drop = FALSE]
    paired <- data$paired[data$paired$cohort == cohort, , drop = FALSE]
    paired <- paired[match(features[-1L], paired$comparator), , drop = FALSE]
    if (anyNA(aucs$feature) || anyNA(paired$comparator)) stop("A manuscript comparator row is missing.")
    p <- paired$p_holm_4_within_cohort
    display_p <- ifelse(p < 0.001, formatC(p, format = "e", digits = 2L), sprintf("%.3f", p))
    data.frame(Cohort = cohort, Measure = labels,
      `AUC (95% CI)` = sprintf("%.3f (%.3f to %.3f)", aucs$auc, aucs$ci95_low, aucs$ci95_high),
      `Delta AUC (95% CI)` = c("Reference", sprintf("%.3f (%.3f to %.3f)",
        paired$delta_primary_minus_comparator, paired$ci95_low, paired$ci95_high)),
      `Holm-adjusted P` = c("Not applicable", display_p), check.names = FALSE, stringsAsFactors = FALSE)
  }))
  escape <- function(value) {
    value <- gsub("&", "&amp;", as.character(value), fixed = TRUE)
    value <- gsub("<", "&lt;", value, fixed = TRUE)
    gsub(">", "&gt;", value, fixed = TRUE)
  }
  html_table <- function(table, widths) {
    if (length(widths) != ncol(table) || sum(widths) != 100) stop("Invalid manuscript table column widths.")
    header <- paste0("<th>", escape(names(table)), "</th>", collapse = "")
    columns <- paste0('<col width="', widths, '%" style="width:', widths, '%">', collapse = "")
    rows <- vapply(seq_len(nrow(table)), function(i) {
      paste0("<tr>", paste0("<td>", escape(unlist(table[i, ], use.names = FALSE)), "</td>", collapse = ""), "</tr>")
    }, character(1))
    c(paste0('<table width="100%" cellpadding="6" cellspacing="0"><colgroup>', columns,
      "</colgroup><thead><tr>", header, "</tr></thead><tbody>"), rows, "</tbody></table>")
  }
  control_n <- vapply(c("GSE117993", "GSE111889"), function(cohort) {
    sum(data$components$cohort == cohort & data$components$diagnosis == "Control")
  }, integer(1))
  lines <- c('<!doctype html><html lang="en"><head><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1"><title>IBD manuscript tables</title>',
    '<style>@page{size:A4 landscape;margin:14mm}body{font:15px/1.5 Arial,sans-serif;color:#17384b;margin:32px auto;max-width:1300px;padding:0 24px}h1{font-size:25px}h2{font-size:20px;margin-top:32px}table{border-collapse:collapse;table-layout:fixed;width:100%;font-size:13px}th,td{padding:10px;text-align:left;vertical-align:top;border-bottom:1px solid #d5dfe4;overflow-wrap:normal;word-break:normal}th{background:#17384b;color:white}tbody tr:nth-child(even){background:#f2f6f8}p{font-size:13px;color:#425863}@media print{body{margin:0;padding:0;font-size:10pt}table{font-size:9pt}th,td{padding:5pt}thead{display:table-header-group}tr{break-inside:avoid;page-break-inside:avoid}h2{break-after:avoid;page-break-after:avoid}}</style></head><body>',
    '<h1>IBD manuscript tables</h1><p>Generated in R from validated computed results. Numbers are formatted for display; full calculated precision remains in the R caches and supplementary workbook.</p>',
    '<h2>Table 1. Datasets, analysis populations, and study roles</h2>', html_table(table1, c(14, 20, 27, 19, 20)),
    sprintf('<p>UC counts precede CD counts. Single-cell counts are donors; bulk counts are participants. TAURUS compartment populations overlap and must not be added. Development: %d biopsies pooled before normalization. HMP2: %d rectal samples pooled. Controls (%d external; %d HMP2) are descriptive only and excluded from UC-CD discrimination.</p>',
      data$development_biopsies, data$context_rectal_samples, control_n[[1L]], control_n[[2L]]),
    '<p>GSE117993 was a pediatric diagnostic cohort. HMP2 included pediatric and adult participants and is context only because of prior project exposure. GSE83687 lacks the required constant-gene annotation; GSE57945 has zero required components. Nonestimability is neither a negative validation nor a claim of zero participants. Material/unit descriptions are declared study metadata.</p>',
    '<h2>Table 2. Discrimination and paired comparisons with the seven-gene ratio</h2>', html_table(table2, c(13, 25, 22, 25, 15)),
    '<p>UC is the positive class. Delta AUC is ratio minus comparator. Intervals are pointwise 95% DeLong intervals. Two-sided paired tests use retrospective Holm correction across four comparisons within each cohort. Development comparisons are exploratory and do not account for score selection; summed IgG is a post hoc addition. HMP2 is context only. Nonsignificance does not establish equivalence. Expanded multiplicity sensitivity is in Table S4.</p>',
    '<p>The score is log2[(G + 0.5)/(A + M + 0.5)], where G = IGHG1 + IGHG2 + IGHG3 + IGHG4, A = IGHA1 + IGHA2, and M = IGHM. Total Ig additionally includes IGHD, IGKC, IGLC1, IGLC2, and IGLC3. IGHG3-to-other-Ig is log2[(IGHG3 + 0.5)/(total Ig - IGHG3 + 0.5)]. Abundance transforms are cohort-specific and remain defined in the analysis modules.</p>',
    '</body></html>')
  writeLines(lines, path, useBytes = TRUE)
  list(table1 = table1, table2 = table2)
}

output_review_zip <- function(directory, relative_files) {
  archive <- file.path(directory, "results", "Publication_Outputs.zip")
  zip::zipr(archive, files = relative_files, root = directory, include_directories = FALSE, mode = "mirror")
  listing <- utils::unzip(archive, list = TRUE)$Name
  if (anyDuplicated(listing) || !setequal(listing, relative_files)) stop("Review ZIP contains unexpected members.")
  check_directory <- tempfile("ibd_output_zip_check_")
  if (!dir.create(check_directory)) stop("Cannot create ZIP verification storage.")
  on.exit(unlink(check_directory, recursive = TRUE), add = TRUE)
  utils::unzip(archive, exdir = check_directory, unzip = "internal")
  original <- unname(tools::md5sum(file.path(directory, relative_files)))
  restored <- unname(tools::md5sum(file.path(check_directory, relative_files)))
  if (anyNA(original) || !identical(original, restored)) stop("Review ZIP byte round-trip failed.")
  archive
}

output_promote <- function(staged_results, destination, staged_cache, cache_path, allowed_files) {
  # Roll back both destinations if a Windows lock or another promotion error occurs.
  backup_results <- tempfile("ibd_previous_results_", tmpdir = dirname(destination))
  backup_cache <- tempfile("ibd_previous_outputs_", tmpdir = dirname(cache_path))
  moved_results <- moved_cache <- installed_results <- installed_cache <- committed <- FALSE
  on.exit({
    if (!committed) {
      if (installed_results) unlink(destination, recursive = TRUE)
      if (installed_cache) unlink(cache_path)
      if (moved_results && !file.rename(backup_results, destination)) warning("Previous results could not be restored: ", backup_results)
      if (moved_cache && !file.rename(backup_cache, cache_path)) warning("Previous output cache could not be restored: ", backup_cache)
    } else {
      if (moved_results) unlink(backup_results, recursive = TRUE)
      if (moved_cache) unlink(backup_cache)
    }
  }, add = TRUE)
  if (file.exists(destination)) {
    if (!dir.exists(destination)) stop("The results destination is not a directory.")
    existing <- gsub("\\", "/", list.files(destination, recursive = TRUE, all.files = TRUE), fixed = TRUE)
    unexpected <- setdiff(existing, allowed_files)
    if (length(unexpected)) {
      # Retire only byte-verified R products from stage07.0, never an original AI asset.
      obsolete <- c("figures/Figure_1_overview.png", "figures/Figure_1_overview.pdf")
      if (!all(unexpected %in% obsolete) || !file.exists(cache_path)) {
        stop("Results contains unrelated files. Move those files before regenerating the declared outputs.")
      }
      previous <- readRDS(cache_path)
      if (!identical(previous$code_revision, "stage07.0") || !isTRUE(previous$validation_passed) ||
          !is.data.frame(previous$inventory) || !all(c("file", "md5") %in% names(previous$inventory))) {
        stop("Cannot establish the provenance of existing Figure 1 files; they were preserved.")
      }
      inventory <- previous$inventory
      indices <- match(paste0("results/", unexpected), inventory$file)
      actual <- unname(tools::md5sum(file.path(destination, unexpected)))
      if (anyDuplicated(inventory$file) || anyNA(indices) || anyNA(actual) ||
          !identical(actual, unname(inventory$md5[indices]))) {
        stop("Existing Figure 1 bytes differ from the previous R output inventory; they were preserved.")
      }
    }
    moved_results <- file.rename(destination, backup_results)
    if (!moved_results) stop("Cannot replace results. Close any open output workbook/figures and retry.")
  }
  if (file.exists(cache_path)) {
    moved_cache <- file.rename(cache_path, backup_cache)
    if (!moved_cache) stop("Cannot replace the previous output cache.")
  }
  installed_results <- file.rename(staged_results, destination)
  if (!installed_results) stop("Cannot promote the verified results directory.")
  installed_cache <- file.rename(staged_cache, cache_path)
  if (!installed_cache) stop("Cannot promote the verified output cache.")
  committed <- TRUE
  invisible(TRUE)
}

run_outputs <- function(context) {
  project <- context$project
  work <- file.path(project, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  report_path <- file.path(work, "outputs_report.txt")
  connection <- file(report_path, "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  counters <- c(PASS = 0L, WARN = 0L, FAIL = 0L)
  emit <- function(lines) {
    writeLines(as.character(lines), connection)
    flush(connection)
    cat(paste(lines, collapse = "\n"), "\n", sep = "")
  }
  record <- function(status, name, detail) {
    counters[[status]] <<- counters[[status]] + 1L
    emit(paste(status, name, detail, sep = " | "))
  }
  status <- 1L
  tryCatch(withCallingHandlers({
    emit(c("IBD publication workbook and figure generation report", "Code revision: stage07.3",
      paste("Started UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      "Inputs: validated newly computed R caches with a consistent dependency chain.",
      "All table assembly, workbook writing, statistical consistency checks, and computational figure rendering are in R.",
      "Reference TSVs are comparison targets only. No model or score is selected in this stage.",
      "Visual review is pending until the generated workbook and figures are inspected."))
    required <- c(openxlsx = "4.2.8.1", zip = "3.0.2", pROC = "1.19.0.1",
      edgeR = "4.4.2", limma = "3.62.2", Matrix = "1.7-1")
    for (package in names(required)) {
      if (!requireNamespace(package, quietly = TRUE) || utils::packageVersion(package) != package_version(required[[package]])) {
        stop("Missing or incompatible package: ", package, ".")
      }
      record("PASS", paste0("package_", package), as.character(utils::packageVersion(package)))
    }
    if (getRversion() != package_version("4.4.2")) record("WARN", "R_version", paste("Observed", as.character(getRversion()), "; validation used R 4.4.2."))
    if (!isTRUE(context$settings$verify_md5)) stop("Output generation requires verified cache provenance.")
    check_validation_logic()
    record("PASS", "validator_self_checks", "Shared exact identifiers/counts, probability and numerical comparison rules.")
    output_reference_flag_self_check()
    record("PASS", "reference_flag_self_check", "Reference True/False spellings map to logical values; invalid/missing flags and changed truth values are rejected.")
    output_jitter_self_check()
    record("PASS", "jitter_self_check", "1000 identifiers retain deterministic offsets within the declared group width; input order reversal preserves identifier-aligned positions.")
    output_workbook_self_check()
    record("PASS", "workbook_self_check", "Six sheets; independent S3 schemas; leading-zero text IDs, booleans, missing values, fractions, and small numbers round-trip.")
    manifest <- context$read_manifest(file.path(project, "config", "inputs.tsv"))
    mapping <- context$read_manifest(file.path(project, "config", "tables.tsv"))
    if (nrow(mapping) != 7L || anyDuplicated(mapping$source_file) ||
        !identical(mapping$worksheet, c("Table_S1", "Table_S2", "Table_S3", "Table_S3", "Table_S4", "Table_S5", "Table_S6")) ||
        !all(vapply(mapping$source_file, context$safe_relative, logical(1)))) stop("Malformed publication table mapping.")
    cache_md5 <- character()
    read_cache <- function(stage, revision, identifiers) {
      path <- file.path(work, paste0(stage, ".rds"))
      if (!file.exists(path)) stop("Missing prerequisite: work/", stage, ".rds.")
      emit(paste("INFO | cache | Reading validated", stage, "results."))
      before <- unname(tools::md5sum(path))
      bundle <- readRDS(path)
      if (is.na(before) || !identical(before, unname(tools::md5sum(path))) ||
          !isTRUE(bundle$validation_passed) || !identical(bundle$code_revision, revision)) {
        stop("Changed, failed, or incompatible cache: ", stage, ".")
      }
      for (identifier in identifiers) {
        item <- manifest[manifest$input_id == identifier, , drop = FALSE]
        if (nrow(item) != 1L || !identical(unname(bundle$source_md5[identifier]), item$expected_md5[[1L]])) {
          stop("Cache source identity disagrees with input definitions: ", stage, " / ", identifier, ".")
        }
      }
      cache_md5[stage] <<- before
      record("PASS", paste0(stage, "_cache"), paste("Validated", revision, "; MD5:", before))
      bundle
    }
    discovery <- read_cache("discovery", "stage02.0", c("discovery_myeloid", "discovery_fibperi"))
    replication <- read_cache("replication", "stage03.0", c("replication_archive", "replication_annotations", "gse214695_sample_metadata"))
    development <- read_cache("development", "stage04.2", c("development_counts", "gse193677_sample_metadata"))
    external <- read_cache("external", "stage05.0", c("external_archive", "gse117993_sample_metadata", "context_counts",
      "context_metadata", "nonestimable_archive", "nonestimable_rpkm"))
    sensitivity <- read_cache("sensitivity", "stage06.0", "discovery_myeloid")
    if (!identical(replication$discovery_cache_md5, unname(cache_md5["discovery"])) ||
        !identical(development$replication_cache_md5, unname(cache_md5["replication"])) ||
        !identical(external$development_cache_md5, unname(cache_md5["development"])) ||
        !identical(unname(sensitivity$cache_md5[c("discovery", "development", "external")]),
          unname(cache_md5[c("discovery", "development", "external")]))) stop("The publication caches belong to inconsistent dependency chains.")
    record("PASS", "cache_dependency_chain", "Discovery -> replication -> development -> external; sensitivity cites the same three source caches.")
    tables <- output_tables(development, replication, external, sensitivity, mapping)
    record("PASS", "supplementary_schemas", "Seven complete computed blocks; original column order and types retained.")
    output_reference_check(tables, mapping, project, record)
    plot_data <- output_plot_data(discovery, replication, development, external, sensitivity, record)
    rm(discovery, replication, development, external, sensitivity)
    invisible(gc())
    record("PASS", "figure_data", "Full participant populations, fixed UC-positive primary/subtype AUCs, 12 paired contrasts, and all 18940 background genes.")
    staging <- tempfile("ibd_outputs_", tmpdir = project)
    if (!dir.create(staging)) stop("Cannot create output staging directory.")
    on.exit(unlink(staging, recursive = TRUE), add = TRUE)
    staged_results <- file.path(staging, "results")
    if (!dir.create(file.path(staged_results, "figures"), recursive = TRUE)) stop("Cannot create staged results directories.")
    workbook <- file.path(staged_results, "Supplementary_Tables.xlsx")
    output_workbook(workbook, tables, mapping)
    difference <- output_roundtrip(workbook, tables, mapping)
    record("PASS", "saved_workbook_package", "Worksheet dimensions match declared content; all internal relationship and content-type targets resolve; unused drawing placeholders removed in R.")
    record("PASS", "saved_workbook_roundtrip", paste("Exactly six sheets and seven blocks; numeric types, values, IDs, flags and missingness checked; maximum absolute difference:",
      format(difference, digits = 8L, scientific = TRUE)))
    manuscript <- output_manuscript_tables(plot_data, file.path(staged_results, "Manuscript_Tables.html"))
    record("PASS", "manuscript_table_assembly", "Two main tables; cohort counts and all comparator statistics formatted from computed results.")
    specifications <- data.frame(stem = c("Figure_2_discrimination", "Figure_3_external_distributions",
      "Figure_4_paired_comparisons", "Figure_5_components", "Figure_S1_transcriptome_context"),
      width = c(12, 10, 10, 10, 10), height = c(6.4, 7, 11, 12, 9), stringsAsFactors = FALSE)
    drawers <- list(output_forest, output_external_distribution,
      output_paired, output_component_distribution, output_background)
    figure_files <- character()
    for (i in seq_len(nrow(specifications))) {
      emit(paste("INFO | figure | Rendering", specifications$stem[[i]], "in R."))
      for (extension in c("png", "pdf")) {
        name <- paste0(specifications$stem[[i]], ".", extension)
        path <- file.path(staged_results, "figures", name)
        output_device(path, extension, specifications$width[[i]], specifications$height[[i]], drawers[[i]], plot_data)
        output_file_check(path)
        figure_files <- c(figure_files, file.path("results", "figures", name))
      }
      record("PASS", specifications$stem[[i]], "PNG at 300 dpi and vector PDF generated; nonempty file signatures checked. Visual inspection remains pending.")
    }
    relative_files <- gsub("\\", "/", c(file.path("results", "Supplementary_Tables.xlsx"),
      file.path("results", "Manuscript_Tables.html"), figure_files), fixed = TRUE)
    for (path in file.path(staging, relative_files)) output_file_check(path)
    archive <- output_review_zip(staging, relative_files)
    record("PASS", "review_zip_roundtrip", "Exactly 12 declared publication files; extracted member bytes agree with the generated originals.")
    paths <- c(relative_files, "results/Publication_Outputs.zip")
    hashes <- unname(tools::md5sum(file.path(staging, paths)))
    if (anyNA(hashes)) stop("Cannot hash publication artifacts.")
    inventory <- data.frame(file = paths, bytes = unname(file.info(file.path(staging, paths))$size), md5 = hashes, stringsAsFactors = FALSE)
    # Require stable prerequisites through the end of rendering, not just at load.
    current_hashes <- unname(tools::md5sum(file.path(work, paste0(names(cache_md5), ".rds"))))
    if (!identical(current_hashes, unname(cache_md5))) stop("A prerequisite cache changed during output generation.")
    if (counters[["WARN"]] > 0L || counters[["FAIL"]] > 0L) stop("Review the reported warnings/failures before promoting publication outputs.")
    staged_cache <- file.path(staging, "outputs.rds")
    bundle <- list(code_revision = "stage07.3", created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      cache_md5 = cache_md5, supplementary = tables, manuscript = manuscript, figure_data = plot_data,
      inventory = inventory, validation_passed = TRUE, visual_review_status = "PENDING",
      excluded_figures = c(Figure_1 = "Existing AI-generated graphical abstract; not generated, copied, or modified by this stage."),
      release_status = "NOT_RUN", validation_counts = counters, session = utils::sessionInfo())
    saveRDS(bundle, staged_cache, compress = "xz")
    output_promote(staged_results, file.path(project, "results"), staged_cache, file.path(work, "outputs.rds"),
      sub("^results/", "", paths))
    record("PASS", "output_promotion", "Verified artifacts and small provenance bundle promoted; all five prerequisite caches retained.")
    emit(c("", "Generated artifact inventory", utils::capture.output(print(inventory, row.names = FALSE)),
      "", paste("Review ZIP:", file.path(project, "results", "Publication_Outputs.zip")),
      paste("Output cache:", file.path(work, "outputs.rds")),
      "Return this report and Publication_Outputs.zip for workbook/figure visual review.",
      "Figure 1 remains the existing AI-generated graphical abstract; this stage does not generate, copy, or modify it.",
      "No source matrix, private configuration, participant-level file or RDS cache is included in the review ZIP.",
      "The manuscript DOCX and supplementary DOCX have not been edited in this stage.",
      "R session information", utils::capture.output(utils::sessionInfo())))
    status <- 0L
  }, warning = function(warning) {
    record("WARN", "R_warning", conditionMessage(warning))
    invokeRestart("muffleWarning")
  }), error = function(error) {
    record("FAIL", "outputs_error", conditionMessage(error))
    status <<- 1L
  })
  emit(c("", paste("Checks: PASS =", counters[["PASS"]], "; WARN =", counters[["WARN"]], "; FAIL =", counters[["FAIL"]]),
    paste("Workbook/figure generation:", if (status == 0L) "PASS" else "BLOCKED"),
    "Visual publication review: PENDING", "Portable clean-run verification and GitHub release: NOT_RUN",
    "A failed run does not certify or refresh older outputs; use this current report.",
    paste("Finished UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Report:", normalizePath(report_path, winslash = "/", mustWork = TRUE))))
  status
}

release_requirements <- function(description) {
  fields <- intersect(c("Depends", "Imports", "LinkingTo"), names(description))
  entries <- trimws(unlist(strsplit(paste(unlist(description[fields]), collapse = ","), ",", fixed = TRUE)))
  entries <- entries[nzchar(entries)]
  result <- data.frame(package = character(), operator = character(), version = character(), stringsAsFactors = FALSE)
  for (entry in entries) {
    parts <- regmatches(entry, regexec("^([A-Za-z][A-Za-z0-9.]*)[[:space:]]*(\\(([><=!]+)[[:space:]]*([^()[:space:]]+)[[:space:]]*\\))?$", entry))[[1L]]
    if (!length(parts) || !parts[[4L]] %in% c("", ">=", "<=", "==", "!=", ">", "<")) stop("Invalid installed dependency declaration: ", entry, ".")
    result[nrow(result) + 1L, ] <- c(parts[[2L]], parts[[4L]], parts[[5L]])
  }
  result
}

release_version_satisfies <- function(actual, operator, required) {
  if (!nzchar(operator)) return(TRUE)
  comparison <- utils::compareVersion(actual, required)
  switch(operator, ">=" = comparison >= 0L, "<=" = comparison <= 0L,
    "==" = comparison == 0L, "!=" = comparison != 0L, ">" = comparison > 0L, "<" = comparison < 0L,
    stop("Unsupported dependency version operator."))
}

release_source_record <- function(description, requirements) {
  record <- description[intersect(c("Package", "Version", "Repository", "biocViews"), names(description))]
  remote_type <- description$RemoteType
  repository <- description$Repository
  if (!is.null(remote_type) && identical(tolower(remote_type), "github")) {
    remote_fields <- c("RemoteType", "RemoteHost", "RemoteUsername", "RemoteRepo", "RemoteRef", "RemoteSha", "RemoteSubdir")
    required <- c("RemoteUsername", "RemoteRepo", "RemoteSha")
    if (!all(required %in% names(description)) || any(!nzchar(unlist(description[required]))) ||
        !grepl("^[0-9a-fA-F]{40}$", description$RemoteSha)) stop("Incomplete commit-pinned GitHub source for ", description$Package, ".")
    record <- c(record, description[intersect(remote_fields, names(description))])
    record$Source <- "GitHub"
  } else if (!is.null(remote_type) && !tolower(remote_type) %in% c("standard", "cran", "bioc", "bioconductor")) {
    stop("Unsupported installed source for ", description$Package, ": RemoteType=", remote_type, ".")
  } else if (!is.null(repository) && identical(repository, "CRAN")) {
    record$Source <- "Repository"
  } else if (!is.null(repository) && grepl("^Bioconductor([[:space:]]+[0-9]+[.][0-9]+)?$", repository)) {
    if (!repository %in% c("Bioconductor", "Bioconductor 3.20")) stop("Installed Bioconductor source is outside release 3.20: ", description$Package, ".")
    record$Source <- "Bioconductor"
  } else {
    stop("Installed metadata do not establish a supported source for ", description$Package,
      "; Repository=", if (is.null(repository)) "<missing>" else repository,
      "; RemoteType=", if (is.null(remote_type)) "<missing>" else remote_type, ". No source was inferred online.")
  }
  record$Requirements <- sort(unique(requirements$package))
  record
}

release_installed_records <- function(packages, libraries, emit = function(lines) invisible(NULL), reader = NULL) {
  # Read only installed DESCRIPTION files. Do not query repository catalogs or
  # infer remote origins from package names or current online versions.
  if (is.null(reader)) reader <- function(package) {
    location <- find.package(package, lib.loc = libraries, quiet = TRUE)
    if (length(location) != 1L || !nzchar(location)) stop("Missing installed dependency: ", package, ".")
    description <- read.dcf(file.path(location, "DESCRIPTION"))
    if (nrow(description) != 1L) stop("Invalid installed DESCRIPTION for ", package, ".")
    as.list(description[1L, , drop = TRUE])
  }
  descriptions <- list()
  get_description <- function(package) {
    if (is.null(descriptions[[package]])) {
      description <- reader(package)
      if (!identical(description$Package, package) || is.null(description$Version)) stop("Invalid installed package identity: ", package, ".")
      descriptions[[package]] <<- description
    }
    descriptions[[package]]
  }
  pending <- unique(packages)
  visited <- character()
  records <- list()
  while (length(pending)) {
    package <- pending[[1L]]
    pending <- pending[-1L]
    if (package %in% visited) next
    visited <- c(visited, package)
    description <- get_description(package)
    if (identical(description$Priority, "base")) next
    requirements <- release_requirements(description)
    for (i in seq_len(nrow(requirements))) {
      dependency <- requirements$package[[i]]
      actual <- if (dependency == "R") as.character(getRversion()) else get_description(dependency)$Version
      if (!release_version_satisfies(actual, requirements$operator[[i]], requirements$version[[i]])) {
        stop("Unsatisfied installed dependency: ", package, " requires ", dependency, " ",
          requirements$operator[[i]], " ", requirements$version[[i]], "; installed ", actual, ".")
      }
    }
    records[[package]] <- release_source_record(description, requirements)
    pending <- unique(c(pending, setdiff(requirements$package, c("R", visited))))
    emit(paste("INFO | installed_dependency |", package, description$Version, "; source:", records[[package]]$Source))
  }
  records[sort(names(records))]
}

release_local_lock_self_check <- function() {
  fixture <- list(
    FixtureA = list(Package = "FixtureA", Version = "1.0", Repository = "CRAN", Depends = "R (>= 4.0)", Imports = "FixtureB (>= 1.2), stats"),
    FixtureB = list(Package = "FixtureB", Version = "1.2", Repository = "Bioconductor 3.20", LinkingTo = "FixtureA (== 1.0)"),
    stats = list(Package = "stats", Version = as.character(getRversion()), Priority = "base"))
  reader <- function(package) {
    if (is.null(fixture[[package]])) stop("Missing fixture dependency.")
    fixture[[package]]
  }
  records <- release_installed_records("FixtureA", character(), reader = reader)
  if (!identical(names(records), c("FixtureA", "FixtureB")) ||
      !identical(records$FixtureB$Source, "Bioconductor")) stop("Local dependency closure self-check failed.")
  rejects <- function(expression) inherits(tryCatch(force(expression), error = function(error) error), "error")
  fixture$FixtureB$Version <- "1.1"
  if (!rejects(release_installed_records("FixtureA", character(), reader = reader))) stop("Unsatisfied dependency was not rejected.")
  fixture$FixtureB <- NULL
  if (!rejects(release_installed_records("FixtureA", character(), reader = reader))) stop("Missing dependency was not rejected.")
  unknown <- list(Package = "FixtureUnknown", Version = "1.0")
  if (!rejects(release_source_record(unknown, release_requirements(unknown)))) stop("Unknown source was not rejected.")
  if (!rejects(release_requirements(list(Imports = "FixtureA (approximately 1.0)")))) stop("Malformed dependency was not rejected.")
  invisible(TRUE)
}

run_release <- function(context) {
  # Administrative preparation only: copy reviewed R products and capture the
  # current dependency metadata. Do not refit analyses, install, or publish.
  project <- context$project
  work <- file.path(project, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  report_path <- file.path(work, "release_report.txt")
  connection <- file(report_path, "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  counters <- c(PASS = 0L, WARN = 0L, FAIL = 0L)
  emit <- function(lines) {
    writeLines(as.character(lines), connection)
    flush(connection)
    cat(paste(lines, collapse = "\n"), "\n", sep = "")
  }
  record <- function(status, name, detail) {
    counters[[status]] <<- counters[[status]] + 1L
    emit(paste(status, name, detail, sep = " | "))
  }
  checksum <- function(path) {
    value <- unname(tools::md5sum(path))
    if (length(value) != 1L || is.na(value)) stop("Cannot hash a release prerequisite.")
    value
  }
  write_table <- function(table, path) utils::write.table(table, path, sep = "\t",
    quote = FALSE, row.names = FALSE, col.names = TRUE, na = "NA", fileEncoding = "UTF-8")
  status <- 1L
  tryCatch(withCallingHandlers({
    emit(c("IBD portable release candidate preparation report", "Code revision: stage08.1",
      paste("Started UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      "Scope: reviewed artifact identity, bundled source metadata, dependency lock, and portable packaging.",
      "No scientific calculation, package installation, environment activation, Git commit, or publication is performed.",
      "Fresh dependency restoration and a complete isolated analysis run remain separate required checks."))
    if (getRversion() != package_version("4.4.2")) stop("Prepare the candidate with the validated R version 4.4.2.")
    if (!isTRUE(context$settings$verify_md5)) stop("Release preparation requires checksum verification.")
    packages <- context$read_manifest(file.path(project, "config", "packages.tsv"))
    if (!identical(names(packages), c("package", "expected_version")) || anyDuplicated(packages$package) || nrow(packages) != 10L) {
      stop("Malformed publication package definitions.")
    }
    for (i in seq_len(nrow(packages))) {
      package <- packages$package[[i]]
      if (!requireNamespace(package, quietly = TRUE) || utils::packageVersion(package) != package_version(packages$expected_version[[i]])) {
        stop("Missing or incompatible package for release preparation: ", package, ".")
      }
      record("PASS", paste0("package_", package), as.character(utils::packageVersion(package)))
    }
    cache_path <- file.path(work, "outputs.rds")
    if (!file.exists(cache_path)) stop("Missing validated work/outputs.rds.")
    output_cache_md5 <- checksum(cache_path)
    output <- readRDS(cache_path)
    if (!identical(output$code_revision, "stage07.3") || !isTRUE(output$validation_passed) ||
        !identical(checksum(cache_path), output_cache_md5)) stop("Changed, failed, or incompatible publication output cache.")
    for (stage in c("discovery", "replication", "development", "external", "sensitivity")) {
      path <- file.path(work, paste0(stage, ".rds"))
      if (!file.exists(path) || !identical(checksum(path), unname(output$cache_md5[stage]))) {
        stop("A prerequisite cache differs from the reviewed output dependency chain: ", stage, ".")
      }
    }
    record("PASS", "reviewed_cache_chain", paste("Validated stage07.3; all five prerequisite caches agree; output cache MD5:", output_cache_md5))
    reviewed <- context$read_manifest(file.path(project, "config", "reviewed_outputs.tsv"))
    if (!identical(names(reviewed), c("file", "bytes", "md5")) || nrow(reviewed) != 12L ||
        anyDuplicated(reviewed$file) || !all(vapply(reviewed$file, context$safe_relative, logical(1))) ||
        !all(startsWith(reviewed$file, "results/")) || any(grepl("Figure_1", reviewed$file, fixed = TRUE))) {
      stop("Malformed reviewed-artifact definitions.")
    }
    inventory <- output$inventory[output$inventory$file != "results/Publication_Outputs.zip", , drop = FALSE]
    if (nrow(inventory) != nrow(reviewed) || anyDuplicated(inventory$file) || !setequal(inventory$file, reviewed$file)) {
      stop("The output cache does not describe exactly the reviewed publication products.")
    }
    inventory <- inventory[match(reviewed$file, inventory$file), , drop = FALSE]
    if (!identical(as.character(inventory$md5), reviewed$md5) ||
        any(inventory$bytes != as.numeric(reviewed$bytes))) stop("Output cache inventory differs from the reviewed artifact set.")
    for (i in seq_len(nrow(reviewed))) {
      path <- file.path(project, reviewed$file[[i]])
      if (!file.exists(path) || file.info(path)$size != as.numeric(reviewed$bytes[[i]]) ||
          !identical(checksum(path), reviewed$md5[[i]])) stop("A reviewed publication artifact changed: ", reviewed$file[[i]], ".")
    }
    record("PASS", "reviewed_publication_artifacts", "All 12 R-generated products match the numerically and visually reviewed stage07.3 return; Figure 1 is excluded.")
    inputs <- context$read_manifest(file.path(project, "config", "inputs.tsv"))
    bundled <- context$read_manifest(file.path(project, "config", "bundled_inputs.tsv"))
    source_identity <- list()
    for (i in seq_len(nrow(inputs))) {
      item <- inputs[i, , drop = FALSE]
      identifier <- item$input_id[[1L]]
      if (!context$safe_relative(item$relative_path[[1L]])) stop("Unsafe declared input path.")
      path <- file.path(context$settings$legacy_root, item$relative_path[[1L]])
      if (identifier %in% names(context$settings$file_overrides)) path <- context$settings$file_overrides[[identifier]]
      if (!file.exists(path)) {
        if (item$required[[1L]] == "TRUE") stop("Missing required native input: ", identifier, ".")
        record("WARN", identifier, "Optional pairing metadata are absent; this file is not an analysis dependency.")
        next
      }
      bytes <- as.numeric(file.info(path)$size)
      md5 <- checksum(path)
      if (is.na(bytes) || bytes != as.numeric(item$expected_bytes[[1L]]) || !identical(md5, item$expected_md5[[1L]])) {
        stop("Native input identity changed: ", identifier, ".")
      }
      source_identity[[length(source_identity) + 1L]] <- data.frame(input_id = identifier, bytes = bytes, md5 = md5, stringsAsFactors = FALSE)
    }
    source_identity <- do.call(rbind, source_identity)
    record("PASS", "native_input_identities", paste("Verified inputs:", nrow(source_identity), "; paths are omitted from the public identity record."))
    if (nrow(bundled) != 4L || !all(bundled$input_id %in% source_identity$input_id)) stop("Incomplete portable source metadata.")
    for (i in seq_len(nrow(bundled))) {
      path <- file.path(project, bundled$bundled_path[[i]])
      if (!identical(checksum(path), bundled$expected_md5[[i]])) stop("Bundled source metadata changed.")
    }
    record("PASS", "portable_metadata", "Four exact source metadata files travel with the code; explicit user overrides retain precedence.")
    staging <- tempfile("ibd_release_", tmpdir = work)
    candidate <- file.path(staging, "IBD_Publication_Code")
    if (!dir.create(candidate, recursive = TRUE)) stop("Cannot create the release candidate staging directory.")
    on.exit(unlink(staging, recursive = TRUE), add = TRUE)
    source_files <- c("run.R", paste0("R/", c("discovery", "replication", "development", "external", "sensitivity", "outputs"), ".R"),
      paste0("config/", c("local.example.R", "inputs.tsv", "tables.tsv", "bundled_inputs.tsv", "packages.tsv", "reviewed_outputs.tsv", "restoration_reference.tsv", "restore_recovery.tsv", "verification_reference.tsv")),
      paste0("reference/", c("tables_reference.zip", "discovery_expected.zip", "replication_expected.zip",
        "development_expected.zip", "external_expected.zip", "sensitivity_expected.zip")),
      bundled$bundled_path, "README.md", ".gitignore", reviewed$file)
    for (relative in source_files) {
      if (!context$safe_relative(relative)) stop("Unsafe candidate member path.")
      from <- file.path(project, relative)
      to <- file.path(candidate, relative)
      dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(from, to, overwrite = FALSE) || !identical(checksum(from), checksum(to))) {
        stop("Candidate byte-copy failed: ", relative, ".")
      }
    }
    record("PASS", "portable_file_copy", "Seven main R files, declared configuration/references, four source metadata files, and reviewed products copied byte for byte; private local.R and caches excluded.")
    repositories <- c(BioCsoft = "https://bioconductor.org/packages/3.20/bioc",
      BioCann = "https://bioconductor.org/packages/3.20/data/annotation",
      BioCexp = "https://bioconductor.org/packages/3.20/data/experiment",
      BioCworkflows = "https://bioconductor.org/packages/3.20/workflows",
      BioCbooks = "https://bioconductor.org/packages/3.20/books", CRAN = "https://cloud.r-project.org")
    previous_libraries <- .libPaths()
    lock_packages <- packages$package
    if (requireNamespace("BiocManager", quietly = TRUE)) lock_packages <- c(lock_packages, "BiocManager")
    release_local_lock_self_check()
    record("PASS", "local_dependency_self_checks", "Recursive/cyclic dependencies and base packages handled; missing/incompatible dependencies, malformed declarations and unknown sources rejected.")
    emit("INFO | renv | Reading recursive installed DESCRIPTION metadata locally; no repository lookup, installation or activation.")
    records <- release_installed_records(lock_packages, previous_libraries, emit)
    installed <- utils::installed.packages(lib.loc = previous_libraries, noCache = TRUE)
    installed <- installed[!duplicated(installed[, "Package"]), , drop = FALSE]
    dependencies <- tools::package_dependencies(packages = lock_packages, db = installed,
      which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE)
    base_packages <- installed[which(installed[, "Priority"] == "base"), "Package"]
    expected <- setdiff(unique(c(lock_packages, unlist(dependencies, use.names = FALSE))), c("R", base_packages))
    if (!setequal(names(records), expected)) stop("Local dependency closure disagrees with the independent installed-package database.")
    record("PASS", "local_dependency_closure", paste("Complete required package set agrees with tools::package_dependencies; packages:", length(records)))
    lock <- list(R = list(Version = as.character(getRversion()), Repositories = as.list(repositories)),
      Bioconductor = list(Version = "3.20"), Packages = records)
    lock_path <- file.path(candidate, "renv.lock")
    emit(paste("INFO | renv | Writing and reopening local dependency lock; package records:", length(records)))
    renv::lockfile_write(lock, file = lock_path, project = candidate)
    if (!identical(.libPaths(), previous_libraries)) stop("Dependency capture changed the active library paths.")
    lock <- renv::lockfile_read(file = lock_path, project = candidate)
    if (!setequal(names(lock$Packages), names(records))) stop("Dependency lock package set changed during serialization.")
    for (package in names(records)) {
      for (field in names(records[[package]])) {
        if (!identical(lock$Packages[[package]][[field]], records[[package]][[field]])) {
          stop("Dependency metadata changed during serialization: ", package, "/", field, ".")
        }
      }
    }
    if (!identical(lock$R$Version, "4.4.2") || !identical(lock$Bioconductor$Version, "3.20") ||
        !all(lock_packages %in% names(lock$Packages))) stop("Incomplete R/Bioconductor dependency lock.")
    software <- do.call(rbind, lapply(names(lock$Packages), function(package) {
      item <- lock$Packages[[package]]
      if (!identical(item$Package, package) || is.null(item$Version) || is.null(item$Source) ||
          !item$Source %in% c("Repository", "Bioconductor", "GitHub")) stop("A dependency lacks a supported remote source: ", package, ".")
      if (utils::packageVersion(package, lib.loc = previous_libraries) != package_version(item$Version)) {
        stop("Dependency lock differs from the installed library: ", package, ".")
      }
      if (package %in% packages$package && package_version(item$Version) !=
          package_version(packages$expected_version[match(package, packages$package)])) stop("A direct dependency version changed.")
      if (item$Source == "GitHub" && (is.null(item$RemoteSha) || !nzchar(item$RemoteSha))) stop("A GitHub package is not commit-pinned.")
      data.frame(package = package, version = item$Version, source = item$Source,
        repository = if (is.null(item$Repository)) "NA" else item$Repository,
        role = if (package %in% packages$package) "direct" else "transitive_or_restore_support", stringsAsFactors = FALSE)
    }))
    validated_packages <- c(output$session$otherPkgs, output$session$loadedOnly)
    for (package in intersect(names(validated_packages), software$package)) {
      if (package_version(validated_packages[[package]]$Version) !=
          package_version(software$version[match(package, software$package)])) stop("A dependency differs from the validated output session: ", package, ".")
    }
    record("PASS", "dependency_lock", paste("Direct packages:", nrow(packages), "; total locked packages:", nrow(software), "; R 4.4.2; Bioconductor 3.20; installed and validated-session versions checked."))
    write_table(software, file.path(candidate, "software_versions.tsv"))
    write_table(source_identity, file.path(candidate, "source_identity.tsv"))
    release_status <- data.frame(item = c("preparation_revision", "scientific_output_revision", "numerical_reproduction",
      "reviewed_artifact_identity", "fresh_dependency_restore", "isolated_complete_run", "GitHub_publication"),
      status = c("stage08.1", "stage07.3", "PASS_AUTHOR_REPORTS", "PASS_REVIEWED_MD5S", "NOT_RUN", "NOT_RUN", "NOT_RUN"), stringsAsFactors = FALSE)
    write_table(release_status, file.path(candidate, "release_status.tsv"))
    members <- list.files(candidate, recursive = TRUE, all.files = TRUE, include.dirs = FALSE)
    members <- gsub("\\", "/", members, fixed = TRUE)
    if (any(grepl("(^|/)(local[.]R|work|library|staging|sandbox|[.]RData|[.]Rhistory)(/|$)|[.]rds$", members))) {
      stop("A private configuration, library or cache leaked into the candidate.")
    }
    private_markers <- unique(gsub("\\", "/", c(context$settings$legacy_root, project,
      Sys.getenv("USERPROFILE", unset = "")), fixed = TRUE))
    private_markers <- private_markers[nchar(private_markers) >= 8L]
    for (relative in members[grepl("[.](R|md|tsv|json)$", members)]) {
      lines <- readLines(file.path(candidate, relative), warn = FALSE, encoding = "UTF-8")
      lines <- gsub("\\", "/", lines, fixed = TRUE)
      if (any(vapply(private_markers, function(marker) any(grepl(marker, lines, fixed = TRUE)), logical(1)))) {
        stop("A machine-specific path leaked into the candidate: ", relative, ".")
      }
    }
    record("PASS", "candidate_privacy_scope", "No private input root, local.R, work cache, installed package library, or transport ZIP is included.")
    file_manifest <- data.frame(file = members, bytes = unname(file.info(file.path(candidate, members))$size),
      md5 = unname(tools::md5sum(file.path(candidate, members))), stringsAsFactors = FALSE)
    write_table(file_manifest, file.path(candidate, "release_file_manifest.tsv"))
    members <- c(members, "release_file_manifest.tsv")
    archive <- file.path(staging, "Release_Candidate.zip")
    zip::zipr(archive, files = paste0("IBD_Publication_Code/", members), root = staging,
      include_directories = FALSE, mode = "mirror")
    verify <- tempfile("ibd_release_zip_")
    if (!dir.create(verify)) stop("Cannot create release ZIP verification storage.")
    on.exit(unlink(verify, recursive = TRUE), add = TRUE)
    listing <- utils::unzip(archive, list = TRUE)$Name
    if (anyDuplicated(listing) || !setequal(listing, paste0("IBD_Publication_Code/", members))) stop("Release ZIP membership differs from the candidate.")
    utils::unzip(archive, exdir = verify, unzip = "internal")
    if (!identical(unname(tools::md5sum(file.path(candidate, members))),
        unname(tools::md5sum(file.path(verify, "IBD_Publication_Code", members))))) stop("Release ZIP bytes changed on extraction.")
    record("PASS", "release_zip_roundtrip", paste("Candidate members:", length(members), "; every extracted member matches its original bytes."))
    retained_archive <- file.path(work, "Release_Candidate.zip")
    if (!file.copy(archive, retained_archive, overwrite = TRUE) || !identical(checksum(archive), checksum(retained_archive))) stop("Cannot retain the release candidate ZIP.")
    emit(c("", "Captured package inventory", utils::capture.output(print(software, row.names = FALSE)),
      paste("Release candidate:", normalizePath(retained_archive, winslash = "/", mustWork = TRUE)),
      paste("Release candidate MD5:", checksum(retained_archive)),
      "Return this report and Release_Candidate.zip for lockfile and portability review.",
      "The dependency lock is a capture of the validated installed environment; restoration has not been tested yet.",
      "The candidate is prepared locally. Nothing was pushed, committed, or published."))
    status <- 0L
  }, warning = function(warning) {
    record("WARN", "R_warning", conditionMessage(warning))
    invokeRestart("muffleWarning")
  }), error = function(error) {
    record("FAIL", "release_error", conditionMessage(error))
    status <<- 1L
  })
  emit(c("", paste("Checks: PASS =", counters[["PASS"]], "; WARN =", counters[["WARN"]], "; FAIL =", counters[["FAIL"]]),
    paste("Portable candidate preparation:", if (status == 0L) "PASS" else "BLOCKED"),
    "Fresh dependency restoration and isolated source-to-result verification: NOT_RUN",
    "GitHub publication: NOT_RUN", "A failed run does not certify an older candidate ZIP.",
    paste("Finished UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Report:", normalizePath(report_path, winslash = "/", mustWork = TRUE))))
  status
}

restore_try_lines <- function(path, reader = readLines) {
  # Transient Windows sharing violations are an unavailable read, not a worker
  # failure. Callers distinguish NULL from a successfully read empty file.
  tryCatch(withCallingHandlers(reader(path, warn = FALSE, encoding = "UTF-8"),
    warning = function(warning) invokeRestart("muffleWarning")), error = function(error) NULL)
}

restore_publish_progress <- function(directory, sequence, lines) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  final <- file.path(directory, sprintf("%06d.ready", sequence))
  temporary <- tempfile("progress_", tmpdir = directory, fileext = ".tmp")
  writeLines(as.character(lines), temporary, useBytes = TRUE)
  if (file.exists(final) || !file.rename(temporary, final)) stop("Cannot publish a closed worker progress event.")
  invisible(final)
}

restore_collect_progress <- function(directory, seen = character()) {
  files <- if (dir.exists(directory)) list.files(directory, pattern = "^[0-9]{6}[.]ready$", full.names = TRUE) else character()
  fresh <- setdiff(files, seen)
  lines <- character()
  for (path in fresh) {
    value <- restore_try_lines(path)
    if (is.null(value)) next
    lines <- c(lines, value)
    seen <- c(seen, path)
  }
  list(lines = lines, seen = seen)
}

restore_snapshot_file <- function(from, to, timeout = 3) {
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("snapshot_", tmpdir = dirname(to), fileext = ".tmp")
  on.exit(unlink(temporary), add = TRUE)
  started <- Sys.time()
  repeat {
    copied <- tryCatch(withCallingHandlers({
      if (!file.exists(from) || !isTRUE(file.copy(from, temporary, overwrite = TRUE))) FALSE else {
        original <- unname(tools::md5sum(from))
        snapshot <- unname(tools::md5sum(temporary))
        if (length(original) != 1L || is.na(original) || !identical(original, snapshot)) FALSE else {
          if (file.exists(to)) unlink(to)
          isTRUE(file.rename(temporary, to)) && identical(original, unname(tools::md5sum(to)))
        }
      }
    }, warning = function(warning) invokeRestart("muffleWarning")), error = function(error) FALSE)
    if (isTRUE(copied)) return(TRUE)
    if (as.numeric(difftime(Sys.time(), started, units = "secs")) >= timeout) return(FALSE)
    Sys.sleep(0.1)
  }
}

restore_io_self_check <- function() {
  unavailable <- function(...) { warning("Simulated sharing violation."); stop("Simulated unavailable read.") }
  if (!is.null(restore_try_lines("fixture", reader = unavailable))) stop("A transient sharing violation was not tolerated.")
  fixture <- tempfile("restore_io_fixture_")
  dir.create(fixture)
  on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
  progress <- file.path(fixture, "progress")
  restore_publish_progress(progress, 1L, c("Fixture first", "Fixture second"))
  writeLines("Uncommitted event", file.path(progress, "unfinished.tmp"))
  observed <- restore_collect_progress(progress)
  if (!identical(observed$lines, c("Fixture first", "Fixture second")) ||
      length(restore_collect_progress(progress, observed$seen)$lines)) stop("Closed progress handoff self-check failed.")
  restore_publish_progress(progress, 2L, "Fixture third")
  if (!identical(restore_collect_progress(progress, observed$seen)$lines, "Fixture third")) stop("Incremental progress self-check failed.")
  original <- file.path(fixture, "original.txt")
  snapshot <- file.path(fixture, "snapshot.txt")
  writeLines(c("Complete fixture", "Last line"), original)
  if (!restore_snapshot_file(original, snapshot) ||
      !identical(unname(tools::md5sum(original)), unname(tools::md5sum(snapshot)))) stop("Diagnostic snapshot self-check failed.")
  if (restore_snapshot_file(file.path(fixture, "missing.txt"), snapshot, timeout = 0)) stop("A missing diagnostic file was accepted.")
  invisible(TRUE)
}

restore_namespace_isolation <- function(packages, library, namespaces = loadedNamespaces(),
    path_for = function(package) getNamespaceInfo(package, "path"), system_library = .Library) {
  rows <- lapply(namespaces, function(package) {
    # The base namespace is part of the R executable, not an ordinary package
    # namespace. getNamespaceInfo(base, "path") is explicitly unsupported.
    if (identical(package, "base")) {
      if (!identical(getNamespace("base"), .BaseNamespaceEnv)) stop("Unexpected base namespace identity.")
      return(data.frame(package = package, role = "R_base_namespace", path = NA_character_, stringsAsFactors = FALSE))
    }
    actual <- normalizePath(path_for(package), winslash = "/", mustWork = TRUE)
    if (package %in% packages) {
      expected <- normalizePath(file.path(library, package), winslash = "/", mustWork = TRUE)
      if (!identical(actual, expected)) stop("A locked namespace was loaded outside the restored library: ", package, ".")
      role <- "locked_private_package"
    } else {
      expected <- normalizePath(file.path(system_library, package), winslash = "/", mustWork = TRUE)
      if (!identical(actual, expected)) stop("An undeclared namespace was loaded outside the R system library: ", package, ".")
      description <- read.dcf(file.path(expected, "DESCRIPTION"))
      if (nrow(description) != 1L || !all(c("Package", "Version", "Priority") %in% colnames(description)) ||
          !identical(unname(description[1L, "Package"]), package) ||
          !identical(unname(description[1L, "Priority"]), "base") ||
          package_version(description[1L, "Version"]) != getRversion()) stop("An undeclared non-base namespace was loaded: ", package, ".")
      role <- "R_base_package"
    }
    data.frame(package = package, role = role, path = actual, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

restore_namespace_self_check <- function() {
  fixture <- tempfile("namespace_fixture_")
  on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
  target <- file.path(fixture, "target")
  system <- file.path(fixture, "system")
  for (path in c(file.path(target, "LockedFixture"), file.path(system, "BaseFixture"),
      file.path(system, "RecommendedFixture"), file.path(fixture, "outside"))) dir.create(path, recursive = TRUE)
  writeLines(c("Package: BaseFixture", paste("Version:", getRversion()), "Priority: base"), file.path(system, "BaseFixture", "DESCRIPTION"))
  writeLines(c("Package: RecommendedFixture", paste("Version:", getRversion()), "Priority: recommended"), file.path(system, "RecommendedFixture", "DESCRIPTION"))
  paths <- c(LockedFixture = file.path(target, "LockedFixture"), BaseFixture = file.path(system, "BaseFixture"),
    RecommendedFixture = file.path(system, "RecommendedFixture"))
  resolver <- function(package) {
    if (identical(package, "base")) stop("The base namespace path must not be queried.")
    unname(paths[[package]])
  }
  observed <- restore_namespace_isolation("LockedFixture", target,
    c("base", "BaseFixture", "LockedFixture"), resolver, system)
  if (!identical(observed$role, c("R_base_namespace", "R_base_package", "locked_private_package")) ||
      !is.na(observed$path[[1L]])) stop("Allowed namespace fixture failed.")
  rejects <- function(expression) inherits(tryCatch(force(expression), error = function(error) error), "error")
  if (!rejects(restore_namespace_isolation("LockedFixture", target, "RecommendedFixture", resolver, system))) stop("An undeclared recommended package was accepted.")
  paths[["LockedFixture"]] <- file.path(fixture, "outside")
  if (!rejects(restore_namespace_isolation("LockedFixture", target, "LockedFixture", resolver, system))) stop("A locked package outside the private target was accepted.")
  paths[["BaseFixture"]] <- file.path(fixture, "outside")
  if (!rejects(restore_namespace_isolation("LockedFixture", target, "BaseFixture", resolver, system))) stop("An external package claiming base priority was accepted.")
  # The metadata graph must be allowed to read R's own base DESCRIPTION files.
  if (length(release_installed_records("stats", c(target, .Library)))) stop("R base metadata were treated as a locked dependency.")
  invisible(TRUE)
}

restore_find_reviewed_attempt <- function(context, current_job, emit = function(lines) invisible(NULL)) {
  project <- context$project
  identity <- context$read_manifest(file.path(project, "config", "restore_recovery.tsv"))
  if (!identical(names(identity), c("file", "bytes", "md5")) ||
      !identical(identity$file, c("restore_stdout.log", "restore_worker_report.txt", "restore_code_manifest.tsv",
        "restore_target_after.tsv", "probe_worker_report.txt", "restore_retained_files.tsv")) ||
      !identical(identity$md5, c("005f48faecc1ae3a7fbb8d1c674f99f7", "77d7194dee71c62c1cd1d6b5804f5dda", "f8d782a0f2b31063e7deae29d37c72df",
        "0535b466f372bd1b466e512bcbfa9012", "7a067fd45e3a7649e52999124916ecee", "88d21718152d5a370247a9a41fda80c2"))) stop("Malformed reviewed recovery identities.")
  versions <- context$read_manifest(file.path(project, "software_versions.tsv"))
  expected_packages <- versions$package
  candidates <- list.dirs(file.path(project, "work"), full.names = TRUE, recursive = FALSE)
  candidates <- candidates[grepl("^fresh_restore_", basename(candidates))]
  candidates <- candidates[normalizePath(candidates, winslash = "/", mustWork = TRUE) != current_job]
  for (candidate in candidates) {
    recovered <- tryCatch(withCallingHandlers({
      paths <- file.path(candidate, identity$file)
      if (!all(file.exists(paths)) || any(file.info(paths)$size != as.numeric(identity$bytes)) ||
          !identical(unname(tools::md5sum(paths)), identity$md5)) stop("Different recovery attempt.")
      request <- readRDS(file.path(candidate, "worker_request.rds"))
      done <- readRDS(file.path(candidate, "restore_done.rds"))
      probe <- readRDS(file.path(candidate, "probe_done.rds"))
      if (!identical(done$status, 0L) || done$counters[["WARN"]] != 0L || done$counters[["FAIL"]] != 0L ||
          !identical(probe$status, 1L) || !identical(request$lock_md5, "beed2a03c5c3cf915cdd31c719934bc5")) stop("Unrecognized recovery state.")
      expected_library <- normalizePath(file.path(candidate, "package_library"), winslash = "/", mustWork = TRUE)
      expected_project <- normalizePath(file.path(candidate, "project"), winslash = "/", mustWork = TRUE)
      if (!identical(request$library, expected_library) || !identical(request$project, expected_project) ||
          !identical(unname(tools::md5sum(file.path(expected_project, "renv.lock"))), request$lock_md5)) stop("Recovery paths or lock changed.")
      manifest <- context$read_manifest(file.path(candidate, "restore_code_manifest.tsv"))
      if (!identical(names(manifest), c("file", "md5")) || anyDuplicated(manifest$file) ||
          !all(vapply(manifest$file, context$safe_relative, logical(1))) ||
          !identical(unname(tools::md5sum(file.path(expected_project, manifest$file))), manifest$md5)) stop("The reviewed private source copy changed.")
      for (name in c("work", "results")) {
        if (length(list.files(file.path(expected_project, name), all.files = TRUE, no.. = TRUE))) stop("The prior attempt contains analysis outputs.")
      }
      retained <- context$read_manifest(file.path(candidate, "restore_retained_files.tsv"))
      if (!identical(names(retained), c("file", "md5")) || anyDuplicated(retained$file) ||
          !all(vapply(retained$file, context$safe_relative, logical(1))) ||
          !identical(unname(tools::md5sum(file.path(expected_library, retained$file))), retained$md5)) stop("A reviewed retained package file changed.")
      installed <- utils::installed.packages(lib.loc = expected_library, noCache = TRUE)
      if (nrow(installed) != 19L || anyDuplicated(installed[, "Package"]) ||
          !setequal(installed[, "Package"], expected_packages)) stop("The reviewed 19-package target is incomplete or changed.")
      for (package in expected_packages) {
        if (utils::packageVersion(package, lib.loc = expected_library) !=
            package_version(versions$version[match(package, versions$package)])) stop("A retained package version changed.")
      }
      list(library = expected_library, identity = identity, packages = sort(expected_packages),
        previous_attempt = basename(candidate))
    }, warning = function(warning) invokeRestart("muffleWarning")), error = function(error) {
      emit(paste("INFO | recovery_candidate |", basename(candidate), ";", conditionMessage(error)))
      NULL
    })
    if (!is.null(recovered)) return(recovered)
  }
  NULL
}

restore_target_inventory <- function(lock, library, path) {
  installed <- utils::installed.packages(lib.loc = library, noCache = TRUE)
  inventory <- data.frame(package = installed[, "Package"], version = installed[, "Version"], stringsAsFactors = FALSE)
  utils::write.table(inventory, path, sep = "\t", quote = FALSE, row.names = FALSE)
  unexpected <- setdiff(inventory$package, names(lock$Packages))
  if (anyDuplicated(inventory$package) || length(unexpected)) stop("Unexpected packages in the private target: ", paste(unexpected, collapse = ", "), ".")
  for (package in inventory$package) {
    if (package_version(inventory$version[match(package, inventory$package)]) !=
        package_version(lock$Packages[[package]]$Version)) stop("A private target package has a different version: ", package, ".")
  }
  inventory
}

restore_complete_target <- function(lock, request, emit) {
  before <- restore_target_inventory(lock, request$library, file.path(request$job, "restore_target_before.tsv"))
  retained_files <- as.character(unlist(lapply(before$package, function(package) {
    paste0(package, "/", list.files(file.path(request$library, package), recursive = TRUE, all.files = TRUE, include.dirs = FALSE))
  }), use.names = FALSE))
  retained <- data.frame(file = retained_files,
    md5 = unname(tools::md5sum(file.path(request$library, retained_files))), stringsAsFactors = FALSE)
  if (anyNA(retained$md5)) stop("Cannot hash the retained private installations.")
  utils::write.table(retained, file.path(request$job, "restore_retained_files.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  pending <- setdiff(names(lock$Packages), before$package)
  emit(paste("INFO | target_inventory | Present:", nrow(before), "; missing:", paste(pending, collapse = ", ")))
  if (isTRUE(request$resumed) && length(pending)) stop("The reviewed target changed after recovery checks; no new installation was requested.")
  while (length(pending)) {
    ready <- pending[vapply(pending, function(package) {
      !any(lock$Packages[[package]]$Requirements %in% pending)
    }, logical(1))]
    if (!length(ready)) stop("Cannot order the missing locked package dependencies.")
    package <- ready[[1L]]
    item <- lock$Packages[[package]]
    if (!identical(item$Source, "Repository") || !identical(item$Repository, "CRAN")) {
      stop("A missing package requires explicit non-CRAN recovery: ", package, ". No source or version was substituted.")
    }
    specification <- paste0(package, "@", item$Version)
    emit(paste("INFO | complete_private_target | Installing exact requested package:", specification,
      "; installed packages will be retained; dependency upgrades are not requested."))
    renv::install(packages = specification, project = request$project, library = request$library,
      repos = c(CRAN = "https://cloud.r-project.org"), dependencies = character(),
      rebuild = TRUE, prompt = FALSE, lock = FALSE)
    if (!file.exists(file.path(request$library, package, "DESCRIPTION")) ||
        utils::packageVersion(package, lib.loc = request$library) != package_version(item$Version)) {
      stop("Explicit installation did not place the locked version in the private target: ", package, ".")
    }
    pending <- setdiff(pending, package)
  }
  after <- restore_target_inventory(lock, request$library, file.path(request$job, "restore_target_after.tsv"))
  if (!identical(unname(tools::md5sum(file.path(request$library, retained$file))), retained$md5)) stop("An existing private installation changed while completing the target.")
  emit(paste("INFO | retained_private_installations |", nrow(before), "existing packages verified unchanged byte for byte; target packages:", nrow(after)))
  invisible(after)
}

restore_worker <- function(request) {
  # This function runs in a separate --vanilla Rscript process. The restore
  # process uses only a copied renv bootstrap; the probe uses the new library.
  report <- file.path(request$job, paste0(request$mode, "_worker_report.txt"))
  connection <- file(report, "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  counters <- c(PASS = 0L, WARN = 0L, FAIL = 0L)
  progress_sequence <- 0L
  emit <- function(lines) {
    writeLines(as.character(lines), connection)
    flush(connection)
    cat(paste(lines, collapse = "\n"), "\n", sep = "")
    utils::flush.console()
    progress_sequence <<- progress_sequence + 1L
    restore_publish_progress(request$progress, progress_sequence, lines)
  }
  record <- function(status, name, detail) {
    counters[[status]] <<- counters[[status]] + 1L
    emit(paste(status, name, detail, sep = " | "))
  }
  checksum <- function(path) unname(tools::md5sum(path))
  status <- 1L
  tryCatch(withCallingHandlers({
    emit(c(paste("IBD isolated dependency worker:", request$mode), "Code revision: stage09.3",
      paste("Started UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE))))
    if (!request$mode %in% c("restore", "probe") || getRversion() != package_version("4.4.2")) stop("Invalid worker mode or R version.")
    if (!identical(checksum(request$lockfile), request$lock_md5)) stop("The accepted dependency lock changed.")
    inherited <- names(Sys.getenv())
    Sys.unsetenv(inherited[startsWith(inherited, "RENV_")])
    Sys.unsetenv(c("R_LIBS", "R_LIBS_USER", "R_LIBS_SITE"))
    Sys.setenv(RENV_PROJECT = request$project, RENV_PATHS_ROOT = file.path(request$job, "renv_state"),
      RENV_PATHS_LIBRARY = request$library, RENV_PATHS_CELLAR = file.path(request$job, "empty_cellar"))
    options(renv.config.cache.enabled = FALSE, renv.config.install.shortcuts = FALSE,
      renv.config.external.libraries = character(), renv.config.user.library = FALSE,
      renv.config.auto.snapshot = FALSE, renv.config.autoloader.enabled = FALSE,
      renv.config.activate.prompt = FALSE, renv.config.ppm.enabled = FALSE,
      renv.config.pak.enabled = FALSE, renv.config.connect.retry = 1L,
      renv.config.connect.timeout = 60L, renv.config.install.jobs = 1L,
      timeout = 180L)
    if (request$mode == "restore") {
      if (!isTRUE(request$resumed) && length(list.files(request$library, all.files = TRUE, no.. = TRUE))) stop("The restore target is not empty.")
      if (dir.exists(file.path(request$job, "renv_state")) &&
          length(list.files(file.path(request$job, "renv_state"), recursive = TRUE, all.files = TRUE, no.. = TRUE))) stop("The private renv state was not empty at first restore.")
      record("PASS", if (isTRUE(request$resumed)) "reviewed_restore_scope" else "empty_restore_scope",
        if (isTRUE(request$resumed)) "All 19 retained packages come from the reviewed completed remote target; new private cache; author library/global cache excluded." else
          "New target library and private cache; reuse of installed packages and the global cache disabled.")
    }
    first_library <- if (request$mode == "restore") request$bootstrap else request$library
    .libPaths(c(first_library, .Library), include.site = FALSE)
    if ("renv" %in% loadedNamespaces()) stop("renv was unexpectedly loaded before worker isolation.")
    if (!requireNamespace("renv", lib.loc = first_library, quietly = TRUE) ||
        utils::packageVersion("renv", lib.loc = first_library) != package_version("1.2.4")) stop("Missing isolated renv 1.2.4.")
    lock <- renv::lockfile_read(file = request$lockfile, project = request$project)
    if (!identical(lock$R$Version, "4.4.2") || !identical(lock$Bioconductor$Version, "3.20") || length(lock$Packages) != 19L) stop("Unexpected dependency lock target.")
    repository_records <- lock$R$Repositories
    if (all(vapply(repository_records, function(item) is.character(item) && length(item) == 1L, logical(1)))) {
      repositories <- unlist(repository_records, use.names = TRUE)
    } else {
      if (!all(vapply(repository_records, function(item) is.list(item) && all(c("Name", "URL") %in% names(item)), logical(1)))) stop("Malformed repository records.")
      repositories <- stats::setNames(vapply(repository_records, function(item) item$URL, character(1)),
        vapply(repository_records, function(item) item$Name, character(1)))
    }
    if (!identical(names(repositories), c("BioCsoft", "BioCann", "BioCexp", "BioCworkflows", "BioCbooks", "CRAN"))) stop("Unexpected repository definitions.")
    options(repos = repositories)
    if (request$mode == "restore") {
      if (!isTRUE(request$resumed)) {
        emit("INFO | restore | Restoring locked versions remotely into the private target. The copied renv is bootstrap support only.")
        arguments <- list(project = request$project, library = request$library,
          lockfile = request$lockfile, packages = names(lock$Packages), rebuild = TRUE, prompt = FALSE)
        if ("retry" %in% names(formals(renv::restore))) arguments$retry <- FALSE
        do.call(renv::restore, arguments)
      } else emit("INFO | restore_resume | Rechecking the reviewed 19-package remote target; no downloads, installation or rebuild are requested.")
      restore_complete_target(lock, request, emit)
    }
    installed <- utils::installed.packages(lib.loc = request$library, noCache = TRUE)
    if (nrow(installed) != 19L || anyDuplicated(installed[, "Package"]) ||
        !setequal(installed[, "Package"], names(lock$Packages))) stop("The restored library does not contain exactly the 19 locked packages.")
    for (package in names(lock$Packages)) {
      if (utils::packageVersion(package, lib.loc = request$library) != package_version(lock$Packages[[package]]$Version)) stop("Restored version differs from the lock: ", package, ".")
    }
    record("PASS", "restored_versions", "All 19 locked versions are physically installed in the separate target library.")
    if (request$mode == "probe") {
      for (package in names(lock$Packages)) {
        loadNamespace(package, lib.loc = request$library)
        actual <- normalizePath(getNamespaceInfo(package, "path"), winslash = "/", mustWork = TRUE)
        expected <- normalizePath(file.path(request$library, package), winslash = "/", mustWork = TRUE)
        if (!identical(actual, expected)) stop("A package was loaded outside the restored library: ", package, ".")
        record("PASS", paste0("isolated_namespace_", package), lock$Packages[[package]]$Version)
      }
      namespace_locations <- restore_namespace_isolation(names(lock$Packages), request$library)
      utils::write.table(namespace_locations, file.path(request$job, "restore_namespace_locations.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
      record("PASS", "namespace_isolation", "Every locked namespace comes from the private target; only R base packages may load from the R system library; the base namespace identity is checked directly.")
      records <- release_installed_records(names(lock$Packages), c(request$library, .Library), emit)
      if (!setequal(names(records), names(lock$Packages))) stop("The restored dependency graph is incomplete.")
      for (package in names(records)) {
        for (field in c("Package", "Version", "Source", "Repository", "Requirements")) {
          if (!identical(records[[package]][[field]], lock$Packages[[package]][[field]])) stop("Restored metadata differ from the accepted lock: ", package, "/", field, ".")
        }
      }
      record("PASS", "restored_dependency_metadata", "Installed versions, source stamps and required dependencies agree with the accepted lock.")
      fixture_path <- tempfile("hdf5_probe_", tmpdir = request$job, fileext = ".h5")
      fixture <- matrix(as.integer(seq_len(6L)), nrow = 2L)
      handle <- hdf5r::H5File$new(fixture_path, mode = "w")
      on.exit(try(handle$close_all(), silent = TRUE), add = TRUE)
      handle[["fixture"]] <- fixture
      observed <- handle[["fixture"]][, ]
      handle$close_all()
      unlink(fixture_path)
      if (!identical(observed, fixture)) stop("HDF5 runtime fixture failed.")
      sparse <- Matrix::sparseMatrix(i = c(1L, 2L), j = c(1L, 2L), x = c(1, 2), dims = c(2L, 2L))
      if (!identical(as.numeric(as.matrix(sparse)), c(1, 0, 0, 2))) stop("Sparse Matrix runtime fixture failed.")
      roc <- pROC::roc(c(0, 0, 1, 1), c(1, 2, 3, 4), levels = c(0, 1), direction = "<", quiet = TRUE)
      if (as.numeric(pROC::auc(roc)) != 1) stop("pROC runtime fixture failed.")
      record("PASS", "compiled_runtime_probes", "HDF5 integer roundtrip, sparse Matrix operation and fixed-direction AUC fixture passed.")
      inventory <- data.frame(package = names(lock$Packages),
        version = vapply(names(lock$Packages), function(package) as.character(utils::packageVersion(package, lib.loc = request$library)), character(1)),
        source = vapply(names(lock$Packages), function(package) records[[package]]$Source, character(1)), stringsAsFactors = FALSE)
      utils::write.table(inventory, file.path(request$job, "restore_installed_packages.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    }
    if (!identical(checksum(request$lockfile), request$lock_md5)) stop("The accepted lock changed during restoration.")
    emit(c("R session information", utils::capture.output(utils::sessionInfo())))
    status <- 0L
  }, warning = function(warning) {
    record("WARN", "R_warning", conditionMessage(warning))
    invokeRestart("muffleWarning")
  }), error = function(error) {
    record("FAIL", "worker_error", conditionMessage(error))
    status <<- 1L
  })
  emit(c(paste("Worker status:", if (status == 0L) "PASS" else "BLOCKED"),
    paste("Finished UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE))))
  list(status = status, counters = counters)
}

run_restore <- function(context) {
  project <- context$project
  work <- file.path(project, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  report_path <- file.path(work, "restore_report.txt")
  connection <- file(report_path, "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  counters <- c(PASS = 0L, WARN = 0L, FAIL = 0L)
  emit <- function(lines) {
    writeLines(as.character(lines), connection)
    flush(connection)
    cat(paste(lines, collapse = "\n"), "\n", sep = "")
  }
  record <- function(status, name, detail) {
    counters[[status]] <<- counters[[status]] + 1L
    emit(paste(status, name, detail, sep = " | "))
  }
  checksum <- function(path) unname(tools::md5sum(path))
  job <- tempfile("fresh_restore_", tmpdir = work)
  if (!dir.create(job, recursive = TRUE)) stop("Cannot create isolated restoration storage.")
  job <- normalizePath(job, winslash = "/", mustWork = TRUE)
  snapshots <- file.path(job, "diagnostic_snapshots")
  dir.create(snapshots)
  parent_libraries <- .libPaths()
  status <- 1L
  active_pid <- NULL
  active_pid_path <- NULL
  stop_child <- function() {
    candidate_pid <- active_pid
    if (is.null(candidate_pid) && !is.null(active_pid_path)) {
      for (attempt in seq_len(30L)) {
        lines <- restore_try_lines(active_pid_path)
        if (length(lines) == 1L && grepl("^[0-9]+$", lines)) { candidate_pid <- as.integer(lines); break }
        Sys.sleep(0.1)
      }
    }
    if (is.null(candidate_pid) || length(candidate_pid) != 1L || is.na(candidate_pid) || candidate_pid <= 0L) return(invisible(NULL))
    # Terminate only the worker launched for this private attempt. On Windows,
    # include its package-installation subprocesses when the user interrupts.
    if (.Platform$OS.type == "windows" && nzchar(Sys.which("taskkill"))) {
      try(system2(Sys.which("taskkill"), c("/PID", as.character(candidate_pid), "/T", "/F"),
        stdout = FALSE, stderr = FALSE), silent = TRUE)
    } else try(tools::pskill(candidate_pid), silent = TRUE)
  }
  on.exit(stop_child(), add = TRUE)
  tryCatch({
    emit(c("IBD fresh dependency restoration report", "Code revision: stage09.3",
      paste("Started UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      "Scope: remote restoration or provenance-checked completion of a private library, then a separate --vanilla namespace/runtime probe.",
      "The author library, validated scientific caches and reviewed results are not restoration targets.",
      "No full scientific analysis or GitHub publication is performed by this stage."))
    restore_io_self_check()
    record("PASS", "worker_io_self_checks", "Unavailable reads tolerated; closed progress events arrive once in order; uncommitted events ignored; diagnostic byte snapshots checked.")
    if (getRversion() != package_version("4.4.2") || !isTRUE(context$settings$verify_md5)) stop("Fresh restoration requires R 4.4.2 and checksum verification.")
    restore_namespace_self_check()
    record("PASS", "namespace_guard_self_checks", "The real base namespace bypasses path lookup; R base metadata are readable; ordinary base packages pass; external locked/undeclared recommended namespaces are rejected.")
    if (!requireNamespace("renv", quietly = TRUE) || utils::packageVersion("renv") != package_version("1.2.4")) stop("Bootstrap renv 1.2.4 is required.")
    if (!requireNamespace("zip", quietly = TRUE)) stop("The validated zip package is required for diagnostic packaging.")
    references <- context$read_manifest(file.path(project, "config", "restoration_reference.tsv"))
    if (!identical(names(references), c("file", "bytes", "md5")) ||
        !identical(references$file, c("renv.lock", "software_versions.tsv", "source_identity.tsv"))) stop("Malformed accepted environment identities.")
    if (!identical(references$md5, c("beed2a03c5c3cf915cdd31c719934bc5", "888f62bb8f7d7f4ac1955caecba6ef00", "c2803e81a937e8fc1e601b268b94e9bb"))) stop("Environment identities differ from the accepted stage08.1 return.")
    for (i in seq_len(nrow(references))) {
      path <- file.path(project, references$file[[i]])
      if (!file.exists(path) || file.info(path)$size != as.numeric(references$bytes[[i]]) ||
          !identical(checksum(path), references$md5[[i]])) stop("Accepted environment file changed: ", references$file[[i]], ".")
    }
    record("PASS", "accepted_environment", "Lock and both inventories are byte-identical to the reviewed stage08.1 author return; 19 locked packages.")
    previous_descriptions <- vapply(c("renv", "zip"), function(package) checksum(file.path(find.package(package), "DESCRIPTION")), character(1))
    isolated <- file.path(job, "project")
    library <- file.path(job, "package_library")
    bootstrap <- file.path(job, "bootstrap")
    for (path in c(isolated, library, bootstrap)) if (!dir.create(path, recursive = TRUE)) stop("Cannot create a private restoration directory.")
    recovery <- restore_find_reviewed_attempt(context, job, emit)
    if (!is.null(recovery)) {
      library <- recovery$library
      utils::write.table(recovery$identity, file.path(job, "restore_recovery_identity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
      record("PASS", "reviewed_complete_restore", paste("Retaining 19 exact-version installations from the reviewed completed remote attempt:",
        recovery$previous_attempt, "; its logs, source manifest, lock, retained file hashes, target versions and empty analysis directories agree."))
    } else {
      if (!isTRUE(getOption("ibd.restore.fresh", FALSE))) stop("The reviewed 19-package target was not found unchanged. No remote installation was started. Preserve work/fresh_restore_* and inspect recovery_candidate messages. Use options(ibd.restore.fresh = TRUE) only to request an entirely new remote installation.")
      emit("INFO | restore_scope | An explicit ibd.restore.fresh request permits a new empty remote target.")
    }
    bundled <- context$read_manifest(file.path(project, "config", "bundled_inputs.tsv"))
    source_files <- c("run.R", paste0("R/", c("discovery", "replication", "development", "external", "sensitivity", "outputs"), ".R"),
      paste0("config/", c("local.example.R", "inputs.tsv", "tables.tsv", "bundled_inputs.tsv", "packages.tsv", "reviewed_outputs.tsv", "restoration_reference.tsv", "restore_recovery.tsv", "verification_reference.tsv")),
      paste0("reference/", c("tables_reference.zip", "discovery_expected.zip", "replication_expected.zip", "development_expected.zip", "external_expected.zip", "sensitivity_expected.zip")),
      bundled$bundled_path, references$file, "README.md", ".gitignore")
    for (relative in source_files) {
      if (!context$safe_relative(relative)) stop("Unsafe isolated source path.")
      from <- file.path(project, relative)
      to <- file.path(isolated, relative)
      dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(from, to) || !identical(checksum(from), checksum(to))) stop("Isolated source byte-copy failed: ", relative, ".")
    }
    for (name in c("work", "results")) dir.create(file.path(isolated, name))
    inputs <- context$read_manifest(file.path(project, "config", "inputs.tsv"))
    overrides <- character()
    for (i in seq_len(nrow(inputs))) {
      identifier <- inputs$input_id[[i]]
      if (identifier %in% bundled$input_id) next
      path <- file.path(context$settings$legacy_root, inputs$relative_path[[i]])
      if (identifier %in% names(context$settings$file_overrides)) path <- context$settings$file_overrides[[identifier]]
      overrides[identifier] <- normalizePath(path, winslash = "/", mustWork = inputs$required[[i]] == "TRUE")
    }
    local_lines <- c("# Private generated input pointers for isolated verification only.",
      paste("legacy_root <-", paste(utils::capture.output(dput(context$settings$legacy_root)), collapse = "\n")),
      "verify_md5 <- TRUE", paste("file_overrides <-", paste(utils::capture.output(dput(overrides)), collapse = "\n")))
    writeLines(local_lines, file.path(isolated, "config", "local.R"), useBytes = TRUE)
    if (length(list.files(file.path(isolated, "work"), all.files = TRUE, no.. = TRUE)) ||
        length(list.files(file.path(isolated, "results"), all.files = TRUE, no.. = TRUE))) stop("The isolated analysis directories are not empty.")
    record("PASS", "empty_analysis_copy", "Portable code and exact metadata copied; work/results are empty; only declared native input paths are shared.")
    renv_path <- find.package("renv")
    if (!file.copy(renv_path, bootstrap, recursive = TRUE) ||
        !identical(checksum(file.path(renv_path, "DESCRIPTION")), checksum(file.path(bootstrap, "renv", "DESCRIPTION")))) stop("Cannot copy renv bootstrap support.")
    code_manifest <- data.frame(file = source_files, md5 = unname(tools::md5sum(file.path(isolated, source_files))), stringsAsFactors = FALSE)
    utils::write.table(code_manifest, file.path(job, "restore_code_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    request <- list(job = normalizePath(job, winslash = "/", mustWork = TRUE),
      project = normalizePath(isolated, winslash = "/", mustWork = TRUE),
      library = normalizePath(library, winslash = "/", mustWork = TRUE),
      bootstrap = normalizePath(bootstrap, winslash = "/", mustWork = TRUE),
      lockfile = normalizePath(file.path(isolated, "renv.lock"), winslash = "/", mustWork = TRUE),
      lock_md5 = references$md5[references$file == "renv.lock"], resumed = !is.null(recovery),
      restoration_mode = if (is.null(recovery)) "EMPTY_REMOTE_TARGET" else "REVIEWED_REMOTE_TARGET_REPROBE")
    request_path <- file.path(job, "worker_request.rds")
    saveRDS(request, request_path)
    launcher <- file.path(job, "dependency_worker.R")
    writeLines(c("arguments <- commandArgs(trailingOnly = TRUE)",
      "request <- readRDS(arguments[[1L]])", "request$mode <- arguments[[2L]]",
      "request$progress <- file.path(request$job, paste0(request$mode, '_progress'))",
      "dir.create(request$progress, showWarnings = FALSE)",
      "pid_temporary <- file.path(request$job, paste0(request$mode, '.pid.tmp'))",
      "writeLines(as.character(Sys.getpid()), pid_temporary)",
      "if (!file.rename(pid_temporary, file.path(request$job, paste0(request$mode, '.pid')))) stop('Cannot publish worker PID.')",
      "result <- tryCatch({", "  module <- new.env(parent = baseenv())",
      "  sys.source(file.path(request$project, 'R', 'outputs.R'), envir = module)",
      "  module$restore_worker(request)",
      "}, error = function(error) {", "  cat('FAIL | launcher_error |', conditionMessage(error), '\\n')",
      "  list(status = 1L, counters = c(PASS = 0L, WARN = 0L, FAIL = 1L))", "})",
      "temporary <- file.path(request$job, paste0(request$mode, '_done.tmp'))", "saveRDS(result, temporary)",
      "if (!file.rename(temporary, file.path(request$job, paste0(request$mode, '_done.rds')))) quit(save = 'no', status = 1L)",
      "quit(save = 'no', status = result$status, runLast = FALSE)"), launcher)
    executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    if (!file.exists(executable)) stop("Cannot find the matching Rscript executable.")
    timeout <- getOption("ibd.restore.timeout", 3600)
    if (!is.numeric(timeout) || length(timeout) != 1L || !is.finite(timeout) || timeout < 60) stop("ibd.restore.timeout must be at least 60 seconds.")
    for (mode in c("restore", "probe")) {
      emit(paste("INFO | worker_start |", mode, "; new --vanilla R process; maximum seconds:", timeout))
      logs <- file.path(job, paste0(mode, c("_stdout.log", "_stderr.log")))
      done <- file.path(job, paste0(mode, "_done.rds"))
      pid_path <- file.path(job, paste0(mode, ".pid"))
      progress_path <- file.path(job, paste0(mode, "_progress"))
      active_pid_path <- pid_path
      launch_status <- system2(executable, args = c("--vanilla", shQuote(launcher), shQuote(request_path), mode),
        stdout = logs[[1L]], stderr = logs[[2L]], wait = FALSE)
      if (!identical(as.integer(launch_status), 0L)) stop("Cannot launch the isolated dependency worker.")
      seen <- character()
      started <- Sys.time()
      heartbeat <- started
      repeat {
        if (is.null(active_pid)) {
          pid_lines <- restore_try_lines(pid_path)
          if (length(pid_lines) == 1L && grepl("^[0-9]+$", pid_lines)) active_pid <- as.integer(pid_lines)
        }
        # Redirected stdout/stderr can be exclusively held by Windows until
        # process exit. Never read or zip those live files.
        progress <- restore_collect_progress(progress_path, seen)
        if (length(progress$lines)) emit(progress$lines)
        seen <- progress$seen
        if (file.exists(done)) break
        elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
        if (is.null(active_pid) && elapsed > 30) stop("The worker did not publish its startup PID within 30 seconds; inspect diagnostic logs.")
        if (elapsed > timeout) { stop_child(); active_pid <- NULL; active_pid_path <- NULL; stop("Worker time limit exceeded; inspect captured network/build logs.") }
        if (as.numeric(difftime(Sys.time(), heartbeat, units = "secs")) >= 15) {
          emit(paste("INFO | worker_running |", mode, "; elapsed seconds:", round(elapsed), "; raw download/build logs will be read after the worker closes them."))
          heartbeat <- Sys.time()
        }
        Sys.sleep(1)
      }
      progress <- restore_collect_progress(progress_path, seen)
      if (length(progress$lines)) emit(progress$lines)
      for (path in logs) {
        if (!restore_snapshot_file(path, file.path(snapshots, basename(path)), timeout = 15)) stop("A completed worker log did not become readable: ", basename(path), ".")
      }
      worker_report <- file.path(job, paste0(mode, "_worker_report.txt"))
      if (file.exists(worker_report) && !restore_snapshot_file(worker_report, file.path(snapshots, basename(worker_report)), timeout = 5)) stop("Cannot snapshot the completed worker report.")
      active_pid <- NULL
      active_pid_path <- NULL
      result <- readRDS(done)
      if (!identical(result$status, 0L)) {
        for (path in file.path(snapshots, basename(logs))) {
          lines <- restore_try_lines(path)
          if (length(lines)) emit(c(paste("Captured diagnostic tail:", basename(path)), utils::tail(lines, 40L)))
        }
        stop("The isolated ", mode, " worker failed; inspect its complete captured log.")
      }
      record("PASS", paste0(mode, "_closed_logs"), "Worker completion received; both output logs copied and checked before diagnostic packaging.")
      if (result$counters[["WARN"]] > 0L) record("WARN", paste0(mode, "_worker_warnings"), paste("Warnings captured:", result$counters[["WARN"]], "; review required."))
      record("PASS", paste0(mode, "_worker"), paste("Fresh process completed; worker checks:", result$counters[["PASS"]], "PASS."))
    }
    if (!identical(.libPaths(), parent_libraries) || !identical(previous_descriptions,
        vapply(c("renv", "zip"), function(package) checksum(file.path(find.package(package), "DESCRIPTION")), character(1)))) stop("The parent library metadata changed.")
    if (!identical(unname(tools::md5sum(file.path(isolated, source_files))), code_manifest$md5)) stop("A copied source/reference file changed during restoration.")
    record("PASS", "parent_and_source_invariance", "Parent library paths/bootstrap metadata and copied code/reference bytes are unchanged.")
    state <- list(code_revision = "stage09.3", created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      validation_passed = counters[["WARN"]] == 0L, request = request,
      source_manifest = code_manifest, analysis_status = "NOT_RUN", validation_counts = counters)
    status <- 0L
  }, interrupt = function(interrupt) {
    stop_child(); active_pid <<- NULL; active_pid_path <<- NULL
    record("FAIL", "restore_interrupted", "Interrupted by the user; only the private restoration attempt was a target.")
  }, error = function(error) {
    stop_child(); active_pid <<- NULL; active_pid_path <<- NULL
    record("FAIL", "restore_error", conditionMessage(error))
  })
  outcome <- if (status != 0L) "BLOCKED" else if (counters[["WARN"]] > 0L) "REVIEW_REQUIRED" else "PASS"
  restoration_mode <- if (exists("request", inherits = FALSE)) request$restoration_mode else "NOT_STARTED"
  summary <- data.frame(item = c("revision", "fresh_restore", "restoration_mode", "accepted_lock_md5", "isolated_complete_analysis", "GitHub_publication"),
    status = c("stage09.3", outcome, restoration_mode, "beed2a03c5c3cf915cdd31c719934bc5", "NOT_RUN", "NOT_RUN"), stringsAsFactors = FALSE)
  utils::write.table(summary, file.path(job, "restore_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  archive <- file.path(work, "Restore_Verification.zip")
  archive_promoted <- FALSE
  tryCatch({
    requested <- list.files(job, pattern = "(_report[.]txt|[.]log|[.]tsv)$", full.names = FALSE)
    sources <- file.path(job, requested)
    if (file.exists(file.path(job, "project", "renv.lock"))) {
      requested <- c(requested, "accepted_renv.lock")
      sources <- c(sources, file.path(job, "project", "renv.lock"))
    }
    captured <- vapply(seq_along(sources), function(i) restore_snapshot_file(sources[[i]], file.path(snapshots, requested[[i]])), logical(1))
    capture <- data.frame(file = requested, status = ifelse(captured, "CAPTURED_BYTE_VERIFIED", "UNAVAILABLE_NOT_ARCHIVED"), stringsAsFactors = FALSE)
    utils::write.table(capture, file.path(snapshots, "diagnostic_capture.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    if (any(!captured)) {
      record(if (status == 0L) "FAIL" else "WARN", "diagnostic_unavailable_files", paste(requested[!captured], collapse = "; "))
      status <- 1L
    }
    review_files <- c(requested[captured], "diagnostic_capture.tsv")
    staged_archive <- file.path(job, "Restore_Verification.complete.zip")
    zip::zipr(staged_archive, files = review_files, root = snapshots, include_directories = FALSE, mode = "mirror")
    listing <- utils::unzip(staged_archive, list = TRUE)$Name
    if (anyDuplicated(listing) || !setequal(listing, review_files)) stop("Diagnostic ZIP member mismatch.")
    verified <- tempfile("restore_review_")
    dir.create(verified)
    on.exit(unlink(verified, recursive = TRUE), add = TRUE)
    utils::unzip(staged_archive, exdir = verified, unzip = "internal")
    if (!identical(unname(tools::md5sum(file.path(snapshots, review_files))), unname(tools::md5sum(file.path(verified, review_files))))) stop("Diagnostic ZIP byte mismatch.")
    archive_md5 <- checksum(staged_archive)
    previous_archive <- tempfile("previous_restore_zip_", tmpdir = work, fileext = ".zip")
    existed <- file.exists(archive)
    if (existed && !file.rename(archive, previous_archive)) stop("Cannot replace the previous diagnostic ZIP; close applications holding it open.")
    if (!file.rename(staged_archive, archive)) {
      if (existed) file.rename(previous_archive, archive)
      stop("Cannot promote the completed diagnostic ZIP.")
    }
    if (!identical(checksum(archive), archive_md5)) {
      unlink(archive)
      if (existed) file.rename(previous_archive, archive)
      stop("Promoted diagnostic ZIP identity changed.")
    }
    if (existed) unlink(previous_archive)
    archive_promoted <- TRUE
    record("PASS", "diagnostic_zip_roundtrip", "Only closed byte-verified snapshots packaged; extracted bytes checked before promoting the completed ZIP. Missing diagnostics, if any, are listed explicitly.")
  }, error = function(error) {
    record("FAIL", "diagnostic_zip_error", conditionMessage(error)); status <<- 1L
  })
  if (status == 0L && counters[["WARN"]] == 0L) {
    state$validation_counts <- counters
    saveRDS(state, file.path(work, "fresh_restore_state.rds"))
  }
  emit(c("", paste("Checks: PASS =", counters[["PASS"]], "; WARN =", counters[["WARN"]], "; FAIL =", counters[["FAIL"]]),
    paste("Fresh dependency restoration:", if (status != 0L) "BLOCKED" else outcome),
    "Full isolated source-to-result analysis: NOT_RUN", "GitHub publication: NOT_RUN",
    paste("Private restoration attempt:", normalizePath(job, winslash = "/", mustWork = TRUE)),
    paste("Review ZIP:", if (archive_promoted) archive else "NOT_CREATED_FOR_THIS_RUN; any previous ZIP is not certified."),
    "Return restore_report.txt and the current Restore_Verification.zip if this run produced it, including after a failed restore.",
    "Keep work/ and results/. A failed attempt does not certify an older fresh_restore_state.rds.",
    paste("Finished UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)), paste("Report:", report_path)))
  if (status == 0L && counters[["WARN"]] > 0L) 1L else status
}

verification_plan <- function() {
  data.frame(stage = c("preflight", "discovery", "replication", "development", "external", "sensitivity", "outputs"),
    revision = c("stage01.3", "stage02.0", "stage03.0", "stage04.2", "stage05.0", "stage06.0", "stage07.3"),
    stringsAsFactors = FALSE)
}

verification_check_report <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  matched <- regmatches(lines, regexec("^Checks: PASS = ([0-9]+) ; WARN = ([0-9]+) ; FAIL = ([0-9]+)$", lines))
  matched <- matched[lengths(matched) == 4L]
  if (length(matched) != 1L) stop("Missing or ambiguous final check counts: ", basename(path), ".")
  counts <- stats::setNames(as.integer(matched[[1L]][-1L]), c("PASS", "WARN", "FAIL"))
  if (anyNA(counts) || counts[["PASS"]] < 1L || counts[["WARN"]] != 0L || counts[["FAIL"]] != 0L) {
    stop("The stage report contains warnings/failures or no completed checks: ", basename(path), ".")
  }
  counts
}

verification_self_check <- function() {
  path <- tempfile("verification_report_fixture_", fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  writeLines("Checks: PASS = 3 ; WARN = 0 ; FAIL = 0", path)
  if (!identical(verification_check_report(path), c(PASS = 3L, WARN = 0L, FAIL = 0L))) stop("Report parser fixture failed.")
  rejects <- function(expression) inherits(tryCatch(force(expression), error = function(error) error), "error")
  for (lines in list("Checks: PASS = 3 ; WARN = 1 ; FAIL = 0", "Checks: PASS = 3 ; WARN = 0 ; FAIL = 1",
      "Checks: PASS = 0 ; WARN = 0 ; FAIL = 0", "Incomplete report",
      rep("Checks: PASS = 3 ; WARN = 0 ; FAIL = 0", 2L))) {
    writeLines(lines, path)
    if (!rejects(verification_check_report(path))) stop("An incomplete or unsuccessful report fixture was accepted.")
  }
  restore_io_self_check()
  restore_namespace_self_check()
  invisible(TRUE)
}

verification_source_files <- function(context) {
  bundled <- context$read_manifest(file.path(context$project, "config", "bundled_inputs.tsv"))
  c("run.R", paste0("R/", c("discovery", "replication", "development", "external", "sensitivity", "outputs"), ".R"),
    paste0("config/", c("local.example.R", "inputs.tsv", "tables.tsv", "bundled_inputs.tsv", "packages.tsv", "reviewed_outputs.tsv",
      "restoration_reference.tsv", "restore_recovery.tsv", "verification_reference.tsv")),
    paste0("reference/", c("tables_reference.zip", "discovery_expected.zip", "replication_expected.zip", "development_expected.zip",
      "external_expected.zip", "sensitivity_expected.zip")), bundled$bundled_path,
    "renv.lock", "software_versions.tsv", "source_identity.tsv", "README.md", ".gitignore")
}

verification_accepted_restore <- function(context) {
  project <- context$project
  state <- readRDS(file.path(project, "work", "fresh_restore_state.rds"))
  if (!identical(state$code_revision, "stage09.3") || !isTRUE(state$validation_passed) ||
      !identical(state$analysis_status, "NOT_RUN") || state$validation_counts[["PASS"]] != 11L ||
      state$validation_counts[["WARN"]] != 0L || state$validation_counts[["FAIL"]] != 0L) stop("A successful reviewed stage09.3 restoration state is required.")
  work <- normalizePath(file.path(project, "work"), winslash = "/", mustWork = TRUE)
  job <- normalizePath(state$request$job, winslash = "/", mustWork = TRUE)
  library <- normalizePath(state$request$library, winslash = "/", mustWork = TRUE)
  isolated <- normalizePath(file.path(job, "project"), winslash = "/", mustWork = TRUE)
  if (!identical(dirname(job), work) || !grepl("^fresh_restore_", basename(job)) ||
      !identical(dirname(dirname(library)), work) || !grepl("^fresh_restore_", basename(dirname(library))) ||
      !identical(basename(library), "package_library") || !identical(state$request$project, isolated) ||
      !identical(state$request$restoration_mode, "REVIEWED_REMOTE_TARGET_REPROBE")) stop("The reviewed private restoration paths changed.")
  references <- context$read_manifest(file.path(project, "config", "verification_reference.tsv"))
  if (!identical(names(references), c("file", "bytes", "md5")) ||
      !identical(references$file, c("restore_summary.tsv", "restore_code_manifest.tsv", "probe_worker_report.txt",
        "restore_namespace_locations.tsv", "restore_retained_files.tsv", "restore_worker_report.txt")) ||
      !identical(references$md5, c("4602433ff4eb12ef8c786280060d33c2", "dd80696b1577b26108dbe1f3ddfba3e8",
        "297ddb26e678f7c4842b2895967cf8f7", "e7b563ec0ef0087eb823eefd3b1b1899", "99493a2de0db1a3e4389b5b73122d49c", "2862652bd5d5e4625ea318fc9ab032eb"))) stop("Malformed reviewed restoration references.")
  paths <- file.path(job, references$file)
  if (!all(file.exists(paths)) || any(file.info(paths)$size != as.numeric(references$bytes)) ||
      !identical(unname(tools::md5sum(paths)), references$md5)) stop("The returned successful restoration diagnostics changed or are missing.")
  for (mode in c("restore", "probe")) {
    result <- readRDS(file.path(job, paste0(mode, "_done.rds")))
    if (!identical(result$status, 0L) || result$counters[["WARN"]] != 0L || result$counters[["FAIL"]] != 0L) stop("The accepted worker completion record is not successful: ", mode, ".")
  }
  manifest <- context$read_manifest(file.path(job, "restore_code_manifest.tsv"))
  if (!identical(manifest, state$source_manifest) || !all(vapply(manifest$file, context$safe_relative, logical(1))) ||
      !identical(unname(tools::md5sum(file.path(isolated, manifest$file))), manifest$md5)) stop("The accepted restoration source copy changed.")
  if (!identical(state$request$lock_md5, "beed2a03c5c3cf915cdd31c719934bc5") ||
      !identical(unname(tools::md5sum(c(file.path(project, "renv.lock"), state$request$lockfile))),
        rep(state$request$lock_md5, 2L))) stop("The accepted R-generated dependency lock changed.")
  for (name in c("work", "results")) {
    if (!dir.exists(file.path(isolated, name)) || length(list.files(file.path(isolated, name), all.files = TRUE, no.. = TRUE))) stop("The restoration project is no longer an unused analysis copy.")
  }
  retained <- context$read_manifest(file.path(job, "restore_retained_files.tsv"))
  if (!identical(names(retained), c("file", "md5")) || anyDuplicated(retained$file) ||
      !all(vapply(retained$file, context$safe_relative, logical(1))) ||
      !identical(unname(tools::md5sum(file.path(library, retained$file))), retained$md5)) stop("The reviewed restored package files changed.")
  versions <- context$read_manifest(file.path(project, "software_versions.tsv"))
  installed <- utils::installed.packages(lib.loc = library, noCache = TRUE)
  if (nrow(installed) != 19L || !setequal(installed[, "Package"], versions$package)) stop("The restored target no longer contains exactly 19 locked packages.")
  for (package in versions$package) {
    if (utils::packageVersion(package, lib.loc = library) != package_version(versions$version[match(package, versions$package)])) stop("A restored package version changed: ", package, ".")
  }
  list(state = state, library = library, retained = retained, references = references)
}

verification_library_guard <- function(request) {
  if (!requireNamespace("renv", lib.loc = request$library, quietly = TRUE)) stop("Cannot load the restored renv namespace.")
  lock <- renv::lockfile_read(file = request$lockfile, project = request$project)
  if (!identical(lock$R$Version, "4.4.2") || !identical(lock$Bioconductor$Version, "3.20") ||
      length(lock$Packages) != 19L || !identical(unname(tools::md5sum(request$lockfile)), request$lock_md5)) stop("The worker dependency lock differs from the reviewed environment.")
  installed <- utils::installed.packages(lib.loc = request$library, noCache = TRUE)
  if (nrow(installed) != 19L || !setequal(installed[, "Package"], names(lock$Packages))) stop("Unexpected private package inventory.")
  for (package in names(lock$Packages)) {
    if (utils::packageVersion(package, lib.loc = request$library) != package_version(lock$Packages[[package]]$Version)) stop("Worker package version differs: ", package, ".")
    loadNamespace(package, lib.loc = request$library)
  }
  locations <- restore_namespace_isolation(names(lock$Packages), request$library)
  records <- release_installed_records(names(lock$Packages), c(request$library, .Library))
  if (!setequal(names(records), names(lock$Packages))) stop("The worker dependency graph differs from the lock.")
  for (package in names(lock$Packages)) {
    for (field in c("Package", "Version", "Source", "Repository", "Requirements")) {
      if (!identical(records[[package]][[field]], lock$Packages[[package]][[field]])) stop("Worker package metadata differ: ", package, "/", field, ".")
    }
  }
  utils::write.table(locations, file.path(request$job, paste0(request$stage, "_namespaces.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  inventory <- data.frame(package = names(lock$Packages),
    version = vapply(names(lock$Packages), function(package) lock$Packages[[package]]$Version, character(1)),
    source = vapply(names(lock$Packages), function(package) records[[package]]$Source, character(1)), stringsAsFactors = FALSE)
  utils::write.table(inventory, file.path(request$job, paste0(request$stage, "_packages.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  lock
}

verification_compare_outputs <- function(request, record) {
  read_tsv <- function(path) utils::read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
    colClasses = "character", check.names = FALSE, stringsAsFactors = FALSE)
  reviewed <- read_tsv(file.path(request$project, "config", "reviewed_outputs.tsv"))
  if (nrow(reviewed) != 12L || !identical(names(reviewed), c("file", "bytes", "md5"))) stop("Unexpected reviewed product definitions.")
  old_paths <- file.path(request$comparison_root, reviewed$file)
  new_paths <- file.path(request$project, reviewed$file)
  if (!all(file.exists(old_paths)) || !identical(unname(tools::md5sum(old_paths)), reviewed$md5) ||
      any(file.info(old_paths)$size != as.numeric(reviewed$bytes))) stop("The reviewed publication products changed before comparison.")
  bundle <- readRDS(file.path(request$project, "work", "outputs.rds"))
  if (!isTRUE(bundle$validation_passed) || !identical(bundle$code_revision, "stage07.3")) stop("The independently computed publication cache did not validate.")
  cache_paths <- file.path(request$project, "work", paste0(names(bundle$cache_md5), ".rds"))
  if (!identical(unname(tools::md5sum(cache_paths)), unname(bundle$cache_md5))) stop("The new publication dependency chain changed.")
  mapping <- read_tsv(file.path(request$project, "config", "tables.tsv"))
  new_difference <- output_roundtrip(file.path(request$project, "results", "Supplementary_Tables.xlsx"), bundle$supplementary, mapping)
  old_difference <- output_roundtrip(file.path(request$comparison_root, "results", "Supplementary_Tables.xlsx"), bundle$supplementary, mapping)
  record("PASS", "independent_workbook_equivalence", paste("Both workbooks agree with the newly computed full-precision values, types, flags, IDs and missingness; six sheets/seven blocks; maximum differences:", format(new_difference, digits = 8L), "/", format(old_difference, digits = 8L)))
  if (!identical(readLines(file.path(request$project, "results", "Manuscript_Tables.html"), warn = FALSE),
      readLines(file.path(request$comparison_root, "results", "Manuscript_Tables.html"), warn = FALSE))) stop("The computed main-table document differs from the reviewed document.")
  record("PASS", "independent_main_tables", "Both main tables and their display values match the reviewed document after line-ending normalization.")
  for (path in new_paths) output_file_check(path)
  publication_zip <- file.path(request$project, "results", "Publication_Outputs.zip")
  members <- utils::unzip(publication_zip, list = TRUE)$Name
  if (anyDuplicated(members) || !setequal(members, reviewed$file)) stop("The new publication ZIP has unexpected members.")
  extracted <- tempfile("verification_publication_")
  dir.create(extracted)
  on.exit(unlink(extracted, recursive = TRUE), add = TRUE)
  utils::unzip(publication_zip, exdir = extracted, unzip = "internal")
  if (!identical(unname(tools::md5sum(file.path(extracted, reviewed$file))), unname(tools::md5sum(new_paths)))) stop("The new publication ZIP differs from its computed products.")
  inventory <- data.frame(file = reviewed$file, bytes = unname(file.info(new_paths)$size),
    md5 = unname(tools::md5sum(new_paths)), reviewed_md5 = reviewed$md5, stringsAsFactors = FALSE)
  inventory$byte_identical <- inventory$md5 == inventory$reviewed_md5
  utils::write.table(inventory, file.path(request$job, "verification_products.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  record("PASS", "independent_publication_package", "Exactly 12 new R-generated products; signatures and ZIP member bytes verified; Figure 1 excluded.")
  png <- endsWith(inventory$file, ".png")
  if (sum(png) != 5L) stop("Unexpected number of computational PNG figures.")
  record("PASS", "figure_identity_audit", paste("PNG byte matches to previously reviewed products:", sum(inventory$byte_identical[png]), "/ 5. Any changed PNG requires visual review; workbook/PDF/ZIP timestamps are not numerical failures."))
  inputs <- read_tsv(file.path(request$project, "config", "inputs.tsv"))
  for (i in seq_len(nrow(inputs))) {
    path <- request$native_paths[[inputs$input_id[[i]]]]
    if (!file.exists(path) || file.info(path)$size != as.numeric(inputs$expected_bytes[[i]]) ||
        !identical(unname(tools::md5sum(path)), inputs$expected_md5[[i]])) stop("A declared native input changed: ", inputs$input_id[[i]], ".")
  }
  record("PASS", "native_input_invariance", "All 14 native input identities agree after the isolated run.")
  invisible(TRUE)
}

verification_worker <- function(request) {
  connection <- file(file.path(request$job, paste0(request$stage, "_worker_report.txt")), "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  counters <- c(PASS = 0L, WARN = 0L, FAIL = 0L)
  sequence <- 0L
  emit <- function(lines) {
    writeLines(as.character(lines), connection); flush(connection)
    cat(paste(lines, collapse = "\n"), "\n", sep = ""); utils::flush.console()
    sequence <<- sequence + 1L
    restore_publish_progress(request$progress, sequence, lines)
  }
  record <- function(status, name, detail) {
    counters[[status]] <<- counters[[status]] + 1L
    emit(paste(status, name, detail, sep = " | "))
  }
  status <- 1L
  stage_counts <- c(PASS = 0L, WARN = 0L, FAIL = 0L)
  tryCatch(withCallingHandlers({
    emit(c(paste("IBD isolated analysis worker:", request$stage), "Code revision: stage10.0",
      paste("Started UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE))))
    plan <- verification_plan()
    if (!request$stage %in% c(plan$stage, "compare") || getRversion() != package_version("4.4.2")) stop("Invalid analysis worker stage or R version.")
    inherited <- names(Sys.getenv())
    Sys.unsetenv(c(inherited[startsWith(inherited, "RENV_")], "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE", "IBD_DATA_ROOT"))
    options(renv.config.auto.snapshot = FALSE, renv.config.autoloader.enabled = FALSE, renv.config.cache.enabled = FALSE,
      renv.config.activate.prompt = FALSE, ibd.worker.return.status = TRUE)
    .libPaths(c(request$library, .Library), include.site = FALSE)
    libraries <- .libPaths()
    setwd(request$project)
    verify_source <- function() {
      if (!identical(unname(tools::md5sum(file.path(request$project, request$source_manifest$file))), request$source_manifest$md5) ||
          !identical(unname(tools::md5sum(file.path(request$project, "config", "local.R"))), request$local_config_md5)) stop("The isolated code or private input pointers changed.")
    }
    verify_source()
    lock <- verification_library_guard(request)
    record("PASS", "isolated_worker_environment", "All 19 locked package versions, source/dependency records and namespace paths verified in this fresh process.")
    if (request$stage == "compare") {
      verification_compare_outputs(request, record)
    } else {
      index <- match(request$stage, plan$stage)
      prior_stages <- if (index > 2L) plan$stage[seq.int(2L, index - 1L)] else character()
      expected_before <- if (length(prior_stages)) paste0(prior_stages, ".rds") else character()
      prior_paths <- if (length(expected_before)) file.path(request$project, "work", expected_before) else character()
      if (!setequal(list.files(file.path(request$project, "work"), pattern = "[.]rds$"), expected_before)) stop("Unexpected caches before the isolated stage.")
      before_md5 <- if (length(prior_paths)) unname(tools::md5sum(prior_paths)) else character()
      if (anyNA(before_md5)) stop("An independently calculated prerequisite cache is missing.")
      if (length(list.files(file.path(request$project, "results"), all.files = TRUE, no.. = TRUE))) stop("Publication products appeared before output generation.")
      emit(paste("INFO | analysis | Running unchanged scientific stage:", request$stage,
        "; raw scientific progress will be captured when its output handles close."))
      result <- source(file.path(request$project, "run.R"), local = new.env(parent = baseenv()), echo = FALSE)$value
      if (!is.numeric(result) || length(result) != 1L || is.na(result) || as.integer(result) != 0L) stop("The scientific stage returned a nonzero or missing status.")
      report <- file.path(request$project, "work", paste0(request$stage, "_report.txt"))
      stage_counts <- verification_check_report(report)
      if (!any(readLines(report, warn = FALSE) == paste("Code revision:", plan$revision[[index]]))) stop("Unexpected scientific module revision in its report.")
      if (length(prior_paths) && !identical(unname(tools::md5sum(prior_paths)), before_md5)) stop("An earlier independent cache was modified.")
      expected_after <- c(expected_before, if (index > 1L) paste0(request$stage, ".rds") else character())
      if (!setequal(list.files(file.path(request$project, "work"), pattern = "[.]rds$"), expected_after)) stop("Unexpected or missing cache after the isolated stage.")
      record("PASS", "scientific_stage_completed", paste(request$stage, "; checks:", stage_counts[["PASS"]], "PASS, 0 WARN, 0 FAIL; earlier independent caches unchanged."))
    }
    verify_source()
    restore_namespace_isolation(names(lock$Packages), request$library)
    if (!identical(.libPaths(), libraries)) stop("Library paths changed during the isolated stage.")
    record("PASS", "worker_invariance", "Source bytes, private input pointers and isolated package paths remain unchanged.")
    emit(c("R session information", utils::capture.output(utils::sessionInfo())))
    status <- if (counters[["WARN"]] == 0L) 0L else 1L
  }, warning = function(warning) {
    record("WARN", "R_warning", conditionMessage(warning)); invokeRestart("muffleWarning")
  }), error = function(error) { record("FAIL", "worker_error", conditionMessage(error)); status <<- 1L })
  emit(c(paste("Worker status:", if (status == 0L) "PASS" else "BLOCKED"),
    paste("Finished UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE))))
  list(status = status, counters = counters, stage_counts = stage_counts)
}

verification_child <- function(request, stage, executable, launcher, timeout, emit) {
  job <- request$job
  snapshots <- file.path(job, "diagnostic_snapshots")
  pid_path <- file.path(job, paste0(stage, ".pid"))
  done <- file.path(job, paste0(stage, "_done.rds"))
  progress_path <- file.path(job, paste0(stage, "_progress"))
  logs <- file.path(job, paste0(stage, c("_stdout.log", "_stderr.log")))
  active_pid <- NULL
  completed <- FALSE
  stop_child <- function() {
    candidate <- active_pid
    if (is.null(candidate)) {
      for (attempt in seq_len(30L)) {
        lines <- restore_try_lines(pid_path)
        if (length(lines) == 1L && grepl("^[0-9]+$", lines)) { candidate <- as.integer(lines); break }
        Sys.sleep(0.1)
      }
    }
    if (length(candidate) != 1L || is.na(candidate) || candidate <= 0L) return(invisible(NULL))
    if (.Platform$OS.type == "windows" && nzchar(Sys.which("taskkill"))) {
      try(system2(Sys.which("taskkill"), c("/PID", as.character(candidate), "/T", "/F"), stdout = FALSE, stderr = FALSE), silent = TRUE)
    } else try(tools::pskill(candidate), silent = TRUE)
  }
  on.exit(if (!completed) stop_child(), add = TRUE)
  emit(paste("INFO | worker_start |", stage, "; fresh --vanilla R process; maximum seconds:", timeout))
  launch_status <- system2(executable,
    args = c("--vanilla", "--default-packages=methods,datasets,utils,grDevices,graphics,stats", shQuote(launcher), paste0("--stage=", stage)),
    stdout = logs[[1L]], stderr = logs[[2L]], wait = FALSE)
  if (!identical(as.integer(launch_status), 0L)) stop("Cannot launch the isolated analysis worker.")
  seen <- character()
  started <- Sys.time()
  heartbeat <- started
  repeat {
    if (is.null(active_pid)) {
      lines <- restore_try_lines(pid_path)
      if (length(lines) == 1L && grepl("^[0-9]+$", lines)) active_pid <- as.integer(lines)
    }
    progress <- restore_collect_progress(progress_path, seen)
    if (length(progress$lines)) emit(progress$lines)
    seen <- progress$seen
    if (file.exists(done)) break
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (is.null(active_pid) && elapsed > 30) stop("The analysis worker did not publish its PID within 30 seconds.")
    if (elapsed > timeout) stop("Analysis worker time limit exceeded: ", stage, ".")
    if (as.numeric(difftime(Sys.time(), heartbeat, units = "secs")) >= 15) {
      emit(paste("INFO | worker_running |", stage, "; elapsed seconds:", round(elapsed),
        "; scientific stdout/stderr are captured after the worker closes them."))
      heartbeat <- Sys.time()
    }
    Sys.sleep(1)
  }
  completed <- TRUE
  progress <- restore_collect_progress(progress_path, seen)
  if (length(progress$lines)) emit(progress$lines)
  for (path in logs) {
    if (!restore_snapshot_file(path, file.path(snapshots, basename(path)), timeout = 15)) stop("Cannot capture a completed analysis log: ", basename(path), ".")
  }
  report <- file.path(job, paste0(stage, "_worker_report.txt"))
  if (file.exists(report) && !restore_snapshot_file(report, file.path(snapshots, basename(report)), timeout = 5)) stop("Cannot capture the completed worker report.")
  result <- readRDS(done)
  if (!identical(result$status, 0L) || result$counters[["WARN"]] != 0L || result$counters[["FAIL"]] != 0L) {
    for (path in file.path(snapshots, basename(logs))) {
      lines <- restore_try_lines(path)
      if (length(lines)) emit(c(paste("Captured diagnostic tail:", basename(path)), utils::tail(lines, 40L)))
    }
    stop("The isolated analysis worker failed: ", stage, ".")
  }
  result$elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  result
}

verification_archive <- function(job, destination) {
  snapshots <- file.path(job, "diagnostic_snapshots")
  names <- list.files(job, pattern = "(_report[.]txt|[.]log|[.]tsv)$", full.names = FALSE)
  paths <- file.path(job, names)
  if (dir.exists(file.path(job, "project", "work"))) {
    reports <- list.files(file.path(job, "project", "work"), pattern = "_report[.]txt$", full.names = FALSE)
    names <- c(names, paste0("analysis_", reports))
    paths <- c(paths, file.path(job, "project", "work", reports))
  }
  publication <- file.path(job, "project", "results", "Publication_Outputs.zip")
  if (file.exists(publication)) { names <- c(names, "Publication_Outputs.zip"); paths <- c(paths, publication) }
  captured <- vapply(seq_along(paths), function(i) restore_snapshot_file(paths[[i]], file.path(snapshots, names[[i]]), timeout = 5), logical(1))
  capture <- data.frame(file = names, status = ifelse(captured, "CAPTURED_BYTE_VERIFIED", "UNAVAILABLE_NOT_ARCHIVED"), stringsAsFactors = FALSE)
  utils::write.table(capture, file.path(snapshots, "diagnostic_capture.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  members <- c(names[captured], "diagnostic_capture.tsv")
  staged <- file.path(job, "Analysis_Verification.complete.zip")
  zip::zipr(staged, files = members, root = snapshots, include_directories = FALSE, mode = "mirror")
  listed <- utils::unzip(staged, list = TRUE)$Name
  if (anyDuplicated(listed) || !setequal(listed, members)) stop("Analysis diagnostic ZIP member mismatch.")
  extracted <- tempfile("analysis_diagnostic_check_")
  dir.create(extracted)
  on.exit(unlink(extracted, recursive = TRUE), add = TRUE)
  utils::unzip(staged, exdir = extracted, unzip = "internal")
  if (!identical(unname(tools::md5sum(file.path(snapshots, members))), unname(tools::md5sum(file.path(extracted, members))))) stop("Analysis diagnostic ZIP byte mismatch.")
  archive_md5 <- unname(tools::md5sum(staged))
  backup <- tempfile("previous_analysis_review_", tmpdir = dirname(destination), fileext = ".zip")
  existed <- file.exists(destination)
  if (existed && !file.rename(destination, backup)) stop("Cannot replace the previous analysis ZIP; close applications holding it open.")
  if (!file.rename(staged, destination)) {
    if (existed) file.rename(backup, destination)
    stop("Cannot promote the completed analysis diagnostic ZIP.")
  }
  if (!identical(unname(tools::md5sum(destination)), archive_md5)) {
    unlink(destination)
    if (existed) file.rename(backup, destination)
    stop("Promoted analysis diagnostic ZIP identity changed.")
  }
  if (existed) unlink(backup)
  list(md5 = archive_md5, unavailable = names[!captured], members = members)
}

run_verify <- function(context) {
  project <- context$project
  work <- file.path(project, "work")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  report_path <- file.path(work, "verification_report.txt")
  connection <- file(report_path, "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  counters <- c(PASS = 0L, WARN = 0L, FAIL = 0L)
  emit <- function(lines) {
    writeLines(as.character(lines), connection); flush(connection)
    cat(paste(lines, collapse = "\n"), "\n", sep = ""); utils::flush.console()
  }
  record <- function(status, name, detail) {
    counters[[status]] <<- counters[[status]] + 1L
    emit(paste(status, name, detail, sep = " | "))
  }
  checksum <- function(paths) unname(tools::md5sum(paths))
  job <- tempfile("clean_analysis_", tmpdir = work)
  if (!dir.create(file.path(job, "diagnostic_snapshots"), recursive = TRUE)) stop("Cannot create clean analysis storage.")
  job <- normalizePath(job, winslash = "/", mustWork = TRUE)
  plan <- verification_plan()
  stages <- data.frame(stage = c(plan$stage, "compare"), status = "NOT_RUN", pass = NA_integer_, warn = NA_integer_,
    fail = NA_integer_, elapsed_seconds = NA_real_, stringsAsFactors = FALSE)
  write_stages <- function() utils::write.table(stages, file.path(job, "verification_stages.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write_stages()
  status <- 1L
  current_stage <- NA_integer_
  archive_promoted <- FALSE
  visual_review <- "NOT_RUN"
  tryCatch(withCallingHandlers({
    emit(c("IBD complete isolated source-to-result verification report", "Code revision: stage10.0",
      paste("Started UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      "Scope: seven unchanged analysis stages in a new empty project, each in a separate --vanilla R process using the reviewed restored library.",
      "Native expression files are shared read-only. No author scientific cache or publication product enters the new analysis directories.",
      "Reviewed publication products are comparison targets only after the full calculation. Figure 1 is excluded.",
      "No package installation, dependency activation, manuscript editing, Git commit or GitHub publication is performed."))
    if (getRversion() != package_version("4.4.2") || !isTRUE(context$settings$verify_md5)) stop("Isolated verification requires R 4.4.2 and input checksums.")
    if (!requireNamespace("zip", quietly = TRUE)) stop("The validated zip package is required for diagnostic packaging.")
    verification_self_check()
    record("PASS", "verification_self_checks", "Incomplete, duplicated, warning/failure and zero-check reports rejected; successful report, namespace and Windows-safe I/O fixtures passed.")
    accepted <- verification_accepted_restore(context)
    record("PASS", "reviewed_restoration", "The successful stage09.3 return, its source copy, 19 versions and all 1540 restored package-file identities agree. No installation requested.")
    old_manifest <- accepted$state$source_manifest
    unchanged <- !old_manifest$file %in% c("run.R", "R/outputs.R", "README.md")
    if (!identical(checksum(file.path(project, old_manifest$file[unchanged])), old_manifest$md5[unchanged])) stop("A scientific module, comparison reference or environment definition changed since restoration.")
    prefix <- function(path) {
      lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
      boundary <- which(lines == "release_requirements <- function(description) {")
      if (length(boundary) != 1L) stop("Cannot identify the frozen scientific output module boundary.")
      lines[seq_len(boundary - 1L)]
    }
    if (!identical(prefix(file.path(project, "R", "outputs.R")),
        prefix(file.path(accepted$state$request$project, "R", "outputs.R")))) stop("The reviewed scientific output-generation code changed.")
    record("PASS", "scientific_code_identity", "The five calculation modules, stage07.3 scientific output code, reference archives, metadata and accepted environment definitions are unchanged.")
    reviewed <- context$read_manifest(file.path(project, "config", "reviewed_outputs.tsv"))
    reviewed_paths <- file.path(project, reviewed$file)
    if (nrow(reviewed) != 12L || !all(file.exists(reviewed_paths)) ||
        !identical(checksum(reviewed_paths), reviewed$md5) || any(file.info(reviewed_paths)$size != as.numeric(reviewed$bytes))) stop("The reviewed author publication products changed or are missing.")
    author_caches <- file.path(work, paste0(c("discovery", "replication", "development", "external", "sensitivity", "outputs"), ".rds"))
    author_caches <- author_caches[file.exists(author_caches)]
    emit("INFO | invariance | Hashing existing author caches without opening fitted objects; these files will not be copied or read by the analysis workers.")
    author_cache_md5 <- checksum(author_caches)
    parent_libraries <- .libPaths()
    source_files <- verification_source_files(context)
    if (anyDuplicated(source_files) || !all(vapply(source_files, context$safe_relative, logical(1)))) stop("Malformed clean-run source whitelist.")
    source_md5 <- checksum(file.path(project, source_files))
    if (anyNA(source_md5)) stop("A required portable source file is missing.")
    isolated <- file.path(job, "project")
    dir.create(isolated)
    for (relative in source_files) {
      to <- file.path(isolated, relative)
      dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(file.path(project, relative), to) || !identical(checksum(to), checksum(file.path(project, relative)))) stop("Clean-run source copy failed: ", relative, ".")
    }
    for (name in c("work", "results")) dir.create(file.path(isolated, name))
    bundled <- context$read_manifest(file.path(project, "config", "bundled_inputs.tsv"))
    inputs <- context$read_manifest(file.path(project, "config", "inputs.tsv"))
    overrides <- character()
    native_paths <- stats::setNames(rep("", nrow(inputs)), inputs$input_id)
    for (i in seq_len(nrow(inputs))) {
      identifier <- inputs$input_id[[i]]
      if (identifier %in% bundled$input_id) {
        path <- file.path(isolated, bundled$bundled_path[match(identifier, bundled$input_id)])
      } else {
        path <- file.path(context$settings$legacy_root, inputs$relative_path[[i]])
        if (identifier %in% names(context$settings$file_overrides)) path <- context$settings$file_overrides[[identifier]]
        path <- normalizePath(path, winslash = "/", mustWork = inputs$required[[i]] == "TRUE")
        overrides[identifier] <- path
      }
      native_paths[identifier] <- normalizePath(path, winslash = "/", mustWork = inputs$required[[i]] == "TRUE")
    }
    writeLines(c("# Private native input pointers for the clean analysis only.",
      paste("legacy_root <-", paste(utils::capture.output(dput(context$settings$legacy_root)), collapse = "\n")),
      "verify_md5 <- TRUE", paste("file_overrides <-", paste(utils::capture.output(dput(overrides)), collapse = "\n"))), file.path(isolated, "config", "local.R"), useBytes = TRUE)
    for (name in c("work", "results")) {
      if (length(list.files(file.path(isolated, name), all.files = TRUE, no.. = TRUE))) stop("The new analysis directories are not empty.")
    }
    manifest <- data.frame(file = source_files, md5 = source_md5, stringsAsFactors = FALSE)
    utils::write.table(manifest, file.path(job, "verification_code_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(accepted$references, file.path(job, "accepted_restoration.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    request <- list(job = job, project = normalizePath(isolated, winslash = "/", mustWork = TRUE),
      library = accepted$library, source_manifest = manifest,
      local_config_md5 = checksum(file.path(isolated, "config", "local.R")),
      lockfile = normalizePath(file.path(isolated, "renv.lock"), winslash = "/", mustWork = TRUE),
      lock_md5 = "beed2a03c5c3cf915cdd31c719934bc5", comparison_root = project, native_paths = native_paths)
    request_path <- file.path(job, "analysis_request.rds")
    saveRDS(request, request_path)
    record("PASS", "empty_source_to_result_copy", "New work/results directories contain no files; only portable code, four exact metadata snapshots, comparison references and private native-input pointers were copied.")
    launcher <- file.path(job, "analysis_worker.R")
    request_literal <- paste(utils::capture.output(dput(request_path)), collapse = "\n")
    writeLines(c(paste0("request <- readRDS(", request_literal, ")"),
      "arguments <- commandArgs(trailingOnly = TRUE)",
      "if (length(arguments) != 1L || !startsWith(arguments, '--stage=')) stop('Invalid worker arguments.')",
      "request$stage <- sub('^--stage=', '', arguments[[1L]])",
      "request$progress <- file.path(request$job, paste0(request$stage, '_progress'))",
      "dir.create(request$progress, showWarnings = FALSE)",
      "pid_temporary <- file.path(request$job, paste0(request$stage, '.pid.tmp'))",
      "writeLines(as.character(Sys.getpid()), pid_temporary)",
      "if (!file.rename(pid_temporary, file.path(request$job, paste0(request$stage, '.pid')))) stop('Cannot publish worker PID.')",
      "result <- tryCatch({", "  module <- new.env(parent = baseenv())",
      "  sys.source(file.path(request$project, 'R', 'outputs.R'), envir = module)",
      "  module$verification_worker(request)", "}, error = function(error) {",
      "  cat('FAIL | launcher_error |', conditionMessage(error), '\\n')",
      "  list(status = 1L, counters = c(PASS = 0L, WARN = 0L, FAIL = 1L))", "})",
      "temporary <- file.path(request$job, paste0(request$stage, '_done.tmp'))", "saveRDS(result, temporary)",
      "if (!file.rename(temporary, file.path(request$job, paste0(request$stage, '_done.rds')))) quit(save = 'no', status = 1L)",
      "quit(save = 'no', status = result$status, runLast = FALSE)"), launcher)
    executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
    if (!file.exists(executable)) stop("Cannot find the matching Rscript executable.")
    timeout <- getOption("ibd.verify.timeout", 7200)
    if (!is.numeric(timeout) || length(timeout) != 1L || !is.finite(timeout) || timeout < 60) stop("ibd.verify.timeout must be at least 60 seconds.")
    for (i in seq_len(nrow(stages))) {
      current_stage <- i
      stages$status[[i]] <- "RUNNING"
      write_stages()
      result <- verification_child(request, stages$stage[[i]], executable, launcher, timeout, emit)
      stages$status[[i]] <- "PASS"
      counts <- if (stages$stage[[i]] == "compare") result$counters else result$stage_counts
      stages$pass[[i]] <- counts[["PASS"]]; stages$warn[[i]] <- counts[["WARN"]]; stages$fail[[i]] <- counts[["FAIL"]]
      stages$elapsed_seconds[[i]] <- result$elapsed_seconds
      write_stages()
      record("PASS", paste0("isolated_", stages$stage[[i]]), paste("Fresh process completed;", counts[["PASS"]], "checks passed; elapsed seconds:", round(result$elapsed_seconds)))
    }
    if (!identical(.libPaths(), parent_libraries) || !identical(checksum(author_caches), author_cache_md5) ||
        !identical(checksum(reviewed_paths), reviewed$md5) || !identical(checksum(file.path(project, source_files)), source_md5)) stop("The author environment, caches, reviewed products or source files changed during the clean run.")
    if (!identical(checksum(file.path(accepted$library, accepted$retained$file)), accepted$retained$md5)) stop("The restored package files changed during the clean run.")
    record("PASS", "author_and_package_invariance", "Author library paths, existing cache bytes, reviewed products, original code and all restored package files remain unchanged.")
    products <- context$read_manifest(file.path(job, "verification_products.tsv"))
    png <- endsWith(products$file, ".png")
    visual_review <- if (sum(png) == 5L && all(products$byte_identical[png] == "TRUE")) "PNG_IDENTITIES_MATCH_REVIEWED" else "REVIEW_REQUIRED"
    new_caches <- file.path(isolated, "work", paste0(plan$stage[-1L], ".rds"))
    cache_inventory <- data.frame(stage = plan$stage[-1L], revision = plan$revision[-1L], md5 = checksum(new_caches), stringsAsFactors = FALSE)
    utils::write.table(cache_inventory, file.path(job, "verification_cache_inventory.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    status <- 0L
  }, warning = function(warning) {
    record("WARN", "R_warning", conditionMessage(warning)); invokeRestart("muffleWarning")
  }), interrupt = function(interrupt) {
    record("FAIL", "verification_interrupted", "Interrupted; only the private analysis worker was stopped. Author caches/results and restored packages are retained.")
  }, error = function(error) { record("FAIL", "verification_error", conditionMessage(error)) })
  if (!is.na(current_stage) && stages$status[[current_stage]] == "RUNNING") stages$status[[current_stage]] <- "BLOCKED"
  write_stages()
  if (counters[["WARN"]] > 0L || counters[["FAIL"]] > 0L) status <- 1L
  summary <- data.frame(item = c("revision", "complete_isolated_analysis", "visual_review", "GitHub_publication"),
    status = c("stage10.0", if (status == 0L) "PASS" else "BLOCKED", visual_review, "NOT_RUN"), stringsAsFactors = FALSE)
  utils::write.table(summary, file.path(job, "verification_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  archive <- file.path(work, "Analysis_Verification.zip")
  tryCatch({
    archived <- verification_archive(job, archive)
    archive_promoted <- TRUE
    if (length(archived$unavailable)) stop("Some diagnostics could not be captured; unavailable members are listed in diagnostic_capture.tsv.")
    record("PASS", "analysis_zip_roundtrip", "Closed diagnostic snapshots and the new R-generated publication ZIP extracted and matched byte for byte; no source matrices, fitted RDS objects or installed packages included.")
  }, error = function(error) { record("FAIL", "verification_zip_error", conditionMessage(error)); status <<- 1L })
  if (status == 0L && counters[["WARN"]] == 0L) {
    state <- list(code_revision = "stage10.0", created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      validation_passed = TRUE, request = request, stages = stages, source_manifest = manifest,
      cache_inventory = cache_inventory, visual_review = visual_review, archive_md5 = archived$md5,
      validation_counts = counters, publication_status = "NOT_RUN")
    saveRDS(state, file.path(work, "clean_run_state.rds"))
  }
  emit(c("", "Isolated stage inventory", utils::capture.output(print(stages, row.names = FALSE)),
    "", paste("Checks: PASS =", counters[["PASS"]], "; WARN =", counters[["WARN"]], "; FAIL =", counters[["FAIL"]]),
    paste("Complete isolated source-to-result verification:", if (status == 0L) "PASS" else "BLOCKED"),
    paste("Computational PNG review:", visual_review), "GitHub publication: NOT_RUN",
    paste("Private analysis project:", file.path(job, "project")),
    paste("Review ZIP:", if (archive_promoted) archive else "NOT_CREATED_FOR_THIS_RUN; any previous ZIP is not certified."),
    "Return verification_report.txt and the current Analysis_Verification.zip, including after failure.",
    "Keep work/ and results/, including the private restoration and clean-analysis folders. Figure 1 and manuscript DOCX files are unchanged.",
    "A failed run does not certify an older clean_run_state.rds or an older verification ZIP.",
    paste("Finished UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)), paste("Report:", report_path)))
  status
}
