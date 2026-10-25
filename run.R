# Entry point for the IBD publication workflow.
# Stages: preflight, discovery, replication, development, external, sensitivity, outputs, release, restore, verify.

local({
  locate_project <- function() {
    source_files <- vapply(sys.frames(), function(frame) {
      value <- frame$ofile
      if (is.character(value) && length(value) == 1L) value else NA_character_
    }, character(1))
    source_files <- source_files[!is.na(source_files)]
    if (length(source_files)) {
      return(dirname(normalizePath(utils::tail(source_files, 1L),
                                   winslash = "/", mustWork = TRUE)))
    }
    arguments <- commandArgs(trailingOnly = FALSE)
    script <- sub("^--file=", "", arguments[startsWith(arguments, "--file=")])
    if (length(script) == 1L) {
      return(dirname(normalizePath(script, winslash = "/", mustWork = TRUE)))
    }
    stop("Cannot locate run.R. Run it using Rscript or source(file.choose()).")
  }

  read_manifest <- function(path) {
    utils::read.delim(path, sep = "\t", header = TRUE, quote = "",
                      comment.char = "", colClasses = "character",
                      check.names = FALSE, na.strings = "NA",
                      stringsAsFactors = FALSE)
  }

  safe_relative <- function(path) {
    normalized <- gsub("\\", "/", path, fixed = TRUE)
    pieces <- strsplit(normalized, "/", fixed = TRUE)[[1L]]
    nzchar(normalized) && !grepl("^(/|[A-Za-z]:)", normalized) &&
      !any(pieces %in% c("..", ""))
  }

  checksum <- function(path) {
    value <- unname(tools::md5sum(path))
    if (length(value) != 1L || is.na(value)) {
      stop("Cannot compute the file checksum: ", path)
    }
    tolower(value)
  }

  absolute_path <- function(path) {
    normalized <- gsub("\\", "/", path.expand(path), fixed = TRUE)
    grepl("^(/|[A-Za-z]:/)", normalized)
  }

  load_settings <- function(project, arguments) {
    settings <- new.env(parent = baseenv())
    settings$legacy_root <- ""
    settings$verify_md5 <- TRUE
    settings$file_overrides <- character()
    local_file <- file.path(project, "config", "local.R")
    if (file.exists(local_file)) sys.source(local_file, envir = settings)
    environment_root <- Sys.getenv("IBD_DATA_ROOT", unset = "")
    if (nzchar(environment_root)) settings$legacy_root <- environment_root
    root_option <- arguments[startsWith(arguments, "--source-root=")]
    if (length(root_option) > 1L) stop("Specify --source-root only once.")
    if (length(root_option)) {
      settings$legacy_root <- sub("^--source-root=", "", root_option)
    }
    if ("--skip-md5" %in% arguments) settings$verify_md5 <- FALSE
    if (!is.character(settings$legacy_root) ||
        length(settings$legacy_root) != 1L ||
        is.na(settings$legacy_root) || !nzchar(settings$legacy_root)) {
      stop("Set legacy_root in config/local.R or supply IBD_DATA_ROOT.")
    }
    if (!is.logical(settings$verify_md5) ||
        length(settings$verify_md5) != 1L || is.na(settings$verify_md5)) {
      stop("verify_md5 must be TRUE or FALSE.")
    }
    if (!absolute_path(settings$legacy_root)) {
      stop("The source directory must be an absolute path.")
    }
    if (!is.character(settings$file_overrides) ||
        (length(settings$file_overrides) &&
         (is.null(names(settings$file_overrides)) ||
          anyNA(names(settings$file_overrides)) ||
          any(!nzchar(names(settings$file_overrides))) ||
          anyDuplicated(names(settings$file_overrides)) ||
          anyNA(settings$file_overrides) ||
          any(!nzchar(settings$file_overrides))))) {
      stop("file_overrides must be a named character vector of input paths.")
    }
    if (length(settings$file_overrides) &&
        !all(absolute_path(settings$file_overrides))) {
      stop("Every file override must be an absolute path.")
    }
    # Exact archived source metadata travel with the code. Explicit user
    # overrides retain precedence; all analysis stages still verify source MD5s.
    bundled <- read_manifest(file.path(project, "config", "bundled_inputs.tsv"))
    inputs <- read_manifest(file.path(project, "config", "inputs.tsv"))
    if (!identical(names(bundled), c("input_id", "bundled_path", "expected_bytes", "expected_md5", "provenance")) ||
        anyDuplicated(bundled$input_id) || !all(bundled$input_id %in% inputs$input_id) ||
        !all(vapply(bundled$bundled_path, safe_relative, logical(1)))) stop("Malformed bundled-input definitions.")
    for (i in seq_len(nrow(bundled))) {
      item <- bundled[i, , drop = FALSE]
      identifier <- item$input_id[[1L]]
      source <- inputs[match(identifier, inputs$input_id), , drop = FALSE]
      path <- file.path(project, item$bundled_path[[1L]])
      if (!identical(item$expected_md5, source$expected_md5) || !identical(item$expected_bytes, source$expected_bytes) ||
          !file.exists(path) || file.info(path)$size != as.numeric(item$expected_bytes) ||
          !identical(checksum(path), item$expected_md5[[1L]])) stop("Bundled source metadata failed identity checks: ", identifier, ".")
      if (!identifier %in% names(settings$file_overrides)) {
        settings$file_overrides[identifier] <- normalizePath(path, winslash = "/", mustWork = TRUE)
      }
    }
    settings
  }

  preflight <- function(project, arguments) {
    work <- file.path(project, "work")
    dir.create(work, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(work)) stop("Cannot create the local work directory.")
    report <- file.path(work, "preflight_report.txt")
    log <- file(report, open = "wt", encoding = "UTF-8")
    on.exit(close(log), add = TRUE)
    counts <- c(PASS = 0L, WARN = 0L, FAIL = 0L)
    emit <- function(lines = "") {
      lines <- as.character(lines)
      writeLines(lines, log)
      flush(log)
      cat(paste(lines, collapse = "\n"), "\n", sep = "")
    }
    record <- function(status, identifier, detail) {
      counts[[status]] <<- counts[[status]] + 1L
      emit(paste(status, identifier, detail, sep = " | "))
    }
    final_status <- 1L
    tryCatch({
      emit(c("IBD publication workflow input readiness report",
             "Code revision: stage01.3",
             paste("Started UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
             paste("Code directory:", project),
             "Scope: input integrity, reference-table structure, and software inventory.",
             "No scientific result, table, figure, or workbook is computed here.", ""))
      settings <- load_settings(project, arguments)
      root <- normalizePath(path.expand(settings$legacy_root),
                            winslash = "/", mustWork = FALSE)
      emit(c(paste("Existing input directory:", root),
             paste("MD5 verification enabled:", settings$verify_md5), ""))
      root_exists <- dir.exists(root)
      record(if (root_exists) "PASS" else "FAIL", "input_directory", root)
      if (!settings$verify_md5) {
        record("WARN", "integrity_scope", "Checksums are disabled; file identity is unverified.")
      }

      input_manifest <- read_manifest(file.path(project, "config", "inputs.tsv"))
      input_columns <- c("input_id", "study", "role", "relative_path", "required",
                         "expected_bytes", "expected_md5", "provenance")
      if (!identical(names(input_manifest), input_columns) ||
          anyNA(input_manifest$input_id) || anyDuplicated(input_manifest$input_id) ||
          any(!input_manifest$required %in% c("TRUE", "FALSE")) ||
          !all(vapply(input_manifest$relative_path, safe_relative, logical(1)))) {
        stop("The input manifest is malformed or contains an unsafe path.")
      }
      unknown_overrides <- setdiff(names(settings$file_overrides), input_manifest$input_id)
      if (length(unknown_overrides)) {
        stop("Unknown input override identifiers: ", paste(unknown_overrides, collapse = ", "))
      }
      emit(c("", "Source input checks"))
      for (i in seq_len(nrow(input_manifest))) {
        item <- input_manifest[i, , drop = FALSE]
        identifier <- item$input_id[[1L]]
        required <- identical(item$required[[1L]], "TRUE")
        path <- file.path(root, item$relative_path[[1L]])
        if (identifier %in% names(settings$file_overrides)) {
          path <- settings$file_overrides[[identifier]]
        }
        path <- normalizePath(path.expand(path), winslash = "/", mustWork = FALSE)
        emit(paste("Checking", identifier, "at", path))
        info <- file.info(path)
        readable <- file.exists(path) && !isTRUE(info$isdir) &&
          !is.na(info$size) && info$size > 0 && file.access(path, 4L) == 0L
        if (!readable) {
          record(if (required) "FAIL" else "WARN", identifier,
                 "The input is missing, empty, unreadable, or a directory.")
          next
        }
        expected_bytes <- suppressWarnings(as.numeric(item$expected_bytes[[1L]]))
        size_matches <- is.na(expected_bytes) || identical(as.numeric(info$size), expected_bytes)
        if (!size_matches) {
          record(if (required) "FAIL" else "WARN", identifier,
                 paste("Byte count differs: observed", format(info$size, scientific = FALSE),
                       "expected", format(expected_bytes, scientific = FALSE)))
        }
        observed_md5 <- if (settings$verify_md5) checksum(path) else "NOT_COMPUTED"
        expected_md5 <- item$expected_md5[[1L]]
        has_expected_md5 <- !is.na(expected_md5) && nzchar(expected_md5)
        emit(paste("Observed bytes:", format(info$size, scientific = FALSE),
                   "; observed MD5:", observed_md5))
        if (settings$verify_md5 && has_expected_md5) {
          if (!identical(observed_md5, tolower(expected_md5))) {
            record(if (required) "FAIL" else "WARN", identifier,
                   paste("Checksum differs; archived MD5:", expected_md5))
          } else if (size_matches) {
            record("PASS", identifier, "Readable and identical to the archived checksum.")
          }
        } else {
          record("WARN", identifier,
                 if (has_expected_md5) "Present; checksum comparison was skipped." else
                   "Present; no archived checksum is available for identity comparison.")
        }
      }

      emit(c("", "Author-provided expected table checks"))
      reference_zip <- file.path(project, "reference", "tables_reference.zip")
      target_map <- read_manifest(file.path(project, "config", "tables.tsv"))
      map_columns <- c("table", "part", "source_file", "worksheet", "expected_rows",
                       "expected_columns", "source_columns", "reference_zip_md5")
      if (!identical(names(target_map), map_columns) || nrow(target_map) != 7L) {
        stop("The table mapping must contain the seven supplied reference blocks.")
      }
      worksheets <- unique(target_map$worksheet)
      expected_worksheets <- paste0("Table_S", seq_len(6L))
      if (!identical(worksheets, expected_worksheets) ||
          any(nchar(worksheets) > 31L) ||
          anyDuplicated(tolower(worksheets)) ||
          anyDuplicated(target_map$source_file)) {
        stop("The table mapping does not define six unique expected worksheets.")
      }
      archive_ok <- identical(checksum(reference_zip), unique(target_map$reference_zip_md5))
      record(if (archive_ok) "PASS" else "FAIL", "reference_archive",
             if (archive_ok) "The uploaded table archive is preserved byte for byte." else
               "The reference archive differs from the author-provided archive.")
      listing <- utils::unzip(reference_zip, list = TRUE)$Name
      listing <- listing[!endsWith(listing, "/")]
      if (!setequal(listing, target_map$source_file) || anyDuplicated(listing)) {
        stop("The reference archive does not contain exactly the seven mapped tables.")
      }
      if (!all(vapply(target_map$source_file, safe_relative, logical(1)))) {
        stop("The reference archive mapping contains an unsafe member path.")
      }
      reference_directory <- tempfile(pattern = "ibd_reference_tables_")
      if (!dir.create(reference_directory)) {
        stop("Cannot create the temporary reference-table directory.")
      }
      on.exit(unlink(reference_directory, recursive = TRUE, force = TRUE), add = TRUE)
      utils::unzip(reference_zip, files = target_map$source_file,
                   exdir = reference_directory, unzip = "internal")
      reference_paths <- file.path(reference_directory, target_map$source_file)
      if (!all(file.exists(reference_paths))) {
        stop("At least one reference table could not be extracted for reading.")
      }
      for (i in seq_len(nrow(target_map))) {
        item <- target_map[i, , drop = FALSE]
        expected <- utils::read.delim(
          reference_paths[[i]], sep = "\t", header = TRUE, quote = "\"",
          comment.char = "", colClasses = "character",
          check.names = FALSE, stringsAsFactors = FALSE
        )
        expected_names <- strsplit(item$source_columns[[1L]], "|", fixed = TRUE)[[1L]]
        structure_ok <- nrow(expected) == as.integer(item$expected_rows[[1L]]) &&
          ncol(expected) == as.integer(item$expected_columns[[1L]]) &&
          identical(names(expected), expected_names)
        record(if (structure_ok) "PASS" else "FAIL", item$source_file[[1L]],
               paste("Rows:", nrow(expected), "; columns:", ncol(expected),
                     "; planned worksheet:", item$worksheet[[1L]]))
      }
      record("PASS", "workbook_structure",
             "Six planned worksheets; both S3 blocks remain on Table_S3 with their own columns.")
      emit("Reference values are comparison targets only; they are not analysis inputs.")

      emit(c("", "Software inventory"))
      target_versions <- c(edgeR = "4.4.2", limma = "3.62.2", Matrix = "1.7-1",
                           hdf5r = "1.3.12", data.table = "1.18.4",
                           statmod = "1.5.2", pROC = "1.19.0.1",
                           openxlsx = NA_character_, renv = NA_character_)
      required_packages <- setdiff(names(target_versions), "renv")
      if (!identical(as.character(getRversion()), "4.4.2")) {
        record("WARN", "R_version",
               paste("Observed", as.character(getRversion()), "; manuscript analyses used 4.4.2."))
      } else record("PASS", "R_version", "Matches the archived R version 4.4.2.")
      for (package in names(target_versions)) {
        available <- tryCatch(requireNamespace(package, quietly = TRUE),
                              error = function(e) FALSE)
        if (!available) {
          record("WARN", paste0("package_", package),
                 if (package %in% required_packages) {
                   "Not loadable; needed for later analysis or workbook generation."
                 } else "Not loadable; dependency locking is not yet available.")
          next
        }
        installed_version <- as.character(utils::packageVersion(package))
        archived_version <- target_versions[[package]]
        version_matches <- is.na(archived_version) ||
          utils::packageVersion(package) == package_version(archived_version)
        record(if (version_matches) "PASS" else "WARN", paste0("package_", package),
               paste("Observed version:", installed_version,
                     "; archived version:", if (is.na(archived_version)) "not specified" else archived_version))
      }
      emit(c("", "R session information", utils::capture.output(utils::sessionInfo())))
      final_status <- if (counts[["FAIL"]] > 0L) 1L else 0L
    }, error = function(error) {
      record("FAIL", "preflight_error", conditionMessage(error))
      final_status <<- 1L
    })
    emit(c("", paste("Checks: PASS =", counts[["PASS"]],
                      "; WARN =", counts[["WARN"]], "; FAIL =", counts[["FAIL"]]),
           paste("Readiness:", if (final_status != 0L) "BLOCKED" else
             if (counts[["WARN"]] > 0L) "REVIEW_REQUIRED" else "INPUT_CHECKS_PASSED"),
           "Scientific reproduction status: NOT_RUN",
           "Excel generation status: NOT_RUN",
           paste("Finished UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
           paste("Report:", normalizePath(report, winslash = "/", mustWork = TRUE))))
    final_status
  }

  arguments <- commandArgs(trailingOnly = TRUE)
  allowed <- arguments %in% c("--skip-md5", "--help") |
    startsWith(arguments, "--source-root=") | startsWith(arguments, "--stage=")
  if (any(!allowed)) stop("Unknown option: ", paste(arguments[!allowed], collapse = ", "))
  stage_option <- arguments[startsWith(arguments, "--stage=")]
  if (length(stage_option) > 1L) stop("Specify --stage only once.")
  stage <- if (length(stage_option)) sub("^--stage=", "", stage_option) else
    getOption("ibd.stage", "preflight")
  options(ibd.stage = NULL)
  if (!is.character(stage) || length(stage) != 1L || is.na(stage) ||
      !stage %in% c("preflight", "discovery", "replication", "development", "external", "sensitivity", "outputs", "release", "restore", "verify")) {
    stop("The stage must be preflight, discovery, replication, development, external, sensitivity, outputs, release, restore, or verify.")
  }
  if ("--help" %in% arguments) {
    cat("Usage: Rscript run.R [--stage=preflight|discovery|replication|development|external|sensitivity|outputs|release|restore|verify] [--source-root=PATH] [--skip-md5]\n",
        "The default action is an input preflight. Analysis stages require verified inputs.\n", sep = "")
  } else {
    project <- locate_project()
    if (identical(stage, "preflight")) {
      status <- preflight(project, arguments)
    } else {
      module_path <- file.path(project, "R", "discovery.R")
      if (!file.exists(module_path)) stop("Extract the complete package; missing R/discovery.R.")
      module <- new.env(parent = baseenv())
      sys.source(module_path, envir = module)
      context <- list(project = project, settings = load_settings(project, arguments),
                      read_manifest = read_manifest, safe_relative = safe_relative)
      if (identical(stage, "discovery")) {
        status <- module$run_discovery(context)
      } else {
        if (stage %in% c("external", "sensitivity", "outputs", "release")) {
          shared_path <- file.path(project, "R", "development.R")
          if (!file.exists(shared_path)) stop("Extract the complete package; missing R/development.R.")
          sys.source(shared_path, envir = module)
        }
        if (identical(stage, "sensitivity")) {
          shared_path <- file.path(project, "R", "external.R")
          if (!file.exists(shared_path)) stop("Extract the complete package; missing R/external.R.")
          sys.source(shared_path, envir = module)
        }
        module_name <- if (stage %in% c("release", "restore", "verify")) "outputs" else stage
        analysis_path <- file.path(project, "R", paste0(module_name, ".R"))
        if (!file.exists(analysis_path)) stop("Extract the complete package; missing R/", stage, ".R.")
        sys.source(analysis_path, envir = module)
        action <- get(paste0("run_", stage), envir = module, inherits = FALSE)
        status <- action(context)
      }
    }
    # Isolated verification workers must finish their post-run checks and
    # publish a closed completion record before the launcher exits.
    if (!interactive() && !isTRUE(getOption("ibd.worker.return.status", FALSE))) {
      quit(save = "no", status = status, runLast = FALSE)
    }
    invisible(status)
  }
})
