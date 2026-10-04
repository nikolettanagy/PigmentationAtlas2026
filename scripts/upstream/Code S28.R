#!/usr/bin/env Rscript

###############################################################################
# PigmentationAtlas
# Module 12.3 – Restart-safe full GWAS–apaQTL colocalization catalogue
#
# Strategy
#   - Reads the Module 12.1 ELIGIBLE work queue.
#   - Runs the validated Module 12.2 engine once per coloc_unit_id.
#   - Skips units already completed successfully unless --force is supplied.
#   - Copies unit-level outputs into a permanent catalogue directory.
#   - Rebuilds catalogue-wide summary, QC and run-status tables after every unit.
#
# Usage
#   Full catalogue:
#     Rscript scripts/module12_3_run_coloc_catalogue.R
#
#   Smoke test with first 10 pending units:
#     Rscript scripts/module12_3_run_coloc_catalogue.R --max-units=10
#
#   Re-run all selected units:
#     Rscript scripts/module12_3_run_coloc_catalogue.R --force
#
#   Start at a specific coloc unit:
#     Rscript scripts/module12_3_run_coloc_catalogue.R --start-from=COLOC_0000100
#
# Notes
#   - The run is intentionally sequential because Module 12.2 writes temporary
#     outputs to a shared single_test directory.
#   - A failed unit is recorded and the catalogue continues by default.
#   - Exit code is 0 if no selected unit failed; otherwise 1.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
})

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

module_start_time <- Sys.time()

###############################################################################
# 1. ARGUMENTS
###############################################################################

args <- commandArgs(trailingOnly = TRUE)

force_rerun <- "--force" %in% args
stop_on_error <- "--stop-on-error" %in% args

extract_arg <- function(prefix, default = NA_character_) {
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0L) {
    return(default)
  }
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

max_units_text <- extract_arg("--max-units=", NA_character_)
max_units <- suppressWarnings(as.integer(max_units_text))
if (is.na(max_units) || max_units < 1L) {
  max_units <- Inf
}

start_from <- extract_arg("--start-from=", NA_character_)
if (!is.na(start_from) && !nzchar(start_from)) {
  start_from <- NA_character_
}

###############################################################################
# 2. CONFIGURATION
###############################################################################

project_dir <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

module12_dir <- file.path(
  project_dir,
  "results",
  "module12",
  "colocalization"
)

work_queue_file <- file.path(
  module12_dir,
  "work_queue",
  "module12_1_coloc_work_queue_eligible.tsv.gz"
)

single_test_dir <- file.path(
  module12_dir,
  "single_test"
)

catalogue_dir <- file.path(
  module12_dir,
  "catalogue"
)

catalogue_units_dir <- file.path(
  catalogue_dir,
  "units"
)

catalogue_logs_dir <- file.path(
  catalogue_dir,
  "logs"
)

catalogue_status_dir <- file.path(
  catalogue_dir,
  "status"
)

catalogue_summary_dir <- file.path(
  catalogue_dir,
  "summary"
)

dir.create(catalogue_units_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(catalogue_logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(catalogue_status_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(catalogue_summary_dir, recursive = TRUE, showWarnings = FALSE)

run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

run_log_file <- file.path(
  catalogue_logs_dir,
  paste0("module12_3_catalogue_", run_timestamp, ".log")
)

run_status_file <- file.path(
  catalogue_status_dir,
  "module12_3_unit_status.tsv"
)

catalogue_summary_file <- file.path(
  catalogue_summary_dir,
  "module12_3_coloc_catalogue_summary.tsv.gz"
)

catalogue_qc_file <- file.path(
  catalogue_summary_dir,
  "module12_3_coloc_catalogue_qc.tsv.gz"
)

catalogue_interpretation_file <- file.path(
  catalogue_summary_dir,
  "module12_3_coloc_interpretation_counts.tsv"
)

module_status_file <- file.path(
  catalogue_summary_dir,
  "module12_3_status.tsv"
)

###############################################################################
# 3. HELPERS
###############################################################################

log_message <- function(...) {
  text <- paste0(..., collapse = "")
  stamped <- paste0(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    text
  )
  cat(stamped, "\n")
  cat(stamped, "\n", file = run_log_file, append = TRUE)
}

stop_module <- function(message_text) {
  log_message("FATAL: ", message_text)
  stop(message_text, call. = FALSE)
}

first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0L) {
    return(NA_character_)
  }
  normalizePath(hit[[1L]], winslash = "/", mustWork = TRUE)
}

locate_module12_2 <- function() {
  candidates <- c(
    file.path(project_dir, "scripts", "module12_2_prepare_and_test_single_coloc.R"),
    file.path(project_dir, "scripts", "module12_2_prepare_and_test_single_coloc_FINAL_FIXED_v4.R"),
    file.path(project_dir, "module12_2_prepare_and_test_single_coloc.R"),
    file.path(project_dir, "module12_2_prepare_and_test_single_coloc_FINAL_FIXED_v4.R")
  )

  found <- first_existing_file(candidates)

  if (is.na(found)) {
    stop_module(
      paste0(
        "Could not locate the validated Module 12.2 script.\n",
        "Expected one of:\n  ",
        paste(candidates, collapse = "\n  ")
      )
    )
  }

  found
}

locate_rscript <- function() {
  from_path <- Sys.which("Rscript")
  if (nzchar(from_path)) {
    return(normalizePath(from_path, winslash = "/", mustWork = TRUE))
  }

  candidates <- c(
    file.path(R.home("bin"), "Rscript.exe"),
    file.path(R.home("bin"), "Rscript")
  )

  found <- first_existing_file(candidates)
  if (is.na(found)) {
    stop_module("Could not locate Rscript.")
  }

  found
}

safe_fread <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }

  tryCatch(
    fread(path, showProgress = FALSE),
    error = function(e) NULL
  )
}

copy_if_exists <- function(source, destination) {
  if (!file.exists(source)) {
    return(FALSE)
  }

  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)

  copied <- file.copy(
    from = source,
    to = destination,
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )

  isTRUE(copied)
}

latest_matching_file <- function(directory, pattern) {
  if (!dir.exists(directory)) {
    return(NA_character_)
  }

  files <- list.files(
    directory,
    pattern = pattern,
    full.names = TRUE
  )

  if (length(files) == 0L) {
    return(NA_character_)
  }

  info <- file.info(files)
  files[[which.max(info$mtime)]]
}

unit_catalogue_status_path <- function(unit_id) {
  file.path(
    catalogue_units_dir,
    unit_id,
    paste0(unit_id, "_catalogue_status.tsv")
  )
}

unit_is_complete <- function(unit_id) {
  status_path <- unit_catalogue_status_path(unit_id)
  status <- safe_fread(status_path)

  if (is.null(status) || nrow(status) == 0L) {
    return(FALSE)
  }

  required <- c("coloc_unit_id", "catalogue_status")
  if (!all(required %in% names(status))) {
    return(FALSE)
  }

  any(
    status$coloc_unit_id == unit_id &
      status$catalogue_status %in% c("PASS", "PASS_WITH_WARNINGS")
  )
}

write_run_status <- function(status_dt) {
  setorder(status_dt, queue_order)
  fwrite(
    status_dt,
    run_status_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}

rebuild_catalogue_tables <- function() {
  unit_dirs <- list.dirs(
    catalogue_units_dir,
    recursive = FALSE,
    full.names = TRUE
  )

  summary_files <- file.path(
    unit_dirs,
    basename(unit_dirs),
    paste0(basename(unit_dirs), "_coloc_summary.tsv")
  )

  qc_files <- file.path(
    unit_dirs,
    basename(unit_dirs),
    paste0(basename(unit_dirs), "_qc.tsv")
  )

  summary_list <- lapply(summary_files[file.exists(summary_files)], safe_fread)
  summary_list <- Filter(Negate(is.null), summary_list)

  qc_list <- lapply(qc_files[file.exists(qc_files)], safe_fread)
  qc_list <- Filter(Negate(is.null), qc_list)

  if (length(summary_list) > 0L) {
    all_summary <- rbindlist(summary_list, fill = TRUE, use.names = TRUE)

    if ("coloc_unit_id" %in% names(all_summary)) {
      setorder(all_summary, coloc_unit_id)
      all_summary <- unique(all_summary, by = "coloc_unit_id")
    }

    fwrite(
      all_summary,
      catalogue_summary_file,
      sep = "\t",
      quote = FALSE,
      na = "NA",
      compress = "gzip"
    )

    if ("interpretation" %in% names(all_summary)) {
      interpretation_counts <- all_summary[
        ,
        .N,
        by = interpretation
      ][order(-N, interpretation)]

      fwrite(
        interpretation_counts,
        catalogue_interpretation_file,
        sep = "\t",
        quote = FALSE,
        na = "NA"
      )
    }
  }

  if (length(qc_list) > 0L) {
    all_qc <- rbindlist(qc_list, fill = TRUE, use.names = TRUE)

    if ("coloc_unit_id" %in% names(all_qc)) {
      setorder(all_qc, coloc_unit_id)
      all_qc <- unique(all_qc, by = "coloc_unit_id")
    }

    fwrite(
      all_qc,
      catalogue_qc_file,
      sep = "\t",
      quote = FALSE,
      na = "NA",
      compress = "gzip"
    )
  }
}

###############################################################################
# 4. INPUT VALIDATION
###############################################################################

if (!file.exists(work_queue_file)) {
  stop_module(
    paste0(
      "Eligible work queue not found:\n",
      work_queue_file
    )
  )
}

module12_2_script <- locate_module12_2()
rscript_executable <- locate_rscript()

queue <- fread(work_queue_file, showProgress = FALSE)

if (!"coloc_unit_id" %in% names(queue)) {
  stop_module("The eligible queue has no coloc_unit_id column.")
}

queue[, coloc_unit_id := as.character(coloc_unit_id)]
queue <- queue[!is.na(coloc_unit_id) & nzchar(coloc_unit_id)]
queue <- unique(queue, by = "coloc_unit_id")
queue[, queue_order := .I]

if (!is.na(start_from)) {
  start_index <- match(start_from, queue$coloc_unit_id)

  if (is.na(start_index)) {
    stop_module(
      paste0(
        "--start-from unit not present in eligible queue: ",
        start_from
      )
    )
  }

  queue <- queue[queue_order >= start_index]
}

if (!force_rerun) {
  completed <- vapply(
    queue$coloc_unit_id,
    unit_is_complete,
    logical(1)
  )
  queue <- queue[!completed]
}

if (is.finite(max_units) && nrow(queue) > max_units) {
  queue <- queue[seq_len(max_units)]
}

###############################################################################
# 5. INITIAL STATUS
###############################################################################

cat(
  "============================================================\n",
  "PigmentationAtlas\n",
  "MODULE 12.3 – FULL COLOCALIZATION CATALOGUE\n",
  "============================================================\n\n",
  sep = ""
)

log_message("Validated Module 12.2 engine: ", module12_2_script)
log_message("Rscript executable: ", rscript_executable)
log_message("Eligible queue file: ", work_queue_file)
log_message("Units selected for this invocation: ", nrow(queue))
log_message("Force rerun: ", force_rerun)
log_message("Stop on error: ", stop_on_error)

if (nrow(queue) == 0L) {
  log_message("No pending units. Rebuilding catalogue summaries only.")
  rebuild_catalogue_tables()

  final_status <- data.table(
    module = "12.3",
    start_time = format(module_start_time, "%Y-%m-%d %H:%M:%S"),
    end_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    selected_units = 0L,
    successful_units = 0L,
    failed_units = 0L,
    skipped_units = 0L,
    final_status = "PASS_NO_PENDING_UNITS"
  )

  fwrite(
    final_status,
    module_status_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  quit(save = "no", status = 0L)
}

existing_status <- safe_fread(run_status_file)

if (is.null(existing_status)) {
  run_status <- data.table(
    queue_order = integer(),
    coloc_unit_id = character(),
    clump_id = character(),
    tissue = character(),
    phenotype_id = character(),
    gene_id = character(),
    started_at = character(),
    ended_at = character(),
    runtime_seconds = numeric(),
    child_exit_code = integer(),
    catalogue_status = character(),
    message = character(),
    unit_directory = character(),
    child_log = character()
  )
} else {
  run_status <- existing_status
}

###############################################################################
# 6. UNIT LOOP
###############################################################################

success_count <- 0L
failure_count <- 0L
warning_count <- 0L

for (i in seq_len(nrow(queue))) {
  unit <- queue[i]
  unit_id <- unit$coloc_unit_id[[1L]]

  unit_started <- Sys.time()

  clump_id <- if ("clump_id" %in% names(unit)) {
    as.character(unit$clump_id[[1L]])
  } else {
    NA_character_
  }

  tissue <- if ("tissue" %in% names(unit)) {
    as.character(unit$tissue[[1L]])
  } else {
    NA_character_
  }

  phenotype_id <- if ("phenotype_id" %in% names(unit)) {
    as.character(unit$phenotype_id[[1L]])
  } else {
    NA_character_
  }

  gene_id <- if ("gene_id" %in% names(unit)) {
    as.character(unit$gene_id[[1L]])
  } else {
    NA_character_
  }

  log_message(
    "[",
    i,
    "/",
    nrow(queue),
    "] Starting ",
    unit_id,
    " | ",
    clump_id,
    " | ",
    tissue,
    " | ",
    phenotype_id
  )

  child_log <- file.path(
    catalogue_logs_dir,
    paste0(unit_id, "_", run_timestamp, ".log")
  )

  command_output <- tryCatch(
    system2(
      command = rscript_executable,
      args = c(
        shQuote(module12_2_script),
        shQuote(unit_id)
      ),
      stdout = child_log,
      stderr = child_log,
      wait = TRUE
    ),
    error = function(e) {
      attr(1L, "system_error") <- conditionMessage(e)
      1L
    }
  )

  child_exit_code <- as.integer(command_output)
  unit_ended <- Sys.time()
  runtime_seconds <- as.numeric(
    difftime(unit_ended, unit_started, units = "secs")
  )

  source_result_dir <- file.path(
    single_test_dir,
    "results",
    paste0(unit_id, "_", clump_id)
  )

  source_summary <- file.path(
    source_result_dir,
    paste0(unit_id, "_coloc_summary.tsv")
  )

  source_snp <- file.path(
    source_result_dir,
    paste0(unit_id, "_coloc_snp_results.tsv.gz")
  )

  source_status <- file.path(
    source_result_dir,
    paste0(unit_id, "_status.tsv")
  )

  source_qc <- file.path(
    single_test_dir,
    "qc",
    paste0(unit_id, "_qc.tsv")
  )

  source_harmonized <- file.path(
    single_test_dir,
    "input",
    paste0(unit_id, "_harmonized_input.tsv.gz")
  )

  source_module_log <- latest_matching_file(
    file.path(single_test_dir, "logs"),
    paste0(
      "^",
      unit_id,
      ".*\\.log$|^module12_2_",
      unit_id,
      "_.*\\.log$"
    )
  )

  unit_dir <- file.path(catalogue_units_dir, unit_id)
  dir.create(unit_dir, recursive = TRUE, showWarnings = FALSE)

  destination_summary <- file.path(
    unit_dir,
    paste0(unit_id, "_coloc_summary.tsv")
  )

  destination_snp <- file.path(
    unit_dir,
    paste0(unit_id, "_coloc_snp_results.tsv.gz")
  )

  destination_status <- file.path(
    unit_dir,
    paste0(unit_id, "_module12_2_status.tsv")
  )

  destination_qc <- file.path(
    unit_dir,
    paste0(unit_id, "_qc.tsv")
  )

  destination_harmonized <- file.path(
    unit_dir,
    paste0(unit_id, "_harmonized_input.tsv.gz")
  )

  destination_module_log <- file.path(
    unit_dir,
    paste0(unit_id, "_module12_2.log")
  )

  copied_summary <- copy_if_exists(source_summary, destination_summary)
  copied_snp <- copy_if_exists(source_snp, destination_snp)
  copied_status <- copy_if_exists(source_status, destination_status)
  copied_qc <- copy_if_exists(source_qc, destination_qc)
  copied_harmonized <- copy_if_exists(
    source_harmonized,
    destination_harmonized
  )

  if (!is.na(source_module_log)) {
    copy_if_exists(source_module_log, destination_module_log)
  }

  source_status_dt <- safe_fread(source_status)
  module12_2_final_status <- NA_character_

  if (
    !is.null(source_status_dt) &&
      nrow(source_status_dt) > 0L &&
      "final_status" %in% names(source_status_dt)
  ) {
    module12_2_final_status <- as.character(
      source_status_dt$final_status[[1L]]
    )
  }

  required_outputs_present <- all(
    copied_summary,
    copied_snp,
    copied_status,
    copied_qc,
    copied_harmonized
  )

  if (
    child_exit_code == 0L &&
      required_outputs_present &&
      module12_2_final_status %in% c("PASS", "PASS_WITH_WARNINGS")
  ) {
    catalogue_status <- module12_2_final_status
    message_text <- "Completed and copied to catalogue."

    success_count <- success_count + 1L
    if (identical(catalogue_status, "PASS_WITH_WARNINGS")) {
      warning_count <- warning_count + 1L
    }
  } else {
    catalogue_status <- "FAIL"
    failure_count <- failure_count + 1L

    missing_outputs <- c(
      summary = !copied_summary,
      snp_results = !copied_snp,
      status = !copied_status,
      qc = !copied_qc,
      harmonized_input = !copied_harmonized
    )

    missing_names <- names(missing_outputs)[missing_outputs]

    message_parts <- c(
      paste0("child_exit_code=", child_exit_code),
      paste0(
        "module12_2_final_status=",
        ifelse(
          is.na(module12_2_final_status),
          "NA",
          module12_2_final_status
        )
      )
    )

    if (length(missing_names) > 0L) {
      message_parts <- c(
        message_parts,
        paste0(
          "missing_outputs=",
          paste(missing_names, collapse = ",")
        )
      )
    }

    system_error <- attr(command_output, "system_error")
    if (!is.null(system_error)) {
      message_parts <- c(
        message_parts,
        paste0("system_error=", system_error)
      )
    }

    message_text <- paste(message_parts, collapse = "; ")
  }

  unit_catalogue_status <- data.table(
    module = "12.3",
    coloc_unit_id = unit_id,
    queue_order = unit$queue_order[[1L]],
    clump_id = clump_id,
    tissue = tissue,
    phenotype_id = phenotype_id,
    gene_id = gene_id,
    started_at = format(unit_started, "%Y-%m-%d %H:%M:%S"),
    ended_at = format(unit_ended, "%Y-%m-%d %H:%M:%S"),
    runtime_seconds = runtime_seconds,
    child_exit_code = child_exit_code,
    module12_2_final_status = module12_2_final_status,
    catalogue_status = catalogue_status,
    message = message_text,
    unit_directory = normalizePath(
      unit_dir,
      winslash = "/",
      mustWork = FALSE
    ),
    child_log = normalizePath(
      child_log,
      winslash = "/",
      mustWork = FALSE
    )
  )

  fwrite(
    unit_catalogue_status,
    unit_catalogue_status_path(unit_id),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  run_status <- run_status[
    coloc_unit_id != unit_id
  ]

  run_status <- rbind(
    run_status,
    unit_catalogue_status[
      ,
      .(
        queue_order,
        coloc_unit_id,
        clump_id,
        tissue,
        phenotype_id,
        gene_id,
        started_at,
        ended_at,
        runtime_seconds,
        child_exit_code,
        catalogue_status,
        message,
        unit_directory,
        child_log
      )
    ],
    fill = TRUE
  )

  write_run_status(run_status)
  rebuild_catalogue_tables()

  log_message(
    "[",
    i,
    "/",
    nrow(queue),
    "] ",
    unit_id,
    " -> ",
    catalogue_status,
    " (",
    round(runtime_seconds, 2),
    " s)"
  )

  if (
    identical(catalogue_status, "FAIL") &&
      stop_on_error
  ) {
    log_message("Stopping because --stop-on-error was supplied.")
    break
  }
}

###############################################################################
# 7. FINAL SUMMARY
###############################################################################

module_end_time <- Sys.time()

selected_status <- run_status[
  coloc_unit_id %in% queue$coloc_unit_id
]

successful_units <- selected_status[
  catalogue_status %in% c("PASS", "PASS_WITH_WARNINGS"),
  .N
]

failed_units <- selected_status[
  catalogue_status == "FAIL",
  .N
]

warning_units <- selected_status[
  catalogue_status == "PASS_WITH_WARNINGS",
  .N
]

final_module_status <- if (failed_units == 0L) {
  if (warning_units > 0L) {
    "PASS_WITH_WARNINGS"
  } else {
    "PASS"
  }
} else if (successful_units > 0L) {
  "PARTIAL_PASS"
} else {
  "FAIL"
}

module_status <- data.table(
  module = "12.3",
  start_time = format(module_start_time, "%Y-%m-%d %H:%M:%S"),
  end_time = format(module_end_time, "%Y-%m-%d %H:%M:%S"),
  runtime_seconds = as.numeric(
    difftime(module_end_time, module_start_time, units = "secs")
  ),
  eligible_units_in_queue = nrow(
    fread(work_queue_file, select = "coloc_unit_id", showProgress = FALSE)
  ),
  selected_units = nrow(queue),
  successful_units = successful_units,
  warning_units = warning_units,
  failed_units = failed_units,
  force_rerun = force_rerun,
  stop_on_error = stop_on_error,
  final_status = final_module_status
)

fwrite(
  module_status,
  module_status_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat(
  "\n============================================================\n",
  "MODULE 12.3 – FINAL SUMMARY\n",
  "============================================================\n\n",
  sep = ""
)

print(module_status)

cat(
  "\nOutputs:\n",
  "  Unit status:\n    ",
  run_status_file,
  "\n",
  "  Catalogue summary:\n    ",
  catalogue_summary_file,
  "\n",
  "  Catalogue QC:\n    ",
  catalogue_qc_file,
  "\n",
  "  Interpretation counts:\n    ",
  catalogue_interpretation_file,
  "\n",
  "  Module status:\n    ",
  module_status_file,
  "\n",
  "  Run log:\n    ",
  run_log_file,
  "\n",
  sep = ""
)

cat(
  "\n============================================================\n",
  "FINAL RESULT: ",
  final_module_status,
  "\n",
  "============================================================\n",
  sep = ""
)

quit(
  save = "no",
  status = if (failed_units == 0L) 0L else 1L
)
