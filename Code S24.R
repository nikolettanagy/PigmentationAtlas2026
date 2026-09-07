###############################################################################
# PigmentationAtlas
# Module 11.3B.3
#
# BUILD THE HARMONIZED apaQTL CATALOGUE
#
# Input:
#   results/module11/harmonization/catalogue/
#   module11_3B_work_queue.tsv
#
# Function source:
#   scripts/module11_3B_2_read_and_qc_apaqtl.R
#
# Main outputs:
#   module11_3B_harmonized_apaQTL.tsv.gz
#   module11_3B_file_summary.tsv
#   module11_3B_failed_files.tsv
#   module11_3B_column_mapping_summary.tsv
#   module11_3B_build_status.tsv
#
# Restart behaviour:
#   Existing valid chunk files are skipped.
###############################################################################


###############################################################################
# 1. CLEAN SESSION
###############################################################################

rm(list = ls())
invisible(gc())

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  width = 220
)

module_start_time <- Sys.time()

cat("\n")
cat("============================================================\n")
cat("PigmentationAtlas\n")
cat("MODULE 11.3B.3 â€“ HARMONIZED apaQTL CATALOGUE BUILDER\n")
cat("Started:", format(module_start_time), "\n")
cat("============================================================\n\n")


###############################################################################
# 2. LOAD PACKAGE
###############################################################################

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop(
    "Package 'data.table' is required.",
    call. = FALSE
  )
}

suppressPackageStartupMessages(
  library(data.table)
)


###############################################################################
# 3. PROJECT PATHS
###############################################################################

project_root <- normalizePath(
  "C:/Users/User/Desktop/PigmentationAtlas",
  winslash = "/",
  mustWork = TRUE
)

setwd(project_root)

reader_script <- file.path(
  project_root,
  "scripts",
  "module11_3B_2_read_and_qc_apaqtl.R"
)

catalogue_dir <- file.path(
  project_root,
  "results",
  "module11",
  "harmonization",
  "catalogue"
)

chunk_dir <- file.path(
  catalogue_dir,
  "chunks"
)

checkpoint_dir <- file.path(
  catalogue_dir,
  "checkpoints"
)

qc_dir <- file.path(
  catalogue_dir,
  "qc"
)

log_dir <- file.path(
  catalogue_dir,
  "logs"
)

dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

work_queue_file <- file.path(
  catalogue_dir,
  "module11_3B_work_queue.tsv"
)

final_catalogue_file <- file.path(
  catalogue_dir,
  "module11_3B_harmonized_apaQTL.tsv.gz"
)

file_summary_output <- file.path(
  catalogue_dir,
  "module11_3B_file_summary.tsv"
)

failed_files_output <- file.path(
  catalogue_dir,
  "module11_3B_failed_files.tsv"
)

column_mapping_output <- file.path(
  catalogue_dir,
  "module11_3B_column_mapping_summary.tsv"
)

progress_output <- file.path(
  checkpoint_dir,
  "module11_3B_progress.tsv"
)

build_status_output <- file.path(
  catalogue_dir,
  "module11_3B_build_status.tsv"
)

session_info_output <- file.path(
  log_dir,
  "module11_3B_3_session_info.txt"
)

log_file <- file.path(
  log_dir,
  paste0(
    "module11_3B_3_build_",
    format(module_start_time, "%Y%m%d_%H%M%S"),
    ".log"
  )
)

###############################################################################
# 4. LOAD THE VALIDATED READER FUNCTION IN AN ISOLATED ENVIRONMENT
###############################################################################

if (!file.exists(reader_script)) {
  stop(
    paste0(
      "Reader script not found:\n",
      reader_script
    ),
    call. = FALSE
  )
}

cat("Loading the validated Module 11.3B.2 reader...\n\n")

reader_environment <- new.env(
  parent = globalenv()
)

source(
  reader_script,
  local = reader_environment,
  echo = FALSE
)

if (
  !exists(
    "read_and_qc_apaqtl",
    envir = reader_environment,
    mode = "function",
    inherits = FALSE
  )
) {
  stop(
    "The function read_and_qc_apaqtl() was not created.",
    call. = FALSE
  )
}

read_and_qc_apaqtl <- get(
  "read_and_qc_apaqtl",
  envir = reader_environment,
  inherits = FALSE
)

rm(reader_environment)
invisible(gc())

cat("\nValidated reader function loaded successfully.\n\n")

###############################################################################
# 5. HELPER FUNCTIONS
###############################################################################

safe_filename <- function(x) {

  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)

  x
}


empty_dt <- function(columns) {

  as.data.table(
    setNames(
      replicate(
        length(columns),
        character(0),
        simplify = FALSE
      ),
      columns
    )
  )
}


write_empty_or_table <- function(x,
                                 output_file,
                                 columns = NULL) {

  if (nrow(x) == 0L && !is.null(columns)) {
    x <- empty_dt(columns)
  }

  fwrite(
    x,
    output_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}


append_binary_file <- function(source_file,
                               destination_file) {

  input_connection <- file(
    source_file,
    open = "rb"
  )

  output_connection <- file(
    destination_file,
    open = "ab"
  )

  on.exit(
    {
      try(close(input_connection), silent = TRUE)
      try(close(output_connection), silent = TRUE)
    },
    add = TRUE
  )

  repeat {

    buffer <- readBin(
      input_connection,
      what = "raw",
      n = 1024L * 1024L
    )

    if (length(buffer) == 0L) {
      break
    }

    writeBin(
      buffer,
      output_connection
    )
  }

  close(input_connection)
  close(output_connection)

  TRUE
}


make_chunk_path <- function(queue_id,
                            clump_id,
                            tissue) {

  file.path(
    chunk_dir,
    paste0(
      safe_filename(queue_id),
      "__",
      safe_filename(clump_id),
      "__",
      safe_filename(tissue),
      ".harmonized.tsv.gz"
    )
  )
}


make_qc_path <- function(queue_id) {

  file.path(
    checkpoint_dir,
    paste0(
      safe_filename(queue_id),
      ".qc.tsv"
    )
  )
}


make_column_map_path <- function(queue_id) {

  file.path(
    checkpoint_dir,
    paste0(
      safe_filename(queue_id),
      ".column_map.tsv"
    )
  )
}


make_failure_path <- function(queue_id) {

  file.path(
    checkpoint_dir,
    paste0(
      safe_filename(queue_id),
      ".failure.tsv"
    )
  )
}


###############################################################################
# 6. READ WORK QUEUE
###############################################################################

if (!file.exists(work_queue_file)) {
  stop(
    paste0(
      "Work queue not found:\n",
      work_queue_file
    ),
    call. = FALSE
  )
}

work_queue <- fread(
  work_queue_file,
  sep = "\t",
  header = TRUE,
  na.strings = c("", "NA", "NaN", "NULL"),
  showProgress = FALSE
)

required_columns <- c(
  "queue_id",
  "clump_id",
  "tissue",
  "canonical_chromosome",
  "canonical_region_start",
  "canonical_region_end",
  "output_file_resolved"
)

missing_columns <- setdiff(
  required_columns,
  names(work_queue)
)

if (length(missing_columns) > 0L) {
  stop(
    paste0(
      "Missing work-queue columns: ",
      paste(missing_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

if (nrow(work_queue) != 516L) {

  warning(
    paste0(
      "Expected 516 work-queue rows, observed ",
      nrow(work_queue),
      "."
    ),
    call. = FALSE
  )
}

cat("Work queue rows:", nrow(work_queue), "\n")
cat("Unique CLUMPs:", uniqueN(work_queue$clump_id), "\n")
cat("Unique tissues:", uniqueN(work_queue$tissue), "\n\n")


###############################################################################
# 7. INITIALIZE LOG
###############################################################################

log_connection <- file(
  log_file,
  open = "wt"
)

sink(
  log_connection,
  type = "output",
  split = TRUE
)

logging_active <- TRUE

on.exit(
  {

    if (
      exists("logging_active", inherits = FALSE) &&
      isTRUE(logging_active)
    ) {

      while (sink.number(type = "output") > 0L) {
        try(sink(type = "output"), silent = TRUE)
      }

      try(
        close(log_connection),
        silent = TRUE
      )
    }
  },
  add = TRUE
)

cat("Log file:\n", log_file, "\n\n")


###############################################################################
# 8. PROCESS ALL FILES
###############################################################################

total_files <- nrow(work_queue)

run_progress <- vector(
  mode = "list",
  length = total_files
)

for (i in seq_len(total_files)) {

  metadata <- work_queue[i]

  queue_id <- metadata$queue_id[1L]
  clump_id <- metadata$clump_id[1L]
  tissue <- metadata$tissue[1L]
  source_file <- metadata$output_file_resolved[1L]

  chunk_file <- make_chunk_path(
    queue_id,
    clump_id,
    tissue
  )

  qc_file <- make_qc_path(queue_id)
  column_map_file <- make_column_map_path(queue_id)
  failure_file <- make_failure_path(queue_id)

  cat(
    sprintf(
      "[%03d/%03d] %s | %s | %s\n",
      i,
      total_files,
      queue_id,
      clump_id,
      tissue
    )
  )

  # Restart-safe skip.
  previously_completed <- (
    file.exists(chunk_file) &&
      file.info(chunk_file)$size > 0 &&
      file.exists(qc_file) &&
      file.info(qc_file)$size > 0
  )

  if (previously_completed) {

    previous_qc <- tryCatch(
      fread(
        qc_file,
        sep = "\t",
        header = TRUE,
        showProgress = FALSE
      ),
      error = function(e) NULL
    )

    if (
      !is.null(previous_qc) &&
      nrow(previous_qc) == 1L &&
      previous_qc$file_status[1L] %chin%
        c("PASS", "PASS_WITH_ROW_QC_FLAGS")
    ) {

      cat("  Status: SKIPPED_ALREADY_COMPLETED\n")

      run_progress[[i]] <- data.table(
        queue_id = queue_id,
        clump_id = clump_id,
        tissue = tissue,
        processing_status = "SKIPPED_ALREADY_COMPLETED",
        chunk_file = chunk_file,
        qc_file = qc_file,
        error_message = NA_character_
      )

      next
    }
  }

  file_result <- tryCatch(

    read_and_qc_apaqtl(
      file_path = source_file,
      metadata = metadata,
      retain_source_columns = FALSE,
      fail_on_missing_pvalue = TRUE,
      fail_on_missing_variant = TRUE
    ),

    error = function(e) e
  )

  if (inherits(file_result, "error")) {

    error_message <- conditionMessage(file_result)

    failure_record <- data.table(
      queue_id = queue_id,
      clump_id = clump_id,
      tissue = tissue,
      source_file = source_file,
      error_message = error_message,
      failure_time = format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    )

    fwrite(
      failure_record,
      failure_file,
      sep = "\t",
      quote = FALSE,
      na = "NA"
    )

    cat("  Status: FAILED\n")
    cat("  Error:", error_message, "\n")

    run_progress[[i]] <- data.table(
      queue_id = queue_id,
      clump_id = clump_id,
      tissue = tissue,
      processing_status = "FAILED",
      chunk_file = NA_character_,
      qc_file = NA_character_,
      error_message = error_message
    )

    next
  }

  fwrite(
    file_result$data,
    chunk_file,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = "gzip"
  )

  fwrite(
    file_result$qc,
    qc_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  column_map_with_metadata <- copy(
    file_result$column_map
  )

  column_map_with_metadata[
    ,
    `:=`(
      queue_id = queue_id,
      clump_id = clump_id,
      tissue = tissue
    )
  ]

  setcolorder(
    column_map_with_metadata,
    c(
      "queue_id",
      "clump_id",
      "tissue",
      "standardized_field",
      "source_column",
      "detected"
    )
  )

  fwrite(
    column_map_with_metadata,
    column_map_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  if (file.exists(failure_file)) {
    unlink(failure_file)
  }

  cat(
    "  Status:",
    file_result$qc$file_status,
    "| rows:",
    format(nrow(file_result$data), big.mark = ","),
    "| runtime:",
    file_result$qc$runtime_seconds,
    "seconds\n"
  )

  run_progress[[i]] <- data.table(
    queue_id = queue_id,
    clump_id = clump_id,
    tissue = tissue,
    processing_status = file_result$qc$file_status,
    chunk_file = chunk_file,
    qc_file = qc_file,
    error_message = NA_character_
  )

  if (i %% 10L == 0L || i == total_files) {

    progress_table <- rbindlist(
      run_progress[
        !vapply(
          run_progress,
          is.null,
          logical(1)
        )
      ],
      use.names = TRUE,
      fill = TRUE
    )

    fwrite(
      progress_table,
      progress_output,
      sep = "\t",
      quote = FALSE,
      na = "NA"
    )

    invisible(gc())

    cat(
      "  Progress checkpoint saved:",
      nrow(progress_table),
      "/",
      total_files,
      "\n"
    )
  }
}


###############################################################################
# 9. FINAL PROGRESS TABLE
###############################################################################

progress_table <- rbindlist(
  run_progress[
    !vapply(
      run_progress,
      is.null,
      logical(1)
    )
  ],
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  progress_table,
  progress_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


###############################################################################
# 10. COMBINE FILE-LEVEL QC TABLES
###############################################################################

qc_files <- list.files(
  checkpoint_dir,
  pattern = "^APAQTL_[0-9]+\\.qc\\.tsv$",
  full.names = TRUE
)

qc_tables <- lapply(
  qc_files,
  function(x) {

    tryCatch(
      fread(
        x,
        sep = "\t",
        header = TRUE,
        showProgress = FALSE
      ),
      error = function(e) NULL
    )
  }
)

qc_tables <- Filter(
  Negate(is.null),
  qc_tables
)

if (length(qc_tables) > 0L) {

  file_summary <- rbindlist(
    qc_tables,
    use.names = TRUE,
    fill = TRUE
  )

  setorder(
    file_summary,
    queue_id
  )

} else {

  file_summary <- data.table()
}

fwrite(
  file_summary,
  file_summary_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


###############################################################################
# 11. COMBINE COLUMN-MAPPING TABLES
###############################################################################

column_map_files <- list.files(
  checkpoint_dir,
  pattern = "^APAQTL_[0-9]+\\.column_map\\.tsv$",
  full.names = TRUE
)

column_map_tables <- lapply(
  column_map_files,
  function(x) {

    tryCatch(
      fread(
        x,
        sep = "\t",
        header = TRUE,
        showProgress = FALSE
      ),
      error = function(e) NULL
    )
  }
)

column_map_tables <- Filter(
  Negate(is.null),
  column_map_tables
)

if (length(column_map_tables) > 0L) {

  column_mapping_summary <- rbindlist(
    column_map_tables,
    use.names = TRUE,
    fill = TRUE
  )

  setorder(
    column_mapping_summary,
    queue_id,
    standardized_field
  )

} else {

  column_mapping_summary <- data.table()
}

fwrite(
  column_mapping_summary,
  column_mapping_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


###############################################################################
# 12. COMBINE FAILURE TABLES
###############################################################################

failure_files <- list.files(
  checkpoint_dir,
  pattern = "^APAQTL_[0-9]+\\.failure\\.tsv$",
  full.names = TRUE
)

failure_tables <- lapply(
  failure_files,
  function(x) {

    tryCatch(
      fread(
        x,
        sep = "\t",
        header = TRUE,
        showProgress = FALSE
      ),
      error = function(e) NULL
    )
  }
)

failure_tables <- Filter(
  Negate(is.null),
  failure_tables
)

if (length(failure_tables) > 0L) {

  failed_files <- rbindlist(
    failure_tables,
    use.names = TRUE,
    fill = TRUE
  )

  setorder(
    failed_files,
    queue_id
  )

} else {

  failed_files <- data.table(
    queue_id = character(),
    clump_id = character(),
    tissue = character(),
    source_file = character(),
    error_message = character(),
    failure_time = character()
  )
}

fwrite(
  failed_files,
  failed_files_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


###############################################################################
# 13. VALIDATE COMPLETED CHUNKS
###############################################################################

successful_queue_ids <- file_summary[
  file_status %chin% c(
    "PASS",
    "PASS_WITH_ROW_QC_FLAGS"
  ),
  queue_id
]

expected_successful_chunks <- work_queue[
  queue_id %chin% successful_queue_ids
]

expected_successful_chunks[
  ,
  chunk_file := mapply(
    make_chunk_path,
    queue_id,
    clump_id,
    tissue,
    USE.NAMES = FALSE
  )
]

expected_successful_chunks[
  ,
  chunk_exists := file.exists(chunk_file)
]

expected_successful_chunks[
  ,
  chunk_nonempty := chunk_exists &
    file.info(chunk_file)$size > 0
]

valid_chunk_table <- expected_successful_chunks[
  chunk_exists == TRUE &
    chunk_nonempty == TRUE
]

setorder(
  valid_chunk_table,
  queue_id
)


###############################################################################
# 14. BUILD THE SINGLE CONCATENATED GZIP CATALOGUE
###############################################################################

if (file.exists(final_catalogue_file)) {
  unlink(final_catalogue_file)
}

cat("\nCombining harmonized chunks...\n")

for (i in seq_len(nrow(valid_chunk_table))) {

  source_chunk <- valid_chunk_table$chunk_file[i]

  if (i == 1L) {

    file.copy(
      source_chunk,
      final_catalogue_file,
      overwrite = TRUE
    )

  } else {

    chunk_dt <- fread(
      source_chunk,
      sep = "\t",
      header = TRUE,
      showProgress = FALSE
    )

    temporary_member <- tempfile(
      pattern = "module11_3B_member_",
      fileext = ".tsv.gz"
    )

    fwrite(
      chunk_dt,
      temporary_member,
      sep = "\t",
      quote = FALSE,
      na = "NA",
      col.names = FALSE,
      compress = "gzip"
    )

    append_binary_file(
      temporary_member,
      final_catalogue_file
    )

    unlink(temporary_member)

    rm(chunk_dt)
    invisible(gc())
  }

  if (i %% 25L == 0L || i == nrow(valid_chunk_table)) {

    cat(
      "  Combined:",
      i,
      "/",
      nrow(valid_chunk_table),
      "chunks\n"
    )
  }
}


###############################################################################
# 15. FINAL BUILD STATUS
###############################################################################

module_end_time <- Sys.time()

runtime_seconds <- as.numeric(
  difftime(
    module_end_time,
    module_start_time,
    units = "secs"
  )
)

processed_file_count <- nrow(file_summary)

pass_file_count <- if (nrow(file_summary) > 0L) {
  sum(file_summary$file_status == "PASS")
} else {
  0L
}

flagged_file_count <- if (nrow(file_summary) > 0L) {
  sum(
    file_summary$file_status ==
      "PASS_WITH_ROW_QC_FLAGS"
  )
} else {
  0L
}

failed_file_count <- nrow(failed_files)

total_harmonized_rows <- if (nrow(file_summary) > 0L) {
  sum(file_summary$source_rows, na.rm = TRUE)
} else {
  0
}

all_expected_files_processed <- (
  processed_file_count == nrow(work_queue)
)

all_chunks_available <- (
  nrow(valid_chunk_table) ==
    pass_file_count + flagged_file_count
)

final_catalogue_exists <- (
  file.exists(final_catalogue_file) &&
    file.info(final_catalogue_file)$size > 0
)

final_status <- if (
  failed_file_count == 0L &&
    all_expected_files_processed &&
    all_chunks_available &&
    final_catalogue_exists
) {

  "PASS"

} else if (
  final_catalogue_exists &&
    nrow(valid_chunk_table) > 0L
) {

  "PARTIAL"

} else {

  "FAIL"
}

build_status <- data.table(
  module = "11.3B.3",
  module_name = "harmonized apaQTL catalogue builder",
  start_time = format(
    module_start_time,
    "%Y-%m-%d %H:%M:%S"
  ),
  end_time = format(
    module_end_time,
    "%Y-%m-%d %H:%M:%S"
  ),
  runtime_seconds = round(
    runtime_seconds,
    3
  ),
  expected_files = nrow(work_queue),
  processed_files = processed_file_count,
  pass_files = pass_file_count,
  pass_with_flags_files = flagged_file_count,
  failed_files = failed_file_count,
  valid_chunks = nrow(valid_chunk_table),
  total_harmonized_rows = total_harmonized_rows,
  final_catalogue_exists = final_catalogue_exists,
  final_catalogue_size_bytes = ifelse(
    final_catalogue_exists,
    as.numeric(
      file.info(final_catalogue_file)$size
    ),
    NA_real_
  ),
  final_status = final_status
)

fwrite(
  build_status,
  build_status_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


###############################################################################
# 16. SESSION INFORMATION
###############################################################################

writeLines(
  c(
    "PigmentationAtlas Module 11.3B.3",
    paste("Start:", format(module_start_time)),
    paste("End:", format(module_end_time)),
    paste("Final status:", final_status),
    "",
    capture.output(sessionInfo())
  ),
  con = session_info_output
)


###############################################################################
# 17. FINAL REPORT
###############################################################################

cat("\n")
cat("============================================================\n")
cat("MODULE 11.3B.3 â€“ BUILD SUMMARY\n")
cat("============================================================\n\n")

print(build_status)

if (nrow(file_summary) > 0L) {

  cat("\nFile-status distribution:\n")

  print(
    file_summary[
      ,
      .N,
      by = file_status
    ][
      order(-N)
    ]
  )

  cat("\nRows by tissue:\n")

  print(
    file_summary[
      ,
      .(
        files = .N,
        harmonized_rows = sum(
          source_rows,
          na.rm = TRUE
        ),
        unique_clumps = uniqueN(clump_id)
      ),
      by = tissue
    ][
      order(tissue)
    ]
  )
}

cat("\nOutputs:\n")
cat("  Final catalogue:\n    ", final_catalogue_file, "\n")
cat("  File summary:\n    ", file_summary_output, "\n")
cat("  Failed files:\n    ", failed_files_output, "\n")
cat("  Column mapping:\n    ", column_mapping_output, "\n")
cat("  Progress checkpoint:\n    ", progress_output, "\n")
cat("  Build status:\n    ", build_status_output, "\n")
cat("  Log:\n    ", log_file, "\n\n")

cat("============================================================\n")

if (identical(final_status, "PASS")) {

  cat("FINAL QC RESULT: PASS\n")
  cat("All ", nrow(work_queue), " apaQTL files were processed successfully.\n", sep = "")
  cat("The harmonized apaQTL catalogue was created.\n")
  cat("Module 11.3B.4 final QC may be started.\n")

} else if (identical(final_status, "PARTIAL")) {

  cat("FINAL QC RESULT: PARTIAL\n")
  cat("A partial catalogue was created.\n")
  cat("Review module11_3B_failed_files.tsv and rerun this script.\n")

} else {

  cat("FINAL QC RESULT: FAIL\n")
  cat("The harmonized catalogue could not be completed.\n")
}

cat("============================================================\n\n")


###############################################################################
# 18. CLOSE LOGGING
###############################################################################

while (sink.number(type = "output") > 0L) {
  sink(type = "output")
}

if (isOpen(log_connection)) {
  close(log_connection)
}

logging_active <- FALSE


###############################################################################
# 19. STOP ONLY ON COMPLETE FAILURE
###############################################################################

if (identical(final_status, "FAIL")) {

  stop(
    paste0(
      "Module 11.3B.3 failed. Review:\n",
      failed_files_output
    ),
    call. = FALSE
  )
}

invisible(
  list(
    build_status = build_status,
    file_summary = file_summary,
    failed_files = failed_files,
    progress = progress_table
  )
)
