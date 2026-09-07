###############################################################################
# PigmentationAtlas
# Module 11.3B.1
#
# INITIALIZATION OF THE HARMONIZED apaQTL CATALOGUE
#
# Purpose:
#   1. Read the canonical CLUMP × tissue region map from Module 11.3A v4
#   2. Validate the region map
#   3. Select successfully extracted apaQTL files
#   4. Resolve and validate file paths
#   5. Create the Module 11.3B work queue
#   6. Write initialization QC reports
#
# Important:
#   - Manifest genomic coordinates are NOT used.
#   - Canonical coordinates come from the validated Module 11.3A v4 region map,
#    which was generated from the Module 11.2 publication-locus manifest.
#   - This module does NOT read the contents of the apaQTL files.
#
# Expected input:
#   results/module11/harmonization/region_mapping_v4/
#   module11_3A_v4_direct_clump_region_map.tsv
#
# Main outputs:
#   results/module11/harmonization/catalogue/
#   module11_3B_work_queue.tsv
#   module11_3B_initialization_qc.tsv
#   module11_3B_missing_files.tsv
#   module11_3B_non_extracted_records.tsv
#   module11_3B_initialization_status.tsv
###############################################################################


###############################################################################
# 1. CLEAN SESSION
###############################################################################

rm(list = ls())
invisible(gc())

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  width = 200
)

module_start_time <- Sys.time()

cat("\n")
cat("============================================================\n")
cat("PigmentationAtlas\n")
cat("MODULE 11.3B.1 – INITIALIZATION\n")
cat("Started:", format(module_start_time), "\n")
cat("============================================================\n\n")


###############################################################################
# 2. REQUIRED PACKAGE
###############################################################################

required_packages <- c("data.table")

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    paste0(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", "),
      "\nInstall with:\ninstall.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages(
  library(data.table)
)

cat("Loaded data.table version:",
    as.character(packageVersion("data.table")), "\n\n")


###############################################################################
# 3. PROJECT PATHS
###############################################################################

project_root <- normalizePath(
  "C:/Users/User/Desktop/PigmentationAtlas",
  winslash = "/",
  mustWork = TRUE
)

setwd(project_root)

region_map_file <- file.path(
  project_root,
  "results",
  "module11",
  "harmonization",
  "region_mapping_v4",
  "module11_3A_v4_direct_clump_region_map.tsv"
)

catalogue_dir <- file.path(
  project_root,
  "results",
  "module11",
  "harmonization",
  "catalogue"
)

qc_dir <- file.path(catalogue_dir, "qc")
log_dir <- file.path(catalogue_dir, "logs")
checkpoint_dir <- file.path(catalogue_dir, "checkpoints")
chunk_dir <- file.path(catalogue_dir, "chunks")

directories_to_create <- c(
  catalogue_dir,
  qc_dir,
  log_dir,
  checkpoint_dir,
  chunk_dir
)

for (directory in directories_to_create) {
  if (!dir.exists(directory)) {
    dir.create(
      directory,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
}

cat("Project root:\n", project_root, "\n\n")
cat("Region map input:\n", region_map_file, "\n\n")
cat("Catalogue directory:\n", catalogue_dir, "\n\n")


###############################################################################
# 4. OUTPUT FILE PATHS
###############################################################################

work_queue_file <- file.path(
  catalogue_dir,
  "module11_3B_work_queue.tsv"
)

initialization_qc_file <- file.path(
  qc_dir,
  "module11_3B_initialization_qc.tsv"
)

missing_files_file <- file.path(
  qc_dir,
  "module11_3B_missing_files.tsv"
)

non_extracted_file <- file.path(
  qc_dir,
  "module11_3B_non_extracted_records.tsv"
)

duplicate_pairs_file <- file.path(
  qc_dir,
  "module11_3B_duplicate_clump_tissue_pairs.tsv"
)

invalid_coordinates_file <- file.path(
  qc_dir,
  "module11_3B_invalid_canonical_coordinates.tsv"
)

initialization_status_file <- file.path(
  catalogue_dir,
  "module11_3B_initialization_status.tsv"
)

session_info_file <- file.path(
  log_dir,
  "module11_3B_initialization_session_info.txt"
)

log_file <- file.path(
  log_dir,
  paste0(
    "module11_3B_initialization_",
    format(module_start_time, "%Y%m%d_%H%M%S"),
    ".log"
  )
)


###############################################################################
# 5. HELPER FUNCTIONS
###############################################################################

normalize_column_names <- function(x) {

  x <- trimws(x)
  x <- tolower(x)

  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x <- gsub("_+", "_", x)

  x
}


find_column <- function(column_names,
                        candidates,
                        required = TRUE,
                        field_label = NULL) {

  hits <- candidates[candidates %in% column_names]

  if (length(hits) == 0L) {

    if (isTRUE(required)) {

      if (is.null(field_label)) {
        field_label <- paste(candidates, collapse = " / ")
      }

      stop(
        paste0(
          "Required column not found: ",
          field_label,
          "\nAvailable columns:\n",
          paste(column_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    return(NA_character_)
  }

  hits[1L]
}


normalize_clump_id <- function(x) {

  x <- trimws(as.character(x))
  x <- toupper(x)

  numeric_part <- suppressWarnings(
    as.integer(gsub("[^0-9]", "", x))
  )

  normalized <- ifelse(
    is.na(numeric_part),
    x,
    sprintf("CLUMP_%04d", numeric_part)
  )

  normalized
}


normalize_tissue <- function(x) {

  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_

  x
}


normalize_status <- function(x) {

  x <- trimws(as.character(x))
  x <- toupper(x)

  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x <- gsub("_+", "_", x)

  x
}


is_absolute_path <- function(x) {

  grepl(
    pattern = "^[A-Za-z]:[/\\\\]|^/|^\\\\\\\\",
    x = x
  )
}


resolve_output_path <- function(path_value,
                                project_root) {

  path_value <- trimws(as.character(path_value))

  if (is.na(path_value) || path_value == "") {
    return(NA_character_)
  }

  path_value <- gsub("\\\\", "/", path_value)

  if (is_absolute_path(path_value)) {
    return(path_value)
  }

  file.path(project_root, path_value)
}


safe_file_size <- function(path) {

  if (is.na(path) || !file.exists(path)) {
    return(NA_real_)
  }

  as.numeric(file.info(path)$size)
}


write_empty_or_table <- function(x,
                                 output_file,
                                 template_columns = NULL) {

  if (nrow(x) == 0L && !is.null(template_columns)) {

    empty_table <- as.data.table(
      setNames(
        replicate(
          length(template_columns),
          character(0),
          simplify = FALSE
        ),
        template_columns
      )
    )

    fwrite(
      empty_table,
      output_file,
      sep = "\t",
      quote = FALSE,
      na = "NA"
    )

  } else {

    fwrite(
      x,
      output_file,
      sep = "\t",
      quote = FALSE,
      na = "NA"
    )
  }
}


qc_row <- function(metric,
                   value,
                   expected = NA_character_,
                   result = NA_character_,
                   interpretation = NA_character_) {

  data.table(
    metric = as.character(metric),
    value = as.character(value),
    expected = as.character(expected),
    result = as.character(result),
    interpretation = as.character(interpretation)
  )
}


###############################################################################
# 6. START LOGGING
###############################################################################

log_connection <- file(log_file, open = "wt")

sink(log_connection, type = "output", split = TRUE)
sink(log_connection, type = "message")

logging_active <- TRUE

on.exit({

  if (exists("logging_active", inherits = FALSE) &&
      isTRUE(logging_active)) {

    try(sink(type = "message"), silent = TRUE)
    try(sink(type = "output"), silent = TRUE)
    try(close(log_connection), silent = TRUE)
  }

}, add = TRUE)

cat("Log file:\n", log_file, "\n\n")


###############################################################################
# 7. CHECK INPUT FILE
###############################################################################

if (!file.exists(region_map_file)) {

  stop(
    paste0(
      "Module 11.3A v4 region map was not found:\n",
      region_map_file
    ),
    call. = FALSE
  )
}

cat("Region map exists: YES\n")
cat("Region map size:",
    format(file.info(region_map_file)$size, big.mark = ","),
    "bytes\n\n")


###############################################################################
# 8. READ REGION MAP
###############################################################################

region_map <- fread(
  region_map_file,
  sep = "\t",
  header = TRUE,
  na.strings = c("", "NA", "NaN", "NULL"),
  showProgress = TRUE
)

if (nrow(region_map) == 0L) {
  stop("The Module 11.3A region map is empty.", call. = FALSE)
}

original_column_names <- names(region_map)
normalized_names <- normalize_column_names(original_column_names)

if (anyDuplicated(normalized_names)) {

  duplicated_normalized <- unique(
    normalized_names[
      duplicated(normalized_names) |
        duplicated(normalized_names, fromLast = TRUE)
    ]
  )

  stop(
    paste0(
      "Column-name normalization created duplicate names:\n",
      paste(duplicated_normalized, collapse = ", ")
    ),
    call. = FALSE
  )
}

setnames(
  region_map,
  old = original_column_names,
  new = normalized_names
)

cat("Region map loaded.\n")
cat("Rows:", format(nrow(region_map), big.mark = ","), "\n")
cat("Columns:", ncol(region_map), "\n\n")

cat("Normalized region-map columns:\n")
print(names(region_map))
cat("\n")


###############################################################################
# 9. DETECT REQUIRED COLUMNS
###############################################################################

clump_col <- find_column(
  names(region_map),
  candidates = c(
    "clump_id",
    "clump",
    "locus_id",
    "locus"
  ),
  field_label = "CLUMP ID"
)

tissue_col <- find_column(
  names(region_map),
  candidates = c(
    "tissue",
    "tissue_name",
    "dataset",
    "condition"
  ),
  field_label = "tissue"
)

status_col <- find_column(
  names(region_map),
  candidates = c(
    "extraction_status",
    "status",
    "apaqtl_extraction_status"
  ),
  field_label = "extraction status"
)

output_file_col <- find_column(
  names(region_map),
  candidates = c(
    "output_file",
    "output_path",
    "extracted_file",
    "apaqtl_file",
    "file"
  ),
  field_label = "apaQTL output file"
)

chromosome_col <- find_column(
  names(region_map),
  candidates = c(
    "canonical_chromosome",
    "master_chromosome",
    "chromosome"
  ),
  field_label = "canonical chromosome"
)

region_start_col <- find_column(
  names(region_map),
  candidates = c(
    "canonical_region_start",
    "master_region_start",
    "region_start"
  ),
  field_label = "canonical region start"
)

region_end_col <- find_column(
  names(region_map),
  candidates = c(
    "canonical_region_end",
    "master_region_end",
    "region_end"
  ),
  field_label = "canonical region end"
)

mapping_status_col <- find_column(
  names(region_map),
  candidates = c(
    "mapping_status",
    "map_status"
  ),
  required = FALSE
)

cat("Detected columns:\n")
cat("  CLUMP ID:              ", clump_col, "\n")
cat("  Tissue:                ", tissue_col, "\n")
cat("  Extraction status:     ", status_col, "\n")
cat("  Output file:           ", output_file_col, "\n")
cat("  Canonical chromosome:  ", chromosome_col, "\n")
cat("  Canonical region start:", region_start_col, "\n")
cat("  Canonical region end:  ", region_end_col, "\n")
cat("  Mapping status:        ",
    ifelse(is.na(mapping_status_col), "<not present>", mapping_status_col),
    "\n\n")


###############################################################################
# 10. CONSTRUCT STANDARDIZED REGION MAP
###############################################################################

region_map_standard <- data.table(
  clump_id = normalize_clump_id(region_map[[clump_col]]),
  tissue = normalize_tissue(region_map[[tissue_col]]),
  extraction_status = normalize_status(region_map[[status_col]]),
  output_file_original = trimws(
    as.character(region_map[[output_file_col]])
  ),
  canonical_chromosome = gsub(
    "^CHR",
    "",
    toupper(trimws(as.character(region_map[[chromosome_col]])))
  ),
  canonical_region_start = suppressWarnings(
    as.integer(region_map[[region_start_col]])
  ),
  canonical_region_end = suppressWarnings(
    as.integer(region_map[[region_end_col]])
  )
)

if (!is.na(mapping_status_col)) {

  region_map_standard[
    ,
    mapping_status := normalize_status(
      region_map[[mapping_status_col]]
    )
  ]

} else {

  region_map_standard[
    ,
    mapping_status := "DIRECTLY_MAPPED"
  ]
}

region_map_standard[
  ,
  source_region_map_row := .I
]

region_map_standard[
  canonical_chromosome %chin% c("23", "X"),
  canonical_chromosome := "X"
]

region_map_standard[
  canonical_chromosome %chin% c("24", "Y"),
  canonical_chromosome := "Y"
]

region_map_standard[
  canonical_chromosome %chin% c("25", "M", "MT", "MITO"),
  canonical_chromosome := "MT"
]

region_map_standard[
  ,
  output_file_resolved := vapply(
    output_file_original,
    resolve_output_path,
    character(1),
    project_root = project_root
  )
]

region_map_standard[
  ,
  output_file_exists := !is.na(output_file_resolved) &
    file.exists(output_file_resolved)
]

region_map_standard[
  ,
  output_file_size_bytes := vapply(
    output_file_resolved,
    safe_file_size,
    numeric(1)
  )
]

region_map_standard[
  ,
  output_file_nonempty := output_file_exists &
    !is.na(output_file_size_bytes) &
    output_file_size_bytes > 0
]


###############################################################################
# 11. BASIC METADATA QC
###############################################################################

region_map_standard[
  ,
  clump_tissue_key := paste(
    clump_id,
    tissue,
    sep = "||"
  )
]

duplicate_pairs <- region_map_standard[
  !is.na(clump_id) &
    !is.na(tissue) &
    duplicated(clump_tissue_key) |
    duplicated(clump_tissue_key, fromLast = TRUE)
][
  order(clump_id, tissue, source_region_map_row)
]

invalid_coordinates <- region_map_standard[
  is.na(canonical_chromosome) |
    canonical_chromosome == "" |
    is.na(canonical_region_start) |
    is.na(canonical_region_end) |
    canonical_region_start < 1L |
    canonical_region_end < canonical_region_start
][
  order(clump_id, tissue)
]

missing_clump_ids <- region_map_standard[
  is.na(clump_id) | clump_id == ""
]

missing_tissues <- region_map_standard[
  is.na(tissue) | tissue == ""
]

cat("Standardized region-map QC:\n")
cat("  Missing CLUMP IDs:", nrow(missing_clump_ids), "\n")
cat("  Missing tissues:", nrow(missing_tissues), "\n")
cat("  Duplicate CLUMP × tissue pairs:", nrow(duplicate_pairs), "\n")
cat("  Invalid canonical coordinates:", nrow(invalid_coordinates), "\n\n")


###############################################################################
# 12. DEFINE EXTRACTION STATUS GROUPS
###############################################################################

cat("Extraction-status distribution:\n")
print(
  region_map_standard[
    ,
    .N,
    by = extraction_status
  ][
    order(-N, extraction_status)
  ]
)
cat("\n")

successful_statuses <- c(
  "EXTRACTED",
  "SUCCESS",
  "SUCCESSFULLY_EXTRACTED"
)

work_queue <- region_map_standard[
  extraction_status %chin% successful_statuses
]

non_extracted_records <- region_map_standard[
  !extraction_status %chin% successful_statuses
]

cat("Rows selected for the work queue:",
    format(nrow(work_queue), big.mark = ","), "\n")

cat("Rows excluded from the work queue:",
    format(nrow(non_extracted_records), big.mark = ","), "\n\n")

if (nrow(work_queue) == 0L) {

  stop(
    paste0(
      "No successfully extracted records were found.\n",
      "Observed statuses: ",
      paste(
        sort(unique(region_map_standard$extraction_status)),
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


###############################################################################
# 13. CHECK WORK-QUEUE FILES
###############################################################################

missing_files <- work_queue[
  is.na(output_file_resolved) |
    output_file_resolved == "" |
    !output_file_exists
][
  order(clump_id, tissue)
]

empty_files <- work_queue[
  output_file_exists == TRUE &
    output_file_nonempty == FALSE
][
  order(clump_id, tissue)
]

valid_work_queue <- work_queue[
  output_file_exists == TRUE &
    output_file_nonempty == TRUE
]

setorder(
  valid_work_queue,
  canonical_chromosome,
  canonical_region_start,
  clump_id,
  tissue
)

valid_work_queue[
  ,
  queue_id := sprintf(
    "APAQTL_%04d",
    seq_len(.N)
  )
]

setcolorder(
  valid_work_queue,
  c(
    "queue_id",
    "clump_id",
    "tissue",
    "canonical_chromosome",
    "canonical_region_start",
    "canonical_region_end",
    "mapping_status",
    "extraction_status",
    "output_file_original",
    "output_file_resolved",
    "output_file_exists",
    "output_file_size_bytes",
    "output_file_nonempty",
    "source_region_map_row",
    "clump_tissue_key"
  )
)

cat("Work-queue file QC:\n")
cat("  Expected extracted records:",
    format(nrow(work_queue), big.mark = ","), "\n")
cat("  Existing non-empty files:",
    format(nrow(valid_work_queue), big.mark = ","), "\n")
cat("  Missing files:",
    format(nrow(missing_files), big.mark = ","), "\n")
cat("  Empty files:",
    format(nrow(empty_files), big.mark = ","), "\n\n")


###############################################################################
# 14. EXPECTED COUNTS FROM MODULE 11.3A
###############################################################################

expected_region_map_rows <- 530L
expected_unique_clumps <- 265L
expected_extracted_files <- 516L
expected_non_extracted_records <- 14L

observed_region_map_rows <- nrow(region_map_standard)

observed_unique_clumps <- uniqueN(
  region_map_standard$clump_id,
  na.rm = TRUE
)

observed_unique_clump_tissue_pairs <- uniqueN(
  region_map_standard$clump_tissue_key,
  na.rm = TRUE
)

observed_extracted_records <- nrow(work_queue)
observed_valid_work_queue <- nrow(valid_work_queue)
observed_non_extracted_records <- nrow(non_extracted_records)


###############################################################################
# 15. INITIALIZATION QC TABLE
###############################################################################

initialization_qc <- rbindlist(
  list(

    qc_row(
      metric = "region_map_file_exists",
      value = file.exists(region_map_file),
      expected = "TRUE",
      result = ifelse(file.exists(region_map_file), "PASS", "FAIL"),
      interpretation = "Module 11.3A v4 region map is available."
    ),

    qc_row(
      metric = "region_map_rows",
      value = observed_region_map_rows,
      expected = expected_region_map_rows,
      result = ifelse(
        observed_region_map_rows == expected_region_map_rows,
        "PASS",
        "REVIEW"
      ),
      interpretation = "Total CLUMP × tissue records in the canonical map."
    ),

    qc_row(
      metric = "unique_clumps",
      value = observed_unique_clumps,
      expected = expected_unique_clumps,
      result = ifelse(
        observed_unique_clumps == expected_unique_clumps,
        "PASS",
        "REVIEW"
      ),
      interpretation = "Unique apaQTL CLUMPs represented in the region map."
    ),

    qc_row(
      metric = "unique_clump_tissue_pairs",
      value = observed_unique_clump_tissue_pairs,
      expected = expected_region_map_rows,
      result = ifelse(
        observed_unique_clump_tissue_pairs == expected_region_map_rows,
        "PASS",
        "FAIL"
      ),
      interpretation = "Each CLUMP × tissue pair should occur once."
    ),

    qc_row(
      metric = "extracted_records",
      value = observed_extracted_records,
      expected = expected_extracted_files,
      result = ifelse(
        observed_extracted_records == expected_extracted_files,
        "PASS",
        "REVIEW"
      ),
      interpretation = "Records expected to have an extracted apaQTL file."
    ),

    qc_row(
      metric = "valid_work_queue_rows",
      value = observed_valid_work_queue,
      expected = expected_extracted_files,
      result = ifelse(
        observed_valid_work_queue == expected_extracted_files,
        "PASS",
        "FAIL"
      ),
      interpretation = "Extracted records with an existing, non-empty file."
    ),

    qc_row(
      metric = "non_extracted_records",
      value = observed_non_extracted_records,
      expected = expected_non_extracted_records,
      result = ifelse(
        observed_non_extracted_records == expected_non_extracted_records,
        "PASS",
        "REVIEW"
      ),
      interpretation = "Records without extracted apaQTL rows."
    ),

    qc_row(
      metric = "missing_output_files",
      value = nrow(missing_files),
      expected = 0,
      result = ifelse(
        nrow(missing_files) == 0L,
        "PASS",
        "FAIL"
      ),
      interpretation = "All work-queue output files must exist."
    ),

    qc_row(
      metric = "empty_output_files",
      value = nrow(empty_files),
      expected = 0,
      result = ifelse(
        nrow(empty_files) == 0L,
        "PASS",
        "FAIL"
      ),
      interpretation = "All work-queue output files must be non-empty."
    ),

    qc_row(
      metric = "duplicate_clump_tissue_rows",
      value = nrow(duplicate_pairs),
      expected = 0,
      result = ifelse(
        nrow(duplicate_pairs) == 0L,
        "PASS",
        "FAIL"
      ),
      interpretation = "Duplicate CLUMP × tissue records are not permitted."
    ),

    qc_row(
      metric = "invalid_canonical_coordinate_rows",
      value = nrow(invalid_coordinates),
      expected = 0,
      result = ifelse(
        nrow(invalid_coordinates) == 0L,
        "PASS",
        "FAIL"
      ),
      interpretation = paste(
        "Canonical coordinates must be complete and region_end",
        "must be greater than or equal to region_start."
      )
    ),

    qc_row(
      metric = "missing_clump_id_rows",
      value = nrow(missing_clump_ids),
      expected = 0,
      result = ifelse(
        nrow(missing_clump_ids) == 0L,
        "PASS",
        "FAIL"
      ),
      interpretation = "Every region-map record must have a CLUMP ID."
    ),

    qc_row(
      metric = "missing_tissue_rows",
      value = nrow(missing_tissues),
      expected = 0,
      result = ifelse(
        nrow(missing_tissues) == 0L,
        "PASS",
        "FAIL"
      ),
      interpretation = "Every region-map record must have a tissue."
    )
  ),
  use.names = TRUE,
  fill = TRUE
)


###############################################################################
# 16. DETERMINE FINAL INITIALIZATION STATUS
###############################################################################

hard_failure_metrics <- initialization_qc[
  result == "FAIL",
  metric
]

review_metrics <- initialization_qc[
  result == "REVIEW",
  metric
]

if (length(hard_failure_metrics) > 0L) {

  final_status <- "FAIL"

} else if (length(review_metrics) > 0L) {

  final_status <- "PASS_WITH_REVIEW"

} else {

  final_status <- "PASS"
}

module_end_time <- Sys.time()
runtime_seconds <- as.numeric(
  difftime(
    module_end_time,
    module_start_time,
    units = "secs"
  )
)

initialization_status <- data.table(
  module = "11.3B.1",
  module_name = "apaQTL catalogue initialization",
  pipeline_version = "1.0.0",
  start_time = format(
    module_start_time,
    "%Y-%m-%d %H:%M:%S"
  ),
  end_time = format(
    module_end_time,
    "%Y-%m-%d %H:%M:%S"
  ),
  runtime_seconds = round(runtime_seconds, 3),
  final_status = final_status,
  hard_failure_count = length(hard_failure_metrics),
  review_count = length(review_metrics),
  region_map_rows = observed_region_map_rows,
  unique_clumps = observed_unique_clumps,
  extracted_records = observed_extracted_records,
  valid_work_queue_rows = observed_valid_work_queue,
  missing_files = nrow(missing_files),
  empty_files = nrow(empty_files)
)


###############################################################################
# 17. WRITE OUTPUT TABLES
###############################################################################

fwrite(
  valid_work_queue,
  work_queue_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  initialization_qc,
  initialization_qc_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

write_empty_or_table(
  missing_files,
  missing_files_file,
  template_columns = names(work_queue)
)

write_empty_or_table(
  non_extracted_records,
  non_extracted_file,
  template_columns = names(region_map_standard)
)

write_empty_or_table(
  duplicate_pairs,
  duplicate_pairs_file,
  template_columns = names(region_map_standard)
)

write_empty_or_table(
  invalid_coordinates,
  invalid_coordinates_file,
  template_columns = names(region_map_standard)
)

fwrite(
  initialization_status,
  initialization_status_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


###############################################################################
# 18. SAVE SESSION INFORMATION
###############################################################################

session_lines <- capture.output(
  sessionInfo()
)

writeLines(
  c(
    "PigmentationAtlas Module 11.3B.1",
    paste("Start:", format(module_start_time)),
    paste("End:", format(module_end_time)),
    paste("Final status:", final_status),
    "",
    session_lines
  ),
  con = session_info_file
)


###############################################################################
# 19. FINAL REPORT
###############################################################################

cat("\n")
cat("============================================================\n")
cat("MODULE 11.3B.1 – INITIALIZATION SUMMARY\n")
cat("============================================================\n\n")

print(initialization_qc)

cat("\nCanonical chromosome distribution in valid work queue:\n")

chromosome_order <- c(
  as.character(1:22),
  "X",
  "Y",
  "MT"
)

chromosome_distribution <- valid_work_queue[
  ,
  .(
    files = .N,
    unique_clumps = uniqueN(clump_id)
  ),
  by = canonical_chromosome
]

chromosome_distribution[
  ,
  chromosome_sort := match(
    canonical_chromosome,
    chromosome_order
  )
]

setorder(
  chromosome_distribution,
  chromosome_sort,
  canonical_chromosome
)

chromosome_distribution[
  ,
  chromosome_sort := NULL
]

print(chromosome_distribution)

cat("\nTissue distribution in valid work queue:\n")

print(
  valid_work_queue[
    ,
    .(
      files = .N,
      unique_clumps = uniqueN(clump_id)
    ),
    by = tissue
  ][
    order(tissue)
  ]
)

cat("\n")
cat("Final initialization status:", final_status, "\n")
cat("Runtime:", round(runtime_seconds, 2), "seconds\n\n")

cat("Outputs written:\n")
cat("  Work queue:\n    ", work_queue_file, "\n")
cat("  Initialization QC:\n    ", initialization_qc_file, "\n")
cat("  Missing files:\n    ", missing_files_file, "\n")
cat("  Non-extracted records:\n    ", non_extracted_file, "\n")
cat("  Duplicate pairs:\n    ", duplicate_pairs_file, "\n")
cat("  Invalid coordinates:\n    ", invalid_coordinates_file, "\n")
cat("  Initialization status:\n    ", initialization_status_file, "\n")
cat("  Session information:\n    ", session_info_file, "\n")
cat("  Full log:\n    ", log_file, "\n\n")

cat("============================================================\n")

if (identical(final_status, "PASS")) {

  cat("FINAL QC RESULT: PASS\n")
  cat("Module 11.3B.2 may be started.\n")

} else if (identical(final_status, "PASS_WITH_REVIEW")) {

  cat("FINAL QC RESULT: PASS WITH REVIEW\n")
  cat("No structural failure was detected, but expected counts changed.\n")
  cat("Review the QC table before starting Module 11.3B.2.\n")

} else {

  cat("FINAL QC RESULT: FAIL\n")
  cat("Module 11.3B.2 must NOT be started.\n")
  cat("Failed metrics:\n")
  cat(
    paste0(
      "  - ",
      hard_failure_metrics,
      collapse = "\n"
    ),
    "\n"
  )
}

cat("============================================================\n\n")


###############################################################################
# 20. CLOSE LOGGING CLEANLY
###############################################################################

if (sink.number(type = "message") > 0)
  sink(type = "message")

if (sink.number(type = "output") > 0)
  sink(type = "output")

try(close(log_connection), silent = TRUE)

logging_active <- FALSE


###############################################################################
# 21. STOP ON HARD FAILURE
###############################################################################

if (identical(final_status, "FAIL")) {

  stop(
    paste0(
      "Module 11.3B.1 initialization failed. ",
      "See QC and log outputs in:\n",
      catalogue_dir
    ),
    call. = FALSE
  )
}

invisible(
  list(
    region_map = region_map_standard,
    work_queue = valid_work_queue,
    initialization_qc = initialization_qc,
    initialization_status = initialization_status
  )
)