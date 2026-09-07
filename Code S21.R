###############################################################################
# PigmentationAtlas
# Module 11.3A v4
#
# ROBUST DIRECT CLUMP × apaQTL REGION MAPPING
#
# Purpose:
#   1. Read canonical CLUMP metadata
#   2. Discover regional apaQTL files in both skin tissues
#   3. Parse CLUMP ID, chromosome, start and end from filenames
#   4. Match files to canonical CLUMPs using:
#        CLUMP ID + chromosome + coordinate overlap
#   5. Reject chromosome-mismatched or ambiguous assignments
#   6. Create the complete CLUMP × tissue region map
#
# Main output:
#   results/module11/harmonization/region_mapping_v4/
#   module11_3A_v4_direct_clump_region_map.tsv
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
cat("MODULE 11.3A v4 – ROBUST DIRECT CLUMP REGION MAPPING\n")
cat("Started:", format(module_start_time), "\n")
cat("============================================================\n\n")


###############################################################################
# 2. REQUIRED PACKAGE
###############################################################################

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop(
    "The data.table package is required.\n",
    "Install it with:\ninstall.packages('data.table')",
    call. = FALSE
  )
}

suppressPackageStartupMessages(
  library(data.table)
)

cat(
  "Loaded data.table version:",
  as.character(packageVersion("data.table")),
  "\n\n"
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

apaqtl_root <- file.path(
  project_root,
  "results",
  "module11",
  "eqtl_regions"
)

output_dir <- file.path(
  project_root,
  "results",
  "module11",
  "harmonization",
  "region_mapping_v4"
)

qc_dir <- file.path(output_dir, "qc")
log_dir <- file.path(output_dir, "logs")

for (directory in c(output_dir, qc_dir, log_dir)) {
  if (!dir.exists(directory)) {
    dir.create(
      directory,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
}

output_map_file <- file.path(
  output_dir,
  "module11_3A_v4_direct_clump_region_map.tsv"
)

file_inventory_output <- file.path(
  qc_dir,
  "module11_3A_v4_apaqtl_file_inventory.tsv"
)

mapping_qc_output <- file.path(
  qc_dir,
  "module11_3A_v4_mapping_qc.tsv"
)

unmatched_files_output <- file.path(
  qc_dir,
  "module11_3A_v4_unmatched_apaqtl_files.tsv"
)

ambiguous_files_output <- file.path(
  qc_dir,
  "module11_3A_v4_ambiguous_apaqtl_files.tsv"
)

chromosome_mismatch_output <- file.path(
  qc_dir,
  "module11_3A_v4_chromosome_mismatches.tsv"
)

duplicate_assignments_output <- file.path(
  qc_dir,
  "module11_3A_v4_duplicate_assignments.tsv"
)

status_output <- file.path(
  output_dir,
  "module11_3A_v4_status.tsv"
)

session_info_output <- file.path(
  log_dir,
  "module11_3A_v4_session_info.txt"
)

log_file <- file.path(
  log_dir,
  paste0(
    "module11_3A_v4_",
    format(module_start_time, "%Y%m%d_%H%M%S"),
    ".log"
  )
)


###############################################################################
# 4. HELPER FUNCTIONS
###############################################################################

normalize_column_names <- function(x) {

  x <- trimws(as.character(x))
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

  ifelse(
    is.na(numeric_part),
    x,
    sprintf("CLUMP_%04d", numeric_part)
  )
}


normalize_chromosome <- function(x) {

  x <- toupper(trimws(as.character(x)))
  x <- sub("^CHR", "", x)

  x[x %in% c("23", "X")] <- "X"
  x[x %in% c("24", "Y")] <- "Y"
  x[x %in% c("25", "M", "MT", "MITO")] <- "MT"
  x[x == ""] <- NA_character_

  x
}


chromosome_sort_value <- function(x) {

  x <- normalize_chromosome(x)

  value <- suppressWarnings(as.integer(x))

  value[x == "X"] <- 23L
  value[x == "Y"] <- 24L
  value[x == "MT"] <- 25L

  value
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
                   observed,
                   expected,
                   status,
                   interpretation) {

  data.table(
    metric = as.character(metric),
    observed = as.character(observed),
    expected = as.character(expected),
    status = as.character(status),
    interpretation = as.character(interpretation)
  )
}


###############################################################################
# 5. LOGGING
###############################################################################

log_connection <- file(log_file, open = "wt")

sink(
  log_connection,
  type = "output",
  split = TRUE
)

logging_active <- TRUE

on.exit({

  if (exists("logging_active", inherits = FALSE) &&
      isTRUE(logging_active)) {

    while (sink.number(type = "output") > 0L) {
      try(sink(type = "output"), silent = TRUE)
    }

    try(close(log_connection), silent = TRUE)
  }

}, add = TRUE)

cat("Project root:\n", project_root, "\n\n")
cat("apaQTL root:\n", apaqtl_root, "\n\n")
cat("Output directory:\n", output_dir, "\n\n")
cat("Log file:\n", log_file, "\n\n")

###############################################################################
# 6. READ MODULE 11.2 MANIFEST
###############################################################################

manifest_file <- file.path(
  project_root,
  "results",
  "module11",
  "eqtl_regions",
  "summary",
  "module11_2_apaqtl_region_manifest.tsv"
)

if (!file.exists(manifest_file)) {

  stop(
    paste(
      "Cannot find:",
      manifest_file
    ),
    call.=FALSE
  )

}

region_manifest <- fread(manifest_file)

cat(
  "Module 11.2 manifest:",
  nrow(region_manifest),
  "rows\n"
)

###############################################################################
# 7. BUILD CANONICAL CLUMPS FROM MODULE 11.2 MANIFEST
###############################################################################

# Normalize manifest column names
original_manifest_names <- names(region_manifest)
normalized_manifest_names <- normalize_column_names(
  original_manifest_names
)

if (anyDuplicated(normalized_manifest_names)) {
  stop(
    "Column-name normalization created duplicate manifest column names.",
    call. = FALSE
  )
}

setnames(
  region_manifest,
  old = original_manifest_names,
  new = normalized_manifest_names
)

cat("Manifest column names:\n")
print(names(region_manifest))
cat("\n")

# Detect required columns
manifest_clump_col <- find_column(
  names(region_manifest),
  candidates = c(
    "clump_id",
    "locus_id",
    "publication_locus_id",
    "clump",
    "locus"
  ),
  field_label = "publication locus / CLUMP ID"
)

manifest_chr_col <- find_column(
  names(region_manifest),
  candidates = c(
    "chromosome",
    "chr",
    "chrom",
    "canonical_chromosome"
  ),
  field_label = "chromosome"
)

manifest_start_col <- find_column(
  names(region_manifest),
  candidates = c(
    "region_start",
    "locus_start",
    "window_start",
    "start",
    "canonical_region_start"
  ),
  field_label = "region start"
)

manifest_end_col <- find_column(
  names(region_manifest),
  candidates = c(
    "region_end",
    "locus_end",
    "window_end",
    "end",
    "canonical_region_end"
  ),
  field_label = "region end"
)

cat("Detected manifest columns:\n")
cat("  CLUMP/locus ID:", manifest_clump_col, "\n")
cat("  Chromosome:    ", manifest_chr_col, "\n")
cat("  Region start:  ", manifest_start_col, "\n")
cat("  Region end:    ", manifest_end_col, "\n\n")

# Build standardized canonical table
canonical_clumps <- region_manifest[
  ,
  .(
    clump_id = normalize_clump_id(
      get(manifest_clump_col)
    ),
    canonical_chromosome = normalize_chromosome(
      get(manifest_chr_col)
    ),
    canonical_region_start = suppressWarnings(
      as.integer(get(manifest_start_col))
    ),
    canonical_region_end = suppressWarnings(
      as.integer(get(manifest_end_col))
    )
  )
]

# One row per publication locus
canonical_clumps <- unique(
  canonical_clumps,
  by = c(
    "clump_id",
    "canonical_chromosome",
    "canonical_region_start",
    "canonical_region_end"
  )
)

canonical_clumps[
  ,
  canonical_key := paste(
    clump_id,
    canonical_chromosome,
    sep = "||"
  )
]

canonical_clumps[
  ,
  source_master_row := seq_len(.N)
]

# Validate rows
invalid_master_rows <- canonical_clumps[
  is.na(clump_id) |
    clump_id == "" |
    is.na(canonical_chromosome) |
    canonical_chromosome == "" |
    is.na(canonical_region_start) |
    is.na(canonical_region_end) |
    canonical_region_start < 1L |
    canonical_region_end < canonical_region_start
]

duplicate_canonical_keys <- canonical_clumps[
  duplicated(canonical_key) |
    duplicated(canonical_key, fromLast = TRUE)
][
  order(
    clump_id,
    canonical_chromosome,
    canonical_region_start
  )
]

cat("Canonical CLUMP QC:\n")
cat("  Rows:", nrow(canonical_clumps), "\n")
cat(
  "  Unique CLUMP IDs:",
  uniqueN(canonical_clumps$clump_id),
  "\n"
)
cat(
  "  Unique CLUMP + chromosome keys:",
  uniqueN(canonical_clumps$canonical_key),
  "\n"
)
cat("  Invalid rows:", nrow(invalid_master_rows), "\n")
cat(
  "  Duplicate CLUMP + chromosome keys:",
  nrow(duplicate_canonical_keys),
  "\n\n"
)

if (nrow(invalid_master_rows) > 0L) {
  print(invalid_master_rows)

  stop(
    "Invalid publication-locus metadata rows were detected.",
    call. = FALSE
  )
}

if (nrow(duplicate_canonical_keys) > 0L) {
  print(duplicate_canonical_keys)

  stop(
    paste0(
      "The same CLUMP ID occurs more than once on the same chromosome.\n",
      "CLUMP ID + chromosome must identify one canonical region."
    ),
    call. = FALSE
  )
}

if (nrow(canonical_clumps) != 265L) {
  stop(
    paste0(
      "Expected 265 canonical publication loci, but found ",
      nrow(canonical_clumps),
      "."
    ),
    call. = FALSE
  )
}

cat(
  "Canonical publication loci successfully created: ",
  nrow(canonical_clumps),
  "\n\n",
  sep = ""
)

###############################################################################
# 8. DEFINE TISSUES
###############################################################################

expected_tissues <- c(
  "Skin_Not_Sun_Exposed_Suprapubic",
  "Skin_Sun_Exposed_Lower_leg"
)

missing_tissue_directories <- expected_tissues[
  !dir.exists(file.path(apaqtl_root, expected_tissues))
]

if (length(missing_tissue_directories) > 0L) {

  stop(
    paste0(
      "Missing expected tissue directories:\n",
      paste(missing_tissue_directories, collapse = "\n")
    ),
    call. = FALSE
  )
}

cat("Expected tissues:\n")
print(expected_tissues)
cat("\n")


###############################################################################
# 9. DISCOVER apaQTL FILES
###############################################################################

apaqtl_files <- unlist(
  lapply(
    expected_tissues,
    function(tissue_name) {

      list.files(
        file.path(apaqtl_root, tissue_name),
        pattern = "\\.apaqtl\\.tsv\\.gz$",
        recursive = FALSE,
        full.names = TRUE,
        ignore.case = TRUE
      )
    }
  ),
  use.names = FALSE
)

apaqtl_files <- normalizePath(
  apaqtl_files,
  winslash = "/",
  mustWork = TRUE
)

apaqtl_files <- unique(apaqtl_files)

cat(
  "Regional apaQTL files discovered:",
  length(apaqtl_files),
  "\n\n"
)

if (length(apaqtl_files) == 0L) {
  stop("No regional apaQTL files were found.", call. = FALSE)
}


###############################################################################
# 10. PARSE FILE METADATA
###############################################################################

file_inventory <- data.table(
  output_file = apaqtl_files
)

file_inventory[
  ,
  filename := basename(output_file)
]

file_inventory[
  ,
  tissue := basename(dirname(output_file))
]

filename_pattern <- paste0(
  "^",
  "(CLUMP_[0-9]+)",
  "_chr",
  "([^_]+)",
  "_([0-9]+)",
  "_([0-9]+)",
  "\\.apaqtl\\.tsv\\.gz$"
)

parsed_values <- tstrsplit(
  file_inventory$filename,
  filename_pattern,
  perl = TRUE
)

regex_capture <- regexec(
  filename_pattern,
  file_inventory$filename,
  perl = TRUE
)

regex_matches <- regmatches(
  file_inventory$filename,
  regex_capture
)

extract_capture <- function(match_vector, position) {

  vapply(
    match_vector,
    function(x) {

      if (length(x) >= position + 1L) {
        x[position + 1L]
      } else {
        NA_character_
      }
    },
    character(1)
  )
}

file_inventory[
  ,
  file_clump_id := normalize_clump_id(
    extract_capture(regex_matches, 1L)
  )
]

file_inventory[
  ,
  file_chromosome := normalize_chromosome(
    extract_capture(regex_matches, 2L)
  )
]

file_inventory[
  ,
  file_region_start := suppressWarnings(
    as.integer(extract_capture(regex_matches, 3L))
  )
]

file_inventory[
  ,
  file_region_end := suppressWarnings(
    as.integer(extract_capture(regex_matches, 4L))
  )
]

file_inventory[
  ,
  filename_parse_status := fifelse(
    !is.na(file_clump_id) &
      !is.na(file_chromosome) &
      !is.na(file_region_start) &
      !is.na(file_region_end) &
      file_region_start >= 1L &
      file_region_end >= file_region_start,
    "PASS",
    "FAIL"
  )
]

file_inventory[
  ,
  output_file_exists := file.exists(output_file)
]

file_inventory[
  ,
  output_file_size_bytes := vapply(
    output_file,
    safe_file_size,
    numeric(1)
  )
]

file_inventory[
  ,
  output_file_nonempty :=
    output_file_exists &
    !is.na(output_file_size_bytes) &
    output_file_size_bytes > 0
]

file_inventory[
  ,
  file_key := paste(
    file_clump_id,
    file_chromosome,
    tissue,
    sep = "||"
  )
]

invalid_filenames <- file_inventory[
  filename_parse_status == "FAIL"
]

duplicate_files <- file_inventory[
  duplicated(file_key) |
    duplicated(file_key, fromLast = TRUE)
][
  order(file_clump_id, file_chromosome, tissue, filename)
]

cat("File inventory QC:\n")
cat("  Parsed successfully:",
    file_inventory[filename_parse_status == "PASS", .N], "\n")
cat("  Filename parse failures:", nrow(invalid_filenames), "\n")
cat("  Duplicate CLUMP + chromosome + tissue keys:",
    nrow(duplicate_files), "\n")
cat("  Missing files:",
    file_inventory[output_file_exists == FALSE, .N], "\n")
cat("  Empty files:",
    file_inventory[output_file_nonempty == FALSE, .N], "\n\n")

if (nrow(invalid_filenames) > 0L) {

  print(invalid_filenames)

  stop(
    "One or more apaQTL filenames could not be parsed.",
    call. = FALSE
  )
}

if (nrow(duplicate_files) > 0L) {

  print(duplicate_files)

  stop(
    paste0(
      "Duplicate regional apaQTL files were found for the same ",
      "CLUMP + chromosome + tissue combination."
    ),
    call. = FALSE
  )
}


###############################################################################
# 11. CREATE CANDIDATE MATCHES
###############################################################################

candidate_matches <- merge(
  file_inventory,
  canonical_clumps,
  by.x = c(
    "file_clump_id",
    "file_chromosome"
  ),
  by.y = c(
    "clump_id",
    "canonical_chromosome"
  ),
  all.x = TRUE,
  allow.cartesian = TRUE,
  sort = FALSE
)

setnames(
  candidate_matches,
  old = c(
    "file_clump_id",
    "file_chromosome"
  ),
  new = c(
    "clump_id",
    "canonical_chromosome"
  )
)

candidate_matches[
  ,
  coordinate_overlap :=
    !is.na(canonical_region_start) &
    !is.na(canonical_region_end) &
    file_region_start <= canonical_region_end &
    file_region_end >= canonical_region_start
]

candidate_matches[
  ,
  overlap_start := pmax(
    file_region_start,
    canonical_region_start,
    na.rm = TRUE
  )
]

candidate_matches[
  ,
  overlap_end := pmin(
    file_region_end,
    canonical_region_end,
    na.rm = TRUE
  )
]

candidate_matches[
  ,
  overlap_bp := fifelse(
    coordinate_overlap,
    overlap_end - overlap_start + 1L,
    0L
  )
]

candidate_matches[
  ,
  file_width := file_region_end - file_region_start + 1L
]

candidate_matches[
  ,
  canonical_width :=
    canonical_region_end - canonical_region_start + 1L
]

candidate_matches[
  ,
  reciprocal_file_overlap := fifelse(
    coordinate_overlap & file_width > 0L,
    overlap_bp / file_width,
    0
  )
]

candidate_matches[
  ,
  reciprocal_canonical_overlap := fifelse(
    coordinate_overlap & canonical_width > 0L,
    overlap_bp / canonical_width,
    0
  )
]

candidate_matches[
  ,
  coordinate_distance :=
    abs(file_region_start - canonical_region_start) +
    abs(file_region_end - canonical_region_end)
]


###############################################################################
# 12. IDENTIFY VALID, UNMATCHED AND AMBIGUOUS FILES
###############################################################################

valid_candidates <- candidate_matches[
  coordinate_overlap == TRUE
]

candidate_counts <- valid_candidates[
  ,
  .(
    valid_candidate_count = .N
  ),
  by = output_file
]

file_inventory <- merge(
  file_inventory,
  candidate_counts,
  by = "output_file",
  all.x = TRUE,
  sort = FALSE
)

file_inventory[
  is.na(valid_candidate_count),
  valid_candidate_count := 0L
]

unmatched_files <- file_inventory[
  valid_candidate_count == 0L
]

ambiguous_files <- file_inventory[
  valid_candidate_count > 1L
]

cat("Candidate matching QC:\n")
cat("  Files with exactly one valid match:",
    file_inventory[valid_candidate_count == 1L, .N], "\n")
cat("  Files with no valid match:", nrow(unmatched_files), "\n")
cat("  Files with multiple valid matches:", nrow(ambiguous_files), "\n\n")

if (nrow(unmatched_files) > 0L) {

  cat("Unmatched files:\n")
  print(
    unmatched_files[
      ,
      .(
        filename,
        tissue,
        file_clump_id,
        file_chromosome,
        file_region_start,
        file_region_end
      )
    ]
  )
  cat("\n")
}

if (nrow(ambiguous_files) > 0L) {

  cat("Ambiguous files:\n")
  print(
    ambiguous_files[
      ,
      .(
        filename,
        tissue,
        file_clump_id,
        file_chromosome,
        file_region_start,
        file_region_end,
        valid_candidate_count
      )
    ]
  )
  cat("\n")
}


###############################################################################
# 13. SELECT UNIQUE DIRECT MATCHES
###############################################################################

selected_matches <- valid_candidates[
  output_file %chin%
    file_inventory[valid_candidate_count == 1L, output_file]
]

selected_matches[
  ,
  mapping_status := "DIRECT_CLUMP_CHROMOSOME_OVERLAP_MATCH"
]

selected_matches[
  ,
  extraction_status := fifelse(
    output_file_exists &
      output_file_nonempty,
    "EXTRACTED",
    "FILE_INVALID"
  )
]

selected_matches[
  ,
  chromosome_matches_canonical :=
    normalize_chromosome(canonical_chromosome) ==
    normalize_chromosome(
      sub(
        "^.*_chr([^_]+)_.*$",
        "\\1",
        filename
      )
    )
]

chromosome_mismatches <- selected_matches[
  is.na(chromosome_matches_canonical) |
    chromosome_matches_canonical == FALSE
]

if (nrow(chromosome_mismatches) > 0L) {

  cat("CRITICAL chromosome mismatches:\n")
  print(chromosome_mismatches)
  cat("\n")
}


###############################################################################
# 14. BUILD COMPLETE CLUMP × TISSUE GRID
###############################################################################

canonical_grid <- CJ(
  canonical_row = seq_len(nrow(canonical_clumps)),
  tissue = expected_tissues,
  unique = TRUE
)

canonical_grid <- merge(
  canonical_grid,
  canonical_clumps[
    ,
    .(
      canonical_row = .I,
      clump_id,
      canonical_chromosome,
      canonical_region_start,
      canonical_region_end,
      source_master_row,
      canonical_key
    )
  ],
  by = "canonical_row",
  all.x = TRUE,
  sort = FALSE
)

selected_for_join <- selected_matches[
  ,
  .(
    clump_id,
    canonical_chromosome,
    tissue,
    canonical_region_start_matched = canonical_region_start,
    canonical_region_end_matched = canonical_region_end,
    output_file,
    filename,
    file_region_start,
    file_region_end,
    output_file_exists,
    output_file_size_bytes,
    output_file_nonempty,
    coordinate_overlap,
    overlap_bp,
    reciprocal_file_overlap,
    reciprocal_canonical_overlap,
    coordinate_distance,
    mapping_status,
    extraction_status,
    chromosome_matches_canonical
  )
]

region_map <- merge(
  canonical_grid,
  selected_for_join,
  by = c(
    "clump_id",
    "canonical_chromosome",
    "tissue"
  ),
  all.x = TRUE,
  sort = FALSE
)

region_map[
  is.na(output_file),
  `:=`(
    mapping_status = "NO_APAQTL_FILE",
    extraction_status = "NO_APAQTL_FILE",
    output_file_exists = FALSE,
    output_file_nonempty = FALSE,
    coordinate_overlap = NA,
    chromosome_matches_canonical = NA
  )
]

region_map[
  ,
  canonical_coordinates_match :=
    is.na(canonical_region_start_matched) |
    (
      canonical_region_start ==
        canonical_region_start_matched &
      canonical_region_end ==
        canonical_region_end_matched
    )
]

region_map[
  ,
  output_file := fifelse(
    is.na(output_file),
    NA_character_,
    output_file
  )
]

region_map[
  ,
  output_path := output_file
]

region_map[
  ,
  clump_tissue_key := paste(
    clump_id,
    tissue,
    sep = "||"
  )
]

region_map[
  ,
  chromosome_sort := chromosome_sort_value(
    canonical_chromosome
  )
]

setorder(
  region_map,
  chromosome_sort,
  canonical_region_start,
  clump_id,
  tissue
)

region_map[
  ,
  region_map_row := seq_len(.N)
]

region_map[
  ,
  chromosome_sort := NULL
]


###############################################################################
# 15. FINAL STRUCTURAL QC
###############################################################################

duplicate_assignments <- region_map[
  duplicated(clump_tissue_key) |
    duplicated(clump_tissue_key, fromLast = TRUE)
]

wrong_chromosome_assignments <- region_map[
  extraction_status == "EXTRACTED" &
    (
      is.na(chromosome_matches_canonical) |
      chromosome_matches_canonical == FALSE
    )
]

nonoverlapping_assignments <- region_map[
  extraction_status == "EXTRACTED" &
    (
      is.na(coordinate_overlap) |
      coordinate_overlap == FALSE
    )
]

canonical_coordinate_conflicts <- region_map[
  canonical_coordinates_match == FALSE
]

extracted_rows <- region_map[
  extraction_status == "EXTRACTED"
]

non_extracted_rows <- region_map[
  extraction_status != "EXTRACTED"
]

expected_clumps <- 265L
expected_tissues_count <- 2L
expected_total_rows <- expected_clumps * expected_tissues_count
expected_extracted_files <- 516L
expected_non_extracted_rows <- 14L

mapping_qc <- rbindlist(
  list(

    qc_row(
      metric = "canonical_clumps",
      observed = nrow(canonical_clumps),
      expected = expected_clumps,
      status = ifelse(
        nrow(canonical_clumps) == expected_clumps,
        "PASS",
        "REVIEW"
      ),
      interpretation =
        "Number of canonical publication loci in the Module 11.2 manifest."
    ),

    qc_row(
      metric = "expected_tissues",
      observed = length(expected_tissues),
      expected = expected_tissues_count,
      status = ifelse(
        length(expected_tissues) == expected_tissues_count,
        "PASS",
        "FAIL"
      ),
      interpretation =
        "The catalogue uses the two GTEx skin tissues."
    ),

    qc_row(
      metric = "complete_clump_tissue_rows",
      observed = nrow(region_map),
      expected = expected_total_rows,
      status = ifelse(
        nrow(region_map) == expected_total_rows,
        "PASS",
        "FAIL"
      ),
      interpretation =
        "Complete canonical CLUMP × tissue grid."
    ),

    qc_row(
      metric = "discovered_apaqtl_files",
      observed = nrow(file_inventory),
      expected = expected_extracted_files,
      status = ifelse(
        nrow(file_inventory) == expected_extracted_files,
        "PASS",
        "REVIEW"
      ),
      interpretation =
        "Regional apaQTL files found in the two tissue directories."
    ),

    qc_row(
      metric = "uniquely_mapped_files",
      observed = nrow(selected_matches),
      expected = nrow(file_inventory),
      status = ifelse(
        nrow(selected_matches) == nrow(file_inventory),
        "PASS",
        "FAIL"
      ),
      interpretation =
        paste(
          "Every discovered file must map uniquely using",
          "CLUMP ID + chromosome + coordinate overlap."
        )
    ),

    qc_row(
      metric = "extracted_records",
      observed = nrow(extracted_rows),
      expected = expected_extracted_files,
      status = ifelse(
        nrow(extracted_rows) == expected_extracted_files,
        "PASS",
        "REVIEW"
      ),
      interpretation =
        "Mapped records with an existing non-empty apaQTL file."
    ),

    qc_row(
      metric = "non_extracted_records",
      observed = nrow(non_extracted_rows),
      expected = expected_non_extracted_rows,
      status = ifelse(
        nrow(non_extracted_rows) == expected_non_extracted_rows,
        "PASS",
        "REVIEW"
      ),
      interpretation =
        "Canonical CLUMP × tissue records without an apaQTL file."
    ),

    qc_row(
      metric = "filename_parse_failures",
      observed = nrow(invalid_filenames),
      expected = 0,
      status = ifelse(
        nrow(invalid_filenames) == 0L,
        "PASS",
        "FAIL"
      ),
      interpretation =
        "All regional filenames must be parsed successfully."
    ),

    qc_row(
      metric = "duplicate_file_keys",
      observed = nrow(duplicate_files),
      expected = 0,
      status = ifelse(
        nrow(duplicate_files) == 0L,
        "PASS",
        "FAIL"
      ),
      interpretation =
        "No duplicate CLUMP + chromosome + tissue files."
    ),

    qc_row(
      metric = "unmatched_files",
      observed = nrow(unmatched_files),
      expected = 0,
      status = ifelse(
        nrow(unmatched_files) == 0L,
        "PASS",
        "FAIL"
      ),
      interpretation =
        "Every apaQTL file must overlap its canonical CLUMP region."
    ),

    qc_row(
      metric = "ambiguous_files",
      observed = nrow(ambiguous_files),
      expected = 0,
      status = ifelse(
        nrow(ambiguous_files) == 0L,
        "PASS",
        "FAIL"
      ),
      interpretation =
        "Every apaQTL file must map to exactly one canonical region."
    ),

    qc_row(
      metric = "chromosome_mismatch_assignments",
      observed = nrow(wrong_chromosome_assignments),
      expected = 0,
      status = ifelse(
        nrow(wrong_chromosome_assignments) == 0L,
        "PASS",
        "FAIL"
      ),
      interpretation =
        "File chromosome must equal canonical chromosome."
    ),

    qc_row(
      metric = "nonoverlapping_assignments",
      observed = nrow(nonoverlapping_assignments),
      expected = 0,
      status = ifelse(
        nrow(nonoverlapping_assignments) == 0L,
        "PASS",
        "FAIL"
      ),
      interpretation =
        "Extracted file region must overlap the canonical region."
    ),

    qc_row(
      metric = "canonical_coordinate_conflicts",
      observed = nrow(canonical_coordinate_conflicts),
      expected = 0,
      status = ifelse(
        nrow(canonical_coordinate_conflicts) == 0L,
        "PASS",
        "FAIL"
      ),
      interpretation =
        "Joined canonical coordinates must remain unchanged."
    ),

    qc_row(
      metric = "duplicate_clump_tissue_assignments",
      observed = nrow(duplicate_assignments),
      expected = 0,
      status = ifelse(
        nrow(duplicate_assignments) == 0L,
        "PASS",
        "FAIL"
      ),
      interpretation =
        "Each CLUMP × tissue pair must occur exactly once."
    )
  ),
  use.names = TRUE,
  fill = TRUE
)


###############################################################################
# 16. DETERMINE FINAL STATUS
###############################################################################

hard_failures <- mapping_qc[
  status == "FAIL",
  metric
]

review_metrics <- mapping_qc[
  status == "REVIEW",
  metric
]

if (length(hard_failures) > 0L) {

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

status_table <- data.table(
  module = "11.3A_v4",
  module_name = "Robust direct CLUMP region mapping",
  pipeline_version = "4.0.0",
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
  hard_failure_count = length(hard_failures),
  review_count = length(review_metrics),
  canonical_clumps = nrow(canonical_clumps),
  complete_region_map_rows = nrow(region_map),
  discovered_apaqtl_files = nrow(file_inventory),
  uniquely_mapped_files = nrow(selected_matches),
  extracted_records = nrow(extracted_rows),
  non_extracted_records = nrow(non_extracted_rows),
  unmatched_files = nrow(unmatched_files),
  ambiguous_files = nrow(ambiguous_files),
  chromosome_mismatches =
    nrow(wrong_chromosome_assignments)
)


###############################################################################
# 17. PREPARE FINAL OUTPUT COLUMNS
###############################################################################

final_region_map <- region_map[
  ,
  .(
    clump_id,
    tissue,
    canonical_chromosome,
    canonical_region_start,
    canonical_region_end,
    mapping_status,
    extraction_status,
    output_file,
    output_path,
    output_file_exists,
    output_file_size_bytes,
    output_file_nonempty,
    file_region_start,
    file_region_end,
    coordinate_overlap,
    overlap_bp,
    reciprocal_file_overlap,
    reciprocal_canonical_overlap,
    coordinate_distance,
    chromosome_matches_canonical,
    canonical_coordinates_match,
    source_master_row,
    canonical_key,
    clump_tissue_key,
    region_map_row
  )
]


###############################################################################
# 18. WRITE OUTPUTS
###############################################################################

fwrite(
  final_region_map,
  output_map_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  file_inventory,
  file_inventory_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  mapping_qc,
  mapping_qc_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

write_empty_or_table(
  unmatched_files,
  unmatched_files_output,
  template_columns = names(file_inventory)
)

write_empty_or_table(
  ambiguous_files,
  ambiguous_files_output,
  template_columns = names(file_inventory)
)

write_empty_or_table(
  wrong_chromosome_assignments,
  chromosome_mismatch_output,
  template_columns = names(region_map)
)

write_empty_or_table(
  duplicate_assignments,
  duplicate_assignments_output,
  template_columns = names(region_map)
)

fwrite(
  status_table,
  status_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


###############################################################################
# 19. SESSION INFORMATION
###############################################################################

session_lines <- capture.output(
  sessionInfo()
)

writeLines(
  c(
    "PigmentationAtlas Module 11.3A v4",
    paste("Start:", format(module_start_time)),
    paste("End:", format(module_end_time)),
    paste("Final status:", final_status),
    "",
    session_lines
  ),
  con = session_info_output
)


###############################################################################
# 20. FINAL REPORT
###############################################################################

cat("\n")
cat("============================================================\n")
cat("MODULE 11.3A v4 – FINAL SUMMARY\n")
cat("============================================================\n\n")

print(mapping_qc)

cat("\nCanonical chromosome distribution:\n")

print(
  final_region_map[
    ,
    .(
      records = .N,
      unique_clumps = uniqueN(clump_id),
      extracted_records =
        sum(extraction_status == "EXTRACTED")
    ),
    by = canonical_chromosome
  ][
    ,
    chromosome_sort :=
      chromosome_sort_value(canonical_chromosome)
  ][
    order(chromosome_sort)
  ][
    ,
    chromosome_sort := NULL
  ]
)

cat("\nTissue distribution:\n")

print(
  final_region_map[
    ,
    .(
      records = .N,
      unique_clumps = uniqueN(clump_id),
      extracted_records =
        sum(extraction_status == "EXTRACTED"),
      non_extracted_records =
        sum(extraction_status != "EXTRACTED")
    ),
    by = tissue
  ][
    order(tissue)
  ]
)

cat("\nFinal status:", final_status, "\n")
cat("Runtime:", round(runtime_seconds, 2), "seconds\n\n")

cat("Main output:\n")
cat(output_map_file, "\n\n")

cat("QC directory:\n")
cat(qc_dir, "\n\n")

if (identical(final_status, "PASS")) {

  cat("FINAL QC RESULT: PASS\n")
  cat("The Module 11.3A v4 region map is valid.\n")

} else if (identical(final_status, "PASS_WITH_REVIEW")) {

  cat("FINAL QC RESULT: PASS WITH REVIEW\n")
  cat("No structural error was detected, but expected counts differ.\n")
  cat("Review the QC table before continuing.\n")

} else {

  cat("FINAL QC RESULT: FAIL\n")
  cat("Do NOT start Module 11.3B.1.\n")
  cat("Failed metrics:\n")

  cat(
    paste0(
      "  - ",
      hard_failures,
      collapse = "\n"
    ),
    "\n"
  )
}

cat("============================================================\n\n")


###############################################################################
# 21. CLOSE LOGGING
###############################################################################

while (sink.number(type = "output") > 0L) {
  try(sink(type = "output"), silent = TRUE)
}

while (sink.number(type = "output") > 0L) {
  try(sink(type = "output"), silent = TRUE)
}

try(close(log_connection), silent = TRUE)

logging_active <- FALSE


###############################################################################
# 22. STOP ON HARD FAILURE
###############################################################################

if (identical(final_status, "FAIL")) {

  stop(
    paste0(
      "Module 11.3A v4 failed. ",
      "Review QC outputs under:\n",
      qc_dir
    ),
    call. = FALSE
  )
}

invisible(
  list(
    canonical_clumps = canonical_clumps,
    file_inventory = file_inventory,
    region_map = final_region_map,
    mapping_qc = mapping_qc,
    status = status_table
  )
)
