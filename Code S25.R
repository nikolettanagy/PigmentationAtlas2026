###############################################################################
# PigmentationAtlas
# Module 11.3B.4
#
# FINAL QC OF THE HARMONIZED apaQTL CATALOGUE
#
# Purpose:
#   - validate all 516 harmonized chunks
#   - identify the exact reasons for PASS_WITH_ROW_QC_FLAGS
#   - verify row counts and catalogue integrity
#   - summarize row-level QC flags
#   - verify the final concatenated gzip catalogue
#
# Memory strategy:
#   Each chunk is read and summarized separately.
###############################################################################


###############################################################################
# 1. CLEAN SESSION
###############################################################################

rm(list = ls())
invisible(gc())

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  width = 240
)

module_start_time <- Sys.time()

cat("\n")
cat("============================================================\n")
cat("PigmentationAtlas\n")
cat("MODULE 11.3B.4 – FINAL apaQTL CATALOGUE QC\n")
cat("Started:", format(module_start_time), "\n")
cat("============================================================\n\n")


###############################################################################
# 2. PACKAGE
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

qc_dir <- file.path(
  catalogue_dir,
  "qc",
  "final_qc"
)

log_dir <- file.path(
  catalogue_dir,
  "logs"
)

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

work_queue_file <- file.path(
  catalogue_dir,
  "module11_3B_work_queue.tsv"
)

file_summary_file <- file.path(
  catalogue_dir,
  "module11_3B_file_summary.tsv"
)

build_status_file <- file.path(
  catalogue_dir,
  "module11_3B_build_status.tsv"
)

final_catalogue_file <- file.path(
  catalogue_dir,
  "module11_3B_harmonized_apaQTL.tsv.gz"
)

file_qc_output <- file.path(
  qc_dir,
  "module11_3B_4_file_level_qc.tsv"
)

flag_reason_output <- file.path(
  qc_dir,
  "module11_3B_4_flag_reason_summary.tsv"
)

row_flag_output <- file.path(
  qc_dir,
  "module11_3B_4_row_flag_summary.tsv"
)

tissue_summary_output <- file.path(
  qc_dir,
  "module11_3B_4_tissue_summary.tsv"
)

column_consistency_output <- file.path(
  qc_dir,
  "module11_3B_4_column_consistency.tsv"
)

catalogue_integrity_output <- file.path(
  qc_dir,
  "module11_3B_4_catalogue_integrity.tsv"
)

final_status_output <- file.path(
  catalogue_dir,
  "module11_3B_4_final_qc_status.tsv"
)

session_info_output <- file.path(
  log_dir,
  "module11_3B_4_session_info.txt"
)

log_file <- file.path(
  log_dir,
  paste0(
    "module11_3B_4_final_qc_",
    format(module_start_time, "%Y%m%d_%H%M%S"),
    ".log"
  )
)


###############################################################################
# 4. HELPER FUNCTIONS
###############################################################################

safe_filename <- function(x) {

  x <- gsub(
    "[^A-Za-z0-9._-]+",
    "_",
    as.character(x)
  )

  gsub("_+", "_", x)
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


is_missing_character <- function(x) {

  is.na(x) |
    trimws(as.character(x)) == ""
}


safe_sum <- function(x) {

  sum(x, na.rm = TRUE)
}


count_gzip_lines <- function(file_path,
                             block_size = 100000L) {

  connection <- gzfile(
    file_path,
    open = "rt"
  )

  on.exit(
    try(close(connection), silent = TRUE),
    add = TRUE
  )

  total_lines <- 0

  repeat {

    block <- readLines(
      connection,
      n = block_size,
      warn = FALSE
    )

    block_length <- length(block)

    if (block_length == 0L) {
      break
    }

    total_lines <- total_lines + block_length
  }

  total_lines
}


count_true_if_present <- function(dt,
                                  column_name) {

  if (!column_name %in% names(dt)) {
    return(NA_integer_)
  }

  sum(
    dt[[column_name]] %in% TRUE,
    na.rm = TRUE
  )
}


count_false_if_present <- function(dt,
                                   column_name) {

  if (!column_name %in% names(dt)) {
    return(NA_integer_)
  }

  sum(
    dt[[column_name]] %in% FALSE,
    na.rm = TRUE
  )
}


###############################################################################
# 5. CHECK REQUIRED INPUTS
###############################################################################

required_files <- c(
  work_queue_file,
  file_summary_file,
  build_status_file,
  final_catalogue_file
)

missing_required_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_required_files) > 0L) {

  stop(
    paste0(
      "Required input file(s) missing:\n",
      paste(
        missing_required_files,
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
}

cat("All required Module 11.3B.3 outputs are present.\n\n")


###############################################################################
# 6. READ CONTROL TABLES
###############################################################################

work_queue <- fread(
  work_queue_file,
  sep = "\t",
  header = TRUE,
  na.strings = c("", "NA", "NaN", "NULL"),
  showProgress = FALSE
)

file_summary <- fread(
  file_summary_file,
  sep = "\t",
  header = TRUE,
  na.strings = c("", "NA", "NaN", "NULL"),
  showProgress = FALSE
)

build_status <- fread(
  build_status_file,
  sep = "\t",
  header = TRUE,
  na.strings = c("", "NA", "NaN", "NULL"),
  showProgress = FALSE
)

cat("Work-queue rows:", nrow(work_queue), "\n")
cat("File-summary rows:", nrow(file_summary), "\n")
cat("Build status:", build_status$final_status[1L], "\n\n")


###############################################################################
# 7. ADD EXPECTED CHUNK PATHS
###############################################################################

work_queue[
  ,
  chunk_file := mapply(
    make_chunk_path,
    queue_id,
    clump_id,
    tissue,
    USE.NAMES = FALSE
  )
]

work_queue[
  ,
  chunk_exists := file.exists(chunk_file)
]

work_queue[
  ,
  chunk_size_bytes := fifelse(
    chunk_exists,
    as.numeric(file.info(chunk_file)$size),
    NA_real_
  )
]

work_queue[
  ,
  chunk_nonempty := (
    chunk_exists &
      !is.na(chunk_size_bytes) &
      chunk_size_bytes > 0
  )
]


###############################################################################
# 8. INITIALIZE LOGGING
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
        try(
          sink(type = "output"),
          silent = TRUE
        )
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
# 9. EXPECTED HARMONIZED COLUMNS
###############################################################################

expected_columns <- c(
  "queue_id",
  "clump_id",
  "tissue",
  "canonical_chromosome",
  "canonical_region_start",
  "canonical_region_end",
  "source_file",
  "variant_id",
  "rsid",
  "chromosome",
  "position",
  "ref_allele",
  "alt_allele",
  "phenotype_id",
  "gene_id",
  "p_value",
  "beta",
  "standard_error",
  "maf",
  "allele_frequency",
  "sample_size",
  "molecular_trait_chromosome",
  "molecular_trait_position",
  "distance",
  "chromosome_matches_canonical",
  "position_within_canonical_region",
  "duplicate_record"
)


###############################################################################
# 10. PROCESS EVERY HARMONIZED CHUNK
###############################################################################

total_files <- nrow(work_queue)

file_qc_list <- vector(
  mode = "list",
  length = total_files
)

column_qc_list <- vector(
  mode = "list",
  length = total_files
)

cat("Starting chunk-level final QC...\n\n")

for (i in seq_len(total_files)) {

  metadata <- work_queue[i]

  queue_id <- metadata$queue_id[1L]
  clump_id <- metadata$clump_id[1L]
  tissue <- metadata$tissue[1L]
  chunk_file <- metadata$chunk_file[1L]

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

  if (
    !file.exists(chunk_file) ||
      file.info(chunk_file)$size <= 0
  ) {

    file_qc_list[[i]] <- data.table(
      queue_id = queue_id,
      clump_id = clump_id,
      tissue = tissue,
      chunk_file = chunk_file,
      chunk_size_bytes = NA_real_,
      rows = NA_integer_,
      columns = NA_integer_,
      missing_required_columns = NA_integer_,
      missing_variant_id = NA_integer_,
      missing_chromosome = NA_integer_,
      missing_position = NA_integer_,
      missing_ref_allele = NA_integer_,
      missing_alt_allele = NA_integer_,
      missing_phenotype_id = NA_integer_,
      missing_gene_id = NA_integer_,
      missing_p_value = NA_integer_,
      invalid_p_value = NA_integer_,
      missing_beta = NA_integer_,
      missing_standard_error = NA_integer_,
      invalid_standard_error = NA_integer_,
      chromosome_mismatch_rows = NA_integer_,
      outside_region_rows = NA_integer_,
      duplicate_rows = NA_integer_,
      invalid_maf_rows = NA_integer_,
      invalid_allele_frequency_rows = NA_integer_,
      metadata_queue_mismatch_rows = NA_integer_,
      metadata_clump_mismatch_rows = NA_integer_,
      metadata_tissue_mismatch_rows = NA_integer_,
      qc_flagged_rows = NA_integer_,
      final_file_qc = "FAIL_MISSING_CHUNK"
    )

    cat("  Result: FAIL_MISSING_CHUNK\n")
    next
  }

  dt <- tryCatch(

    fread(
      chunk_file,
      sep = "\t",
      header = TRUE,
      na.strings = c("", "NA", "NaN", "NULL"),
      showProgress = FALSE
    ),

    error = function(e) e
  )

  if (inherits(dt, "error")) {

    file_qc_list[[i]] <- data.table(
      queue_id = queue_id,
      clump_id = clump_id,
      tissue = tissue,
      chunk_file = chunk_file,
      chunk_size_bytes = as.numeric(
        file.info(chunk_file)$size
      ),
      rows = NA_integer_,
      columns = NA_integer_,
      missing_required_columns = NA_integer_,
      missing_variant_id = NA_integer_,
      missing_chromosome = NA_integer_,
      missing_position = NA_integer_,
      missing_ref_allele = NA_integer_,
      missing_alt_allele = NA_integer_,
      missing_phenotype_id = NA_integer_,
      missing_gene_id = NA_integer_,
      missing_p_value = NA_integer_,
      invalid_p_value = NA_integer_,
      missing_beta = NA_integer_,
      missing_standard_error = NA_integer_,
      invalid_standard_error = NA_integer_,
      chromosome_mismatch_rows = NA_integer_,
      outside_region_rows = NA_integer_,
      duplicate_rows = NA_integer_,
      invalid_maf_rows = NA_integer_,
      invalid_allele_frequency_rows = NA_integer_,
      metadata_queue_mismatch_rows = NA_integer_,
      metadata_clump_mismatch_rows = NA_integer_,
      metadata_tissue_mismatch_rows = NA_integer_,
      qc_flagged_rows = NA_integer_,
      final_file_qc = "FAIL_UNREADABLE_CHUNK"
    )

    cat(
      "  Result: FAIL_UNREADABLE_CHUNK |",
      conditionMessage(dt),
      "\n"
    )

    next
  }

  actual_columns <- names(dt)

  missing_required_columns <- setdiff(
    expected_columns,
    actual_columns
  )

  column_qc_list[[i]] <- data.table(
    queue_id = queue_id,
    clump_id = clump_id,
    tissue = tissue,
    expected_column = expected_columns,
    present = expected_columns %in% actual_columns
  )

  n_rows <- nrow(dt)

  missing_variant_id <- if ("variant_id" %in% actual_columns) {
    sum(is_missing_character(dt$variant_id))
  } else {
    n_rows
  }

  missing_chromosome <- if ("chromosome" %in% actual_columns) {
    sum(is_missing_character(dt$chromosome))
  } else {
    n_rows
  }

  missing_position <- if ("position" %in% actual_columns) {
    sum(is.na(dt$position))
  } else {
    n_rows
  }

  missing_ref_allele <- if ("ref_allele" %in% actual_columns) {
    sum(is_missing_character(dt$ref_allele))
  } else {
    n_rows
  }

  missing_alt_allele <- if ("alt_allele" %in% actual_columns) {
    sum(is_missing_character(dt$alt_allele))
  } else {
    n_rows
  }

  missing_phenotype_id <- if ("phenotype_id" %in% actual_columns) {
    sum(is_missing_character(dt$phenotype_id))
  } else {
    n_rows
  }

  missing_gene_id <- if ("gene_id" %in% actual_columns) {
    sum(is_missing_character(dt$gene_id))
  } else {
    n_rows
  }

  missing_p_value <- if ("p_value" %in% actual_columns) {
    sum(is.na(dt$p_value))
  } else {
    n_rows
  }

  invalid_p_value <- if ("p_value" %in% actual_columns) {
    sum(
      !is.na(dt$p_value) &
        (
          !is.finite(dt$p_value) |
            dt$p_value < 0 |
            dt$p_value > 1
        )
    )
  } else {
    n_rows
  }

  missing_beta <- if ("beta" %in% actual_columns) {
    sum(is.na(dt$beta))
  } else {
    n_rows
  }

  missing_standard_error <- if ("standard_error" %in% actual_columns) {
    sum(is.na(dt$standard_error))
  } else {
    n_rows
  }

  invalid_standard_error <- if ("standard_error" %in% actual_columns) {
    sum(
      !is.na(dt$standard_error) &
        (
          !is.finite(dt$standard_error) |
            dt$standard_error <= 0
        )
    )
  } else {
    n_rows
  }

  chromosome_mismatch_rows <- count_false_if_present(
    dt,
    "chromosome_matches_canonical"
  )

  outside_region_rows <- count_false_if_present(
    dt,
    "position_within_canonical_region"
  )

  duplicate_rows <- count_true_if_present(
    dt,
    "duplicate_record"
  )

  invalid_maf_rows <- if ("maf" %in% actual_columns) {

    sum(
      !is.na(dt$maf) &
        (
          !is.finite(dt$maf) |
            dt$maf < 0 |
            dt$maf > 0.5
        )
    )

  } else {

    NA_integer_
  }

  invalid_allele_frequency_rows <- if (
    "allele_frequency" %in% actual_columns
  ) {

    sum(
      !is.na(dt$allele_frequency) &
        (
          !is.finite(dt$allele_frequency) |
            dt$allele_frequency < 0 |
            dt$allele_frequency > 1
        )
    )

  } else {

    NA_integer_
  }

  metadata_queue_mismatch_rows <- if (
    "queue_id" %in% actual_columns
  ) {

    sum(
      is.na(dt$queue_id) |
        dt$queue_id != queue_id
    )

  } else {

    n_rows
  }

  metadata_clump_mismatch_rows <- if (
    "clump_id" %in% actual_columns
  ) {

    sum(
      is.na(dt$clump_id) |
        dt$clump_id != clump_id
    )

  } else {

    n_rows
  }

  metadata_tissue_mismatch_rows <- if (
    "tissue" %in% actual_columns
  ) {

    sum(
      is.na(dt$tissue) |
        dt$tissue != tissue
    )

  } else {

    n_rows
  }

  row_flag_matrix <- cbind(
    chromosome_mismatch = if (
      "chromosome_matches_canonical" %in% actual_columns
    ) {
      dt$chromosome_matches_canonical %in% FALSE
    } else {
      rep(TRUE, n_rows)
    },

    outside_region = if (
      "position_within_canonical_region" %in% actual_columns
    ) {
      dt$position_within_canonical_region %in% FALSE
    } else {
      rep(TRUE, n_rows)
    },

    duplicate_record = if (
      "duplicate_record" %in% actual_columns
    ) {
      dt$duplicate_record %in% TRUE
    } else {
      rep(TRUE, n_rows)
    },

    invalid_p_value = if ("p_value" %in% actual_columns) {
      is.na(dt$p_value) |
        !is.finite(dt$p_value) |
        dt$p_value < 0 |
        dt$p_value > 1
    } else {
      rep(TRUE, n_rows)
    },

    missing_beta = if ("beta" %in% actual_columns) {
      is.na(dt$beta)
    } else {
      rep(TRUE, n_rows)
    },

    invalid_standard_error = if (
      "standard_error" %in% actual_columns
    ) {
      is.na(dt$standard_error) |
        !is.finite(dt$standard_error) |
        dt$standard_error <= 0
    } else {
      rep(TRUE, n_rows)
    },

    invalid_allele_frequency = if (
      "allele_frequency" %in% actual_columns
    ) {
      !is.na(dt$allele_frequency) &
        (
          !is.finite(dt$allele_frequency) |
            dt$allele_frequency < 0 |
            dt$allele_frequency > 1
        )
    } else {
      rep(FALSE, n_rows)
    }
  )

  qc_flagged_rows <- sum(
    rowSums(row_flag_matrix) > 0
  )

  critical_fail_count <- sum(
    length(missing_required_columns),
    missing_variant_id,
    missing_chromosome,
    missing_position,
    missing_phenotype_id,
    missing_gene_id,
    missing_p_value,
    invalid_p_value,
    metadata_queue_mismatch_rows,
    metadata_clump_mismatch_rows,
    metadata_tissue_mismatch_rows,
    na.rm = TRUE
  )

  warning_count <- sum(
    missing_ref_allele,
    missing_alt_allele,
    missing_beta,
    missing_standard_error,
    invalid_standard_error,
    chromosome_mismatch_rows,
    outside_region_rows,
    duplicate_rows,
    invalid_maf_rows,
    invalid_allele_frequency_rows,
    na.rm = TRUE
  )

  final_file_qc <- if (critical_fail_count > 0L) {

    "FAIL"

  } else if (warning_count > 0L) {

    "PASS_WITH_FLAGS"

  } else {

    "PASS"
  }

  file_qc_list[[i]] <- data.table(
    queue_id = queue_id,
    clump_id = clump_id,
    tissue = tissue,
    chunk_file = chunk_file,
    chunk_size_bytes = as.numeric(
      file.info(chunk_file)$size
    ),
    rows = n_rows,
    columns = ncol(dt),
    missing_required_columns = length(
      missing_required_columns
    ),
    missing_variant_id = missing_variant_id,
    missing_chromosome = missing_chromosome,
    missing_position = missing_position,
    missing_ref_allele = missing_ref_allele,
    missing_alt_allele = missing_alt_allele,
    missing_phenotype_id = missing_phenotype_id,
    missing_gene_id = missing_gene_id,
    missing_p_value = missing_p_value,
    invalid_p_value = invalid_p_value,
    missing_beta = missing_beta,
    missing_standard_error = missing_standard_error,
    invalid_standard_error = invalid_standard_error,
    chromosome_mismatch_rows = chromosome_mismatch_rows,
    outside_region_rows = outside_region_rows,
    duplicate_rows = duplicate_rows,
    invalid_maf_rows = invalid_maf_rows,
    invalid_allele_frequency_rows =
      invalid_allele_frequency_rows,
    metadata_queue_mismatch_rows =
      metadata_queue_mismatch_rows,
    metadata_clump_mismatch_rows =
      metadata_clump_mismatch_rows,
    metadata_tissue_mismatch_rows =
      metadata_tissue_mismatch_rows,
    qc_flagged_rows = qc_flagged_rows,
    final_file_qc = final_file_qc
  )

  cat(
    "  Result:",
    final_file_qc,
    "| rows:",
    format(n_rows, big.mark = ","),
    "| flagged rows:",
    format(qc_flagged_rows, big.mark = ","),
    "\n"
  )

  rm(dt, row_flag_matrix)
  invisible(gc())
}


###############################################################################
# 11. COMBINE FILE-LEVEL RESULTS
###############################################################################

file_qc <- rbindlist(
  file_qc_list,
  use.names = TRUE,
  fill = TRUE
)

setorder(
  file_qc,
  queue_id
)

fwrite(
  file_qc,
  file_qc_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

column_consistency <- rbindlist(
  column_qc_list[
    !vapply(
      column_qc_list,
      is.null,
      logical(1)
    )
  ],
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  column_consistency,
  column_consistency_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


###############################################################################
# 12. ROW-FLAG SUMMARY
###############################################################################

flag_columns <- c(
  "missing_variant_id",
  "missing_chromosome",
  "missing_position",
  "missing_ref_allele",
  "missing_alt_allele",
  "missing_phenotype_id",
  "missing_gene_id",
  "missing_p_value",
  "invalid_p_value",
  "missing_beta",
  "missing_standard_error",
  "invalid_standard_error",
  "chromosome_mismatch_rows",
  "outside_region_rows",
  "duplicate_rows",
  "invalid_maf_rows",
  "invalid_allele_frequency_rows",
  "metadata_queue_mismatch_rows",
  "metadata_clump_mismatch_rows",
  "metadata_tissue_mismatch_rows"
)

row_flag_summary <- rbindlist(
  lapply(
    flag_columns,
    function(flag_name) {

      values <- file_qc[[flag_name]]

      data.table(
        flag = flag_name,
        affected_files = sum(
          !is.na(values) &
            values > 0
        ),
        affected_rows = sum(
          values,
          na.rm = TRUE
        )
      )
    }
  )
)

row_flag_summary[
  ,
  flag_severity := fifelse(
    flag %in% c(
      "missing_variant_id",
      "missing_chromosome",
      "missing_position",
      "missing_phenotype_id",
      "missing_gene_id",
      "missing_p_value",
      "invalid_p_value",
      "metadata_queue_mismatch_rows",
      "metadata_clump_mismatch_rows",
      "metadata_tissue_mismatch_rows"
    ),
    "CRITICAL",
    "WARNING"
  )
]

setcolorder(
  row_flag_summary,
  c(
    "flag",
    "flag_severity",
    "affected_files",
    "affected_rows"
  )
)

setorder(
  row_flag_summary,
  -affected_rows,
  flag
)

fwrite(
  row_flag_summary,
  row_flag_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


###############################################################################
# 13. EXPLAIN ORIGINAL PASS_WITH_ROW_QC_FLAGS
###############################################################################

file_summary_subset <- file_summary[
  ,
  .(
    queue_id,
    original_file_status = file_status,
    original_missing_rsid = missing_rsid,
    original_missing_variant_id = missing_variant_id,
    original_missing_chromosome = missing_chromosome,
    original_missing_position = missing_position,
    original_missing_phenotype_id = missing_phenotype_id,
    original_missing_gene_id = missing_gene_id,
    original_missing_p_value = missing_p_value,
    original_invalid_p_value = invalid_p_value,
    original_missing_beta = missing_beta,
    original_missing_standard_error =
      missing_standard_error,
    original_chromosome_mismatch_rows =
      chromosome_mismatch_rows,
    original_outside_region_rows =
      outside_canonical_region_rows,
    original_duplicate_rows = duplicate_rows,
    original_invalid_maf_rows = invalid_maf_rows,
    original_invalid_allele_frequency_rows =
      invalid_allele_frequency_rows
  )
]

flag_reason_table <- merge(
  file_summary_subset,
  file_qc,
  by = "queue_id",
  all.x = TRUE
)

flag_reason_table[
  ,
  flag_reason := paste(
    c(
      if (
        !is.na(original_missing_rsid) &&
          original_missing_rsid > 0
      ) {
        "MISSING_RSID"
      } else {
        character()
      },

      if (
        !is.na(original_missing_variant_id) &&
          original_missing_variant_id > 0
      ) {
        "MISSING_VARIANT_ID"
      } else {
        character()
      },

      if (
        !is.na(original_missing_p_value) &&
          original_missing_p_value > 0
      ) {
        "MISSING_P_VALUE"
      } else {
        character()
      },

      if (
        !is.na(original_missing_beta) &&
          original_missing_beta > 0
      ) {
        "MISSING_BETA"
      } else {
        character()
      },

      if (
        !is.na(original_missing_standard_error) &&
          original_missing_standard_error > 0
      ) {
        "MISSING_STANDARD_ERROR"
      } else {
        character()
      },

      if (
        !is.na(original_chromosome_mismatch_rows) &&
          original_chromosome_mismatch_rows > 0
      ) {
        "CHROMOSOME_MISMATCH"
      } else {
        character()
      },

      if (
        !is.na(original_outside_region_rows) &&
          original_outside_region_rows > 0
      ) {
        "OUTSIDE_CANONICAL_REGION"
      } else {
        character()
      },

      if (
        !is.na(original_duplicate_rows) &&
          original_duplicate_rows > 0
      ) {
        "DUPLICATE_ROWS"
      } else {
        character()
      },

      if (
        !is.na(original_invalid_maf_rows) &&
          original_invalid_maf_rows > 0
      ) {
        "INVALID_MAF"
      } else {
        character()
      },

      if (
        !is.na(original_invalid_allele_frequency_rows) &&
          original_invalid_allele_frequency_rows > 0
      ) {
        "INVALID_ALLELE_FREQUENCY"
      } else {
        character()
      }
    ),
    collapse = ";"
  ),
  by = seq_len(nrow(flag_reason_table))
]

flag_reason_table[
  flag_reason == "",
  flag_reason := "NO_EXPLICIT_ROW_ERROR"
]

flag_reason_summary <- flag_reason_table[
  ,
  .(
    files = .N,
    rows = sum(rows, na.rm = TRUE)
  ),
  by = .(
    original_file_status,
    flag_reason
  )
][
  order(
    original_file_status,
    -files
  )
]

fwrite(
  flag_reason_summary,
  flag_reason_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


###############################################################################
# 14. TISSUE SUMMARY
###############################################################################

tissue_summary <- file_qc[
  ,
  .(
    files = .N,
    total_rows = sum(rows, na.rm = TRUE),
    total_size_bytes = sum(
      chunk_size_bytes,
      na.rm = TRUE
    ),
    pass_files = sum(
      final_file_qc == "PASS",
      na.rm = TRUE
    ),
    pass_with_flags_files = sum(
      final_file_qc == "PASS_WITH_FLAGS",
      na.rm = TRUE
    ),
    failed_files = sum(
      grepl("^FAIL", final_file_qc),
      na.rm = TRUE
    ),
    flagged_rows = sum(
      qc_flagged_rows,
      na.rm = TRUE
    ),
    unique_clumps = uniqueN(clump_id)
  ),
  by = tissue
][
  order(tissue)
]

fwrite(
  tissue_summary,
  tissue_summary_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


###############################################################################
# 15. FINAL CATALOGUE INTEGRITY CHECK
###############################################################################

cat("\nCounting lines in final gzip catalogue...\n")

final_catalogue_line_count <- count_gzip_lines(
  final_catalogue_file
)

expected_data_rows <- sum(
  file_qc$rows,
  na.rm = TRUE
)

expected_total_lines <- expected_data_rows + 1

catalogue_header <- fread(
  final_catalogue_file,
  sep = "\t",
  header = TRUE,
  nrows = 0,
  showProgress = FALSE
)

final_catalogue_columns <- names(catalogue_header)

catalogue_integrity <- data.table(
  metric = c(
    "expected_work_queue_files",
    "observed_file_qc_rows",
    "existing_nonempty_chunks",
    "expected_data_rows_from_chunks",
    "final_catalogue_total_lines",
    "expected_total_lines_header_plus_data",
    "final_catalogue_data_rows",
    "final_catalogue_column_count",
    "expected_column_count",
    "missing_expected_columns",
    "build_status_pass",
    "final_catalogue_exists",
    "final_catalogue_nonempty"
  ),

  observed = c(
    as.character(nrow(work_queue)),
    as.character(nrow(file_qc)),
    as.character(
      sum(
        work_queue$chunk_nonempty,
        na.rm = TRUE
      )
    ),
    as.character(expected_data_rows),
    as.character(final_catalogue_line_count),
    as.character(expected_total_lines),
    as.character(
      final_catalogue_line_count - 1
    ),
    as.character(
      length(final_catalogue_columns)
    ),
    as.character(length(expected_columns)),
    as.character(
      length(
        setdiff(
          expected_columns,
          final_catalogue_columns
        )
      )
    ),
    as.character(
      identical(
        build_status$final_status[1L],
        "PASS"
      )
    ),
    as.character(
      file.exists(final_catalogue_file)
    ),
    as.character(
      file.exists(final_catalogue_file) &&
        file.info(final_catalogue_file)$size > 0
    )
  ),

  expected = c(
    "516",
    "516",
    "516",
    as.character(expected_data_rows),
    as.character(expected_total_lines),
    as.character(expected_total_lines),
    as.character(expected_data_rows),
    as.character(length(expected_columns)),
    as.character(length(expected_columns)),
    "0",
    "TRUE",
    "TRUE",
    "TRUE"
  )
)

catalogue_integrity[
  ,
  status := "FAIL"
]

catalogue_integrity[
  metric == "final_catalogue_column_count",
  status := fifelse(
    length(setdiff(expected_columns, final_catalogue_columns)) == 0L,
    "PASS",
    "FAIL"
  )
]

catalogue_integrity[
  metric == "missing_expected_columns",
  status := fifelse(
    observed == "0",
    "PASS",
    "FAIL"
  )
]

catalogue_integrity[
  !metric %in% c(
    "final_catalogue_column_count",
    "missing_expected_columns"
  ),
  status := fifelse(
    observed == expected,
    "PASS",
    "FAIL"
  )
]

fwrite(
  catalogue_integrity,
  catalogue_integrity_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


###############################################################################
# 16. DETERMINE FINAL STATUS
###############################################################################

critical_flags <- row_flag_summary[
  flag_severity == "CRITICAL",
  sum(affected_rows, na.rm = TRUE)
]

warning_flags <- row_flag_summary[
  flag_severity == "WARNING",
  sum(affected_rows, na.rm = TRUE)
]

failed_file_count <- sum(
  grepl("^FAIL", file_qc$final_file_qc),
  na.rm = TRUE
)

integrity_fail_count <- sum(
  catalogue_integrity$status == "FAIL"
)

missing_column_count <- sum(
  column_consistency$present == FALSE,
  na.rm = TRUE
)

module_end_time <- Sys.time()

runtime_seconds <- as.numeric(
  difftime(
    module_end_time,
    module_start_time,
    units = "secs"
  )
)

final_status <- if (
  failed_file_count > 0L ||
    critical_flags > 0L ||
    integrity_fail_count > 0L ||
    missing_column_count > 0L
) {

  "FAIL"

} else if (warning_flags > 0L) {

  "PASS_WITH_WARNINGS"

} else {

  "PASS"
}

final_qc_status <- data.table(
  module = "11.3B.4",
  module_name = "final harmonized apaQTL catalogue QC",
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
  checked_files = nrow(file_qc),
  failed_files = failed_file_count,
  total_rows = expected_data_rows,
  files_with_original_qc_flags = sum(
    file_summary$file_status ==
      "PASS_WITH_ROW_QC_FLAGS"
  ),
  files_with_final_flags = sum(
    file_qc$final_file_qc ==
      "PASS_WITH_FLAGS"
  ),
  critical_flagged_rows = critical_flags,
  warning_flagged_rows = warning_flags,
  integrity_failures = integrity_fail_count,
  missing_expected_column_instances =
    missing_column_count,
  final_catalogue_lines =
    final_catalogue_line_count,
  final_catalogue_data_rows =
    final_catalogue_line_count - 1,
  final_status = final_status
)

fwrite(
  final_qc_status,
  final_status_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


###############################################################################
# 17. SESSION INFORMATION
###############################################################################

writeLines(
  c(
    "PigmentationAtlas Module 11.3B.4",
    paste("Start:", format(module_start_time)),
    paste("End:", format(module_end_time)),
    paste("Final status:", final_status),
    "",
    capture.output(sessionInfo())
  ),
  con = session_info_output
)


###############################################################################
# 18. FINAL REPORT
###############################################################################

cat("\n")
cat("============================================================\n")
cat("MODULE 11.3B.4 – FINAL QC SUMMARY\n")
cat("============================================================\n\n")

print(final_qc_status)

cat("\nFinal file-QC distribution:\n")

print(
  file_qc[
    ,
    .N,
    by = final_file_qc
  ][
    order(-N)
  ]
)

cat("\nOriginal Module 11.3B.3 file-status distribution:\n")

print(
  file_summary[
    ,
    .N,
    by = file_status
  ][
    order(-N)
  ]
)

cat("\nOriginal QC flag explanations:\n")

print(flag_reason_summary)

cat("\nRow-level QC flags:\n")

print(row_flag_summary)

cat("\nTissue summary:\n")

print(tissue_summary)

cat("\nCatalogue-integrity checks:\n")

print(catalogue_integrity)

cat("\nOutputs:\n")

cat(
  "  File-level QC:\n    ",
  file_qc_output,
  "\n"
)

cat(
  "  Original flag reasons:\n    ",
  flag_reason_output,
  "\n"
)

cat(
  "  Row-level flag summary:\n    ",
  row_flag_output,
  "\n"
)

cat(
  "  Tissue summary:\n    ",
  tissue_summary_output,
  "\n"
)

cat(
  "  Column consistency:\n    ",
  column_consistency_output,
  "\n"
)

cat(
  "  Catalogue integrity:\n    ",
  catalogue_integrity_output,
  "\n"
)

cat(
  "  Final QC status:\n    ",
  final_status_output,
  "\n"
)

cat(
  "  Log:\n    ",
  log_file,
  "\n\n"
)

cat("============================================================\n")

if (identical(final_status, "PASS")) {

  cat("FINAL QC RESULT: PASS\n")
  cat("The harmonized apaQTL catalogue passed all final checks.\n")
  cat("Module 11.3B is complete.\n")

} else if (
  identical(
    final_status,
    "PASS_WITH_WARNINGS"
  )
) {

  cat("FINAL QC RESULT: PASS WITH WARNINGS\n")
  cat("No critical integrity problem was detected.\n")
  cat("Review the warning-level row flags before downstream analysis.\n")

} else {

  cat("FINAL QC RESULT: FAIL\n")
  cat("Critical catalogue-integrity or row-level problems were detected.\n")
  cat("Review the Module 11.3B.4 QC tables before continuing.\n")
}

cat("============================================================\n\n")


###############################################################################
# 19. CLOSE LOGGING SAFELY
###############################################################################

while (sink.number(type = "output") > 0L) {
  try(
    sink(type = "output"),
    silent = TRUE
  )
}

try(
  close(log_connection),
  silent = TRUE
)

logging_active <- FALSE


###############################################################################
# 20. STOP ONLY IF CRITICAL QC FAILED
###############################################################################

if (identical(final_status, "FAIL")) {

  stop(
    paste0(
      "Module 11.3B.4 failed. Review:\n",
      final_status_output
    ),
    call. = FALSE
  )
}

invisible(
  list(
    final_status = final_qc_status,
    file_qc = file_qc,
    flag_reason_summary = flag_reason_summary,
    row_flag_summary = row_flag_summary,
    tissue_summary = tissue_summary,
    catalogue_integrity = catalogue_integrity
  )
)