#!/usr/bin/env Rscript

###############################################################################
# PigmentationAtlas
# Module 12.2 – Prepare and Test a Single GWAS–apaQTL Colocalization Unit
#
# Purpose
#   Select one ELIGIBLE unit from the Module 12.1 work queue, harmonize GWAS
#   and apaQTL variants, prepare coloc.abf() inputs, run a single test
#   colocalization, and write detailed summary, SNP-level and QC outputs.
#
# Usage
#   Default: first ELIGIBLE unit
#     Rscript scripts/module12_2_prepare_and_test_single_coloc.R
#
#   By coloc unit ID:
#     Rscript scripts/module12_2_prepare_and_test_single_coloc.R COLOC_0000001
#
#   By eligible queue row number:
#     Rscript scripts/module12_2_prepare_and_test_single_coloc.R 25
#
# Important
#   - The analysis unit is PLINK CLUMP × tissue × phenotype_id × gene_id.
#   - This module is intended as a single-unit validation step before the
#     complete restart-safe Module 12.3 catalogue run.
#   - Alleles are harmonized to the GWAS effect allele.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  warn = 1
)

module_start_time <- Sys.time()

###############################################################################
# 1. CONFIGURATION
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

work_queue_dir <- file.path(
  module12_dir,
  "work_queue"
)

single_test_dir <- file.path(
  module12_dir,
  "single_test"
)

input_dir <- file.path(
  single_test_dir,
  "input"
)

results_dir <- file.path(
  single_test_dir,
  "results"
)

qc_dir <- file.path(
  single_test_dir,
  "qc"
)

log_dir <- file.path(
  single_test_dir,
  "logs"
)

dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

eligible_queue_file <- file.path(
  work_queue_dir,
  "module12_1_coloc_work_queue_eligible.tsv.gz"
)

module12_1_status_file <- file.path(
  work_queue_dir,
  "module12_1_status.tsv"
)

# coloc prior probabilities
prior_p1 <- 1e-4
prior_p2 <- 1e-4
prior_p12 <- 1e-5

# Minimum requirements after full allele harmonization
minimum_harmonized_variants <- 50L
minimum_unique_variants <- 50L

# Dataset type. Pigmentation GWAS is treated as quantitative by default.
# Change to "cc" only if the GWAS phenotype is a case-control trait and
# suitable case fraction is available.
gwas_type <- "cc"
apaqtl_type <- "quant"

# Optional sample sizes. If NA, Module 12.2 will try to detect them from files.
# For coloc with beta/varbeta, N is strongly recommended but not mandatory.
default_gwas_sample_size <- 419469
default_apaqtl_sample_size <- NA_real_

# Case fraction required only when gwas_type == "cc"
default_gwas_case_fraction <- 0.3737566

###############################################################################
# 2. HELPER FUNCTIONS
###############################################################################

stop_module <- function(message_text) {
  cat("\nFATAL ERROR:\n", message_text, "\n", sep = "")
  stop(message_text, call. = FALSE)
}

first_matching_column <- function(column_names, candidates) {
  hit <- candidates[candidates %in% column_names]

  if (length(hit) == 0L) {
    return(NA_character_)
  }

  hit[[1L]]
}

normalize_chr <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x <- toupper(x)
  x[x == "23"] <- "X"
  x[x == "24"] <- "Y"
  x[x %in% c("M", "MTDNA")] <- "MT"
  x
}

normalize_position <- function(x) {
  suppressWarnings(as.integer(as.character(x)))
}

normalize_allele <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[!nzchar(x)] <- NA_character_
  x
}

clean_identifier <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(trimws(x))] <- NA_character_
  x
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

is_valid_p <- function(x) {
  x <- safe_numeric(x)
  !is.na(x) & is.finite(x) & x > 0 & x <= 1
}

is_valid_frequency <- function(x) {
  x <- safe_numeric(x)
  !is.na(x) & is.finite(x) & x > 0 & x < 1
}

is_valid_se <- function(x) {
  x <- safe_numeric(x)
  !is.na(x) & is.finite(x) & x > 0
}

complement_allele <- function(x) {
  x <- normalize_allele(x)

  result <- rep(NA_character_, length(x))
  result[x == "A"] <- "T"
  result[x == "T"] <- "A"
  result[x == "C"] <- "G"
  result[x == "G"] <- "C"

  result
}

is_palindromic_pair <- function(a1, a2) {
  a1 <- normalize_allele(a1)
  a2 <- normalize_allele(a2)

  pair <- paste(
    pmin(a1, a2),
    pmax(a1, a2),
    sep = "/"
  )

  pair %in% c("A/T", "C/G")
}

build_position_key <- function(chr, pos) {
  chr <- normalize_chr(chr)
  pos <- normalize_position(pos)

  valid <- !is.na(chr) & nzchar(chr) & !is.na(pos)

  result <- rep(NA_character_, length(chr))
  result[valid] <- paste(chr[valid], pos[valid], sep = ":")
  result
}

detect_gwas_columns <- function(column_names) {
  list(
    chr = first_matching_column(
      column_names,
      c("chromosome", "chr", "canonical_chromosome", "CHR")
    ),
    position = first_matching_column(
      column_names,
      c("base_pair_location", "position", "pos", "canonical_position", "BP")
    ),
    effect = first_matching_column(
      column_names,
      c("effect_allele", "alt_allele", "A1")
    ),
    other = first_matching_column(
      column_names,
      c("other_allele", "ref_allele", "A2")
    ),
    beta = first_matching_column(
      column_names,
      c("beta", "effect_size", "BETA")
    ),
    se = first_matching_column(
      column_names,
      c("standard_error", "se", "SE")
    ),
    p = first_matching_column(
      column_names,
      c("p_value", "pvalue", "pval", "P")
    ),
    af = first_matching_column(
      column_names,
      c(
        "effect_allele_frequency",
        "allele_frequency",
        "eaf",
        "EAF",
        "af"
      )
    ),
    maf = first_matching_column(
      column_names,
      c("maf", "minor_allele_frequency")
    ),
    rsid = first_matching_column(
      column_names,
      c("rsid", "rs_id", "SNP")
    ),
    variant = first_matching_column(
      column_names,
      c("variant_id", "variant", "ID")
    ),
    n = first_matching_column(
      column_names,
      c("n", "N", "sample_size", "samplesize")
    )
  )
}

detect_apaqtl_columns <- function(column_names) {
  list(
    phenotype = first_matching_column(
      column_names,
      c("phenotype_id", "molecular_trait_id", "phenotype")
    ),
    gene = first_matching_column(
      column_names,
      c("gene_id", "gene", "target_gene_id")
    ),
    chr = first_matching_column(
      column_names,
      c(
        "chromosome",
        "chr",
        "variant_chromosome",
        "canonical_chromosome"
      )
    ),
    position = first_matching_column(
      column_names,
      c(
        "position",
        "base_pair_location",
        "pos",
        "variant_position",
        "canonical_position"
      )
    ),
    ref = first_matching_column(
      column_names,
      c(
        "ref_allele",
        "reference_allele",
        "other_allele",
        "allele0"
      )
    ),
    alt = first_matching_column(
      column_names,
      c(
        "alt_allele",
        "alternate_allele",
        "effect_allele",
        "allele1"
      )
    ),
    beta = first_matching_column(
      column_names,
      c("beta", "slope", "effect_size")
    ),
    se = first_matching_column(
      column_names,
      c("standard_error", "se", "slope_se")
    ),
    p = first_matching_column(
      column_names,
      c("p_value", "pval_nominal", "pvalue", "pval")
    ),
    af = first_matching_column(
      column_names,
      c(
        "allele_frequency",
        "effect_allele_frequency",
        "af"
      )
    ),
    maf = first_matching_column(
      column_names,
      c("maf", "minor_allele_frequency")
    ),
    rsid = first_matching_column(
      column_names,
      c("rsid", "rs_id", "SNP")
    ),
    variant = first_matching_column(
      column_names,
      c("variant_id", "variant", "snp_id")
    ),
    n = first_matching_column(
      column_names,
      c("n", "N", "sample_size", "samplesize")
    )
  )
}

derive_maf <- function(af, maf) {
  af <- safe_numeric(af)
  maf <- safe_numeric(maf)

  result <- maf

  missing_maf <- is.na(result) & is_valid_frequency(af)
  result[missing_maf] <- pmin(
    af[missing_maf],
    1 - af[missing_maf]
  )

  result[
    !is.na(result) &
      (!is.finite(result) | result <= 0 | result > 0.5)
  ] <- NA_real_

  result
}

resolve_apaqtl_sample_size <- function(tissue_name) {
  normalized_tissue <- gsub(
    "[^a-z0-9]+",
    "_",
    tolower(as.character(tissue_name))
  )

  if (grepl("not_sun_exposed.*suprapubic", normalized_tissue)) {
    return(517)
  }

  if (grepl("sun_exposed.*lower_leg", normalized_tissue)) {
    return(605)
  }

  NA_real_
}

resolve_sample_size <- function(
  input_table,
  sample_size_column = NA_character_,
  default_sample_size = NA_real_
) {
  detected_values <- numeric(0)

  if (
    !is.na(sample_size_column) &&
      sample_size_column %in% names(input_table)
  ) {
    detected_values <- safe_numeric(
      input_table[[sample_size_column]]
    )

    detected_values <- detected_values[
      !is.na(detected_values) &
        is.finite(detected_values) &
        detected_values > 0
    ]
  }

  if (length(detected_values) > 0L) {
    # A scalar N is required by coloc. The median is robust when per-variant
    # sample sizes differ slightly because of variant-level missingness.
    return(as.numeric(round(stats::median(detected_values))))
  }

  default_sample_size <- safe_numeric(default_sample_size)

  if (
    length(default_sample_size) >= 1L &&
      !is.na(default_sample_size[[1L]]) &&
      is.finite(default_sample_size[[1L]]) &&
      default_sample_size[[1L]] > 0
  ) {
    return(as.numeric(default_sample_size[[1L]]))
  }

  NA_real_
}

select_unit <- function(queue, selector) {
  if (nrow(queue) == 0L) {
    stop_module("The Module 12.1 eligible queue is empty.")
  }

  if (is.null(selector) || !nzchar(selector)) {
    return(queue[1L])
  }

  if (grepl("^COLOC_[0-9]+$", selector)) {
    hit <- queue[coloc_unit_id == selector]

    if (nrow(hit) != 1L) {
      stop_module(
        paste0(
          "coloc_unit_id not found or not unique: ",
          selector
        )
      )
    }

    return(hit)
  }

  selector_index <- suppressWarnings(as.integer(selector))

  if (
    is.na(selector_index) ||
      selector_index < 1L ||
      selector_index > nrow(queue)
  ) {
    stop_module(
      paste0(
        "Invalid selector. Use a valid COLOC ID or an eligible queue row ",
        "between 1 and ",
        nrow(queue),
        "."
      )
    )
  }

  queue[selector_index]
}

###############################################################################
# 3. VALIDATE INPUTS
###############################################################################

if (!file.exists(eligible_queue_file)) {
  stop_module(
    paste0(
      "Module 12.1 eligible work queue does not exist:\n",
      eligible_queue_file
    )
  )
}

if (!file.exists(module12_1_status_file)) {
  stop_module(
    paste0(
      "Module 12.1 status file does not exist:\n",
      module12_1_status_file
    )
  )
}

status_table <- fread(
  module12_1_status_file,
  sep = "\t",
  header = TRUE,
  showProgress = FALSE
)

status_column <- first_matching_column(
  names(status_table),
  c("final_status", "status")
)

if (is.na(status_column)) {
  stop_module("Module 12.1 status file lacks a status column.")
}

module12_1_status <- as.character(
  status_table[[status_column]][1L]
)

if (!module12_1_status %in% c("PASS", "PASS_WITH_WARNINGS")) {
  stop_module(
    paste0(
      "Module 12.1 is not approved for downstream analysis. Status: ",
      module12_1_status
    )
  )
}

eligible_queue <- fread(
  eligible_queue_file,
  sep = "\t",
  header = TRUE,
  showProgress = FALSE
)

required_queue_columns <- c(
  "coloc_unit_id",
  "clump_id",
  "tissue",
  "phenotype_id",
  "gene_id",
  "apaqtl_chunk_file",
  "gwas_file",
  "queue_status"
)

missing_queue_columns <- setdiff(
  required_queue_columns,
  names(eligible_queue)
)

if (length(missing_queue_columns) > 0L) {
  stop_module(
    paste0(
      "Eligible queue is missing required columns: ",
      paste(missing_queue_columns, collapse = ", ")
    )
  )
}

eligible_queue <- eligible_queue[
  queue_status == "ELIGIBLE"
]

selector <- commandArgs(trailingOnly = TRUE)

selector <- if (length(selector) == 0L) {
  NULL
} else {
  selector[[1L]]
}

selected_unit <- select_unit(
  eligible_queue,
  selector
)

coloc_unit_id <- selected_unit$coloc_unit_id[[1L]]
clump_id <- selected_unit$clump_id[[1L]]
tissue <- selected_unit$tissue[[1L]]
phenotype_id <- selected_unit$phenotype_id[[1L]]
gene_id <- selected_unit$gene_id[[1L]]
apaqtl_file <- selected_unit$apaqtl_chunk_file[[1L]]
gwas_file <- selected_unit$gwas_file[[1L]]

default_apaqtl_sample_size <- resolve_apaqtl_sample_size(tissue)

if (is.na(default_apaqtl_sample_size)) {
  stop_module(
    paste0(
      "No apaQTL sample size mapping is defined for tissue: ",
      tissue
    )
  )
}

if (!file.exists(apaqtl_file)) {
  stop_module(
    paste0(
      "apaQTL file not found:\n",
      apaqtl_file
    )
  )
}

if (!file.exists(gwas_file)) {
  stop_module(
    paste0(
      "GWAS file not found:\n",
      gwas_file
    )
  )
}

unit_output_dir <- file.path(
  results_dir,
  paste0(
    coloc_unit_id,
    "_",
    clump_id
  )
)

dir.create(
  unit_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

prepared_input_output <- file.path(
  input_dir,
  paste0(
    coloc_unit_id,
    "_harmonized_input.tsv.gz"
  )
)

summary_output <- file.path(
  unit_output_dir,
  paste0(
    coloc_unit_id,
    "_coloc_summary.tsv"
  )
)

snp_output <- file.path(
  unit_output_dir,
  paste0(
    coloc_unit_id,
    "_coloc_snp_results.tsv.gz"
  )
)

qc_output <- file.path(
  qc_dir,
  paste0(
    coloc_unit_id,
    "_qc.tsv"
  )
)

status_output <- file.path(
  unit_output_dir,
  paste0(
    coloc_unit_id,
    "_status.tsv"
  )
)

log_file <- file.path(
  log_dir,
  paste0(
    "module12_2_",
    coloc_unit_id,
    "_",
    format(module_start_time, "%Y%m%d_%H%M%S"),
    ".log"
  )
)

session_info_output <- file.path(
  log_dir,
  paste0(
    "module12_2_",
    coloc_unit_id,
    "_session_info.txt"
  )
)

log_connection <- file(
  log_file,
  open = "wt",
  encoding = "UTF-8"
)

sink(log_connection, type = "output", split = TRUE)
sink(log_connection, type = "message", append = TRUE)

on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(type = "output"), silent = TRUE)
  try(close(log_connection), silent = TRUE)
}, add = TRUE)

cat("============================================================\n")
cat("PigmentationAtlas\n")
cat("MODULE 12.2 – SINGLE COLOCALIZATION TEST\n")
cat("============================================================\n\n")

cat("coloc_unit_id: ", coloc_unit_id, "\n", sep = "")
cat("CLUMP:         ", clump_id, "\n", sep = "")
cat("Tissue:        ", tissue, "\n", sep = "")
cat("Phenotype:     ", phenotype_id, "\n", sep = "")
cat("Gene:          ", gene_id, "\n\n", sep = "")

###############################################################################
# 4. READ AND STANDARDIZE GWAS
###############################################################################

gwas_raw <- fread(
  gwas_file,
  sep = "\t",
  header = TRUE,
  showProgress = FALSE
)

gwas_columns <- detect_gwas_columns(
  names(gwas_raw)
)

required_gwas_columns <- c(
  gwas_columns$chr,
  gwas_columns$position,
  gwas_columns$effect,
  gwas_columns$other,
  gwas_columns$beta,
  gwas_columns$se,
  gwas_columns$p
)

if (any(is.na(required_gwas_columns))) {
  stop_module(
    paste0(
      "GWAS file lacks one or more required columns. Detected mapping:\n",
      paste(
        names(gwas_columns),
        unlist(gwas_columns),
        sep = "=",
        collapse = "\n"
      )
    )
  )
}

gwas <- data.table(
  chromosome = normalize_chr(
    gwas_raw[[gwas_columns$chr]]
  ),
  position = normalize_position(
    gwas_raw[[gwas_columns$position]]
  ),
  gwas_effect_allele = normalize_allele(
    gwas_raw[[gwas_columns$effect]]
  ),
  gwas_other_allele = normalize_allele(
    gwas_raw[[gwas_columns$other]]
  ),
  gwas_beta = safe_numeric(
    gwas_raw[[gwas_columns$beta]]
  ),
  gwas_se = safe_numeric(
    gwas_raw[[gwas_columns$se]]
  ),
  gwas_p_value = safe_numeric(
    gwas_raw[[gwas_columns$p]]
  )
)

gwas[
  ,
  gwas_rsid := if (!is.na(gwas_columns$rsid)) {
    clean_identifier(gwas_raw[[gwas_columns$rsid]])
  } else {
    NA_character_
  }
]

gwas[
  ,
  gwas_variant_id := if (!is.na(gwas_columns$variant)) {
    clean_identifier(gwas_raw[[gwas_columns$variant]])
  } else {
    NA_character_
  }
]

gwas[
  ,
  gwas_af := if (!is.na(gwas_columns$af)) {
    safe_numeric(gwas_raw[[gwas_columns$af]])
  } else {
    NA_real_
  }
]

gwas[
  ,
  gwas_maf := derive_maf(
    gwas_af,
    if (!is.na(gwas_columns$maf)) {
      safe_numeric(gwas_raw[[gwas_columns$maf]])
    } else {
      NA_real_
    }
  )
]

gwas[
  ,
  position_key := build_position_key(
    chromosome,
    position
  )
]

gwas[
  ,
  gwas_varbeta := gwas_se^2
]

gwas_valid <- gwas[
  !is.na(position_key) &
    !is.na(gwas_effect_allele) &
    !is.na(gwas_other_allele) &
    is.finite(gwas_beta) &
    is_valid_se(gwas_se) &
    is_valid_p(gwas_p_value)
]

# Prefer the most significant valid row for duplicated positions.
setorder(
  gwas_valid,
  position_key,
  gwas_p_value
)

gwas_position_counts <- gwas_valid[, .N, by = position_key]
gwas_duplicate_positions <- gwas_position_counts[N > 1L, .N]
gwas_duplicate_rows <- gwas_position_counts[N > 1L, sum(N - 1L)]

gwas_valid <- gwas_valid[
  ,
  .SD[1L],
  by = position_key
]

###############################################################################
# 5. READ AND STANDARDIZE apaQTL
###############################################################################

apaqtl_raw <- fread(
  apaqtl_file,
  sep = "\t",
  header = TRUE,
  showProgress = FALSE
)

apaqtl_columns <- detect_apaqtl_columns(
  names(apaqtl_raw)
)

required_apaqtl_columns <- c(
  apaqtl_columns$phenotype,
  apaqtl_columns$gene,
  apaqtl_columns$chr,
  apaqtl_columns$position,
  apaqtl_columns$ref,
  apaqtl_columns$alt,
  apaqtl_columns$beta,
  apaqtl_columns$se,
  apaqtl_columns$p
)

if (any(is.na(required_apaqtl_columns))) {
  stop_module(
    paste0(
      "apaQTL file lacks one or more required columns. Detected mapping:\n",
      paste(
        names(apaqtl_columns),
        unlist(apaqtl_columns),
        sep = "=",
        collapse = "\n"
      )
    )
  )
}

apaqtl_phenotype_values <- clean_identifier(
  apaqtl_raw[[apaqtl_columns$phenotype]]
)

apaqtl_gene_values <- clean_identifier(
  apaqtl_raw[[apaqtl_columns$gene]]
)

apaqtl_keep_rows <-
  apaqtl_phenotype_values == phenotype_id &
  apaqtl_gene_values == gene_id

apaqtl_keep_rows[is.na(apaqtl_keep_rows)] <- FALSE

apaqtl_subset <- apaqtl_raw[apaqtl_keep_rows]

if (nrow(apaqtl_subset) == 0L) {
  stop_module(
    paste0(
      "No rows matched phenotype_id=",
      phenotype_id,
      " and gene_id=",
      gene_id,
      " in the selected apaQTL chunk."
    )
  )
}

apaqtl <- data.table(
  chromosome = normalize_chr(
    apaqtl_subset[[apaqtl_columns$chr]]
  ),
  position = normalize_position(
    apaqtl_subset[[apaqtl_columns$position]]
  ),
  apaqtl_ref_allele = normalize_allele(
    apaqtl_subset[[apaqtl_columns$ref]]
  ),
  apaqtl_alt_allele = normalize_allele(
    apaqtl_subset[[apaqtl_columns$alt]]
  ),
  apaqtl_beta = safe_numeric(
    apaqtl_subset[[apaqtl_columns$beta]]
  ),
  apaqtl_se = safe_numeric(
    apaqtl_subset[[apaqtl_columns$se]]
  ),
  apaqtl_p_value = safe_numeric(
    apaqtl_subset[[apaqtl_columns$p]]
  )
)

apaqtl[
  ,
  apaqtl_rsid := if (!is.na(apaqtl_columns$rsid)) {
    clean_identifier(apaqtl_subset[[apaqtl_columns$rsid]])
  } else {
    NA_character_
  }
]

apaqtl[
  ,
  apaqtl_variant_id := if (!is.na(apaqtl_columns$variant)) {
    clean_identifier(apaqtl_subset[[apaqtl_columns$variant]])
  } else {
    NA_character_
  }
]

apaqtl[
  ,
  apaqtl_af := if (!is.na(apaqtl_columns$af)) {
    safe_numeric(apaqtl_subset[[apaqtl_columns$af]])
  } else {
    NA_real_
  }
]

apaqtl[
  ,
  apaqtl_maf := derive_maf(
    apaqtl_af,
    if (!is.na(apaqtl_columns$maf)) {
      safe_numeric(apaqtl_subset[[apaqtl_columns$maf]])
    } else {
      NA_real_
    }
  )
]

apaqtl[
  ,
  position_key := build_position_key(
    chromosome,
    position
  )
]

apaqtl[
  ,
  apaqtl_varbeta := apaqtl_se^2
]

apaqtl_valid <- apaqtl[
  !is.na(position_key) &
    !is.na(apaqtl_ref_allele) &
    !is.na(apaqtl_alt_allele) &
    is.finite(apaqtl_beta) &
    is_valid_se(apaqtl_se) &
    is_valid_p(apaqtl_p_value)
]

setorder(
  apaqtl_valid,
  position_key,
  apaqtl_p_value
)

apaqtl_position_counts <- apaqtl_valid[, .N, by = position_key]
apaqtl_duplicate_positions <- apaqtl_position_counts[N > 1L, .N]
apaqtl_duplicate_rows <- apaqtl_position_counts[N > 1L, sum(N - 1L)]

apaqtl_valid <- apaqtl_valid[
  ,
  .SD[1L],
  by = position_key
]

###############################################################################
# 6. MERGE AND HARMONIZE ALLELES
###############################################################################

merged <- merge(
  gwas_valid,
  apaqtl_valid,
  by = c(
    "position_key",
    "chromosome",
    "position"
  ),
  all = FALSE,
  sort = FALSE
)

if (nrow(merged) == 0L) {
  stop_module(
    "No GWAS–apaQTL variants overlapped by chromosome and position."
  )
}

merged[
  ,
  allele_relationship := fcase(
    apaqtl_alt_allele == gwas_effect_allele &
      apaqtl_ref_allele == gwas_other_allele,
    "DIRECT",

    apaqtl_alt_allele == gwas_other_allele &
      apaqtl_ref_allele == gwas_effect_allele,
    "SWAPPED",

    complement_allele(apaqtl_alt_allele) == gwas_effect_allele &
      complement_allele(apaqtl_ref_allele) == gwas_other_allele,
    "COMPLEMENT_DIRECT",

    complement_allele(apaqtl_alt_allele) == gwas_other_allele &
      complement_allele(apaqtl_ref_allele) == gwas_effect_allele,
    "COMPLEMENT_SWAPPED",

    default = "INCOMPATIBLE"
  )
]

merged[
  ,
  palindromic := is_palindromic_pair(
    gwas_effect_allele,
    gwas_other_allele
  )
]

# Palindromic variants are retained only when both allele frequencies exist
# and clearly support orientation. Otherwise they are excluded conservatively.
merged[
  ,
  palindromic_frequency_resolvable := (
    palindromic == FALSE |
      (
        is_valid_frequency(gwas_af) &
          is_valid_frequency(apaqtl_af) &
          (
            abs(gwas_af - apaqtl_af) <= 0.10 |
              abs(gwas_af - (1 - apaqtl_af)) <= 0.10
          )
      )
  )
]

merged[
  ,
  harmonization_status := fcase(
    allele_relationship == "INCOMPATIBLE",
    "EXCLUDE_INCOMPATIBLE_ALLELES",

    palindromic == TRUE &
      palindromic_frequency_resolvable == FALSE,
    "EXCLUDE_AMBIGUOUS_PALINDROMIC",

    allele_relationship %in% c(
      "DIRECT",
      "COMPLEMENT_DIRECT"
    ),
    "KEEP_DIRECT",

    allele_relationship %in% c(
      "SWAPPED",
      "COMPLEMENT_SWAPPED"
    ),
    "KEEP_FLIPPED",

    default = "EXCLUDE_UNRESOLVED"
  )
]

merged[
  ,
  apaqtl_beta_harmonized := fcase(
    harmonization_status == "KEEP_DIRECT",
    apaqtl_beta,

    harmonization_status == "KEEP_FLIPPED",
    -apaqtl_beta,

    default = NA_real_
  )
]

merged[
  ,
  apaqtl_af_harmonized := fcase(
    harmonization_status == "KEEP_DIRECT",
    apaqtl_af,

    harmonization_status == "KEEP_FLIPPED" &
      is_valid_frequency(apaqtl_af),
    1 - apaqtl_af,

    default = NA_real_
  )
]

harmonized <- merged[
  harmonization_status %in% c(
    "KEEP_DIRECT",
    "KEEP_FLIPPED"
  )
]

if (nrow(harmonized) < minimum_harmonized_variants) {
  stop_module(
    paste0(
      "Only ",
      nrow(harmonized),
      " variants remained after allele harmonization; minimum required is ",
      minimum_harmonized_variants,
      "."
    )
  )
}

if (uniqueN(harmonized$position_key) < minimum_unique_variants) {
  stop_module(
    paste0(
      "Only ",
      uniqueN(harmonized$position_key),
      " unique variants remained after harmonization; minimum required is ",
      minimum_unique_variants,
      "."
    )
  )
}

###############################################################################
# 7. RESOLVE SAMPLE SIZES AND MAF
###############################################################################

gwas_n <- resolve_sample_size(
  gwas_raw,
  gwas_columns$n,
  default_gwas_sample_size
)

apaqtl_n <- resolve_sample_size(
  apaqtl_subset,
  apaqtl_columns$n,
  default_apaqtl_sample_size
)

harmonized[
  ,
  coloc_maf_gwas := derive_maf(
    gwas_af,
    gwas_maf
  )
]

harmonized[
  ,
  coloc_maf_apaqtl := derive_maf(
    apaqtl_af_harmonized,
    apaqtl_maf
  )
]

# coloc requires MAF only for certain input configurations. Since beta and
# varbeta are provided, rows without MAF can still be analyzed.
variant_names <- ifelse(
  !is.na(harmonized$gwas_rsid) &
    grepl("^rs[0-9]+$", harmonized$gwas_rsid),
  harmonized$gwas_rsid,
  harmonized$position_key
)

variant_names <- make.unique(
  variant_names,
  sep = "_dup"
)

harmonized[, coloc_snp := variant_names]

###############################################################################
# 8. PREPARE COLOC DATASETS
###############################################################################

dataset_gwas <- list(
  beta = harmonized$gwas_beta,
  varbeta = harmonized$gwas_se^2,
  snp = harmonized$coloc_snp,
  position = harmonized$position,
  type = "cc",
  N = gwas_n,
  s = default_gwas_case_fraction,
  MAF = harmonized$coloc_maf_gwas
)

dataset_apaqtl <- list(
  beta = harmonized$apaqtl_beta_harmonized,
  varbeta = harmonized$apaqtl_se^2,
  snp = harmonized$coloc_snp,
  position = harmonized$position,
  type = "quant",
  N = apaqtl_n,
  MAF = harmonized$coloc_maf_apaqtl
)

if (any(!is.na(harmonized$coloc_maf_gwas))) {
  dataset_gwas$MAF <- harmonized$coloc_maf_gwas
}

if (any(!is.na(harmonized$coloc_maf_apaqtl))) {
  dataset_apaqtl$MAF <- harmonized$coloc_maf_apaqtl
}

if (!is.na(gwas_n)) {
  dataset_gwas$N <- gwas_n
}

if (!is.na(apaqtl_n)) {
  dataset_apaqtl$N <- apaqtl_n
}

if (identical(gwas_type, "cc")) {
  if (
    is.na(default_gwas_case_fraction) ||
      default_gwas_case_fraction <= 0 ||
      default_gwas_case_fraction >= 1
  ) {
    stop_module(
      "gwas_type is 'cc', but no valid case fraction was configured."
    )
  }

  dataset_gwas$s <- default_gwas_case_fraction
}

check_gwas <- check_dataset(
  dataset_gwas,
  suffix = "GWAS"
)

check_apaqtl <- check_dataset(
  dataset_apaqtl,
  suffix = "apaQTL"
)

###############################################################################
# 9. RUN COLOCALIZATION
###############################################################################

coloc_result <- coloc.abf(
  dataset1 = dataset_gwas,
  dataset2 = dataset_apaqtl,
  p1 = prior_p1,
  p2 = prior_p2,
  p12 = prior_p12
)

posterior_summary <- as.list(
  coloc_result$summary
)

summary_table <- data.table(
  module = "12.2",
  coloc_unit_id = coloc_unit_id,
  clump_id = clump_id,
  tissue = tissue,
  phenotype_id = phenotype_id,
  gene_id = gene_id,
  gwas_file = normalizePath(
    gwas_file,
    winslash = "/",
    mustWork = TRUE
  ),
  apaqtl_file = normalizePath(
    apaqtl_file,
    winslash = "/",
    mustWork = TRUE
  ),
  n_gwas_valid = nrow(gwas_valid),
  n_apaqtl_valid = nrow(apaqtl_valid),
  n_position_overlaps = nrow(merged),
  n_harmonized_variants = nrow(harmonized),
  gwas_sample_size = gwas_n,
  apaqtl_sample_size = apaqtl_n,
  gwas_type = gwas_type,
  apaqtl_type = apaqtl_type,
  prior_p1 = prior_p1,
  prior_p2 = prior_p2,
  prior_p12 = prior_p12,
  nsnps = as.numeric(posterior_summary[["nsnps"]]),
  PP.H0 = as.numeric(posterior_summary[["PP.H0.abf"]]),
  PP.H1 = as.numeric(posterior_summary[["PP.H1.abf"]]),
  PP.H2 = as.numeric(posterior_summary[["PP.H2.abf"]]),
  PP.H3 = as.numeric(posterior_summary[["PP.H3.abf"]]),
  PP.H4 = as.numeric(posterior_summary[["PP.H4.abf"]])
)

summary_table[
  ,
  PP.H4_over_H3_H4 := ifelse(
    (PP.H3 + PP.H4) > 0,
    PP.H4 / (PP.H3 + PP.H4),
    NA_real_
  )
]

summary_table[
  ,
  posterior_sum := PP.H0 + PP.H1 + PP.H2 + PP.H3 + PP.H4
]

summary_table[
  ,
  interpretation := fcase(
    PP.H4 >= 0.80,
    "STRONG_COLOCALIZATION",

    PP.H4 >= 0.50,
    "MODERATE_COLOCALIZATION",

    PP.H3 >= 0.80,
    "DISTINCT_ASSOCIATION_SIGNALS",

    PP.H2 >= 0.80,
    "APAQTL_SIGNAL_ONLY",

    PP.H1 >= 0.80,
    "GWAS_SIGNAL_ONLY",

    PP.H0 >= 0.80,
    "NO_ASSOCIATION_SIGNAL",

    default = "INCONCLUSIVE"
  )
]

snp_results <- as.data.table(
  coloc_result$results
)

setnames(
  snp_results,
  old = intersect(
    c("snp", "position"),
    names(snp_results)
  ),
  new = intersect(
    c("coloc_snp", "coloc_position"),
    c("coloc_snp", "coloc_position")
  )
)

snp_results <- merge(
  harmonized,
  snp_results,
  by = "coloc_snp",
  all.x = TRUE,
  sort = FALSE
)

setorder(
  snp_results,
  -SNP.PP.H4
)

summary_table[
  ,
  top_coloc_variant := if (
    nrow(snp_results) > 0L &&
      "SNP.PP.H4" %in% names(snp_results)
  ) {
    snp_results$coloc_snp[[1L]]
  } else {
    NA_character_
  }
]

summary_table[
  ,
  top_coloc_variant_PP.H4 := if (
    nrow(snp_results) > 0L &&
      "SNP.PP.H4" %in% names(snp_results)
  ) {
    snp_results$SNP.PP.H4[[1L]]
  } else {
    NA_real_
  }
]

###############################################################################
# 10. QC
###############################################################################

harmonization_distribution <- merged[
  ,
  .N,
  by = .(
    allele_relationship,
    harmonization_status
  )
][
  order(-N)
]

critical_qc_failures <- c(
  nrow(harmonized) < minimum_harmonized_variants,
  uniqueN(harmonized$position_key) < minimum_unique_variants,
  any(!is.finite(summary_table$posterior_sum)),
  abs(summary_table$posterior_sum - 1) > 1e-6,
  any(
    unlist(
      summary_table[
        ,
        .(PP.H0, PP.H1, PP.H2, PP.H3, PP.H4)
      ]
    ) < 0
  ),
  any(
    unlist(
      summary_table[
        ,
        .(PP.H0, PP.H1, PP.H2, PP.H3, PP.H4)
      ]
    ) > 1
  )
)

critical_failure_count <- sum(
  critical_qc_failures,
  na.rm = TRUE
)

warning_count <- sum(
  c(
    gwas_duplicate_positions > 0L,
    apaqtl_duplicate_positions > 0L,
    any(
      merged$harmonization_status ==
        "EXCLUDE_AMBIGUOUS_PALINDROMIC"
    ),
    any(
      merged$harmonization_status ==
        "EXCLUDE_INCOMPATIBLE_ALLELES"
    ),
    is.na(gwas_n),
    is.na(apaqtl_n),
    all(is.na(harmonized$coloc_maf_gwas)),
    all(is.na(harmonized$coloc_maf_apaqtl))
  ),
  na.rm = TRUE
)

final_status <- if (
  critical_failure_count > 0L
) {
  "FAIL"
} else if (
  warning_count > 0L
) {
  "PASS_WITH_WARNINGS"
} else {
  "PASS"
}

module_end_time <- Sys.time()

runtime_seconds <- as.numeric(
  difftime(
    module_end_time,
    module_start_time,
    units = "secs"
  )
)

qc_table <- data.table(
  module = "12.2",
  coloc_unit_id = coloc_unit_id,
  start_time = format(
    module_start_time,
    "%Y-%m-%d %H:%M:%S"
  ),
  end_time = format(
    module_end_time,
    "%Y-%m-%d %H:%M:%S"
  ),
  runtime_seconds = runtime_seconds,
  raw_gwas_rows = nrow(gwas_raw),
  valid_gwas_rows = nrow(gwas_valid),
  duplicated_gwas_positions = gwas_duplicate_positions,
  duplicated_gwas_rows_beyond_first = gwas_duplicate_rows,
  raw_apaqtl_rows_for_unit = nrow(apaqtl_subset),
  valid_apaqtl_rows = nrow(apaqtl_valid),
  duplicated_apaqtl_positions = apaqtl_duplicate_positions,
  duplicated_apaqtl_rows_beyond_first = apaqtl_duplicate_rows,
  position_overlap_rows = nrow(merged),
  direct_rows = sum(
    merged$allele_relationship == "DIRECT"
  ),
  swapped_rows = sum(
    merged$allele_relationship == "SWAPPED"
  ),
  complement_direct_rows = sum(
    merged$allele_relationship == "COMPLEMENT_DIRECT"
  ),
  complement_swapped_rows = sum(
    merged$allele_relationship == "COMPLEMENT_SWAPPED"
  ),
  incompatible_allele_rows = sum(
    merged$allele_relationship == "INCOMPATIBLE"
  ),
  ambiguous_palindromic_rows = sum(
    merged$harmonization_status ==
      "EXCLUDE_AMBIGUOUS_PALINDROMIC"
  ),
  harmonized_rows = nrow(harmonized),
  gwas_sample_size_available = !is.na(gwas_n),
  apaqtl_sample_size_available = !is.na(apaqtl_n),
  gwas_maf_available_rows = sum(
    !is.na(harmonized$coloc_maf_gwas)
  ),
  apaqtl_maf_available_rows = sum(
    !is.na(harmonized$coloc_maf_apaqtl)
  ),
  posterior_sum = summary_table$posterior_sum,
  critical_failures = critical_failure_count,
  warnings = warning_count,
  final_status = final_status
)

status_table <- data.table(
  module = "12.2",
  coloc_unit_id = coloc_unit_id,
  clump_id = clump_id,
  tissue = tissue,
  phenotype_id = phenotype_id,
  gene_id = gene_id,
  final_status = final_status,
  n_harmonized_variants = nrow(harmonized),
  PP.H4 = summary_table$PP.H4,
  interpretation = summary_table$interpretation,
  summary_file = normalizePath(
    summary_output,
    winslash = "/",
    mustWork = FALSE
  ),
  snp_results_file = normalizePath(
    snp_output,
    winslash = "/",
    mustWork = FALSE
  ),
  created_at = format(
    module_end_time,
    "%Y-%m-%d %H:%M:%S"
  )
)

###############################################################################
# 11. WRITE OUTPUTS
###############################################################################

fwrite(
  harmonized,
  prepared_input_output,
  sep = "\t",
  quote = FALSE,
  na = "NA",
  compress = "gzip"
)

fwrite(
  summary_table,
  summary_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  snp_results,
  snp_output,
  sep = "\t",
  quote = FALSE,
  na = "NA",
  compress = "gzip"
)

fwrite(
  qc_table,
  qc_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  status_table,
  status_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

harmonization_distribution_output <- file.path(
  qc_dir,
  paste0(
    coloc_unit_id,
    "_harmonization_distribution.tsv"
  )
)

fwrite(
  harmonization_distribution,
  harmonization_distribution_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

writeLines(
  capture.output(sessionInfo()),
  session_info_output,
  useBytes = TRUE
)

###############################################################################
# 12. FINAL REPORT
###############################################################################

cat("\n============================================================\n")
cat("MODULE 12.2 – FINAL SUMMARY\n")
cat("============================================================\n\n")

print(summary_table)

cat("\nQC:\n")
print(qc_table)

cat("\nAllele harmonization distribution:\n")
print(harmonization_distribution)

cat("\nOutputs:\n")
cat("  Harmonized input:\n    ", prepared_input_output, "\n", sep = "")
cat("  Coloc summary:\n    ", summary_output, "\n", sep = "")
cat("  SNP-level results:\n    ", snp_output, "\n", sep = "")
cat("  QC:\n    ", qc_output, "\n", sep = "")
cat("  Status:\n    ", status_output, "\n", sep = "")
cat("  Log:\n    ", log_file, "\n\n", sep = "")

cat("============================================================\n")

if (identical(final_status, "PASS")) {
  cat("FINAL RESULT: PASS\n")
  cat("The selected colocalization unit completed successfully.\n")
  cat("Module 12.3 full catalogue development may proceed.\n")
} else if (identical(final_status, "PASS_WITH_WARNINGS")) {
  cat("FINAL RESULT: PASS WITH WARNINGS\n")
  cat("The selected colocalization unit completed successfully, but QC\n")
  cat("warnings should be reviewed before the full catalogue run.\n")
} else {
  cat("FINAL RESULT: FAIL\n")
  cat("The selected colocalization unit failed QC.\n")
}

cat("============================================================\n")

if (identical(final_status, "FAIL")) {
  stop(
    paste0(
      "Module 12.2 failed. Review:\n",
      qc_output
    ),
    call. = FALSE
  )
}

invisible(
  list(
    final_status = status_table,
    summary = summary_table,
    snp_results = snp_results,
    qc = qc_table,
    harmonized_input = harmonized
  )
)

