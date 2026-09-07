#!/usr/bin/env Rscript

###############################################################################
# PigmentationAtlas
# Module 12.1 – Initialize GWAS–apaQTL Colocalization Work Queue
#
# Purpose
#   Build a restart-safe, QC-audited work queue for downstream colocalization.
#   The analysis unit is:
#
#       PLINK CLUMP × tissue × phenotype_id × gene_id
#
# Inputs
#   1. Module 11.3B harmonized apaQTL file-level QC table
#   2. Module 11.3B work queue
#   3. Harmonized per-file apaQTL chunks
#   4. CLUMP-specific GWAS workspace files
#   5. Module 11.3B.4 final QC status
#
# Main outputs
#   results/module12/colocalization/work_queue/
#     module12_1_coloc_work_queue.tsv.gz
#     module12_1_coloc_work_queue_eligible.tsv.gz
#     module12_1_coloc_work_queue_excluded.tsv.gz
#     module12_1_file_mapping.tsv
#     module12_1_qc_summary.tsv
#     module12_1_status.tsv
#
# Notes
#   - Additional columns are permitted.
#   - Missing mandatory columns are fatal.
#   - rsID is not required; chromosome-position-allele keys are preferred.
#   - Allele orientation is not harmonized in this module. That belongs to
#     Module 12.2.
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
# 1. CONFIGURATION
###############################################################################

project_dir <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

module11_dir <- file.path(
  project_dir,
  "results",
  "module11",
  "harmonization"
)

catalogue_dir <- file.path(
  module11_dir,
  "catalogue"
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

qc_dir <- file.path(
  work_queue_dir,
  "qc"
)

log_dir <- file.path(
  module12_dir,
  "logs"
)

dir.create(work_queue_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

module11_work_queue_file <- file.path(
  module11_dir,
  "catalogue",
  "module11_3B_work_queue.tsv"
)

if (!file.exists(module11_work_queue_file)) {
  alternative_work_queue <- file.path(
    module11_dir,
    "module11_3B_work_queue.tsv"
  )

  if (file.exists(alternative_work_queue)) {
    module11_work_queue_file <- alternative_work_queue
  }
}

module11_file_qc_file <- file.path(
  catalogue_dir,
  "qc",
  "final_qc",
  "module11_3B_4_file_level_qc.tsv"
)

module11_final_status_file <- file.path(
  catalogue_dir,
  "module11_3B_4_final_qc_status.tsv"
)

clump_root_dir <- file.path(
  project_dir,
  "results",
  "clumps"
)

# Conservative default threshold. Units below this threshold are retained in
# the excluded table and can later be reconsidered.
minimum_shared_variants <- 50L

# A signal-level threshold is deliberately not imposed here. Module 12.2 will
# calculate appropriate p-value and variance checks after allele harmonization.
minimum_apaqtl_rows <- 1L
minimum_gwas_rows <- 1L

###############################################################################
# 2. OUTPUT PATHS
###############################################################################

all_queue_output <- file.path(
  work_queue_dir,
  "module12_1_coloc_work_queue.tsv.gz"
)

eligible_queue_output <- file.path(
  work_queue_dir,
  "module12_1_coloc_work_queue_eligible.tsv.gz"
)

excluded_queue_output <- file.path(
  work_queue_dir,
  "module12_1_coloc_work_queue_excluded.tsv.gz"
)

file_mapping_output <- file.path(
  work_queue_dir,
  "module12_1_file_mapping.tsv"
)

qc_summary_output <- file.path(
  qc_dir,
  "module12_1_qc_summary.tsv"
)

exclusion_reason_output <- file.path(
  qc_dir,
  "module12_1_exclusion_reason_summary.tsv"
)

tissue_summary_output <- file.path(
  qc_dir,
  "module12_1_tissue_summary.tsv"
)

clump_summary_output <- file.path(
  qc_dir,
  "module12_1_clump_summary.tsv"
)

status_output <- file.path(
  work_queue_dir,
  "module12_1_status.tsv"
)

session_info_output <- file.path(
  log_dir,
  "module12_1_session_info.txt"
)

log_file <- file.path(
  log_dir,
  paste0(
    "module12_1_initialize_coloc_work_queue_",
    format(module_start_time, "%Y%m%d_%H%M%S"),
    ".log"
  )
)

###############################################################################
# 3. LOGGING
###############################################################################

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
cat("MODULE 12.1 – INITIALIZE COLOCALIZATION WORK QUEUE\n")
cat("============================================================\n\n")

cat("Project directory:\n  ", project_dir, "\n\n", sep = "")
cat("Start time:\n  ", format(module_start_time), "\n\n", sep = "")

###############################################################################
# 4. HELPER FUNCTIONS
###############################################################################

stop_module <- function(message_text) {
  cat("\nFATAL ERROR:\n", message_text, "\n", sep = "")
  stop(message_text, call. = FALSE)
}

first_existing_file <- function(paths) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  paths <- paths[file.exists(paths)]

  if (length(paths) == 0L) {
    return(NA_character_)
  }

  normalizePath(
    paths[[1L]],
    winslash = "/",
    mustWork = TRUE
  )
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

normalize_allele <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[!nzchar(x)] <- NA_character_
  x
}

normalize_position <- function(x) {
  suppressWarnings(as.integer(as.character(x)))
}

clean_identifier <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(trimws(x))] <- NA_character_
  x
}

resolve_path <- function(x) {
  x <- clean_identifier(x)

  if (is.na(x)) {
    return(NA_character_)
  }

  candidate_paths <- c(
    x,
    file.path(project_dir, x)
  )

  first_existing_file(candidate_paths)
}

build_position_key <- function(chr, pos) {
  chr <- normalize_chr(chr)
  pos <- normalize_position(pos)

  valid <- !is.na(chr) & nzchar(chr) & !is.na(pos)

  result <- rep(NA_character_, length(chr))
  result[valid] <- paste(chr[valid], pos[valid], sep = ":")
  result
}

build_unordered_allele_key <- function(chr, pos, allele1, allele2) {
  chr <- normalize_chr(chr)
  pos <- normalize_position(pos)
  allele1 <- normalize_allele(allele1)
  allele2 <- normalize_allele(allele2)

  valid <- (
    !is.na(chr) &
      nzchar(chr) &
      !is.na(pos) &
      !is.na(allele1) &
      !is.na(allele2)
  )

  lower_allele <- ifelse(
    valid,
    pmin(allele1, allele2),
    NA_character_
  )

  upper_allele <- ifelse(
    valid,
    pmax(allele1, allele2),
    NA_character_
  )

  result <- rep(NA_character_, length(chr))

  result[valid] <- paste(
    chr[valid],
    pos[valid],
    lower_allele[valid],
    upper_allele[valid],
    sep = ":"
  )

  result
}

read_header <- function(path) {
  fread(
    path,
    sep = "\t",
    header = TRUE,
    nrows = 0L,
    showProgress = FALSE
  )
}

read_selected_columns <- function(path, columns) {
  fread(
    path,
    sep = "\t",
    header = TRUE,
    select = unique(columns),
    showProgress = FALSE
  )
}

find_gwas_file_for_clump <- function(clump_id, clump_dir = NA_character_) {
  search_dirs <- unique(
    c(
      clump_dir,
      file.path(clump_root_dir, clump_id),
      list.dirs(
        clump_root_dir,
        recursive = FALSE,
        full.names = TRUE
      )[
        grepl(
          paste0("^", clump_id, "(_|$)"),
          basename(
            list.dirs(
              clump_root_dir,
              recursive = FALSE,
              full.names = TRUE
            )
          )
        )
      ]
    )
  )

  search_dirs <- search_dirs[
    !is.na(search_dirs) &
      nzchar(search_dirs) &
      dir.exists(search_dirs)
  ]

  if (length(search_dirs) == 0L) {
    return(NA_character_)
  }

  candidate_files <- unique(
    unlist(
      lapply(
        search_dirs,
        function(search_dir) {
          list.files(
            search_dir,
            pattern = "\\.(tsv|txt)(\\.gz)?$",
            recursive = TRUE,
            full.names = TRUE,
            ignore.case = TRUE
          )
        }
      ),
      use.names = FALSE
    )
  )

  if (length(candidate_files) == 0L) {
    return(NA_character_)
  }

  base_names <- tolower(basename(candidate_files))

  exclude_pattern <- paste(
    c(
      "metadata",
      "candidate",
      "annotation",
      "gene",
      "qc",
      "summary",
      "coloc",
      "apaqtl",
      "finemap",
      "credible",
      "susie",
      "log"
    ),
    collapse = "|"
  )

  candidate_files <- candidate_files[
    !grepl(exclude_pattern, base_names)
  ]

  if (length(candidate_files) == 0L) {
    return(NA_character_)
  }

  base_names <- tolower(basename(candidate_files))

  score <- integer(length(candidate_files))
  score <- score + ifelse(grepl("gwas", base_names), 100L, 0L)
  score <- score + ifelse(grepl("region", base_names), 30L, 0L)
  score <- score + ifelse(grepl("clump", base_names), 20L, 0L)
  score <- score + ifelse(grepl("harmon", base_names), 10L, 0L)
  score <- score + ifelse(grepl("\\.gz$", base_names), 5L, 0L)

  candidate_files <- candidate_files[
    order(
      -score,
      nchar(candidate_files),
      candidate_files
    )
  ]

  for (candidate_file in candidate_files) {
    header_result <- try(
      names(read_header(candidate_file)),
      silent = TRUE
    )

    if (inherits(header_result, "try-error")) {
      next
    }

    chr_column <- first_matching_column(
      header_result,
      c(
        "chromosome",
        "chr",
        "canonical_chromosome",
        "CHR"
      )
    )

    position_column <- first_matching_column(
      header_result,
      c(
        "base_pair_location",
        "position",
        "pos",
        "canonical_position",
        "BP"
      )
    )

    p_column <- first_matching_column(
      header_result,
      c(
        "p_value",
        "pvalue",
        "pval",
        "P"
      )
    )

    if (
      !is.na(chr_column) &&
        !is.na(position_column) &&
        !is.na(p_column)
    ) {
      return(
        normalizePath(
          candidate_file,
          winslash = "/",
          mustWork = TRUE
        )
      )
    }
  }

  NA_character_
}

detect_apaqtl_columns <- function(column_names) {
  list(
    clump = first_matching_column(
      column_names,
      c("clump_id", "metadata_clump_id")
    ),
    tissue = first_matching_column(
      column_names,
      c("tissue", "metadata_tissue")
    ),
    phenotype = first_matching_column(
      column_names,
      c(
        "phenotype_id",
        "molecular_trait_id",
        "phenotype"
      )
    ),
    gene = first_matching_column(
      column_names,
      c(
        "gene_id",
        "gene",
        "target_gene_id"
      )
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
    variant = first_matching_column(
      column_names,
      c(
        "variant_id",
        "variant",
        "snp_id"
      )
    ),
    rsid = first_matching_column(
      column_names,
      c("rsid", "rs_id", "SNP")
    ),
    p = first_matching_column(
      column_names,
      c(
        "p_value",
        "pval_nominal",
        "pvalue",
        "pval"
      )
    ),
    beta = first_matching_column(
      column_names,
      c(
        "beta",
        "slope",
        "effect_size"
      )
    ),
    se = first_matching_column(
      column_names,
      c(
        "standard_error",
        "se",
        "slope_se"
      )
    ),
    maf = first_matching_column(
      column_names,
      c(
        "maf",
        "minor_allele_frequency"
      )
    ),
    af = first_matching_column(
      column_names,
      c(
        "allele_frequency",
        "effect_allele_frequency",
        "af"
      )
    )
  )
}

detect_gwas_columns <- function(column_names) {
  list(
    chr = first_matching_column(
      column_names,
      c(
        "chromosome",
        "chr",
        "canonical_chromosome",
        "CHR"
      )
    ),
    position = first_matching_column(
      column_names,
      c(
        "base_pair_location",
        "position",
        "pos",
        "canonical_position",
        "BP"
      )
    ),
    effect = first_matching_column(
      column_names,
      c(
        "effect_allele",
        "alt_allele",
        "A1"
      )
    ),
    other = first_matching_column(
      column_names,
      c(
        "other_allele",
        "ref_allele",
        "A2"
      )
    ),
    variant = first_matching_column(
      column_names,
      c(
        "variant_id",
        "variant",
        "ID"
      )
    ),
    rsid = first_matching_column(
      column_names,
      c(
        "rsid",
        "rs_id",
        "SNP"
      )
    ),
    p = first_matching_column(
      column_names,
      c(
        "p_value",
        "pvalue",
        "pval",
        "P"
      )
    ),
    beta = first_matching_column(
      column_names,
      c(
        "beta",
        "effect_size",
        "BETA"
      )
    ),
    se = first_matching_column(
      column_names,
      c(
        "standard_error",
        "se",
        "SE"
      )
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
    )
  )
}

###############################################################################
# 5. VALIDATE MODULE 11 INPUTS
###############################################################################

required_input_files <- c(
  module11_work_queue_file,
  module11_file_qc_file,
  module11_final_status_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0L) {
  stop_module(
    paste(
      "Missing required Module 11 input file(s):",
      paste(missing_input_files, collapse = "\n"),
      sep = "\n"
    )
  )
}

if (!dir.exists(clump_root_dir)) {
  stop_module(
    paste0(
      "CLUMP workspace directory does not exist:\n",
      clump_root_dir
    )
  )
}

final_status_table <- fread(
  module11_final_status_file,
  sep = "\t",
  header = TRUE,
  showProgress = FALSE
)

status_column <- first_matching_column(
  names(final_status_table),
  c("final_status", "status")
)

if (is.na(status_column)) {
  stop_module(
    "Module 11.3B.4 final status table has no final_status/status column."
  )
}

module11_final_status <- as.character(
  final_status_table[[status_column]][1L]
)

if (
  !module11_final_status %in%
    c("PASS", "PASS_WITH_WARNINGS")
) {
  stop_module(
    paste0(
      "Module 11.3B.4 is not approved for downstream analysis. ",
      "Observed status: ",
      module11_final_status
    )
  )
}

cat(
  "Module 11.3B.4 status: ",
  module11_final_status,
  "\n\n",
  sep = ""
)

###############################################################################
# 6. READ MODULE 11 WORK QUEUE AND FILE QC
###############################################################################

module11_work_queue <- fread(
  module11_work_queue_file,
  sep = "\t",
  header = TRUE,
  showProgress = FALSE
)

file_qc <- fread(
  module11_file_qc_file,
  sep = "\t",
  header = TRUE,
  showProgress = FALSE
)

if (nrow(module11_work_queue) != 516L) {
  stop_module(
    paste0(
      "Expected 516 Module 11 work-queue rows, observed ",
      nrow(module11_work_queue),
      "."
    )
  )
}

queue_clump_column <- first_matching_column(
  names(module11_work_queue),
  c("clump_id", "metadata_clump_id")
)

queue_tissue_column <- first_matching_column(
  names(module11_work_queue),
  c("tissue", "metadata_tissue")
)

queue_id_column <- first_matching_column(
  names(module11_work_queue),
  c("queue_id", "clump_tissue_key")
)

if (
  is.na(queue_clump_column) ||
    is.na(queue_tissue_column)
) {
  stop_module(
    "Module 11 work queue lacks clump_id and/or tissue."
  )
}

if (is.na(queue_id_column)) {
  module11_work_queue[
    ,
    module12_source_queue_id := sprintf(
      "SOURCE_%04d",
      .I
    )
  ]
  queue_id_column <- "module12_source_queue_id"
}

###############################################################################
# 7. RESOLVE HARMONIZED CHUNK PATHS
###############################################################################

file_qc_queue_column <- first_matching_column(
  names(file_qc),
  c("queue_id", "clump_tissue_key")
)

file_qc_clump_column <- first_matching_column(
  names(file_qc),
  c("clump_id", "metadata_clump_id")
)

file_qc_tissue_column <- first_matching_column(
  names(file_qc),
  c("tissue", "metadata_tissue")
)

file_qc_status_column <- first_matching_column(
  names(file_qc),
  c("final_file_qc", "file_status", "status")
)

harmonized_path_candidates <- c(
  "harmonized_output_file",
  "harmonized_file",
  "chunk_file",
  "chunk_output_file",
  "output_file",
  "resolved_output_file",
  "catalogue_chunk_file",
  "final_output_file"
)

file_qc_harmonized_path_column <- first_matching_column(
  names(file_qc),
  harmonized_path_candidates
)

queue_harmonized_path_column <- first_matching_column(
  names(module11_work_queue),
  harmonized_path_candidates
)

file_mapping <- data.table(
  source_queue_id = as.character(
    module11_work_queue[[queue_id_column]]
  ),
  clump_id = as.character(
    module11_work_queue[[queue_clump_column]]
  ),
  tissue = as.character(
    module11_work_queue[[queue_tissue_column]]
  )
)

file_mapping[
  ,
  harmonized_chunk_file := NA_character_
]

if (
  !is.na(file_qc_harmonized_path_column) &&
    !is.na(file_qc_queue_column)
) {
  qc_path_map <- unique(
    file_qc[
      ,
      .(
        source_queue_id = as.character(
          get(file_qc_queue_column)
        ),
        harmonized_chunk_candidate = as.character(
          get(file_qc_harmonized_path_column)
        ),
        module11_file_qc = if (
          !is.na(file_qc_status_column)
        ) {
          as.character(get(file_qc_status_column))
        } else {
          NA_character_
        }
      )
    ]
  )

  file_mapping <- merge(
    file_mapping,
    qc_path_map,
    by = "source_queue_id",
    all.x = TRUE,
    sort = FALSE
  )

  file_mapping[
    ,
    harmonized_chunk_file := vapply(
      harmonized_chunk_candidate,
      resolve_path,
      character(1L)
    )
  ]
}

if (
  !is.na(queue_harmonized_path_column)
) {
  queue_candidate_paths <- vapply(
    module11_work_queue[[queue_harmonized_path_column]],
    resolve_path,
    character(1L)
  )

  file_mapping[
    is.na(harmonized_chunk_file),
    harmonized_chunk_file := queue_candidate_paths[
      is.na(harmonized_chunk_file)
    ]
  ]
}

# Fallback: discover likely harmonized chunks and match by queue ID or CLUMP/tissue.
unresolved_chunk_rows <- which(
  is.na(file_mapping$harmonized_chunk_file)
)

if (length(unresolved_chunk_rows) > 0L) {
  cat(
    "Resolving ",
    length(unresolved_chunk_rows),
    " harmonized chunk path(s) by filename search...\n",
    sep = ""
  )

  catalogue_candidates <- list.files(
    catalogue_dir,
    pattern = "\\.(tsv|txt)(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  catalogue_candidates <- catalogue_candidates[
    !grepl(
      "final_qc|work_queue|status|summary|integrity|consistency",
      tolower(catalogue_candidates)
    )
  ]

  candidate_basenames <- basename(catalogue_candidates)

  for (mapping_index in unresolved_chunk_rows) {
    source_queue_id <- file_mapping$source_queue_id[mapping_index]
    clump_id <- file_mapping$clump_id[mapping_index]
    tissue <- file_mapping$tissue[mapping_index]

    tissue_token <- gsub("[^A-Za-z0-9]+", "_", tissue)
    tissue_token <- gsub("_+", "_", tissue_token)

    hit_index <- which(
      grepl(source_queue_id, candidate_basenames, fixed = TRUE) |
        (
          grepl(clump_id, candidate_basenames, fixed = TRUE) &
            grepl(
              substr(tolower(tissue_token), 1L, 12L),
              tolower(candidate_basenames),
              fixed = TRUE
            )
        )
    )

    if (length(hit_index) > 0L) {
      file_mapping$harmonized_chunk_file[mapping_index] <-
        normalizePath(
          catalogue_candidates[hit_index[[1L]]],
          winslash = "/",
          mustWork = TRUE
        )
    }
  }
}

file_mapping[
  ,
  chunk_exists := (
    !is.na(harmonized_chunk_file) &
      file.exists(harmonized_chunk_file)
  )
]

if (any(!file_mapping$chunk_exists)) {
  missing_chunk_table <- file_mapping[
    chunk_exists == FALSE
  ]

  fwrite(
    missing_chunk_table,
    file_mapping_output,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  stop_module(
    paste0(
      "Could not resolve ",
      nrow(missing_chunk_table),
      " harmonized apaQTL chunk file(s). ",
      "Review:\n",
      file_mapping_output
    )
  )
}

###############################################################################
# 8. RESOLVE CLUMP WORKSPACES AND GWAS FILES
###############################################################################

clump_directories <- list.dirs(
  clump_root_dir,
  recursive = FALSE,
  full.names = TRUE
)

clump_directory_map <- data.table(
  clump_dir = normalizePath(
    clump_directories,
    winslash = "/",
    mustWork = TRUE
  )
)

clump_directory_map[
  ,
  directory_name := basename(clump_dir)
]

unique_clumps <- unique(file_mapping$clump_id)

gwas_map <- data.table(
  clump_id = unique_clumps
)

gwas_map[
  ,
  clump_dir := vapply(
    clump_id,
    function(current_clump) {
      hits <- clump_directory_map[
        grepl(
          paste0("^", current_clump, "(_|$)"),
          directory_name
        ),
        clump_dir
      ]

      if (length(hits) == 0L) {
        return(NA_character_)
      }

      hits[[1L]]
    },
    character(1L)
  )
]

cat(
  "Resolving GWAS files for ",
  nrow(gwas_map),
  " CLUMPs...\n",
  sep = ""
)

gwas_map[
  ,
  gwas_file := vapply(
    seq_len(.N),
    function(index_value) {
      find_gwas_file_for_clump(
        clump_id[index_value],
        clump_dir[index_value]
      )
    },
    character(1L)
  )
]

gwas_map[
  ,
  gwas_file_exists := (
    !is.na(gwas_file) &
      file.exists(gwas_file)
  )
]

file_mapping <- merge(
  file_mapping,
  gwas_map,
  by = "clump_id",
  all.x = TRUE,
  sort = FALSE
)

fwrite(
  file_mapping,
  file_mapping_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

missing_gwas_clumps <- unique(
  file_mapping[
    gwas_file_exists == FALSE,
    clump_id
  ]
)

if (length(missing_gwas_clumps) > 0L) {
  cat(
    "WARNING: No valid GWAS file was found for ",
    length(missing_gwas_clumps),
    " CLUMP(s).\n",
    sep = ""
  )
}

###############################################################################
# 9. PROCESS EACH CLUMP × TISSUE HARMONIZED CHUNK
###############################################################################

queue_parts <- vector(
  "list",
  nrow(file_mapping)
)

for (file_index in seq_len(nrow(file_mapping))) {
  mapping_row <- file_mapping[file_index]

  cat(
    sprintf(
      "[%03d/%03d] %s | %s\n",
      file_index,
      nrow(file_mapping),
      mapping_row$clump_id,
      mapping_row$tissue
    )
  )

  apaqtl_header <- try(
    read_header(mapping_row$harmonized_chunk_file),
    silent = TRUE
  )

  if (inherits(apaqtl_header, "try-error")) {
    queue_parts[[file_index]] <- data.table(
      source_queue_id = mapping_row$source_queue_id,
      clump_id = mapping_row$clump_id,
      tissue = mapping_row$tissue,
      phenotype_id = NA_character_,
      gene_id = NA_character_,
      apaqtl_chunk_file = mapping_row$harmonized_chunk_file,
      gwas_file = mapping_row$gwas_file,
      n_apaqtl_rows = NA_integer_,
      n_apaqtl_unique_positions = NA_integer_,
      n_apaqtl_unique_allele_keys = NA_integer_,
      n_gwas_rows = NA_integer_,
      n_gwas_unique_positions = NA_integer_,
      n_gwas_unique_allele_keys = NA_integer_,
      n_shared_positions = NA_integer_,
      n_shared_allele_keys = NA_integer_,
      overlap_key_type = NA_character_,
      queue_status = "EXCLUDED",
      exclusion_reason = "APAQTL_CHUNK_READ_ERROR"
    )
    next
  }

  apaqtl_columns <- detect_apaqtl_columns(
    names(apaqtl_header)
  )

  mandatory_apaqtl_columns <- c(
    apaqtl_columns$phenotype,
    apaqtl_columns$gene,
    apaqtl_columns$chr,
    apaqtl_columns$position,
    apaqtl_columns$p
  )

  if (any(is.na(mandatory_apaqtl_columns))) {
    queue_parts[[file_index]] <- data.table(
      source_queue_id = mapping_row$source_queue_id,
      clump_id = mapping_row$clump_id,
      tissue = mapping_row$tissue,
      phenotype_id = NA_character_,
      gene_id = NA_character_,
      apaqtl_chunk_file = mapping_row$harmonized_chunk_file,
      gwas_file = mapping_row$gwas_file,
      n_apaqtl_rows = NA_integer_,
      n_apaqtl_unique_positions = NA_integer_,
      n_apaqtl_unique_allele_keys = NA_integer_,
      n_gwas_rows = NA_integer_,
      n_gwas_unique_positions = NA_integer_,
      n_gwas_unique_allele_keys = NA_integer_,
      n_shared_positions = NA_integer_,
      n_shared_allele_keys = NA_integer_,
      overlap_key_type = NA_character_,
      queue_status = "EXCLUDED",
      exclusion_reason = paste0(
        "MISSING_APAQTL_COLUMNS:",
        paste(
          names(apaqtl_columns)[is.na(unlist(apaqtl_columns))],
          collapse = ","
        )
      )
    )
    next
  }

  apaqtl_select <- unique(
    na.omit(
      unlist(
        apaqtl_columns,
        use.names = FALSE
      )
    )
  )

  apaqtl_data <- try(
    read_selected_columns(
      mapping_row$harmonized_chunk_file,
      apaqtl_select
    ),
    silent = TRUE
  )

  if (inherits(apaqtl_data, "try-error")) {
    queue_parts[[file_index]] <- data.table(
      source_queue_id = mapping_row$source_queue_id,
      clump_id = mapping_row$clump_id,
      tissue = mapping_row$tissue,
      phenotype_id = NA_character_,
      gene_id = NA_character_,
      apaqtl_chunk_file = mapping_row$harmonized_chunk_file,
      gwas_file = mapping_row$gwas_file,
      n_apaqtl_rows = NA_integer_,
      n_apaqtl_unique_positions = NA_integer_,
      n_apaqtl_unique_allele_keys = NA_integer_,
      n_gwas_rows = NA_integer_,
      n_gwas_unique_positions = NA_integer_,
      n_gwas_unique_allele_keys = NA_integer_,
      n_shared_positions = NA_integer_,
      n_shared_allele_keys = NA_integer_,
      overlap_key_type = NA_character_,
      queue_status = "EXCLUDED",
      exclusion_reason = "APAQTL_CHUNK_READ_ERROR"
    )
    next
  }

  apaqtl_data[
    ,
    module12_phenotype_id := clean_identifier(
      get(apaqtl_columns$phenotype)
    )
  ]

  apaqtl_data[
    ,
    module12_gene_id := clean_identifier(
      get(apaqtl_columns$gene)
    )
  ]

  apaqtl_data[
    ,
    module12_position_key := build_position_key(
      get(apaqtl_columns$chr),
      get(apaqtl_columns$position)
    )
  ]

  if (
    !is.na(apaqtl_columns$ref) &&
      !is.na(apaqtl_columns$alt)
  ) {
    apaqtl_data[
      ,
      module12_allele_key := build_unordered_allele_key(
        get(apaqtl_columns$chr),
        get(apaqtl_columns$position),
        get(apaqtl_columns$ref),
        get(apaqtl_columns$alt)
      )
    ]
  } else {
    apaqtl_data[
      ,
      module12_allele_key := NA_character_
    ]
  }

  apaqtl_data[
    ,
    module12_valid_p := {
      current_p <- suppressWarnings(
        as.numeric(get(apaqtl_columns$p))
      )

      !is.na(current_p) &
        is.finite(current_p) &
        current_p > 0 &
        current_p <= 1
    }
  ]

  valid_apaqtl <- apaqtl_data[
    !is.na(module12_phenotype_id) &
      !is.na(module12_gene_id) &
      !is.na(module12_position_key) &
      module12_valid_p == TRUE
  ]

  if (nrow(valid_apaqtl) == 0L) {
    queue_parts[[file_index]] <- data.table(
      source_queue_id = mapping_row$source_queue_id,
      clump_id = mapping_row$clump_id,
      tissue = mapping_row$tissue,
      phenotype_id = NA_character_,
      gene_id = NA_character_,
      apaqtl_chunk_file = mapping_row$harmonized_chunk_file,
      gwas_file = mapping_row$gwas_file,
      n_apaqtl_rows = 0L,
      n_apaqtl_unique_positions = 0L,
      n_apaqtl_unique_allele_keys = 0L,
      n_gwas_rows = NA_integer_,
      n_gwas_unique_positions = NA_integer_,
      n_gwas_unique_allele_keys = NA_integer_,
      n_shared_positions = 0L,
      n_shared_allele_keys = 0L,
      overlap_key_type = NA_character_,
      queue_status = "EXCLUDED",
      exclusion_reason = "NO_VALID_APAQTL_ROWS"
    )
    next
  }

  if (
    is.na(mapping_row$gwas_file) ||
      !file.exists(mapping_row$gwas_file)
  ) {
    apaqtl_groups <- valid_apaqtl[
      ,
      .(
        n_apaqtl_rows = .N,
        n_apaqtl_unique_positions = uniqueN(
          module12_position_key,
          na.rm = TRUE
        ),
        n_apaqtl_unique_allele_keys = uniqueN(
          module12_allele_key,
          na.rm = TRUE
        )
      ),
      by = .(
        phenotype_id = module12_phenotype_id,
        gene_id = module12_gene_id
      )
    ]

    apaqtl_groups[
      ,
      `:=`(
        source_queue_id = mapping_row$source_queue_id,
        clump_id = mapping_row$clump_id,
        tissue = mapping_row$tissue,
        apaqtl_chunk_file = mapping_row$harmonized_chunk_file,
        gwas_file = mapping_row$gwas_file,
        n_gwas_rows = NA_integer_,
        n_gwas_unique_positions = NA_integer_,
        n_gwas_unique_allele_keys = NA_integer_,
        n_shared_positions = NA_integer_,
        n_shared_allele_keys = NA_integer_,
        overlap_key_type = NA_character_,
        queue_status = "EXCLUDED",
        exclusion_reason = "GWAS_FILE_NOT_FOUND"
      )
    ]

    queue_parts[[file_index]] <- apaqtl_groups
    next
  }

  gwas_header <- try(
    read_header(mapping_row$gwas_file),
    silent = TRUE
  )

  if (inherits(gwas_header, "try-error")) {
    queue_parts[[file_index]] <- data.table(
      source_queue_id = mapping_row$source_queue_id,
      clump_id = mapping_row$clump_id,
      tissue = mapping_row$tissue,
      phenotype_id = NA_character_,
      gene_id = NA_character_,
      apaqtl_chunk_file = mapping_row$harmonized_chunk_file,
      gwas_file = mapping_row$gwas_file,
      n_apaqtl_rows = nrow(valid_apaqtl),
      n_apaqtl_unique_positions = uniqueN(
        valid_apaqtl$module12_position_key,
        na.rm = TRUE
      ),
      n_apaqtl_unique_allele_keys = uniqueN(
        valid_apaqtl$module12_allele_key,
        na.rm = TRUE
      ),
      n_gwas_rows = NA_integer_,
      n_gwas_unique_positions = NA_integer_,
      n_gwas_unique_allele_keys = NA_integer_,
      n_shared_positions = NA_integer_,
      n_shared_allele_keys = NA_integer_,
      overlap_key_type = NA_character_,
      queue_status = "EXCLUDED",
      exclusion_reason = "GWAS_READ_ERROR"
    )
    next
  }

  gwas_columns <- detect_gwas_columns(
    names(gwas_header)
  )

  mandatory_gwas_columns <- c(
    gwas_columns$chr,
    gwas_columns$position,
    gwas_columns$p
  )

  if (any(is.na(mandatory_gwas_columns))) {
    queue_parts[[file_index]] <- data.table(
      source_queue_id = mapping_row$source_queue_id,
      clump_id = mapping_row$clump_id,
      tissue = mapping_row$tissue,
      phenotype_id = NA_character_,
      gene_id = NA_character_,
      apaqtl_chunk_file = mapping_row$harmonized_chunk_file,
      gwas_file = mapping_row$gwas_file,
      n_apaqtl_rows = nrow(valid_apaqtl),
      n_apaqtl_unique_positions = uniqueN(
        valid_apaqtl$module12_position_key,
        na.rm = TRUE
      ),
      n_apaqtl_unique_allele_keys = uniqueN(
        valid_apaqtl$module12_allele_key,
        na.rm = TRUE
      ),
      n_gwas_rows = NA_integer_,
      n_gwas_unique_positions = NA_integer_,
      n_gwas_unique_allele_keys = NA_integer_,
      n_shared_positions = NA_integer_,
      n_shared_allele_keys = NA_integer_,
      overlap_key_type = NA_character_,
      queue_status = "EXCLUDED",
      exclusion_reason = "MISSING_GWAS_COLUMNS"
    )
    next
  }

  gwas_select <- unique(
    na.omit(
      unlist(
        gwas_columns,
        use.names = FALSE
      )
    )
  )

  gwas_data <- try(
    read_selected_columns(
      mapping_row$gwas_file,
      gwas_select
    ),
    silent = TRUE
  )

  if (inherits(gwas_data, "try-error")) {
    queue_parts[[file_index]] <- data.table(
      source_queue_id = mapping_row$source_queue_id,
      clump_id = mapping_row$clump_id,
      tissue = mapping_row$tissue,
      phenotype_id = NA_character_,
      gene_id = NA_character_,
      apaqtl_chunk_file = mapping_row$harmonized_chunk_file,
      gwas_file = mapping_row$gwas_file,
      n_apaqtl_rows = nrow(valid_apaqtl),
      n_apaqtl_unique_positions = uniqueN(
        valid_apaqtl$module12_position_key,
        na.rm = TRUE
      ),
      n_apaqtl_unique_allele_keys = uniqueN(
        valid_apaqtl$module12_allele_key,
        na.rm = TRUE
      ),
      n_gwas_rows = NA_integer_,
      n_gwas_unique_positions = NA_integer_,
      n_gwas_unique_allele_keys = NA_integer_,
      n_shared_positions = NA_integer_,
      n_shared_allele_keys = NA_integer_,
      overlap_key_type = NA_character_,
      queue_status = "EXCLUDED",
      exclusion_reason = "GWAS_READ_ERROR"
    )
    next
  }

  gwas_data[
    ,
    module12_position_key := build_position_key(
      get(gwas_columns$chr),
      get(gwas_columns$position)
    )
  ]

  if (
    !is.na(gwas_columns$effect) &&
      !is.na(gwas_columns$other)
  ) {
    gwas_data[
      ,
      module12_allele_key := build_unordered_allele_key(
        get(gwas_columns$chr),
        get(gwas_columns$position),
        get(gwas_columns$effect),
        get(gwas_columns$other)
      )
    ]
  } else {
    gwas_data[
      ,
      module12_allele_key := NA_character_
    ]
  }

  gwas_data[
    ,
    module12_valid_p := {
      current_p <- suppressWarnings(
        as.numeric(get(gwas_columns$p))
      )

      !is.na(current_p) &
        is.finite(current_p) &
        current_p > 0 &
        current_p <= 1
    }
  ]

  valid_gwas <- gwas_data[
    !is.na(module12_position_key) &
      module12_valid_p == TRUE
  ]

  gwas_position_keys <- unique(
    valid_gwas$module12_position_key[
      !is.na(valid_gwas$module12_position_key)
    ]
  )

  gwas_allele_keys <- unique(
    valid_gwas$module12_allele_key[
      !is.na(valid_gwas$module12_allele_key)
    ]
  )

  grouped_queue <- valid_apaqtl[
    ,
    {
      current_position_keys <- unique(
        module12_position_key[
          !is.na(module12_position_key)
        ]
      )

      current_allele_keys <- unique(
        module12_allele_key[
          !is.na(module12_allele_key)
        ]
      )

      shared_positions <- length(
        intersect(
          current_position_keys,
          gwas_position_keys
        )
      )

      shared_allele_keys <- length(
        intersect(
          current_allele_keys,
          gwas_allele_keys
        )
      )

      use_allele_key <- (
        length(current_allele_keys) > 0L &&
          length(gwas_allele_keys) > 0L
      )

      effective_shared_count <- if (use_allele_key) {
        shared_allele_keys
      } else {
        shared_positions
      }

      effective_key_type <- if (use_allele_key) {
        "CHR_POS_UNORDERED_ALLELES"
      } else {
        "CHR_POS"
      }

      exclusion_reasons <- character()

      if (.N < minimum_apaqtl_rows) {
        exclusion_reasons <- c(
          exclusion_reasons,
          "INSUFFICIENT_APAQTL_ROWS"
        )
      }

      if (nrow(valid_gwas) < minimum_gwas_rows) {
        exclusion_reasons <- c(
          exclusion_reasons,
          "INSUFFICIENT_GWAS_ROWS"
        )
      }

      if (effective_shared_count < minimum_shared_variants) {
        exclusion_reasons <- c(
          exclusion_reasons,
          "INSUFFICIENT_SHARED_VARIANTS"
        )
      }

      current_status <- if (
        length(exclusion_reasons) == 0L
      ) {
        "ELIGIBLE"
      } else {
        "EXCLUDED"
      }

      current_reason <- if (
        length(exclusion_reasons) == 0L
      ) {
        "NONE"
      } else {
        paste(
          unique(exclusion_reasons),
          collapse = ";"
        )
      }

      list(
        n_apaqtl_rows = .N,
        n_apaqtl_unique_positions = length(
          current_position_keys
        ),
        n_apaqtl_unique_allele_keys = length(
          current_allele_keys
        ),
        n_gwas_rows = nrow(valid_gwas),
        n_gwas_unique_positions = length(
          gwas_position_keys
        ),
        n_gwas_unique_allele_keys = length(
          gwas_allele_keys
        ),
        n_shared_positions = shared_positions,
        n_shared_allele_keys = shared_allele_keys,
        effective_shared_variants = effective_shared_count,
        overlap_key_type = effective_key_type,
        queue_status = current_status,
        exclusion_reason = current_reason
      )
    },
    by = .(
      phenotype_id = module12_phenotype_id,
      gene_id = module12_gene_id
    )
  ]

  grouped_queue[
    ,
    `:=`(
      source_queue_id = mapping_row$source_queue_id,
      clump_id = mapping_row$clump_id,
      tissue = mapping_row$tissue,
      apaqtl_chunk_file = mapping_row$harmonized_chunk_file,
      gwas_file = mapping_row$gwas_file
    )
  ]

  setcolorder(
    grouped_queue,
    c(
      "source_queue_id",
      "clump_id",
      "tissue",
      "phenotype_id",
      "gene_id",
      "apaqtl_chunk_file",
      "gwas_file",
      "n_apaqtl_rows",
      "n_apaqtl_unique_positions",
      "n_apaqtl_unique_allele_keys",
      "n_gwas_rows",
      "n_gwas_unique_positions",
      "n_gwas_unique_allele_keys",
      "n_shared_positions",
      "n_shared_allele_keys",
      "effective_shared_variants",
      "overlap_key_type",
      "queue_status",
      "exclusion_reason"
    )
  )

  queue_parts[[file_index]] <- grouped_queue
}

###############################################################################
# 10. COMBINE AND FINALIZE WORK QUEUE
###############################################################################

coloc_queue <- rbindlist(
  queue_parts,
  use.names = TRUE,
  fill = TRUE
)

if (nrow(coloc_queue) == 0L) {
  stop_module(
    "No colocalization units were generated."
  )
}

coloc_queue[
  ,
  coloc_unit_id := sprintf(
    "COLOC_%07d",
    seq_len(.N)
  )
]

coloc_queue[
  ,
  clump_tissue_phenotype_gene_key := paste(
    clump_id,
    tissue,
    phenotype_id,
    gene_id,
    sep = "|"
  )
]

duplicate_unit_rows <- coloc_queue[
  duplicated(
    clump_tissue_phenotype_gene_key
  ) |
    duplicated(
      clump_tissue_phenotype_gene_key,
      fromLast = TRUE
    )
]

if (nrow(duplicate_unit_rows) > 0L) {
  duplicate_output <- file.path(
    qc_dir,
    "module12_1_duplicate_units.tsv"
  )

  fwrite(
    duplicate_unit_rows,
    duplicate_output,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  stop_module(
    paste0(
      "Duplicate CLUMP × tissue × phenotype × gene units detected. Review:\n",
      duplicate_output
    )
  )
}

setorder(
  coloc_queue,
  clump_id,
  tissue,
  gene_id,
  phenotype_id
)

coloc_queue[
  ,
  queue_order := seq_len(.N)
]

setcolorder(
  coloc_queue,
  c(
    "queue_order",
    "coloc_unit_id",
    "source_queue_id",
    "clump_id",
    "tissue",
    "phenotype_id",
    "gene_id",
    "clump_tissue_phenotype_gene_key",
    setdiff(
      names(coloc_queue),
      c(
        "queue_order",
        "coloc_unit_id",
        "source_queue_id",
        "clump_id",
        "tissue",
        "phenotype_id",
        "gene_id",
        "clump_tissue_phenotype_gene_key"
      )
    )
  )
)

eligible_queue <- coloc_queue[
  queue_status == "ELIGIBLE"
]

excluded_queue <- coloc_queue[
  queue_status == "EXCLUDED"
]

###############################################################################
# 11. QC SUMMARIES
###############################################################################

exclusion_reason_summary <- excluded_queue[
  ,
  .N,
  by = exclusion_reason
][
  order(-N, exclusion_reason)
]

tissue_summary <- coloc_queue[
  ,
  .(
    coloc_units = .N,
    eligible_units = sum(
      queue_status == "ELIGIBLE",
      na.rm = TRUE
    ),
    excluded_units = sum(
      queue_status == "EXCLUDED",
      na.rm = TRUE
    ),
    unique_clumps = uniqueN(clump_id),
    unique_genes = uniqueN(
      gene_id,
      na.rm = TRUE
    ),
    unique_phenotypes = uniqueN(
      phenotype_id,
      na.rm = TRUE
    ),
    median_shared_variants = median(
      effective_shared_variants,
      na.rm = TRUE
    ),
    minimum_shared_variants = suppressWarnings(
      min(
        effective_shared_variants,
        na.rm = TRUE
      )
    ),
    maximum_shared_variants = suppressWarnings(
      max(
        effective_shared_variants,
        na.rm = TRUE
      )
    )
  ),
  by = tissue
]

clump_summary <- coloc_queue[
  ,
  .(
    coloc_units = .N,
    eligible_units = sum(
      queue_status == "ELIGIBLE",
      na.rm = TRUE
    ),
    excluded_units = sum(
      queue_status == "EXCLUDED",
      na.rm = TRUE
    ),
    tissues = uniqueN(tissue),
    unique_genes = uniqueN(
      gene_id,
      na.rm = TRUE
    ),
    unique_phenotypes = uniqueN(
      phenotype_id,
      na.rm = TRUE
    )
  ),
  by = clump_id
][
  order(clump_id)
]

module_end_time <- Sys.time()

runtime_seconds <- as.numeric(
  difftime(
    module_end_time,
    module_start_time,
    units = "secs"
  )
)

missing_chunk_count <- sum(
  !file_mapping$chunk_exists,
  na.rm = TRUE
)

missing_gwas_clump_count <- uniqueN(
  file_mapping[
    gwas_file_exists == FALSE,
    clump_id
  ]
)

eligible_count <- nrow(eligible_queue)
excluded_count <- nrow(excluded_queue)

final_status <- if (
  missing_chunk_count > 0L ||
    nrow(coloc_queue) == 0L
) {
  "FAIL"
} else if (
  missing_gwas_clump_count > 0L ||
    excluded_count > 0L
) {
  "PASS_WITH_WARNINGS"
} else {
  "PASS"
}

qc_summary <- data.table(
  module = "12.1",
  module_name = "initialize GWAS-apaQTL colocalization work queue",
  start_time = format(
    module_start_time,
    "%Y-%m-%d %H:%M:%S"
  ),
  end_time = format(
    module_end_time,
    "%Y-%m-%d %H:%M:%S"
  ),
  runtime_seconds = runtime_seconds,
  module11_final_status = module11_final_status,
  expected_clump_tissue_files = 516L,
  mapped_apaqtl_chunks = sum(
    file_mapping$chunk_exists,
    na.rm = TRUE
  ),
  unique_clumps = uniqueN(
    file_mapping$clump_id
  ),
  unique_tissues = uniqueN(
    file_mapping$tissue
  ),
  clumps_with_gwas_file = uniqueN(
    file_mapping[
      gwas_file_exists == TRUE,
      clump_id
    ]
  ),
  clumps_without_gwas_file = missing_gwas_clump_count,
  total_coloc_units = nrow(coloc_queue),
  eligible_coloc_units = eligible_count,
  excluded_coloc_units = excluded_count,
  minimum_shared_variants = minimum_shared_variants,
  duplicate_units = nrow(duplicate_unit_rows),
  final_status = final_status
)

status_table <- data.table(
  module = "12.1",
  final_status = final_status,
  work_queue_file = normalizePath(
    all_queue_output,
    winslash = "/",
    mustWork = FALSE
  ),
  eligible_queue_file = normalizePath(
    eligible_queue_output,
    winslash = "/",
    mustWork = FALSE
  ),
  excluded_queue_file = normalizePath(
    excluded_queue_output,
    winslash = "/",
    mustWork = FALSE
  ),
  total_units = nrow(coloc_queue),
  eligible_units = eligible_count,
  excluded_units = excluded_count,
  created_at = format(
    module_end_time,
    "%Y-%m-%d %H:%M:%S"
  )
)

###############################################################################
# 12. WRITE OUTPUTS
###############################################################################

fwrite(
  coloc_queue,
  all_queue_output,
  sep = "\t",
  quote = FALSE,
  na = "NA",
  compress = "gzip"
)

fwrite(
  eligible_queue,
  eligible_queue_output,
  sep = "\t",
  quote = FALSE,
  na = "NA",
  compress = "gzip"
)

fwrite(
  excluded_queue,
  excluded_queue_output,
  sep = "\t",
  quote = FALSE,
  na = "NA",
  compress = "gzip"
)

fwrite(
  file_mapping,
  file_mapping_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  qc_summary,
  qc_summary_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  exclusion_reason_summary,
  exclusion_reason_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  tissue_summary,
  tissue_summary_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  clump_summary,
  clump_summary_output,
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

writeLines(
  capture.output(sessionInfo()),
  session_info_output,
  useBytes = TRUE
)

###############################################################################
# 13. FINAL REPORT
###############################################################################

cat("\n============================================================\n")
cat("MODULE 12.1 – FINAL SUMMARY\n")
cat("============================================================\n\n")

print(qc_summary)

cat("\nQueue-status distribution:\n")
print(
  coloc_queue[
    ,
    .N,
    by = queue_status
  ][
    order(-N)
  ]
)

cat("\nExclusion-reason distribution:\n")
print(exclusion_reason_summary)

cat("\nTissue summary:\n")
print(tissue_summary)

cat("\nOutputs:\n")
cat("  Complete work queue:\n    ", all_queue_output, "\n", sep = "")
cat("  Eligible work queue:\n    ", eligible_queue_output, "\n", sep = "")
cat("  Excluded work queue:\n    ", excluded_queue_output, "\n", sep = "")
cat("  File mapping:\n    ", file_mapping_output, "\n", sep = "")
cat("  QC summary:\n    ", qc_summary_output, "\n", sep = "")
cat("  Status:\n    ", status_output, "\n", sep = "")
cat("  Log:\n    ", log_file, "\n\n", sep = "")

cat("============================================================\n")

if (identical(final_status, "PASS")) {
  cat("FINAL RESULT: PASS\n")
  cat("The colocalization work queue was initialized successfully.\n")
  cat("Module 12.2 may be started.\n")
} else if (identical(final_status, "PASS_WITH_WARNINGS")) {
  cat("FINAL RESULT: PASS WITH WARNINGS\n")
  cat("The work queue was created, but excluded units or missing GWAS files\n")
  cat("must be reviewed before the full colocalization run.\n")
  cat("Eligible units may proceed to Module 12.2 testing.\n")
} else {
  cat("FINAL RESULT: FAIL\n")
  cat("The colocalization work queue could not be initialized safely.\n")
}

cat("============================================================\n")

if (identical(final_status, "FAIL")) {
  stop(
    paste0(
      "Module 12.1 failed. Review:\n",
      qc_summary_output
    ),
    call. = FALSE
  )
}

invisible(
  list(
    final_status = status_table,
    work_queue = coloc_queue,
    eligible_queue = eligible_queue,
    excluded_queue = excluded_queue,
    file_mapping = file_mapping,
    qc_summary = qc_summary
  )
)
