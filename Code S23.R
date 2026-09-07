###############################################################################
# PigmentationAtlas
# Module 11.3B.2
#
# READ AND QC A SINGLE apaQTL FILE
#
# Purpose:
#   1. Define read_and_qc_apaqtl()
#   2. Detect apaQTL source columns
#   3. Standardize variant, phenotype, position and association statistics
#   4. Add canonical CLUMP × tissue metadata
#   5. Perform file-level and row-level QC
#   6. Test the function on the first work-queue file
#
# Input:
#   results/module11/harmonization/catalogue/
#   module11_3B_work_queue.tsv
#
# Test outputs:
#   results/module11/harmonization/catalogue/qc/
#   module11_3B_2_test_file_qc.tsv
#   module11_3B_2_test_column_map.tsv
#   module11_3B_2_test_duplicate_rows.tsv
#   module11_3B_2_test_harmonized_preview.tsv
#
# Important:
#   This module tests one file only.
#   The complete 516-file catalogue is built in Module 11.3B.3.
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
cat("MODULE 11.3B.2 – SINGLE-FILE READER AND QC\n")
cat("Started:", format(module_start_time), "\n")
cat("============================================================\n\n")


###############################################################################
# 2. REQUIRED PACKAGE
###############################################################################

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop(
    "The data.table package is required. Install it with install.packages('data.table').",
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

qc_dir <- file.path(catalogue_dir, "qc")
log_dir <- file.path(catalogue_dir, "logs")

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

work_queue_file <- file.path(
  catalogue_dir,
  "module11_3B_work_queue.tsv"
)

test_file_qc_output <- file.path(
  qc_dir,
  "module11_3B_2_test_file_qc.tsv"
)

test_column_map_output <- file.path(
  qc_dir,
  "module11_3B_2_test_column_map.tsv"
)

test_duplicate_output <- file.path(
  qc_dir,
  "module11_3B_2_test_duplicate_rows.tsv"
)

test_preview_output <- file.path(
  qc_dir,
  "module11_3B_2_test_harmonized_preview.tsv"
)

test_full_output <- file.path(
  qc_dir,
  "module11_3B_2_test_harmonized_full.tsv.gz"
)

test_status_output <- file.path(
  catalogue_dir,
  "module11_3B_2_test_status.tsv"
)

session_info_output <- file.path(
  log_dir,
  "module11_3B_2_session_info.txt"
)


###############################################################################
# 4. GENERAL HELPER FUNCTIONS
###############################################################################

normalize_column_names <- function(x) {

  x <- trimws(as.character(x))
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x <- gsub("_+", "_", x)

  x
}


normalize_chromosome <- function(x) {

  x <- toupper(trimws(as.character(x)))
  x <- gsub("^CHR", "", x)

  x[x %chin% c("23")] <- "X"
  x[x %chin% c("24")] <- "Y"
  x[x %chin% c("25", "M", "MITO")] <- "MT"
  x[x == ""] <- NA_character_

  x
}


normalize_allele <- function(x) {

  x <- toupper(trimws(as.character(x)))
  x[x == ""] <- NA_character_

  x
}


find_first_column <- function(column_names,
                              candidates,
                              required = FALSE,
                              field_name = NULL) {

  hits <- candidates[candidates %in% column_names]

  if (length(hits) > 0L) {
    return(hits[1L])
  }

  if (isTRUE(required)) {

    if (is.null(field_name)) {
      field_name <- paste(candidates, collapse = " / ")
    }

    stop(
      paste0(
        "Required field could not be identified: ",
        field_name,
        "\nAvailable columns:\n",
        paste(column_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  NA_character_
}


safe_numeric <- function(x) {

  suppressWarnings(
    as.numeric(as.character(x))
  )
}


count_missing <- function(x) {

  sum(
    is.na(x) |
      trimws(as.character(x)) == ""
  )
}


empty_table_with_columns <- function(column_names) {

  as.data.table(
    setNames(
      replicate(
        length(column_names),
        character(0),
        simplify = FALSE
      ),
      column_names
    )
  )
}


###############################################################################
# 5. VARIANT-ID PARSING
###############################################################################

parse_variant_id <- function(variant_id) {

  variant_id <- as.character(variant_id)

  n <- length(variant_id)

  result <- data.table(
    parsed_chromosome = rep(NA_character_, n),
    parsed_position = rep(NA_real_, n),
    parsed_ref = rep(NA_character_, n),
    parsed_alt = rep(NA_character_, n),
    variant_parse_status = rep("UNPARSED", n)
  )

  if (n == 0L) {
    return(result)
  }

  cleaned <- trimws(variant_id)

  # Pattern examples:
  # chr11_69049973_T_G_b38
  # 11_69049973_T_G_b38
  # chr11:69049973:T:G
  # 11:69049973:T:G

  pattern_four_fields <- paste0(
    "^(?:chr)?",
    "([0-9]+|X|Y|MT|M)",
    "[:_]",
    "([0-9]+)",
    "[:_]",
    "([A-Za-z]+)",
    "[:_]",
    "([A-Za-z]+)"
  )

  matched_four <- regexec(
    pattern_four_fields,
    cleaned,
    ignore.case = TRUE,
    perl = TRUE
  )

  pieces_four <- regmatches(cleaned, matched_four)

  four_ok <- lengths(pieces_four) >= 5L

  if (any(four_ok)) {

    result[
      four_ok,
      parsed_chromosome := normalize_chromosome(
        vapply(
          pieces_four[four_ok],
          `[`,
          character(1),
          2L
        )
      )
    ]

    result[
      four_ok,
      parsed_position := safe_numeric(
        vapply(
          pieces_four[four_ok],
          `[`,
          character(1),
          3L
        )
      )
    ]

    result[
      four_ok,
      parsed_ref := normalize_allele(
        vapply(
          pieces_four[four_ok],
          `[`,
          character(1),
          4L
        )
      )
    ]

    result[
      four_ok,
      parsed_alt := normalize_allele(
        vapply(
          pieces_four[four_ok],
          `[`,
          character(1),
          5L
        )
      )
    ]

    result[
      four_ok,
      variant_parse_status := "PARSED_CHR_POS_REF_ALT"
    ]
  }

  # Simpler pattern:
  # chr11_69049973
  # chr11:69049973
  # 11_69049973

  still_unparsed <- result$variant_parse_status == "UNPARSED"

  if (any(still_unparsed)) {

    pattern_two_fields <- paste0(
      "^(?:chr)?",
      "([0-9]+|X|Y|MT|M)",
      "[:_]",
      "([0-9]+)"
    )

    matched_two <- regexec(
      pattern_two_fields,
      cleaned[still_unparsed],
      ignore.case = TRUE,
      perl = TRUE
    )

    pieces_two <- regmatches(
      cleaned[still_unparsed],
      matched_two
    )

    two_ok <- lengths(pieces_two) >= 3L

    target_indices <- which(still_unparsed)[two_ok]

    if (length(target_indices) > 0L) {

      result[
        target_indices,
        parsed_chromosome := normalize_chromosome(
          vapply(
            pieces_two[two_ok],
            `[`,
            character(1),
            2L
          )
        )
      ]

      result[
        target_indices,
        parsed_position := safe_numeric(
          vapply(
            pieces_two[two_ok],
            `[`,
            character(1),
            3L
          )
        )
      ]

      result[
        target_indices,
        variant_parse_status := "PARSED_CHR_POS"
      ]
    }
  }

  result
}


###############################################################################
# 6. apaQTL COLUMN DICTIONARY
###############################################################################

apaqtl_column_candidates <- list(

  variant_id = c(
    "variant_id",
    "variant",
    "variantid",
    "snp",
    "snp_id",
    "snpid",
    "marker",
    "markername",
    "id"
  ),

  rsid = c(
    "rsid",
    "rs_id",
    "dbsnp",
    "dbsnp_id"
  ),

  chromosome = c(
    "chromosome",
    "chrom",
    "chr",
    "variant_chromosome",
    "variant_chr"
  ),

  position = c(
    "position",
    "pos",
    "bp",
    "base_pair_location",
    "base_pair_position",
    "variant_position",
    "variant_pos"
  ),

  ref_allele = c(
    "ref",
    "reference_allele",
    "ref_allele",
    "allele0",
    "a0",
    "other_allele",
    "non_effect_allele"
  ),

  alt_allele = c(
    "alt",
    "alternate_allele",
    "alternative_allele",
    "alt_allele",
    "allele1",
    "a1",
    "effect_allele"
  ),

  phenotype_id = c(
    "phenotype_id",
    "phenotype",
    "phenotypeid",
    "molecular_trait_id",
    "molecular_trait",
    "gene_id",
    "gene",
    "geneid",
    "transcript_id",
    "transcript",
    "intron_id",
    "event_id",
    "cluster_id"
  ),

  gene_id = c(
    "gene_id",
    "gene",
    "geneid",
    "ensembl_gene_id",
    "ensembl_id"
  ),

  p_value = c(
    "p_value",
    "pvalue",
    "pval",
    "pval_nominal",
    "nominal_p_value",
    "p_value_nominal",
    "p"
  ),

  beta = c(
    "beta",
    "slope",
    "effect",
    "effect_size",
    "estimate",
    "nes"
  ),

  standard_error = c(
    "standard_error",
    "standarderror",
    "se",
    "slope_se",
    "beta_se",
    "stderr"
  ),

  maf = c(
    "maf",
    "minor_allele_frequency"
  ),

  allele_frequency = c(
    "allele_frequency",
    "effect_allele_frequency",
    "eaf",
    "af",
    "a1_frequency",
    "alt_af"
  ),

  sample_size = c(
    "sample_size",
    "samplesize",
    "n",
    "n_samples",
    "n_total"
  ),

  molecular_trait_chromosome = c(
    "molecular_trait_chromosome",
    "phenotype_chromosome",
    "gene_chromosome",
    "gene_chr"
  ),

  molecular_trait_position = c(
    "molecular_trait_position",
    "phenotype_position",
    "gene_position",
    "tss",
    "tss_position",
    "phenotype_pos"
  ),

  distance = c(
    "distance",
    "tss_distance",
    "distance_to_tss",
    "variant_phenotype_distance"
  )
)


###############################################################################
# 7. COLUMN-DETECTION FUNCTION
###############################################################################

detect_apaqtl_columns <- function(column_names) {

  mapping <- lapply(
    names(apaqtl_column_candidates),
    function(field) {

      find_first_column(
        column_names = column_names,
        candidates = apaqtl_column_candidates[[field]],
        required = FALSE,
        field_name = field
      )
    }
  )

  names(mapping) <- names(apaqtl_column_candidates)

  mapping
}


make_column_map_table <- function(mapping) {

  data.table(
    standardized_field = names(mapping),
    source_column = unlist(mapping, use.names = FALSE),
    detected = !is.na(unlist(mapping, use.names = FALSE))
  )
}


###############################################################################
# 8. SAFE COLUMN EXTRACTION
###############################################################################

extract_character_column <- function(dt,
                                     source_column,
                                     default = NA_character_) {

  if (is.na(source_column) || !source_column %in% names(dt)) {
    return(rep(default, nrow(dt)))
  }

  as.character(dt[[source_column]])
}


extract_numeric_column <- function(dt,
                                   source_column,
                                   default = NA_real_) {

  if (is.na(source_column) || !source_column %in% names(dt)) {
    return(rep(default, nrow(dt)))
  }

  safe_numeric(dt[[source_column]])
}


###############################################################################
# 9. MAIN FUNCTION
###############################################################################

read_and_qc_apaqtl <- function(file_path,
                               metadata,
                               retain_source_columns = TRUE,
                               fail_on_missing_pvalue = TRUE,
                               fail_on_missing_variant = TRUE) {

  file_start_time <- Sys.time()

  required_metadata <- c(
    "queue_id",
    "clump_id",
    "tissue",
    "canonical_chromosome",
    "canonical_region_start",
    "canonical_region_end"
  )

  missing_metadata <- setdiff(
    required_metadata,
    names(metadata)
  )

  if (length(missing_metadata) > 0L) {

    stop(
      paste0(
        "Required metadata field(s) missing: ",
        paste(missing_metadata, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (length(file_path) != 1L ||
      is.na(file_path) ||
      trimws(file_path) == "") {

    stop("A single valid file path must be supplied.", call. = FALSE)
  }

  if (!file.exists(file_path)) {

    stop(
      paste0(
        "apaQTL file does not exist:\n",
        file_path
      ),
      call. = FALSE
    )
  }

  file_size_bytes <- as.numeric(
    file.info(file_path)$size
  )

  if (is.na(file_size_bytes) || file_size_bytes <= 0) {

    stop(
      paste0(
        "apaQTL file is empty:\n",
        file_path
      ),
      call. = FALSE
    )
  }

  source_dt <- tryCatch(

    fread(
      file_path,
      sep = "auto",
      header = TRUE,
      na.strings = c(
        "",
        "NA",
        "NaN",
        "NULL",
        "."
      ),
      showProgress = FALSE
    ),

    error = function(e) {

      stop(
        paste0(
          "Failed to read apaQTL file:\n",
          file_path,
          "\nOriginal error:\n",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  if (nrow(source_dt) == 0L) {

    stop(
      paste0(
        "apaQTL file contains a header but no data rows:\n",
        file_path
      ),
      call. = FALSE
    )
  }

  source_column_names_original <- names(source_dt)
  source_column_names_normalized <- normalize_column_names(
    source_column_names_original
  )

  if (anyDuplicated(source_column_names_normalized)) {

    duplicate_names <- unique(
      source_column_names_normalized[
        duplicated(source_column_names_normalized) |
          duplicated(
            source_column_names_normalized,
            fromLast = TRUE
          )
      ]
    )

    stop(
      paste0(
        "Column-name normalization produced duplicate columns:\n",
        paste(duplicate_names, collapse = ", "),
        "\nFile:\n",
        file_path
      ),
      call. = FALSE
    )
  }

  setnames(
    source_dt,
    old = source_column_names_original,
    new = source_column_names_normalized
  )

  column_mapping <- detect_apaqtl_columns(
    names(source_dt)
  )

  column_map_table <- make_column_map_table(
    column_mapping
  )

  variant_source_available <- any(
    !is.na(
      unlist(
        column_mapping[
          c(
            "variant_id",
            "rsid",
            "chromosome",
            "position"
          )
        ]
      )
    )
  )

  if (isTRUE(fail_on_missing_variant) &&
      !variant_source_available) {

    stop(
      paste0(
        "No usable variant identifier or genomic position column was found.\n",
        "Available columns:\n",
        paste(names(source_dt), collapse = ", "),
        "\nFile:\n",
        file_path
      ),
      call. = FALSE
    )
  }

  if (isTRUE(fail_on_missing_pvalue) &&
      is.na(column_mapping$p_value)) {

    stop(
      paste0(
        "No p-value column was identified.\n",
        "Available columns:\n",
        paste(names(source_dt), collapse = ", "),
        "\nFile:\n",
        file_path
      ),
      call. = FALSE
    )
  }

  variant_id_raw <- extract_character_column(
    source_dt,
    column_mapping$variant_id
  )

  rsid_raw <- extract_character_column(
    source_dt,
    column_mapping$rsid
  )

  chromosome_direct <- normalize_chromosome(
    extract_character_column(
      source_dt,
      column_mapping$chromosome
    )
  )

  position_direct <- extract_numeric_column(
    source_dt,
    column_mapping$position
  )

  ref_direct <- normalize_allele(
    extract_character_column(
      source_dt,
      column_mapping$ref_allele
    )
  )

  alt_direct <- normalize_allele(
    extract_character_column(
      source_dt,
      column_mapping$alt_allele
    )
  )

  parsed_variant <- parse_variant_id(
    variant_id_raw
  )

  chromosome_final <- fifelse(
    !is.na(chromosome_direct),
    chromosome_direct,
    parsed_variant$parsed_chromosome
  )

  position_final <- fifelse(
    !is.na(position_direct),
    position_direct,
    parsed_variant$parsed_position
  )

  ref_final <- fifelse(
    !is.na(ref_direct),
    ref_direct,
    parsed_variant$parsed_ref
  )

  alt_final <- fifelse(
    !is.na(alt_direct),
    alt_direct,
    parsed_variant$parsed_alt
  )

  reconstructed_variant_id <- fifelse(
    !is.na(chromosome_final) &
      !is.na(position_final) &
      !is.na(ref_final) &
      !is.na(alt_final),
    paste(
      paste0("chr", chromosome_final),
      as.integer(position_final),
      ref_final,
      alt_final,
      sep = "_"
    ),
    NA_character_
  )

  variant_id_final <- fifelse(
    !is.na(variant_id_raw) &
      trimws(variant_id_raw) != "",
    variant_id_raw,
    reconstructed_variant_id
  )

  phenotype_id <- extract_character_column(
    source_dt,
    column_mapping$phenotype_id
  )

  gene_id <- extract_character_column(
    source_dt,
    column_mapping$gene_id
  )

  phenotype_id <- fifelse(
    !is.na(phenotype_id) &
      trimws(phenotype_id) != "",
    phenotype_id,
    gene_id
  )

  p_value <- extract_numeric_column(
    source_dt,
    column_mapping$p_value
  )

  beta <- extract_numeric_column(
    source_dt,
    column_mapping$beta
  )

  standard_error <- extract_numeric_column(
    source_dt,
    column_mapping$standard_error
  )

  maf <- extract_numeric_column(
    source_dt,
    column_mapping$maf
  )

  allele_frequency <- extract_numeric_column(
    source_dt,
    column_mapping$allele_frequency
  )

  sample_size <- extract_numeric_column(
    source_dt,
    column_mapping$sample_size
  )

  molecular_trait_chromosome <- normalize_chromosome(
    extract_character_column(
      source_dt,
      column_mapping$molecular_trait_chromosome
    )
  )

  molecular_trait_position <- extract_numeric_column(
    source_dt,
    column_mapping$molecular_trait_position
  )

  distance <- extract_numeric_column(
    source_dt,
    column_mapping$distance
  )

  canonical_chromosome <- normalize_chromosome(
    metadata$canonical_chromosome[1L]
  )

  canonical_region_start <- safe_numeric(
    metadata$canonical_region_start[1L]
  )

  canonical_region_end <- safe_numeric(
    metadata$canonical_region_end[1L]
  )

  chromosome_matches_canonical <- !is.na(chromosome_final) &
    chromosome_final == canonical_chromosome

  position_within_canonical_region <- !is.na(position_final) &
    position_final >= canonical_region_start &
    position_final <= canonical_region_end

  pvalue_valid <- !is.na(p_value) &
    p_value >= 0 &
    p_value <= 1

  allele_frequency_valid <- is.na(allele_frequency) |
    (
      allele_frequency >= 0 &
        allele_frequency <= 1
    )

  maf_valid <- is.na(maf) |
    (
      maf >= 0 &
        maf <= 0.5
    )

  harmonized_dt <- data.table(

    queue_id = rep(
      as.character(metadata$queue_id[1L]),
      nrow(source_dt)
    ),

    clump_id = rep(
      as.character(metadata$clump_id[1L]),
      nrow(source_dt)
    ),

    tissue = rep(
      as.character(metadata$tissue[1L]),
      nrow(source_dt)
    ),

    canonical_chromosome = rep(
      canonical_chromosome,
      nrow(source_dt)
    ),

    canonical_region_start = rep(
      canonical_region_start,
      nrow(source_dt)
    ),

    canonical_region_end = rep(
      canonical_region_end,
      nrow(source_dt)
    ),

    variant_id = variant_id_final,
    rsid = rsid_raw,
    chromosome = chromosome_final,
    position = position_final,
    ref_allele = ref_final,
    alt_allele = alt_final,

    phenotype_id = phenotype_id,
    gene_id = gene_id,

    p_value = p_value,
    beta = beta,
    standard_error = standard_error,

    maf = maf,
    allele_frequency = allele_frequency,
    sample_size = sample_size,

    molecular_trait_chromosome =
      molecular_trait_chromosome,

    molecular_trait_position =
      molecular_trait_position,

    distance = distance,

    variant_parse_status =
      parsed_variant$variant_parse_status,

    chromosome_matches_canonical =
      chromosome_matches_canonical,

    position_within_canonical_region =
      position_within_canonical_region,

    pvalue_valid = pvalue_valid,
    maf_valid = maf_valid,

    allele_frequency_valid =
      allele_frequency_valid,

    source_file = rep(
      normalizePath(
        file_path,
        winslash = "/",
        mustWork = TRUE
      ),
      nrow(source_dt)
    ),

    source_row = seq_len(nrow(source_dt))
  )

  # Stable record key.
  harmonized_dt[
    ,
    record_key := paste(
      clump_id,
      tissue,
      fifelse(
        !is.na(variant_id),
        variant_id,
        paste(chromosome, position, sep = ":")
      ),
      fifelse(
        !is.na(phenotype_id),
        phenotype_id,
        "<NO_PHENOTYPE>"
      ),
      sep = "||"
    )
  ]

  # Duplicate definition:
  # same CLUMP, tissue, variant and molecular phenotype.
  harmonized_dt[
    ,
    duplicate_record := duplicated(record_key) |
      duplicated(record_key, fromLast = TRUE)
  ]

  duplicate_rows <- harmonized_dt[
    duplicate_record == TRUE
  ]

  # Add source columns using a source_ prefix.
  if (isTRUE(retain_source_columns)) {

    source_copy <- copy(source_dt)

    setnames(
      source_copy,
      old = names(source_copy),
      new = paste0(
        "source_",
        names(source_copy)
      )
    )

    harmonized_dt <- cbind(
      harmonized_dt,
      source_copy
    )
  }

  file_end_time <- Sys.time()

  runtime_seconds <- as.numeric(
    difftime(
      file_end_time,
      file_start_time,
      units = "secs"
    )
  )

  critical_problem_count <- sum(
    is.na(harmonized_dt$variant_id) |
      is.na(harmonized_dt$p_value) |
      harmonized_dt$pvalue_valid == FALSE
  )

  coordinate_problem_count <- sum(
    is.na(harmonized_dt$chromosome) |
      is.na(harmonized_dt$position) |
      harmonized_dt$chromosome_matches_canonical == FALSE |
      harmonized_dt$position_within_canonical_region == FALSE
  )

  file_status <- if (
    nrow(harmonized_dt) == 0L
  ) {

    "FAIL_EMPTY"

  } else if (
    all(is.na(harmonized_dt$p_value))
  ) {

    "FAIL_NO_PVALUES"

  } else if (
    all(
      is.na(harmonized_dt$variant_id) &
        is.na(harmonized_dt$position)
    )
  ) {

    "FAIL_NO_VARIANTS"

  } else if (
    critical_problem_count > 0L ||
      coordinate_problem_count > 0L
  ) {

    "PASS_WITH_ROW_QC_FLAGS"

  } else {

    "PASS"
  }

  file_qc <- data.table(

    queue_id = as.character(
      metadata$queue_id[1L]
    ),

    clump_id = as.character(
      metadata$clump_id[1L]
    ),

    tissue = as.character(
      metadata$tissue[1L]
    ),

    source_file = normalizePath(
      file_path,
      winslash = "/",
      mustWork = TRUE
    ),

    file_size_bytes = file_size_bytes,

    source_rows = nrow(source_dt),
    source_columns = ncol(source_dt),

    unique_variants = uniqueN(
      harmonized_dt$variant_id,
      na.rm = TRUE
    ),

    unique_positions = uniqueN(
      paste(
        harmonized_dt$chromosome,
        harmonized_dt$position,
        sep = ":"
      ),
      na.rm = TRUE
    ),

    unique_phenotypes = uniqueN(
      harmonized_dt$phenotype_id,
      na.rm = TRUE
    ),

    unique_genes = uniqueN(
      harmonized_dt$gene_id,
      na.rm = TRUE
    ),

    missing_variant_id = sum(
      is.na(harmonized_dt$variant_id)
    ),

    missing_rsid = sum(
      is.na(harmonized_dt$rsid)
    ),

    missing_chromosome = sum(
      is.na(harmonized_dt$chromosome)
    ),

    missing_position = sum(
      is.na(harmonized_dt$position)
    ),

    missing_phenotype_id = sum(
      is.na(harmonized_dt$phenotype_id)
    ),

    missing_gene_id = sum(
      is.na(harmonized_dt$gene_id)
    ),

    missing_p_value = sum(
      is.na(harmonized_dt$p_value)
    ),

    invalid_p_value = sum(
      !is.na(harmonized_dt$p_value) &
        harmonized_dt$pvalue_valid == FALSE
    ),

    missing_beta = sum(
      is.na(harmonized_dt$beta)
    ),

    missing_standard_error = sum(
      is.na(harmonized_dt$standard_error)
    ),

    chromosome_mismatch_rows = sum(
      !is.na(harmonized_dt$chromosome) &
        harmonized_dt$chromosome_matches_canonical == FALSE
    ),

    outside_canonical_region_rows = sum(
      !is.na(harmonized_dt$position) &
        harmonized_dt$position_within_canonical_region == FALSE
    ),

    duplicate_rows = nrow(duplicate_rows),

    invalid_maf_rows = sum(
      harmonized_dt$maf_valid == FALSE,
      na.rm = TRUE
    ),

    invalid_allele_frequency_rows = sum(
      harmonized_dt$allele_frequency_valid == FALSE,
      na.rm = TRUE
    ),

    runtime_seconds = round(
      runtime_seconds,
      4
    ),

    file_status = file_status
  )

  list(
    data = harmonized_dt,
    qc = file_qc,
    column_map = column_map_table,
    duplicate_rows = duplicate_rows
  )
}


###############################################################################
# 10. READ WORK QUEUE
###############################################################################

if (!file.exists(work_queue_file)) {

  stop(
    paste0(
      "Module 11.3B work queue not found:\n",
      work_queue_file,
      "\nRun Module 11.3B.1 first."
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

required_work_queue_columns <- c(
  "queue_id",
  "clump_id",
  "tissue",
  "canonical_chromosome",
  "canonical_region_start",
  "canonical_region_end",
  "output_file_resolved"
)

missing_work_queue_columns <- setdiff(
  required_work_queue_columns,
  names(work_queue)
)

if (length(missing_work_queue_columns) > 0L) {

  stop(
    paste0(
      "Required work-queue column(s) missing: ",
      paste(
        missing_work_queue_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

if (nrow(work_queue) == 0L) {
  stop("The work queue is empty.", call. = FALSE)
}

cat("Work queue loaded.\n")
cat("Rows:", nrow(work_queue), "\n")
cat("Unique CLUMPs:", uniqueN(work_queue$clump_id), "\n\n")


###############################################################################
# 11. SELECT TEST FILE
###############################################################################

test_metadata <- work_queue[1L]

test_file <- test_metadata$output_file_resolved[1L]

cat("Test queue entry:\n")
print(
  test_metadata[
    ,
    .(
      queue_id,
      clump_id,
      tissue,
      canonical_chromosome,
      canonical_region_start,
      canonical_region_end,
      output_file_resolved
    )
  ]
)

cat("\n")


###############################################################################
# 12. RUN SINGLE-FILE TEST
###############################################################################

test_result <- tryCatch(

  read_and_qc_apaqtl(
    file_path = test_file,
    metadata = test_metadata,
    retain_source_columns = TRUE,
    fail_on_missing_pvalue = TRUE,
    fail_on_missing_variant = TRUE
  ),

  error = function(e) {

    error_message <- conditionMessage(e)

    cat("\n")
    cat("SINGLE-FILE TEST FAILED\n")
    cat(error_message, "\n\n")

    test_status <- data.table(
      module = "11.3B.2",
      test_queue_id = test_metadata$queue_id,
      test_clump_id = test_metadata$clump_id,
      test_tissue = test_metadata$tissue,
      test_file = test_file,
      final_status = "FAIL",
      error_message = error_message
    )

    fwrite(
      test_status,
      test_status_output,
      sep = "\t",
      quote = FALSE,
      na = "NA"
    )

    stop(
      paste0(
        "Module 11.3B.2 single-file test failed.\n",
        "See:\n",
        test_status_output
      ),
      call. = FALSE
    )
  }
)


###############################################################################
# 13. WRITE TEST OUTPUTS
###############################################################################

fwrite(
  test_result$qc,
  test_file_qc_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  test_result$column_map,
  test_column_map_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

if (nrow(test_result$duplicate_rows) > 0L) {

  fwrite(
    test_result$duplicate_rows,
    test_duplicate_output,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

} else {

  fwrite(
    empty_table_with_columns(
      names(test_result$data)
    ),
    test_duplicate_output,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}

preview_rows <- min(
  100L,
  nrow(test_result$data)
)

fwrite(
  test_result$data[
    seq_len(preview_rows)
  ],
  test_preview_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  test_result$data,
  test_full_output,
  sep = "\t",
  quote = FALSE,
  na = "NA",
  compress = "gzip"
)


###############################################################################
# 14. TEST STATUS
###############################################################################

module_end_time <- Sys.time()

module_runtime_seconds <- as.numeric(
  difftime(
    module_end_time,
    module_start_time,
    units = "secs"
  )
)

test_pass <- test_result$qc$file_status %chin% c(
  "PASS",
  "PASS_WITH_ROW_QC_FLAGS"
)

test_status <- data.table(
  module = "11.3B.2",
  module_name = "single-file apaQTL reader and QC",
  test_queue_id = test_metadata$queue_id,
  test_clump_id = test_metadata$clump_id,
  test_tissue = test_metadata$tissue,
  test_file = test_file,
  source_rows = test_result$qc$source_rows,
  harmonized_rows = nrow(test_result$data),
  detected_source_columns =
    sum(test_result$column_map$detected),
  file_status = test_result$qc$file_status,
  final_status = ifelse(
    test_pass,
    "PASS",
    "FAIL"
  ),
  runtime_seconds = round(
    module_runtime_seconds,
    3
  )
)

fwrite(
  test_status,
  test_status_output,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


###############################################################################
# 15. SAVE SESSION INFORMATION
###############################################################################

writeLines(
  c(
    "PigmentationAtlas Module 11.3B.2",
    paste("Start:", format(module_start_time)),
    paste("End:", format(module_end_time)),
    paste("Final status:", test_status$final_status),
    "",
    capture.output(sessionInfo())
  ),
  con = session_info_output
)


###############################################################################
# 16. FINAL REPORT
###############################################################################

cat("\n")
cat("============================================================\n")
cat("MODULE 11.3B.2 – SINGLE-FILE TEST SUMMARY\n")
cat("============================================================\n\n")

cat("Detected source-column mapping:\n")
print(test_result$column_map)

cat("\nFile-level QC:\n")
print(test_result$qc)

cat("\nFirst harmonized rows:\n")

preview_columns <- intersect(
  c(
    "queue_id",
    "clump_id",
    "tissue",
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
    "chromosome_matches_canonical",
    "position_within_canonical_region",
    "duplicate_record"
  ),
  names(test_result$data)
)

print(
  test_result$data[
    seq_len(min(10L, .N)),
    ..preview_columns
  ]
)

cat("\nOutputs written:\n")
cat("  Test file QC:\n    ", test_file_qc_output, "\n")
cat("  Column mapping:\n    ", test_column_map_output, "\n")
cat("  Duplicate rows:\n    ", test_duplicate_output, "\n")
cat("  Harmonized preview:\n    ", test_preview_output, "\n")
cat("  Harmonized test file:\n    ", test_full_output, "\n")
cat("  Test status:\n    ", test_status_output, "\n")
cat("  Session information:\n    ", session_info_output, "\n\n")

cat("============================================================\n")

if (identical(test_status$final_status, "PASS")) {

  cat("FINAL QC RESULT: PASS\n")

  if (
    identical(
      test_result$qc$file_status,
      "PASS_WITH_ROW_QC_FLAGS"
    )
  ) {

    cat(
      "The reader worked, but one or more rows received QC flags.\n"
    )

    cat(
      "Review the file-level QC before running all 516 files.\n"
    )

  } else {

    cat(
      "The single-file reader and harmonizer worked without row-level QC problems.\n"
    )
  }

  cat("Module 11.3B.3 may be prepared.\n")

} else {

  cat("FINAL QC RESULT: FAIL\n")
  cat("Do not start Module 11.3B.3.\n")
}

cat("============================================================\n\n")

invisible(test_result)