#!/usr/bin/env Rscript

###############################################################################
# PigmentationAtlas
# Module 12.5C – Gene-level Candidate Prioritization
#
# Purpose
#   Collapse record-level colocalization results to one row per clump × gene,
#   while preserving the strongest evidence and summarizing support across
#   tissues/phenotypes. The resulting table is the recommended input for
#   Module 12.6 gene reclassification.
#
# Primary input
#   results/module12/colocalization/biological_annotation/tables/
#   module12_5B_candidate_gene_prioritization_GENCODEv50.tsv.gz
#
# Optional cross-check input
#   results/module12/colocalization/prioritization/tables/
#   module12_4_ranked_coloc_catalogue/
#   module12_4_ranked_coloc_catalogue.tsv
#
# Outputs
#   results/module12/colocalization/biological_annotation/gene_level/
#     module12_5C_gene_level_candidate_prioritization_GENCODEv50.tsv.gz
#     module12_5C_best_gene_per_clump_GENCODEv50.tsv.gz
#     module12_5C_gene_level_summary.tsv
#     module12_5C_priority_class_summary.tsv
#     module12_5C_qc_summary.tsv
#     module12_5C_status.tsv
#     module12_5C_sessionInfo.txt
###############################################################################

options(stringsAsFactors = FALSE, scipen = 999)

# =============================================================================
# 0. PACKAGES
# =============================================================================

required_packages <- c("data.table")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(data.table)
})

# =============================================================================
# 1. HELPERS
# =============================================================================

section <- function(title) {
  cat("\n", paste(rep("=", 78), collapse = ""), "\n", sep = "")
  cat(title, "\n")
  cat(paste(rep("=", 78), collapse = ""), "\n", sep = "")
}

stopf <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

first_existing <- function(x, candidates, required = TRUE) {
  hit <- candidates[candidates %in% names(x)]

  if (length(hit) == 0L) {
    if (required) {
      stopf(
        "None of the expected columns were found: %s",
        paste(candidates, collapse = ", ")
      )
    }
    return(NA_character_)
  }

  hit[1L]
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_sum_logical <- function(x) {
  x <- as.logical(x)
  if (length(x) == 0L || all(is.na(x))) return(NA_integer_)
  sum(x, na.rm = TRUE)
}

collapse_unique <- function(x, sep = ";") {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) return(NA_character_)
  paste(sort(unique(x)), collapse = sep)
}

normalize_gene_id <- function(x) {
  sub("\\.[0-9]+$", "", as.character(x))
}

as_bool <- function(x) {
  if (is.logical(x)) return(x)

  y <- toupper(trimws(as.character(x)))

  out <- rep(NA, length(y))
  out[y %in% c("TRUE", "T", "1", "YES", "Y")] <- TRUE
  out[y %in% c("FALSE", "F", "0", "NO", "N")] <- FALSE
  out
}

priority_rank_from_class <- function(x) {
  ranks <- c(
    "STRONG_H4" = 1L,
    "MODERATE_H4" = 2L,
    "SUGGESTIVE_H4" = 3L,
    "H3_DOMINANT" = 4L,
    "INCONCLUSIVE" = 5L
  )

  unname(ranks[as.character(x)])
}

# =============================================================================
# 2. PATHS
# =============================================================================

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

input_12_5B <- file.path(
  project_root,
  "results",
  "module12",
  "colocalization",
  "biological_annotation",
  "tables",
  "module12_5B_candidate_gene_prioritization_GENCODEv50.tsv.gz"
)

input_12_4 <- file.path(
  project_root,
  "results",
  "module12",
  "colocalization",
  "prioritization",
  "tables",
  "module12_4_ranked_coloc_catalogue",
  "module12_4_ranked_coloc_catalogue.tsv"
)

output_dir <- file.path(
  project_root,
  "results",
  "module12",
  "colocalization",
  "biological_annotation",
  "gene_level"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_gene_level <- file.path(
  output_dir,
  "module12_5C_gene_level_candidate_prioritization_GENCODEv50.tsv.gz"
)

output_best_gene <- file.path(
  output_dir,
  "module12_5C_best_gene_per_clump_GENCODEv50.tsv.gz"
)

output_summary <- file.path(
  output_dir,
  "module12_5C_gene_level_summary.tsv"
)

output_priority_summary <- file.path(
  output_dir,
  "module12_5C_priority_class_summary.tsv"
)

output_qc <- file.path(
  output_dir,
  "module12_5C_qc_summary.tsv"
)

output_status <- file.path(
  output_dir,
  "module12_5C_status.tsv"
)

output_session <- file.path(
  output_dir,
  "module12_5C_sessionInfo.txt"
)

# =============================================================================
# 3. INPUT VALIDATION
# =============================================================================

section("MODULE 12.5C – INPUT VALIDATION")

cat("Project root:\n ", project_root, "\n\n", sep = "")
cat("Module 12.5B input:\n ", input_12_5B, "\n\n", sep = "")
cat("Module 12.4 cross-check input:\n ", input_12_4, "\n\n", sep = "")
cat("Output directory:\n ", output_dir, "\n", sep = "")

if (!file.exists(input_12_5B)) {
  stopf("Required Module 12.5B input file does not exist:\n%s", input_12_5B)
}

# =============================================================================
# 4. READ MODULE 12.5B
# =============================================================================

section("READING MODULE 12.5B")

prior <- fread(input_12_5B)

cat("Rows:    ", nrow(prior), "\n", sep = "")
cat("Columns: ", ncol(prior), "\n", sep = "")

required_core <- c(
  "clump_id",
  "gene_id",
  "priority_class",
  "priority_rank",
  "PP.H3",
  "PP.H4"
)

missing_core <- setdiff(required_core, names(prior))

if (length(missing_core) > 0L) {
  stopf(
    "Required columns are missing from Module 12.5B: %s",
    paste(missing_core, collapse = ", ")
  )
}

# Preserve original versioned gene ID and create a normalized ID for grouping.
prior[, gene_id_original := as.character(gene_id)]
prior[, gene_id := normalize_gene_id(gene_id)]

# Standardize key numeric columns.
numeric_candidates <- intersect(
  c(
    "priority_rank",
    "global_rank",
    "PP.H0",
    "PP.H1",
    "PP.H2",
    "PP.H3",
    "PP.H4",
    "PP.H4_over_H3_H4",
    "h4_minus_h3",
    "nsnps",
    "n_harmonized_variants",
    "top_coloc_variant_PP.H4"
  ),
  names(prior)
)

prior[, (numeric_candidates) := lapply(.SD, as.numeric), .SDcols = numeric_candidates]

# Reconstruct priority rank if needed.
bad_priority_rank <- is.na(prior$priority_rank)

if (any(bad_priority_rank)) {
  replacement <- priority_rank_from_class(prior$priority_class)
  prior[bad_priority_rank, priority_rank := replacement[bad_priority_rank]]
}

if (any(is.na(prior$priority_rank))) {
  warning("Some rows still have missing priority_rank after class-based recovery.")
}

# Locate optional annotation columns adaptively.
gene_symbol_col <- first_existing(
  prior,
  c(
    "gene_symbol",
    "gene_name",
    "external_gene_name",
    "gencode_gene_name",
    "symbol"
  ),
  required = FALSE
)

gene_biotype_col <- first_existing(
  prior,
  c(
    "gene_biotype",
    "gene_type",
    "gencode_gene_type",
    "biotype"
  ),
  required = FALSE
)

biotype_class_col <- first_existing(
  prior,
  c(
    "gene_biotype_class",
    "biotype_class",
    "candidate_biotype_class",
    "broad_biotype_class"
  ),
  required = FALSE
)

chromosome_col <- first_existing(
  prior,
  c("chromosome", "chr", "gene_chromosome", "seqname"),
  required = FALSE
)

gene_start_col <- first_existing(
  prior,
  c("gene_start", "start", "gene_start_bp"),
  required = FALSE
)

gene_end_col <- first_existing(
  prior,
  c("gene_end", "end", "gene_end_bp"),
  required = FALSE
)

gene_strand_col <- first_existing(
  prior,
  c("gene_strand", "strand"),
  required = FALSE
)

cat("Gene symbol column:   ",
    ifelse(is.na(gene_symbol_col), "<not found>", gene_symbol_col),
    "\n",
    sep = "")

cat("Gene biotype column:  ",
    ifelse(is.na(gene_biotype_col), "<not found>", gene_biotype_col),
    "\n",
    sep = "")

cat("Biotype class column: ",
    ifelse(is.na(biotype_class_col), "<not found>", biotype_class_col),
    "\n",
    sep = "")

# =============================================================================
# 5. RECORD-LEVEL QC
# =============================================================================

section("RECORD-LEVEL QC")

record_qc <- prior[
  ,
  .N,
  by = .(clump_id, gene_id)
]

cat("Unique clumps:                ", uniqueN(prior$clump_id), "\n", sep = "")
cat("Unique genes:                 ", uniqueN(prior$gene_id), "\n", sep = "")
cat("Unique clump × gene pairs:    ", nrow(record_qc), "\n", sep = "")
cat("Median records/pair:          ", median(record_qc$N), "\n", sep = "")
cat("Mean records/pair:            ", round(mean(record_qc$N), 2), "\n", sep = "")
cat("Maximum records/pair:         ", max(record_qc$N), "\n", sep = "")

# Optional consistency check against Module 12.4.
crosscheck_12_4 <- data.table(
  metric = character(),
  value = character()
)

if (file.exists(input_12_4)) {
  ranked_12_4 <- fread(
    input_12_4,
    select = intersect(
      c("clump_id", "gene_id", "coloc_unit_id"),
      names(fread(input_12_4, nrows = 0L))
    )
  )

  ranked_12_4[, gene_id := normalize_gene_id(gene_id)]

  crosscheck_12_4 <- data.table(
    metric = c(
      "module12_4_rows",
      "module12_4_unique_clumps",
      "module12_4_unique_clump_gene",
      "module12_5B_rows",
      "module12_5B_unique_clumps",
      "module12_5B_unique_clump_gene"
    ),
    value = as.character(c(
      nrow(ranked_12_4),
      uniqueN(ranked_12_4$clump_id),
      uniqueN(ranked_12_4, by = c("clump_id", "gene_id")),
      nrow(prior),
      uniqueN(prior$clump_id),
      uniqueN(prior, by = c("clump_id", "gene_id"))
    ))
  )

  cat("\nModule 12.4/12.5B structural cross-check completed.\n")
} else {
  cat("\nModule 12.4 file not found; cross-check skipped.\n")
}

# =============================================================================
# 6. SELECT REPRESENTATIVE RECORD PER CLUMP × GENE
# =============================================================================

section("SELECTING REPRESENTATIVE RECORD PER CLUMP × GENE")

# Evidence hierarchy:
#   1. Lower priority_rank is better.
#   2. Higher PP.H4 is better.
#   3. Higher conditional H4 support is better when available.
#   4. Larger H4-H3 difference is better when available.
#   5. Lower global_rank is better when available.
#   6. Stable deterministic tie-breaks.

sort_columns <- c(
  "clump_id",
  "gene_id",
  "priority_rank",
  "PP.H4"
)

if ("PP.H4_over_H3_H4" %in% names(prior)) {
  sort_columns <- c(sort_columns, "PP.H4_over_H3_H4")
}

if ("h4_minus_h3" %in% names(prior)) {
  sort_columns <- c(sort_columns, "h4_minus_h3")
}

if ("global_rank" %in% names(prior)) {
  sort_columns <- c(sort_columns, "global_rank")
}

if ("tissue" %in% names(prior)) {
  sort_columns <- c(sort_columns, "tissue")
}

sort_order <- c(
  1L,   # clump_id ascending
  1L,   # gene_id ascending
  1L,   # priority_rank ascending
  -1L   # PP.H4 descending
)

if ("PP.H4_over_H3_H4" %in% names(prior)) {
  sort_order <- c(sort_order, -1L)
}

if ("h4_minus_h3" %in% names(prior)) {
  sort_order <- c(sort_order, -1L)
}

if ("global_rank" %in% names(prior)) {
  sort_order <- c(sort_order, 1L)
}

if ("tissue" %in% names(prior)) {
  sort_order <- c(sort_order, 1L)
}

setorderv(
  prior,
  cols = sort_columns,
  order = sort_order,
  na.last = TRUE
)

prior[, representative_record := seq_len(.N) == 1L, by = .(clump_id, gene_id)]

representative <- prior[representative_record == TRUE]

if (nrow(representative) != nrow(record_qc)) {
  stopf(
    "Representative record count (%d) does not equal unique clump × gene count (%d).",
    nrow(representative),
    nrow(record_qc)
  )
}

cat("Representative records selected: ", nrow(representative), "\n", sep = "")

# =============================================================================
# 7. BUILD GENE-LEVEL TABLE
# =============================================================================

section("BUILDING GENE-LEVEL TABLE")

# Summaries calculated across all record-level observations for each gene/clump.
support_summary <- prior[
  ,
  .(
    n_coloc_records = .N,
    n_tissues = if ("tissue" %in% names(.SD)) uniqueN(tissue, na.rm = TRUE) else NA_integer_,
    tissues = if ("tissue" %in% names(.SD)) collapse_unique(tissue) else NA_character_,
    n_phenotypes = if ("phenotype_id" %in% names(.SD)) uniqueN(phenotype_id, na.rm = TRUE) else NA_integer_,
    phenotype_ids = if ("phenotype_id" %in% names(.SD)) collapse_unique(phenotype_id) else NA_character_,
    n_strong_h4_records = sum(priority_class == "STRONG_H4", na.rm = TRUE),
    n_moderate_h4_records = sum(priority_class == "MODERATE_H4", na.rm = TRUE),
    n_suggestive_h4_records = sum(priority_class == "SUGGESTIVE_H4", na.rm = TRUE),
    n_h3_dominant_records = sum(priority_class == "H3_DOMINANT", na.rm = TRUE),
    n_inconclusive_records = sum(priority_class == "INCONCLUSIVE", na.rm = TRUE),
    best_priority_rank_across_records = safe_min(priority_rank),
    maximum_PP.H4 = safe_max(PP.H4),
    maximum_PP.H3 = safe_max(PP.H3),
    maximum_PP.H4_over_H3_H4 =
      if ("PP.H4_over_H3_H4" %in% names(.SD)) safe_max(PP.H4_over_H3_H4) else NA_real_,
    maximum_h4_minus_h3 =
      if ("h4_minus_h3" %in% names(.SD)) safe_max(h4_minus_h3) else NA_real_,
    n_high_confidence_h4 =
      if ("high_confidence_h4" %in% names(.SD)) safe_sum_logical(as_bool(high_confidence_h4)) else NA_integer_,
    n_adequate_variant_count =
      if ("adequate_variant_count" %in% names(.SD)) safe_sum_logical(as_bool(adequate_variant_count)) else NA_integer_,
    top_coloc_variants =
      if ("top_coloc_variant" %in% names(.SD)) collapse_unique(top_coloc_variant) else NA_character_
  ),
  by = .(clump_id, gene_id)
]

# Keep all columns from the representative record except the helper flag.
representative[, representative_record := NULL]

gene_level <- merge(
  representative,
  support_summary,
  by = c("clump_id", "gene_id"),
  all.x = TRUE,
  sort = FALSE
)

# Rename representative-record evidence fields for clarity while preserving
# familiar names used by downstream scripts.
rename_map <- c(
  "tissue" = "best_tissue",
  "phenotype_id" = "best_phenotype_id",
  "coloc_unit_id" = "representative_coloc_unit_id",
  "interpretation" = "best_interpretation",
  "evidence_label" = "best_evidence_label"
)

for (old_name in names(rename_map)) {
  new_name <- unname(rename_map[[old_name]])
  if (old_name %in% names(gene_level) && !new_name %in% names(gene_level)) {
    setnames(gene_level, old_name, new_name)
  }
}

# Add a clearly named gene-level eligibility flag.
gene_level[
  ,
  h4_supported_gene := priority_class %in% c(
    "STRONG_H4",
    "MODERATE_H4",
    "SUGGESTIVE_H4"
  )
]

# Deterministic ordering.
gene_level_sort_columns <- c("clump_id", "priority_rank", "PP.H4", "gene_id")
gene_level_sort_order <- c(1L, 1L, -1L, 1L)

if ("PP.H4_over_H3_H4" %in% names(gene_level)) {
  gene_level_sort_columns <- c(
    "clump_id",
    "priority_rank",
    "PP.H4",
    "PP.H4_over_H3_H4",
    "gene_id"
  )
  gene_level_sort_order <- c(1L, 1L, -1L, -1L, 1L)
}

setorderv(
  gene_level,
  cols = gene_level_sort_columns,
  order = gene_level_sort_order,
  na.last = TRUE
)

# Verify uniqueness.
duplicate_pairs <- gene_level[
  ,
  .N,
  by = .(clump_id, gene_id)
][N > 1L]

if (nrow(duplicate_pairs) > 0L) {
  stopf(
    "Gene-level output still contains %d duplicated clump × gene pairs.",
    nrow(duplicate_pairs)
  )
}

cat("Gene-level rows:              ", nrow(gene_level), "\n", sep = "")
cat("Unique clump × gene pairs:    ",
    uniqueN(gene_level, by = c("clump_id", "gene_id")),
    "\n",
    sep = "")
cat("Duplicated pairs:             ", nrow(duplicate_pairs), "\n", sep = "")
cat("H4-supported gene-level rows: ", sum(gene_level$h4_supported_gene), "\n", sep = "")

# =============================================================================
# 8. SELECT BEST H4-SUPPORTED GENE PER CLUMP
# =============================================================================

section("SELECTING BEST H4-SUPPORTED GENE PER CLUMP")

eligible_gene_level <- gene_level[h4_supported_gene == TRUE]

if (nrow(eligible_gene_level) == 0L) {
  warning("No H4-supported gene-level candidates were found.")
  best_gene <- eligible_gene_level[0L]
} else {
  best_sort_columns <- c(
    "clump_id",
    "priority_rank",
    "PP.H4"
  )

  best_sort_order <- c(
    1L,
    1L,
    -1L
  )

  if ("PP.H4_over_H3_H4" %in% names(eligible_gene_level)) {
    best_sort_columns <- c(
      best_sort_columns,
      "PP.H4_over_H3_H4"
    )
    best_sort_order <- c(best_sort_order, -1L)
  }

  if ("h4_minus_h3" %in% names(eligible_gene_level)) {
    best_sort_columns <- c(
      best_sort_columns,
      "h4_minus_h3"
    )
    best_sort_order <- c(best_sort_order, -1L)
  }

  if ("global_rank" %in% names(eligible_gene_level)) {
    best_sort_columns <- c(
      best_sort_columns,
      "global_rank"
    )
    best_sort_order <- c(best_sort_order, 1L)
  }

  best_sort_columns <- c(best_sort_columns, "gene_id")
  best_sort_order <- c(best_sort_order, 1L)

  setorderv(
    eligible_gene_level,
    cols = best_sort_columns,
    order = best_sort_order,
    na.last = TRUE
  )

  best_gene <- eligible_gene_level[, .SD[1L], by = clump_id]
}

if (anyDuplicated(best_gene$clump_id) > 0L) {
  stopf("Best-gene table contains duplicated clump IDs.")
}

cat("Clumps with H4-supported genes: ", uniqueN(eligible_gene_level$clump_id), "\n", sep = "")
cat("Best genes selected:             ", nrow(best_gene), "\n", sep = "")

# Add within-clump gene rank for all gene-level rows.
rank_columns <- c("priority_rank", "PP.H4")
rank_order <- c(1L, -1L)

if ("PP.H4_over_H3_H4" %in% names(gene_level)) {
  rank_columns <- c(rank_columns, "PP.H4_over_H3_H4")
  rank_order <- c(rank_order, -1L)
}

if ("h4_minus_h3" %in% names(gene_level)) {
  rank_columns <- c(rank_columns, "h4_minus_h3")
  rank_order <- c(rank_order, -1L)
}

if ("global_rank" %in% names(gene_level)) {
  rank_columns <- c(rank_columns, "global_rank")
  rank_order <- c(rank_order, 1L)
}

rank_columns <- c(rank_columns, "gene_id")
rank_order <- c(rank_order, 1L)

setorderv(
  gene_level,
  cols = c("clump_id", rank_columns),
  order = c(1L, rank_order),
  na.last = TRUE
)

gene_level[, gene_level_rank_within_clump := seq_len(.N), by = clump_id]

# =============================================================================
# 9. SUMMARY TABLES
# =============================================================================

section("GENERATING SUMMARY TABLES")

gene_level_summary <- data.table(
  metric = c(
    "input_record_rows",
    "input_unique_clumps",
    "input_unique_genes",
    "input_unique_clump_gene_pairs",
    "gene_level_rows",
    "gene_level_unique_clumps",
    "gene_level_unique_genes",
    "gene_level_duplicated_clump_gene_pairs",
    "h4_supported_gene_level_rows",
    "clumps_with_h4_supported_genes",
    "best_genes_selected",
    "median_records_per_clump_gene",
    "mean_records_per_clump_gene",
    "maximum_records_per_clump_gene"
  ),
  value = c(
    nrow(prior),
    uniqueN(prior$clump_id),
    uniqueN(prior$gene_id),
    nrow(record_qc),
    nrow(gene_level),
    uniqueN(gene_level$clump_id),
    uniqueN(gene_level$gene_id),
    nrow(duplicate_pairs),
    sum(gene_level$h4_supported_gene),
    uniqueN(eligible_gene_level$clump_id),
    nrow(best_gene),
    median(record_qc$N),
    round(mean(record_qc$N), 4),
    max(record_qc$N)
  )
)

priority_class_summary <- gene_level[
  ,
  .(
    n_gene_level_candidates = .N,
    n_clumps = uniqueN(clump_id),
    percent_of_gene_level_candidates = round(100 * .N / nrow(gene_level), 2)
  ),
  by = .(priority_rank, priority_class)
][
  order(priority_rank, priority_class)
]

qc_summary <- rbindlist(
  list(
    data.table(
      metric = c(
        "required_columns_present",
        "gene_level_pairs_unique",
        "best_gene_clumps_unique",
        "all_best_genes_h4_supported",
        "all_priority_ranks_nonmissing",
        "all_PP_H4_in_unit_interval"
      ),
      value = as.character(c(
        length(missing_core) == 0L,
        nrow(duplicate_pairs) == 0L,
        anyDuplicated(best_gene$clump_id) == 0L,
        all(best_gene$h4_supported_gene),
        all(!is.na(gene_level$priority_rank)),
        all(
          is.na(gene_level$PP.H4) |
            (gene_level$PP.H4 >= 0 & gene_level$PP.H4 <= 1)
        )
      ))
    ),
    crosscheck_12_4
  ),
  use.names = TRUE,
  fill = TRUE
)

module_status <- data.table(
  module = "12.5C",
  module_name = "Gene-level Candidate Prioritization",
  status = "COMPLETED",
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  input_file = input_12_5B,
  output_gene_level = output_gene_level,
  output_best_gene = output_best_gene,
  input_rows = nrow(prior),
  gene_level_rows = nrow(gene_level),
  clumps_with_h4_supported_genes = uniqueN(eligible_gene_level$clump_id),
  best_genes_selected = nrow(best_gene)
)

# =============================================================================
# 10. WRITE OUTPUTS
# =============================================================================

section("WRITING OUTPUTS")

fwrite(
  gene_level,
  output_gene_level,
  sep = "\t",
  quote = FALSE,
  na = "NA",
  compress = "gzip"
)

fwrite(
  best_gene,
  output_best_gene,
  sep = "\t",
  quote = FALSE,
  na = "NA",
  compress = "gzip"
)

fwrite(
  gene_level_summary,
  output_summary,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  priority_class_summary,
  output_priority_summary,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  qc_summary,
  output_qc,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  module_status,
  output_status,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

capture.output(sessionInfo(), file = output_session)

cat("Gene-level table:\n ", output_gene_level, "\n\n", sep = "")
cat("Best gene per clump:\n ", output_best_gene, "\n\n", sep = "")
cat("Summary:\n ", output_summary, "\n\n", sep = "")
cat("Priority summary:\n ", output_priority_summary, "\n\n", sep = "")
cat("QC summary:\n ", output_qc, "\n\n", sep = "")
cat("Status:\n ", output_status, "\n\n", sep = "")
cat("Session info:\n ", output_session, "\n", sep = "")

# =============================================================================
# 11. FINAL REPORT
# =============================================================================

section("MODULE 12.5C COMPLETED")

cat("Input record rows:             ", nrow(prior), "\n", sep = "")
cat("Unique clump × gene pairs:     ", nrow(record_qc), "\n", sep = "")
cat("Gene-level output rows:        ", nrow(gene_level), "\n", sep = "")
cat("Duplicated gene-level pairs:   ", nrow(duplicate_pairs), "\n", sep = "")
cat("H4-supported gene candidates:  ", sum(gene_level$h4_supported_gene), "\n", sep = "")
cat("Clumps with H4 support:        ", uniqueN(eligible_gene_level$clump_id), "\n", sep = "")
cat("Best genes selected:           ", nrow(best_gene), "\n", sep = "")

cat("\nPriority classes at gene level:\n")
print(priority_class_summary)

cat("\nRecommended Module 12.6 input:\n")
cat(output_best_gene, "\n")
