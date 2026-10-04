#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

cat("============================================================\n")
cat("PigmentationAtlas\n")
cat("MODULE 12.5 – BIOLOGICAL ANNOTATION & CANDIDATE PRIORITIZATION\n")
cat("============================================================\n\n")

start_time <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply PigmentationAtlas project root")
project_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)

input_file <- file.path(
  project_dir,
  "results", "module12_rebuild_staging", "colocalization", "prioritization",
  "tables", "module12_4_ranked_coloc_catalogue.tsv.gz"
)

outdir <- file.path(
  project_dir,
  "results", "module12_rebuild_staging", "colocalization", "biological_annotation"
)
tables_dir <- file.path(outdir, "tables")
logs_dir <- file.path(outdir, "logs")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  logs_dir,
  paste0(
    "module12_5_",
    format(start_time, "%Y%m%d_%H%M%S"),
    ".log"
  )
)

log_message <- function(...) {
  msg <- paste0(..., collapse = "")
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

fail <- function(msg) {
  log_message("FATAL: ", msg)
  cat("\nFINAL RESULT: FAIL\n")
  quit(save = "no", status = 1L)
}

if (!file.exists(input_file)) {
  fail(paste0("Missing input: ", input_file))
}

dt <- fread(input_file, showProgress = FALSE)

required <- c(
  "coloc_unit_id", "clump_id", "tissue", "phenotype_id",
  "gene_id", "nsnps", "PP.H3", "PP.H4",
  "priority_class", "global_rank"
)

missing_required <- setdiff(required, names(dt))
if (length(missing_required) > 0L) {
  fail(
    paste0(
      "Missing required columns: ",
      paste(missing_required, collapse = ", ")
    )
  )
}

log_message("Input coloc units: ", nrow(dt))
log_message("Unique genes: ", uniqueN(dt$gene_id))

# ------------------------------------------------------------------
# 1. Stable identifiers and safe annotation placeholders
# ------------------------------------------------------------------

dt[, gene_id_versioned := as.character(gene_id)]
dt[, gene_id_stable := sub("\\.[0-9]+$", "", gene_id_versioned)]

# Use existing annotation columns if they are ever added upstream.
if ("gene_name" %in% names(dt)) {
  dt[, gene_symbol := as.character(gene_name)]
} else if ("gene_symbol" %in% names(dt)) {
  dt[, gene_symbol := as.character(gene_symbol)]
} else if ("gene" %in% names(dt)) {
  dt[, gene_symbol := as.character(gene)]
} else {
  dt[, gene_symbol := rep(NA_character_, .N)]
}

if (!"gene_biotype" %in% names(dt)) {
  dt[, gene_biotype := rep(NA_character_, .N)]
} else {
  dt[, gene_biotype := as.character(gene_biotype)]
}

# ------------------------------------------------------------------
# 2. Candidate class
# ------------------------------------------------------------------

dt[, candidate_class := fifelse(
  is.na(gene_biotype) | gene_biotype == "",
  "Unknown",
  fifelse(
    gene_biotype == "protein_coding",
    "Protein_coding",
    fifelse(
      grepl(
        "lncRNA|lincRNA|antisense|sense_intronic|sense_overlapping|processed_transcript|macro_lncRNA|bidirectional_promoter_lncRNA",
        gene_biotype,
        ignore.case = TRUE
      ),
      "Regulatory_RNA",
      fifelse(
        grepl("pseudogene", gene_biotype, ignore.case = TRUE),
        "Pseudogene",
        "Other"
      )
    )
  )
)]

# ------------------------------------------------------------------
# 3. Evidence and candidate score
# ------------------------------------------------------------------

dt[, evidence_level := fifelse(
  PP.H4 >= 0.80, "Strong",
  fifelse(
    PP.H4 >= 0.50, "Moderate",
    fifelse(PP.H4 >= 0.20, "Suggestive", "Low")
  )
)]

dt[, candidate_score :=
     100 * PP.H4 +
     10 * fifelse(PP.H4 > PP.H3, 1, 0) +
     5 * fifelse(nsnps >= 100, 1, 0) +
     5 * fifelse(nsnps >= 500, 1, 0) +
     5 * fifelse(candidate_class == "Regulatory_RNA", 1, 0)
]

dt[, candidate_score := round(candidate_score, 3)]

setorder(
  dt,
  priority_rank,
  -candidate_score,
  -PP.H4,
  global_rank
)

# ------------------------------------------------------------------
# 4. Full candidate catalogue
# ------------------------------------------------------------------

candidate_file <- file.path(
  tables_dir,
  "module12_5_candidate_gene_prioritization.tsv.gz"
)

fwrite(
  dt,
  candidate_file,
  sep = "\t",
  quote = FALSE,
  na = "NA",
  compress = "gzip"
)

# ------------------------------------------------------------------
# 5. Gene-level summary
# ------------------------------------------------------------------

safe_max <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

safe_min <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

gene_summary <- dt[
  ,
  .(
    gene_symbol = {
      x <- unique(na.omit(gene_symbol))
      if (length(x) == 0L) NA_character_ else x[[1L]]
    },
    gene_biotype = {
      x <- unique(na.omit(gene_biotype))
      if (length(x) == 0L) NA_character_ else x[[1L]]
    },
    candidate_class = {
      x <- unique(na.omit(candidate_class))
      if (length(x) == 0L) "Unknown" else x[[1L]]
    },
    n_coloc_units = .N,
    n_clumps = uniqueN(clump_id),
    n_tissues = uniqueN(tissue),
    n_phenotypes = uniqueN(phenotype_id),
    best_PP.H4 = safe_max(PP.H4),
    best_PP.H3 = safe_max(PP.H3),
    best_candidate_score = safe_max(candidate_score),
    best_priority_rank = safe_min(priority_rank),
    best_priority_class = priority_class[which.min(priority_rank)][1L],
    best_coloc_unit_id = coloc_unit_id[
      order(priority_rank, -candidate_score, -PP.H4)
    ][1L],
    best_clump_id = clump_id[
      order(priority_rank, -candidate_score, -PP.H4)
    ][1L],
    best_tissue = tissue[
      order(priority_rank, -candidate_score, -PP.H4)
    ][1L],
    best_phenotype_id = phenotype_id[
      order(priority_rank, -candidate_score, -PP.H4)
    ][1L]
  ),
  by = .(gene_id_stable)
]

setorder(
  gene_summary,
  best_priority_rank,
  -best_candidate_score,
  -best_PP.H4
)

gene_summary[, gene_rank := seq_len(.N)]
setcolorder(
  gene_summary,
  c("gene_rank", setdiff(names(gene_summary), "gene_rank"))
)

gene_summary_file <- file.path(
  tables_dir,
  "module12_5_gene_priority_summary.tsv"
)

fwrite(
  gene_summary,
  gene_summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------------
# 6. Manuscript-focused subsets
# ------------------------------------------------------------------

strong_units <- dt[priority_class == "STRONG_H4"]
moderate_units <- dt[priority_class == "MODERATE_H4"]
suggestive_units <- dt[priority_class == "SUGGESTIVE_H4"]

fwrite(
  strong_units,
  file.path(tables_dir, "module12_5_strong_H4_candidates.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

fwrite(
  dt[priority_class %chin% c("STRONG_H4", "MODERATE_H4")],
  file.path(tables_dir, "module12_5_strong_and_moderate_candidates.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

top_gene_candidates <- gene_summary[
  best_priority_class %chin% c("STRONG_H4", "MODERATE_H4")
]

fwrite(
  top_gene_candidates,
  file.path(tables_dir, "module12_5_top_gene_candidates.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

# ------------------------------------------------------------------
# 7. Summary tables
# ------------------------------------------------------------------

candidate_class_counts <- dt[
  ,
  .(
    n_units = .N,
    n_genes = uniqueN(gene_id_stable),
    n_clumps = uniqueN(clump_id),
    median_PP.H4 = median(PP.H4, na.rm = TRUE),
    max_PP.H4 = max(PP.H4, na.rm = TRUE)
  ),
  by = .(priority_class, candidate_class)
][order(priority_class, candidate_class)]

fwrite(
  candidate_class_counts,
  file.path(tables_dir, "module12_5_candidate_class_counts.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

annotation_completeness <- data.table(
  metric = c(
    "input_coloc_units",
    "unique_genes",
    "units_with_gene_symbol",
    "units_with_gene_biotype",
    "genes_with_gene_symbol",
    "genes_with_gene_biotype"
  ),
  value = c(
    nrow(dt),
    uniqueN(dt$gene_id_stable),
    sum(!is.na(dt$gene_symbol) & dt$gene_symbol != ""),
    sum(!is.na(dt$gene_biotype) & dt$gene_biotype != ""),
    uniqueN(dt[!is.na(gene_symbol) & gene_symbol != "", gene_id_stable]),
    uniqueN(dt[!is.na(gene_biotype) & gene_biotype != "", gene_id_stable])
  )
)

fwrite(
  annotation_completeness,
  file.path(tables_dir, "module12_5_annotation_completeness.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

# ------------------------------------------------------------------
# 8. QC and final status
# ------------------------------------------------------------------

critical_failures <- 0L
warnings <- 0L

if (nrow(dt) != 7840L) warnings <- warnings + 1L
if (anyDuplicated(dt$coloc_unit_id) > 0L) critical_failures <- critical_failures + 1L
if (nrow(strong_units) != 25L) warnings <- warnings + 1L
if (all(is.na(dt$gene_biotype))) warnings <- warnings + 1L

final_status <- if (critical_failures > 0L) {
  "FAIL"
} else if (warnings > 0L) {
  "PASS_WITH_WARNINGS"
} else {
  "PASS"
}

end_time <- Sys.time()

status <- data.table(
  module = "12.5",
  start_time = format(start_time, "%Y-%m-%d %H:%M:%S"),
  end_time = format(end_time, "%Y-%m-%d %H:%M:%S"),
  runtime_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),
  input_coloc_units = nrow(dt),
  unique_genes = uniqueN(dt$gene_id_stable),
  strong_h4_units = nrow(strong_units),
  moderate_h4_units = nrow(moderate_units),
  suggestive_h4_units = nrow(suggestive_units),
  annotated_gene_symbols = uniqueN(
    dt[!is.na(gene_symbol) & gene_symbol != "", gene_id_stable]
  ),
  annotated_gene_biotypes = uniqueN(
    dt[!is.na(gene_biotype) & gene_biotype != "", gene_id_stable]
  ),
  critical_failures = critical_failures,
  warnings = warnings,
  final_status = final_status
)

fwrite(
  status,
  file.path(tables_dir, "module12_5_status.tsv"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("\n============================================================\n")
cat("MODULE 12.5 – FINAL SUMMARY\n")
cat("============================================================\n\n")
print(status)

cat("\nMain outputs:\n")
cat("  ", candidate_file, "\n", sep = "")
cat("  ", gene_summary_file, "\n", sep = "")
cat("  ", file.path(tables_dir, "module12_5_strong_H4_candidates.tsv"), "\n", sep = "")
cat("  ", file.path(tables_dir, "module12_5_top_gene_candidates.tsv"), "\n", sep = "")

if (all(is.na(dt$gene_biotype))) {
  cat(
    "\nNOTE: gene biotype and symbol annotation are not present in the\n",
    "Module 12.4 catalogue. Statistical prioritization was completed,\n",
    "but GENCODE annotation must be merged in a subsequent annotation step.\n",
    sep = ""
  )
}

cat("\n============================================================\n")
cat("FINAL RESULT: ", final_status, "\n", sep = "")
cat("============================================================\n")

quit(
  save = "no",
  status = if (identical(final_status, "FAIL")) 1L else 0L
)
