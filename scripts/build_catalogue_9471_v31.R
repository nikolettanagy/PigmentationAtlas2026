#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("Missing ", flag, call. = FALSE)
  args[[i + 1L]]
}
root <- normalizePath(getarg("--project"), winslash = "/", mustWork = TRUE)
output <- normalizePath(getarg("--output"), winslash = "/", mustWork = FALSE)
if (file.exists(output) || dir.exists(output))
  stop("Output already exists; choose a fresh directory.", call. = FALSE)
stage_root <- file.path(root, "results/module12_rebuild_staging/colocalization")
queue_path <- file.path(stage_root, "work_queue/module12_1_coloc_work_queue_eligible.tsv.gz")
staged_results <- file.path(stage_root, "single_test/results")
original_units <- file.path(root, "results/module12/colocalization/catalogue/units")
if (!file.exists(queue_path) || !dir.exists(staged_results) || !dir.exists(original_units))
  stop("Required staging queue, staged results, or original per-unit directory missing.", call. = FALSE)
queue <- fread(queue_path, showProgress = FALSE)
if (!"coloc_unit_id" %in% names(queue)) stop("Queue lacks coloc_unit_id.", call. = FALSE)
ids <- as.character(queue$coloc_unit_id)
if (length(ids) != 9471L || uniqueN(ids) != 9471L || anyNA(ids))
  stop("Staging queue must contain exactly 9471 unique IDs.", call. = FALSE)
if (any(!grepl("^COLOC_[A-Za-z0-9_]+$", ids)))
  stop("Unexpected coloc_unit_id in queue.", call. = FALSE)

# Locate per-unit results only. Existing aggregate catalogue files are never read.
staged <- list.files(staged_results, pattern = "_coloc_summary\\.tsv$", recursive = TRUE,
                     full.names = TRUE)
staged_id <- sub("_coloc_summary\\.tsv$", "", basename(staged))
if (anyDuplicated(staged_id[staged_id %in% ids]))
  stop("More than one staged summary for a queued ID; resolve ambiguity first.", call. = FALSE)
staged_map <- setNames(staged, staged_id)
old_paths <- file.path(original_units, ids, paste0(ids, "_coloc_summary.tsv"))
names(old_paths) <- ids
selected <- staged_map[ids]
fallback <- is.na(selected)
selected[fallback] <- old_paths[ids[fallback]]
missing <- ids[!file.exists(selected)]
if (length(missing)) stop("Missing per-unit summaries (first IDs): ",
                           paste(head(missing, 20L), collapse = ", "), call. = FALSE)

read_one <- function(path, id) {
  x <- fread(path, showProgress = FALSE)
  if (nrow(x) != 1L) stop("Expected one row in ", path, call. = FALSE)
  if ("coloc_unit_id" %in% names(x)) {
    if (is.na(x$coloc_unit_id[[1L]]) || as.character(x$coloc_unit_id[[1L]]) != id)
      stop("ID mismatch: ", path, call. = FALSE)
  } else x[, coloc_unit_id := id]
  x
}
records <- lapply(seq_along(ids), function(i) read_one(selected[[i]], ids[[i]]))
columns <- names(records[[1L]])
bad <- which(!vapply(records, function(x) setequal(names(x), columns), logical(1)))
if (length(bad)) stop("Different summary columns in ", selected[[bad[[1L]]]], call. = FALSE)
for (i in seq_along(records)) setcolorder(records[[i]], columns)
summary <- rbindlist(records, use.names = TRUE)
setorder(summary, coloc_unit_id)
if (nrow(summary) != 9471L || uniqueN(summary$coloc_unit_id) != 9471L)
  stop("Catalogue row or ID count mismatch.", call. = FALSE)
pp_candidates <- list(paste0("PP.H", 0:4), paste0("PP.H", 0:4, ".abf"))
pp_index <- which(vapply(pp_candidates, function(x) all(x %in% names(summary)), logical(1)))
if (!length(pp_index)) stop("Missing complete H0–H4 posterior columns.", call. = FALSE)
pp <- pp_candidates[[pp_index[[1L]]]]
if (all(pp %in% names(summary))) {
  p <- as.matrix(summary[, ..pp])
  storage.mode(p) <- "double"
  if (anyNA(p) || any(!is.finite(p)) || any(p < 0 | p > 1) ||
      any(abs(rowSums(p) - 1) > 1e-6))
    stop("Invalid posterior vector.", call. = FALSE)
}

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
tmp <- tempfile(".catalogue_9471_v31_", tmpdir = dirname(output))
dir.create(tmp)
fwrite(summary, file.path(tmp, "module12_3_coloc_catalogue_summary.tsv.gz"),
       sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
sources <- data.table(coloc_unit_id = ids, source = ifelse(fallback, "original_unit", "staged_unit"),
                      source_path = unname(selected),
                      source_md5 = unname(tools::md5sum(selected)))
fwrite(sources, file.path(tmp, "source_manifest.tsv.gz"), sep = "\t", compress = "gzip")
interpretation <- intersect(c("interpretation", "coloc_interpretation"), names(summary))
if (length(interpretation)) {
  labels <- summary[[interpretation[[1L]]]]
  counts <- as.data.table(table(labels, useNA = "ifany"))
  setnames(counts, c("interpretation", "N"))
  setorder(counts, -N, interpretation)
  fwrite(counts, file.path(tmp, "interpretation_counts.tsv"), sep = "\t")
}
audit <- data.table(result = "PASS", queue_rows = length(ids),
                    independent_summary_rows = nrow(summary),
                    staged_units = sum(!fallback), original_per_unit_fallback = sum(fallback),
                    old_aggregate_read = FALSE,
                    posterior_checked = TRUE)
fwrite(audit, file.path(tmp, "audit.tsv"), sep = "\t")
if (!file.rename(tmp, output))
  stop("Could not finalize output; staged files remain at ", tmp, call. = FALSE)
print(audit)
cat("Output: ", output, "\n", sep = "")
