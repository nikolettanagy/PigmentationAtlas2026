#!/usr/bin/env Rscript
# Audit 9471 coloc units recalculated with one GWAS case fraction; compare 7840 legacy units.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply PigmentationAtlas project root")
suppressPackageStartupMessages(library(data.table))
root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
base <- file.path(root, "results/module12_rebuild_staging/colocalization")
queue <- fread(file.path(base, "work_queue/module12_1_coloc_work_queue_eligible.tsv.gz"))
oldbase <- file.path(root, "results/module12/colocalization")
oldq <- fread(file.path(oldbase, "work_queue/module12_1_coloc_work_queue_eligible.tsv.gz"))
oldsum <- fread(file.path(oldbase,
  "catalogue/summary/module12_3_coloc_catalogue_summary.tsv.gz"))
stopifnot(nrow(queue) == 9471L, nrow(oldq) == 7840L,
          nrow(oldsum) == 7840L,
          !anyDuplicated(queue$coloc_unit_id), !anyDuplicated(oldq$coloc_unit_id),
          !anyDuplicated(oldsum$coloc_unit_id))
key <- "clump_tissue_phenotype_gene_key"
if (anyDuplicated(queue[[key]]) || anyDuplicated(oldq[[key]]) ||
    !setequal(oldq$coloc_unit_id, oldsum$coloc_unit_id))
  stop("Original or staged identity mapping is inconsistent")
old_index <- match(queue[[key]], oldq[[key]])
old_sum_index <- match(oldq$coloc_unit_id, oldsum$coloc_unit_id)
old_h4 <- rep(NA_real_, nrow(queue))
old_h4[!is.na(old_index)] <- oldsum$PP.H4[old_sum_index[old_index[!is.na(old_index)]]]
new_audit <- fread(file.path(base, "new_units_audit.tsv"))
remaining <- fread(file.path(base, "remaining_units_run_status.tsv"))
if (nrow(new_audit) != 1631L || any(!new_audit$result_verified) ||
    nrow(remaining) != 7840L || anyDuplicated(remaining$coloc_unit_id) ||
    any(!remaining$status %in% c("PASS", "PASS_WITH_WARNINGS")))
  stop("Run the 1631 + 7840 unit audits successfully before aggregating")
new_summary <- vector("list", nrow(queue))
new_qc <- vector("list", nrow(queue))
for (i in seq_len(nrow(queue))) {
  id <- queue$coloc_unit_id[[i]]
  clump <- queue$clump_id[[i]]
  dir <- file.path(base, "single_test/results", paste0(id, "_", clump))
  path <- file.path(dir, paste0(id, "_coloc_summary.tsv"))
  qcpath <- file.path(base, "single_test/qc", paste0(id, "_qc.tsv"))
  statuspath <- file.path(dir, paste0(id, "_status.tsv"))
  if (!file.exists(path) || !file.exists(qcpath) || !file.exists(statuspath))
    stop("Unit outputs missing: ", id)
  s <- fread(path)
  q <- fread(qcpath)
  st <- fread(statuspath)
  if (nrow(s) != 1L || nrow(q) != 1L || nrow(st) != 1L ||
      s$coloc_unit_id[[1L]] != id || q$coloc_unit_id[[1L]] != id ||
      !st$final_status[[1L]] %in% c("PASS", "PASS_WITH_WARNINGS") ||
      s$clump_id[[1L]] != clump ||
      s$gene_id[[1L]] != queue$gene_id[[i]] ||
      s$phenotype_id[[1L]] != queue$phenotype_id[[i]] ||
      !is.finite(s$PP.H4[[1L]]) || s$nsnps[[1L]] < 50L)
    stop("Unit result or biological identity invalid: ", id)
  new_summary[[i]] <- s
  new_qc[[i]] <- q
  if (i %% 1000L == 0L) cat("Checked", i, "of", nrow(queue), "\n")
}
summary <- rbindlist(new_summary, use.names = TRUE, fill = TRUE)
qc <- rbindlist(new_qc, use.names = TRUE, fill = TRUE)
if (nrow(summary) != 9471L || anyDuplicated(summary$coloc_unit_id))
  stop("Unreconciled full catalogue")
old_class <- function(h4, h3) fifelse(h4 >= .8, "STRONG_H4",
  fifelse(h4 >= .5, "MODERATE_H4", fifelse(h4 >= .2, "SUGGESTIVE_H4",
  fifelse(h3 >= .8, "H3_DOMINANT", "INCONCLUSIVE"))))
old_h3 <- rep(NA_real_, nrow(queue))
old_h3[!is.na(old_index)] <- oldsum$PP.H3[old_sum_index[old_index[!is.na(old_index)]]]
comparison <- data.table(coloc_unit_id = queue$coloc_unit_id,
  clump_id = queue$clump_id, biological_key = queue[[key]],
  previously_eligible = !is.na(old_index), old_PP.H4 = old_h4,
  new_PP.H4 = as.numeric(summary$PP.H4),
  old_class = old_class(old_h4, old_h3),
  new_class = old_class(as.numeric(summary$PP.H4), as.numeric(summary$PP.H3)))
comparison[previously_eligible == FALSE, old_class := NA_character_]
comparison[, delta_PP.H4 := new_PP.H4 - old_PP.H4]
out <- file.path(base, "unified_catalogue/summary")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(summary, file.path(out, "module12_3_coloc_catalogue_summary.tsv.gz"),
       sep = "\t", compress = "gzip", quote = FALSE, na = "NA")
fwrite(qc, file.path(out, "module12_3_coloc_catalogue_qc.tsv.gz"),
       sep = "\t", compress = "gzip", quote = FALSE, na = "NA")
fwrite(comparison, file.path(out, "old_vs_unified_by_biological_key.tsv.gz"),
       sep = "\t", compress = "gzip", quote = FALSE, na = "NA")
legacy <- comparison[previously_eligible == TRUE]
cat("PASS: 9471 unified units (7840 recomputed + 1631 prior verified).\n")
cat("Former 7840: maximum absolute PP.H4 difference:",
    max(abs(legacy$delta_PP.H4)), "median:", median(abs(legacy$delta_PP.H4)), "\n")
cat("Former 7840 with changed priority class:",
    sum(legacy$old_class != legacy$new_class), "\n")
cat("Unified strong H4:", sum(summary$PP.H4 >= .8), "moderate H4:",
    sum(summary$PP.H4 >= .5 & summary$PP.H4 < .8), "\n")
cat("Full results:", out, "\n")
