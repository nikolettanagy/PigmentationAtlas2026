#!/usr/bin/env Rscript
# Compare original and isolated 12.6 outputs by CLUMP ID.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply PigmentationAtlas project root")
suppressPackageStartupMessages(library(data.table))
root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
old_path <- file.path(root, "results/module12/reclassification",
                      "module12_6_clump_gene_reclassification.tsv.gz")
new_dir <- file.path(root, "results/module12_rebuild_staging/reclassification")
new_path <- file.path(new_dir, "module12_6_clump_gene_reclassification.tsv.gz")
old <- fread(old_path)
new <- fread(new_path)
if (nrow(old) != 1325L || nrow(new) != 1325L ||
    anyDuplicated(old$clump_id) || anyDuplicated(new$clump_id) ||
    !setequal(old$clump_id, new$clump_id)) stop("Clump master changed")
setkey(old, clump_id)
setkey(new, clump_id)
old <- old[new$clump_id]
norm <- function(z) {z <- as.character(z); z[is.na(z)] <- "<missing>"; z}
if (!identical(norm(old$baseline_gene_id), norm(new$baseline_gene_id)))
  stop("Conventional baseline gene assignment changed")
result <- data.table(
  clump_id = new$clump_id,
  old_evaluable = old$evaluable,
  new_evaluable = new$evaluable,
  old_reclassified = old$interpretation_change,
  new_reclassified = new$interpretation_change,
  baseline_gene_id = new$baseline_gene_id,
  old_prioritized_gene_id = old$prioritized_gene_id,
  new_prioritized_gene_id = new$prioritized_gene_id,
  old_priority_class = old$priority_class,
  new_priority_class = new$priority_class,
  old_PP.H4 = old$PP.H4,
  new_PP.H4 = new$PP.H4)
result[, change := fifelse(!old_evaluable & new_evaluable, "NEWLY_EVALUABLE",
                   fifelse(old_evaluable & !new_evaluable, "LOST_EVALUABILITY",
                   fifelse(old_evaluable & new_evaluable &
                           norm(old_prioritized_gene_id) != norm(new_prioritized_gene_id),
                           "BEST_GENE_CHANGED",
                   fifelse(old_evaluable & new_evaluable,
                           "SAME_BEST_GENE", "NEITHER_EVALUABLE"))))]
summary <- result[, .(clumps = .N,
                      old_reclassified = sum(old_reclassified == TRUE, na.rm = TRUE),
                      new_reclassified = sum(new_reclassified == TRUE, na.rm = TRUE)),
                  by = change][order(change)]
fwrite(result, file.path(new_dir, "old_vs_staged_reclassification_by_clump.tsv"), sep = "\t")
fwrite(summary, file.path(new_dir, "old_vs_staged_reclassification_summary.tsv"), sep = "\t")
cat("Old: evaluable", sum(old$evaluable), "reclassified",
    sum(old$interpretation_change == TRUE, na.rm = TRUE), "\n")
cat("New: evaluable", sum(new$evaluable), "reclassified",
    sum(new$interpretation_change == TRUE, na.rm = TRUE), "\n")
print(summary)
cat("Details:", file.path(new_dir, "old_vs_staged_reclassification_by_clump.tsv"), "\n")
