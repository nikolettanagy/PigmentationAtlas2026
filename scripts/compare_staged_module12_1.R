#!/usr/bin/env Rscript
# Usage: Rscript scripts/compare_staged_module12_1.R PROJECT_ROOT
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the full PigmentationAtlas project root.")
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
old_path <- file.path(root, "results/module12/colocalization/work_queue/module12_1_coloc_work_queue.tsv.gz")
new_dir <- file.path(root, "results/module12_rebuild_staging/colocalization/work_queue")
new_path <- file.path(new_dir, "module12_1_coloc_work_queue.tsv.gz")
if (!file.exists(old_path) || !file.exists(new_path)) stop("Old or staged work queue is missing.")
old <- data.table::fread(old_path)
new <- data.table::fread(new_path)
if (nrow(old) != 9584L || sum(old$queue_status == "ELIGIBLE") != 7840L)
  stop("Original queue does not match the audited 9584 / 7840 baseline.")
if (anyNA(old$queue_status) || anyNA(new$queue_status) ||
    any(!old$queue_status %in% c("ELIGIBLE", "EXCLUDED")) ||
    any(!new$queue_status %in% c("ELIGIBLE", "EXCLUDED")))
  stop("Unexpected queue status in original or staged work queue.")
key <- "clump_tissue_phenotype_gene_key"
stopifnot(key %in% names(old), key %in% names(new),
          !anyDuplicated(old[[key]]), !anyDuplicated(new[[key]]))
old_keys <- old[[key]]; new_keys <- new[[key]]
shared <- intersect(old_keys, new_keys)
old_shared <- old[match(shared, old_keys)]
new_shared <- new[match(shared, new_keys)]
new_only <- new[!get(key) %in% old_keys]
old_only <- old[!get(key) %in% new_keys]
newly_eligible <- new_only[queue_status == "ELIGIBLE"]
lost_eligible <- old_only[queue_status == "ELIGIBLE"]
changed_status <- sum(old_shared$queue_status != new_shared$queue_status)
shared_newly_eligible <- sum(old_shared$queue_status == "EXCLUDED" &
                             new_shared$queue_status == "ELIGIBLE")
shared_lost_eligible <- sum(old_shared$queue_status == "ELIGIBLE" &
                            new_shared$queue_status == "EXCLUDED")
old_counts <- table(factor(old$queue_status, levels = c("ELIGIBLE", "EXCLUDED")))
new_counts <- table(factor(new$queue_status, levels = c("ELIGIBLE", "EXCLUDED")))
cat("Old queue:", nrow(old), "total,", old_counts[[1L]], "eligible,",
    old_counts[[2L]], "excluded\n")
cat("New queue:", nrow(new), "total,", new_counts[[1L]], "eligible,",
    new_counts[[2L]], "excluded\n")
cat("Shared unit keys:", length(shared), "changed eligibility:", changed_status, "\n")
cat("Shared newly eligible:", shared_newly_eligible,
    "shared lost eligibility:", shared_lost_eligible, "\n")
cat("New-only unit keys:", nrow(new_only), "of which eligible:", nrow(newly_eligible), "\n")
cat("Old-only unit keys:", nrow(old_only), "of which formerly eligible:", nrow(lost_eligible), "\n")
report <- data.table::data.table(metric = c(
  "old_total", "old_eligible", "old_excluded", "new_total", "new_eligible",
  "new_excluded", "shared_keys", "shared_changed_status", "new_only",
  "new_only_eligible", "old_only", "old_only_eligible",
  "shared_newly_eligible", "shared_lost_eligible"),
  value = c(nrow(old), old_counts[[1L]], old_counts[[2L]], nrow(new),
            new_counts[[1L]], new_counts[[2L]], length(shared), changed_status,
            nrow(new_only), nrow(newly_eligible), nrow(old_only), nrow(lost_eligible),
            shared_newly_eligible, shared_lost_eligible))
data.table::fwrite(report, file.path(new_dir, "old_vs_staged_summary.tsv"), sep = "\t")
data.table::fwrite(newly_eligible[, .(source_queue_id, clump_id, tissue,
                                       phenotype_id, gene_id, effective_shared_variants)],
                   file.path(new_dir, "newly_eligible_units.tsv"), sep = "\t")
if (nrow(lost_eligible) || shared_lost_eligible) {
  cat("REVIEW: some previously eligible units disappeared or changed status.\n")
} else {
  cat("Previously eligible unit keys retained their status.\n")
}
cat("Comparison files saved under:", new_dir, "\n")
