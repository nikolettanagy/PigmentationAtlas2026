#!/usr/bin/env Rscript
# Usage: Rscript scripts/audit_new_coloc_units.R PROJECT_ROOT
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply PigmentationAtlas project root")
suppressPackageStartupMessages(library(data.table))
root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
base <- file.path(root, "results/module12_rebuild_staging/colocalization")
qdir <- file.path(base, "work_queue")
queue <- fread(file.path(qdir, "module12_1_coloc_work_queue_eligible.tsv.gz"))
new <- fread(file.path(qdir, "newly_eligible_units.tsv"))
key <- c("clump_id", "tissue", "phenotype_id", "gene_id")
targets <- merge(queue, unique(new[, ..key]), by = key, sort = FALSE)
if (nrow(targets) != 1631L || anyDuplicated(targets$coloc_unit_id))
  stop("New eligible unit list is incomplete or duplicated")
status_path <- file.path(base, "new_units_run_status.tsv")
if (!file.exists(status_path)) stop("No run status available")
run <- fread(status_path)
if (anyDuplicated(run$coloc_unit_id) || any(!run$coloc_unit_id %in% targets$coloc_unit_id))
  stop("Run status has duplicate or unexpected IDs")
rows <- vector("list", nrow(targets))
for (i in seq_len(nrow(targets))) {
  t <- targets[i]
  id <- t$coloc_unit_id[[1L]]
  result_dir <- file.path(base, "single_test/results", paste0(id, "_", t$clump_id[[1L]]))
  sum_file <- file.path(result_dir, paste0(id, "_coloc_summary.tsv"))
  stat_file <- file.path(result_dir, paste0(id, "_status.tsv"))
  r <- run[coloc_unit_id == id]
  state <- if (nrow(r)) r$status[[1L]] else "MISSING"
  posterior <- NA_real_
  variants <- NA_integer_
  result_ok <- FALSE
  if (file.exists(sum_file) && file.exists(stat_file)) {
    s <- tryCatch(fread(sum_file), error = function(e) NULL)
    qc <- tryCatch(fread(stat_file), error = function(e) NULL)
    if (!is.null(s) && nrow(s) == 1L && !is.null(qc) && nrow(qc) == 1L &&
        s$coloc_unit_id[[1L]] == id && qc$final_status[[1L]] == state) {
      posterior <- as.numeric(s$PP.H4[[1L]])
      variants <- as.integer(s$nsnps[[1L]])
      result_ok <- is.finite(posterior) && posterior >= 0 && posterior <= 1 &&
        !is.na(variants) && variants >= 50L &&
        isTRUE(all.equal(as.numeric(r$PP.H4[[1L]]), posterior, tolerance = 1e-6))
    }
  }
  rows[[i]] <- data.table(coloc_unit_id = id, clump_id = t$clump_id[[1L]],
    status = state, result_verified = result_ok, PP.H4 = posterior, nsnps = variants)
}
audit <- rbindlist(rows)
output <- file.path(base, "new_units_audit.tsv")
fwrite(audit, output, sep = "\t")
cat("Target units:", nrow(audit), "\n")
print(audit[, .N, by = .(status, result_verified)][order(status, result_verified)])
cat("Strong H4 (>=0.8):", sum(audit$PP.H4 >= .8, na.rm = TRUE), "\n")
cat("Moderate H4 (>=0.5 and <0.8):",
    sum(audit$PP.H4 >= .5 & audit$PP.H4 < .8, na.rm = TRUE), "\n")
cat("Report:", output, "\n")
if (any(!audit$result_verified)) quit(status = 1L)
cat("PASS: all 1631 newly eligible units have verified coloc output.\n")
