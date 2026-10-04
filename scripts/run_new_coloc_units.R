#!/usr/bin/env Rscript
# Run only newly eligible 12.1 units in the isolated 12.2 output directory.
# Usage: Rscript scripts/run_new_coloc_units.R PROJECT_ROOT [--max-units=10]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L)
  stop("Usage: Rscript scripts/run_new_coloc_units.R PROJECT_ROOT [--max-units=N]")
suppressPackageStartupMessages(library(data.table))
root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
limit <- Inf
if (length(args) == 2L) {
  if (!grepl("^--max-units=[1-9][0-9]*$", args[[2L]])) stop("Invalid --max-units")
  limit <- as.integer(sub("^--max-units=", "", args[[2L]]))
}
base <- file.path(root, "results/module12_rebuild_staging/colocalization")
qdir <- file.path(base, "work_queue")
queue <- fread(file.path(qdir, "module12_1_coloc_work_queue_eligible.tsv.gz"))
new <- fread(file.path(qdir, "newly_eligible_units.tsv"))
key <- c("clump_id", "tissue", "phenotype_id", "gene_id")
stopifnot(all(key %in% names(queue)), all(key %in% names(new)),
          !anyDuplicated(queue$clump_tissue_phenotype_gene_key))
targets <- merge(queue, unique(new[, ..key]), by = key, sort = FALSE)
setorder(targets, queue_order)
if (nrow(targets) != 1631L) stop("Expected 1631 new eligible units; found ", nrow(targets))
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows")
                     "Rscript.exe" else "Rscript")
engine <- normalizePath("scripts/staged_module12_2_smoke.R", winslash = "/", mustWork = TRUE)
if (!file.exists(rscript)) stop("Rscript executable missing: ", rscript)
status_path <- file.path(base, "new_units_run_status.tsv")
logs <- file.path(base, "new_units_logs")
dir.create(logs, recursive = TRUE, showWarnings = FALSE)
status <- if (file.exists(status_path)) fread(status_path) else data.table(
  coloc_unit_id = character(), clump_id = character(), tissue = character(),
  phenotype_id = character(), gene_id = character(), status = character(),
  exit_code = integer(), PP.H4 = numeric(), nsnps = integer(), log_file = character())
selected <- 0L
failures <- 0L
for (i in seq_len(nrow(targets))) {
  unit <- targets[i]
  id <- unit$coloc_unit_id[[1L]]
  result_dir <- file.path(base, "single_test/results", paste0(id, "_", unit$clump_id[[1L]]))
  summary_path <- file.path(result_dir, paste0(id, "_coloc_summary.tsv"))
  unit_status <- file.path(result_dir, paste0(id, "_status.tsv"))
  prior <- status[coloc_unit_id == id & status %in% c("PASS", "PASS_WITH_WARNINGS")]
  if (nrow(prior) && file.exists(summary_path) && file.exists(unit_status)) next
  if (selected >= limit) break
  selected <- selected + 1L
  log <- file.path(logs, paste0(id, ".log"))
  cat(sprintf("[%d/%d] %s %s\n", selected, nrow(targets), id, unit$clump_id[[1L]]))
  exit_code <- tryCatch(system2(rscript,
    args = c(shQuote(engine), shQuote(root), shQuote(id)),
    stdout = log, stderr = log, wait = TRUE), error = function(e) 1L)
  final <- "FAIL"
  h4 <- NA_real_
  nsnps <- NA_integer_
  if (exit_code == 0L && file.exists(summary_path) && file.exists(unit_status)) {
    s <- tryCatch(fread(summary_path), error = function(e) NULL)
    qc <- tryCatch(fread(unit_status), error = function(e) NULL)
    if (!is.null(s) && nrow(s) == 1L && !is.null(qc) && nrow(qc) == 1L &&
        s$coloc_unit_id[[1L]] == id &&
        qc$final_status[[1L]] %in% c("PASS", "PASS_WITH_WARNINGS")) {
      final <- qc$final_status[[1L]]
      h4 <- as.numeric(s$PP.H4[[1L]])
      nsnps <- as.integer(s$nsnps[[1L]])
    }
  }
  status <- status[coloc_unit_id != id]
  status <- rbind(status, data.table(coloc_unit_id = id,
    clump_id = unit$clump_id[[1L]], tissue = unit$tissue[[1L]],
    phenotype_id = unit$phenotype_id[[1L]], gene_id = unit$gene_id[[1L]],
    status = final, exit_code = as.integer(exit_code), PP.H4 = h4,
    nsnps = nsnps, log_file = log), fill = TRUE)
  fwrite(status, status_path, sep = "\t")
  cat("  ", final, " PP.H4=", h4, " nsnps=", nsnps, "\n", sep = "")
  if (final == "FAIL") failures <- failures + 1L
}
cat("New units targeted:", nrow(targets), "processed now:", selected,
    "failures now:", failures, "recorded successes:",
    sum(status$status %in% c("PASS", "PASS_WITH_WARNINGS")), "\n")
cat("Run status:", status_path, "\n")
if (failures) quit(status = 1L)
