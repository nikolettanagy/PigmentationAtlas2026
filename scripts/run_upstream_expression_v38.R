#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("Missing ", flag, call. = FALSE)
  args[[i + 1L]]
}
root <- normalizePath(opt("--workspace"), winslash = "/", mustWork = TRUE)
source_root <- normalizePath(opt("--input-root"), winslash = "/", mustWork = TRUE)
rscript <- opt("--rscript")
if (!file.exists(rscript)) stop("Rscript missing: ", rscript, call. = FALSE)
if (!requireNamespace("arrow", quietly = TRUE)) stop("R package arrow is required.", call. = FALSE)
if (root == source_root || startsWith(paste0(root, "/"), paste0(source_root, "/")))
  stop("Workspace must be separate from the original project.", call. = FALSE)
setwd(root)
Sys.setenv(PIGMENTATIONATLAS_ROOT = root)
log_dir <- file.path(root, "logs", "expression_v38")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
status_file <- file.path(log_dir, "status.tsv")
if (file.exists(status_file)) stop("Existing phase-2 status; refusing overwrite.", call. = FALSE)
manifest <- fread("config/expression_datasets.tsv")
enabled <- tolower(as.character(manifest$enabled)) %in% c("true", "1", "yes")
if (sum(enabled) != 2L) stop("Expected exactly two enabled expression datasets.", call. = FALSE)
for (rel in manifest$file_path[enabled]) {
  if (grepl("^([A-Za-z]:|/)|\\.\\.", rel)) stop("Unsafe input path: ", rel, call. = FALSE)
  src <- file.path(source_root, rel)
  dst <- file.path(root, rel)
  if (!file.exists(src)) stop("Missing source expression input: ", src, call. = FALSE)
  if (file.exists(dst)) stop("Refusing existing workspace input: ", dst, call. = FALSE)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(src, dst) || file.info(src)$size != file.info(dst)$size)
    stop("Copy failed or size mismatch: ", rel, call. = FALSE)
  cat("COPIED", rel, "\n")
}
run <- function(stage, script, expected) {
  if (file.exists(expected)) stop("Refusing existing result: ", expected, call. = FALSE)
  out <- file.path(log_dir, paste0(stage, ".stdout.log"))
  err <- file.path(log_dir, paste0(stage, ".stderr.log"))
  code <- tryCatch(suppressWarnings(system2(rscript, args = shQuote(script),
                                           stdout = out, stderr = err, wait = TRUE)),
                   error = function(e) { cat(conditionMessage(e), file = err); 127L })
  ok <- as.integer(code) == 0L && file.exists(expected) &&
    !is.na(file.info(expected)$size) && file.info(expected)$size > 0
  status <- data.table(stage = stage, result = if (ok) "PASS" else "FAIL",
                       exit_code = as.integer(code), output = expected)
  fwrite(status, status_file, sep = "\t", append = file.exists(status_file),
         col.names = !file.exists(status_file))
  if (!ok) {
    for (p in c(out, err)) if (file.exists(p)) cat(paste(tail(readLines(p, warn = FALSE), 30L), collapse = "\n"), "\n")
    stop(stage, " failed; see ", log_dir, call. = FALSE)
  }
  cat("PASS", stage, expected, "\n")
}
run("S7_expression_database", "scripts/upstream/Code S7.R",
    "results/expression/database/expression_database_index.tsv")
index <- fread("results/expression/database/expression_database_index.tsv")
if (nrow(index) != 2L || uniqueN(index$dataset_id) != 2L)
  stop("S7 did not build both expression datasets.", call. = FALSE)
run("S8_expression_evidence", "scripts/upstream/Code S8.R",
    "results/expression/expression_summary.tsv")
summary <- fread("results/expression/expression_summary.tsv")
if (!nrow(summary)) stop("S8 summary is empty.", call. = FALSE)
cat("PASS expression phase; datasets=2, summary rows=", nrow(summary), "\n", sep = "")
