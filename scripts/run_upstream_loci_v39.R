#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value for ", flag, call. = FALSE)
  args[[i + 1L]]
}
root <- normalizePath(opt("--workspace"), winslash = "/", mustWork = TRUE)
rscript <- opt("--rscript", file.path(R.home("bin"), "Rscript.exe"))
if (!file.exists(rscript)) stop("Rscript missing: ", rscript, call. = FALSE)
start <- opt("--start-at", "S9")
stages <- paste0("S", 9:13)
if (!start %in% stages) stop("--start-at must be S9, S10, S11, S12, or S13", call. = FALSE)
setwd(root)
Sys.setenv(PIGMENTATIONATLAS_ROOT = root)
log_dir <- file.path(root, "logs", "loci_v39")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
status_file <- file.path(log_dir, "status.tsv")
if (!file.exists("results/expression/expression_summary.tsv")) stop("S8 output missing.", call. = FALSE)
if (!file.exists("results/GCST90691600_genomic_loci_v3.tsv")) stop("S3 output missing.", call. = FALSE)
if (!file.exists("results/GCST90691600_harmonized.tsv.gz")) stop("S1 output missing.", call. = FALSE)
if (!file.exists("results/GCST90691600_clump_uniqueID_r2_0.1_kb1000.clumps")) stop("PLINK clumps missing.", call. = FALSE)
if (nrow(fread("results/expression/expression_summary.tsv")) != 1983L) stop("S8 row count mismatch.", call. = FALSE)
clumps <- fread("results/GCST90691600_clump_uniqueID_r2_0.1_kb1000.clumps")
if (nrow(clumps) != 1325L) stop("PLINK clump count mismatch.", call. = FALSE)
loci <- fread("results/GCST90691600_genomic_loci_v3.tsv")
if (!nrow(loci) || anyDuplicated(loci$locus_id)) stop("Invalid S3 loci.", call. = FALSE)
checks <- list(
  S9 = function() { p <- "results/loci/GCST90691600_locus_extraction_summary.tsv"; file.exists(p) && nrow(fread(p)) == nrow(loci) },
  S10 = function() { p <- "results/clump_master_metadata.tsv"; if (!file.exists(p)) return(FALSE); x <- fread(p); nrow(x) == 1325L && !anyDuplicated(x$clump_id) },
  S11 = function() { p <- "results/clump_workspace_manifest.tsv"; if (!file.exists(p)) return(FALSE); x <- fread(p); nrow(x) == 1325L && all(x$workspace_status == "ready") },
  S12 = function() { p <- "results/clump_gwas_extraction_manifest.tsv"; if (!file.exists(p)) return(FALSE); x <- fread(p); nrow(x) == 1325L && all(x$extraction_status == "ready") },
  S13 = function() { p <- "results/clump_gwas_validation_manifest.tsv"; if (!file.exists(p)) return(FALSE); x <- fread(p); nrow(x) == 1325L && all(x$validation_status == "ready") }
)
outputs <- c(S9="results/loci/GCST90691600_locus_extraction_summary.tsv", S10="results/clump_master_metadata.tsv", S11="results/clump_workspace_manifest.tsv", S12="results/clump_gwas_extraction_manifest.tsv", S13="results/clump_gwas_validation_manifest.tsv")
start_index <- match(start, stages)
if (start_index > 1L) for (stage in stages[seq_len(start_index - 1L)]) {
  if (!isTRUE(checks[[stage]]())) stop("Cannot resume: previous stage invalid: ", stage, call. = FALSE)
  cat("VERIFIED", stage, outputs[[stage]], "\n")
}
for (stage in stages[start_index:length(stages)]) {
  if (file.exists(outputs[[stage]])) stop("Existing output for ", stage, "; review before rerunning: ", outputs[[stage]], call. = FALSE)
  script <- file.path("scripts", "upstream", paste0("Code ", stage, ".R"))
  if (!file.exists(script)) stop("Missing script: ", script, call. = FALSE)
  stdout <- file.path(log_dir, paste0(stage, ".stdout.log"))
  stderr <- file.path(log_dir, paste0(stage, ".stderr.log"))
  code <- tryCatch(suppressWarnings(system2(rscript, shQuote(script), stdout=stdout, stderr=stderr, wait=TRUE)), error=function(e) {cat(conditionMessage(e), file=stderr); 127L})
  passed <- as.integer(code) == 0L && isTRUE(tryCatch(checks[[stage]](), error=function(e) FALSE))
  row <- data.table(stage=stage, result=if (passed) "PASS" else "FAIL", exit_code=as.integer(code), output=outputs[[stage]])
  fwrite(row, status_file, sep="\t", append=file.exists(status_file), col.names=!file.exists(status_file))
  if (!passed) {
    for (p in c(stdout, stderr)) if (file.exists(p)) cat(paste(tail(readLines(p, warn=FALSE), 35L), collapse="\n"), "\n")
    stop(stage, " failed; see ", log_dir, call. = FALSE)
  }
  cat("PASS", stage, outputs[[stage]], "\n")
}
cat("PASS locus and clump preparation; loci=", nrow(loci), ", clumps=1325\n", sep="")
