#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value for ", flag, call. = FALSE)
  args[[i + 1L]]
}
root <- normalizePath(opt("--workspace"), winslash = "/", mustWork = TRUE)
start <- opt("--start-at", "S1")
if (!start %in% c("S1", "S2")) stop("--start-at must be S1 or S2.", call. = FALSE)
rscript <- opt("--rscript", file.path(R.home("bin"), "Rscript.exe"))
if (!file.exists(rscript)) stop("Rscript executable missing: ", rscript, call. = FALSE)
setwd(root)
Sys.setenv(PIGMENTATIONATLAS_ROOT = root)
log_dir <- file.path(root, "logs", "phase1_v36")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
status_file <- file.path(log_dir, "status.tsv")
if (file.exists(status_file)) stop("Existing v36 status; use a new audit directory.", call. = FALSE)
record <- function(stage, result, code, output) {
  line <- data.frame(stage = stage, result = result, exit_code = code,
                     output = output, time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))
  write.table(line, status_file, sep = "\t", row.names = FALSE,
              col.names = !file.exists(status_file), append = file.exists(status_file),
              quote = FALSE)
  cat(result, stage, output, "\n")
}
tail_log <- function(path) {
  if (file.exists(path)) {
    lines <- readLines(path, warn = FALSE)
    if (length(lines)) cat(paste(tail(lines, 30L), collapse = "\n"), "\n")
  }
}
run <- function(stage, command, arguments, expected) {
  if (file.exists(expected)) stop("Refusing existing output: ", expected, call. = FALSE)
  stdout <- file.path(log_dir, paste0(stage, ".stdout.log"))
  stderr <- file.path(log_dir, paste0(stage, ".stderr.log"))
  code <- tryCatch(
    suppressWarnings(system2(command, args = arguments, stdout = stdout,
                             stderr = stderr, wait = TRUE)),
    error = function(e) { cat(conditionMessage(e), "\n", file = stderr); 127L }
  )
  info <- file.info(expected)
  if (!identical(as.integer(code), 0L) || is.na(info$size) || info$size < 1) {
    record(stage, "FAIL", as.integer(code), expected)
    tail_log(stdout)
    tail_log(stderr)
    stop(stage, " failed; see ", log_dir, call. = FALSE)
  }
  record(stage, "PASS", 0L, expected)
}
run_r <- function(stage, script, expected) {
  if (!file.exists(script)) stop("Missing script: ", script, call. = FALSE)
  run(stage, rscript, shQuote(script), expected)
}
gwas <- "results/GCST90691600_harmonized.tsv.gz"
qc_file <- "results/GCST90691600_QC_summary.tsv"
if (start == "S1") {
  run_r("S1_GWAS", "scripts/upstream/Code S1.R", gwas)
}
if (!file.exists(gwas) || !file.exists(qc_file)) stop("S1 output/QC missing.", call. = FALSE)
qc <- read.delim(qc_file, check.names = FALSE)
if (nrow(qc) != 1L || qc$n_input_rows != 21874448 ||
    qc$n_clean_rows != 21831407 || qc$n_genomewide_significant != 33455) {
  stop("S1 QC does not match expected values.", call. = FALSE)
}
record("S1_QC", "PASS", 0L, gwas)
run_r("S2_clumping_input", "scripts/upstream/Code S2.R",
      "results/GCST90691600_clumping_input_uniqueID.tsv")
prefix <- "reference/1000G_GRCh38/1000G_GRCh38_EUR_uniqueID"
clump <- "results/GCST90691600_clump_uniqueID_r2_0.1_kb1000"
plink <- "tools/plink2/plink2.exe"
if (!file.exists(plink)) stop("Missing PLINK executable.", call. = FALSE)
run("PLINK_clumping", plink,
    c("--pfile", prefix, "vzs",
      "--clump", "results/GCST90691600_clumping_input_uniqueID.tsv",
      "--clump-id-field", "ID", "--clump-p-field", "P",
      "--clump-p1", "5e-8", "--clump-p2", "0.05",
      "--clump-r2", "0.1", "--clump-kb", "1000", "--out", clump),
    paste0(clump, ".clumps"))
run_r("intervals", "scripts/upstream/03_build_genomewide_loci.R",
      "results/GCST90691600_clump_intervals.tsv")
run_r("S3_loci", "scripts/upstream/Code S3.R",
      "results/GCST90691600_genomic_loci_v3.tsv")
run_r("S4_GENCODE", "scripts/upstream/Code S4.R",
      "data/annotations/processed/GENCODE_v50_genes_GRCh38.tsv.gz")
run_r("S5_overlaps", "scripts/upstream/Code S5.R",
      "results/GCST90691600_locus_gene_overlaps_v3.tsv")
run_r("S6_candidates", "scripts/upstream/Code S6.R",
      "results/GCST90691600_candidate_gene_prioritization_v2.tsv")
cat("PASS phase 1 complete. Status:", status_file, "\n")
