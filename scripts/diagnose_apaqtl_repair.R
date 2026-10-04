#!/usr/bin/env Rscript
# Read-only diagnosis: Rscript scripts/diagnose_apaqtl_repair.R PROJECT_ROOT QUEUE_ID
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply project root and queue ID.")
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
source("scripts/module11_3B_2_read_and_qc_apaqtl.R", local = TRUE)
cat_dir <- file.path(root, "results/module11/harmonization/catalogue")
queue <- data.table::fread(file.path(cat_dir, "module11_3B_work_queue.tsv"))
qc <- data.table::fread(file.path(cat_dir, "qc/final_qc/module11_3B_4_file_level_qc.tsv"))
q <- queue[queue_id == args[2]]
oldqc <- qc[queue_id == args[2]]
stopifnot(nrow(q) == 1L, nrow(oldqc) == 1L)
metadata <- as.list(q[, .(queue_id, clump_id, tissue, canonical_chromosome,
                          canonical_region_start, canonical_region_end)])
x <- read_and_qc_apaqtl(q$output_file_resolved, metadata,
                        retain_source_columns = FALSE)
d <- x$data
old <- data.table::fread(oldqc$chunk_file,
                         select = c("variant_id", "phenotype_id", "position",
                                    "chromosome", "canonical_chromosome"))
cat("Queue:", args[2], "canonical chr", q$canonical_chromosome, "\n")
cat("Rows: saved", nrow(old), "expected QC", oldqc$rows,
    "recomputed", nrow(d), "\n")
cat("Chromosome mismatches: saved", oldqc$chromosome_mismatch_rows,
    "recomputed", sum(!d$chromosome_matches_canonical, na.rm = TRUE), "\n")
cat("Missing chromosome-match flag:", sum(is.na(d$chromosome_matches_canonical)), "\n")
cat("Outside region: saved", oldqc$outside_region_rows,
    "recomputed", sum(!d$position_within_canonical_region, na.rm = TRUE), "\n")
cat("Missing in-region flag:", sum(is.na(d$position_within_canonical_region)), "\n")
cat("Reader status:", as.character(x$qc$file_status[[1L]]), "\n")
cat("Saved first three:\n")
print(old[1:min(3L, nrow(old)), .(chromosome, canonical_chromosome,
                                   position, variant_id)])
cat("Recomputed first three:\n")
print(d[1:min(3L, nrow(d)), .(chromosome, canonical_chromosome,
                               position, variant_id)])
cat("No files were changed.\n")
