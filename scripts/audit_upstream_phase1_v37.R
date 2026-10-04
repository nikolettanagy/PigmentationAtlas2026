#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript audit_upstream_phase1_v37.R WORKSPACE")
root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
qc <- fread(file.path(root, "results/GCST90691600_QC_summary.tsv"))
clumps <- fread(file.path(root, "results/GCST90691600_clump_uniqueID_r2_0.1_kb1000.clumps"))
intervals <- fread(file.path(root, "results/GCST90691600_clump_intervals.tsv"))
genes <- fread(file.path(root, "data/annotations/processed/GENCODE_v50_genes_GRCh38.tsv.gz"))
candidates <- fread(file.path(root, "results/GCST90691600_candidate_gene_prioritization_v2.tsv"))
checks <- data.table(
  metric = c("gwas_input_rows", "gwas_clean_rows", "significant_variants",
             "plink_clumps", "clump_intervals", "gencode_genes", "candidate_rows"),
  observed = c(qc$n_input_rows[[1L]], qc$n_clean_rows[[1L]],
               qc$n_genomewide_significant[[1L]], nrow(clumps),
               nrow(intervals), nrow(genes), nrow(candidates)),
  expected = c(21874448, 21831407, 33455, 1325, 1325, NA, NA)
)
checks[, status := ifelse(is.na(expected), "RECORDED",
                          ifelse(observed == expected, "PASS", "FAIL"))]
print(checks)
if (any(checks$status == "FAIL")) quit(save = "no", status = 1L)
cat("PASS phase 1 fixed-count audit; descriptive counts recorded.\n")
