#!/usr/bin/env Rscript
# Usage: Rscript scripts/preflight_full_project.R C:/path/to/PigmentationAtlas
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the full PigmentationAtlas project directory.")
root <- normalizePath(args[1], mustWork = TRUE)
check <- function(kind, rel) {
  path <- file.path(root, rel)
  present <- if (kind == "dir") dir.exists(path) else file.exists(path)
  size <- if (present && kind == "file") file.info(path)$size else NA_real_
  cat(sprintf("%-7s %-13s %-80s %s\n",
              if (present) "PRESENT" else "MISSING", kind, rel,
              if (is.na(size)) "" else paste0(round(size / 1048576, 1), " MiB")))
  present
}
cat("PigmentationAtlas project root:", root, "\n")
files <- c(
  "data/gwas/GCST90691600.h.tsv",
  "results/GCST90691600_harmonized.tsv.gz",
  "results/GCST90691600_clumping_input_uniqueID.tsv",
  "results/GCST90691600_clump_uniqueID_r2_0.1_kb1000.clumps",
  "results/GCST90691600_clump_intervals.tsv",
  "data/annotations/processed/GENCODE_v50_genes_GRCh38.tsv.gz",
  "reference/1000G_GRCh38/1000G_GRCh38_EUR_uniqueID.pgen",
  "reference/1000G_GRCh38/1000G_GRCh38_EUR_uniqueID.pvar.zst",
  "reference/1000G_GRCh38/1000G_GRCh38_EUR_uniqueID.psam",
  "results/GCST90691600_candidate_gene_prioritization_v2.tsv",
  "config/expression_datasets.tsv",
  "scripts/Code S23.R", "scripts/Code S24.R",
  "scripts/nem kellőek/03_build_genomewide_loci.R",
  "scripts/module11_3B_2_read_and_qc_apaqtl.R"
)
present <- vapply(files, function(p) check("file", p), logical(1))
cat("\nFound", sum(present), "of", length(files), "key project files.\n")
cat("A missing final helper is expected in the current script collection; S24\n",
    "sources it by name. This preflight checks presence, not file validity.\n", sep = "")
