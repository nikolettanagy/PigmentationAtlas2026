#!/usr/bin/env Rscript
# Usage: Rscript scripts/capture_release_environment.R PROJECT_ROOT
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply PigmentationAtlas project root")
suppressPackageStartupMessages(library(data.table))
root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
out <- file.path(root, "results/module12_rebuild_staging/release_audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
packages <- c("data.table", "coloc", "susieR", "ggplot2", "patchwork",
              "ggalluvial", "ggrepel", "scales", "openxlsx", "arrow",
              "rtracklayer", "renv")
versions <- vapply(packages, function(x) {
  if (requireNamespace(x, quietly = TRUE)) as.character(packageVersion(x))
  else NA_character_
}, character(1L))
fwrite(data.table(package = packages, installed_version = unname(versions)),
       file.path(out, "R_package_versions.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
fwrite(data.table(component = c("R", "platform", "OS", "GWAS_case_fraction",
                                "expected_unified_coloc_units"),
                  value = c(as.character(getRversion()), R.version$platform,
                            paste(Sys.info()[["sysname"]], Sys.info()[["release"]]),
                            sprintf("%.16g", 156778 / 419469), "9471")),
       file.path(out, "run_parameters.tsv"), sep = "\t")
paths <- c(
  "data/gwas/GCST90691600.h.tsv",
  "data/annotations/processed/GENCODE_v50_genes_GRCh38.tsv.gz",
  "reference/1000G_GRCh38/1000G_GRCh38_EUR_uniqueID.pgen",
  "reference/1000G_GRCh38/1000G_GRCh38_EUR_uniqueID.pvar.zst",
  "results/module11/harmonization/catalogue/module11_3B_work_queue.tsv",
  "results/module12_rebuild_staging/colocalization/work_queue/module12_1_coloc_work_queue_eligible.tsv.gz",
  "results/module12_rebuild_staging/colocalization/unified_catalogue/summary/module12_3_coloc_catalogue_summary.tsv.gz",
  "results/module12_rebuild_staging/reclassification/module12_6_clump_gene_reclassification.tsv.gz"
)
full <- file.path(root, paths)
info <- file.info(full)
small <- file.exists(full) & !is.na(info$size) & info$size <= 100 * 1024^2
hash <- rep(NA_character_, length(paths))
hash[small] <- unname(tools::md5sum(full[small]))
manifest <- data.table(relative_path = paths, exists = file.exists(full),
                       bytes = info$size, modified_local = as.character(info$mtime),
                       md5_if_at_most_100MiB = hash)
fwrite(manifest, file.path(out, "input_output_manifest.tsv"), sep = "\t", na = "NA")
cat("PASS: release environment and input/output metadata saved under:", out, "\n")
cat("R:", as.character(getRversion()), "coloc:", versions[["coloc"]], "\n")
cat("Manifest files present:", sum(manifest$exists), "/", nrow(manifest), "\n")
