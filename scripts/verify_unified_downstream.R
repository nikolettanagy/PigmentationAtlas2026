#!/usr/bin/env Rscript
# Prove whether staged 12.4–12.6 outputs built from the provisional merged
# catalogue are identical in all scientific source fields to the unified run.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply PigmentationAtlas project root")
suppressPackageStartupMessages(library(data.table))
root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
base <- file.path(root, "results/module12_rebuild_staging/colocalization")
hybrid <- fread(file.path(base, "catalogue/summary/module12_3_coloc_catalogue_summary.tsv.gz"))
unified <- fread(file.path(base, "unified_catalogue/summary/module12_3_coloc_catalogue_summary.tsv.gz"))
fields <- c("clump_id", "tissue", "phenotype_id", "gene_id", "nsnps",
            "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4",
            "PP.H4_over_H3_H4", "posterior_sum", "top_coloc_variant",
            "top_coloc_variant_PP.H4")
posterior_fields <- c("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4",
                      "PP.H4_over_H3_H4", "posterior_sum",
                      "top_coloc_variant_PP.H4")
numeric_tolerance <- 1e-12
if (nrow(hybrid) != 9471L || nrow(unified) != 9471L ||
    anyDuplicated(hybrid$coloc_unit_id) || anyDuplicated(unified$coloc_unit_id) ||
    !setequal(hybrid$coloc_unit_id, unified$coloc_unit_id) ||
    !all(fields %in% names(hybrid)) || !all(fields %in% names(unified)))
  stop("The two catalogues cannot be compared by unit ID")
setkey(hybrid, coloc_unit_id)
setkey(unified, coloc_unit_id)
hybrid <- hybrid[unified$coloc_unit_id]
results <- rbindlist(lapply(fields, function(f) {
  a <- hybrid[[f]]; b <- unified[[f]]
  if (is.numeric(a) && is.numeric(b)) {
    comparable <- !is.na(a) & !is.na(b)
    max_delta <- if (!any(comparable)) 0 else
      max(abs(a[comparable] - b[comparable]))
    exact_mismatch <- sum((is.na(a) != is.na(b)) |
                          (!is.na(a) & !is.na(b) & a != b))
    mismatch <- sum((is.na(a) != is.na(b)) |
                    (!is.na(a) & !is.na(b) &
                     abs(a - b) > if (f %in% posterior_fields) numeric_tolerance else 0))
  } else {
    max_delta <- NA_real_
    mismatch <- sum((is.na(a) != is.na(b)) |
                    (!is.na(a) & !is.na(b) & as.character(a) != as.character(b)))
    exact_mismatch <- mismatch
  }
  data.table(field = f, exact_numeric_differences = exact_mismatch,
             material_mismatches = mismatch, maximum_absolute_difference = max_delta)
}))
report <- file.path(base, "unified_catalogue/summary/unified_vs_staged_source_fields.tsv")
fwrite(results, report, sep = "\t")
print(results)
cat("Report:", report, "\n")
if (any(results$material_mismatches != 0L))
  stop("Scientific fields differ beyond 1e-12 posterior tolerance; rerun downstream modules")
cat("PASS: all 9471 downstream scientific fields agree; posterior tolerance 1e-12.\n")
