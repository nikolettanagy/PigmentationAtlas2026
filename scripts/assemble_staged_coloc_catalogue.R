#!/usr/bin/env Rscript
# Re-key the old 7840 unchanged units by biological key and append 1631 new runs.
# All outputs stay under results/module12_rebuild_staging.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply PigmentationAtlas project root")
suppressPackageStartupMessages(library(data.table))
root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
oldbase <- file.path(root, "results/module12/colocalization")
newbase <- file.path(root, "results/module12_rebuild_staging/colocalization")
oldq <- fread(file.path(oldbase, "work_queue/module12_1_coloc_work_queue_eligible.tsv.gz"))
newq <- fread(file.path(newbase, "work_queue/module12_1_coloc_work_queue_eligible.tsv.gz"))
key <- "clump_tissue_phenotype_gene_key"
if (nrow(oldq) != 7840L || nrow(newq) != 9471L ||
    anyDuplicated(oldq[[key]]) || anyDuplicated(newq[[key]]))
  stop("Old or staged eligible queue does not match verified counts")
old_summary <- fread(file.path(oldbase,
  "catalogue/summary/module12_3_coloc_catalogue_summary.tsv.gz"))
old_qc <- fread(file.path(oldbase,
  "catalogue/summary/module12_3_coloc_catalogue_qc.tsv.gz"))
if (nrow(old_summary) != 7840L || anyDuplicated(old_summary$coloc_unit_id) ||
    nrow(old_qc) != 7840L || anyDuplicated(old_qc$coloc_unit_id))
  stop("Original summary or QC does not contain 7840 unique units")
if (!setequal(oldq$coloc_unit_id, old_summary$coloc_unit_id) ||
    !setequal(oldq$coloc_unit_id, old_qc$coloc_unit_id))
  stop("Original summary/QC unit IDs differ from original eligible queue")
old_summary <- old_summary[match(oldq$coloc_unit_id, old_summary$coloc_unit_id)]
old_qc <- old_qc[match(oldq$coloc_unit_id, old_qc$coloc_unit_id)]
biological_cols <- c("clump_id", "tissue", "phenotype_id", "gene_id")
for (col in biological_cols) {
  if (col %in% names(old_summary) &&
      !identical(as.character(old_summary[[col]]), as.character(oldq[[col]])))
    stop("Original summary biological key mismatch: ", col)
}
idx <- match(oldq[[key]], newq[[key]])
if (anyNA(idx)) stop("An originally eligible biological unit disappeared")
for (col in c(biological_cols, "apaqtl_chunk_file", "gwas_file")) {
  if (!identical(as.character(oldq[[col]]), as.character(newq[[col]][idx])))
    stop("Previously eligible unit input changed: ", col)
}
old_ids <- oldq$coloc_unit_id
new_ids <- newq$coloc_unit_id[idx]
old_summary[, coloc_unit_id := new_ids]
old_qc[, coloc_unit_id := new_ids]
new_only <- newq[!get(key) %in% oldq[[key]]]
if (nrow(new_only) != 1631L) stop("Expected exactly 1631 new eligible units")
audit <- fread(file.path(newbase, "new_units_audit.tsv"))
if (nrow(audit) != 1631L || any(!audit$result_verified) ||
    !setequal(audit$coloc_unit_id, new_only$coloc_unit_id))
  stop("Run the full new-unit audit successfully before assembling")
sum_new <- vector("list", 1631L)
qc_new <- vector("list", 1631L)
for (i in seq_len(1631L)) {
  unit <- new_only[i]
  id <- unit$coloc_unit_id[[1L]]
  folder <- file.path(newbase, "single_test/results",
                      paste0(id, "_", unit$clump_id[[1L]]))
  sum_new[[i]] <- fread(file.path(folder, paste0(id, "_coloc_summary.tsv")))
  qc_new[[i]] <- fread(file.path(newbase, "single_test/qc", paste0(id, "_qc.tsv")))
  if (nrow(sum_new[[i]]) != 1L || sum_new[[i]]$coloc_unit_id[[1L]] != id ||
      nrow(qc_new[[i]]) != 1L || qc_new[[i]]$coloc_unit_id[[1L]] != id)
    stop("Unexpected new unit summary/QC: ", id)
}
all_summary <- rbindlist(c(list(old_summary), sum_new), fill = TRUE, use.names = TRUE)
all_qc <- rbindlist(c(list(old_qc), qc_new), fill = TRUE, use.names = TRUE)
if (nrow(all_summary) != 9471L || anyDuplicated(all_summary$coloc_unit_id) ||
    nrow(all_qc) != 9471L || anyDuplicated(all_qc$coloc_unit_id) ||
    !setequal(all_summary$coloc_unit_id, newq$coloc_unit_id) ||
    !setequal(all_qc$coloc_unit_id, newq$coloc_unit_id))
  stop("Combined catalogue does not cover the staged eligible queue")
setorder(all_summary, coloc_unit_id)
setorder(all_qc, coloc_unit_id)
out <- file.path(newbase, "catalogue/summary")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(all_summary, file.path(out, "module12_3_coloc_catalogue_summary.tsv.gz"),
       sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
fwrite(all_qc, file.path(out, "module12_3_coloc_catalogue_qc.tsv.gz"),
       sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
fwrite(data.table(old_coloc_unit_id = old_ids, staged_coloc_unit_id = new_ids,
                  biological_key = oldq[[key]]),
       file.path(out, "legacy_to_staged_unit_ids.tsv"), sep = "\t")
cat("PASS: 7840 unchanged units re-keyed; 1631 new units appended.\n")
cat("Staged catalogue: 9471 unique eligible units.\n")
cat("Old GWAS case fraction was 0.3737566; new engine used 156778/419469.\n")
cat("Combined catalogue assembled at:", out, "\n")
