#!/usr/bin/env Rscript
# Finish the isolated 12.1 work queue using the 516 v11 checkpoints.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript scripts/finish_staged_module12_1.R PROJECT_ROOT")
suppressPackageStartupMessages(library(data.table))
root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
queue_dir <- file.path(root, "results/module12_rebuild_staging/colocalization/work_queue")
mapping_path <- file.path(queue_dir, "module12_1_file_mapping.tsv")
cache_dir <- file.path(queue_dir, "chunk_checkpoints_v11")
old_script <- "scripts/staged_module12_1.R"
if (!file.exists(mapping_path) || !file.exists(old_script))
  stop("Staged file mapping or v11 script missing. Run from the inner candidate folder.")
mapping <- fread(mapping_path)
if (nrow(mapping) != 516L) stop("Expected 516 mapped source files; found ", nrow(mapping))
script_md5 <- unname(tools::md5sum(old_script))
accepted_script_md5 <- c(script_md5, "a1725f3adc9f511c42e24d568a4b38e4")
parts <- vector("list", 516L)
for (i in seq_len(516L)) {
  path <- file.path(cache_dir, sprintf("part_%04d.rds", i))
  if (!file.exists(path)) stop("Checkpoint missing: ", path)
  obj <- readRDS(path)
  paths <- c(mapping$harmonized_chunk_file[[i]], mapping$gwas_file[[i]])
  info <- file.info(paths)
  sig <- list(paths = paths, bytes = info$size,
              modified = as.numeric(info$mtime), script = script_md5)
  if (!obj$signature$script %in% accepted_script_md5 ||
      !identical(obj$signature[names(sig) != "script"],
                 sig[names(sig) != "script"]))
    stop("Input changed since checkpoint ", i)
  if (!is.data.table(obj$data) || nrow(obj$data) < 1L)
    stop("Empty or invalid checkpoint ", i)
  if (any(obj$data$source_queue_id != mapping$source_queue_id[[i]]))
    stop("Source ID mismatch in checkpoint ", i)
  parts[[i]] <- obj$data
  if (i %% 50L == 0L) cat("Loaded and validated", i, "/ 516 checkpoints\n")
}
cat("Combining checkpoint results...\n")
queue <- rbindlist(parts, use.names = TRUE, fill = TRUE)
rm(parts)
gc(verbose = FALSE)
queue[, clump_tissue_phenotype_gene_key := paste(clump_id, tissue,
                                                  phenotype_id, gene_id, sep = "|")]
if (anyDuplicated(queue$clump_tissue_phenotype_gene_key)) {
  duplicates <- queue[duplicated(clump_tissue_phenotype_gene_key) |
                        duplicated(clump_tissue_phenotype_gene_key, fromLast = TRUE)]
  fwrite(duplicates, file.path(queue_dir, "module12_1_duplicate_units_from_checkpoints.tsv"),
         sep = "\t")
  stop("Duplicate analysis units; inspect module12_1_duplicate_units_from_checkpoints.tsv")
}
setorder(queue, clump_id, tissue, gene_id, phenotype_id)
queue[, queue_order := seq_len(.N)]
queue[, coloc_unit_id := sprintf("COLOC_%07d", seq_len(.N))]
setcolorder(queue, c("queue_order", "coloc_unit_id", "source_queue_id",
                     "clump_id", "tissue", "phenotype_id", "gene_id",
                     "clump_tissue_phenotype_gene_key",
                     setdiff(names(queue), c("queue_order", "coloc_unit_id",
                       "source_queue_id", "clump_id", "tissue", "phenotype_id",
                       "gene_id", "clump_tissue_phenotype_gene_key"))))
if (anyNA(queue$queue_status) ||
    any(!queue$queue_status %in% c("ELIGIBLE", "EXCLUDED")))
  stop("Invalid queue status in checkpoints.")
out <- file.path(queue_dir, "module12_1_coloc_work_queue.tsv.gz")
tmp <- tempfile(pattern = "module12_1_", tmpdir = queue_dir, fileext = ".tsv.gz")
fwrite(queue, tmp, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
if (!file.rename(tmp, out)) stop("Could not move completed queue into place.")
eligible <- queue[queue_status == "ELIGIBLE"]
excluded <- queue[queue_status == "EXCLUDED"]
fwrite(eligible, file.path(queue_dir, "module12_1_coloc_work_queue_eligible.tsv.gz"),
       sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
fwrite(excluded, file.path(queue_dir, "module12_1_coloc_work_queue_excluded.tsv.gz"),
       sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
fwrite(data.table(module = "12.1", final_status = "PASS_WITH_WARNINGS",
                  work_queue_file = out,
                  eligible_queue_file = file.path(queue_dir,
                    "module12_1_coloc_work_queue_eligible.tsv.gz"),
                  excluded_queue_file = file.path(queue_dir,
                    "module12_1_coloc_work_queue_excluded.tsv.gz"),
                  total_units = nrow(queue), eligible_units = nrow(eligible),
                  excluded_units = nrow(excluded),
                  created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
       file.path(queue_dir, "module12_1_status.tsv"), sep = "\t")
cat("PASS: staged 12.1 queue saved:", out, "\n")
cat("Total:", nrow(queue), "Eligible:", sum(queue$queue_status == "ELIGIBLE"),
    "Excluded:", sum(queue$queue_status == "EXCLUDED"), "\n")
