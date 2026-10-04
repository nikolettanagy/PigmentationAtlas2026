#!/usr/bin/env Rscript
# Usage: Rscript scripts/stage_apaqtl_repair.R C:/path/to/PigmentationAtlas representatives
# After reviewing the representative report, use "all" for every flagged chunk.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !args[2] %in% c("representatives", "all")) {
  stop("Usage: Rscript scripts/stage_apaqtl_repair.R PROJECT_ROOT representatives|all")
}
project_root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
mode <- args[2]
source("scripts/module11_3B_2_read_and_qc_apaqtl.R", local = TRUE)
cat_dir <- file.path(project_root, "results/module11/harmonization/catalogue")
queue <- data.table::fread(file.path(cat_dir, "module11_3B_work_queue.tsv"))
qc <- data.table::fread(file.path(cat_dir, "qc/final_qc/module11_3B_4_file_level_qc.tsv"))
stopifnot(!anyDuplicated(queue$queue_id), !anyDuplicated(qc$queue_id))
bad <- qc[chromosome_mismatch_rows > 0L,
          .(queue_id, clump_id, tissue, rows, chromosome_mismatch_rows,
            chunk_file, old_file_qc = final_file_qc)]
bad <- merge(bad, queue[, .(queue_id, canonical_chromosome,
                             canonical_region_start, canonical_region_end,
                             output_file_resolved)], by = "queue_id", sort = FALSE)
if (nrow(bad) != sum(qc$chromosome_mismatch_rows > 0L))
  stop("An affected QC record has no matching work-queue record.")
data.table::setorder(bad, canonical_chromosome, queue_id)
if (mode == "representatives")
  bad <- bad[, .SD[1L], by = canonical_chromosome]
stage <- file.path(project_root, "results/module11/catalogue_repair_staging")
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
report <- list()
for (i in seq_len(nrow(bad))) {
  row <- bad[i]
  source_file <- as.character(row$output_file_resolved)
  old_file <- as.character(row$chunk_file)
  if (!file.exists(source_file) || !file.exists(old_file))
    stop("Missing source or old chunk for ", row$queue_id)
  new_file <- file.path(stage, basename(old_file))
  if (file.exists(new_file)) {
    cat("Already staged:", row$queue_id, "(preserving existing file)\n")
    next
  }
  metadata <- as.list(row[, .(queue_id, clump_id, tissue, canonical_chromosome,
                               canonical_region_start, canonical_region_end)])
  x <- read_and_qc_apaqtl(source_file, metadata,
                          retain_source_columns = FALSE,
                          fail_on_missing_pvalue = TRUE,
                          fail_on_missing_variant = TRUE)
  d <- x$data
  if (nrow(d) == 0L || anyNA(d$chromosome_matches_canonical) ||
      any(!d$chromosome_matches_canonical) ||
      anyNA(d$position_within_canonical_region) ||
      any(!d$position_within_canonical_region))
    stop("Rebuilt chromosome, region, or nonempty validation failed: ", row$queue_id)
  variant_chr <- sub("^chr([^_]+)_.*$", "\\1", as.character(d$variant_id))
  if (any(variant_chr != as.character(row$canonical_chromosome)))
    stop("Rebuilt variant IDs disagree with canonical chromosome: ", row$queue_id)
  old <- data.table::fread(old_file,
                           select = c("variant_id", "phenotype_id", "position"))
  keys <- function(z) sort(paste(z$variant_id, z$phenotype_id, z$position,
                                 sep = "\t"))
  old_keys <- keys(old)
  new_keys <- keys(d)
  shared_keys <- sum(new_keys %in% old_keys)
  tmp <- tempfile(pattern = "repair_", tmpdir = stage, fileext = ".tsv.gz")
  data.table::fwrite(d, tmp, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
  if (!file.rename(tmp, new_file)) {
    unlink(tmp)
    stop("Could not move staged chunk into place: ", row$queue_id)
  }
  report[[length(report) + 1L]] <- data.table::data.table(
    queue_id = row$queue_id, chromosome = row$canonical_chromosome,
    old_rows = nrow(old), new_rows = nrow(d),
    shared_variant_phenotype_position_keys = shared_keys,
    old_mismatches = row$chromosome_mismatch_rows,
    new_mismatches = sum(!d$chromosome_matches_canonical),
    source_file = source_file, old_chunk = old_file, staged_chunk = new_file,
    old_md5 = unname(tools::md5sum(old_file)),
    staged_md5 = unname(tools::md5sum(new_file)))
  cat("PASS:", row$queue_id, "chr", row$canonical_chromosome,
      "old rows", nrow(old), "new rows", nrow(d),
      "shared keys", shared_keys, "new mismatches 0\n")
}
if (length(report)) {
  report_file <- file.path(stage, paste0("repair_report_", mode, ".tsv"))
  data.table::fwrite(data.table::rbindlist(report), report_file, sep = "\t")
  cat("Report:", report_file, "\n")
}
cat("Original catalogue and original chunks were not modified.\n")
