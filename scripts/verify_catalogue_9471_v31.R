#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("Missing ", flag, call. = FALSE)
  args[[i + 1L]]
}
root <- normalizePath(arg("--project"), winslash = "/", mustWork = TRUE)
base <- file.path(root, "results/module12_rebuild_staging/colocalization")
output_arg <- match("--output", args)
independent <- if (is.na(output_arg)) {
  file.path(base, "catalogue_9471_independent")
} else {
  if (output_arg == length(args)) stop("Missing --output value.", call. = FALSE)
  normalizePath(args[[output_arg + 1L]], winslash = "/", mustWork = TRUE)
}
fresh <- file.path(independent, "module12_3_coloc_catalogue_summary.tsv.gz")
reference <- file.path(base, "unified_catalogue/summary/module12_3_coloc_catalogue_summary.tsv.gz")
manifest <- file.path(independent, "source_manifest.tsv.gz")
if (!all(file.exists(c(fresh, reference, manifest))))
  stop("Independent summary, unified reference, or manifest missing.", call. = FALSE)
a <- fread(fresh, showProgress = FALSE)
b <- fread(reference, showProgress = FALSE)
m <- fread(manifest, showProgress = FALSE)
if (!all(c("coloc_unit_id") %in% names(a)) || !"coloc_unit_id" %in% names(b))
  stop("Missing ID column.", call. = FALSE)
if (nrow(a) != 9471L || nrow(b) != 9471L || nrow(m) != 9471L ||
    uniqueN(a$coloc_unit_id) != 9471L || uniqueN(b$coloc_unit_id) != 9471L ||
    uniqueN(m$coloc_unit_id) != 9471L ||
    !setequal(a$coloc_unit_id, b$coloc_unit_id) ||
    !setequal(a$coloc_unit_id, m$coloc_unit_id))
  stop("9471-unit ID membership failed.", call. = FALSE)
if (any(m$source != "staged_unit")) stop("Non-staging source in manifest.", call. = FALSE)
setkey(a, coloc_unit_id)
setkey(b, coloc_unit_id)
b <- b[a$coloc_unit_id]
common <- intersect(setdiff(names(a), "coloc_unit_id"), names(b))
only_a <- setdiff(names(a), names(b))
only_b <- setdiff(names(b), names(a))
differences <- data.table(column = character(), mismatches = integer(),
                          max_abs_difference = numeric())
for (column in common) {
  x <- a[[column]]
  y <- b[[column]]
  if (is.numeric(x) && is.numeric(y)) {
    unequal_na <- xor(is.na(x), is.na(y))
    valid <- !is.na(x) & !is.na(y)
    d <- abs(x[valid] - y[valid])
    bad <- sum(unequal_na) + sum(!is.finite(d) | d > 1e-12)
    maxd <- if (length(d) && all(is.finite(d))) max(d) else NA_real_
  } else {
    bad <- sum(is.na(x) != is.na(y) | (!is.na(x) & !is.na(y) & as.character(x) != as.character(y)))
    maxd <- NA_real_
  }
  if (bad) differences <- rbind(differences, data.table(
    column = column, mismatches = as.integer(bad), max_abs_difference = maxd))
}
posterior <- grep("PP[._]?H[0-4]|posterior", names(a), value = TRUE, ignore.case = TRUE)
result <- if (nrow(differences) == 0L && length(only_a) == 0L && length(only_b) == 0L)
  "PASS" else "FAIL"
cat("Result:", result, "\n")
cat("Independent rows:", nrow(a), " Unified rows:", nrow(b),
    " Staged manifest rows:", nrow(m), "\n")
cat("Compared columns:", length(common), " Posterior columns:",
    paste(posterior, collapse = ", "), "\n")
cat("Only independent columns:", paste(only_a, collapse = ", "), "\n")
cat("Only unified columns:", paste(only_b, collapse = ", "), "\n")
if (nrow(differences)) print(differences)
if (length(intersect(c("interpretation", "coloc_interpretation"), names(a)))) {
  col <- intersect(c("interpretation", "coloc_interpretation"), names(a))[[1L]]
  print(as.data.table(table(a[[col]], useNA = "ifany")))
}
quit(save = "no", status = if (result == "PASS") 0L else 1L)
