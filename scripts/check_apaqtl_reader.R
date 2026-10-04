#!/usr/bin/env Rscript
# Run from the first_release_candidate root.
source("scripts/module11_3B_2_read_and_qc_apaqtl.R", local = TRUE)
stopifnot(exists("read_and_qc_apaqtl", mode = "function"))
file <- paste0("examples/CLUMP_0001/apaqtl/",
               "Skin_Not_Sun_Exposed_Suprapubic/",
               "CLUMP_0001_chr1_1122649_1461195.apaqtl.tsv.gz")
metadata <- list(queue_id = "APAQTL_TEST_CLUMP_0001",
                 clump_id = "CLUMP_0001",
                 tissue = "Skin_Not_Sun_Exposed_Suprapubic",
                 canonical_chromosome = "1",
                 canonical_region_start = 1122649,
                 canonical_region_end = 1461195)
result <- read_and_qc_apaqtl(file_path = file, metadata = metadata,
                             retain_source_columns = FALSE,
                             fail_on_missing_pvalue = TRUE,
                             fail_on_missing_variant = TRUE)
stopifnot(is.list(result),
          all(c("data", "qc", "column_map", "duplicate_rows") %in% names(result)),
          nrow(result$data) > 0L,
          nrow(result$qc) == 1L)
cat("PASS: sourceable Module 11.3B.2 reader processed", nrow(result$data),
    "apaQTL rows from CLUMP_0001.\n")
cat("Reader status:", as.character(result$qc$file_status[[1L]]), "\n")
