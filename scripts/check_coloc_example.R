#!/usr/bin/env Rscript
# Run from the first_release_candidate directory: Rscript scripts/check_coloc_example.R
if (!requireNamespace("coloc", quietly = TRUE)) {
  stop("Install coloc first: install.packages('coloc')", call. = FALSE)
}
input <- "examples/CLUMP_0001/coloc/COLOC_0000001_harmonized_input.tsv.gz"
expected_file <- "examples/CLUMP_0001/coloc/expected_summary.tsv"
stopifnot(file.exists(input), file.exists(expected_file))
x <- read.delim(gzfile(input), check.names = FALSE)
expected <- read.delim(expected_file, check.names = FALSE)
stopifnot(nrow(x) == 1524L, nrow(expected) == 1L,
          !anyNA(x[, c("gwas_beta", "gwas_se", "apaqtl_beta_harmonized",
                        "apaqtl_se", "coloc_snp", "position")]))
gwas_file <- "examples/CLUMP_0001/gwas/gwas_region.tsv.gz"
apaqtl_file <- paste0("examples/CLUMP_0001/apaqtl/",
                     "Skin_Not_Sun_Exposed_Suprapubic/",
                     "CLUMP_0001_chr1_1122649_1461195.apaqtl.tsv.gz")
gwas_raw <- read.delim(gzfile(gwas_file), check.names = FALSE)
apaqtl_raw <- read.delim(gzfile(apaqtl_file), check.names = FALSE)
apaqtl_raw <- apaqtl_raw[apaqtl_raw$phenotype_id == expected$phenotype_id, ]
# Join by exact variant identifiers: multiple different alleles can occupy
# the same genomic coordinate in the raw files.
gi <- match(x$gwas_variant_id, gwas_raw$variant_id)
qi <- match(x$apaqtl_variant_id, apaqtl_raw$variant_id)
stopifnot(!anyNA(gi), !anyNA(qi), !anyDuplicated(gi), !anyDuplicated(qi))
g <- gwas_raw[gi, ]
q <- apaqtl_raw[qi, ]
near <- function(a, b) all(is.finite(a) & is.finite(b) & abs(a - b) < 1e-9)
stopifnot(near(x$position, g$base_pair_location),
          near(x$position, q$variant_position),
          identical(as.character(x$gwas_effect_allele), as.character(g$effect_allele)),
          identical(as.character(x$gwas_other_allele), as.character(g$other_allele)),
          identical(as.character(x$apaqtl_ref_allele), as.character(q$ref)),
          identical(as.character(x$apaqtl_alt_allele), as.character(q$alt)),
          near(x$gwas_beta, g$beta), near(x$gwas_se, g$standard_error),
          near(x$apaqtl_beta, q$slope), near(x$apaqtl_se, q$slope_se),
          all(x$harmonization_status == "KEEP_DIRECT"),
          near(x$apaqtl_beta_harmonized, q$slope))
cat("PASS: all 1524 harmonized variants match raw GWAS and apaQTL records.\n")
official_cases <- 156778
official_controls <- 262691
stopifnot(official_cases + official_controls == expected$gwas_sample_size)
official_fraction <- official_cases / (official_cases + official_controls)
gwas <- list(beta = x$gwas_beta, varbeta = x$gwas_se^2,
             snp = x$coloc_snp, position = x$position,
             type = "cc", N = official_cases + official_controls,
             s = official_fraction,
             MAF = x$coloc_maf_gwas)
apaqtl <- list(beta = x$apaqtl_beta_harmonized,
               varbeta = x$apaqtl_se^2, snp = x$coloc_snp,
               position = x$position, type = "quant", N = 517,
               MAF = x$coloc_maf_apaqtl)
result <- coloc::coloc.abf(gwas, apaqtl, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
names <- paste0("PP.H", 0:4)
actual <- as.numeric(result$summary[paste0(names, ".abf")])
target <- as.numeric(expected[1, names])
if (anyNA(actual) || any(abs(actual - target) > 0.005)) {
  stop(sprintf("Posterior mismatch. Actual: %s; expected: %s",
               paste(signif(actual, 5), collapse = ", "),
               paste(signif(target, 5), collapse = ", ")), call. = FALSE)
}
cat("PASS: COLOC_0000001, 1524 variants; PP.H0–H4 within 0.005.\n")
cat("GWAS case fraction from Catalog counts: ",
    format(official_fraction, digits = 10), "\n", sep = "")
cat("coloc version: ", as.character(packageVersion("coloc")), "\n", sep = "")
