library(data.table)

input_file <- "data/gwas/GCST90691600.h.tsv"
output_file <- "results/GCST90691600_harmonized.tsv.gz"
qc_file <- "results/GCST90691600_QC_summary.tsv"

cat("Reading GWAS file:\n", input_file, "\n")

gwas <- fread(
  input_file,
  sep = "\t",
  showProgress = TRUE
)

cat("Rows read:", nrow(gwas), "\n")

required_columns <- c(
  "chromosome",
  "base_pair_location",
  "effect_allele",
  "other_allele",
  "beta",
  "standard_error",
  "effect_allele_frequency",
  "p_value",
  "rsid",
  "variant_id"
)

missing_columns <- setdiff(required_columns, names(gwas))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

gwas_clean <- gwas[
  !is.na(chromosome) &
  !is.na(base_pair_location) &
  !is.na(effect_allele) &
  !is.na(other_allele) &
  !is.na(beta) &
  !is.na(standard_error) &
  !is.na(effect_allele_frequency) &
  !is.na(p_value) &
  !is.na(rsid) &
  standard_error > 0 &
  p_value > 0 &
  p_value <= 1 &
  effect_allele_frequency > 0 &
  effect_allele_frequency < 1
]

gwas_clean[, chromosome := as.character(chromosome)]
gwas_clean[, effect_allele := toupper(effect_allele)]
gwas_clean[, other_allele := toupper(other_allele)]

gwas_clean[, maf := pmin(
  effect_allele_frequency,
  1 - effect_allele_frequency
)]

gwas_clean[, z_score := beta / standard_error]

gwas_clean <- unique(
  gwas_clean,
  by = "rsid"
)

setorder(
  gwas_clean,
  chromosome,
  base_pair_location
)

qc_summary <- data.table(
  input_file = input_file,
  n_input_rows = nrow(gwas),
  n_clean_rows = nrow(gwas_clean),
  n_removed_rows = nrow(gwas) - nrow(gwas_clean),
  n_unique_rsid = uniqueN(gwas_clean$rsid),
  n_genomewide_significant = sum(
    gwas_clean$p_value < 5e-8,
    na.rm = TRUE
  ),
  minimum_p_value = min(
    gwas_clean$p_value,
    na.rm = TRUE
  ),
  median_maf = median(
    gwas_clean$maf,
    na.rm = TRUE
  )
)

fwrite(
  gwas_clean,
  output_file,
  sep = "\t"
)

fwrite(
  qc_summary,
  qc_file,
  sep = "\t"
)

cat("\nQC summary:\n")
print(qc_summary)

cat("\nSaved harmonized GWAS to:\n", output_file, "\n")
cat("Saved QC summary to:\n", qc_file, "\n")