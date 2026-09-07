
############################################################
# PigmentationAtlas
# Module 06
# Initialize locus workspaces and extract GWAS regions
############################################################

suppressPackageStartupMessages({
  library(data.table)
})

trait_id <- "GCST90691600"

gwas_file <- file.path(
  "results",
  paste0(trait_id, "_harmonized.tsv.gz")
)

loci_file <- file.path(
  "results",
  paste0(trait_id, "_genomic_loci_v3.tsv")
)

output_root <- file.path("results", "loci")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

log_message <- function(...) {
  cat(
    format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"),
    ...,
    "\n"
  )
}

required_gwas_columns <- c(
  "chromosome",
  "base_pair_location",
  "effect_allele",
  "other_allele",
  "beta",
  "standard_error",
  "effect_allele_frequency",
  "p_value",
  "rsid",
  "variant_id",
  "maf",
  "z_score"
)

required_locus_columns <- c(
  "locus_id",
  "chr",
  "start",
  "end",
  "lead_id",
  "lead_pos",
  "lead_p"
)

log_message("Reading GWAS:", gwas_file)
gwas <- fread(gwas_file)

log_message("Reading loci:", loci_file)
loci <- fread(loci_file)

missing_gwas <- setdiff(required_gwas_columns, names(gwas))
missing_loci <- setdiff(required_locus_columns, names(loci))

if (length(missing_gwas) > 0L) {
  stop(
    "Missing GWAS columns: ",
    paste(missing_gwas, collapse = ", ")
  )
}

if (length(missing_loci) > 0L) {
  stop(
    "Missing locus columns: ",
    paste(missing_loci, collapse = ", ")
  )
}

if (anyDuplicated(loci$locus_id)) {
  stop("Duplicated locus_id values found.")
}

if (any(loci$start > loci$end, na.rm = TRUE)) {
  stop("At least one locus has start > end.")
}

setkey(gwas, chromosome, base_pair_location)

summary_list <- vector("list", nrow(loci))

log_message(
  "Starting extraction for",
  nrow(loci),
  "loci"
)

for (i in seq_len(nrow(loci))) {

  locus_row <- loci[i]

  locus_id <- locus_row$locus_id
  locus_chr <- as.integer(locus_row$chr)
  locus_start <- as.integer(locus_row$start)
  locus_end <- as.integer(locus_row$end)

  locus_dir <- file.path(output_root, locus_id)

  dir.create(
    locus_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  region <- gwas[
    chromosome == locus_chr &
      base_pair_location >= locus_start &
      base_pair_location <= locus_end
  ]

  setorder(region, base_pair_location, variant_id)

  lead_present <- locus_row$lead_pos %in% region$base_pair_location

  top_variant <- if (nrow(region) > 0L) {
    region[which.min(p_value)]
  } else {
    NULL
  }

  metadata <- copy(locus_row)

  metadata[, `:=`(
    trait_id = trait_id,
    analysis_start = locus_start,
    analysis_end = locus_end,
    analysis_size_bp = locus_end - locus_start + 1L,
    n_gwas_variants = nrow(region),
    n_unique_positions = uniqueN(region$base_pair_location),
    n_unique_variant_ids = uniqueN(region$variant_id),
    minimum_gwas_position = if (nrow(region) > 0L) {
      min(region$base_pair_location)
    } else {
      NA_integer_
    },
    maximum_gwas_position = if (nrow(region) > 0L) {
      max(region$base_pair_location)
    } else {
      NA_integer_
    },
    lead_position_present = lead_present,
    minimum_p_in_region = if (nrow(region) > 0L) {
      min(region$p_value, na.rm = TRUE)
    } else {
      NA_real_
    },
    top_variant_in_region = if (nrow(region) > 0L) {
      top_variant$variant_id
    } else {
      NA_character_
    },
    top_rsid_in_region = if (nrow(region) > 0L) {
      top_variant$rsid
    } else {
      NA_character_
    },
    top_position_in_region = if (nrow(region) > 0L) {
      top_variant$base_pair_location
    } else {
      NA_integer_
    }
  )]

  fwrite(
    metadata,
    file.path(locus_dir, "metadata.tsv"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  fwrite(
    region,
    file.path(locus_dir, "gwas.tsv.gz"),
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = "gzip"
  )

  summary_list[[i]] <- metadata[, .(
    trait_id,
    locus_id,
    chr,
    analysis_start,
    analysis_end,
    analysis_size_bp,
    lead_id,
    lead_pos,
    lead_p,
    n_clumps,
    n_gwas_variants,
    n_unique_positions,
    n_unique_variant_ids,
    minimum_gwas_position,
    maximum_gwas_position,
    minimum_p_in_region,
    top_variant_in_region,
    top_rsid_in_region,
    top_position_in_region,
    lead_position_present
  )]

  if (
    i == 1L ||
    i %% 10L == 0L ||
    i == nrow(loci)
  ) {
    log_message(
      "Processed",
      i,
      "of",
      nrow(loci),
      "loci"
    )
  }
}

extraction_summary <- rbindlist(
  summary_list,
  use.names = TRUE,
  fill = TRUE
)

setorder(extraction_summary, chr, analysis_start)

summary_file <- file.path(
  output_root,
  paste0(
    trait_id,
    "_locus_extraction_summary.tsv"
  )
)

fwrite(
  extraction_summary,
  summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

qc_summary <- data.table(
  trait_id = trait_id,
  n_loci = nrow(extraction_summary),
  total_extracted_variants =
    sum(extraction_summary$n_gwas_variants),
  minimum_variants_per_locus =
    min(extraction_summary$n_gwas_variants),
  median_variants_per_locus =
    median(extraction_summary$n_gwas_variants),
  maximum_variants_per_locus =
    max(extraction_summary$n_gwas_variants),
  loci_with_zero_variants =
    sum(extraction_summary$n_gwas_variants == 0L),
  loci_missing_lead_position =
    sum(!extraction_summary$lead_position_present),
  loci_where_top_position_differs_from_lead =
    sum(
      extraction_summary$top_position_in_region !=
        extraction_summary$lead_pos,
      na.rm = TRUE
    )
)

qc_file <- file.path(
  output_root,
  paste0(
    trait_id,
    "_locus_extraction_QC.tsv"
  )
)

fwrite(
  qc_summary,
  qc_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

log_message("Module 06 completed.")
log_message("Summary:", summary_file)
log_message("QC:", qc_file)

print(qc_summary)

