###############################################################################
# PigmentationAtlas
# Module 10.3 — Production QC, duplicate management and publication manifest
#
# Purpose:
#   Post-process Module 10.1.1 SuSiE outputs without rerunning fine-mapping.
#   - separates the regional GWAS lead from the top-PIP variant
#   - selects one representative from overlapping/duplicate loci
#   - assigns publication-level QC statuses
#   - writes production and publication manifests
#   - writes duplicate, exclusion and annotation tables
#
# Required inputs:
#   results/susie_manifest.tsv
#   results/clumps/CLUMP_xxxx/gwas/gwas_region.tsv.gz
#   results/clumps/CLUMP_xxxx/susie/susie_pip.tsv
#
# Main outputs:
#   results/susie_manifest_production.tsv
#   results/susie_manifest_publication.tsv
#   results/susie_qc_report.tsv
#   results/susie_duplicate_groups.tsv
#   results/susie_production_summary.tsv
#   results/susie_publication_figure_manifest.tsv
#   results/clumps/CLUMP_xxxx/susie/publication_annotation.tsv
#   logs/module10_3_production_qc.log
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
})

###############################################################################
# 1. USER CONFIGURATION
###############################################################################

project_dir <- "C:/Users/User/Desktop/PigmentationAtlas"
setwd(project_dir)

# Input/output files
input_manifest_file <- "results/susie_manifest.tsv"
production_manifest_file <- "results/susie_manifest_production.tsv"
publication_manifest_file <- "results/susie_manifest_publication.tsv"
qc_report_file <- "results/susie_qc_report.tsv"
duplicate_report_file <- "results/susie_duplicate_groups.tsv"
summary_file <- "results/susie_production_summary.tsv"
figure_manifest_file <- "results/susie_publication_figure_manifest.tsv"
log_file <- "logs/module10_3_production_qc.log"

# Statistical/QC thresholds
genomewide_significance_threshold <- 5e-8
estimate_s_warning_threshold <- 0.10
low_max_pip_warning <- 0.10
many_credible_sets_warning <- 10L

# Duplicate definition used only when duplicate_group is absent or incomplete.
# Two same-chromosome loci are grouped if they have the same regional GWAS lead
# OR if reciprocal interval overlap reaches this threshold.
duplicate_reciprocal_overlap_threshold <- 0.80

# Publication inclusion policy
include_pass <- TRUE
include_pass_with_warning <- TRUE
exclude_nonconverged <- TRUE
exclude_failed <- TRUE
exclude_duplicate_nonrepresentatives <- TRUE

# Figure file names produced by Module 10.2
figure_pdf_suffix <- "_susie_finemapping.pdf"
figure_png_suffix <- "_susie_finemapping.png"

# If TRUE, overwrite existing Module 10.3 outputs.
overwrite_outputs <- TRUE

###############################################################################
# 2. HELPERS
###############################################################################

dir.create("logs", recursive = TRUE, showWarnings = FALSE)

now_text <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_message <- function(...) {
  msg <- paste0("[", now_text(), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

safe_read <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) stop("Required file not found: ", path)
    return(NULL)
  }
  fread(path)
}

first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0L) NA_character_ else hit[1L]
}

pick_column <- function(dt, candidates, required = TRUE, label = "column") {
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) == 0L) {
    if (required) {
      stop("Could not identify ", label, ". Tried: ", paste(candidates, collapse = ", "))
    }
    return(NA_character_)
  }
  hit[1L]
}

as_numeric_safe <- function(x) suppressWarnings(as.numeric(x))
as_integer_safe <- function(x) suppressWarnings(as.integer(x))

clean_character <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

collapse_unique <- function(x, sep = ";") {
  x <- clean_character(x)
  x <- unique(x[!is.na(x)])
  if (length(x) == 0L) NA_character_ else paste(x, collapse = sep)
}

append_flag <- function(existing, flag) {
  flags <- unlist(strsplit(ifelse(is.na(existing), "", existing), ";", fixed = TRUE))
  flags <- trimws(c(flags, flag))
  flags <- unique(flags[nzchar(flags) & !is.na(flags)])
  if (length(flags) == 0L) NA_character_ else paste(flags, collapse = ";")
}

format_p <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 3))
}

format_num <- function(x, digits = 4L) {
  ifelse(is.na(x), "NA", formatC(x, format = "fg", digits = digits))
}

resolve_clump_files <- function(clump_id) {
  clump_dir <- file.path("results", "clumps", clump_id)
  list(
    clump_dir = clump_dir,
    gwas_file = first_existing_file(c(
      file.path(clump_dir, "gwas", "gwas_region.tsv.gz"),
      file.path(clump_dir, "gwas", "gwas_region.tsv"),
      file.path(clump_dir, "gwas_region.tsv.gz"),
      file.path(clump_dir, "gwas_region.tsv")
    )),
    pip_file = first_existing_file(c(
      file.path(clump_dir, "susie", "susie_pip.tsv"),
      file.path(clump_dir, "susie_pip.tsv")
    )),
    annotation_file = file.path(clump_dir, "susie", "publication_annotation.tsv"),
    figure_pdf = file.path(clump_dir, "figures", paste0(clump_id, figure_pdf_suffix)),
    figure_png = file.path(clump_dir, "figures", paste0(clump_id, figure_png_suffix))
  )
}

extract_gwas_lead <- function(gwas_file) {
  gwas <- safe_read(gwas_file, required = TRUE)
  if (nrow(gwas) == 0L) stop("GWAS region file is empty: ", gwas_file)

  p_col <- pick_column(gwas, c("p_value", "pval", "p", "P", "p_value_gwas"), TRUE, "GWAS p-value column")
  rsid_col <- pick_column(gwas, c("rsid", "rsID", "SNP", "variant_rsid"), FALSE, "GWAS rsID column")
  chr_col <- pick_column(gwas, c("chromosome", "chr", "CHR", "#CHROM"), FALSE, "chromosome column")
  pos_col <- pick_column(gwas, c("base_pair_location", "position", "pos", "BP"), FALSE, "position column")
  variant_col <- pick_column(gwas, c("variant_id", "ID", "reference_variant_id"), FALSE, "variant ID column")

  gwas[, .p_internal := as_numeric_safe(get(p_col))]
  usable <- gwas[is.finite(.p_internal) & .p_internal >= 0 & .p_internal <= 1]
  if (nrow(usable) == 0L) stop("No valid p-values in GWAS file: ", gwas_file)

  lead <- usable[which.min(.p_internal)]

  list(
    rsid = if (!is.na(rsid_col)) clean_character(lead[[rsid_col]])[1L] else NA_character_,
    variant_id = if (!is.na(variant_col)) clean_character(lead[[variant_col]])[1L] else NA_character_,
    chromosome = if (!is.na(chr_col)) clean_character(lead[[chr_col]])[1L] else NA_character_,
    position = if (!is.na(pos_col)) as_integer_safe(lead[[pos_col]])[1L] else NA_integer_,
    p_value = lead$.p_internal[1L],
    n_region_variants = nrow(gwas),
    n_genomewide_significant = sum(usable$.p_internal <= genomewide_significance_threshold, na.rm = TRUE)
  )
}

extract_top_pip <- function(pip_file) {
  pip <- safe_read(pip_file, required = TRUE)
  if (nrow(pip) == 0L) stop("SuSiE PIP file is empty: ", pip_file)

  pip_col <- pick_column(pip, c("pip", "PIP"), TRUE, "PIP column")
  p_col <- pick_column(pip, c("p_value", "pval", "p", "P"), FALSE, "PIP-table GWAS p-value column")
  rsid_col <- pick_column(pip, c("rsid", "rsID", "SNP"), FALSE, "PIP-table rsID column")
  variant_col <- pick_column(pip, c("reference_variant_id", "variant_id", "ID"), FALSE, "PIP-table variant ID column")
  chr_col <- pick_column(pip, c("chromosome", "chr", "CHR"), FALSE, "PIP-table chromosome column")
  pos_col <- pick_column(pip, c("position", "base_pair_location", "pos", "BP"), FALSE, "PIP-table position column")
  cs_col <- pick_column(pip, c("credible_sets", "credible_set"), FALSE, "credible set column")

  pip[, .pip_internal := as_numeric_safe(get(pip_col))]
  usable <- pip[is.finite(.pip_internal)]
  if (nrow(usable) == 0L) stop("No valid PIP values in: ", pip_file)

  top <- usable[order(-.pip_internal)][1L]

  list(
    rsid = if (!is.na(rsid_col)) clean_character(top[[rsid_col]])[1L] else NA_character_,
    variant_id = if (!is.na(variant_col)) clean_character(top[[variant_col]])[1L] else NA_character_,
    chromosome = if (!is.na(chr_col)) clean_character(top[[chr_col]])[1L] else NA_character_,
    position = if (!is.na(pos_col)) as_integer_safe(top[[pos_col]])[1L] else NA_integer_,
    p_value = if (!is.na(p_col)) as_numeric_safe(top[[p_col]])[1L] else NA_real_,
    pip = top$.pip_internal[1L],
    credible_sets = if (!is.na(cs_col)) clean_character(top[[cs_col]])[1L] else NA_character_
  )
}

# Union-find implementation for deterministic duplicate grouping.
find_root <- function(parent, i) {
  while (parent[i] != i) {
    parent[i] <- parent[parent[i]]
    i <- parent[i]
  }
  list(parent = parent, root = i)
}

union_roots <- function(parent, a, b) {
  fa <- find_root(parent, a); parent <- fa$parent; ra <- fa$root
  fb <- find_root(parent, b); parent <- fb$parent; rb <- fb$root
  if (ra != rb) parent[rb] <- ra
  parent
}

build_duplicate_groups <- function(dt) {
  n <- nrow(dt)
  if (n < 2L) return(dt)

  parent <- seq_len(n)

  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      chr_i <- clean_character(dt$chromosome[i])
      chr_j <- clean_character(dt$chromosome[j])
      if (is.na(chr_i) || is.na(chr_j) || chr_i != chr_j) next

      same_gwas_lead <- FALSE
      if (!is.na(dt$gwas_lead_rsid[i]) && !is.na(dt$gwas_lead_rsid[j])) {
        same_gwas_lead <- dt$gwas_lead_rsid[i] == dt$gwas_lead_rsid[j]
      } else if (!is.na(dt$gwas_lead_variant[i]) && !is.na(dt$gwas_lead_variant[j])) {
        same_gwas_lead <- dt$gwas_lead_variant[i] == dt$gwas_lead_variant[j]
      }

      reciprocal_overlap <- 0
      if (all(is.finite(c(dt$region_start[i], dt$region_end[i], dt$region_start[j], dt$region_end[j])))) {
        overlap <- max(0, min(dt$region_end[i], dt$region_end[j]) - max(dt$region_start[i], dt$region_start[j]) + 1)
        len_i <- dt$region_end[i] - dt$region_start[i] + 1
        len_j <- dt$region_end[j] - dt$region_start[j] + 1
        reciprocal_overlap <- min(overlap / len_i, overlap / len_j)
      }

      same_existing_group <- !is.na(dt$duplicate_group[i]) &&
        !is.na(dt$duplicate_group[j]) &&
        dt$duplicate_group[i] == dt$duplicate_group[j]

      if (same_gwas_lead || same_existing_group ||
          reciprocal_overlap >= duplicate_reciprocal_overlap_threshold) {
        parent <- union_roots(parent, i, j)
      }
    }
  }

  roots <- integer(n)
  for (i in seq_len(n)) {
    fr <- find_root(parent, i)
    parent <- fr$parent
    roots[i] <- fr$root
  }

  counts <- table(roots)
  duplicate_roots <- as.integer(names(counts[counts > 1L]))

  dt[, `:=`(
    duplicate_flag = FALSE,
    duplicate_group = NA_character_
  )]

  if (length(duplicate_roots) > 0L) {
    ordered_roots <- duplicate_roots[order(vapply(duplicate_roots, function(r) {
      min(dt$clump_id[roots == r])
    }, character(1L)))]

    for (k in seq_along(ordered_roots)) {
      members <- which(roots == ordered_roots[k])
      dt[members, `:=`(
        duplicate_flag = TRUE,
        duplicate_group = sprintf("DUPGRP_%04d", k)
      )]
    }
  }

  dt
}

select_duplicate_representatives <- function(dt) {
  dt[, `:=`(
    representative_locus = TRUE,
    representative_clump_id = clump_id,
    duplicate_selection_reason = NA_character_
  )]

  groups <- unique(na.omit(dt$duplicate_group))
  if (length(groups) == 0L) return(dt)

  for (grp in groups) {
    idx <- which(dt$duplicate_group == grp)
    candidates <- copy(dt[idx])

    candidates[, `:=`(
      rank_converged = fifelse(!is.na(converged) & converged, 0L, 1L),
      rank_estimate_s = fifelse(is.finite(estimate_s), estimate_s, Inf),
      rank_pip = fifelse(is.finite(top_pip), -top_pip, Inf),
      rank_n_variants = fifelse(is.finite(n_susie_variants), -n_susie_variants, Inf),
      rank_width = fifelse(is.finite(region_width), -region_width, Inf)
    )]

    setorder(
      candidates,
      rank_converged,
      rank_estimate_s,
      rank_pip,
      rank_n_variants,
      rank_width,
      clump_id
    )

    representative_id <- candidates$clump_id[1L]
    dt[idx, `:=`(
      representative_locus = clump_id == representative_id,
      representative_clump_id = representative_id
    )]

    reason <- paste0(
      "priority=converged_then_lower_estimate_s_then_higher_PIP_then_more_variants_then_wider_region;",
      "representative=", representative_id
    )
    dt[idx, duplicate_selection_reason := reason]
  }

  dt
}

assign_publication_status <- function(dt) {
  dt[, `:=`(
    production_qc_flags = clean_character(qc_flags),
    publication_status = NA_character_,
    publication_reason = NA_character_,
    include_in_publication = FALSE
  )]

  for (i in seq_len(nrow(dt))) {
    flags <- character()

    status_i <- clean_character(dt$status[i])
    converged_i <- isTRUE(dt$converged[i])
    duplicate_nonrep <- isTRUE(dt$duplicate_flag[i]) && !isTRUE(dt$representative_locus[i])

    if (duplicate_nonrep && exclude_duplicate_nonrepresentatives) {
      pub_status <- "EXCLUDED_DUPLICATE"
      flags <- c(flags, "duplicate_nonrepresentative")
    } else if ((is.na(status_i) || status_i %in% c("failed", "missing", "error")) && exclude_failed) {
      pub_status <- "EXCLUDED_FAILED"
      flags <- c(flags, ifelse(is.na(status_i), "missing_susie_status", paste0("susie_status_", status_i)))
    } else if (!converged_i && exclude_nonconverged) {
      pub_status <- "EXCLUDED_NONCONVERGED"
      flags <- c(flags, "nonconverged")
    } else {
      if (!is.na(dt$estimate_s[i]) && dt$estimate_s[i] > estimate_s_warning_threshold) {
        flags <- c(flags, "gwas_ld_inconsistency")
      }
      if (!is.na(dt$top_pip[i]) && dt$top_pip[i] < low_max_pip_warning) {
        flags <- c(flags, "low_max_pip")
      }
      if (!is.na(dt$n_credible_sets[i]) && dt$n_credible_sets[i] > many_credible_sets_warning) {
        flags <- c(flags, "many_credible_sets")
      }
      if (is.na(dt$n_credible_sets[i]) || dt$n_credible_sets[i] == 0L) {
        flags <- c(flags, "no_credible_set")
      }
      if (!is.na(dt$gwas_lead_p[i]) && dt$gwas_lead_p[i] > genomewide_significance_threshold) {
        flags <- c(flags, "region_not_genomewide_significant")
      }
      if (!is.na(dt$top_pip_p[i]) && dt$top_pip_p[i] > genomewide_significance_threshold &&
          !is.na(dt$gwas_lead_p[i]) && dt$gwas_lead_p[i] <= genomewide_significance_threshold) {
        flags <- c(flags, "top_pip_variant_not_significant_despite_significant_region")
      }
      if (!is.na(dt$gwas_lead_rsid[i]) && !is.na(dt$top_pip_rsid[i]) &&
          dt$gwas_lead_rsid[i] != dt$top_pip_rsid[i]) {
        flags <- c(flags, "top_pip_differs_from_gwas_lead")
      }

      # Preserve informative Module 10.1.1 flags except the old misleading name.
      old_flags <- unlist(strsplit(ifelse(is.na(dt$qc_flags[i]), "", dt$qc_flags[i]), ";", fixed = TRUE))
      old_flags <- trimws(old_flags)
      old_flags <- old_flags[nzchar(old_flags)]
      old_flags <- setdiff(old_flags, c(
        "max_pip_variant_not_genomewide_significant",
        "overlapping_or_repeated_locus",
        "nonconverged"
      ))
      flags <- c(flags, old_flags)

      flags <- unique(flags)
      pub_status <- if (length(flags) == 0L) "PASS" else "PASS_WITH_WARNING"
    }

    flags <- unique(flags)
    reason <- if (length(flags) == 0L) NA_character_ else paste(flags, collapse = ";")

    include <- (pub_status == "PASS" && include_pass) ||
      (pub_status == "PASS_WITH_WARNING" && include_pass_with_warning)

    dt[i, `:=`(
      production_qc_flags = reason,
      publication_status = pub_status,
      publication_reason = reason,
      include_in_publication = include
    )]
  }

  dt
}

write_annotation <- function(row, paths) {
  dir.create(dirname(paths$annotation_file), recursive = TRUE, showWarnings = FALSE)

  annotation <- data.table(
    clump_id = row$clump_id,
    publication_status = row$publication_status,
    include_in_publication = row$include_in_publication,
    publication_reason = row$publication_reason,
    chromosome = row$chromosome,
    region_start = row$region_start,
    region_end = row$region_end,
    gwas_lead_rsid = row$gwas_lead_rsid,
    gwas_lead_variant = row$gwas_lead_variant,
    gwas_lead_position = row$gwas_lead_position,
    gwas_lead_p = row$gwas_lead_p,
    top_pip_rsid = row$top_pip_rsid,
    top_pip_variant = row$top_pip_variant,
    top_pip_position = row$top_pip_position,
    top_pip = row$top_pip,
    top_pip_p = row$top_pip_p,
    converged = row$converged,
    n_iterations = row$n_iterations,
    estimate_s = row$estimate_s,
    n_credible_sets = row$n_credible_sets,
    duplicate_flag = row$duplicate_flag,
    duplicate_group = row$duplicate_group,
    representative_locus = row$representative_locus,
    representative_clump_id = row$representative_clump_id,
    figure_pdf = paths$figure_pdf,
    figure_pdf_exists = file.exists(paths$figure_pdf),
    figure_png = paths$figure_png,
    figure_png_exists = file.exists(paths$figure_png),
    generated_at = now_text()
  )

  fwrite(annotation, paths$annotation_file, sep = "\t", na = "NA")
  annotation
}

###############################################################################
# 3. LOAD MODULE 10 MANIFEST
###############################################################################

if (!file.exists(input_manifest_file)) {
  stop("Input SuSiE manifest not found: ", input_manifest_file)
}

if (!overwrite_outputs) {
  existing <- c(production_manifest_file, publication_manifest_file, qc_report_file,
                duplicate_report_file, summary_file, figure_manifest_file)
  if (any(file.exists(existing))) {
    stop("Module 10.3 outputs already exist and overwrite_outputs=FALSE.")
  }
}

cat("", file = log_file)
log_message("Module 10.3 started.")

manifest <- fread(input_manifest_file)
if (!"clump_id" %in% names(manifest)) stop("Input manifest has no clump_id column.")
if (anyDuplicated(manifest$clump_id)) stop("Input manifest contains duplicated clump_id values.")

# Ensure columns needed downstream exist.
required_defaults <- list(
  status = NA_character_,
  reason = NA_character_,
  chromosome = NA_character_,
  region_start = NA_integer_,
  region_end = NA_integer_,
  n_susie_variants = NA_integer_,
  n_credible_sets = NA_integer_,
  max_pip = NA_real_,
  lead_variant = NA_character_,
  lead_rsid = NA_character_,
  lead_position = NA_integer_,
  lead_p_value = NA_real_,
  converged = NA,
  n_iterations = NA_integer_,
  estimate_s = NA_real_,
  estimate_s_method = NA_character_,
  runtime_seconds = NA_real_,
  qc_status = NA_character_,
  qc_flags = NA_character_,
  duplicate_group = NA_character_,
  duplicate_flag = FALSE
)
for (nm in names(required_defaults)) {
  if (!nm %in% names(manifest)) manifest[, (nm) := required_defaults[[nm]]]
}

manifest[, `:=`(
  chromosome = clean_character(chromosome),
  region_start = as_integer_safe(region_start),
  region_end = as_integer_safe(region_end),
  n_susie_variants = as_integer_safe(n_susie_variants),
  n_credible_sets = as_integer_safe(n_credible_sets),
  max_pip = as_numeric_safe(max_pip),
  converged = as.logical(converged),
  n_iterations = as_integer_safe(n_iterations),
  estimate_s = as_numeric_safe(estimate_s),
  runtime_seconds = as_numeric_safe(runtime_seconds),
  qc_flags = clean_character(qc_flags),
  duplicate_group = clean_character(duplicate_group),
  duplicate_flag = as.logical(duplicate_flag)
)]

log_message("Loaded ", nrow(manifest), " loci from ", input_manifest_file, ".")

###############################################################################
# 4. RECONSTRUCT GWAS LEAD AND TOP-PIP VARIANT
###############################################################################

reconstructed <- vector("list", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
  clump_id <- manifest$clump_id[i]
  paths <- resolve_clump_files(clump_id)
  log_message("Inspecting ", clump_id, " | ", i, "/", nrow(manifest))

  row <- data.table(
    clump_id = clump_id,
    gwas_file = paths$gwas_file,
    pip_file = paths$pip_file,
    gwas_lead_rsid = NA_character_,
    gwas_lead_variant = NA_character_,
    gwas_lead_chromosome = NA_character_,
    gwas_lead_position = NA_integer_,
    gwas_lead_p = NA_real_,
    n_gwas_region_recounted = NA_integer_,
    n_gwas_significant = NA_integer_,
    top_pip_rsid = clean_character(manifest$lead_rsid[i]),
    top_pip_variant = clean_character(manifest$lead_variant[i]),
    top_pip_chromosome = clean_character(manifest$chromosome[i]),
    top_pip_position = as_integer_safe(manifest$lead_position[i]),
    top_pip = as_numeric_safe(manifest$max_pip[i]),
    top_pip_p = as_numeric_safe(manifest$lead_p_value[i]),
    top_pip_credible_sets = NA_character_,
    reconstruction_status = "ready",
    reconstruction_error = NA_character_
  )

  errors <- character()

  if (is.na(paths$gwas_file)) {
    errors <- c(errors, "missing_gwas_region_file")
  } else {
    gwas_lead <- tryCatch(extract_gwas_lead(paths$gwas_file), error = function(e) e)
    if (inherits(gwas_lead, "error")) {
      errors <- c(errors, paste0("gwas_lead_error:", conditionMessage(gwas_lead)))
    } else {
      set(row, i = 1L, j = "gwas_lead_rsid", value = gwas_lead$rsid)
      set(row, i = 1L, j = "gwas_lead_variant", value = gwas_lead$variant_id)
      set(row, i = 1L, j = "gwas_lead_chromosome", value = gwas_lead$chromosome)
      set(row, i = 1L, j = "gwas_lead_position", value = gwas_lead$position)
      set(row, i = 1L, j = "gwas_lead_p", value = gwas_lead$p_value)
      set(row, i = 1L, j = "n_gwas_region_recounted", value = gwas_lead$n_region_variants)
      set(row, i = 1L, j = "n_gwas_significant", value = gwas_lead$n_genomewide_significant)
    }
  }

  if (is.na(paths$pip_file)) {
    errors <- c(errors, "missing_susie_pip_file")
  } else {
    top_pip_info <- tryCatch(extract_top_pip(paths$pip_file), error = function(e) e)
    if (inherits(top_pip_info, "error")) {
      errors <- c(errors, paste0("top_pip_error:", conditionMessage(top_pip_info)))
    } else {
      # Use set() rather than := so the existing top_pip column cannot mask
      # the local top_pip_info object during data.table evaluation.
      set(row, i = 1L, j = "top_pip_rsid", value = top_pip_info$rsid)
      set(row, i = 1L, j = "top_pip_variant", value = top_pip_info$variant_id)
      set(row, i = 1L, j = "top_pip_chromosome", value = top_pip_info$chromosome)
      set(row, i = 1L, j = "top_pip_position", value = top_pip_info$position)
      set(row, i = 1L, j = "top_pip", value = top_pip_info$pip)
      set(row, i = 1L, j = "top_pip_p", value = top_pip_info$p_value)
      set(row, i = 1L, j = "top_pip_credible_sets", value = top_pip_info$credible_sets)
    }
  }

  if (length(errors) > 0L) {
    row[, `:=`(
      reconstruction_status = "warning",
      reconstruction_error = paste(errors, collapse = ";")
    )]
  }

  reconstructed[[i]] <- row
}

reconstructed <- rbindlist(reconstructed, use.names = TRUE, fill = TRUE)
production <- merge(manifest, reconstructed, by = "clump_id", all.x = TRUE, sort = FALSE)
production[, .input_order := match(clump_id, manifest$clump_id)]
setorder(production, .input_order)
production[, .input_order := NULL]

# Use reconstructed chromosome where the original manifest value is absent.
production[is.na(chromosome), chromosome := gwas_lead_chromosome]
production[, region_width := fifelse(
  is.finite(region_start) & is.finite(region_end),
  as.integer(region_end - region_start + 1L),
  NA_integer_
)]

# Retain backward compatibility but make the old column semantics explicit.
production[, `:=`(
  legacy_lead_rsid_is_top_pip = TRUE,
  legacy_lead_variant_is_top_pip = TRUE
)]

###############################################################################
# 5. DUPLICATE MANAGEMENT
###############################################################################

production <- build_duplicate_groups(production)
production <- select_duplicate_representatives(production)

###############################################################################
# 6. PUBLICATION-LEVEL QC
###############################################################################

production <- assign_publication_status(production)

# Add compact status rank for sorting.
production[, publication_status_rank := fcase(
  publication_status == "PASS", 1L,
  publication_status == "PASS_WITH_WARNING", 2L,
  publication_status == "EXCLUDED_NONCONVERGED", 3L,
  publication_status == "EXCLUDED_DUPLICATE", 4L,
  publication_status == "EXCLUDED_FAILED", 5L,
  default = 9L
)]

###############################################################################
# 7. ANNOTATION AND FIGURE MANIFEST
###############################################################################

annotation_rows <- vector("list", nrow(production))
for (i in seq_len(nrow(production))) {
  paths <- resolve_clump_files(production$clump_id[i])
  annotation_rows[[i]] <- write_annotation(production[i], paths)
}
figure_manifest <- rbindlist(annotation_rows, use.names = TRUE, fill = TRUE)

###############################################################################
# 8. REPORT TABLES
###############################################################################

publication_manifest <- production[include_in_publication == TRUE]
setorder(publication_manifest, chromosome, region_start, clump_id)

qc_report <- production[, .(
  publication_status_rank,
  clump_id,
  chromosome,
  region_start,
  region_end,
  gwas_lead_rsid,
  gwas_lead_p,
  top_pip_rsid,
  top_pip,
  top_pip_p,
  converged,
  n_iterations,
  estimate_s,
  n_credible_sets,
  duplicate_flag,
  duplicate_group,
  representative_locus,
  representative_clump_id,
  publication_status,
  publication_reason,
  include_in_publication,
  reconstruction_status,
  reconstruction_error
)]
setorder(qc_report, publication_status_rank, clump_id)
qc_report[, publication_status_rank := NULL]

duplicate_groups <- production[duplicate_flag == TRUE, .(
  group_size = .N,
  representative_clump_id = unique(representative_clump_id),
  member_clumps = paste(clump_id, collapse = ";"),
  excluded_members = collapse_unique(clump_id[!representative_locus]),
  chromosome = collapse_unique(chromosome),
  group_start = suppressWarnings(min(region_start, na.rm = TRUE)),
  group_end = suppressWarnings(max(region_end, na.rm = TRUE)),
  gwas_leads = collapse_unique(gwas_lead_rsid),
  strongest_gwas_p = suppressWarnings(min(gwas_lead_p, na.rm = TRUE)),
  selection_rule = unique(duplicate_selection_reason)
), by = duplicate_group]

if (nrow(duplicate_groups) > 0L) {
  duplicate_groups[!is.finite(group_start), group_start := NA_integer_]
  duplicate_groups[!is.finite(group_end), group_end := NA_integer_]
  duplicate_groups[!is.finite(strongest_gwas_p), strongest_gwas_p := NA_real_]
}

status_summary <- production[, .N, by = publication_status][order(publication_status)]
summary_metrics <- data.table(
  metric = c(
    "run_completed_at",
    "input_loci",
    "publication_loci",
    "excluded_loci",
    "duplicate_groups",
    "duplicate_loci",
    "converged_loci",
    "nonconverged_loci",
    "median_estimate_s_all",
    "median_estimate_s_included",
    "median_max_pip_included",
    "median_credible_sets_included",
    "median_runtime_seconds_all"
  ),
  value = as.character(c(
    now_text(),
    nrow(production),
    nrow(publication_manifest),
    sum(!production$include_in_publication),
    uniqueN(na.omit(production$duplicate_group)),
    sum(production$duplicate_flag, na.rm = TRUE),
    sum(production$converged %in% TRUE, na.rm = TRUE),
    sum(production$converged %in% FALSE, na.rm = TRUE),
    median(production$estimate_s, na.rm = TRUE),
    median(publication_manifest$estimate_s, na.rm = TRUE),
    median(publication_manifest$top_pip, na.rm = TRUE),
    median(publication_manifest$n_credible_sets, na.rm = TRUE),
    median(production$runtime_seconds, na.rm = TRUE)
  ))
)

status_summary_long <- status_summary[, .(
  metric = paste0("status_", publication_status),
  value = as.character(N)
)]
summary_table <- rbindlist(list(summary_metrics, status_summary_long), use.names = TRUE, fill = TRUE)

###############################################################################
# 9. WRITE OUTPUTS
###############################################################################

fwrite(production, production_manifest_file, sep = "\t", na = "NA")
fwrite(publication_manifest, publication_manifest_file, sep = "\t", na = "NA")
fwrite(qc_report, qc_report_file, sep = "\t", na = "NA")
fwrite(duplicate_groups, duplicate_report_file, sep = "\t", na = "NA")
fwrite(summary_table, summary_file, sep = "\t", na = "NA")
fwrite(figure_manifest, figure_manifest_file, sep = "\t", na = "NA")

###############################################################################
# 10. CONSOLE SUMMARY
###############################################################################

log_message("Module 10.3 finished.")
log_message("Production loci: ", nrow(production))
log_message("Publication loci: ", nrow(publication_manifest))
log_message("Duplicate groups: ", uniqueN(na.omit(production$duplicate_group)))
log_message(
  "Status counts: ",
  paste(status_summary[, paste0(publication_status, "=", N)], collapse = "; ")
)

cat("\n============================================================\n")
cat("MODULE 10.3 COMPLETE\n")
cat("============================================================\n")
print(status_summary)
cat("\nInput loci:       ", nrow(production), "\n", sep = "")
cat("Publication loci: ", nrow(publication_manifest), "\n", sep = "")
cat("Excluded loci:    ", sum(!production$include_in_publication), "\n", sep = "")
cat("Duplicate groups: ", uniqueN(na.omit(production$duplicate_group)), "\n", sep = "")
cat("\nProduction manifest: ", production_manifest_file, "\n", sep = "")
cat("Publication manifest: ", publication_manifest_file, "\n", sep = "")
cat("QC report:            ", qc_report_file, "\n", sep = "")
cat("Duplicate report:     ", duplicate_report_file, "\n", sep = "")
cat("Summary:              ", summary_file, "\n", sep = "")
cat("Figure manifest:      ", figure_manifest_file, "\n", sep = "")
cat("Log:                  ", log_file, "\n", sep = "")
cat("============================================================\n")
