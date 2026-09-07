###############################################################################
# PigmentationAtlas
# Module 10.1 — Robust Bayesian fine-mapping with SuSiE-RSS
#
# Inputs:
#   results/clump_ld_manifest.tsv
#   results/clumps/CLUMP_xxxx/gwas/gwas_region.tsv.gz
#   results/clumps/CLUMP_xxxx/ld/ld_variant_map.tsv
#   results/clumps/CLUMP_xxxx/ld/ld_matrix.unphased.vcor1.bin
#
# Outputs per clump:
#   results/clumps/CLUMP_xxxx/susie/susie_fit.rds
#   results/clumps/CLUMP_xxxx/susie/susie_pip.tsv
#   results/clumps/CLUMP_xxxx/susie/susie_credible_sets.tsv
#   results/clumps/CLUMP_xxxx/susie/susie_summary.tsv
#   results/clumps/CLUMP_xxxx/susie/susie_diagnostics.tsv
#
# Global outputs:
#   results/susie_manifest.tsv
#   logs/module10_susie.log
#
# Compatible with:
#   susieR 0.14.2
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(susieR)
})

###############################################################################
# 1. USER CONFIGURATION
###############################################################################

project_dir <- "C:/Users/User/Desktop/PigmentationAtlas"
setwd(project_dir)

# "pilot" = first pilot_n ready clumps
# "full"  = all ready clumps
run_mode <- "pilot"
pilot_n  <- 5L

# GWAS sample size
gwas_sample_size <- 419469L

# SuSiE settings
susie_L                         <- 10L
susie_coverage                  <- 0.95
susie_min_abs_corr              <- 0.5
susie_prior_variance            <- 50
susie_estimate_residual_variance <- FALSE
susie_max_iter                  <- 500L
susie_tol                       <- 1e-4

# Existing outputs
overwrite_existing <- TRUE

# Minimum number of usable variants
min_variants_for_susie <- 10L

# QC thresholds
genomewide_significance_threshold <- 5e-8
max_credible_sets_warning <- 10L
low_max_pip_warning <- 0.10
estimate_s_warning_threshold <- 0.10
duplicate_reciprocal_overlap_threshold <- 0.80

# Numerical tolerances
symmetry_tolerance      <- 1e-8
diagonal_tolerance      <- 1e-6
correlation_tolerance   <- 1e-6
psd_tolerance           <- -1e-6

# Matrix validation
# Full eigenvalue decomposition can be expensive for large matrices.
# It is performed only below this variant count.
max_variants_for_full_eigencheck <- 5000L

# LD/GWAS consistency diagnostic
# estimate_s_rss() is useful but may also be expensive.
run_estimate_s_rss <- TRUE
max_variants_for_estimate_s <- 5000L

# Input/output paths
ld_manifest_file <- "results/clump_ld_manifest.tsv"
susie_manifest_file <- "results/susie_manifest.tsv"
log_dir <- "logs"
log_file <- file.path(log_dir, "module10_susie.log")

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 2. HELPERS
###############################################################################

timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

log_message <- function(..., append = TRUE) {
  msg <- paste0("[", timestamp(), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = append)
}

safe_error_text <- function(e) {
  paste(class(e)[1], conditionMessage(e), sep = ": ")
}

first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0L) NA_character_ else hit[1L]
}

normalize_variant_id <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x <- sub("_b3[78]$", "", x, ignore.case = TRUE)
  x <- gsub(":", "_", x, fixed = TRUE)
  x
}

parse_reference_variant_id <- function(x) {
  x <- as.character(x)
  parts <- tstrsplit(x, ":", fixed = TRUE)

  if (length(parts) < 4L) {
    stop("reference_variant_id cannot be parsed as CHR:POS:REF:ALT.")
  }

  data.table(
    reference_chromosome = parts[[1L]],
    reference_position = suppressWarnings(as.integer(parts[[2L]])),
    reference_allele = toupper(parts[[3L]]),
    alternate_allele = toupper(parts[[4L]])
  )
}

read_plink_square_binary <- function(path, n_variants) {
  if (!file.exists(path)) {
    stop("LD binary file does not exist: ", path)
  }

  expected_values <- as.double(n_variants) * as.double(n_variants)
  expected_bytes <- expected_values * 4

  actual_bytes <- file.info(path)$size

  if (is.na(actual_bytes)) {
    stop("Cannot determine LD file size: ", path)
  }

  if (actual_bytes != expected_bytes) {
    stop(
      "Unexpected LD file size. Expected ", format(expected_bytes, scientific = FALSE),
      " bytes for a ", n_variants, " x ", n_variants,
      " single-precision square matrix, observed ",
      format(actual_bytes, scientific = FALSE), " bytes."
    )
  }

  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)

  values <- readBin(
    con,
    what = "numeric",
    n = expected_values,
    size = 4L,
    endian = "little"
  )

  if (length(values) != expected_values) {
    stop(
      "Incomplete LD matrix read: expected ", expected_values,
      " values, obtained ", length(values), "."
    )
  }

  # R fills matrices column-wise. For a symmetric square matrix, row-major and
  # column-major representations yield the same reconstructed matrix.
  matrix(values, nrow = n_variants, ncol = n_variants)
}

validate_ld_matrix <- function(R) {
  if (!is.matrix(R)) stop("R is not a matrix.")
  if (nrow(R) != ncol(R)) stop("R is not square.")
  if (any(!is.finite(R))) stop("R contains non-finite values.")

  max_asymmetry <- max(abs(R - t(R)))
  max_diagonal_deviation <- max(abs(diag(R) - 1))
  min_correlation <- min(R)
  max_correlation <- max(R)

  if (max_asymmetry > symmetry_tolerance) {
    stop("LD matrix is not symmetric; max asymmetry = ", signif(max_asymmetry, 5))
  }

  if (max_diagonal_deviation > diagonal_tolerance) {
    stop(
      "LD matrix diagonal differs from 1; max deviation = ",
      signif(max_diagonal_deviation, 5)
    )
  }

  if (min_correlation < -1 - correlation_tolerance ||
      max_correlation > 1 + correlation_tolerance) {
    stop(
      "LD matrix contains correlations outside [-1,1]: min = ",
      signif(min_correlation, 5), ", max = ", signif(max_correlation, 5)
    )
  }

  # Remove tiny floating-point asymmetry and diagonal drift.
  R <- (R + t(R)) / 2
  diag(R) <- 1

  min_eigenvalue <- NA_real_
  psd_check <- "not_run_large_matrix"

  if (nrow(R) <= max_variants_for_full_eigencheck) {
    eig <- eigen(R, symmetric = TRUE)
    min_eigenvalue <- min(eig$values)

    if (min_eigenvalue < psd_tolerance) {
      stop(
        "LD matrix is not positive semidefinite within tolerance; minimum eigenvalue = ",
        signif(min_eigenvalue, 6)
      )
    }

    # Float32 PLINK matrices can contain tiny negative eigenvalues caused by
    # numerical rounding. Clip only these tolerance-level values to zero,
    # reconstruct the matrix, and renormalize it to a correlation matrix.
    if (min_eigenvalue < 0) {
      eig$values[eig$values < 0] <- 0
      R <- eig$vectors %*% (eig$values * t(eig$vectors))
      scale_vec <- sqrt(pmax(diag(R), .Machine$double.eps))
      R <- R / tcrossprod(scale_vec)
      R <- (R + t(R)) / 2
      diag(R) <- 1
      psd_check <- "corrected_tiny_negative_eigenvalues"
    } else {
      psd_check <- "passed"
    }
  }

  list(
    R = R,
    max_asymmetry = max_asymmetry,
    max_diagonal_deviation = max_diagonal_deviation,
    min_correlation = min_correlation,
    max_correlation = max_correlation,
    min_eigenvalue = min_eigenvalue,
    psd_check = psd_check
  )
}

resolve_input_files <- function(clump_id) {
  clump_dir <- file.path("results", "clumps", clump_id)

  list(
    clump_dir = clump_dir,
    gwas_file = first_existing_file(c(
      file.path(clump_dir, "gwas", "gwas_region.tsv.gz"),
      file.path(clump_dir, "gwas", "gwas_region.tsv"),
      file.path(clump_dir, "gwas_region.tsv.gz"),
      file.path(clump_dir, "gwas_region.tsv")
    )),
    ld_map_file = first_existing_file(c(
      file.path(clump_dir, "ld", "ld_variant_map.tsv"),
      file.path(clump_dir, "ld_variant_map.tsv")
    )),
    ld_matrix_file = first_existing_file(c(
      file.path(clump_dir, "ld", "ld_matrix.unphased.vcor1.bin"),
      file.path(clump_dir, "ld_matrix.unphased.vcor1.bin")
    )),
    ld_vars_file = first_existing_file(c(
      file.path(clump_dir, "ld", "ld_matrix.unphased.vcor1.bin.vars"),
      file.path(clump_dir, "ld_matrix.unphased.vcor1.bin.vars")
    ))
  )
}

orient_gwas_to_reference <- function(gwas, ld_map) {
  required_gwas <- c(
    "variant_id", "effect_allele", "other_allele",
    "beta", "standard_error"
  )
  missing_gwas <- setdiff(required_gwas, names(gwas))
  if (length(missing_gwas) > 0L) {
    stop("Missing GWAS columns: ", paste(missing_gwas, collapse = ", "))
  }

  required_map <- c(
    "variant_id", "reference_variant_id", "ld_order"
  )
  missing_map <- setdiff(required_map, names(ld_map))
  if (length(missing_map) > 0L) {
    stop("Missing LD map columns: ", paste(missing_map, collapse = ", "))
  }

  gwas <- copy(gwas)
  ld_map <- copy(ld_map)

  gwas[, variant_id_norm := normalize_variant_id(variant_id)]
  ld_map[, variant_id_norm := normalize_variant_id(variant_id)]

  if (anyDuplicated(gwas$variant_id_norm)) {
    duplicate_ids <- unique(gwas$variant_id_norm[duplicated(gwas$variant_id_norm)])
    stop(
      "Duplicated normalized GWAS variant IDs; examples: ",
      paste(head(duplicate_ids, 5L), collapse = ", ")
    )
  }

  if (anyDuplicated(ld_map$variant_id_norm)) {
    duplicate_ids <- unique(ld_map$variant_id_norm[duplicated(ld_map$variant_id_norm)])
    stop(
      "Duplicated normalized LD-map variant IDs; examples: ",
      paste(head(duplicate_ids, 5L), collapse = ", ")
    )
  }

  ref_parts <- parse_reference_variant_id(ld_map$reference_variant_id)
  ld_map <- cbind(ld_map, ref_parts)

  merged <- merge(
    ld_map,
    gwas,
    by = "variant_id_norm",
    all.x = TRUE,
    sort = FALSE,
    suffixes = c("_ldmap", "_gwas")
  )

  setorder(merged, ld_order)

  if (nrow(merged) != nrow(ld_map)) {
    stop("Unexpected row-count change while merging LD map with GWAS.")
  }

  missing_gwas_rows <- merged[is.na(beta) | is.na(standard_error)]
  if (nrow(missing_gwas_rows) > 0L) {
    stop(
      nrow(missing_gwas_rows),
      " LD variants have no usable GWAS beta/SE."
    )
  }

  merged[, `:=`(
    effect_allele = toupper(effect_allele_gwas),
    other_allele = toupper(other_allele_gwas),
    reference_allele = toupper(reference_allele),
    alternate_allele = toupper(alternate_allele)
  )]

  merged[, allele_orientation := fifelse(
    effect_allele == reference_allele & other_allele == alternate_allele,
    "effect_is_reference",
    fifelse(
      effect_allele == alternate_allele & other_allele == reference_allele,
      "effect_is_alternate",
      "mismatch"
    )
  )]

  mismatch <- merged[allele_orientation == "mismatch"]

  if (nrow(mismatch) > 0L) {
    stop(
      nrow(mismatch),
      " allele mismatches between GWAS and reference panel; examples: ",
      paste(
        head(
          paste0(
            mismatch$variant_id_norm, " GWAS=",
            mismatch$effect_allele, "/", mismatch$other_allele,
            " REF/ALT=", mismatch$reference_allele, "/",
            mismatch$alternate_allele
          ),
          5L
        ),
        collapse = "; "
      )
    )
  }

  if (any(!is.finite(merged$beta)) ||
      any(!is.finite(merged$standard_error)) ||
      any(merged$standard_error <= 0)) {
    stop("GWAS beta/standard_error contains invalid values.")
  }

  merged[, z_original := beta / standard_error]

  # PLINK ref-based LD correlates REF-allele dosage.
  # If GWAS effect allele is ALT, reverse the sign.
  merged[, z_ref := fifelse(
    allele_orientation == "effect_is_reference",
    z_original,
    -z_original
  )]

  if (any(!is.finite(merged$z_ref))) {
    stop("Reference-oriented z-score contains non-finite values.")
  }

  merged
}

validate_variant_order <- function(aligned, ld_vars_file) {
  if (is.na(ld_vars_file) || !file.exists(ld_vars_file)) {
    return(list(
      checked = FALSE,
      status = "vars_file_not_found"
    ))
  }

  vars <- fread(ld_vars_file, header = FALSE)
  if (ncol(vars) < 1L) stop("LD .vars file is empty or malformed.")

  vars <- as.character(vars[[1L]])

  if (length(vars) != nrow(aligned)) {
    stop(
      "LD .vars length (", length(vars),
      ") differs from LD map length (", nrow(aligned), ")."
    )
  }

  expected <- as.character(aligned$reference_variant_id)

  if (!identical(vars, expected)) {
    mismatch_i <- which(vars != expected)
    stop(
      "LD variant order mismatch. First mismatch at row ",
      mismatch_i[1L], ": .vars=", vars[mismatch_i[1L]],
      ", map=", expected[mismatch_i[1L]]
    )
  }

  list(
    checked = TRUE,
    status = "passed"
  )
}

extract_credible_sets <- function(fit, aligned, R) {
  cs_obj <- susieR::susie_get_cs(
    fit,
    Xcorr = R,
    coverage = susie_coverage,
    min_abs_corr = susie_min_abs_corr
  )

  if (is.null(cs_obj$cs) || length(cs_obj$cs) == 0L) {
    return(data.table())
  }

  cs_rows <- vector("list", length(cs_obj$cs))

  for (i in seq_along(cs_obj$cs)) {
    indices <- as.integer(cs_obj$cs[[i]])

    cs_name <- names(cs_obj$cs)[i]
    if (is.null(cs_name) || is.na(cs_name) || cs_name == "") {
      cs_name <- paste0("CS", i)
    }

    purity_min <- NA_real_
    purity_mean <- NA_real_
    purity_median <- NA_real_

    if (!is.null(cs_obj$purity) && nrow(cs_obj$purity) >= i) {
      purity_names <- names(cs_obj$purity)

      min_col <- intersect(
        c("min.abs.corr", "min_abs_corr", "min"),
        purity_names
      )
      mean_col <- intersect(
        c("mean.abs.corr", "mean_abs_corr", "mean"),
        purity_names
      )
      median_col <- intersect(
        c("median.abs.corr", "median_abs_corr", "median"),
        purity_names
      )

      if (length(min_col) > 0L) {
        purity_min <- as.numeric(cs_obj$purity[i, min_col[1L]])
      }
      if (length(mean_col) > 0L) {
        purity_mean <- as.numeric(cs_obj$purity[i, mean_col[1L]])
      }
      if (length(median_col) > 0L) {
        purity_median <- as.numeric(cs_obj$purity[i, median_col[1L]])
      }
    }

    component <- suppressWarnings(as.integer(gsub("[^0-9]", "", cs_name)))
    if (is.na(component) || component < 1L || component > nrow(fit$alpha)) {
      component <- i
    }

    alpha_component <- fit$alpha[component, indices]

    cs_rows[[i]] <- data.table(
      credible_set = cs_name,
      component = component,
      cs_rank = seq_along(indices),
      ld_order = aligned$ld_order[indices],
      variant_id = aligned$variant_id_gwas[indices],
      rsid = if ("rsid_gwas" %in% names(aligned)) aligned$rsid_gwas[indices] else NA_character_,
      reference_variant_id = aligned$reference_variant_id[indices],
      chromosome = aligned$reference_chromosome[indices],
      position = aligned$reference_position[indices],
      reference_allele = aligned$reference_allele[indices],
      alternate_allele = aligned$alternate_allele[indices],
      effect_allele = aligned$effect_allele[indices],
      other_allele = aligned$other_allele[indices],
      beta = aligned$beta[indices],
      standard_error = aligned$standard_error[indices],
      z_original = aligned$z_original[indices],
      z_ref = aligned$z_ref[indices],
      pip = fit$pip[indices],
      alpha = alpha_component,
      cs_size = length(indices),
      purity_min_abs_corr = purity_min,
      purity_mean_abs_corr = purity_mean,
      purity_median_abs_corr = purity_median
    )
  }

  rbindlist(cs_rows, use.names = TRUE, fill = TRUE)
}

make_pip_table <- function(fit, aligned, cs_table) {
  pip <- data.table(
    ld_order = aligned$ld_order,
    variant_id = aligned$variant_id_gwas,
    rsid = if ("rsid_gwas" %in% names(aligned)) aligned$rsid_gwas else NA_character_,
    reference_variant_id = aligned$reference_variant_id,
    chromosome = aligned$reference_chromosome,
    position = aligned$reference_position,
    reference_allele = aligned$reference_allele,
    alternate_allele = aligned$alternate_allele,
    effect_allele = aligned$effect_allele,
    other_allele = aligned$other_allele,
    allele_orientation = aligned$allele_orientation,
    beta = aligned$beta,
    standard_error = aligned$standard_error,
    z_original = aligned$z_original,
    z_ref = aligned$z_ref,
    p_value = if ("p_value" %in% names(aligned)) aligned$p_value else NA_real_,
    effect_allele_frequency = if ("effect_allele_frequency" %in% names(aligned)) {
      aligned$effect_allele_frequency
    } else {
      NA_real_
    },
    pip = fit$pip
  )

  pip[, pip_rank := frank(-pip, ties.method = "min")]

  if (nrow(cs_table) > 0L) {
    cs_membership <- cs_table[
      ,
      .(
        credible_sets = paste(unique(credible_set), collapse = ";"),
        n_credible_sets = uniqueN(credible_set)
      ),
      by = reference_variant_id
    ]

    pip <- merge(
      pip,
      cs_membership,
      by = "reference_variant_id",
      all.x = TRUE,
      sort = FALSE
    )

    setorder(pip, ld_order)
  } else {
    pip[, `:=`(
      credible_sets = NA_character_,
      n_credible_sets = 0L
    )]
  }

  pip
}

empty_manifest_row <- function(clump_id) {
  data.table(
    clump_id = clump_id,
    status = NA_character_,
    reason = NA_character_,
    chromosome = NA_character_,
    region_start = NA_integer_,
    region_end = NA_integer_,
    n_gwas_region = NA_integer_,
    n_ld_variants = NA_integer_,
    n_susie_variants = NA_integer_,
    n_credible_sets = NA_integer_,
    max_pip = NA_real_,
    lead_variant = NA_character_,
    lead_rsid = NA_character_,
    lead_position = NA_integer_,
    lead_p_value = NA_real_,
    n_pip_ge_0_1 = NA_integer_,
    n_pip_ge_0_5 = NA_integer_,
    converged = NA,
    n_iterations = NA_integer_,
    estimate_s = NA_real_,
    ld_order_check = NA_character_,
    ld_psd_check = NA_character_,
    ld_min_eigenvalue = NA_real_,
    max_ld_asymmetry = NA_real_,
    max_ld_diagonal_deviation = NA_real_,
    runtime_seconds = NA_real_,
    susie_version = as.character(packageVersion("susieR")),
    qc_status = NA_character_,
    qc_flags = NA_character_,
    duplicate_group = NA_character_,
    duplicate_flag = FALSE,
    completed_at = NA_character_
  )
}

###############################################################################
# 3. LOAD AND SELECT CLUMPS
###############################################################################

if (!file.exists(ld_manifest_file)) {
  stop("LD manifest not found: ", ld_manifest_file)
}

ld_manifest <- fread(ld_manifest_file)

if (!"clump_id" %in% names(ld_manifest)) {
  stop("clump_ld_manifest.tsv has no clump_id column.")
}

status_col <- intersect(c("status", "ld_status"), names(ld_manifest))
if (length(status_col) == 0L) {
  stop("LD manifest has no status/ld_status column.")
}
status_col <- status_col[1L]

selected <- ld_manifest[
  get(status_col) == "ready"
]

if (nrow(selected) == 0L) {
  stop("No ready clumps found in LD manifest.")
}

if (run_mode == "pilot") {
  selected <- head(selected, pilot_n)
} else if (run_mode != "full") {
  stop("run_mode must be 'pilot' or 'full'.")
}

log_message(
  "Module 10.1 started | mode=", run_mode,
  " | selected clumps=", nrow(selected),
  " | sample size=", gwas_sample_size,
  append = FALSE
)

###############################################################################
# 4. RUN ONE CLUMP
###############################################################################

run_one_clump <- function(clump_id) {
  start_time <- Sys.time()
  result <- empty_manifest_row(clump_id)

  log_message("START ", clump_id)

  tryCatch({
    paths <- resolve_input_files(clump_id)

    missing_paths <- names(paths)[
      names(paths) %in% c("gwas_file", "ld_map_file", "ld_matrix_file") &
        is.na(unlist(paths))
    ]

    if (length(missing_paths) > 0L) {
      stop("Missing required input(s): ", paste(missing_paths, collapse = ", "))
    }

    susie_dir <- file.path(paths$clump_dir, "susie")
    dir.create(susie_dir, recursive = TRUE, showWarnings = FALSE)

    fit_file <- file.path(susie_dir, "susie_fit.rds")
    pip_file <- file.path(susie_dir, "susie_pip.tsv")
    cs_file <- file.path(susie_dir, "susie_credible_sets.tsv")
    summary_file <- file.path(susie_dir, "susie_summary.tsv")
    diagnostics_file <- file.path(susie_dir, "susie_diagnostics.tsv")

    expected_outputs <- c(
      fit_file, pip_file, cs_file, summary_file, diagnostics_file
    )

    if (!overwrite_existing && all(file.exists(expected_outputs))) {
      old_summary <- fread(summary_file)

      result <- copy(old_summary)
      result[, `:=`(
        status = "reused",
        reason = "all_outputs_already_exist",
        runtime_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
        completed_at = timestamp()
      )]

      log_message("REUSED ", clump_id)
      return(result)
    }

    gwas <- fread(paths$gwas_file)
    ld_map <- fread(paths$ld_map_file)

    result[, `:=`(
      n_gwas_region = nrow(gwas),
      n_ld_variants = nrow(ld_map)
    )]

    if (nrow(ld_map) < min_variants_for_susie) {
      stop(
        "Too few LD variants: ", nrow(ld_map),
        " < ", min_variants_for_susie
      )
    }

    aligned <- orient_gwas_to_reference(gwas, ld_map)

    if (!identical(aligned$ld_order, seq_len(nrow(aligned)))) {
      stop("ld_order is not a complete consecutive sequence starting at 1.")
    }

    order_check <- validate_variant_order(aligned, paths$ld_vars_file)

    R_raw <- read_plink_square_binary(
      paths$ld_matrix_file,
      n_variants = nrow(aligned)
    )

    ld_validation <- validate_ld_matrix(R_raw)
    R <- ld_validation$R
    rm(R_raw)
    invisible(gc())

    if (length(aligned$z_ref) != nrow(R)) {
      stop("z-score vector length differs from LD matrix dimension.")
    }

    estimate_s <- NA_real_

    if (run_estimate_s_rss &&
        nrow(R) <= max_variants_for_estimate_s) {
      estimate_s <- tryCatch(
        susieR::estimate_s_rss(
          z = aligned$z_ref,
          R = R,
          n = gwas_sample_size,
          method = "null-mle"
        ),
        error = function(e) {
          log_message(
            "WARNING ", clump_id,
            " estimate_s_rss failed: ", safe_error_text(e)
          )
          NA_real_
        }
      )
    }

    fit <- susieR::susie_rss(
      z = aligned$z_ref,
      R = R,
      n = gwas_sample_size,
      L = min(susie_L, nrow(R)),
      estimate_residual_variance = susie_estimate_residual_variance,
      prior_variance = susie_prior_variance,
      check_prior = TRUE,
      coverage = susie_coverage,
      min_abs_corr = susie_min_abs_corr,
      max_iter = susie_max_iter,
      tol = susie_tol,
      verbose = FALSE
    )

    if (length(fit$pip) != nrow(aligned)) {
      stop("SuSiE PIP vector length differs from input variant count.")
    }

    cs_table <- extract_credible_sets(fit, aligned, R)
    pip_table <- make_pip_table(fit, aligned, cs_table)

    lead <- pip_table[which.max(pip)]

    chromosome_value <- unique(aligned$reference_chromosome)
    chromosome_value <- if (length(chromosome_value) == 1L) {
      chromosome_value
    } else {
      paste(chromosome_value, collapse = ";")
    }

    converged_value <- if (!is.null(fit$converged)) {
      isTRUE(fit$converged)
    } else {
      NA
    }

    n_iterations_value <- if (!is.null(fit$niter)) {
      as.integer(fit$niter)
    } else {
      NA_integer_
    }

    result[, `:=`(
      status = "ready",
      reason = NA_character_,
      chromosome = chromosome_value,
      region_start = min(aligned$reference_position, na.rm = TRUE),
      region_end = max(aligned$reference_position, na.rm = TRUE),
      n_susie_variants = nrow(aligned),
      n_credible_sets = if (nrow(cs_table) == 0L) 0L else uniqueN(cs_table$credible_set),
      max_pip = lead$pip,
      lead_variant = lead$reference_variant_id,
      lead_rsid = lead$rsid,
      lead_position = lead$position,
      lead_p_value = lead$p_value,
      n_pip_ge_0_1 = sum(pip_table$pip >= 0.1, na.rm = TRUE),
      n_pip_ge_0_5 = sum(pip_table$pip >= 0.5, na.rm = TRUE),
      converged = converged_value,
      n_iterations = n_iterations_value,
      estimate_s = estimate_s,
      ld_order_check = order_check$status,
      ld_psd_check = ld_validation$psd_check,
      ld_min_eigenvalue = ld_validation$min_eigenvalue,
      max_ld_asymmetry = ld_validation$max_asymmetry,
      max_ld_diagonal_deviation = ld_validation$max_diagonal_deviation,
      runtime_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
      completed_at = timestamp()
    )]

    qc_flags <- character()
    if (!isTRUE(converged_value)) qc_flags <- c(qc_flags, "nonconverged")
    if (is.na(lead$p_value) || lead$p_value > genomewide_significance_threshold) {
      qc_flags <- c(qc_flags, "max_pip_variant_not_genomewide_significant")
    }
    if (lead$pip < low_max_pip_warning) qc_flags <- c(qc_flags, "low_max_pip")
    if (result$n_credible_sets > max_credible_sets_warning) {
      qc_flags <- c(qc_flags, "many_credible_sets")
    }
    if (!is.na(estimate_s) && estimate_s > estimate_s_warning_threshold) {
      qc_flags <- c(qc_flags, "gwas_ld_inconsistency")
    }
    if (result$n_credible_sets == 0L) qc_flags <- c(qc_flags, "no_credible_set")

    qc_flags <- unique(qc_flags)
    qc_status <- if (length(qc_flags) == 0L) "pass" else "warning"
    final_status <- if (!isTRUE(converged_value)) "warning_nonconverged" else "ready"

    result[, `:=`(
      status = final_status,
      qc_status = qc_status,
      qc_flags = if (length(qc_flags) == 0L) NA_character_ else paste(qc_flags, collapse = ";")
    )]

    diagnostics <- data.table(
      clump_id = clump_id,
      metric = c(
        "n_gwas_region",
        "n_ld_variants",
        "n_susie_variants",
        "max_abs_z",
        "estimate_s_rss",
        "ld_max_asymmetry",
        "ld_max_diagonal_deviation",
        "ld_min_correlation",
        "ld_max_correlation",
        "ld_min_eigenvalue",
        "ld_psd_check",
        "ld_order_check",
        "susie_converged",
        "susie_iterations",
        "susie_L",
        "susie_coverage",
        "susie_min_abs_corr",
        "estimate_residual_variance",
        "sample_size"
      ),
      value = as.character(c(
        nrow(gwas),
        nrow(ld_map),
        nrow(aligned),
        max(abs(aligned$z_ref)),
        estimate_s,
        ld_validation$max_asymmetry,
        ld_validation$max_diagonal_deviation,
        ld_validation$min_correlation,
        ld_validation$max_correlation,
        ld_validation$min_eigenvalue,
        ld_validation$psd_check,
        order_check$status,
        converged_value,
        n_iterations_value,
        min(susie_L, nrow(R)),
        susie_coverage,
        susie_min_abs_corr,
        susie_estimate_residual_variance,
        gwas_sample_size
      ))
    )

    saveRDS(fit, fit_file, compress = TRUE)
    fwrite(pip_table, pip_file, sep = "\t", na = "NA")
    fwrite(cs_table, cs_file, sep = "\t", na = "NA")
    fwrite(result, summary_file, sep = "\t", na = "NA")
    fwrite(diagnostics, diagnostics_file, sep = "\t", na = "NA")

    rm(fit, R, gwas, ld_map, aligned, pip_table, cs_table)
    invisible(gc())

    log_message(
      toupper(result$status), " ", clump_id,
      " | variants=", result$n_susie_variants,
      " | CS=", result$n_credible_sets,
      " | max PIP=", signif(result$max_pip, 4),
      " | runtime=", round(result$runtime_seconds, 1), " sec"
    )

    result
  }, error = function(e) {
    result[, `:=`(
      status = "failed",
      reason = safe_error_text(e),
      runtime_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
      completed_at = timestamp()
    )]

    log_message(
      "FAILED ", clump_id,
      " | ", result$reason
    )

    result
  })
}

###############################################################################
# 5. EXECUTE
###############################################################################

results_list <- vector("list", nrow(selected))

for (i in seq_len(nrow(selected))) {
  clump_id <- selected$clump_id[i]

  log_message(
    "Progress ", i, "/", nrow(selected),
    " | ", clump_id
  )

  results_list[[i]] <- run_one_clump(clump_id)

  # Save progress after every clump.
  current_manifest <- rbindlist(
    results_list[seq_len(i)],
    use.names = TRUE,
    fill = TRUE
  )

  fwrite(
    current_manifest,
    susie_manifest_file,
    sep = "\t",
    na = "NA"
  )
}

susie_manifest <- rbindlist(
  results_list,
  use.names = TRUE,
  fill = TRUE
)

# Detect strongly overlapping or repeated fine-mapping regions. These are
# flagged for review but are not merged automatically.
valid_idx <- which(
  susie_manifest$status %in% c("ready", "warning_nonconverged", "reused") &
  !is.na(susie_manifest$chromosome) &
  !is.na(susie_manifest$region_start) &
  !is.na(susie_manifest$region_end)
)

if (length(valid_idx) > 1L) {
  duplicate_edges <- list()
  edge_n <- 0L

  for (a in seq_len(length(valid_idx) - 1L)) {
    i <- valid_idx[a]
    for (b in (a + 1L):length(valid_idx)) {
      j <- valid_idx[b]
      if (susie_manifest$chromosome[i] != susie_manifest$chromosome[j]) next

      overlap <- max(0L, min(susie_manifest$region_end[i], susie_manifest$region_end[j]) -
                         max(susie_manifest$region_start[i], susie_manifest$region_start[j]) + 1L)
      len_i <- susie_manifest$region_end[i] - susie_manifest$region_start[i] + 1L
      len_j <- susie_manifest$region_end[j] - susie_manifest$region_start[j] + 1L
      reciprocal_overlap <- if (overlap > 0L) min(overlap / len_i, overlap / len_j) else 0
      same_lead <- !is.na(susie_manifest$lead_variant[i]) &&
                   susie_manifest$lead_variant[i] == susie_manifest$lead_variant[j]

      if (reciprocal_overlap >= duplicate_reciprocal_overlap_threshold || same_lead) {
        edge_n <- edge_n + 1L
        duplicate_edges[[edge_n]] <- c(i, j)
      }
    }
  }

  if (length(duplicate_edges) > 0L) {
    parent <- seq_len(nrow(susie_manifest))
    find_root <- function(x) {
      while (parent[x] != x) {
        parent[x] <<- parent[parent[x]]
        x <- parent[x]
      }
      x
    }
    union_root <- function(x, y) {
      rx <- find_root(x); ry <- find_root(y)
      if (rx != ry) parent[ry] <<- rx
    }
    for (edge in duplicate_edges) union_root(edge[1], edge[2])

    roots <- vapply(valid_idx, find_root, integer(1))
    root_counts <- table(roots)
    duplicated_roots <- as.integer(names(root_counts[root_counts > 1L]))
    group_counter <- 0L

    for (root in duplicated_roots) {
      group_counter <- group_counter + 1L
      members <- valid_idx[roots == root]
      group_name <- sprintf("DUPGRP_%04d", group_counter)
      susie_manifest[members, `:=`(
        duplicate_group = group_name,
        duplicate_flag = TRUE,
        qc_status = "warning",
        qc_flags = fifelse(
          is.na(qc_flags) | qc_flags == "",
          "overlapping_or_repeated_locus",
          paste0(qc_flags, ";overlapping_or_repeated_locus")
        )
      )]
    }
  }
}

fwrite(
  susie_manifest,
  susie_manifest_file,
  sep = "\t",
  na = "NA"
)

###############################################################################
# 6. FINAL SUMMARY
###############################################################################

status_counts <- susie_manifest[, .N, by = status][order(status)]

log_message("Module 10.1 finished.")
log_message(
  "Status counts: ",
  paste0(status_counts$status, "=", status_counts$N, collapse = "; ")
)
log_message(
  "Total credible sets: ",
  sum(susie_manifest$n_credible_sets, na.rm = TRUE)
)
log_message(
  "Median max PIP among successful clumps: ",
  signif(
    median(
      susie_manifest[status %in% c("ready", "warning_nonconverged", "reused")]$max_pip,
      na.rm = TRUE
    ),
    4
  )
)

cat("\n")
cat("============================================================\n")
cat("MODULE 10.1 COMPLETE\n")
cat("============================================================\n")
print(status_counts)
cat("\nManifest:", susie_manifest_file, "\n")
cat("Log:     ", log_file, "\n")
cat("============================================================\n")
