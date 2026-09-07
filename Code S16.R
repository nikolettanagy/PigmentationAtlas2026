###############################################################################
# PigmentationAtlas
# Module 10.2 — Automated regional SuSiE fine-mapping figures
#
# Reads Module 10.1 outputs and creates one publication-ready PDF and PNG
# per clump with three panels:
#   A. GWAS regional association (-log10 P)
#   B. SuSiE posterior inclusion probabilities
#   C. Credible-set membership
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
})

###############################################################################
# 1. USER CONFIGURATION
###############################################################################

project_dir <- "C:/Users/User/Desktop/PigmentationAtlas"
setwd(project_dir)

run_mode <- "pilot"       # "pilot" or "full"
pilot_n <- 5L

overwrite_existing <- TRUE

susie_manifest_file <- "results/susie_manifest.tsv"
figure_manifest_file <- "results/susie_figure_manifest.tsv"
log_dir <- "logs"
log_file <- file.path(log_dir, "module10_2_figures.log")

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# Figure dimensions
pdf_width <- 9.0
pdf_height <- 8.0
png_width_inches <- 9.0
png_height_inches <- 8.0
png_resolution <- 300L

# Plotting thresholds
pip_reference_lines <- c(0.1, 0.5, 0.9)
genomewide_significance_threshold <- 5e-8

###############################################################################
# 2. HELPERS
###############################################################################

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_message <- function(..., append = TRUE) {
  msg <- paste0("[", timestamp(), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = append)
}

safe_error_text <- function(e) paste(class(e)[1], conditionMessage(e), sep = ": ")

first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0L) NA_character_ else hit[1L]
}

resolve_files <- function(clump_id) {
  clump_dir <- file.path("results", "clumps", clump_id)
  susie_dir <- file.path(clump_dir, "susie")

  list(
    clump_dir = clump_dir,
    susie_dir = susie_dir,
    pip_file = first_existing_file(c(
      file.path(susie_dir, "susie_pip.tsv")
    )),
    cs_file = first_existing_file(c(
      file.path(susie_dir, "susie_credible_sets.tsv")
    )),
    summary_file = first_existing_file(c(
      file.path(susie_dir, "susie_summary.tsv")
    ))
  )
}

format_bp_axis <- function(x) {
  format(round(x / 1e6, 3), trim = TRUE, scientific = FALSE)
}

safe_neglog10 <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  p[p <= 0] <- .Machine$double.xmin
  -log10(p)
}

credible_set_palette <- function(n) {
  if (n <= 0L) return(character())
  grDevices::hcl.colors(n, palette = "Dark 3")
}

open_pdf <- function(path) {
  grDevices::pdf(
    path,
    width = pdf_width,
    height = pdf_height,
    useDingbats = FALSE,
    onefile = TRUE
  )
}

open_png <- function(path) {
  grDevices::png(
    path,
    width = png_width_inches,
    height = png_height_inches,
    units = "in",
    res = png_resolution,
    type = "cairo-png"
  )
}

plot_locus <- function(pip, cs, clump_id, qc_status, qc_flags) {
  pip <- copy(pip)
  setorder(pip, position)

  if (!"p_value" %in% names(pip)) pip[, p_value := NA_real_]
  if (!"pip" %in% names(pip)) stop("PIP table has no pip column.")
  if (!"position" %in% names(pip)) stop("PIP table has no position column.")

  pip[, neglog10p := safe_neglog10(p_value)]
  x <- pip$position / 1e6

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  layout(matrix(1:3, ncol = 1), heights = c(1.25, 1.0, 0.85))
  par(
    mar = c(0.8, 5.0, 2.6, 1.2),
    oma = c(4.0, 0.5, 3.0, 0.5),
    mgp = c(2.8, 0.8, 0),
    tcl = -0.3,
    las = 1,
    family = "sans"
  )

  # Panel A: GWAS association
  finite_gwas <- is.finite(pip$neglog10p)
  ymax_gwas <- if (any(finite_gwas)) max(pip$neglog10p[finite_gwas], 8, na.rm = TRUE) else 8

  plot(
    x, pip$neglog10p,
    pch = 16,
    cex = 0.55,
    xlab = "",
    ylab = expression(-log[10](italic(P))),
    xaxt = "n",
    ylim = c(0, ymax_gwas * 1.06)
  )
  abline(h = -log10(genomewide_significance_threshold), lty = 2)

  if (any(finite_gwas)) {
    lead_gwas_i <- which.max(pip$neglog10p)
    points(x[lead_gwas_i], pip$neglog10p[lead_gwas_i], pch = 21, cex = 1.05, lwd = 1.2)
    lead_label <- if (!is.na(pip$rsid[lead_gwas_i]) && pip$rsid[lead_gwas_i] != "") {
      pip$rsid[lead_gwas_i]
    } else {
      pip$reference_variant_id[lead_gwas_i]
    }
    text(
      x[lead_gwas_i], pip$neglog10p[lead_gwas_i],
      labels = lead_label,
      pos = 3,
      cex = 0.72,
      xpd = NA
    )
  }
  mtext("A  GWAS regional association", side = 3, adj = 0, line = 0.5, font = 2)

  # Panel B: PIP
  par(mar = c(0.8, 5.0, 2.4, 1.2))
  plot(
    x, pip$pip,
    type = "h",
    lwd = 1.2,
    xlab = "",
    ylab = "PIP",
    xaxt = "n",
    ylim = c(0, 1.04)
  )
  points(x, pip$pip, pch = 16, cex = 0.48)
  abline(h = pip_reference_lines, lty = c(3, 2, 3))

  lead_pip_i <- which.max(pip$pip)
  points(x[lead_pip_i], pip$pip[lead_pip_i], pch = 21, cex = 1.1, lwd = 1.2)
  pip_label <- if (!is.na(pip$rsid[lead_pip_i]) && pip$rsid[lead_pip_i] != "") {
    pip$rsid[lead_pip_i]
  } else {
    pip$reference_variant_id[lead_pip_i]
  }
  text(
    x[lead_pip_i], pip$pip[lead_pip_i],
    labels = paste0(pip_label, "  PIP=", formatC(pip$pip[lead_pip_i], digits = 3, format = "f")),
    pos = 3,
    cex = 0.72,
    xpd = NA
  )
  mtext("B  SuSiE posterior inclusion probabilities", side = 3, adj = 0, line = 0.4, font = 2)

  # Panel C: credible sets
  par(mar = c(3.5, 5.0, 2.4, 1.2))

  if (nrow(cs) == 0L) {
    plot(
      range(x), c(0, 1),
      type = "n",
      xlab = "Genomic position (Mb)",
      ylab = "Credible set",
      yaxt = "n"
    )
    text(mean(range(x)), 0.5, "No credible set passed the purity filter")
  } else {
    cs <- copy(cs)
    cs[, credible_set := as.character(credible_set)]
    cs_levels <- unique(cs$credible_set)
    cs[, cs_index := match(credible_set, cs_levels)]
    cs_colors <- credible_set_palette(length(cs_levels))

    plot(
      range(x), c(0.5, length(cs_levels) + 0.5),
      type = "n",
      xlab = "Genomic position (Mb)",
      ylab = "Credible set",
      yaxt = "n"
    )
    axis(2, at = seq_along(cs_levels), labels = cs_levels, las = 1)

    for (k in seq_along(cs_levels)) {
      d <- cs[credible_set == cs_levels[k]]
      points(
        d$position / 1e6,
        rep(k, nrow(d)),
        pch = 16,
        cex = pmax(0.75, 1.7 * sqrt(d$pip)),
        col = cs_colors[k]
      )
      if (nrow(d) > 1L) {
        segments(min(d$position) / 1e6, k, max(d$position) / 1e6, k, col = cs_colors[k])
      }
    }
  }

  mtext("C  95% credible-set membership", side = 3, adj = 0, line = 0.4, font = 2)

  chromosome <- unique(pip$chromosome)
  chromosome <- paste(chromosome[!is.na(chromosome)], collapse = ";")
  title_text <- paste0(
    clump_id,
    " | chr", chromosome, ":",
    format(min(pip$position, na.rm = TRUE), big.mark = ",", scientific = FALSE),
    "–",
    format(max(pip$position, na.rm = TRUE), big.mark = ",", scientific = FALSE)
  )

  qc_text <- paste0(
    "QC: ", qc_status,
    if (!is.na(qc_flags) && qc_flags != "") paste0(" — ", qc_flags) else ""
  )

  mtext(title_text, side = 3, outer = TRUE, line = 1.2, font = 2, cex = 1.05)
  mtext(qc_text, side = 3, outer = TRUE, line = 0.1, cex = 0.72)
}

###############################################################################
# 3. LOAD MANIFEST
###############################################################################

if (!file.exists(susie_manifest_file)) {
  stop("SuSiE manifest not found: ", susie_manifest_file)
}

manifest <- fread(susie_manifest_file)

required_manifest <- c("clump_id", "status")
missing_manifest <- setdiff(required_manifest, names(manifest))
if (length(missing_manifest) > 0L) {
  stop("Missing manifest columns: ", paste(missing_manifest, collapse = ", "))
}

selected <- manifest[status %in% c("ready", "warning_nonconverged", "reused")]

if (run_mode == "pilot") {
  selected <- head(selected, pilot_n)
} else if (run_mode != "full") {
  stop("run_mode must be 'pilot' or 'full'.")
}

if (nrow(selected) == 0L) stop("No eligible clumps found for plotting.")

log_message(
  "Module 10.2 started | mode=", run_mode,
  " | selected clumps=", nrow(selected),
  append = FALSE
)

###############################################################################
# 4. CREATE FIGURES
###############################################################################

results <- vector("list", nrow(selected))

for (i in seq_len(nrow(selected))) {
  clump_id <- selected$clump_id[i]
  started <- Sys.time()

  log_message("START ", clump_id, " | ", i, "/", nrow(selected))

  results[[i]] <- tryCatch({
    paths <- resolve_files(clump_id)

    if (is.na(paths$pip_file)) stop("susie_pip.tsv not found.")
    if (is.na(paths$cs_file)) stop("susie_credible_sets.tsv not found.")

    pip <- fread(paths$pip_file)
    cs <- fread(paths$cs_file)

    figure_dir <- file.path(paths$clump_dir, "figures")
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

    pdf_file <- file.path(figure_dir, paste0(clump_id, "_susie_finemapping.pdf"))
    png_file <- file.path(figure_dir, paste0(clump_id, "_susie_finemapping.png"))

    if (!overwrite_existing && file.exists(pdf_file) && file.exists(png_file)) {
      log_message("REUSED ", clump_id)
      data.table(
        clump_id = clump_id,
        status = "reused",
        reason = NA_character_,
        pdf_file = pdf_file,
        png_file = png_file,
        n_variants = nrow(pip),
        n_credible_sets = if (nrow(cs) == 0L) 0L else uniqueN(cs$credible_set),
        runtime_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
        completed_at = timestamp()
      )
    } else {
      qc_status <- if ("qc_status" %in% names(selected)) selected$qc_status[i] else selected$status[i]
      qc_flags <- if ("qc_flags" %in% names(selected)) selected$qc_flags[i] else NA_character_

      open_pdf(pdf_file)
      plot_locus(pip, cs, clump_id, qc_status, qc_flags)
      dev.off()

      open_png(png_file)
      plot_locus(pip, cs, clump_id, qc_status, qc_flags)
      dev.off()

      log_message(
        "READY ", clump_id,
        " | variants=", nrow(pip),
        " | CS=", if (nrow(cs) == 0L) 0L else uniqueN(cs$credible_set)
      )

      data.table(
        clump_id = clump_id,
        status = "ready",
        reason = NA_character_,
        pdf_file = pdf_file,
        png_file = png_file,
        n_variants = nrow(pip),
        n_credible_sets = if (nrow(cs) == 0L) 0L else uniqueN(cs$credible_set),
        runtime_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
        completed_at = timestamp()
      )
    }
  }, error = function(e) {
    log_message("FAILED ", clump_id, " | ", safe_error_text(e))
    data.table(
      clump_id = clump_id,
      status = "failed",
      reason = safe_error_text(e),
      pdf_file = NA_character_,
      png_file = NA_character_,
      n_variants = NA_integer_,
      n_credible_sets = NA_integer_,
      runtime_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
      completed_at = timestamp()
    )
  })

  current <- rbindlist(results[seq_len(i)], use.names = TRUE, fill = TRUE)
  fwrite(current, figure_manifest_file, sep = "\t", na = "NA")
}

figure_manifest <- rbindlist(results, use.names = TRUE, fill = TRUE)
fwrite(figure_manifest, figure_manifest_file, sep = "\t", na = "NA")

status_counts <- figure_manifest[, .N, by = status][order(status)]

log_message("Module 10.2 finished.")
log_message(
  "Status counts: ",
  paste0(status_counts$status, "=", status_counts$N, collapse = "; ")
)

cat("\n")
cat("============================================================\n")
cat("MODULE 10.2 COMPLETE\n")
cat("============================================================\n")
print(status_counts)
cat("\nManifest:", figure_manifest_file, "\n")
cat("Log:     ", log_file, "\n")
cat("============================================================\n")
