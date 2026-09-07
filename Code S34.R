#!/usr/bin/env Rscript

# PigmentationAtlas
# Module 13 / Figure 2: Genome-wide colocalization summary
#
# Usage from the PigmentationAtlas project root:
#   Rscript scripts/module13_figure2_colocalization_summary.R
#
# Optional command-line arguments:
#   1. input catalogue
#   2. output directory
# Example:
#   Rscript scripts/module13_figure2_colocalization_summary.R \
#     results/module12_5B_candidate_gene_prioritization_GENCODEv50.tsv.gz \
#     figures/main

options(stringsAsFactors = FALSE)

required_packages <- c("data.table", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing R package(s): ", paste(missing_packages, collapse = ", "),
    "\nInstall them with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "), "))",
    call. = FALSE
  )
}

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)

args <- commandArgs(trailingOnly = TRUE)

input_candidates <- c(
  if (length(args) >= 1L) args[[1L]] else character(),
  "results/module12_5B_candidate_gene_prioritization_GENCODEv50.tsv.gz",
  "results/module12_5B/module12_5B_candidate_gene_prioritization_GENCODEv50.tsv.gz",
  "module12_5B_candidate_gene_prioritization_GENCODEv50.tsv.gz"
)
input_candidates <- unique(input_candidates[nzchar(input_candidates)])
input_file <- input_candidates[file.exists(input_candidates)][1L]

if (is.na(input_file)) {
  stop(
    "Input catalogue not found. Checked:\n  ",
    paste(input_candidates, collapse = "\n  "),
    call. = FALSE
  )
}

output_dir <- if (length(args) >= 2L) args[[2L]] else "figures/main"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading: ", input_file)
dt <- fread(input_file)

required_columns <- c(
  "priority_class", "PP.H3", "PP.H4", "gene_id_stable",
  "gene_symbol", "candidate_class", "tissue"
)
missing_columns <- setdiff(required_columns, names(dt))
if (length(missing_columns) > 0L) {
  stop("Missing required column(s): ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

# -----------------------------------------------------------------------------
# Data preparation
# -----------------------------------------------------------------------------
priority_levels <- c(
  "STRONG_H4", "MODERATE_H4", "SUGGESTIVE_H4",
  "H3_DOMINANT", "INCONCLUSIVE"
)
priority_labels <- c(
  STRONG_H4 = "Strong H4",
  MODERATE_H4 = "Moderate H4",
  SUGGESTIVE_H4 = "Suggestive H4",
  H3_DOMINANT = "H3 dominant",
  INCONCLUSIVE = "Inconclusive"
)

priority_colours <- c(
  "Strong H4" = "#7F0000",
  "Moderate H4" = "#D7301F",
  "Suggestive H4" = "#FC8D59",
  "H3 dominant" = "#4575B4",
  "Inconclusive" = "#BDBDBD"
)

dt[, priority_class := factor(priority_class, levels = priority_levels)]
dt[, priority_label := factor(
  priority_labels[as.character(priority_class)],
  levels = unname(priority_labels[priority_levels])
)]

dt[, PP.H3 := as.numeric(PP.H3)]
dt[, PP.H4 := as.numeric(PP.H4)]
dt <- dt[is.finite(PP.H3) & is.finite(PP.H4)]

n_total <- nrow(dt)
if (n_total == 0L) stop("No valid PP.H3/PP.H4 rows remained after filtering.", call. = FALSE)

category_summary <- dt[, .(
  n_units = .N,
  n_genes = uniqueN(gene_id_stable)
), by = priority_label]
category_summary[, percentage := 100 * n_units / sum(n_units)]
category_summary[, count_label := paste0(
  format(n_units, big.mark = ",", scientific = FALSE),
  "  (", sprintf("%.2f", percentage), "%)"
)]

# Preserve all five classes even if a future dataset has a zero-count category.
category_template <- data.table(
  priority_label = factor(unname(priority_labels[priority_levels]),
                          levels = unname(priority_labels[priority_levels]))
)
category_summary <- merge(category_template, category_summary,
                          by = "priority_label", all.x = TRUE, sort = FALSE)
category_summary[is.na(n_units), `:=`(n_units = 0L, n_genes = 0L,
                                     percentage = 0,
                                     count_label = "0  (0.00%)")]

# -----------------------------------------------------------------------------
# Shared publication theme
# -----------------------------------------------------------------------------
theme_pa <- theme_classic(base_size = 10.5, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 11.5, hjust = 0),
    plot.subtitle = element_text(size = 9.5, colour = "grey30", margin = margin(b = 6)),
    axis.title = element_text(face = "bold", size = 9.5),
    axis.text = element_text(size = 8.5, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8.5),
    legend.key.height = unit(4, "mm"),
    plot.margin = margin(7, 8, 7, 8)
  )

# -----------------------------------------------------------------------------
# Panel A: categorical selectivity
# -----------------------------------------------------------------------------
panel_a <- ggplot(category_summary,
                  aes(x = priority_label, y = n_units, fill = priority_label)) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.25) +
  geom_text(
    aes(label = count_label),
    hjust = -0.08, size = 3.0, fontface = "bold"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = priority_colours, drop = FALSE) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.22))
  ) +
  labs(
    title = "A  Genome-wide colocalization classification",
    subtitle = paste0(format(n_total, big.mark = ","), " locus–gene–tissue units"),
    x = NULL,
    y = "Number of colocalization units"
  ) +
  guides(fill = "none") +
  theme_pa +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

# -----------------------------------------------------------------------------
# Panel B: PP.H4 distribution on a scale that shows the low-probability mass
# -----------------------------------------------------------------------------
# A square-root y transformation preserves zero and reveals the long tail.
panel_b <- ggplot(dt, aes(x = PP.H4, fill = priority_label)) +
  geom_histogram(
    binwidth = 0.02, boundary = 0,
    colour = "white", linewidth = 0.12
  ) +
  geom_vline(xintercept = c(0.20, 0.50, 0.80),
             linetype = c("dotted", "dashed", "solid"),
             linewidth = 0.45, colour = "grey20") +
  annotate("text", x = c(0.20, 0.50, 0.80), y = Inf,
           label = c("0.20", "0.50", "0.80"),
           vjust = 1.4, hjust = -0.15, size = 2.7) +
  scale_fill_manual(values = priority_colours, drop = FALSE) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  scale_y_sqrt(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "B  Distribution of shared-signal posterior probability",
    subtitle = "Vertical lines denote suggestive, moderate and strong PP.H4 thresholds",
    x = "Posterior probability of a shared causal signal (PP.H4)",
    y = "Number of units (square-root scale)"
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_pa +
  theme(legend.position = "bottom")

# -----------------------------------------------------------------------------
# Panel C: direct competition between distinct and shared causal hypotheses
# -----------------------------------------------------------------------------
plot_dt <- copy(dt)
plot_dt[, point_alpha := fifelse(priority_label == "Inconclusive", 0.16, 0.70)]
plot_dt[, point_size := fifelse(priority_label %in% c("Strong H4", "Moderate H4"), 1.55, 0.85)]

panel_c <- ggplot(plot_dt, aes(x = PP.H3, y = PP.H4, colour = priority_label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.45, colour = "grey35") +
  geom_point(aes(alpha = point_alpha, size = point_size), stroke = 0) +
  geom_hline(yintercept = c(0.20, 0.50, 0.80),
             linetype = c("dotted", "dashed", "solid"),
             linewidth = 0.35, colour = "grey45") +
  scale_colour_manual(values = priority_colours, drop = FALSE) +
  scale_alpha_identity() +
  scale_size_identity() +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  coord_equal() +
  labs(
    title = "C  Shared versus distinct causal architecture",
    subtitle = "Points above the diagonal favour H4 over H3",
    x = "Distinct causal signals (PP.H3)",
    y = "Shared causal signal (PP.H4)"
  ) +
  guides(colour = "none") +
  theme_pa

# -----------------------------------------------------------------------------
# Assemble and export
# -----------------------------------------------------------------------------
figure_2 <- panel_a / (panel_b | panel_c) +
  plot_layout(heights = c(0.82, 1.18), widths = c(1.35, 1)) +
  plot_annotation(
    title = "PigmentationAtlas genome-wide colocalization landscape",
    subtitle = paste0(
      "Only ", category_summary[priority_label == "Strong H4", n_units],
      " of ", format(n_total, big.mark = ","),
      " tested units reached Strong H4 support"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, family = "Arial"),
      plot.subtitle = element_text(size = 10.5, colour = "grey25", family = "Arial",
                                   margin = margin(b = 8)),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

pdf_file <- file.path(output_dir, "Figure_2_PigmentationAtlas_colocalization_summary.pdf")
png_file <- file.path(output_dir, "Figure_2_PigmentationAtlas_colocalization_summary.png")
tiff_file <- file.path(output_dir, "Figure_2_PigmentationAtlas_colocalization_summary.tiff")
summary_file <- file.path(output_dir, "Figure_2_colocalization_category_counts.tsv")

ggsave(pdf_file, figure_2, width = 180, height = 205, units = "mm", device = cairo_pdf)
ggsave(png_file, figure_2, width = 180, height = 205, units = "mm", dpi = 600, bg = "white")
ggsave(tiff_file, figure_2, width = 180, height = 205, units = "mm",
       dpi = 600, compression = "lzw", bg = "white")
fwrite(category_summary, summary_file, sep = "\t")

message("\nFigure 2 completed successfully.")
message("PDF:  ", pdf_file)
message("PNG:  ", png_file)
message("TIFF: ", tiff_file)
message("Data: ", summary_file)
message("Strong H4 units: ", category_summary[priority_label == "Strong H4", n_units])
message("Total units: ", n_total)
