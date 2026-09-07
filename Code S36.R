###############################################################################
# PigmentationAtlas
# Module 13.2 – Figure 4: Genome-wide Gene Reclassification Landscape
#
# Purpose:
#   Generate a publication-quality multi-panel Figure 4 from Module 12.6.
#
# Panels:
#   A. Genome-wide distribution of evaluable clumps by PP.H4 and outcome
#   B. Flow from conventional assignment method to reclassification outcome
#   C. Highest-confidence representative reclassification events
#   D. Chromosome-level distribution of reclassified and unchanged clumps
#
# Outputs:
#   figures/Figure4_Genomewide_Reclassification_Landscape.pdf
#   figures/Figure4_Genomewide_Reclassification_Landscape.png
#   figures/Figure4_Genomewide_Reclassification_Landscape.tiff
#
# Author: PigmentationAtlas
###############################################################################

options(stringsAsFactors = FALSE)
options(scipen = 999)

###############################################################################
# 1. Packages
###############################################################################

required_packages <- c(
  "data.table",
  "ggplot2",
  "patchwork",
  "ggalluvial",
  "ggrepel",
  "scales"
)

install_missing_packages <- function(packages) {
  missing <- packages[
    !vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  ]

  if (length(missing) > 0L) {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

install_missing_packages(required_packages)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggalluvial)
  library(ggrepel)
  library(scales)
})

###############################################################################
# 2. Project paths
###############################################################################

project_dir <- normalizePath(
  "C:/Users/User/Desktop/PigmentationAtlas",
  winslash = "/",
  mustWork = TRUE
)

input_dir <- file.path(
  project_dir,
  "results",
  "module12",
  "reclassification"
)

figure_dir <- file.path(project_dir, "figures")
figure_data_dir <- file.path(figure_dir, "figure_data")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_data_dir, recursive = TRUE, showWarnings = FALSE)

reclassification_file <- file.path(
  input_dir,
  "module12_6_clump_gene_reclassification.tsv.gz"
)

###############################################################################
# 3. Helper functions
###############################################################################

message_section <- function(text) {
  cat(
    "\n",
    paste(rep("=", 78), collapse = ""),
    "\n",
    text,
    "\n",
    paste(rep("=", 78), collapse = ""),
    "\n",
    sep = ""
  )
}

assert_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " was not found:\n", path, call. = FALSE)
  }
}

assert_columns <- function(data, columns, label) {
  missing_columns <- setdiff(columns, names(data))

  if (length(missing_columns) > 0L) {
    stop(
      label,
      " is missing required columns:\n",
      paste(missing_columns, collapse = ", "),
      "\n\nAvailable columns:\n",
      paste(names(data), collapse = ", "),
      call. = FALSE
    )
  }
}

standardize_chromosome <- function(x) {
  x <- sub("^chr", "", as.character(x), ignore.case = TRUE)
  factor(x, levels = as.character(1:22), ordered = TRUE)
}

clean_assignment_label <- function(x) {
  labels <- c(
    "OVERLAPPING_PROTEIN_CODING" = "Overlapping coding gene",
    "NEAREST_PROTEIN_CODING" = "Nearest coding gene"
  )

  result <- unname(labels[as.character(x)])
  result[is.na(result)] <- as.character(x)[is.na(result)]
  result
}

clean_outcome_label <- function(x) {
  labels <- c(
    "UNCHANGED_PROTEIN_CODING" = "Unchanged coding gene",
    "CODING_TO_OTHER_PROTEIN_CODING" = "Different coding gene",
    "CODING_TO_PSEUDOGENE" = "Pseudogene",
    "CODING_TO_REGULATORY_RNA" = "Regulatory RNA",
    "CODING_TO_OTHER_OR_UNKNOWN" = "Other / unknown"
  )

  result <- unname(labels[as.character(x)])
  result[is.na(result)] <- as.character(x)[is.na(result)]
  result
}

safe_gene_label <- function(symbol, gene_id) {
  symbol <- trimws(as.character(symbol))
  gene_id <- trimws(as.character(gene_id))

  fifelse(
    !is.na(symbol) & symbol != "",
    symbol,
    fifelse(!is.na(gene_id) & gene_id != "", gene_id, "Unknown")
  )
}

###############################################################################
# 4. Visual settings
###############################################################################

outcome_levels <- c(
  "Unchanged coding gene",
  "Different coding gene",
  "Pseudogene",
  "Regulatory RNA",
  "Other / unknown"
)

outcome_colors <- c(
  "Unchanged coding gene" = "#00BFC4",
  "Different coding gene" = "#F8766D",
  "Pseudogene" = "#B79F00",
  "Regulatory RNA" = "#C77CFF",
  "Other / unknown" = "#7F7F7F"
)

base_theme <- theme_classic(base_size = 11) +
  theme(
    text = element_text(family = "Arial", colour = "black"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    plot.margin = margin(7, 9, 7, 7)
  )

###############################################################################
# 5. Read and validate data
###############################################################################

message_section("FIGURE 4 – INPUT VALIDATION")

assert_file_exists(reclassification_file, "Module 12.6 reclassification table")

reclassification <- fread(reclassification_file, showProgress = TRUE)

required_columns <- c(
  "clump_id",
  "parent_locus",
  "chromosome",
  "index_position",
  "baseline_gene_id",
  "baseline_gene_symbol",
  "baseline_assignment",
  "prioritized_gene_id",
  "prioritized_gene_symbol",
  "priority_class",
  "PP.H4",
  "evaluable",
  "interpretation_change",
  "reclassification_type"
)

assert_columns(reclassification, required_columns, "Reclassification table")

reclassification[, chromosome := standardize_chromosome(chromosome)]
reclassification[, index_position := as.numeric(index_position)]
reclassification[, PP.H4 := as.numeric(PP.H4)]

plot_data <- reclassification[evaluable == TRUE]

if (nrow(plot_data) == 0L) {
  stop("No evaluable clumps were found in the reclassification table.", call. = FALSE)
}

plot_data[, outcome := clean_outcome_label(reclassification_type)]
plot_data[, outcome := factor(outcome, levels = outcome_levels)]
plot_data[, assignment_method := clean_assignment_label(baseline_assignment)]
plot_data[, baseline_label := safe_gene_label(baseline_gene_symbol, baseline_gene_id)]
plot_data[, prioritized_label := safe_gene_label(prioritized_gene_symbol, prioritized_gene_id)]
plot_data[, transition_label := paste0(baseline_label, " 192 ", prioritized_label)]

cat("Evaluable clumps:", nrow(plot_data), "\n")
cat("Reclassified clumps:", sum(plot_data$interpretation_change == TRUE), "\n")
cat("Unchanged clumps:", sum(plot_data$interpretation_change == FALSE), "\n")

###############################################################################
# 6. Panel A – Genome-wide PP.H4 landscape
###############################################################################

message_section("BUILDING PANEL A")

# GRCh38 autosomal chromosome lengths.
chr_lengths <- data.table(
  chromosome = factor(as.character(1:22), levels = as.character(1:22), ordered = TRUE),
  chr_length = c(
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
    159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
    114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
    58617616, 64444167, 46709983, 50818468
  )
)

chr_lengths[, offset := shift(cumsum(chr_length), fill = 0)]
chr_lengths[, midpoint := offset + chr_length / 2]

panel_a_data <- merge(
  plot_data,
  chr_lengths[, .(chromosome, offset, midpoint)],
  by = "chromosome",
  all.x = TRUE,
  sort = FALSE
)

panel_a_data[, genome_position := offset + index_position]

label_data <- panel_a_data[
  interpretation_change == TRUE & !is.na(PP.H4)
][
  order(-PP.H4)
][
  1:min(.N, 8L)
]

panel_a <- ggplot(
  panel_a_data,
  aes(x = genome_position, y = PP.H4, colour = outcome)
) +
  geom_hline(yintercept = c(0.5, 0.8), linetype = c("dashed", "solid"), linewidth = 0.35, colour = "grey45") +
  geom_point(size = 2.4, alpha = 0.85) +
  ggrepel::geom_text_repel(
    data = label_data,
    aes(label = prioritized_label),
    size = 2.7,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    breaks = chr_lengths$midpoint,
    labels = as.character(chr_lengths$chromosome),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = c(0, 1.04),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_colour_manual(values = outcome_colors, drop = TRUE) +
  labs(
    title = "Genome-wide distribution of evaluable colocalization-supported clumps",
    x = "Chromosome",
    y = "Shared-causal posterior probability (PP.H4)",
    colour = "Reclassification outcome"
  ) +
  base_theme +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.text.x = element_text(size = 8)
  )

###############################################################################
# 7. Panel B – Assignment-to-outcome alluvial plot
###############################################################################

message_section("BUILDING PANEL B")

panel_b_data <- plot_data[
  ,
  .N,
  by = .(assignment_method, outcome)
][N > 0]

panel_b_data[, outcome := factor(outcome, levels = outcome_levels)]

panel_b <- ggplot(
  panel_b_data,
  aes(
    axis1 = assignment_method,
    axis2 = outcome,
    y = N
  )
) +
  ggalluvial::geom_alluvium(
    aes(fill = outcome),
    width = 0.18,
    alpha = 0.82,
    knot.pos = 0.45
  ) +
  ggalluvial::geom_stratum(
    width = 0.18,
    fill = "white",
    colour = "black",
    linewidth = 0.35
  ) +
  ggplot2::geom_text(
  stat = "stratum",
  aes(label = after_stat(stratum)),
  size = 3
  ) +
  scale_x_discrete(
    limits = c("Conventional assignment", "Colocalization outcome"),
    expand = c(0.12, 0.12)
  ) +
  scale_fill_manual(values = outcome_colors, drop = TRUE) +
  labs(
    title = "From conventional gene assignment to colocalization-supported outcome",
    x = NULL,
    y = "Number of evaluable clumps",
    fill = "Outcome"
  ) +
  base_theme +
  theme(
    legend.position = "none",
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(face = "bold", size = 9)
  )

###############################################################################
# 8. Panel C – Representative high-confidence reclassification events
###############################################################################

message_section("BUILDING PANEL C")

n_top <- 15L

panel_c_data <- plot_data[
  interpretation_change == TRUE & !is.na(PP.H4)
][
  order(-PP.H4, priority_class, clump_id)
][
  1:min(.N, n_top)
]

panel_c_data[, display_label := paste0(transition_label, "  [", clump_id, "]")]
panel_c_data[, display_label := factor(display_label, levels = rev(display_label))]

panel_c <- ggplot(
  panel_c_data,
  aes(x = PP.H4, y = display_label, colour = outcome)
) +
  geom_segment(
    aes(x = 0, xend = PP.H4, yend = display_label),
    colour = "grey78",
    linewidth = 0.55
  ) +
  geom_point(size = 3.2) +
  geom_vline(xintercept = c(0.5, 0.8), linetype = c("dashed", "solid"), linewidth = 0.35, colour = "grey45") +
  scale_x_continuous(
    limits = c(0, 1.02),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_colour_manual(values = outcome_colors, drop = TRUE) +
  labs(
    title = "Highest-confidence representative reclassification events",
    subtitle = "Conventional gene 192 colocalization-prioritized gene",
    x = "Shared-causal posterior probability (PP.H4)",
    y = NULL,
    colour = "Outcome"
  ) +
  base_theme +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 7.4)
  )

###############################################################################
# 9. Panel D – Chromosome-level distribution
###############################################################################

message_section("BUILDING PANEL D")

panel_d_data <- plot_data[
  ,
  .N,
  by = .(chromosome, outcome)
]

panel_d_data[, chromosome := factor(chromosome, levels = as.character(1:22), ordered = TRUE)]
panel_d_data[, outcome := factor(outcome, levels = outcome_levels)]

panel_d <- ggplot(
  panel_d_data,
  aes(x = chromosome, y = N, fill = outcome)
) +
  geom_col(width = 0.78) +
  scale_fill_manual(values = outcome_colors, drop = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Chromosome-level distribution of evaluable clumps",
    x = "Chromosome",
    y = "Number of evaluable clumps",
    fill = "Reclassification outcome"
  ) +
  base_theme +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8)
  )

###############################################################################
# 10. Assemble Figure 4
###############################################################################

message_section("ASSEMBLING FIGURE 4")

n_evaluable <- nrow(plot_data)
n_reclassified <- sum(plot_data$interpretation_change == TRUE, na.rm = TRUE)
percent_reclassified <- 100 * n_reclassified / n_evaluable

figure_title <- "Genome-wide landscape of colocalization-driven gene reclassification"
figure_subtitle <- sprintf(
  "Among %d evaluable pigmentation GWAS clumps, %d (%.2f%%) were assigned to a different candidate gene.",
  n_evaluable,
  n_reclassified,
  percent_reclassified
)

figure4 <- (
  panel_a +
    plot_annotation(tag_levels = "A")
) / (
  panel_b + panel_c + panel_d +
    plot_layout(widths = c(1.00, 1.30, 1.00))
) +
  plot_layout(heights = c(0.86, 1.14), guides = "collect") +
  plot_annotation(
    title = figure_title,
    subtitle = figure_subtitle,
    tag_levels = "A",
    theme = theme(
      text = element_text(family = "Arial", colour = "black"),
      plot.title = element_text(size = 17, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 10.5, hjust = 0),
      plot.tag = element_text(size = 14, face = "bold")
    )
  ) &
  theme(legend.position = "bottom")

###############################################################################
# 11. Write outputs
###############################################################################

message_section("WRITING FIGURE 4 OUTPUTS")

pdf_file <- file.path(figure_dir, "Figure4_Genomewide_Reclassification_Landscape.pdf")
png_file <- file.path(figure_dir, "Figure4_Genomewide_Reclassification_Landscape.png")
tiff_file <- file.path(figure_dir, "Figure4_Genomewide_Reclassification_Landscape.tiff")

ggsave(
  pdf_file,
  figure4,
  width = 15.5,
  height = 11.3,
  units = "in",
  device = cairo_pdf
)

ggsave(
  png_file,
  figure4,
  width = 15.5,
  height = 11.3,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  tiff_file,
  figure4,
  width = 15.5,
  height = 11.3,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

fwrite(
  panel_a_data,
  file.path(figure_data_dir, "Figure4_panelA_genomewide_data.tsv"),
  sep = "\t"
)

fwrite(
  panel_b_data,
  file.path(figure_data_dir, "Figure4_panelB_alluvial_data.tsv"),
  sep = "\t"
)

fwrite(
  panel_c_data,
  file.path(figure_data_dir, "Figure4_panelC_top_reclassification_events.tsv"),
  sep = "\t"
)

fwrite(
  panel_d_data,
  file.path(figure_data_dir, "Figure4_panelD_chromosome_counts.tsv"),
  sep = "\t"
)

status <- data.table(
  module = "Module 13.2",
  figure = "Figure 4",
  analysis = "Genome-wide Gene Reclassification Landscape",
  status = "SUCCESS",
  evaluable_clumps = n_evaluable,
  reclassified_clumps = n_reclassified,
  reclassified_percent = percent_reclassified,
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

fwrite(
  status,
  file.path(figure_data_dir, "Figure4_status.tsv"),
  sep = "\t"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(figure_data_dir, "Figure4_sessionInfo.txt")
)

###############################################################################
# 12. Console report
###############################################################################

message_section("FIGURE 4 COMPLETED")

cat("Evaluable clumps:    ", n_evaluable, "\n")
cat("Reclassified clumps: ", n_reclassified, "\n")
cat("Reclassified percent:", sprintf("%.2f%%", percent_reclassified), "\n\n")
cat("PDF:\n", pdf_file, "\n\n")
cat("PNG:\n", png_file, "\n\n")
cat("TIFF:\n", tiff_file, "\n")
