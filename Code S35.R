###############################################################################
# PigmentationAtlas
# Module 13.1 – Figure 3: Colocalization-driven Gene Reclassification
#
# Purpose:
#   Generate a publication-quality, fully reproducible multi-panel Figure 3
#   directly from Module 12.6 outputs.
#
# Panels:
#   A. Analysis workflow
#   B. Main reclassification result
#   C. Reclassification categories
#   D. Reclassification by colocalization priority class
#
# Outputs:
#   figures/Figure3_Gene_Reclassification.pdf
#   figures/Figure3_Gene_Reclassification.png
#   figures/Figure3_Gene_Reclassification.tiff
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
  "scales"
)

install_missing_packages <- function(packages) {
  missing <- packages[
    !vapply(
      packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(missing) > 0L) {
    install.packages(
      missing,
      repos = "https://cloud.r-project.org"
    )
  }
}

install_missing_packages(required_packages)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
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

results_dir <- file.path(
  project_dir,
  "results",
  "module12",
  "reclassification"
)

figure_dir <- file.path(
  project_dir,
  "figures"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

reclassification_file <- file.path(
  results_dir,
  "module12_6_clump_gene_reclassification.tsv.gz"
)

type_summary_file <- file.path(
  results_dir,
  "module12_6_reclassification_type_summary.tsv"
)

priority_summary_file <- file.path(
  results_dir,
  "module12_6_reclassification_by_priority_class.tsv"
)

overall_summary_file <- file.path(
  results_dir,
  "module12_6_reclassification_overall_summary.tsv"
)

status_file <- file.path(
  results_dir,
  "module12_6_status.tsv"
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
    stop(
      label,
      " was not found:\n",
      path,
      call. = FALSE
    )
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

clean_priority_label <- function(x) {
  x <- toupper(gsub("[ -]", "_", trimws(as.character(x))))

  labels <- c(
    "STRONG_H4" = "Strong H4",
    "MODERATE_H4" = "Moderate H4",
    "SUGGESTIVE_H4" = "Suggestive H4"
  )

  result <- unname(labels[x])
  result[is.na(result)] <- x[is.na(result)]
  result
}

clean_type_label <- function(x) {
  labels <- c(
    "UNCHANGED_PROTEIN_CODING" = "Unchanged protein-coding gene",
    "CODING_TO_OTHER_PROTEIN_CODING" = "Different protein-coding gene",
    "CODING_TO_PSEUDOGENE" = "Pseudogene",
    "CODING_TO_REGULATORY_RNA" = "Regulatory RNA",
    "CODING_TO_OTHER_OR_UNKNOWN" = "Other or unknown"
  )

  result <- unname(labels[as.character(x)])
  result[is.na(result)] <- as.character(x)[is.na(result)]
  result
}

###############################################################################
# 4. Input validation and reading
###############################################################################

message_section("FIGURE 3 – INPUT VALIDATION")

assert_file_exists(reclassification_file, "Reclassification table")
assert_file_exists(type_summary_file, "Reclassification type summary")
assert_file_exists(priority_summary_file, "Priority-class summary")
assert_file_exists(overall_summary_file, "Overall summary")
assert_file_exists(status_file, "Module 12.6 status")

reclassification <- fread(reclassification_file)
type_summary <- fread(type_summary_file)
priority_summary <- fread(priority_summary_file)
overall_summary <- fread(overall_summary_file)
status <- fread(status_file)

assert_columns(
  reclassification,
  c(
    "clump_id",
    "evaluable",
    "interpretation_change",
    "reclassification_type",
    "priority_class"
  ),
  "Reclassification table"
)

assert_columns(
  type_summary,
  c("reclassification_type", "N"),
  "Reclassification type summary"
)

assert_columns(
  priority_summary,
  c(
    "priority_class",
    "n_clumps",
    "n_reclassified",
    "percent_reclassified"
  ),
  "Priority-class summary"
)

assert_columns(
  overall_summary,
  c("metric", "value"),
  "Overall summary"
)

###############################################################################
# 5. Core numbers and QC
###############################################################################

message_section("CALCULATING FIGURE STATISTICS")

n_total <- nrow(reclassification)
n_evaluable <- sum(reclassification$evaluable, na.rm = TRUE)
n_reclassified <- sum(
  reclassification$interpretation_change == TRUE,
  na.rm = TRUE
)
n_unchanged <- sum(
  reclassification$interpretation_change == FALSE,
  na.rm = TRUE
)

percent_reclassified <- 100 * n_reclassified / n_evaluable
percent_unchanged <- 100 * n_unchanged / n_evaluable

if (n_evaluable != n_reclassified + n_unchanged) {
  stop(
    "QC failure: evaluable clumps do not equal reclassified + unchanged.",
    call. = FALSE
  )
}

if (anyDuplicated(reclassification$clump_id)) {
  stop(
    "QC failure: duplicated clump_id values in reclassification table.",
    call. = FALSE
  )
}

cat("Total PLINK clumps:          ", n_total, "\n")
cat("Evaluable clumps:            ", n_evaluable, "\n")
cat("Reclassified clumps:         ", n_reclassified, "\n")
cat("Unchanged clumps:            ", n_unchanged, "\n")
cat("Percent reclassified:        ", sprintf("%.2f%%", percent_reclassified), "\n")

###############################################################################
# 6. Shared theme
###############################################################################

base_family <- "sans"

figure_theme <- theme_classic(base_size = 9, base_family = base_family) +
  theme(
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0
    ),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    plot.margin = margin(6, 8, 6, 8)
  )

###############################################################################
# 7. Panel A – Workflow
###############################################################################

workflow_boxes <- data.table(
  x = c(1, 2.5, 4, 5.5),
  y = 1,
  label = c(
    "GWAS lead\nvariant",
    "Conventional assignment\n(overlapping or nearest\nprotein-coding gene)",
    "Colocalization\n(H4-supported\ncandidate gene)",
    "Gene\nreclassification"
  )
)

workflow_arrows <- data.table(
  x = c(1.45, 2.95, 4.45),
  xend = c(2.05, 3.55, 5.05),
  y = 1,
  yend = 1
)

panel_a <- ggplot() +
  geom_label(
    data = workflow_boxes,
    aes(x = x, y = y, label = label),
    size = 3.0,
    label.size = 0.35,
    label.padding = unit(0.22, "lines"),
    family = base_family,
    lineheight = 0.95
  ) +
  geom_segment(
    data = workflow_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.45,
    arrow = arrow(length = unit(0.12, "inches"))
  ) +
  annotate(
    "text",
    x = 3.25,
    y = 1.65,
    label = "Conventional and colocalization-supported gene assignments were compared",
    size = 3.1,
    fontface = "bold",
    family = base_family
  ) +
  coord_cartesian(xlim = c(0.35, 6.15), ylim = c(0.35, 1.9), clip = "off") +
  theme_void(base_family = base_family) +
  theme(plot.margin = margin(4, 8, 2, 8))

###############################################################################
# 8. Panel B – Main result
###############################################################################

main_result <- data.table(
  classification = factor(
    c("Reclassified", "Unchanged"),
    levels = c("Reclassified", "Unchanged")
  ),
  n = c(n_reclassified, n_unchanged),
  percent = c(percent_reclassified, percent_unchanged)
)

main_result[
  ,
  label := paste0(
    n,
    " (",
    sprintf("%.2f", percent),
    "%)"
  )
]

panel_b <- ggplot(
  main_result,
  aes(x = "", y = n, fill = classification)
) +
  geom_col(width = 0.72) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.2,
    fontface = "bold",
    family = base_family
  ) +
  annotate(
    "text",
    x = 1,
    y = n_evaluable + n_evaluable * 0.12,
    label = paste0(
      n_reclassified,
      " / ",
      n_evaluable,
      " clumps reassigned"
    ),
    size = 3.5,
    fontface = "bold",
    family = base_family
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18)),
    breaks = pretty_breaks(n = 4)
  ) +
  labs(
    title = "Overall gene reclassification",
    x = NULL,
    y = "Number of evaluable clumps"
  ) +
  guides(fill = guide_legend(reverse = FALSE)) +
  figure_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  )

###############################################################################
# 9. Panel C – Reclassification categories
###############################################################################

category_data <- type_summary[
  reclassification_type != "NOT_EVALUABLE"
]

category_data[
  ,
  category_label := clean_type_label(reclassification_type)
]

category_order <- c(
  "Different protein-coding gene",
  "Unchanged protein-coding gene",
  "Pseudogene",
  "Regulatory RNA",
  "Other or unknown"
)

category_data[
  ,
  category_label := factor(
    category_label,
    levels = rev(category_order)
  )
]

category_data[
  ,
  percent_evaluable := 100 * N / n_evaluable
]

panel_c <- ggplot(
  category_data,
  aes(x = N, y = category_label, fill = category_label)
) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_text(
    aes(
      label = paste0(
        N,
        " (",
        sprintf("%.1f", percent_evaluable),
        "%)"
      )
    ),
    hjust = -0.08,
    size = 3.0,
    family = base_family
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.22)),
    breaks = pretty_breaks(n = 5)
  ) +
  labs(
    title = "Reclassification outcome",
    x = "Number of evaluable clumps",
    y = NULL
  ) +
  figure_theme +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "none"
  )

###############################################################################
# 10. Panel D – Reclassification by priority class
###############################################################################

priority_plot_data <- copy(priority_summary)

priority_plot_data[
  ,
  priority_label := clean_priority_label(priority_class)
]

priority_plot_data[
  ,
  priority_label := factor(
    priority_label,
    levels = c("Strong H4", "Moderate H4", "Suggestive H4")
  )
]

priority_plot_data[
  ,
  n_unchanged := n_clumps - n_reclassified
]

priority_long <- melt(
  priority_plot_data,
  id.vars = c("priority_class", "priority_label", "n_clumps", "percent_reclassified"),
  measure.vars = c("n_reclassified", "n_unchanged"),
  variable.name = "classification",
  value.name = "n"
)

priority_long[
  ,
  classification := factor(
    classification,
    levels = c("n_reclassified", "n_unchanged"),
    labels = c("Reclassified", "Unchanged")
  )
]

panel_d <- ggplot(
  priority_long,
  aes(x = priority_label, y = n, fill = classification)
) +
  geom_col(width = 0.7) +
  geom_text(
    data = priority_plot_data,
    aes(
      x = priority_label,
      y = n_clumps,
      label = paste0(
        sprintf("%.1f", percent_reclassified),
        "%"
      )
    ),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 3.0,
    fontface = "bold",
    family = base_family
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18)),
    breaks = pretty_breaks(n = 5)
  ) +
  labs(
    title = "Reclassification by colocalization support",
    x = NULL,
    y = "Number of evaluable clumps"
  ) +
  figure_theme +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "bottom"
  )

###############################################################################
# 11. Assemble Figure 3
###############################################################################

message_section("ASSEMBLING FIGURE 3")

figure3 <- (
  panel_a /
    (panel_b | panel_c | panel_d)
) +
  plot_layout(
    heights = c(0.72, 2.15),
    widths = c(1.0, 1.3, 1.25),
    guides = "collect"
  ) +
  plot_annotation(
    title = "Colocalization-driven reassignment of candidate genes",
    subtitle = paste0(
      "Among ",
      n_evaluable,
      " clumps with H4-supported candidate genes, ",
      n_reclassified,
      " (",
      sprintf("%.2f", percent_reclassified),
      "%) differed from the conventional protein-coding gene assignment."
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(
        size = 13,
        face = "bold",
        family = base_family
      ),
      plot.subtitle = element_text(
        size = 9.5,
        family = base_family,
        margin = margin(b = 6)
      ),
      plot.tag = element_text(
        size = 12,
        face = "bold",
        family = base_family
      )
    )
  ) &
  theme(
    text = element_text(family = base_family),
    legend.position = "bottom"
  )

###############################################################################
# 12. Save outputs
###############################################################################

message_section("WRITING FIGURE OUTPUTS")

pdf_file <- file.path(
  figure_dir,
  "Figure3_Gene_Reclassification.pdf"
)

png_file <- file.path(
  figure_dir,
  "Figure3_Gene_Reclassification.png"
)

tiff_file <- file.path(
  figure_dir,
  "Figure3_Gene_Reclassification.tiff"
)

width_in <- 12.0
height_in <- 7.1

ggsave(
  filename = pdf_file,
  plot = figure3,
  width = width_in,
  height = height_in,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = png_file,
  plot = figure3,
  width = width_in,
  height = height_in,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = tiff_file,
  plot = figure3,
  width = width_in,
  height = height_in,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

###############################################################################
# 13. Figure data and metadata
###############################################################################

figure_data_file <- file.path(
  figure_dir,
  "Figure3_Gene_Reclassification_data.tsv"
)

figure_data <- rbindlist(
  list(
    data.table(
      panel = "B",
      category = as.character(main_result$classification),
      n = main_result$n,
      percent = main_result$percent
    ),
    data.table(
      panel = "C",
      category = as.character(category_data$category_label),
      n = category_data$N,
      percent = category_data$percent_evaluable
    ),
    data.table(
      panel = "D",
      category = as.character(priority_plot_data$priority_label),
      n = priority_plot_data$n_reclassified,
      percent = priority_plot_data$percent_reclassified
    )
  ),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  figure_data,
  figure_data_file,
  sep = "\t"
)

status_output <- data.table(
  module = "Module 13.1",
  analysis = "Figure 3 – Gene Reclassification",
  status = "SUCCESS",
  total_clumps = n_total,
  evaluable_clumps = n_evaluable,
  reclassified_clumps = n_reclassified,
  unchanged_clumps = n_unchanged,
  percent_reclassified = percent_reclassified,
  pdf_file = pdf_file,
  png_file = png_file,
  tiff_file = tiff_file,
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

fwrite(
  status_output,
  file.path(
    figure_dir,
    "Figure3_Gene_Reclassification_status.tsv"
  ),
  sep = "\t"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    figure_dir,
    "Figure3_Gene_Reclassification_sessionInfo.txt"
  )
)

###############################################################################
# 14. Console report
###############################################################################

message_section("MODULE 13.1 COMPLETED")

cat("Figure 3 PDF:\n", pdf_file, "\n\n")
cat("Figure 3 PNG:\n", png_file, "\n\n")
cat("Figure 3 TIFF:\n", tiff_file, "\n\n")
cat("Figure data:\n", figure_data_file, "\n\n")
cat("Main result: ", n_reclassified, " / ", n_evaluable,
    " clumps reclassified (", sprintf("%.2f%%", percent_reclassified), ")\n",
    sep = "")
