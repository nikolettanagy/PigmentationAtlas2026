###############################################################################
# PigmentationAtlas
# Figure 5 – Conceptual biological summary of colocalization-driven
# candidate-gene reassignment
#
# Purpose:
#   Create a publication-quality synthesis figure that summarizes how
#   PigmentationAtlas changes conventional proximity-based interpretation
#   of pigmentation GWAS clumps.
#
# Inputs:
#   results/module12/reclassification/
#     module12_6_clump_gene_reclassification.tsv.gz
#
# Outputs:
#   figures/Figure5_Conceptual_Biological_Summary.pdf
#   figures/Figure5_Conceptual_Biological_Summary.png
#   figures/Figure5_Conceptual_Biological_Summary.tiff
#   figures/figure_data/Figure5_*.tsv
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
# 2. Paths
###############################################################################

project_dir <- normalizePath(
  "C:/Users/User/Desktop/PigmentationAtlas",
  winslash = "/",
  mustWork = TRUE
)

input_file <- file.path(
  project_dir,
  "results",
  "module12",
  "reclassification",
  "module12_6_clump_gene_reclassification.tsv.gz"
)

figure_dir <- file.path(project_dir, "figures")
figure_data_dir <- file.path(figure_dir, "figure_data")

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_data_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_pdf <- file.path(
  figure_dir,
  "Figure5_Conceptual_Biological_Summary.pdf"
)

output_png <- file.path(
  figure_dir,
  "Figure5_Conceptual_Biological_Summary.png"
)

output_tiff <- file.path(
  figure_dir,
  "Figure5_Conceptual_Biological_Summary.tiff"
)

###############################################################################
# 3. Helpers
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

priority_order <- c(
  "STRONG_H4",
  "MODERATE_H4",
  "SUGGESTIVE_H4"
)

priority_labels <- c(
  "STRONG_H4" = "Strong H4",
  "MODERATE_H4" = "Moderate H4",
  "SUGGESTIVE_H4" = "Suggestive H4"
)

outcome_order <- c(
  "UNCHANGED_PROTEIN_CODING",
  "CODING_TO_OTHER_PROTEIN_CODING",
  "CODING_TO_PSEUDOGENE",
  "CODING_TO_REGULATORY_RNA",
  "CODING_TO_OTHER_OR_UNKNOWN"
)

outcome_labels <- c(
  "UNCHANGED_PROTEIN_CODING" = "Unchanged coding gene",
  "CODING_TO_OTHER_PROTEIN_CODING" = "Different coding gene",
  "CODING_TO_PSEUDOGENE" = "Pseudogene",
  "CODING_TO_REGULATORY_RNA" = "Regulatory RNA",
  "CODING_TO_OTHER_OR_UNKNOWN" = "Other / unknown"
)

# Palette chosen to match Figures 3–4.
outcome_colours <- c(
  "UNCHANGED_PROTEIN_CODING" = "#16B8BE",
  "CODING_TO_OTHER_PROTEIN_CODING" = "#F8766D",
  "CODING_TO_PSEUDOGENE" = "#B79F00",
  "CODING_TO_REGULATORY_RNA" = "#C77CFF",
  "CODING_TO_OTHER_OR_UNKNOWN" = "#7F7F7F"
)

theme_pigmentation <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = base_size + 1,
        margin = margin(b = 8)
      ),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      plot.margin = margin(8, 10, 8, 10)
    )
}

###############################################################################
# 4. Read input
###############################################################################

message_section("READING MODULE 12.6 RECLASSIFICATION TABLE")

assert_file_exists(
  input_file,
  "Module 12.6 reclassification table"
)

reclass <- fread(
  input_file,
  showProgress = TRUE
)

assert_columns(
  reclass,
  c(
    "clump_id",
    "baseline_gene_symbol",
    "prioritized_gene_symbol",
    "priority_class",
    "PP.H4",
    "evaluable",
    "reclassification_type",
    "interpretation_change"
  ),
  "Module 12.6 reclassification table"
)

reclass[
  ,
  evaluable := as.logical(evaluable)
]

reclass[
  ,
  interpretation_change := as.logical(interpretation_change)
]

evaluated <- reclass[
  evaluable == TRUE
]

if (nrow(evaluated) == 0L) {
  stop(
    "No evaluable clumps were found in Module 12.6 output.",
    call. = FALSE
  )
}

evaluated[
  ,
  priority_class := toupper(
    gsub("[ -]", "_", trimws(priority_class))
  )
]

evaluated[
  ,
  reclassification_type := toupper(
    gsub("[ -]", "_", trimws(reclassification_type))
  )
]

n_evaluable <- nrow(evaluated)
n_reclassified <- sum(
  evaluated$interpretation_change == TRUE,
  na.rm = TRUE
)
n_unchanged <- sum(
  evaluated$interpretation_change == FALSE,
  na.rm = TRUE
)
pct_reclassified <- 100 * n_reclassified / n_evaluable

cat("Evaluable clumps:   ", n_evaluable, "\n")
cat("Reclassified clumps:", n_reclassified, "\n")
cat("Unchanged clumps:   ", n_unchanged, "\n")
cat("Reclassified:       ", sprintf("%.2f%%", pct_reclassified), "\n")

###############################################################################
# 5. Figure data
###############################################################################

message_section("PREPARING FIGURE DATA")

outcome_summary <- evaluated[
  ,
  .N,
  by = reclassification_type
]

outcome_summary <- merge(
  data.table(reclassification_type = outcome_order),
  outcome_summary,
  by = "reclassification_type",
  all.x = TRUE
)

outcome_summary[
  is.na(N),
  N := 0L
]

outcome_summary[
  ,
  outcome_label := outcome_labels[reclassification_type]
]

outcome_summary[
  ,
  percent := 100 * N / n_evaluable
]

outcome_summary[
  ,
  outcome_label := factor(
    outcome_label,
    levels = rev(outcome_labels[outcome_order])
  )
]

priority_summary <- evaluated[
  priority_class %in% priority_order,
  .N,
  by = .(
    priority_class,
    reclassification_type
  )
]

priority_grid <- CJ(
  priority_class = priority_order,
  reclassification_type = outcome_order,
  unique = TRUE
)

priority_summary <- merge(
  priority_grid,
  priority_summary,
  by = c(
    "priority_class",
    "reclassification_type"
  ),
  all.x = TRUE
)

priority_summary[
  is.na(N),
  N := 0L
]

priority_summary[
  ,
  priority_label := priority_labels[priority_class]
]

priority_summary[
  ,
  outcome_label := outcome_labels[reclassification_type]
]

priority_summary[
  ,
  priority_label := factor(
    priority_label,
    levels = priority_labels[priority_order]
  )
]

priority_summary[
  ,
  outcome_label := factor(
    outcome_label,
    levels = outcome_labels[outcome_order]
  )
]

priority_totals <- priority_summary[
  ,
  .(total = sum(N)),
  by = priority_label
]

priority_summary <- merge(
  priority_summary,
  priority_totals,
  by = "priority_label",
  all.x = TRUE
)

priority_summary[
  ,
  percent_within_class := fifelse(
    total > 0,
    100 * N / total,
    0
  )
]

noncoding_summary <- outcome_summary[
  reclassification_type %in% c(
    "CODING_TO_PSEUDOGENE",
    "CODING_TO_REGULATORY_RNA",
    "CODING_TO_OTHER_OR_UNKNOWN"
  ),
  .(
    noncoding_or_other = sum(N)
  )
]

n_noncoding_or_other <- noncoding_summary$noncoding_or_other[[1L]]

fwrite(
  outcome_summary,
  file.path(
    figure_data_dir,
    "Figure5_outcome_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  priority_summary,
  file.path(
    figure_data_dir,
    "Figure5_priority_by_outcome.tsv"
  ),
  sep = "\t"
)

###############################################################################
# 6. Panel A – Conceptual comparison
###############################################################################

message_section("BUILDING PANEL A")

panel_a <- ggplot() +
  annotate(
    "rect",
    xmin = 0.4,
    xmax = 4.7,
    ymin = 1.25,
    ymax = 2.75,
    fill = "#F3F3F3",
    colour = "#808080",
    linewidth = 0.4
  ) +
  annotate(
    "rect",
    xmin = 5.3,
    xmax = 9.6,
    ymin = 1.25,
    ymax = 2.75,
    fill = "#EAF7F7",
    colour = "#16B8BE",
    linewidth = 0.5
  ) +
  annotate(
    "text",
    x = 2.55,
    y = 2.48,
    label = "Conventional interpretation",
    fontface = "bold",
    size = 4
  ) +
  annotate(
    "text",
    x = 7.45,
    y = 2.48,
    label = "PigmentationAtlas interpretation",
    fontface = "bold",
    size = 4
  ) +
  annotate(
    "label",
    x = 1.15,
    y = 1.85,
    label = "GWAS\nlead variant",
    size = 3.4,
    label.size = 0.25,
    fill = "white"
  ) +
  annotate(
    "segment",
    x = 1.8,
    xend = 2.55,
    y = 1.85,
    yend = 1.85,
    arrow = arrow(length = unit(0.12, "inches")),
    linewidth = 0.45
  ) +
  annotate(
    "label",
    x = 3.35,
    y = 1.85,
    label = "Overlapping or nearest\nprotein-coding gene",
    size = 3.25,
    label.size = 0.25,
    fill = "white"
  ) +
  annotate(
    "label",
    x = 6.05,
    y = 1.85,
    label = "GWAS clump +\nQTL evidence",
    size = 3.25,
    label.size = 0.25,
    fill = "white"
  ) +
  annotate(
    "segment",
    x = 6.75,
    xend = 7.5,
    y = 1.85,
    yend = 1.85,
    arrow = arrow(length = unit(0.12, "inches")),
    linewidth = 0.45
  ) +
  annotate(
    "label",
    x = 8.3,
    y = 1.85,
    label = "Colocalization-supported\ntarget gene",
    size = 3.25,
    label.size = 0.25,
    fill = "white"
  ) +
  annotate(
    "text",
    x = 5,
    y = 0.75,
    label = "Proximity-based annotation is replaced by evidence-based target-gene prioritization",
    fontface = "italic",
    size = 3.6
  ) +
  coord_cartesian(
    xlim = c(0, 10),
    ylim = c(0.4, 3.05),
    clip = "off"
  ) +
  labs(
    title = "From genomic proximity to colocalization-supported target genes"
  ) +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 11,
      margin = margin(b = 8)
    ),
    plot.margin = margin(8, 10, 8, 10)
  )

###############################################################################
# 7. Panel B – Outcome composition
###############################################################################

message_section("BUILDING PANEL B")

panel_b <- ggplot(
  outcome_summary[
    N > 0
  ],
  aes(
    x = N,
    y = outcome_label,
    fill = reclassification_type
  )
) +
  geom_col(
    width = 0.68,
    show.legend = FALSE
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%d (%.1f%%)",
        N,
        percent
      )
    ),
    hjust = -0.08,
    size = 3.3
  ) +
  scale_fill_manual(
    values = outcome_colours,
    drop = FALSE
  ) +
  scale_x_continuous(
    expand = expansion(
      mult = c(0, 0.22)
    )
  ) +
  labs(
    title = "Biological outcome of candidate-gene reassignment",
    x = "Number of evaluable clumps",
    y = NULL
  ) +
  theme_pigmentation(base_size = 10)

###############################################################################
# 8. Panel C – Evidence class × outcome
###############################################################################

message_section("BUILDING PANEL C")

panel_c <- ggplot(
  priority_summary[
    N > 0
  ],
  aes(
    x = priority_label,
    y = N,
    fill = reclassification_type
  )
) +
  geom_col(
    width = 0.68
  ) +
  geom_text(
    data = priority_totals,
    aes(
      x = priority_label,
      y = total,
      label = total
    ),
    inherit.aes = FALSE,
    vjust = -0.45,
    fontface = "bold",
    size = 3.2
  ) +
  scale_fill_manual(
    values = outcome_colours,
    labels = outcome_labels,
    breaks = outcome_order,
    drop = FALSE,
    name = "Outcome"
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.12)
    )
  ) +
  labs(
    title = "Reclassification is observed across all H4 evidence classes",
    x = NULL,
    y = "Number of evaluable clumps"
  ) +
  theme_pigmentation(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    )
  )

###############################################################################
# 9. Panel D – Take-home synthesis
###############################################################################

message_section("BUILDING PANEL D")

summary_text <- paste0(
  "Among ", n_evaluable,
  " clumps with H4-supported candidate genes,\n",
  n_reclassified, " (",
  sprintf("%.2f", pct_reclassified),
  "%) differed from the conventional assignment."
)

detail_text <- paste0(
  outcome_summary[
    reclassification_type == "CODING_TO_OTHER_PROTEIN_CODING",
    N
  ],
  " reassigned to another coding gene\n",
  outcome_summary[
    reclassification_type == "CODING_TO_PSEUDOGENE",
    N
  ],
  " reassigned to a pseudogene\n",
  outcome_summary[
    reclassification_type == "CODING_TO_REGULATORY_RNA",
    N
  ],
  " reassigned to a regulatory RNA"
)

panel_d <- ggplot() +
  annotate(
    "rect",
    xmin = 0.4,
    xmax = 9.6,
    ymin = 0.5,
    ymax = 4.6,
    fill = "#FAFAFA",
    colour = "#A0A0A0",
    linewidth = 0.45
  ) +
  annotate(
    "text",
    x = 5,
    y = 4.05,
    label = "Key biological interpretation",
    fontface = "bold",
    size = 4.4
  ) +
  annotate(
    "text",
    x = 5,
    y = 3.05,
    label = paste0(
      n_reclassified,
      " / ",
      n_evaluable
    ),
    fontface = "bold",
    size = 9
  ) +
  annotate(
    "text",
    x = 5,
    y = 2.45,
    label = paste0(
      sprintf("%.2f", pct_reclassified),
      "% reassigned"
    ),
    fontface = "bold",
    size = 5.2,
    colour = outcome_colours[
      "CODING_TO_OTHER_PROTEIN_CODING"
    ]
  ) +
  annotate(
    "text",
    x = 5,
    y = 1.65,
    label = summary_text,
    size = 3.35,
    lineheight = 1.08
  ) +
  annotate(
    "text",
    x = 5,
    y = 0.92,
    label = detail_text,
    size = 3.0,
    lineheight = 1.08
  ) +
  coord_cartesian(
    xlim = c(0, 10),
    ylim = c(0.25, 4.85),
    clip = "off"
  ) +
  labs(
    title = "PigmentationAtlas reveals widespread reinterpretation of GWAS loci"
  ) +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 11,
      margin = margin(b = 8)
    ),
    plot.margin = margin(8, 10, 8, 10)
  )

###############################################################################
# 10. Assemble figure
###############################################################################

message_section("ASSEMBLING FIGURE 5")

figure5 <- (
  panel_a
) / (
  panel_b | panel_c | panel_d
) +
  plot_layout(
    heights = c(0.9, 1.25),
    widths = c(1.05, 1.05, 1.1),
    guides = "collect"
  ) +
  plot_annotation(
    title = "PigmentationAtlas transforms proximity-based GWAS annotation into evidence-based target-gene prioritization",
    subtitle = paste0(
      "Colocalization altered the conventionally assigned candidate gene in ",
      n_reclassified,
      " of ",
      n_evaluable,
      " evaluable pigmentation GWAS clumps (",
      sprintf("%.2f", pct_reclassified),
      "%)."
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 18,
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = 10.5,
        margin = margin(b = 10)
      ),
      plot.tag = element_text(
        face = "plain",
        size = 15
      )
    )
  ) &
  theme(
    legend.position = "bottom"
  )

###############################################################################
# 11. Save outputs
###############################################################################

message_section("WRITING FIGURE OUTPUTS")

ggsave(
  filename = output_pdf,
  plot = figure5,
  width = 14,
  height = 9,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = output_png,
  plot = figure5,
  width = 14,
  height = 9,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = output_tiff,
  plot = figure5,
  width = 14,
  height = 9,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

status <- data.table(
  figure = "Figure 5",
  analysis = "Conceptual biological summary",
  status = "SUCCESS",
  evaluable_clumps = n_evaluable,
  reclassified_clumps = n_reclassified,
  unchanged_clumps = n_unchanged,
  reclassified_percent = pct_reclassified,
  input_file = input_file,
  pdf_file = output_pdf,
  png_file = output_png,
  tiff_file = output_tiff,
  generated_at = format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S"
  )
)

fwrite(
  status,
  file.path(
    figure_data_dir,
    "Figure5_status.tsv"
  ),
  sep = "\t"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    figure_data_dir,
    "Figure5_sessionInfo.txt"
  )
)

###############################################################################
# 12. Console report
###############################################################################

message_section("FIGURE 5 COMPLETED")

cat("Input:\n", input_file, "\n\n")
cat("Evaluable clumps:   ", n_evaluable, "\n")
cat("Reclassified clumps:", n_reclassified, "\n")
cat("Reclassified:       ", sprintf("%.2f%%", pct_reclassified), "\n\n")

cat("Outputs:\n")
cat(output_pdf, "\n")
cat(output_png, "\n")
cat(output_tiff, "\n")
