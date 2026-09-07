###############################################################
# PigmentationAtlas
# Supplementary Figure S6
# Additional genome-wide reclassification analyses
###############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

###############################################################
# Project paths
###############################################################

project_dir <- normalizePath(
  ".",
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

output_dir <- file.path(
  project_dir,
  "figures",
  "supplementary"
)

source_data_dir <- file.path(
  project_dir,
  "results",
  "figure_source_data"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  source_data_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

###############################################################
# Read data
###############################################################

dt <- fread(
  input_file,
  showProgress = FALSE
)

required_columns <- c(
  "clump_id",
  "parent_locus",
  "baseline_gene_symbol",
  "baseline_gene_biotype",
  "baseline_assignment",
  "prioritized_gene_symbol",
  "prioritized_gene_biotype",
  "candidate_class",
  "priority_class",
  "PP.H4",
  "PP.H3",
  "prioritized_biotype_class",
  "evaluable",
  "same_gene",
  "reclassified",
  "reclassification_type",
  "interpretation_change"
)

missing_columns <- setdiff(
  required_columns,
  names(dt)
)

if (length(missing_columns) > 0L) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

dt[, PP.H4 := as.numeric(PP.H4)]
dt[, PP.H3 := as.numeric(PP.H3)]

cat("\nInput dimensions:\n")
print(dim(dt))

cat("\nEvaluation status:\n")
print(dt[, .N, by = evaluable])

###############################################################
# Harmonized labels
###############################################################

dt[, evaluation_outcome :=
  fifelse(
    evaluable == FALSE | is.na(evaluable),
    "Not evaluable",
    fifelse(
      reclassified == TRUE,
      "Reclassified",
      "Unchanged"
    )
  )
]

evaluation_levels <- c(
  "Unchanged",
  "Reclassified",
  "Not evaluable"
)

dt[, evaluation_outcome := factor(
  evaluation_outcome,
  levels = evaluation_levels
)]

###############################################################
# Reclassification type labels
###############################################################

reclassification_label_map <- c(
  "CODING_TO_OTHER_PROTEIN_CODING" =
    "Coding to another\nprotein-coding gene",

  "CODING_TO_REGULATORY_RNA" =
    "Coding to\nregulatory RNA",

  "CODING_TO_PSEUDOGENE" =
    "Coding to\npseudogene",

  "CODING_TO_OTHER_NONCODING" =
    "Coding to other\nnon-coding gene",

  "UNCHANGED" =
    "Unchanged",

  "SAME_GENE" =
    "Unchanged",

  "NOT_EVALUABLE" =
    "Not evaluable"
)

dt[, reclassification_type_label :=
  reclassification_label_map[reclassification_type]
]

dt[
  is.na(reclassification_type_label),
  reclassification_type_label :=
    gsub("_", " ", reclassification_type)
]

###############################################################
# Priority class labels
###############################################################

priority_label_map <- c(
  "STRONG_H4" = "Strong H4",
  "MODERATE_H4" = "Moderate H4",
  "SUGGESTIVE_H4" = "Suggestive H4",
  "H3_DOMINANT" = "H3 dominant",
  "INCONCLUSIVE" = "Inconclusive"
)

dt[, priority_class_label :=
  priority_label_map[priority_class]
]

dt[
  is.na(priority_class_label) &
    !is.na(priority_class) &
    priority_class != "",
  priority_class_label :=
    tools::toTitleCase(
      tolower(
        gsub("_", " ", priority_class)
      )
    )
]

priority_levels <- c(
  "Strong H4",
  "Moderate H4",
  "Suggestive H4",
  "H3 dominant",
  "Inconclusive"
)

dt[, priority_class_label := factor(
  priority_class_label,
  levels = priority_levels
)]

###############################################################
# Colours
###############################################################

evaluation_colors <- c(
  "Unchanged" = "#4D4D4D",
  "Reclassified" = "#D55E00",
  "Not evaluable" = "#BDBDBD"
)

priority_colors <- c(
  "Strong H4" = "#0072B2",
  "Moderate H4" = "#56B4E9",
  "Suggestive H4" = "#009E73",
  "H3 dominant" = "#999999",
  "Inconclusive" = "#D0D0D0"
)

type_colors <- c(
  "Coding to another\nprotein-coding gene" = "#4C78A8",
  "Coding to\nregulatory RNA" = "#E45756",
  "Coding to\npseudogene" = "#B279A2",
  "Coding to other\nnon-coding gene" = "#72B7B2",
  "Unchanged" = "#777777",
  "Not evaluable" = "#CCCCCC"
)

base_theme <- theme_classic(
  base_size = 12
) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 10,
      margin = margin(b = 7)
    ),
    axis.title = element_text(
      size = 11
    ),
    axis.text = element_text(
      size = 9.5
    ),
    legend.title = element_text(
      face = "bold"
    )
  )

###############################################################
# Panel A
# Overall assignment outcome
###############################################################

panelA_data <- dt[
  ,
  .N,
  by = evaluation_outcome
][
  ,
  percentage := 100 * N / sum(N)
]

panelA_data[, label := paste0(
  comma(N),
  "\n",
  sprintf("%.1f%%", percentage)
)]

pA <- ggplot(
  panelA_data,
  aes(
    x = evaluation_outcome,
    y = N,
    fill = evaluation_outcome
  )
) +
  geom_col(
    width = 0.68
  ) +
  geom_text(
    aes(label = label),
    vjust = -0.35,
    size = 3.8
  ) +
  scale_fill_manual(
    values = evaluation_colors,
    guide = "none"
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.14)
    ),
    labels = comma
  ) +
  labs(
    title = "A  Genome-wide assignment outcome",
    subtitle = "Colocalization-based evaluation of conventional gene assignments",
    x = NULL,
    y = "Number of clumps"
  ) +
  base_theme

###############################################################
# Panel B
# Reclassification types
###############################################################

panelB_data <- dt[
  reclassified == TRUE &
    !is.na(reclassification_type_label),
  .N,
  by = reclassification_type_label
][
  order(-N)
]

panelB_data[, percentage := 100 * N / sum(N)]

panelB_data[, reclassification_type_label :=
  factor(
    reclassification_type_label,
    levels = rev(reclassification_type_label)
  )
]

panelB_data[, label := paste0(
  comma(N),
  " (",
  sprintf("%.1f%%", percentage),
  ")"
)]

pB <- ggplot(
  panelB_data,
  aes(
    x = N,
    y = reclassification_type_label,
    fill = reclassification_type_label
  )
) +
  geom_col(
    width = 0.68
  ) +
  geom_text(
    aes(label = label),
    hjust = -0.08,
    size = 3.5
  ) +
  scale_fill_manual(
    values = type_colors,
    guide = "none"
  ) +
  scale_x_continuous(
    expand = expansion(
      mult = c(0, 0.24)
    ),
    labels = comma
  ) +
  labs(
    title = "B  Types of gene reassignment",
    subtitle = "Composition of reclassified clumps",
    x = "Number of reclassified clumps",
    y = NULL
  ) +
  base_theme

###############################################################
# Panel C
# Reclassification by priority class
###############################################################

panelC_data <- dt[
  evaluable == TRUE &
    !is.na(priority_class_label),
  .N,
  by = .(
    priority_class_label,
    evaluation_outcome
  )
]

panelC_data[
  ,
  total_priority := sum(N),
  by = priority_class_label
]

panelC_data[
  ,
  percentage := 100 * N / total_priority
]

pC <- ggplot(
  panelC_data,
  aes(
    x = priority_class_label,
    y = N,
    fill = evaluation_outcome
  )
) +
  geom_col(
    position = "stack",
    width = 0.72
  ) +
  scale_fill_manual(
    values = evaluation_colors,
    drop = FALSE,
    name = "Assignment outcome"
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(
      mult = c(0, 0.08)
    )
  ) +
  labs(
    title = "C  Reclassification by evidence strength",
    subtitle = "Assignment outcome across colocalization priority classes",
    x = NULL,
    y = "Number of evaluable clumps"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    legend.position = "bottom"
  )

###############################################################
# Panel D
# PP.H4 distribution by outcome
###############################################################

panelD_data <- dt[
  evaluable == TRUE &
    is.finite(PP.H4) &
    evaluation_outcome %in% c(
      "Unchanged",
      "Reclassified"
    )
]

panelD_summary <- panelD_data[
  ,
  .(
    n = .N,
    median_PP_H4 = median(PP.H4),
    mean_PP_H4 = mean(PP.H4)
  ),
  by = evaluation_outcome
]

cat("\nPP.H4 summary by assignment outcome:\n")
print(panelD_summary)

pD <- ggplot(
  panelD_data,
  aes(
    x = evaluation_outcome,
    y = PP.H4,
    fill = evaluation_outcome
  )
) +
  geom_violin(
    trim = TRUE,
    scale = "width",
    alpha = 0.80,
    colour = "grey25",
    linewidth = 0.40
  ) +
  geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.42
  ) +
  geom_hline(
    yintercept = c(
      0.20,
      0.50,
      0.80
    ),
    linetype = "dotted",
    linewidth = 0.40
  ) +
  scale_fill_manual(
    values = evaluation_colors,
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  labs(
    title = "D  Shared-causal support by assignment outcome",
    subtitle = "PP.H4 distributions among evaluable clumps",
    x = NULL,
    y = "Shared-causal posterior probability (PP.H4)"
  ) +
  base_theme

###############################################################
# Optional statistical comparison
###############################################################

if (
  all(
    c(
      "Unchanged",
      "Reclassified"
    ) %in% unique(
      as.character(
        panelD_data$evaluation_outcome
      )
    )
  )
) {

  wilcox_result <- wilcox.test(
    PP.H4 ~ evaluation_outcome,
    data = panelD_data,
    exact = FALSE
  )

  wilcox_label <- if (
    wilcox_result$p.value < 0.001
  ) {
    "Wilcoxon P < 0.001"
  } else {
    paste0(
      "Wilcoxon P = ",
      formatC(
        wilcox_result$p.value,
        digits = 3,
        format = "fg"
      )
    )
  }

  pD <- pD +
    annotate(
      "text",
      x = 1.5,
      y = 0.97,
      label = wilcox_label,
      size = 3.5
    )
}

###############################################################
# Combined figure
###############################################################

combined_figure <- (
  pA |
    pB
) /
  (
    pC |
      pD
  ) +
  plot_layout(
    widths = c(1, 1),
    heights = c(1, 1.05)
  ) +
  plot_annotation(
    title = "Additional genome-wide reclassification analyses",
    subtitle = paste0(
      "Comparison of conventional protein-coding gene assignments ",
      "with colocalization-prioritized candidate genes"
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 19
      ),
      plot.subtitle = element_text(
        size = 11,
        margin = margin(b = 12)
      )
    )
  )

###############################################################
# Source-data tables
###############################################################

fwrite(
  panelA_data,
  file.path(
    source_data_dir,
    "Supplementary_Figure_S6_overall_assignment_outcome.tsv"
  ),
  sep = "\t"
)

fwrite(
  panelB_data,
  file.path(
    source_data_dir,
    "Supplementary_Figure_S6_reclassification_types.tsv"
  ),
  sep = "\t"
)

fwrite(
  panelC_data,
  file.path(
    source_data_dir,
    "Supplementary_Figure_S6_reclassification_by_priority_class.tsv"
  ),
  sep = "\t"
)

fwrite(
  panelD_summary,
  file.path(
    source_data_dir,
    "Supplementary_Figure_S6_PP_H4_by_assignment_outcome.tsv"
  ),
  sep = "\t"
)

fwrite(
  dt[
    ,
    .(
      clump_id,
      parent_locus,
      chromosome,
      index_variant,
      baseline_gene_symbol,
      baseline_gene_biotype,
      baseline_assignment,
      prioritized_gene_symbol,
      prioritized_gene_biotype,
      candidate_class,
      priority_class,
      PP.H4,
      PP.H3,
      evaluable,
      same_gene,
      reclassified,
      reclassification_type,
      interpretation_change
    )
  ],
  file.path(
    source_data_dir,
    "Supplementary_Figure_S6_reclassification_source_data.tsv.gz"
  ),
  sep = "\t"
)

###############################################################
# Save figure
###############################################################

output_base <- file.path(
  output_dir,
  "Supplementary_Figure_S6_genomewide_reclassification_analyses"
)

ggsave(
  filename = paste0(
    output_base,
    ".pdf"
  ),
  plot = combined_figure,
  width = 14,
  height = 11,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = paste0(
    output_base,
    ".png"
  ),
  plot = combined_figure,
  width = 14,
  height = 11,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = paste0(
    output_base,
    ".tiff"
  ),
  plot = combined_figure,
  width = 14,
  height = 11,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

###############################################################
# Console report
###############################################################

cat("\n=============================================\n")
cat("Supplementary Figure S6 completed.\n")
cat("Total clumps: ", nrow(dt), "\n", sep = "")
cat(
  "Evaluable clumps: ",
  sum(dt$evaluable == TRUE, na.rm = TRUE),
  "\n",
  sep = ""
)
cat(
  "Reclassified clumps: ",
  sum(dt$reclassified == TRUE, na.rm = TRUE),
  "\n",
  sep = ""
)
cat(
  "Unchanged clumps: ",
  sum(
    dt$evaluable == TRUE &
      dt$reclassified == FALSE,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)
cat("Outputs:\n")
cat(paste0(output_base, ".pdf\n"))
cat(paste0(output_base, ".png\n"))
cat(paste0(output_base, ".tiff\n"))
cat("=============================================\n")