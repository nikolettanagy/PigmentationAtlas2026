###############################################################
# PigmentationAtlas
# Supplementary Figure S5
# Extended H3/H4 posterior probability analyses
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
  "colocalization",
  "catalogue",
  "summary",
  "module12_3_coloc_catalogue_summary.tsv.gz"
)

output_dir <- file.path(
  project_dir,
  "figures",
  "supplementary"
)

results_dir <- file.path(
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
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

###############################################################
# Read input
###############################################################

dt <- fread(
  input_file,
  showProgress = FALSE
)

required_columns <- c(
  "coloc_unit_id",
  "clump_id",
  "tissue",
  "phenotype_id",
  "gene_id",
  paste0("PP.H", 0:4)
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

posterior_columns <- paste0("PP.H", 0:4)

for (column_name in posterior_columns) {
  set(
    dt,
    j = column_name,
    value = as.numeric(dt[[column_name]])
  )
}

dt <- dt[
  complete.cases(dt[, ..posterior_columns])
]

for (column_name in posterior_columns) {
  dt <- dt[
    is.finite(get(column_name))
  ]
}

cat("\nInput dimensions after QC:\n")
print(dim(dt))

###############################################################
# Posterior sum QC
###############################################################

dt[, posterior_sum_calculated := rowSums(
  .SD
), .SDcols = posterior_columns]

cat("\nPosterior probability sum QC:\n")
cat(
  "Median posterior sum: ",
  median(dt$posterior_sum_calculated),
  "\n",
  sep = ""
)

cat(
  "Maximum absolute deviation from 1: ",
  max(abs(dt$posterior_sum_calculated - 1)),
  "\n",
  sep = ""
)

###############################################################
# Derived H3/H4 measures
###############################################################

dt[, H3_H4_sum := PP.H3 + PP.H4]

dt[, conditional_H4 := fifelse(
  H3_H4_sum > 0,
  PP.H4 / H3_H4_sum,
  NA_real_
)]

# Numerical protection
dt[
  conditional_H4 < 0,
  conditional_H4 := 0
]

dt[
  conditional_H4 > 1,
  conditional_H4 := 1
]

###############################################################
# Posterior entropy
#
# Raw entropy:
#   -sum(p * log(p))
#
# Normalized by log(5), producing a range of 0–1:
#   0 = posterior concentrated in one hypothesis
#   1 = equal support across all five hypotheses
###############################################################

posterior_matrix <- as.matrix(
  dt[, ..posterior_columns]
)

posterior_matrix[
  posterior_matrix <= 0
] <- NA_real_

entropy_raw <- -rowSums(
  posterior_matrix * log(posterior_matrix),
  na.rm = TRUE
)

dt[, posterior_entropy := entropy_raw / log(5)]

dt[
  posterior_entropy < 0,
  posterior_entropy := 0
]

dt[
  posterior_entropy > 1,
  posterior_entropy := 1
]

###############################################################
# Evidence classes
###############################################################

dt[, evidence_class := fifelse(
  PP.H4 >= 0.80,
  "Strong H4",
  fifelse(
    PP.H4 >= 0.50,
    "Moderate H4",
    fifelse(
      PP.H4 >= 0.20,
      "Suggestive H4",
      fifelse(
        PP.H3 > PP.H4,
        "H3 dominant",
        "Inconclusive"
      )
    )
  )
)]

evidence_levels <- c(
  "Strong H4",
  "Moderate H4",
  "Suggestive H4",
  "H3 dominant",
  "Inconclusive"
)

dt[, evidence_class := factor(
  evidence_class,
  levels = evidence_levels
)]

###############################################################
# Tissue labels
###############################################################

dt[, tissue_short := fifelse(
  tissue == "Skin_Sun_Exposed_Lower_leg",
  "Sun-exposed skin",
  fifelse(
    tissue == "Skin_Not_Sun_Exposed_Suprapubic",
    "Non-sun-exposed skin",
    gsub("_", " ", tissue)
  )
)]

###############################################################
# Global graphical settings
###############################################################

threshold_values <- c(
  0.20,
  0.50,
  0.80
)

evidence_colors <- c(
  "Strong H4" = "#0877B9",
  "Moderate H4" = "#56B4E9",
  "Suggestive H4" = "#009E73",
  "H3 dominant" = "#999999",
  "Inconclusive" = "#C7C7C7"
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
      margin = margin(
        b = 7
      )
    ),
    axis.title = element_text(
      size = 11
    ),
    axis.text = element_text(
      size = 9.5
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )

###############################################################
# Panel A
# Two-dimensional PP.H3 versus PP.H4 density
###############################################################

pA <- ggplot(
  dt,
  aes(
    x = PP.H3,
    y = PP.H4
  )
) +
  geom_bin_2d(
    bins = 55
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.55
  ) +
  geom_vline(
    xintercept = threshold_values,
    linetype = "dotted",
    linewidth = 0.30
  ) +
  geom_hline(
    yintercept = threshold_values,
    linetype = "dotted",
    linewidth = 0.30
  ) +
  scale_fill_viridis_c(
    option = "C",
    trans = "log10",
    name = "Analyses\nper bin",
    labels = comma
  ) +
  coord_equal(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, 0.2)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    title = "A  H3–H4 posterior density",
    subtitle = "Dashed diagonal indicates equal support for H3 and H4",
    x = "PP.H3: distinct causal variants",
    y = "PP.H4: shared causal variant"
  ) +
  base_theme

###############################################################
# Panel B
# Conditional H4 probability
###############################################################

conditional_dt <- dt[
  is.finite(conditional_H4)
]

conditional_median <- median(
  conditional_dt$conditional_H4,
  na.rm = TRUE
)

pB <- ggplot(
  conditional_dt,
  aes(
    x = conditional_H4
  )
) +
  geom_histogram(
    bins = 60,
    fill = "grey45",
    colour = "white",
    linewidth = 0.15
  ) +
  geom_vline(
    xintercept = 0.50,
    linetype = "dashed",
    linewidth = 0.65
  ) +
  geom_vline(
    xintercept = conditional_median,
    linetype = "dotted",
    linewidth = 0.70
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0(
      "Median = ",
      formatC(
        conditional_median,
        digits = 3,
        format = "f"
      )
    ),
    hjust = 1.08,
    vjust = 1.45,
    size = 3.4
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(
      mult = c(0, 0.01)
    )
  ) +
  labs(
    title = "B  Conditional support for H4",
    subtitle = "Relative support for shared versus distinct causal variants",
    x = expression(
      frac(
        "PP.H4",
        "PP.H3 + PP.H4"
      )
    ),
    y = "Number of analyses"
  ) +
  base_theme

###############################################################
# Panel C
# Posterior entropy
###############################################################

entropy_median <- median(
  dt$posterior_entropy,
  na.rm = TRUE
)

pC <- ggplot(
  dt,
  aes(
    x = posterior_entropy,
    fill = evidence_class
  )
) +
  geom_histogram(
    bins = 55,
    position = "stack",
    colour = "white",
    linewidth = 0.10
  ) +
  geom_vline(
    xintercept = entropy_median,
    linetype = "dashed",
    linewidth = 0.65
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0(
      "Median normalized entropy = ",
      formatC(
        entropy_median,
        digits = 3,
        format = "f"
      )
    ),
    hjust = 1.05,
    vjust = 1.45,
    size = 3.3
  ) +
  scale_fill_manual(
    values = evidence_colors,
    drop = FALSE,
    name = "Evidence class"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(
      mult = c(0, 0.01)
    )
  ) +
  labs(
    title = "C  Posterior uncertainty",
    subtitle = paste0(
      "Normalized entropy: 0 = concentrated evidence; ",
      "1 = equal support across H0–H4"
    ),
    x = "Normalized posterior entropy",
    y = "Number of analyses"
  ) +
  base_theme

###############################################################
# Panel D
# Top PP.H4 analyses
###############################################################

top_n <- 30L

top_results <- copy(
  dt[
    order(
      -PP.H4,
      -conditional_H4,
      PP.H3
    )
  ][
    seq_len(
      min(
        top_n,
        .N
      )
    )
  ]
)

# Prefer gene symbol where available; otherwise use gene_id
gene_symbol_candidates <- c(
  "gene_name",
  "gene_symbol",
  "symbol"
)

available_symbol_column <- intersect(
  gene_symbol_candidates,
  names(top_results)
)

if (length(available_symbol_column) > 0L) {

  selected_symbol_column <- available_symbol_column[1]

  top_results[, gene_display := as.character(
    get(selected_symbol_column)
  )]

  top_results[
    is.na(gene_display) |
      gene_display == "",
    gene_display := gene_id
  ]

} else {

  top_results[, gene_display := gene_id]

}

# Make labels unique where the same gene occurs more than once
top_results[, repeated_gene := duplicated(
  gene_display
) | duplicated(
  gene_display,
  fromLast = TRUE
)]

top_results[
  repeated_gene == TRUE,
  gene_display := paste0(
    gene_display,
    " | ",
    tissue_short,
    " | ",
    clump_id
  )
]

top_results[, result_label := factor(
  gene_display,
  levels = rev(gene_display)
)]

pD <- ggplot(
  top_results,
  aes(
    x = PP.H4,
    y = result_label,
    colour = evidence_class
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = PP.H4,
      y = result_label,
      yend = result_label
    ),
    colour = "grey78",
    linewidth = 0.55
  ) +
  geom_point(
    size = 2.7
  ) +
  geom_vline(
    xintercept = threshold_values,
    linetype = "dotted",
    linewidth = 0.35
  ) +
  scale_colour_manual(
    values = evidence_colors,
    drop = FALSE,
    name = "Evidence class"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  labs(
    title = paste0(
      "D  Top ",
      nrow(top_results),
      " shared-causal signals"
    ),
    subtitle = "Gene–clump–tissue analyses ranked by PP.H4",
    x = "Shared-causal posterior probability (PP.H4)",
    y = NULL
  ) +
  base_theme +
  theme(
    axis.text.y = element_text(
      size = 7.5
    ),
    legend.position = "none"
  )

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
    heights = c(1, 1.25),
    guides = "collect"
  ) +
  plot_annotation(
    title = "Extended H3/H4 posterior probability analyses",
    subtitle = paste0(
      "Relative support, uncertainty, and strongest shared-causal ",
      "signals across 7,840 PigmentationAtlas colocalization analyses"
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 19
      ),
      plot.subtitle = element_text(
        size = 11,
        margin = margin(
          b = 12
        )
      )
    )
  ) &
  theme(
    legend.position = "bottom"
  )

###############################################################
# Source-data tables
###############################################################

global_summary <- data.table(
  metric = c(
    "Number of analyses",
    "Median PP.H3",
    "Median PP.H4",
    "Median conditional H4",
    "Median normalized entropy",
    "H4 dominant conditional probability count",
    "H3 dominant conditional probability count"
  ),
  value = c(
    nrow(dt),
    median(dt$PP.H3, na.rm = TRUE),
    median(dt$PP.H4, na.rm = TRUE),
    conditional_median,
    entropy_median,
    sum(dt$conditional_H4 > 0.50, na.rm = TRUE),
    sum(dt$conditional_H4 < 0.50, na.rm = TRUE)
  )
)

fwrite(
  global_summary,
  file.path(
    results_dir,
    "Supplementary_Figure_S5_global_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  top_results[
    ,
    .(
      coloc_unit_id,
      clump_id,
      tissue,
      phenotype_id,
      gene_id,
      gene_display,
      PP.H3,
      PP.H4,
      conditional_H4,
      posterior_entropy,
      evidence_class
    )
  ],
  file.path(
    results_dir,
    "Supplementary_Figure_S5_top_colocalization_results.tsv"
  ),
  sep = "\t"
)

fwrite(
  dt[
    ,
    .(
      coloc_unit_id,
      clump_id,
      tissue,
      phenotype_id,
      gene_id,
      PP.H3,
      PP.H4,
      conditional_H4,
      posterior_entropy,
      evidence_class
    )
  ],
  file.path(
    results_dir,
    "Supplementary_Figure_S5_extended_H3_H4_source_data.tsv.gz"
  ),
  sep = "\t"
)

###############################################################
# Save figure
###############################################################

output_base <- file.path(
  output_dir,
  "Supplementary_Figure_S5_extended_H3_H4_posterior_analyses"
)

ggsave(
  filename = paste0(
    output_base,
    ".pdf"
  ),
  plot = combined_figure,
  width = 14,
  height = 12,
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
  height = 12,
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
  height = 12,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

###############################################################
# Console report
###############################################################

cat("\n=============================================\n")
cat("Supplementary Figure S5 completed.\n")
cat("Analyses: ", format(nrow(dt), big.mark = ","), "\n", sep = "")
cat(
  "Median PP.H3: ",
  formatC(
    median(dt$PP.H3),
    digits = 4,
    format = "f"
  ),
  "\n",
  sep = ""
)
cat(
  "Median PP.H4: ",
  formatC(
    median(dt$PP.H4),
    digits = 4,
    format = "f"
  ),
  "\n",
  sep = ""
)
cat(
  "Median conditional H4: ",
  formatC(
    conditional_median,
    digits = 4,
    format = "f"
  ),
  "\n",
  sep = ""
)
cat(
  "Median normalized entropy: ",
  formatC(
    entropy_median,
    digits = 4,
    format = "f"
  ),
  "\n",
  sep = ""
)
cat("Outputs:\n")
cat(paste0(output_base, ".pdf\n"))
cat(paste0(output_base, ".png\n"))
cat(paste0(output_base, ".tiff\n"))
cat("=============================================\n")