###############################################################
# PigmentationAtlas
# Supplementary Figure S4
# Tissue-specific comparison of colocalization evidence
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

project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)

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

dir.create(
  output_dir,
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
  "tissue",
  "phenotype_id",
  "gene_id",
  "PP.H3",
  "PP.H4"
)

missing_columns <- setdiff(required_columns, names(dt))

if (length(missing_columns) > 0L) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

dt[, PP.H3 := as.numeric(PP.H3)]
dt[, PP.H4 := as.numeric(PP.H4)]

dt <- dt[
  is.finite(PP.H3) &
  is.finite(PP.H4) &
  !is.na(tissue)
]

cat("\nInput dimensions:\n")
print(dim(dt))

cat("\nDetected tissues:\n")
print(dt[, .N, by = tissue][order(-N)])

###############################################################
# Tissue labels
###############################################################

tissue_lookup <- data.table(
  tissue = c(
    "Skin_Sun_Exposed_Lower_leg",
    "Skin_Not_Sun_Exposed_Suprapubic"
  ),
  tissue_short = c(
    "Sun-exposed skin",
    "Non-sun-exposed skin"
  )
)

dt <- merge(
  dt,
  tissue_lookup,
  by = "tissue",
  all.x = TRUE,
  sort = FALSE
)

dt[
  is.na(tissue_short),
  tissue_short := gsub("_", " ", tissue)
]

tissue_levels <- c(
  "Sun-exposed skin",
  "Non-sun-exposed skin"
)

dt[, tissue_short := factor(
  tissue_short,
  levels = tissue_levels
)]

###############################################################
# Evidence classes
###############################################################

dt[, evidence_class :=
  fifelse(
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
  )
]

dt[, evidence_class := factor(
  evidence_class,
  levels = c(
    "Strong H4",
    "Moderate H4",
    "Suggestive H4",
    "H3 dominant",
    "Inconclusive"
  )
)]

###############################################################
# Summary statistics by tissue
###############################################################

tissue_summary <- dt[
  ,
  .(
    analyses = .N,
    median_PP_H4 = median(PP.H4, na.rm = TRUE),
    mean_PP_H4 = mean(PP.H4, na.rm = TRUE),
    suggestive_or_higher = sum(PP.H4 >= 0.20),
    moderate_or_higher = sum(PP.H4 >= 0.50),
    strong = sum(PP.H4 >= 0.80)
  ),
  by = .(tissue_short)
]

cat("\nTissue-level summary:\n")
print(tissue_summary)

fwrite(
  tissue_summary,
  file.path(
    project_dir,
    "results",
    "Supplementary_Figure_S4_tissue_summary.tsv"
  ),
  sep = "\t"
)

###############################################################
# Panel A
# PP.H4 distribution by tissue
###############################################################

pA <- ggplot(
  dt,
  aes(
    x = tissue_short,
    y = PP.H4
  )
) +
  geom_violin(
    trim = TRUE,
    scale = "width",
    fill = "grey85",
    colour = "grey25",
    linewidth = 0.45
  ) +
  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.45
  ) +
  geom_hline(
    yintercept = c(0.20, 0.50, 0.80),
    linetype = "dotted",
    linewidth = 0.45
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "A  Distribution of shared-causal evidence",
    subtitle = "PP.H4 across all gene–clump combinations",
    x = NULL,
    y = "Shared-causal posterior probability (PP.H4)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 10
    ),
    axis.text.x = element_text(
      size = 10
    )
  )

###############################################################
# Prepare paired tissue comparison
###############################################################

pair_keys <- c(
  "clump_id",
  "phenotype_id",
  "gene_id"
)

# Defensive aggregation:
# if duplicate rows occur within a tissue, retain the maximum PP.H4
paired_long <- dt[
  ,
  .(
    PP.H4 = max(PP.H4, na.rm = TRUE)
  ),
  by = c(pair_keys, "tissue_short")
]

paired_wide <- dcast(
  paired_long,
  clump_id + phenotype_id + gene_id ~ tissue_short,
  value.var = "PP.H4"
)

expected_pair_columns <- tissue_levels

missing_tissue_columns <- setdiff(
  expected_pair_columns,
  names(paired_wide)
)

if (length(missing_tissue_columns) > 0L) {
  stop(
    "Could not construct paired comparison. Missing tissue columns: ",
    paste(missing_tissue_columns, collapse = ", ")
  )
}

setnames(
  paired_wide,
  old = tissue_levels,
  new = c("PP_H4_sun", "PP_H4_non_sun")
)

paired <- paired_wide[
  is.finite(PP_H4_sun) &
  is.finite(PP_H4_non_sun)
]

paired[, delta_PP_H4 := PP_H4_sun - PP_H4_non_sun]
paired[, abs_delta_PP_H4 := abs(delta_PP_H4)]

cat("\nNumber of paired gene–clump comparisons: ",
    format(nrow(paired), big.mark = ","), "\n", sep = "")

if (nrow(paired) < 3L) {
  stop("Too few paired observations for tissue comparison.")
}

spearman_result <- cor.test(
  paired$PP_H4_sun,
  paired$PP_H4_non_sun,
  method = "spearman",
  exact = FALSE
)

spearman_rho <- unname(spearman_result$estimate)
spearman_p <- spearman_result$p.value

correlation_label <- paste0(
  "Spearman \u03c1 = ",
  formatC(spearman_rho, digits = 2, format = "f"),
  "\nP ",
  if (spearman_p < 0.001) {
    "< 0.001"
  } else {
    paste0(
      "= ",
      formatC(spearman_p, digits = 2, format = "fg")
    )
  },
  "\nPaired comparisons = ",
  format(nrow(paired), big.mark = ",")
)

###############################################################
# Panel B
# Sun versus non-sun paired PP.H4
###############################################################

pB <- ggplot(
  paired,
  aes(
    x = PP_H4_non_sun,
    y = PP_H4_sun
  )
) +
  geom_point(
    shape = 16,
    size = 1.15,
    alpha = 0.22
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.55
  ) +
  geom_vline(
    xintercept = c(0.20, 0.50, 0.80),
    linetype = "dotted",
    linewidth = 0.35
  ) +
  geom_hline(
    yintercept = c(0.20, 0.50, 0.80),
    linetype = "dotted",
    linewidth = 0.35
  ) +
  annotate(
    "text",
    x = 0.97,
    y = 0.03,
    label = correlation_label,
    hjust = 1,
    vjust = 0,
    size = 3.4
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
    title = "B  Concordance between skin tissues",
    subtitle = "Matched gene–clump comparisons",
    x = "PP.H4 in non-sun-exposed skin",
    y = "PP.H4 in sun-exposed skin"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 10
    )
  )

###############################################################
# Panel C
# Difference in PP.H4 between tissues
###############################################################

median_delta <- median(
  paired$delta_PP_H4,
  na.rm = TRUE
)

delta_label <- paste0(
  "Median \u0394PP.H4 = ",
  formatC(
    median_delta,
    digits = 3,
    format = "f"
  )
)

pC <- ggplot(
  paired,
  aes(x = delta_PP_H4)
) +
  geom_histogram(
    bins = 60,
    fill = "grey45",
    colour = "white",
    linewidth = 0.15
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_vline(
    xintercept = median_delta,
    linetype = "dotted",
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = delta_label,
    hjust = 1.08,
    vjust = 1.5,
    size = 3.5
  ) +
  scale_x_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.25)
  ) +
  labs(
    title = "C  Tissue-specific difference in PP.H4",
    subtitle = "Positive values indicate stronger evidence in sun-exposed skin",
    x = expression(
      Delta * "PP.H4 (sun-exposed \u2212 non-sun-exposed)"
    ),
    y = "Number of paired comparisons"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 10
    )
  )

###############################################################
# Combined figure
###############################################################

figure_title <- paste0(
  "Tissue-specific comparison of colocalization evidence"
)

figure_subtitle <- paste0(
  "Shared-causal posterior probabilities across matched ",
  "PigmentationAtlas gene–clump analyses"
)

fig <- (
  pA |
  pB
) /
  pC +
  plot_layout(
    heights = c(1, 0.82)
  ) +
  plot_annotation(
    title = figure_title,
    subtitle = figure_subtitle,
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
# Save outputs
###############################################################

output_base <- file.path(
  output_dir,
  "Supplementary_Figure_S4_tissue_specific_colocalization"
)

ggsave(
  filename = paste0(output_base, ".pdf"),
  plot = fig,
  width = 13,
  height = 10,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = paste0(output_base, ".png"),
  plot = fig,
  width = 13,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = paste0(output_base, ".tiff"),
  plot = fig,
  width = 13,
  height = 10,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

fwrite(
  paired[
    order(-abs_delta_PP_H4)
  ],
  file.path(
    project_dir,
    "results",
    "Supplementary_Figure_S4_paired_tissue_comparison.tsv.gz"
  ),
  sep = "\t"
)

cat("\n========================================\n")
cat("Supplementary Figure S4 completed.\n")
cat("Paired comparisons: ", nrow(paired), "\n", sep = "")
cat(
  "Spearman rho: ",
  formatC(spearman_rho, digits = 4, format = "f"),
  "\n",
  sep = ""
)
cat(
  "Median delta PP.H4: ",
  formatC(median_delta, digits = 4, format = "f"),
  "\n",
  sep = ""
)
cat("Outputs:\n")
cat(paste0(output_base, ".pdf\n"))
cat(paste0(output_base, ".png\n"))
cat(paste0(output_base, ".tiff\n"))
cat("========================================\n")