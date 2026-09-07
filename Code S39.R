###############################################################################
# PigmentationAtlas
# Supplementary Figure S2
#
# Extended posterior probability analyses across all completed
# molecular QTL colocalization tests
#
# Panels:
# A. Distribution of PP.H0–PP.H4
# B. Relationship between PP.H3 and PP.H4
# C. PP.H4 distribution by evidence class
# D. Ranked PP.H4 profile
#
# Outputs:
#   figures/Supplementary_Figure_S2_posterior_probabilities.pdf
#   figures/Supplementary_Figure_S2_posterior_probabilities.png
#   figures/Supplementary_Figure_S2_posterior_probabilities.tiff
#   results/Supplementary_Figure_S2_summary_statistics.tsv
###############################################################################

rm(list = ls())
gc()

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

###############################################################################
# 1. Packages
###############################################################################

required_packages <- c(
  "data.table",
  "ggplot2",
  "patchwork",
  "scales"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)

###############################################################################
# 2. Project directories
###############################################################################

project_dir <- normalizePath(
  "C:/Users/User/Desktop/PigmentationAtlas",
  winslash = "/",
  mustWork = TRUE
)

results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################################
# 3. Helper functions
###############################################################################

find_coloc_file <- function(results_directory) {

  # First try the most likely final output locations.
  candidate_paths <- c(
    file.path(
      results_directory,
      "module12",
      "module12_4_all_colocalization_results.tsv.gz"
    ),
    file.path(
      results_directory,
      "module12",
      "module12_4_all_colocalization_results.tsv"
    ),
    file.path(
      results_directory,
      "module12",
      "colocalization",
      "module12_4_all_colocalization_results.tsv.gz"
    ),
    file.path(
      results_directory,
      "module12",
      "colocalization",
      "module12_4_all_colocalization_results.tsv"
    ),
    file.path(
      results_directory,
      "module12",
      "colocalization",
      "all_colocalization_results.tsv.gz"
    ),
    file.path(
      results_directory,
      "module12",
      "colocalization",
      "all_colocalization_results.tsv"
    ),
    file.path(
      results_directory,
      "all_colocalization_results.tsv.gz"
    ),
    file.path(
      results_directory,
      "all_colocalization_results.tsv"
    )
  )

  existing_candidates <- candidate_paths[file.exists(candidate_paths)]

  if (length(existing_candidates) > 0) {
    selected <- existing_candidates[1]

    message(
      "Selected colocalization file: ",
      normalizePath(selected, winslash = "/")
    )

    return(selected)
  }

  # Recursive fallback search.
  recursive_hits <- list.files(
    path = results_directory,
    pattern = "(all.*coloc|coloc.*result|posterior).*\\.(tsv|txt)(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  # Exclude files that are clearly summaries, figures, annotations or
  # gene-level collapsed outputs.
  recursive_hits <- recursive_hits[
    !grepl(
      paste0(
        "summary|figure|plot|legend|annotation|biological_annotation|",
        "gene_level|candidate_gene|reclassification|top_"
      ),
      recursive_hits,
      ignore.case = TRUE
    )
  ]

  if (length(recursive_hits) == 0) {

    stop(
      paste(
        "No complete colocalization result file was found.",
        "Run the PowerShell search command printed below the script",
        "and insert the exact file path into coloc_file."
      )
    )
  }

  file_info <- file.info(recursive_hits)

  # Prefer the largest matching file, because the complete result table
  # should normally be larger than intermediate summaries.
  recursive_hits <- recursive_hits[
    order(file_info$size, decreasing = TRUE)
  ]

  selected <- recursive_hits[1]

  message(
    "Recursively selected colocalization file: ",
    normalizePath(selected, winslash = "/")
  )

  selected
}


first_existing_column <- function(data, candidates) {

  hit <- candidates[candidates %in% names(data)]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}


convert_numeric <- function(x) {

  if (is.numeric(x)) {
    return(x)
  }

  suppressWarnings(
    as.numeric(
      gsub(",", ".", as.character(x), fixed = TRUE)
    )
  )
}


clean_class_name <- function(x) {

  x <- toupper(trimws(as.character(x)))
  x <- gsub("[ .-]+", "_", x)

  x[x %in% c(
    "STRONG",
    "STRONGH4",
    "STRONG_H_4",
    "H4_STRONG"
  )] <- "STRONG_H4"

  x[x %in% c(
    "MODERATE",
    "MODERATEH4",
    "MODERATE_H_4",
    "H4_MODERATE"
  )] <- "MODERATE_H4"

  x[x %in% c(
    "SUGGESTIVE",
    "SUGGESTIVEH4",
    "SUGGESTIVE_H_4",
    "H4_SUGGESTIVE"
  )] <- "SUGGESTIVE_H4"

  x[x %in% c(
    "H3",
    "H3DOMINANT",
    "H3_DOMINANT_SIGNAL",
    "DISTINCT_CAUSAL_VARIANTS"
  )] <- "H3_DOMINANT"

  x[x %in% c(
    "NONE",
    "OTHER",
    "UNRESOLVED",
    "NO_EVIDENCE"
  )] <- "INCONCLUSIVE"

  x
}


theme_pigmentation_atlas <- function(base_size = 10) {

  theme_classic(base_size = base_size) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = base_size + 1,
        hjust = 0
      ),

      plot.subtitle = element_text(
        size = base_size - 1,
        margin = margin(b = 7)
      ),

      axis.title = element_text(
        face = "bold"
      ),

      axis.text = element_text(
        colour = "black"
      ),

      strip.background = element_rect(
        fill = "grey92",
        colour = "black",
        linewidth = 0.4
      ),

      strip.text = element_text(
        face = "bold",
        colour = "black"
      ),

      legend.title = element_text(
        face = "bold"
      ),

      legend.position = "right",

      plot.margin = margin(
        t = 8,
        r = 10,
        b = 8,
        l = 8
      )
    )
}

###############################################################################
# 4. Locate and read the complete colocalization result table
###############################################################################

coloc_file <- file.path(
  project_dir,
  "results",
  "module12",
  "colocalization",
  "catalogue",
  "summary",
  "module12_3_coloc_catalogue_summary.tsv.gz"
)

if (!file.exists(coloc_file)) {
  stop("Colocalization master file not found: ", coloc_file)
}

message(
  "Using colocalization master file: ",
  normalizePath(coloc_file, winslash = "/")
)

coloc <- fread(
  coloc_file,
  showProgress = FALSE
)

coloc <- fread(
  coloc_file,
  na.strings = c("", "NA", "NaN", "NULL", "."),
  showProgress = TRUE
)

cat("\nColocalization table dimensions:\n")
print(dim(coloc))

cat("\nDetected columns:\n")
print(names(coloc))

###############################################################################
# 5. Detect posterior probability columns
###############################################################################

pp_h0_col <- first_existing_column(
  coloc,
  c(
    "PP.H0.abf",
    "PP.H0",
    "PPH0",
    "pp_h0",
    "posterior_h0"
  )
)

pp_h1_col <- first_existing_column(
  coloc,
  c(
    "PP.H1.abf",
    "PP.H1",
    "PPH1",
    "pp_h1",
    "posterior_h1"
  )
)

pp_h2_col <- first_existing_column(
  coloc,
  c(
    "PP.H2.abf",
    "PP.H2",
    "PPH2",
    "pp_h2",
    "posterior_h2"
  )
)

pp_h3_col <- first_existing_column(
  coloc,
  c(
    "PP.H3.abf",
    "PP.H3",
    "PPH3",
    "pp_h3",
    "posterior_h3"
  )
)

pp_h4_col <- first_existing_column(
  coloc,
  c(
    "PP.H4.abf",
    "PP.H4",
    "PPH4",
    "pp_h4",
    "posterior_h4"
  )
)

detected_pp_columns <- c(
  H0 = pp_h0_col,
  H1 = pp_h1_col,
  H2 = pp_h2_col,
  H3 = pp_h3_col,
  H4 = pp_h4_col
)

if (any(is.na(detected_pp_columns))) {

  missing_hypotheses <- names(
    detected_pp_columns
  )[is.na(detected_pp_columns)]

  stop(
    paste0(
      "The following posterior probability columns were not detected: ",
      paste(missing_hypotheses, collapse = ", "),
      ". Check the printed column names and extend the candidate list."
    )
  )
}

cat("\nDetected posterior columns:\n")
print(detected_pp_columns)

###############################################################################
# 6. Standardize posterior probabilities
###############################################################################

coloc[, PP_H0 := convert_numeric(get(pp_h0_col))]
coloc[, PP_H1 := convert_numeric(get(pp_h1_col))]
coloc[, PP_H2 := convert_numeric(get(pp_h2_col))]
coloc[, PP_H3 := convert_numeric(get(pp_h3_col))]
coloc[, PP_H4 := convert_numeric(get(pp_h4_col))]

posterior_columns <- c(
  "PP_H0",
  "PP_H1",
  "PP_H2",
  "PP_H3",
  "PP_H4"
)

# Remove rows without valid H3/H4 probabilities.
coloc_plot <- coloc[
  is.finite(PP_H3) &
    is.finite(PP_H4)
]

if (nrow(coloc_plot) == 0) {
  stop("No valid PP.H3/PP.H4 values were found.")
}

# Check probability limits.
for (column_name in posterior_columns) {

  invalid_n <- sum(
    coloc_plot[[column_name]] < 0 |
      coloc_plot[[column_name]] > 1,
    na.rm = TRUE
  )

  if (invalid_n > 0) {
    stop(
      column_name,
      " contains ",
      invalid_n,
      " values outside the interval [0, 1]."
    )
  }
}

###############################################################################
# 7. Posterior sum quality control
###############################################################################

coloc_plot[, posterior_sum :=
  PP_H0 + PP_H1 + PP_H2 + PP_H3 + PP_H4
]

posterior_sum_deviation <- abs(
  coloc_plot$posterior_sum - 1
)

cat("\nPosterior probability sum QC:\n")
cat(
  "Median sum:",
  median(coloc_plot$posterior_sum, na.rm = TRUE),
  "\n"
)
cat(
  "Maximum absolute deviation from 1:",
  max(posterior_sum_deviation, na.rm = TRUE),
  "\n"
)

if (
  max(posterior_sum_deviation, na.rm = TRUE) > 0.01
) {

  warning(
    paste(
      "Some posterior probability sums deviate from 1 by more than 0.01.",
      "Check whether values were rounded in the source table."
    )
  )
}

###############################################################################
# 8. Determine evidence classes
#
# Threshold scheme used in the PigmentationAtlas analysis:
# Strong H4:     PP.H4 >= 0.80
# Moderate H4:   0.50 <= PP.H4 < 0.80
# Suggestive H4: 0.20 <= PP.H4 < 0.50 and PP.H4 > PP.H3
# H3 dominant:   PP.H3 >= 0.20 and PP.H3 > PP.H4
# Inconclusive:  all remaining analyses
###############################################################################

class_col <- first_existing_column(
  coloc_plot,
  c(
    "priority_class",
    "evidence_class",
    "coloc_class",
    "classification",
    "final_class",
    "posterior_class"
  )
)

if (!is.na(class_col)) {

  coloc_plot[, evidence_class := clean_class_name(
    get(class_col)
  )]

  allowed_classes <- c(
    "STRONG_H4",
    "MODERATE_H4",
    "SUGGESTIVE_H4",
    "H3_DOMINANT",
    "INCONCLUSIVE"
  )

  unknown_classes <- setdiff(
    unique(coloc_plot$evidence_class),
    allowed_classes
  )

  unknown_classes <- unknown_classes[
    !is.na(unknown_classes)
  ]

  if (length(unknown_classes) > 0) {

    message(
      "Unrecognized source classes were reclassified from PP.H3/PP.H4: ",
      paste(unknown_classes, collapse = ", ")
    )

    coloc_plot[, evidence_class := NA_character_]
  }
}

if (
  is.na(class_col) ||
    all(is.na(coloc_plot$evidence_class))
) {

  coloc_plot[
    ,
    evidence_class := fifelse(
      PP_H4 >= 0.80,
      "STRONG_H4",

      fifelse(
        PP_H4 >= 0.50,
        "MODERATE_H4",

        fifelse(
          PP_H4 >= 0.20 & PP_H4 > PP_H3,
          "SUGGESTIVE_H4",

          fifelse(
            PP_H3 >= 0.20 & PP_H3 > PP_H4,
            "H3_DOMINANT",
            "INCONCLUSIVE"
          )
        )
      )
    )
  ]
}

class_levels <- c(
  "STRONG_H4",
  "MODERATE_H4",
  "SUGGESTIVE_H4",
  "H3_DOMINANT",
  "INCONCLUSIVE"
)

class_labels <- c(
  "Strong H4",
  "Moderate H4",
  "Suggestive H4",
  "H3 dominant",
  "Inconclusive"
)

coloc_plot[, evidence_class := factor(
  evidence_class,
  levels = class_levels,
  labels = class_labels
)]

###############################################################################
# 9. Class counts and consistency check
###############################################################################

class_counts <- coloc_plot[
  ,
  .(n = .N),
  by = evidence_class
]

class_counts[, percentage :=
  100 * n / sum(n)
]

cat("\nEvidence-class counts:\n")
print(class_counts)

expected_total <- 7840L

if (nrow(coloc_plot) != expected_total) {

  warning(
    paste0(
      "The input contains ",
      format(nrow(coloc_plot), big.mark = ","),
      " valid analyses rather than the expected ",
      format(expected_total, big.mark = ","),
      ". Confirm that the complete final colocalization table was selected."
    )
  )
}

###############################################################################
# 10. Data for Panel A
###############################################################################

posterior_long <- melt(
  coloc_plot,
  measure.vars = posterior_columns,
  variable.name = "hypothesis",
  value.name = "posterior_probability"
)

posterior_long[, hypothesis := factor(
  hypothesis,
  levels = posterior_columns,
  labels = c(
    "H0: no association",
    "H1: GWAS only",
    "H2: molecular QTL only",
    "H3: distinct causal variants",
    "H4: shared causal variant"
  )
)]

###############################################################################
# 11. Data for Panel C
###############################################################################

class_summary <- coloc_plot[
  ,
  .(
    n = .N,
    median_PP_H4 = median(PP_H4, na.rm = TRUE),
    q1_PP_H4 = quantile(PP_H4, 0.25, na.rm = TRUE),
    q3_PP_H4 = quantile(PP_H4, 0.75, na.rm = TRUE),
    minimum_PP_H4 = min(PP_H4, na.rm = TRUE),
    maximum_PP_H4 = max(PP_H4, na.rm = TRUE)
  ),
  by = evidence_class
]

###############################################################################
# 12. Data for Panel D
###############################################################################

ranked_h4 <- copy(coloc_plot)

setorder(
  ranked_h4,
  -PP_H4
)

ranked_h4[, rank := .I]

ranked_h4[, rank_fraction :=
  rank / .N
]

###############################################################################
# 13. Plot A: PP.H0–PP.H4 distributions
###############################################################################

plot_A <- ggplot(
  posterior_long,
  aes(
    x = posterior_probability
  )
) +

  geom_histogram(
    bins = 40,
    boundary = 0,
    closed = "left",
    colour = "white",
    linewidth = 0.15
  ) +

  facet_wrap(
    ~ hypothesis,
    ncol = 1,
    scales = "free_y"
  ) +

  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.01))
  ) +

  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.08))
  ) +

  labs(
    title = "A  Posterior probability distributions",
    subtitle = "Complete distributions of the five coloc hypotheses",
    x = "Posterior probability",
    y = "Number of analyses"
  ) +

  theme_pigmentation_atlas(base_size = 9) +

  theme(
    strip.text = element_text(
      size = 8.2
    )
  )

###############################################################################
# 14. Plot B: PP.H3 versus PP.H4
###############################################################################

plot_B <- ggplot(
  coloc_plot,
  aes(
    x = PP_H3,
    y = PP_H4,
    shape = evidence_class
  )
) +

  geom_point(
    alpha = 0.38,
    size = 1.15,
    stroke = 0.15
  ) +

  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +

  geom_hline(
    yintercept = c(0.20, 0.50, 0.80),
    linetype = "dotted",
    linewidth = 0.35
  ) +

  geom_vline(
    xintercept = 0.20,
    linetype = "dotted",
    linewidth = 0.35
  ) +

  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.01, 0.01))
  ) +

  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.01, 0.01))
  ) +

  labs(
    title = "B  Shared versus distinct causal-variant evidence",
    subtitle = "Each point represents one GWAS–molecular QTL comparison",
    x = "PP.H3: distinct causal variants",
    y = "PP.H4: shared causal variant",
    shape = "Evidence class"
  ) +

  theme_pigmentation_atlas() +

  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  )

###############################################################################
# 15. Plot C: PP.H4 by evidence class
###############################################################################

plot_C <- ggplot(
  coloc_plot,
  aes(
    x = evidence_class,
    y = PP_H4
  )
) +

  geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.35
  ) +

  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.4
  ) +

  geom_hline(
    yintercept = c(0.20, 0.50, 0.80),
    linetype = "dotted",
    linewidth = 0.35
  ) +

  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.02))
  ) +

  labs(
    title = "C  PP.H4 by evidence class",
    subtitle = "Distributions reflect the classification thresholds and H3/H4 dominance",
    x = NULL,
    y = "PP.H4"
  ) +

  theme_pigmentation_atlas() +

  theme(
    axis.text.x = element_text(
      angle = 28,
      hjust = 1,
      vjust = 1
    ),
    legend.position = "none"
  )

###############################################################################
# 16. Plot D: ranked PP.H4 profile
###############################################################################

plot_D <- ggplot(
  ranked_h4,
  aes(
    x = rank_fraction,
    y = PP_H4
  )
) +

  geom_line(
    linewidth = 0.65
  ) +

  geom_hline(
    yintercept = c(0.20, 0.50, 0.80),
    linetype = "dotted",
    linewidth = 0.4
  ) +

  annotate(
    "text",
    x = 0.98,
    y = 0.815,
    label = "Strong H4 threshold",
    hjust = 1,
    vjust = 0,
    size = 2.8
  ) +

  annotate(
    "text",
    x = 0.98,
    y = 0.515,
    label = "Moderate H4 threshold",
    hjust = 1,
    vjust = 0,
    size = 2.8
  ) +

  annotate(
    "text",
    x = 0.98,
    y = 0.215,
    label = "Suggestive H4 threshold",
    hjust = 1,
    vjust = 0,
    size = 2.8
  ) +

  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.01))
  ) +

  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.01))
  ) +

  labs(
    title = "D  Ranked PP.H4 profile",
    subtitle = "Analyses ordered from the highest to the lowest PP.H4",
    x = "Fraction of ranked colocalization analyses",
    y = "PP.H4"
  ) +

  theme_pigmentation_atlas()

###############################################################################
# 17. Combine panels
###############################################################################

figure_S2 <- (
  plot_A | plot_B
) / (
  plot_C | plot_D
) +

  plot_layout(
    widths = c(0.92, 1.08),
    heights = c(1.08, 0.92),
    guides = "collect"
  ) +

  plot_annotation(
    title = paste(
      "Extended posterior probability analyses across",
      "PigmentationAtlas colocalization tests"
    ),

    subtitle = paste0(
      "Complete posterior distributions and H3–H4 relationships across ",
      format(nrow(coloc_plot), big.mark = ","),
      " GWAS–molecular QTL comparisons"
    ),

    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 15,
        hjust = 0
      ),

      plot.subtitle = element_text(
        size = 10.5,
        hjust = 0,
        margin = margin(b = 10)
      ),

      plot.margin = margin(
        t = 12,
        r = 12,
        b = 10,
        l = 12
      )
    )
  ) &

  theme(
    legend.position = "bottom"
  )

###############################################################################
# 18. Display
###############################################################################

print(figure_S2)

###############################################################################
# 19. Export figure
###############################################################################

pdf_file <- file.path(
  figures_dir,
  "Supplementary_Figure_S2_posterior_probabilities.pdf"
)

png_file <- file.path(
  figures_dir,
  "Supplementary_Figure_S2_posterior_probabilities.png"
)

tiff_file <- file.path(
  figures_dir,
  "Supplementary_Figure_S2_posterior_probabilities.tiff"
)

ggsave(
  filename = pdf_file,
  plot = figure_S2,
  width = 13,
  height = 10,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = png_file,
  plot = figure_S2,
  width = 13,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = tiff_file,
  plot = figure_S2,
  width = 13,
  height = 10,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

###############################################################################
# 20. Export summary statistics
###############################################################################

posterior_summary <- rbindlist(
  lapply(
    posterior_columns,
    function(column_name) {

      values <- coloc_plot[[column_name]]

      data.table(
        section = "posterior_probability",
        variable = column_name,
        category = NA_character_,
        n = sum(!is.na(values)),
        mean = mean(values, na.rm = TRUE),
        median = median(values, na.rm = TRUE),
        q1 = quantile(values, 0.25, na.rm = TRUE),
        q3 = quantile(values, 0.75, na.rm = TRUE),
        minimum = min(values, na.rm = TRUE),
        maximum = max(values, na.rm = TRUE)
      )
    }
  )
)

class_summary_export <- class_summary[
  ,
  .(
    section = "evidence_class",
    variable = "PP_H4",
    category = as.character(evidence_class),
    n,
    mean = NA_real_,
    median = median_PP_H4,
    q1 = q1_PP_H4,
    q3 = q3_PP_H4,
    minimum = minimum_PP_H4,
    maximum = maximum_PP_H4
  )
]

summary_export <- rbindlist(
  list(
    posterior_summary,
    class_summary_export
  ),
  use.names = TRUE,
  fill = TRUE
)

summary_file <- file.path(
  results_dir,
  "Supplementary_Figure_S2_summary_statistics.tsv"
)

fwrite(
  summary_export,
  summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

###############################################################################
# 21. Export class counts
###############################################################################

class_count_file <- file.path(
  results_dir,
  "Supplementary_Figure_S2_evidence_class_counts.tsv"
)

fwrite(
  class_counts,
  class_count_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

###############################################################################
# 22. Console summary
###############################################################################

cat("\n")
cat("============================================================\n")
cat("Supplementary Figure S2 completed successfully\n")
cat("============================================================\n")
cat("Input:   ", normalizePath(coloc_file, winslash = "/"), "\n")
cat("Rows:    ", format(nrow(coloc_plot), big.mark = ","), "\n")
cat("PDF:     ", pdf_file, "\n")
cat("PNG:     ", png_file, "\n")
cat("TIFF:    ", tiff_file, "\n")
cat("Summary: ", summary_file, "\n")
cat("Classes: ", class_count_file, "\n")
cat("============================================================\n\n")

print(class_counts)