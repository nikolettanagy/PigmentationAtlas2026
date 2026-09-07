#!/usr/bin/env Rscript

# PigmentationAtlas Publication Engine
# Module 13 v1.1
# Figure 3: Candidate-gene biotype landscape and prioritized regulatory RNAs

options(stringsAsFactors = FALSE)

script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0L) {
    dirname(normalizePath(
      sub("^--file=", "", file_arg[[1]]),
      winslash = "/",
      mustWork = FALSE
    ))
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
})

source(file.path(script_dir, "00_packages.R"))
source(file.path(script_dir, "01_palette.R"))
source(file.path(script_dir, "02_theme.R"))
source(file.path(script_dir, "03_helpers.R"))

project_root <- pa_project_root()
pa_ensure_dirs(project_root)

input_file <- file.path(
  project_root,
  "results",
  "module13",
  "publication_dataset.tsv.gz"
)

if (!file.exists(input_file)) {
  stop(
    "Publication-ready dataset not found:\n  ",
    input_file,
    "\nRun Module13_00_prepare_publication_data.R first.",
    call. = FALSE
  )
}

output_dir <- file.path(project_root, "figures", "main")
results_dir <- file.path(project_root, "results", "module13")

message("Reading publication-ready dataset:")
message("  ", normalizePath(input_file, winslash = "/", mustWork = TRUE))

dt <- data.table::fread(input_file)

pa_validate_columns(
  dt,
  required = c(
    "gene_id_stable",
    "gene_symbol",
    "candidate_class",
    "priority_class",
    "PP_H3",
    "PP_H4",
    "clump_id",
    "tissue",
    "phenotype_id"
  ),
  object_name = "Module 13 publication dataset"
)

dt[, candidate_group := data.table::fcase(
  candidate_class == "Protein coding", "Protein coding",
  candidate_class == "Regulatory RNA", "Regulatory RNA",
  candidate_class == "Pseudogene", "Pseudogene",
  default = "Other / unknown"
)]

candidate_levels <- c(
  "Protein coding",
  "Regulatory RNA",
  "Pseudogene",
  "Other / unknown"
)

candidate_colours <- c(
  "Protein coding" = "#4C78A8",
  "Regulatory RNA" = "#E45756",
  "Pseudogene" = "#72B7B2",
  "Other / unknown" = "#BAB0AC"
)

dt[, candidate_group := factor(
  candidate_group,
  levels = candidate_levels,
  ordered = TRUE
)]

dt[, priority_class := factor(
  priority_class,
  levels = PA_CLASS_LEVELS,
  ordered = TRUE
)]

# -------------------------------------------------------------------------
# Panel A: all candidate-gene colocalization units
# -------------------------------------------------------------------------

overall_summary <- dt[
  ,
  .(count = .N),
  by = candidate_group
][
  ,
  percent := 100 * count / sum(count)
]

overall_summary[
  ,
  candidate_order := match(
    as.character(candidate_group),
    candidate_levels
  )
]
data.table::setorder(overall_summary, candidate_order)
overall_summary[, candidate_order := NULL]

panel_a <- ggplot2::ggplot(
  overall_summary,
  ggplot2::aes(
    x = forcats::fct_rev(candidate_group),
    y = count,
    fill = candidate_group
  )
) +
  ggplot2::geom_col(width = 0.68) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf(
        "%s (%.1f%%)",
        scales::comma(count),
        percent
      )
    ),
    hjust = -0.06,
    size = 2.7,
    family = "Arial"
  ) +
  ggplot2::coord_flip(clip = "off") +
  ggplot2::scale_fill_manual(
    values = candidate_colours,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::comma,
    expand = ggplot2::expansion(mult = c(0, 0.25))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Candidate-gene colocalization units",
    subtitle = "All evaluated gene–clump–tissue combinations"
  ) +
  ggplot2::theme(
    legend.position = "none"
  )

# -------------------------------------------------------------------------
# Panel B: candidate classes among H4-supported units
# Strong + Moderate + Suggestive H4
# -------------------------------------------------------------------------

h4_levels <- c(
  "Strong H4",
  "Moderate H4",
  "Suggestive H4"
)

supported_dt <- dt[
  priority_class %in% h4_levels
]

supported_summary <- supported_dt[
  ,
  .(count = .N),
  by = .(candidate_group, priority_class)
]

supported_totals <- supported_dt[
  ,
  .(total = .N),
  by = candidate_group
]

supported_summary <- merge(
  supported_summary,
  supported_totals,
  by = "candidate_group",
  all.x = TRUE
)

supported_summary[
  ,
  percent_within_group := 100 * count / total
]

panel_b <- ggplot2::ggplot(
  supported_summary,
  ggplot2::aes(
    x = candidate_group,
    y = count,
    fill = priority_class
  )
) +
  ggplot2::geom_col(width = 0.68) +
  ggplot2::geom_text(
    data = supported_totals,
    ggplot2::aes(
      x = candidate_group,
      y = total,
      label = scales::comma(total)
    ),
    inherit.aes = FALSE,
    vjust = -0.45,
    size = 2.7,
    family = "Arial"
  ) +
  ggplot2::scale_fill_manual(
    values = PA_COLOURS[h4_levels],
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::comma,
    expand = ggplot2::expansion(mult = c(0, 0.14))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "H4-supported units",
    fill = "Evidence class",
    subtitle = "Strong, moderate, or suggestive shared-causal evidence"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 22,
      hjust = 1
    ),
    legend.position = "right"
  )

# -------------------------------------------------------------------------
# Panel C: strongest regulatory-RNA candidates
# Select the best PP.H4 result for each unique regulatory-RNA gene
# -------------------------------------------------------------------------

regulatory_best <- dt[
  candidate_group == "Regulatory RNA" &
    is.finite(PP_H4)
][
  order(-PP_H4, PP_H3),
  .SD[1],
  by = gene_id_stable
]

regulatory_best[
  ,
  display_gene := data.table::fifelse(
    !is.na(gene_symbol) &
      nzchar(gene_symbol) &
      gene_symbol != gene_id_stable,
    gene_symbol,
    gene_id_stable
  )
]

regulatory_top <- regulatory_best[
  order(-PP_H4)
][
  seq_len(min(.N, 12L))
]

regulatory_top[
  ,
  display_gene := factor(
    display_gene,
    levels = rev(display_gene)
  )
]

panel_c <- ggplot2::ggplot(
  regulatory_top,
  ggplot2::aes(
    x = PP_H4,
    y = display_gene,
    colour = priority_class
  )
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend = PP_H4,
      yend = display_gene
    ),
    colour = "grey78",
    linewidth = 0.45
  ) +
  ggplot2::geom_point(
    size = 2.6,
    alpha = 0.90
  ) +
  ggplot2::geom_vline(
    xintercept = 0.20,
    linetype = "dotted",
    linewidth = 0.38,
    colour = "grey45"
  ) +
  ggplot2::geom_vline(
    xintercept = 0.50,
    linetype = "dashed",
    linewidth = 0.42,
    colour = "grey35"
  ) +
  ggplot2::scale_colour_manual(
    values = PA_COLOURS[PA_CLASS_LEVELS],
    drop = FALSE
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, max(0.55, max(regulatory_top$PP_H4) * 1.08)),
    breaks = seq(0, 0.6, 0.1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    x = "Best shared-causal posterior probability (PP.H4)",
    y = NULL,
    colour = "Classification",
    subtitle = "Highest PP.H4 result per regulatory-RNA gene"
  ) +
  ggplot2::theme(
    legend.position = "right",
    panel.grid.major.y = ggplot2::element_blank()
  )

# -------------------------------------------------------------------------
# Assemble and export
# -------------------------------------------------------------------------

figure_3 <- (
  panel_a /
    (panel_b | panel_c)
) +
  patchwork::plot_layout(
    heights = c(0.78, 1.22),
    widths = c(1.02, 0.98),
    guides = "collect"
  ) +
  patchwork::plot_annotation(
    title = "Candidate-gene biotype landscape and prioritized regulatory RNAs",
    tag_levels = "A",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        family = "Arial",
        face = "bold",
        size = 13,
        hjust = 0
      ),
      plot.tag = ggplot2::element_text(
        family = "Arial",
        face = "bold",
        size = 12
      )
    )
  ) &
  ggplot2::theme(
    legend.position = "right"
  )

pa_save_figure(
  plot = figure_3,
  filename_stem = "Figure3_candidate_gene_biotype_landscape",
  output_dir = output_dir,
  width_mm = 180,
  height_mm = 170,
  dpi = 600
)

data.table::fwrite(
  overall_summary[, .(
    candidate_class = as.character(candidate_group),
    count,
    percent = round(percent, 4)
  )],
  file.path(
    results_dir,
    "Figure3A_candidate_class_statistics.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

data.table::fwrite(
  supported_summary[, .(
    candidate_class = as.character(candidate_group),
    priority_class = as.character(priority_class),
    count,
    total_in_candidate_class = total,
    percent_within_candidate_class = round(
      percent_within_group,
      4
    )
  )],
  file.path(
    results_dir,
    "Figure3B_h4_supported_statistics.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

data.table::fwrite(
  regulatory_top[, .(
    gene_id_stable,
    gene_symbol,
    clump_id,
    tissue,
    phenotype_id,
    PP_H3,
    PP_H4,
    priority_class = as.character(priority_class)
  )],
  file.path(
    results_dir,
    "Figure3C_top_regulatory_RNA_candidates.tsv"
  ),
  sep = "\t",
  quote = FALSE
)

caption_text <- paste0(
  "Figure 3. Candidate-gene biotype landscape and prioritized regulatory ",
  "RNAs. (A) Distribution of all evaluated gene–clump–tissue colocalization ",
  "units across broad candidate-gene classes. Bars are labelled with counts ",
  "and percentages. (B) Candidate-gene composition of units with strong, ",
  "moderate, or suggestive evidence for a shared causal signal. Numbers above ",
  "bars indicate the total number of H4-supported units in each candidate ",
  "class. (C) Regulatory-RNA genes ranked by their best observed PP.H4 value; ",
  "one highest-scoring colocalization result is shown per gene."
)

pa_write_caption(
  caption_text,
  file.path(
    results_dir,
    "Figure3_caption.txt"
  )
)

pa_log_session(
  file.path(
    project_root,
    "logs",
    "module13",
    "Figure3_sessionInfo.txt"
  )
)

message("Figure 3 completed successfully.")
