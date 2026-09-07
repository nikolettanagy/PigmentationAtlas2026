###############################################################################
# PigmentationAtlas
# Supplementary Figure S1
#
# Quality control and analytical summary of the PigmentationAtlas pipeline
#
# Output:
#   figures/Supplementary_Figure_S1_pipeline_audit.pdf
#   figures/Supplementary_Figure_S1_pipeline_audit.png
#   figures/Supplementary_Figure_S1_pipeline_audit.tiff
#
# R version: 4.6.0
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
  "ggplot2",
  "scales"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(ggplot2)
library(scales)

###############################################################################
# 2. Project directories
###############################################################################

project_dir <- normalizePath(
  "C:/Users/User/Desktop/PigmentationAtlas",
  winslash = "/",
  mustWork = TRUE
)

figures_dir <- file.path(project_dir, "figures")

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################################
# 3. Validated PigmentationAtlas summary values
#
# IMPORTANT:
# Check these values against the final Supplementary Tables before submission.
###############################################################################

n_harmonized_variants <- 21831407
n_gws_variants        <- 33455

n_ld_clumps           <- 1325
n_genomic_loci        <- 194
n_empty_loci          <- 0
n_missing_lead_snps   <- 0

n_candidate_genes     <- 1603

n_coloc_analyses      <- 7840

n_strong_h4           <- 25
n_moderate_h4         <- 51
n_suggestive_h4       <- 207
n_h3_dominant         <- 281
n_inconclusive        <- 7276

n_h4_supported_clumps <- 113
n_reassigned          <- 97
n_unchanged           <- 16

reassignment_percent <- 100 * n_reassigned / n_h4_supported_clumps

###############################################################################
# 4. Internal consistency checks
###############################################################################

stopifnot(
  n_strong_h4 +
    n_moderate_h4 +
    n_suggestive_h4 +
    n_h3_dominant +
    n_inconclusive ==
    n_coloc_analyses
)

stopifnot(
  n_reassigned + n_unchanged ==
    n_h4_supported_clumps
)

stopifnot(
  n_empty_loci == 0,
  n_missing_lead_snps == 0
)

###############################################################################
# 5. Helper functions
###############################################################################

format_count <- function(x) {
  comma(x, accuracy = 1)
}

format_million <- function(x) {
  paste0(
    format(
      round(x / 1e6, 1),
      nsmall = 1,
      trim = TRUE
    ),
    " million"
  )
}

###############################################################################
# 6. Pipeline audit data
###############################################################################

audit <- data.frame(
  step = 1:7,

  stage = c(
    "GWAS harmonization",
    "Genome-wide signal selection",
    "LD clumping",
    "Genomic locus construction",
    "Candidate-gene annotation",
    "Molecular QTL colocalization",
    "Effector-gene prioritization"
  ),

  main_value = c(
    format_million(n_harmonized_variants),
    format_count(n_gws_variants),
    format_count(n_ld_clumps),
    format_count(n_genomic_loci),
    format_count(n_candidate_genes),
    format_count(n_coloc_analyses),
    format_count(n_h4_supported_clumps)
  ),

  main_label = c(
    "harmonized variants",
    "genome-wide significant variants",
    "independent PLINK clumps",
    "merged genomic loci",
    "candidate genes",
    "completed colocalization analyses",
    "H4-supported evaluable clumps"
  ),

  qc_text = c(
    "Variant identifiers, alleles, positions and association statistics harmonized",
    "Variants retained at P < 5 × 10\u207b\u2078",
    "Independent association signals defined using the EUR LD reference panel",
    paste0(
      "Complete locus construction: ",
      n_empty_loci,
      " empty loci; ",
      n_missing_lead_snps,
      " missing lead variants"
    ),
    "GENCODE v50 gene annotation and biotype assignment",
    paste0(
      "H4 evidence: ",
      n_strong_h4,
      " strong, ",
      n_moderate_h4,
      " moderate and ",
      n_suggestive_h4,
      " suggestive analyses"
    ),
    paste0(
      n_reassigned,
      " reassigned; ",
      n_unchanged,
      " unchanged (",
      sprintf("%.1f", reassignment_percent),
      "% reassignment)"
    )
  ),

  status = rep("PASSED", 7),

  stringsAsFactors = FALSE
)

# Vertical positions: first analytical step at the top.
audit$y <- rev(seq_len(nrow(audit)))

###############################################################################
# 7. Layout parameters
###############################################################################

box_xmin <- 0.8
box_xmax <- 9.2

status_xmin <- 7.85
status_xmax <- 8.90

box_height <- 0.72

audit$ymin <- audit$y - box_height / 2
audit$ymax <- audit$y + box_height / 2

###############################################################################
# 8. Arrow data
###############################################################################

arrows <- data.frame(
  x = rep(5.0, nrow(audit) - 1),
  xend = rep(5.0, nrow(audit) - 1),

  y = audit$y[-nrow(audit)] - box_height / 2,
  yend = audit$y[-1] + box_height / 2
)

###############################################################################
# 9. Final outcome summary
###############################################################################

outcome <- data.frame(
  x = c(3.15, 5.00, 6.85),

  value = c(
    n_reassigned,
    n_unchanged,
    reassignment_percent
  ),

  value_text = c(
    format_count(n_reassigned),
    format_count(n_unchanged),
    paste0(sprintf("%.1f", reassignment_percent), "%")
  ),

  label = c(
    "Reassigned clumps",
    "Unchanged clumps",
    "Reassignment rate"
  )
)

###############################################################################
# 10. Create the figure
###############################################################################

figure_S1 <- ggplot() +

  # Connecting arrows
  geom_segment(
    data = arrows,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    linewidth = 0.55,
    arrow = arrow(
      length = unit(0.16, "cm"),
      type = "closed"
    )
  ) +

  # Main stage boxes
  geom_rect(
    data = audit,
    aes(
      xmin = box_xmin,
      xmax = box_xmax,
      ymin = ymin,
      ymax = ymax
    ),
    fill = "white",
    colour = "black",
    linewidth = 0.55
  ) +

  # Left analytical-stage band
  geom_rect(
    data = audit,
    aes(
      xmin = box_xmin,
      xmax = 2.65,
      ymin = ymin,
      ymax = ymax
    ),
    fill = "grey92",
    colour = "black",
    linewidth = 0.45
  ) +

  # QC status boxes
  geom_rect(
    data = audit,
    aes(
      xmin = status_xmin,
      xmax = status_xmax,
      ymin = y - 0.18,
      ymax = y + 0.18
    ),
    fill = "grey88",
    colour = "black",
    linewidth = 0.4
  ) +

  # Stage names
  geom_text(
    data = audit,
    aes(
      x = 1.72,
      y = y,
      label = stage
    ),
    fontface = "bold",
    size = 3.25,
    lineheight = 0.95
  ) +

  # Main numerical values
  geom_text(
    data = audit,
    aes(
      x = 3.05,
      y = y + 0.13,
      label = main_value
    ),
    hjust = 0,
    fontface = "bold",
    size = 4.1
  ) +

  # Main numerical labels
  geom_text(
    data = audit,
    aes(
      x = 3.05,
      y = y - 0.16,
      label = main_label
    ),
    hjust = 0,
    size = 3.0
  ) +

  # QC descriptions
  geom_text(
    data = audit,
    aes(
      x = 5.40,
      y = y,
      label = qc_text
    ),
    hjust = 0,
    size = 2.65,
    lineheight = 0.93
  ) +

  # PASS labels
  geom_text(
    data = audit,
    aes(
      x = (status_xmin + status_xmax) / 2,
      y = y,
      label = "\u2713  PASSED"
    ),
    fontface = "bold",
    size = 2.75
  ) +

  # Section heading above final outcome
  annotate(
    "text",
    x = 5,
    y = -0.22,
    label = "FINAL COLOCALIZATION-DRIVEN GENE CLASSIFICATION",
    fontface = "bold",
    size = 3.65
  ) +

  # Outcome boxes
  geom_rect(
    data = outcome,
    aes(
      xmin = x - 0.72,
      xmax = x + 0.72,
      ymin = -1.15,
      ymax = -0.42
    ),
    fill = "white",
    colour = "black",
    linewidth = 0.55
  ) +

  geom_text(
    data = outcome,
    aes(
      x = x,
      y = -0.69,
      label = value_text
    ),
    fontface = "bold",
    size = 4.5
  ) +

  geom_text(
    data = outcome,
    aes(
      x = x,
      y = -0.98,
      label = label
    ),
    size = 2.9
  ) +

  # Overall annotations
  annotate(
    "text",
    x = 0.8,
    y = 7.78,
    label = "ANALYTICAL STEP",
    hjust = 0,
    fontface = "bold",
    size = 3.0
  ) +

  annotate(
    "text",
    x = 3.05,
    y = 7.78,
    label = "PIPELINE OUTPUT",
    hjust = 0,
    fontface = "bold",
    size = 3.0
  ) +

  annotate(
    "text",
    x = 5.40,
    y = 7.78,
    label = "QUALITY-CONTROL AUDIT",
    hjust = 0,
    fontface = "bold",
    size = 3.0
  ) +

  annotate(
    "text",
    x = 8.37,
    y = 7.78,
    label = "STATUS",
    fontface = "bold",
    size = 3.0
  ) +

  coord_cartesian(
    xlim = c(0.55, 9.45),
    ylim = c(-1.35, 8.05),
    clip = "off"
  ) +

  labs(
    title = "Quality control and analytical summary of the PigmentationAtlas pipeline",
    subtitle = paste(
      "Sequential audit of variant harmonization, signal definition,",
      "gene annotation, molecular QTL integration and effector-gene prioritization"
    ),
    caption = paste0(
      "All seven analytical stages completed successfully. ",
      "Among ",
      format_count(n_h4_supported_clumps),
      " evaluable H4-supported clumps, ",
      format_count(n_reassigned),
      " (",
      sprintf("%.1f", reassignment_percent),
      "%) were assigned to a gene different from the proximity-based baseline."
    )
  ) +

  theme_void(base_size = 11) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0,
      margin = margin(b = 6)
    ),

    plot.subtitle = element_text(
      size = 10.5,
      hjust = 0,
      margin = margin(b = 17)
    ),

    plot.caption = element_text(
      size = 8.5,
      hjust = 0,
      lineheight = 1.05,
      margin = margin(t = 12)
    ),

    plot.margin = margin(
      t = 18,
      r = 18,
      b = 15,
      l = 18
    )
  )

###############################################################################
# 11. Display
###############################################################################

print(figure_S1)

###############################################################################
# 12. Export
###############################################################################

pdf_file <- file.path(
  figures_dir,
  "Supplementary_Figure_S1_pipeline_audit.pdf"
)

png_file <- file.path(
  figures_dir,
  "Supplementary_Figure_S1_pipeline_audit.png"
)

tiff_file <- file.path(
  figures_dir,
  "Supplementary_Figure_S1_pipeline_audit.tiff"
)

ggsave(
  filename = pdf_file,
  plot = figure_S1,
  width = 12,
  height = 10,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = png_file,
  plot = figure_S1,
  width = 12,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = tiff_file,
  plot = figure_S1,
  width = 12,
  height = 10,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

###############################################################################
# 13. Export the numerical audit table
###############################################################################

audit_table <- audit[
  ,
  c(
    "step",
    "stage",
    "main_value",
    "main_label",
    "qc_text",
    "status"
  )
]

audit_table_file <- file.path(
  figures_dir,
  "Supplementary_Figure_S1_pipeline_audit_values.tsv"
)

write.table(
  audit_table,
  file = audit_table_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

###############################################################################
# 14. Console summary
###############################################################################

cat("\n")
cat("============================================================\n")
cat("Supplementary Figure S1 completed successfully\n")
cat("============================================================\n")
cat("PDF:   ", pdf_file, "\n")
cat("PNG:   ", png_file, "\n")
cat("TIFF:  ", tiff_file, "\n")
cat("Values:", audit_table_file, "\n")
cat("============================================================\n\n")