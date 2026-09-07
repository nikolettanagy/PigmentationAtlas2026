###############################################################################
# PigmentationAtlas
# Figure 1 – Analytical workflow
###############################################################################

options(stringsAsFactors = FALSE)
options(scipen = 999)

required_packages <- c("ggplot2", "patchwork", "data.table")

install_missing_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0L) {
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

install_missing_packages(required_packages)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

project_dir <- normalizePath(
  "C:/Users/User/Desktop/PigmentationAtlas",
  winslash = "/",
  mustWork = TRUE
)

figure_dir <- file.path(project_dir, "figures")
figure_data_dir <- file.path(figure_dir, "figure_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_data_dir, recursive = TRUE, showWarnings = FALSE)

output_pdf <- file.path(figure_dir, "Figure1_PigmentationAtlas_Workflow.pdf")
output_png <- file.path(figure_dir, "Figure1_PigmentationAtlas_Workflow.png")
output_tiff <- file.path(figure_dir, "Figure1_PigmentationAtlas_Workflow.tiff")

col_input <- "#D9D9D9"
col_gwas <- "#4C78A8"
col_locus <- "#72B7B2"
col_inference <- "#54A24B"
col_priority <- "#F2A541"
col_output <- "#E45756"
col_outline <- "#4A4A4A"
col_arrow <- "#333333"

add_box <- function(p, xmin, xmax, ymin, ymax, label, fill,
                    text_size = 3.5, fontface = "plain",
                    border = col_outline, linewidth = 0.45) {
  p +
    annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             fill = fill, colour = border, linewidth = linewidth) +
    annotate("text", x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
             label = label, size = text_size, fontface = fontface,
             lineheight = 0.95)
}

add_arrow <- function(p, x, xend, y, yend, linewidth = 0.55) {
  p + annotate(
    "segment", x = x, xend = xend, y = y, yend = yend,
    colour = col_arrow, linewidth = linewidth,
    arrow = arrow(length = unit(0.12, "inches"), type = "closed")
  )
}

add_stage_band <- function(p, xmin, xmax, label, fill) {
  p +
    annotate("rect", xmin = xmin, xmax = xmax, ymin = 0.20, ymax = 0.62,
             fill = fill, colour = NA) +
    annotate("text", x = (xmin + xmax) / 2, y = 0.41,
             label = label, size = 3.2, fontface = "bold")
}

p <- ggplot() +
  coord_cartesian(xlim = c(0, 16), ylim = c(0, 10), clip = "off") +
  theme_void(base_size = 10) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(8, 12, 8, 12))

p <- p +
  annotate("text", x = 1.6, y = 9.25, label = "INPUTS", fontface = "bold", size = 4.2) +
  annotate("text", x = 5.0, y = 9.25, label = "GWAS PROCESSING", fontface = "bold", size = 4.2) +
  annotate("text", x = 9.0, y = 9.25, label = "STATISTICAL INFERENCE", fontface = "bold", size = 4.2) +
  annotate("text", x = 14.0, y = 9.25, label = "BIOLOGICAL INTERPRETATION", fontface = "bold", size = 4.2)

p <- add_box(p, 0.45, 2.75, 7.35, 8.35,
             "Pigmentation GWAS\nsummary statistics", col_input, 3.6, "bold")
p <- add_box(p, 0.45, 2.75, 5.80, 6.80,
             "QTL resources\nGTEx and eQTL datasets", col_input, 3.4)
p <- add_box(p, 0.45, 2.75, 4.25, 5.25,
             "Reference resources\n1000G EUR • GENCODE v50", col_input, 3.25)

p <- add_box(p, 3.45, 5.60, 7.35, 8.35,
             "Harmonization\nand quality control", col_gwas, 3.35, "bold")
p <- add_box(p, 3.45, 5.60, 5.80, 6.80,
             "LD clumping\nindependent signals", col_gwas, 3.35)
p <- add_box(p, 3.45, 5.60, 4.25, 5.25,
             "Genomic-locus\nconstruction", col_locus, 3.35)
p <- add_box(p, 3.45, 5.60, 2.70, 3.70,
             "Clump workspaces\nand candidate genes", col_locus, 3.25)

p <- add_box(p, 6.65, 8.85, 7.35, 8.35,
             "GWAS–QTL\nvariant matching", col_inference, 3.35, "bold")
p <- add_box(p, 9.10, 11.30, 7.35, 8.35,
             "Bayesian\ncolocalization", col_inference, 3.35, "bold")
p <- add_box(p, 6.65, 8.85, 5.35, 6.35,
             "Posterior probabilities\nPP.H3 and PP.H4", col_inference, 3.25)
p <- add_box(p, 9.10, 11.30, 5.35, 6.35,
             "Evidence classes\nStrong • Moderate • Suggestive", col_priority, 3.05)
p <- add_box(p, 6.65, 8.85, 3.35, 4.35,
             "Biological annotation\nand record-level catalogue", col_priority, 3.05)
p <- add_box(p, 9.10, 11.30, 3.35, 4.35,
             "Gene-level aggregation\none clump × one gene", col_priority, 3.05)

p <- add_box(p, 12.20, 15.45, 7.35, 8.35,
             "Best H4-supported\ncandidate gene per clump", col_output, 3.35, "bold")
p <- add_box(p, 12.20, 15.45, 5.55, 6.55,
             "Conventional assignment\nnearest or overlapping coding gene", "#F5B7B1", 3.10)
p <- add_box(p, 12.20, 15.45, 3.75, 4.75,
             "Colocalization-driven\ngene reclassification", col_output, 3.35, "bold")
p <- add_box(p, 12.20, 15.45, 1.95, 2.95,
             "PigmentationAtlas\ncandidate-gene landscape", "#F8D7DA", 3.35, "bold")

p <- add_arrow(p, 2.75, 3.45, 7.85, 7.85)
p <- add_arrow(p, 2.75, 3.45, 6.30, 6.30)
p <- add_arrow(p, 2.75, 3.45, 4.75, 4.75)
p <- add_arrow(p, 4.53, 4.53, 7.35, 6.80)
p <- add_arrow(p, 4.53, 4.53, 5.80, 5.25)
p <- add_arrow(p, 4.53, 4.53, 4.25, 3.70)
p <- add_arrow(p, 5.60, 6.65, 3.20, 7.85)
p <- add_arrow(p, 8.85, 9.10, 7.85, 7.85)
p <- add_arrow(p, 7.75, 7.75, 7.35, 6.35)
p <- add_arrow(p, 10.20, 10.20, 7.35, 6.35)
p <- add_arrow(p, 7.75, 7.75, 5.35, 4.35)
p <- add_arrow(p, 10.20, 10.20, 5.35, 4.35)
p <- add_arrow(p, 8.85, 9.10, 3.85, 3.85)
p <- add_arrow(p, 11.30, 12.20, 7.85, 7.85)
p <- add_arrow(p, 11.30, 12.20, 3.85, 4.25)
p <- add_arrow(p, 13.83, 13.83, 7.35, 6.55)
p <- add_arrow(p, 13.83, 13.83, 5.55, 4.75)
p <- add_arrow(p, 13.83, 13.83, 3.75, 2.95)

p <- add_stage_band(p, 0.35, 2.85, "DATA INPUT", "#ECECEC")
p <- add_stage_band(p, 3.15, 5.90, "SIGNAL DEFINITION", "#DCE6F1")
p <- add_stage_band(p, 6.30, 11.65, "COLOCALIZATION AND PRIORITIZATION", "#DDEEDD")
p <- add_stage_band(p, 12.00, 15.65, "GENE REINTERPRETATION", "#F9D6D5")

p <- p + annotate(
  "text", x = 8.0, y = 1.10,
  label = "Primary analytical unit: independent PLINK clump; genomic loci retained for annotation and visualization",
  size = 3.0, fontface = "italic"
)

figure1 <- p + plot_annotation(
  title = "PigmentationAtlas analytical workflow",
  subtitle = "From genome-wide pigmentation GWAS signals to colocalization-supported candidate-gene prioritization and locus reinterpretation",
  theme = theme(
    plot.title = element_text(face = "bold", size = 19, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10.5, margin = margin(b = 10))
  )
)

ggsave(output_pdf, figure1, width = 16, height = 9.5, units = "in", device = cairo_pdf)
ggsave(output_png, figure1, width = 16, height = 9.5, units = "in", dpi = 300, bg = "white")
ggsave(output_tiff, figure1, width = 16, height = 9.5, units = "in", dpi = 600,
       compression = "lzw", bg = "white")

status <- data.table(
  figure = "Figure 1",
  analysis = "PigmentationAtlas analytical workflow",
  status = "SUCCESS",
  pdf_file = output_pdf,
  png_file = output_png,
  tiff_file = output_tiff,
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

fwrite(status, file.path(figure_data_dir, "Figure1_status.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(figure_data_dir, "Figure1_sessionInfo.txt"))

cat("\n", paste(rep("=", 78), collapse = ""), "\n", sep = "")
cat("FIGURE 1 COMPLETED\n")
cat(paste(rep("=", 78), collapse = ""), "\n", sep = "")
cat("PDF:  ", output_pdf, "\n")
cat("PNG:  ", output_png, "\n")
cat("TIFF: ", output_tiff, "\n")
