suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
})

# ============================================================
# PigmentRegAtlas
# Module 04A: Build a processed GENCODE gene reference
# ============================================================

gtf_file <- paste0(
  "data/annotations/gencode/",
  "gencode.v50.primary_assembly.annotation.gtf.gz"
)

output_file <- paste0(
  "data/annotations/processed/",
  "GENCODE_v50_genes_GRCh38.tsv.gz"
)

output_qc <- paste0(
  "results/",
  "GENCODE_v50_gene_reference_QC.tsv"
)

# ------------------------------------------------------------
# Create output directory
# ------------------------------------------------------------

dir.create(
  "data/annotations/processed",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Check input
# ------------------------------------------------------------

if (!file.exists(gtf_file)) {
  stop(
    "GENCODE GTF file not found: ",
    gtf_file
  )
}

cat("Reading GENCODE annotation:\n")
cat(gtf_file, "\n\n")

# ------------------------------------------------------------
# Import GTF
# ------------------------------------------------------------

gtf <- import(gtf_file)

cat("Total GTF records imported:", length(gtf), "\n")

# Keep gene-level records only
genes_gr <- gtf[
  mcols(gtf)$type == "gene"
]

cat("Gene records retained:", length(genes_gr), "\n")

if (length(genes_gr) == 0L) {
  stop("No gene records were detected in the GTF.")
}

# ------------------------------------------------------------
# Detect metadata column names
# ------------------------------------------------------------

metadata_columns <- names(
  mcols(genes_gr)
)

cat("\nGENCODE metadata columns:\n")
print(metadata_columns)

if (!"gene_id" %in% metadata_columns) {
  stop("gene_id column is missing from the GTF.")
}

if (!"gene_name" %in% metadata_columns) {
  stop("gene_name column is missing from the GTF.")
}

gene_type_column <- NULL

if ("gene_type" %in% metadata_columns) {
  gene_type_column <- "gene_type"
}

if (
  is.null(gene_type_column) &&
  "gene_biotype" %in% metadata_columns
) {
  gene_type_column <- "gene_biotype"
}

if (is.null(gene_type_column)) {
  stop(
    "Neither gene_type nor gene_biotype was found."
  )
}

cat(
  "\nUsing gene biotype column:",
  gene_type_column,
  "\n"
)

# ------------------------------------------------------------
# Convert to data.table
# ------------------------------------------------------------

genes <- data.table(
  gene_id_version = as.character(
    mcols(genes_gr)$gene_id
  ),
  gene_name = as.character(
    mcols(genes_gr)$gene_name
  ),
  gene_type = as.character(
    mcols(genes_gr)[[gene_type_column]]
  ),
  chr = as.character(
    seqnames(genes_gr)
  ),
  start = as.integer(
    start(genes_gr)
  ),
  end = as.integer(
    end(genes_gr)
  ),
  strand = as.character(
    strand(genes_gr)
  )
)

# Remove Ensembl version suffix
genes[, gene_id := sub(
  "\\.[0-9]+$",
  "",
  gene_id_version
)]

# Standardize chromosome naming
genes[, chr := sub(
  "^chr",
  "",
  chr,
  ignore.case = TRUE
)]

# ------------------------------------------------------------
# Define broad gene categories
# ------------------------------------------------------------

lncRNA_types <- c(
  "lncRNA",
  "lincRNA",
  "antisense",
  "sense_intronic",
  "sense_overlapping",
  "processed_transcript",
  "3prime_overlapping_ncRNA",
  "macro_lncRNA",
  "bidirectional_promoter_lncRNA",
  "non_coding",
  "TEC"
)

genes[, is_protein_coding :=
  gene_type == "protein_coding"
]

genes[, is_lncRNA :=
  gene_type %in% lncRNA_types |
    grepl(
      "lnc|linc|antisense|non_coding|processed_transcript",
      gene_type,
      ignore.case = TRUE
    )
]

genes[, broad_gene_class := fifelse(
  is_protein_coding,
  "protein_coding",
  fifelse(
    is_lncRNA,
    "lncRNA",
    "other"
  )
)]

# ------------------------------------------------------------
# Calculate TSS
# ------------------------------------------------------------

genes[, tss := fifelse(
  strand == "-",
  end,
  start
)]

# ------------------------------------------------------------
# Restrict to standard chromosomes
# ------------------------------------------------------------

standard_chromosomes <- c(
  as.character(1:22),
  "X",
  "Y",
  "MT",
  "M"
)

genes[, is_standard_chromosome :=
  chr %in% standard_chromosomes
]

# Main atlas currently uses autosomes 1–22
genes[, is_autosomal :=
  chr %in% as.character(1:22)
]

# ------------------------------------------------------------
# Remove incomplete or duplicated records
# ------------------------------------------------------------

genes <- genes[
  !is.na(gene_id) &
  gene_id != "" &
  !is.na(chr) &
  chr != "" &
  !is.na(start) &
  !is.na(end) &
  start > 0 &
  end >= start
]

setorder(
  genes,
  gene_id,
  chr,
  start,
  end
)

duplicate_gene_ids <- genes[
  duplicated(gene_id) |
  duplicated(gene_id, fromLast = TRUE)
]

if (nrow(duplicate_gene_ids) > 0L) {
  warning(
    nrow(duplicate_gene_ids),
    " records have duplicated version-free gene IDs."
  )
}

# One record per version-free Ensembl gene ID
genes <- unique(
  genes,
  by = "gene_id"
)

# Genomic ordering
genes[, chr_order := suppressWarnings(
  as.integer(chr)
)]

genes[
  chr == "X",
  chr_order := 23L
]

genes[
  chr == "Y",
  chr_order := 24L
]

genes[
  chr %in% c("MT", "M"),
  chr_order := 25L
]

setorder(
  genes,
  chr_order,
  start,
  end,
  gene_id
)

genes[, chr_order := NULL]

# Preferred output column order
setcolorder(
  genes,
  c(
    "gene_id",
    "gene_id_version",
    "gene_name",
    "gene_type",
    "broad_gene_class",
    "is_protein_coding",
    "is_lncRNA",
    "chr",
    "start",
    "end",
    "tss",
    "strand",
    "is_standard_chromosome",
    "is_autosomal"
  )
)

# ------------------------------------------------------------
# Save reference
# ------------------------------------------------------------

fwrite(
  genes,
  output_file,
  sep = "\t"
)

# ------------------------------------------------------------
# QC summary
# ------------------------------------------------------------

qc <- data.table(
  parameter = c(
    "gencode_release",
    "genome_build",
    "total_gtf_records",
    "gene_records_before_filtering",
    "genes_in_processed_reference",
    "autosomal_genes",
    "standard_chromosome_genes",
    "protein_coding_genes",
    "lncRNA_genes",
    "other_genes",
    "genes_with_unique_ensembl_id",
    "minimum_gene_length_bp",
    "median_gene_length_bp",
    "maximum_gene_length_bp"
  ),
  value = c(
    "v50",
    "GRCh38",
    length(gtf),
    length(genes_gr),
    nrow(genes),
    sum(genes$is_autosomal),
    sum(genes$is_standard_chromosome),
    sum(genes$is_protein_coding),
    sum(genes$is_lncRNA),
    sum(genes$broad_gene_class == "other"),
    uniqueN(genes$gene_id),
    min(genes$end - genes$start + 1L),
    median(genes$end - genes$start + 1L),
    max(genes$end - genes$start + 1L)
  )
)

fwrite(
  qc,
  output_qc,
  sep = "\t"
)

cat("\nGENCODE reference QC:\n")
print(qc)

cat("\nProcessed gene reference saved to:\n")
cat(output_file, "\n")

cat("\nQC report saved to:\n")
cat(output_qc, "\n")

cat("\nFirst 10 processed genes:\n")
print(
  genes[1:min(10L, .N)]
)

cat("\nTPCN2 record:\n")
print(
  genes[
    gene_name == "TPCN2"
  ]
)

cat("\nENSG00000261070 record:\n")
print(
  genes[
    gene_id == "ENSG00000261070"
  ]
)