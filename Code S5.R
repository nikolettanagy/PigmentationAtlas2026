suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# PigmentRegAtlas
# Module 04B: Annotate genomic loci with GENCODE genes
# ============================================================

locus_file <- "results/GCST90691600_genomic_loci_v3.tsv"

gene_reference_file <- paste0(
  "data/annotations/processed/",
  "GENCODE_v50_genes_GRCh38.tsv.gz"
)

output_gene_overlaps <- paste0(
  "results/",
  "GCST90691600_locus_gene_overlaps_v3.tsv"
)

output_master_annotation <- paste0(
  "results/",
  "GCST90691600_master_locus_annotation_v3.tsv"
)

output_qc <- paste0(
  "results/",
  "GCST90691600_gene_annotation_QC.tsv"
)

collapse_unique <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0L) {
    return(NA_character_)
  }
  paste(sort(x), collapse = ";")
}

cat("Reading loci:\n", locus_file, "\n")

loci <- fread(locus_file)

cat("Loci read:", nrow(loci), "\n")

cat("Reading processed GENCODE reference:\n")
cat(gene_reference_file, "\n")

genes <- fread(gene_reference_file)

cat("Genes read:", nrow(genes), "\n")

required_locus_columns <- c(
  "locus_id",
  "chr",
  "start",
  "end",
  "lead_id",
  "lead_pos",
  "lead_p"
)

missing_locus_columns <- setdiff(
  required_locus_columns,
  names(loci)
)

if (length(missing_locus_columns) > 0L) {
  stop(
    "Missing locus columns: ",
    paste(missing_locus_columns, collapse = ", ")
  )
}

required_gene_columns <- c(
  "gene_id",
  "gene_name",
  "gene_type",
  "broad_gene_class",
  "is_protein_coding",
  "is_lncRNA",
  "chr",
  "start",
  "end",
  "tss",
  "strand"
)

missing_gene_columns <- setdiff(
  required_gene_columns,
  names(genes)
)

if (length(missing_gene_columns) > 0L) {
  stop(
    "Missing gene columns: ",
    paste(missing_gene_columns, collapse = ", ")
  )
}

loci[, locus_index := .I]

# Autosomal atlas
genes <- genes[
  is_autosomal == TRUE
]

# ------------------------------------------------------------
# Gene overlaps
# ------------------------------------------------------------

cat("Calculating locus-gene overlaps...\n")

setkey(
  genes,
  chr,
  start,
  end
)

overlap_list <- vector(
  "list",
  nrow(loci)
)

for (i in seq_len(nrow(loci))) {

  locus_chr <- as.character(loci$chr[i])
  locus_start <- loci$start[i]
  locus_end <- loci$end[i]

  hits <- genes[
    chr == locus_chr &
      start <= locus_end &
      end >= locus_start
  ]

  if (nrow(hits) == 0L) {
    next
  }

  hits[, `:=`(
    locus_index = i,
    locus_id = loci$locus_id[i],
    locus_start = locus_start,
    locus_end = locus_end,
    lead_id = loci$lead_id[i],
    lead_pos = loci$lead_pos[i],
    lead_p = loci$lead_p[i]
  )]

  hits[, overlap_bp := pmax(
    0L,
    pmin(end, locus_end) -
      pmax(start, locus_start) +
      1L
  )]

  hits[, lead_inside_gene :=
    lead_pos >= start &
      lead_pos <= end
  ]

  overlap_list[[i]] <- hits
}

gene_overlaps <- rbindlist(
  overlap_list,
  fill = TRUE
)

if (nrow(gene_overlaps) > 0L) {

  setorder(
    gene_overlaps,
    locus_index,
    start,
    end
  )

  fwrite(
    gene_overlaps,
    output_gene_overlaps,
    sep = "\t"
  )
}

# ------------------------------------------------------------
# Summaries per locus
# ------------------------------------------------------------

if (nrow(gene_overlaps) > 0L) {

  overlap_summary <- gene_overlaps[, .(
    n_overlapping_genes = uniqueN(gene_id),

    overlapping_gene_ids = collapse_unique(
      gene_id
    ),

    overlapping_gene_names = collapse_unique(
      gene_name
    ),

    overlapping_gene_types = collapse_unique(
      gene_type
    ),

    n_overlapping_protein_coding = uniqueN(
      gene_id[is_protein_coding == TRUE]
    ),

    overlapping_protein_coding_ids = collapse_unique(
      gene_id[is_protein_coding == TRUE]
    ),

    overlapping_protein_coding_names = collapse_unique(
      gene_name[is_protein_coding == TRUE]
    ),

    n_overlapping_lncRNA = uniqueN(
      gene_id[is_lncRNA == TRUE]
    ),

    overlapping_lncRNA_ids = collapse_unique(
      gene_id[is_lncRNA == TRUE]
    ),

    overlapping_lncRNA_names = collapse_unique(
      gene_name[is_lncRNA == TRUE]
    ),

    genes_containing_lead = collapse_unique(
      gene_name[lead_inside_gene == TRUE]
    ),

    gene_ids_containing_lead = collapse_unique(
      gene_id[lead_inside_gene == TRUE]
    )

  ), by = .(
    locus_index,
    locus_id
  )]

} else {

  overlap_summary <- data.table(
    locus_index = integer(),
    locus_id = character()
  )
}

# ------------------------------------------------------------
# Nearest gene functions
# ------------------------------------------------------------

find_nearest_gene <- function(
  chromosome,
  position,
  gene_table
) {

  x <- gene_table[
    chr == chromosome
  ]

  if (nrow(x) == 0L) {
    return(NULL)
  }

  x[, distance_bp := fifelse(
    position < start,
    start - position,
    fifelse(
      position > end,
      position - end,
      0L
    )
  )]

  x <- x[
    order(
      distance_bp,
      abs(tss - position)
    )
  ]

  x[1]
}

nearest_all_list <- vector(
  "list",
  nrow(loci)
)

nearest_pc_list <- vector(
  "list",
  nrow(loci)
)

nearest_lnc_list <- vector(
  "list",
  nrow(loci)
)

for (i in seq_len(nrow(loci))) {

  chromosome <- as.character(loci$chr[i])
  position <- loci$lead_pos[i]

  nearest_all <- find_nearest_gene(
    chromosome,
    position,
    genes
  )

  nearest_pc <- find_nearest_gene(
    chromosome,
    position,
    genes[is_protein_coding == TRUE]
  )

  nearest_lnc <- find_nearest_gene(
    chromosome,
    position,
    genes[is_lncRNA == TRUE]
  )

  nearest_all_list[[i]] <- data.table(
    locus_index = i,
    nearest_gene_id = nearest_all$gene_id,
    nearest_gene_name = nearest_all$gene_name,
    nearest_gene_type = nearest_all$gene_type,
    nearest_gene_distance_bp = nearest_all$distance_bp
  )

  nearest_pc_list[[i]] <- data.table(
    locus_index = i,
    nearest_protein_coding_gene_id =
      nearest_pc$gene_id,
    nearest_protein_coding_gene_name =
      nearest_pc$gene_name,
    nearest_protein_coding_gene_type =
      nearest_pc$gene_type,
    nearest_protein_coding_distance_bp =
      nearest_pc$distance_bp
  )

  nearest_lnc_list[[i]] <- data.table(
    locus_index = i,
    nearest_lncRNA_gene_id =
      nearest_lnc$gene_id,
    nearest_lncRNA_gene_name =
      nearest_lnc$gene_name,
    nearest_lncRNA_gene_type =
      nearest_lnc$gene_type,
    nearest_lncRNA_distance_bp =
      nearest_lnc$distance_bp
  )
}

nearest_all_table <- rbindlist(
  nearest_all_list,
  fill = TRUE
)

nearest_pc_table <- rbindlist(
  nearest_pc_list,
  fill = TRUE
)

nearest_lnc_table <- rbindlist(
  nearest_lnc_list,
  fill = TRUE
)

# ------------------------------------------------------------
# Final master annotation
# ------------------------------------------------------------

master <- merge(
  loci,
  overlap_summary,
  by = c(
    "locus_index",
    "locus_id"
  ),
  all.x = TRUE
)

master <- merge(
  master,
  nearest_all_table,
  by = "locus_index",
  all.x = TRUE
)

master <- merge(
  master,
  nearest_pc_table,
  by = "locus_index",
  all.x = TRUE
)

master <- merge(
  master,
  nearest_lnc_table,
  by = "locus_index",
  all.x = TRUE
)

count_columns <- c(
  "n_overlapping_genes",
  "n_overlapping_protein_coding",
  "n_overlapping_lncRNA"
)

for (column_name in count_columns) {
  set(
    master,
    which(is.na(master[[column_name]])),
    column_name,
    0L
  )
}

setorder(
  master,
  locus_index
)

master[, locus_index := NULL]

fwrite(
  master,
  output_master_annotation,
  sep = "\t"
)

# ------------------------------------------------------------
# QC summary
# ------------------------------------------------------------

qc <- data.table(
  parameter = c(
    "number_of_loci",
    "loci_with_overlapping_gene",
    "loci_with_overlapping_protein_coding",
    "loci_with_overlapping_lncRNA",
    "median_overlapping_genes_per_locus",
    "maximum_overlapping_genes_per_locus",
    "median_overlapping_lncRNA_per_locus",
    "maximum_overlapping_lncRNA_per_locus",
    "loci_with_gene_containing_lead"
  ),
  value = c(
    nrow(master),
    sum(master$n_overlapping_genes > 0),
    sum(master$n_overlapping_protein_coding > 0),
    sum(master$n_overlapping_lncRNA > 0),
    median(master$n_overlapping_genes),
    max(master$n_overlapping_genes),
    median(master$n_overlapping_lncRNA),
    max(master$n_overlapping_lncRNA),
    sum(
      !is.na(master$genes_containing_lead) &
        master$genes_containing_lead != ""
    )
  )
)

fwrite(
  qc,
  output_qc,
  sep = "\t"
)

cat("\nModule 04B QC:\n")
print(qc)

cat(
  "\nMaster annotation saved to:\n",
  output_master_annotation,
  "\n"
)

cat(
  "Detailed overlap table saved to:\n",
  output_gene_overlaps,
  "\n"
)

cat(
  "QC saved to:\n",
  output_qc,
  "\n"
)

cat("\nTPCN2 locus:\n")

print(
  master[
    chr == "11" &
      lead_pos >= 68000000 &
      lead_pos <= 70000000
  ]
)
