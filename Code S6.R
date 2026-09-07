suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# PigmentRegAtlas
# Module 04C v2
# Structural candidate-gene prioritization
#
# Candidate genes:
#   1. all genes overlapping the genomic locus
#   2. nearest protein-coding gene to the lead variant
#   3. nearest lncRNA to the lead variant
#
# Molecular evidence such as eQTL, colocalization,
# fine-mapping and regulatory annotation will be integrated
# in later modules.
# ============================================================

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

trait_id <- "GCST90691600"

locus_file <- file.path(
  "results",
  paste0(trait_id, "_genomic_loci_v3.tsv")
)

overlap_file <- file.path(
  "results",
  paste0(trait_id, "_locus_gene_overlaps_v3.tsv")
)

# Adjust only this path if the reference has another filename.
gene_reference_file <- file.path(
  "data",
  "annotations",
  "processed",
  "GENCODE_v50_genes_GRCh38.tsv.gz"
)

output_all <- file.path(
  "results",
  paste0(trait_id, "_candidate_gene_prioritization_v2.tsv")
)

output_top <- file.path(
  "results",
  paste0(trait_id, "_top_candidate_genes_v2.tsv")
)

output_summary <- file.path(
  "results",
  paste0(trait_id, "_candidate_gene_summary_v2.tsv")
)

output_qc <- file.path(
  "results",
  paste0(trait_id, "_candidate_gene_QC_v2.tsv")
)

output_tpcn2 <- file.path(
  "results",
  paste0(trait_id, "_TPCN2_candidate_genes_v2.tsv")
)

output_missing_loci <- file.path(
  "results",
  paste0(trait_id, "_previously_unannotated_loci_v2.tsv")
)

top_n_per_locus <- 10L

# ------------------------------------------------------------
# Input validation
# ------------------------------------------------------------

input_files <- c(
  locus_file,
  overlap_file,
  gene_reference_file
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0L) {
  stop(
    "The following input files were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

cat("Reading input files...\n")

loci <- fread(
  locus_file,
  na.strings = c("", "NA", "NaN")
)

overlaps <- fread(
  overlap_file,
  na.strings = c("", "NA", "NaN")
)

reference <- fread(
  gene_reference_file,
  na.strings = c("", "NA", "NaN")
)

cat("Loci:", nrow(loci), "\n")
cat("Gene-locus overlap records:", nrow(overlaps), "\n")
cat("GENCODE reference genes:", nrow(reference), "\n\n")

# ------------------------------------------------------------
# Standardize chromosome representation
# ------------------------------------------------------------

normalize_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x
}

for (object_name in c("loci", "overlaps", "reference")) {
  object <- get(object_name)

  if (!"chr" %in% names(object)) {
    stop("Missing chr column in object: ", object_name)
  }

  object[, chr := normalize_chr(chr)]

  assign(
    object_name,
    object,
    envir = .GlobalEnv
  )
}

# ------------------------------------------------------------
# Standardize locus coordinate columns
# ------------------------------------------------------------

if (!"locus_start" %in% names(loci)) {
  if ("start" %in% names(loci)) {
    setnames(loci, "start", "locus_start")
  } else {
    stop("No locus_start or start column found in locus file.")
  }
}

if (!"locus_end" %in% names(loci)) {
  if ("end" %in% names(loci)) {
    setnames(loci, "end", "locus_end")
  } else {
    stop("No locus_end or end column found in locus file.")
  }
}

locus_required <- c(
  "locus_id",
  "chr",
  "locus_start",
  "locus_end",
  "lead_id",
  "lead_pos",
  "lead_p"
)

missing_locus_columns <- setdiff(
  locus_required,
  names(loci)
)

if (length(missing_locus_columns) > 0L) {
  stop(
    "Missing locus columns: ",
    paste(missing_locus_columns, collapse = ", ")
  )
}

loci[
  ,
  `:=`(
    locus_id = trimws(as.character(locus_id)),
    locus_start = as.numeric(locus_start),
    locus_end = as.numeric(locus_end),
    lead_pos = as.numeric(lead_pos),
    lead_p = as.numeric(lead_p)
  )
]

# Keep one row per locus
setorder(loci, locus_id, lead_p)

loci <- loci[
  ,
  .SD[1L],
  by = locus_id
]

# ------------------------------------------------------------
# Standardize GENCODE reference
# ------------------------------------------------------------

reference_required <- c(
  "gene_id",
  "gene_name",
  "gene_type",
  "chr",
  "start",
  "end",
  "strand"
)

missing_reference_columns <- setdiff(
  reference_required,
  names(reference)
)

if (length(missing_reference_columns) > 0L) {
  stop(
    "Missing GENCODE reference columns: ",
    paste(missing_reference_columns, collapse = ", ")
  )
}

reference[
  ,
  `:=`(
    gene_id = as.character(gene_id),
    gene_name = as.character(gene_name),
    gene_type = as.character(gene_type),
    start = as.numeric(start),
    end = as.numeric(end),
    strand = as.character(strand)
  )
]

# Remove Ensembl version only if gene_id contains one
reference[
  ,
  gene_id_unversioned := sub(
    "\\.[0-9]+$",
    "",
    gene_id
  )
]

if (!"gene_id_version" %in% names(reference)) {
  reference[, gene_id_version := gene_id]
}

reference[, gene_id := gene_id_unversioned]
reference[, gene_id_unversioned := NULL]

# Derive TSS when absent
if (!"tss" %in% names(reference)) {
  reference[
    ,
    tss := fifelse(
      strand == "-",
      end,
      start
    )
  ]
}

reference[, tss := as.numeric(tss)]

# Derive protein-coding indicator when absent
if (!"is_protein_coding" %in% names(reference)) {
  reference[
    ,
    is_protein_coding := gene_type == "protein_coding"
  ]
}

# Derive lncRNA indicator when absent
if (!"is_lncRNA" %in% names(reference)) {
  lncRNA_types <- c(
    "lncRNA",
    "lincRNA",
    "antisense",
    "sense_intronic",
    "sense_overlapping",
    "processed_transcript",
    "3prime_overlapping_ncRNA",
    "bidirectional_promoter_lncRNA",
    "macro_lncRNA",
    "non_coding"
  )

  reference[
    ,
    is_lncRNA := gene_type %in% lncRNA_types
  ]
}

reference[
  ,
  `:=`(
    is_protein_coding = as.logical(is_protein_coding),
    is_lncRNA = as.logical(is_lncRNA)
  )
]

reference[
  is.na(is_protein_coding),
  is_protein_coding := FALSE
]

reference[
  is.na(is_lncRNA),
  is_lncRNA := FALSE
]

if (!"broad_gene_class" %in% names(reference)) {
  reference[
    ,
    broad_gene_class := fcase(
      is_protein_coding, "protein_coding",
      is_lncRNA, "lncRNA",
      default = "other"
    )
  ]
}

reference <- reference[
  !is.na(chr) &
    !is.na(start) &
    !is.na(end) &
    !is.na(gene_id) &
    gene_id != ""
]

# One row per gene
setorder(
  reference,
  gene_id,
  -is_protein_coding,
  -is_lncRNA
)

reference <- reference[
  ,
  .SD[1L],
  by = gene_id
]

# ------------------------------------------------------------
# Standardize overlap table
# ------------------------------------------------------------

overlap_required <- c(
  "locus_id",
  "gene_id"
)

missing_overlap_columns <- setdiff(
  overlap_required,
  names(overlaps)
)

if (length(missing_overlap_columns) > 0L) {
  stop(
    "Missing overlap-table columns: ",
    paste(missing_overlap_columns, collapse = ", ")
  )
}

overlaps[
  ,
  `:=`(
    locus_id = trimws(as.character(locus_id)),
    gene_id = sub(
      "\\.[0-9]+$",
      "",
      as.character(gene_id)
    )
  )
]

# Only locus-gene membership is required here.
overlap_membership <- unique(
  overlaps[
    !is.na(locus_id) &
      locus_id != "" &
      !is.na(gene_id) &
      gene_id != "",
    .(
      locus_id,
      gene_id,
      overlap_bp = if (
        "overlap_bp" %in% names(overlaps)
      ) {
        as.numeric(overlap_bp)
      } else {
        NA_real_
      }
    )
  ]
)

# In case duplicate transcript-derived rows exist,
# retain the maximum overlap for each locus-gene pair.
overlap_membership <- overlap_membership[
  ,
  .(
    overlap_bp = if (all(is.na(overlap_bp))) {
      NA_real_
    } else {
      max(overlap_bp, na.rm = TRUE)
    }
  ),
  by = .(
    locus_id,
    gene_id
  )
]

# ------------------------------------------------------------
# Identify loci that previously had no overlapping gene
# ------------------------------------------------------------

previously_unannotated <- loci[
  !unique(overlap_membership$locus_id),
  on = "locus_id"
]

fwrite(
  previously_unannotated,
  output_missing_loci,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat(
  "Previously unannotated loci:",
  nrow(previously_unannotated),
  "\n"
)

# ------------------------------------------------------------
# Create candidates from overlapping genes
# ------------------------------------------------------------

overlap_candidates <- merge(
  overlap_membership,
  reference,
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)

overlap_candidates <- merge(
  overlap_candidates,
  loci,
  by = "locus_id",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("_gene", "_locus")
)

# Prefer gene chromosome after merge
if ("chr_gene" %in% names(overlap_candidates)) {
  overlap_candidates[, chr := chr_gene]
} else if ("chr" %in% names(overlap_candidates)) {
  overlap_candidates[, chr := chr]
}

overlap_candidates[
  ,
  candidate_source := "overlapping_gene"
]

# ------------------------------------------------------------
# Find nearest gene of a requested class
# ------------------------------------------------------------

find_nearest_gene <- function(
  locus_table,
  gene_table,
  class_column,
  source_label
) {

  eligible_genes <- gene_table[
    get(class_column) == TRUE
  ]

  if (nrow(eligible_genes) == 0L) {
    stop(
      "No eligible genes found for class: ",
      class_column
    )
  }

  results <- vector(
    "list",
    nrow(locus_table)
  )

  for (i in seq_len(nrow(locus_table))) {

    current_locus <- locus_table[i]

    chromosome_genes <- eligible_genes[
      chr == current_locus$chr
    ]

    if (nrow(chromosome_genes) == 0L) {
      next
    }

    chromosome_genes[
      ,
      temporary_distance := fifelse(
        current_locus$lead_pos < start,
        start - current_locus$lead_pos,
        fifelse(
          current_locus$lead_pos > end,
          current_locus$lead_pos - end,
          0
        )
      )
    ]

    chromosome_genes[
      ,
      temporary_tss_distance :=
        abs(current_locus$lead_pos - tss)
    ]

    setorder(
      chromosome_genes,
      temporary_distance,
      temporary_tss_distance,
      gene_id
    )

    nearest <- chromosome_genes[1L]

    result <- cbind(
      current_locus,
      nearest,
      fill = TRUE
    )

    result[
      ,
      candidate_source := source_label
    ]

    result[
      ,
      overlap_bp := 0
    ]

    result[
      ,
      c(
        "temporary_distance",
        "temporary_tss_distance"
      ) := NULL
    ]

    results[[i]] <- result
  }

  rbindlist(
    results,
    fill = TRUE,
    use.names = TRUE
  )
}

cat("Finding nearest protein-coding genes...\n")

nearest_pc <- find_nearest_gene(
  locus_table = loci,
  gene_table = reference,
  class_column = "is_protein_coding",
  source_label = "nearest_protein_coding"
)

cat("Finding nearest lncRNAs...\n")

nearest_lncRNA <- find_nearest_gene(
  locus_table = loci,
  gene_table = reference,
  class_column = "is_lncRNA",
  source_label = "nearest_lncRNA"
)

# ------------------------------------------------------------
# Standardize columns before combining
# ------------------------------------------------------------

standard_candidate_columns <- c(
  "locus_id",
  "chr",
  "locus_start",
  "locus_end",
  "lead_id",
  "lead_pos",
  "lead_p",
  "gene_id",
  "gene_id_version",
  "gene_name",
  "gene_type",
  "broad_gene_class",
  "is_protein_coding",
  "is_lncRNA",
  "start",
  "end",
  "tss",
  "strand",
  "overlap_bp",
  "candidate_source"
)

standardize_candidate_table <- function(x) {

  for (column_name in standard_candidate_columns) {
    if (!column_name %in% names(x)) {
      x[, (column_name) := NA]
    }
  }

  x[
    ,
    ..standard_candidate_columns
  ]
}

overlap_candidates <- standardize_candidate_table(
  overlap_candidates
)

nearest_pc <- standardize_candidate_table(
  nearest_pc
)

nearest_lncRNA <- standardize_candidate_table(
  nearest_lncRNA
)

# ------------------------------------------------------------
# Combine and collapse duplicate locus-gene entries
# ------------------------------------------------------------

candidates <- rbindlist(
  list(
    overlap_candidates,
    nearest_pc,
    nearest_lncRNA
  ),
  fill = TRUE,
  use.names = TRUE
)

candidates[
  ,
  gene_id := sub(
    "\\.[0-9]+$",
    "",
    as.character(gene_id)
  )
]

candidates <- candidates[
  !is.na(locus_id) &
    locus_id != "" &
    !is.na(gene_id) &
    gene_id != ""
]

# Collapse multiple candidate sources
candidates <- candidates[
  ,
  {
    source_values <- unique(
      unlist(
        strsplit(
          candidate_source,
          ";",
          fixed = TRUE
        )
      )
    )

    source_order <- c(
      "overlapping_gene",
      "nearest_protein_coding",
      "nearest_lncRNA"
    )

    source_values <- source_order[
      source_order %in% source_values
    ]

    selected <- .SD[1L]

    selected[
      ,
      candidate_source := paste(
        source_values,
        collapse = ";"
      )
    ]

    if (
      "overlapping_gene" %in% source_values &&
        any(overlap_bp > 0, na.rm = TRUE)
    ) {
      selected[
        ,
        overlap_bp := max(
          overlap_bp,
          na.rm = TRUE
        )
      ]
    } else {
      selected[
        ,
        overlap_bp := 0
      ]
    }

    selected
  },
  by = .(
    locus_id,
    gene_id
  )
]

# Remove grouping columns duplicated by .SD
duplicate_group_columns <- intersect(
  c("locus_id", "gene_id"),
  names(candidates)[duplicated(names(candidates))]
)

if (length(duplicate_group_columns) > 0L) {
  candidates[
    ,
    (duplicate_group_columns) := NULL
  ]
}

# ------------------------------------------------------------
# Reattach authoritative locus information
# ------------------------------------------------------------

gene_columns <- candidates[
  ,
  !c(
    "chr",
    "locus_start",
    "locus_end",
    "lead_id",
    "lead_pos",
    "lead_p"
  ),
  with = FALSE
]

candidates <- merge(
  gene_columns,
  loci[
    ,
    .(
      locus_id,
      chr,
      locus_start,
      locus_end,
      lead_id,
      lead_pos,
      lead_p
    )
  ],
  by = "locus_id",
  all.x = TRUE,
  sort = FALSE
)

# ------------------------------------------------------------
# Calculate structural evidence
# ------------------------------------------------------------

candidates[
  ,
  `:=`(
    start = as.numeric(start),
    end = as.numeric(end),
    tss = as.numeric(tss),
    lead_pos = as.numeric(lead_pos),
    locus_start = as.numeric(locus_start),
    locus_end = as.numeric(locus_end),
    overlap_bp = as.numeric(overlap_bp)
  )
]

candidates[
  is.na(overlap_bp),
  overlap_bp := 0
]

candidates[
  ,
  lead_inside_gene :=
    lead_pos >= start &
    lead_pos <= end
]

candidates[
  ,
  gene_overlaps_locus :=
    start <= locus_end &
    end >= locus_start
]

candidates[
  ,
  distance_lead_to_gene_bp := fifelse(
    lead_pos < start,
    start - lead_pos,
    fifelse(
      lead_pos > end,
      lead_pos - end,
      0
    )
  )
]

candidates[
  ,
  distance_lead_to_tss_bp :=
    abs(lead_pos - tss)
]

candidates[
  ,
  distance_locus_to_gene_bp := fifelse(
    end < locus_start,
    locus_start - end,
    fifelse(
      start > locus_end,
      start - locus_end,
      0
    )
  )
]

candidates[
  ,
  gene_length_bp := end - start + 1
]

candidates[
  ,
  overlap_fraction_gene := fifelse(
    overlap_bp > 0 &
      gene_length_bp > 0,
    pmin(
      1,
      overlap_bp / gene_length_bp
    ),
    0
  )
]

# ------------------------------------------------------------
# Transparent scoring
# ------------------------------------------------------------

candidates[
  ,
  score_lead_inside_gene := fifelse(
    lead_inside_gene,
    5,
    0
  )
]

candidates[
  ,
  score_gene_distance := fcase(
    distance_lead_to_gene_bp == 0, 4,
    distance_lead_to_gene_bp <= 10000, 3,
    distance_lead_to_gene_bp <= 50000, 2,
    distance_lead_to_gene_bp <= 250000, 1,
    default = 0
  )
]

candidates[
  ,
  score_tss_distance := fcase(
    distance_lead_to_tss_bp <= 10000, 3,
    distance_lead_to_tss_bp <= 50000, 2,
    distance_lead_to_tss_bp <= 250000, 1,
    default = 0
  )
]

candidates[
  ,
  score_gene_class := fifelse(
    is_protein_coding | is_lncRNA,
    1,
    0
  )
]

candidates[
  ,
  score_overlap := fcase(
    overlap_fraction_gene >= 0.90, 1,
    overlap_fraction_gene >= 0.50, 0.5,
    default = 0
  )
]

# Candidate-source membership is recorded but does not receive
# an additional score. This avoids double-counting proximity.
candidates[
  ,
  structural_priority_score :=
    score_lead_inside_gene +
    score_gene_distance +
    score_tss_distance +
    score_gene_class +
    score_overlap
]

# ------------------------------------------------------------
# Evidence labels
# ------------------------------------------------------------

candidates[
  ,
  structural_evidence := {

    evidence <- character()

    if (gene_overlaps_locus) {
      evidence <- c(
        evidence,
        "gene_overlaps_locus"
      )
    }

    if (lead_inside_gene) {
      evidence <- c(
        evidence,
        "lead_inside_gene"
      )
    }

    if (distance_lead_to_gene_bp <= 10000) {
      evidence <- c(
        evidence,
        "gene_within_10kb"
      )
    } else if (
      distance_lead_to_gene_bp <= 50000
    ) {
      evidence <- c(
        evidence,
        "gene_within_50kb"
      )
    } else if (
      distance_lead_to_gene_bp <= 250000
    ) {
      evidence <- c(
        evidence,
        "gene_within_250kb"
      )
    }

    if (distance_lead_to_tss_bp <= 10000) {
      evidence <- c(
        evidence,
        "TSS_within_10kb"
      )
    } else if (
      distance_lead_to_tss_bp <= 50000
    ) {
      evidence <- c(
        evidence,
        "TSS_within_50kb"
      )
    } else if (
      distance_lead_to_tss_bp <= 250000
    ) {
      evidence <- c(
        evidence,
        "TSS_within_250kb"
      )
    }

    if (is_protein_coding) {
      evidence <- c(
        evidence,
        "protein_coding"
      )
    } else if (is_lncRNA) {
      evidence <- c(
        evidence,
        "lncRNA"
      )
    }

    if (length(evidence) == 0L) {
      evidence <- "nearest_gene_only"
    }

    paste(
      unique(evidence),
      collapse = ";"
    )
  },
  by = seq_len(nrow(candidates))
]

# ------------------------------------------------------------
# Ranking within loci
# ------------------------------------------------------------

setorder(
  candidates,
  locus_id,
  -structural_priority_score,
  distance_lead_to_gene_bp,
  distance_lead_to_tss_bp,
  -gene_overlaps_locus,
  gene_id
)

candidates[
  ,
  candidate_rank := seq_len(.N),
  by = locus_id
]

candidates[
  ,
  candidate_tier := fcase(
    candidate_rank == 1L, "Tier 1",
    candidate_rank <= 3L, "Tier 2",
    candidate_rank <= 10L, "Tier 3",
    default = "Tier 4"
  )
]

candidates[
  ,
  `:=`(
    is_top_candidate =
      candidate_rank == 1L,

    is_top_three_candidate =
      candidate_rank <= 3L,

    is_top_ten_candidate =
      candidate_rank <= top_n_per_locus
  )
]

# ------------------------------------------------------------
# Reserve columns for later molecular evidence
# ------------------------------------------------------------

candidates[
  ,
  `:=`(
    eqtl_evidence_score = NA_real_,
    colocalization_evidence_score = NA_real_,
    fine_mapping_evidence_score = NA_real_,
    regulatory_annotation_score = NA_real_,
    expression_specificity_score = NA_real_,
    integrated_priority_score = NA_real_
  )
]

# ------------------------------------------------------------
# Column order
# ------------------------------------------------------------

preferred_order <- c(
  "locus_id",
  "candidate_rank",
  "candidate_tier",
  "is_top_candidate",
  "is_top_three_candidate",
  "is_top_ten_candidate",
  "candidate_source",
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
  "locus_start",
  "locus_end",
  "lead_id",
  "lead_pos",
  "lead_p",
  "gene_overlaps_locus",
  "lead_inside_gene",
  "distance_lead_to_gene_bp",
  "distance_lead_to_tss_bp",
  "distance_locus_to_gene_bp",
  "overlap_bp",
  "gene_length_bp",
  "overlap_fraction_gene",
  "score_lead_inside_gene",
  "score_gene_distance",
  "score_tss_distance",
  "score_gene_class",
  "score_overlap",
  "structural_priority_score",
  "structural_evidence",
  "eqtl_evidence_score",
  "colocalization_evidence_score",
  "fine_mapping_evidence_score",
  "regulatory_annotation_score",
  "expression_specificity_score",
  "integrated_priority_score"
)

setcolorder(
  candidates,
  intersect(
    preferred_order,
    names(candidates)
  )
)

# ------------------------------------------------------------
# Save complete results
# ------------------------------------------------------------

fwrite(
  candidates,
  output_all,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

top_candidates <- candidates[
  candidate_rank <= top_n_per_locus
]

fwrite(
  top_candidates,
  output_top,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# Locus summary
# ------------------------------------------------------------

locus_summary <- candidates[
  order(candidate_rank),
  {
    top_three <- .SD[
      candidate_rank <= 3L
    ]

    top_ten <- .SD[
      candidate_rank <= 10L
    ]

    list(
      chr = chr[1L],
      locus_start = locus_start[1L],
      locus_end = locus_end[1L],
      lead_id = lead_id[1L],
      lead_pos = lead_pos[1L],
      lead_p = lead_p[1L],

      n_candidate_genes = .N,

      n_overlapping_genes = sum(
        gene_overlaps_locus,
        na.rm = TRUE
      ),

      top_gene_id = gene_id[1L],
      top_gene_name = gene_name[1L],
      top_gene_type = gene_type[1L],
      top_gene_source = candidate_source[1L],
      top_gene_score =
        structural_priority_score[1L],
      top_gene_distance_bp =
        distance_lead_to_gene_bp[1L],

      nearest_protein_coding_gene = {
        value <- gene_name[
          grepl(
            "nearest_protein_coding",
            candidate_source,
            fixed = TRUE
          )
        ]

        if (length(value) == 0L) {
          NA_character_
        } else {
          value[1L]
        }
      },

      nearest_lncRNA = {
        value <- gene_name[
          grepl(
            "nearest_lncRNA",
            candidate_source,
            fixed = TRUE
          )
        ]

        if (length(value) == 0L) {
          NA_character_
        } else {
          value[1L]
        }
      },

      top_three_gene_ids = paste(
        top_three$gene_id,
        collapse = ";"
      ),

      top_three_gene_names = paste(
        top_three$gene_name,
        collapse = ";"
      ),

      top_ten_gene_names = paste(
        top_ten$gene_name,
        collapse = ";"
      )
    )
  },
  by = locus_id
]

fwrite(
  locus_summary,
  output_summary,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# TPCN2 locus
# ------------------------------------------------------------

tpcn2_candidates <- candidates[
  locus_id == "LOCUS_0131_chr11"
]

fwrite(
  tpcn2_candidates,
  output_tpcn2,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# QC
# ------------------------------------------------------------

candidate_counts <- candidates[
  ,
  .N,
  by = locus_id
]

annotated_locus_count <- uniqueN(
  candidates$locus_id
)

missing_final_loci <- loci[
  !unique(candidates$locus_id),
  on = "locus_id"
]

top_genes <- candidates[
  candidate_rank == 1L
]

qc <- data.table(
  parameter = c(
    "trait_id",
    "algorithm",
    "number_of_input_loci",
    "number_of_overlap_records",
    "number_of_final_candidate_records",
    "number_of_annotated_loci",
    "number_of_unannotated_loci",
    "previously_unannotated_loci_rescued",
    "median_candidate_genes_per_locus",
    "mean_candidate_genes_per_locus",
    "maximum_candidate_genes_per_locus",
    "top_candidate_protein_coding",
    "top_candidate_lncRNA",
    "top_candidate_other",
    "top_candidate_contains_lead",
    "top_candidate_overlaps_locus",
    "median_top_candidate_distance_bp",
    "maximum_structural_priority_score"
  ),

  value = c(
    trait_id,
    paste0(
      "structural_prioritization_with_",
      "overlap_and_nearest_gene_rescue_v2"
    ),
    nrow(loci),
    nrow(overlap_membership),
    nrow(candidates),
    annotated_locus_count,
    nrow(missing_final_loci),
    nrow(previously_unannotated),
    median(candidate_counts$N),
    round(mean(candidate_counts$N), 2),
    max(candidate_counts$N),

    top_genes[
      is_protein_coding == TRUE,
      .N
    ],

    top_genes[
      is_lncRNA == TRUE,
      .N
    ],

    top_genes[
      !is_protein_coding &
        !is_lncRNA,
      .N
    ],

    top_genes[
      lead_inside_gene == TRUE,
      .N
    ],

    top_genes[
      gene_overlaps_locus == TRUE,
      .N
    ],

    median(
      top_genes$distance_lead_to_gene_bp,
      na.rm = TRUE
    ),

    max(
      candidates$structural_priority_score,
      na.rm = TRUE
    )
  )
)

fwrite(
  qc,
  output_qc,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# Console output
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("PigmentRegAtlas Module 04C v2 completed\n")
cat("============================================================\n\n")

print(qc)

cat("\nPreviously unannotated loci after rescue:\n")

print(
  candidates[
    locus_id %in% previously_unannotated$locus_id,
    .(
      locus_id,
      candidate_rank,
      gene_id,
      gene_name,
      gene_type,
      candidate_source,
      distance_lead_to_gene_bp,
      distance_lead_to_tss_bp,
      structural_priority_score
    )
  ][
    candidate_rank <= 5L
  ]
)

cat("\nTPCN2 locus top candidates:\n")

print(
  tpcn2_candidates[
    1:min(15L, .N),
    .(
      candidate_rank,
      gene_id,
      gene_name,
      gene_type,
      candidate_source,
      lead_inside_gene,
      distance_lead_to_gene_bp,
      distance_lead_to_tss_bp,
      structural_priority_score
    )
  ]
)

cat("\nOutput files:\n")
cat(output_all, "\n")
cat(output_top, "\n")
cat(output_summary, "\n")
cat(output_qc, "\n")
cat(output_tpcn2, "\n")
cat(output_missing_loci, "\n")