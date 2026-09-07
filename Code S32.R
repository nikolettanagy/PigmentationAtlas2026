###############################################################################
# PigmentationAtlas
# Module 12.6 – Colocalization-driven Gene Reclassification
#
# Purpose:
#   Compare conventional locus-to-gene assignment
#   (overlapping or nearest protein-coding gene)
#   with the best colocalization-supported candidate gene.
#
# Main outputs:
#   results/module12/reclassification/
#
# Author: PigmentationAtlas
###############################################################################

options(stringsAsFactors = FALSE)
options(scipen = 999)

###############################################################################
# 1. Packages
###############################################################################

required_packages <- c(
  "data.table"
)

install_missing_packages <- function(packages) {
  missing <- packages[
    !vapply(
      packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(missing) > 0L) {
    install.packages(
      missing,
      repos = "https://cloud.r-project.org"
    )
  }
}

install_missing_packages(required_packages)

suppressPackageStartupMessages({
  library(data.table)
})

###############################################################################
# 2. Project paths
###############################################################################

project_dir <- normalizePath(
  "C:/Users/User/Desktop/PigmentationAtlas",
  winslash = "/",
  mustWork = TRUE
)

results_dir <- file.path(project_dir, "results")

output_dir <- file.path(
  results_dir,
  "module12",
  "reclassification"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

clump_metadata_file <- file.path(
  results_dir,
  "clump_master_metadata.tsv"
)

prioritization_file <- file.path(
  results_dir,
  "module12",
  "colocalization",
  "biological_annotation",
  "gene_level",
  "module12_5C_best_gene_per_clump_GENCODEv50.tsv.gz"
)

###############################################################################
# 3. Helper functions
###############################################################################

message_section <- function(text) {
  cat(
    "\n",
    paste(rep("=", 78), collapse = ""),
    "\n",
    text,
    "\n",
    paste(rep("=", 78), collapse = ""),
    "\n",
    sep = ""
  )
}

assert_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      label,
      " was not found:\n",
      path,
      call. = FALSE
    )
  }
}

assert_columns <- function(data, columns, label) {
  missing_columns <- setdiff(columns, names(data))

  if (length(missing_columns) > 0L) {
    stop(
      label,
      " is missing required columns:\n",
      paste(missing_columns, collapse = ", "),
      "\n\nAvailable columns:\n",
      paste(names(data), collapse = ", "),
      call. = FALSE
    )
  }
}

first_existing_column <- function(data, candidates) {
  hit <- candidates[candidates %in% names(data)]

  if (length(hit) == 0L) {
    return(NA_character_)
  }

  hit[[1L]]
}

standardize_chromosome <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x
}

stable_ensembl_id <- function(x) {
  sub("\\..*$", "", as.character(x))
}

`%+%` <- function(x, y) {
  paste0(x, y)
}

normalize_biotype <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[ -]", "_", x)
  x
}

is_protein_coding <- function(x) {
  normalize_biotype(x) %in% c(
    "protein_coding",
    "protein_coding_gene"
  )
}

classify_candidate_biotype <- function(biotype, candidate_class = NA_character_) {
  b <- normalize_biotype(biotype)
  c <- tolower(trimws(as.character(candidate_class)))

  result <- rep("OTHER_OR_UNKNOWN", length(b))

  result[grepl("protein", b) & grepl("coding", b)] <-
    "PROTEIN_CODING"

  result[
    grepl(
      "lncrna|linc_rna|antisense|sense_intronic|sense_overlapping|" %+%
        "processed_transcript|macro_lncRNA|bidirectional_promoter_lncRNA",
      b,
      ignore.case = TRUE
    )
  ] <- "REGULATORY_RNA"

  result[
    grepl(
      "pseudogene",
      b,
      ignore.case = TRUE
    )
  ] <- "PSEUDOGENE"

  result[
    grepl("regulatory", c, ignore.case = TRUE) |
      grepl("rna", c, ignore.case = TRUE)
  ] <- "REGULATORY_RNA"

  result[
    grepl("pseudogene", c, ignore.case = TRUE)
  ] <- "PSEUDOGENE"

  result[
    grepl("protein", c, ignore.case = TRUE) &
      grepl("coding", c, ignore.case = TRUE)
  ] <- "PROTEIN_CODING"

  result
}

priority_class_rank <- function(x) {
  normalized <- toupper(
    gsub("[ -]", "_", trimws(as.character(x)))
  )

  ranks <- c(
    "STRONG_H4" = 1L,
    "MODERATE_H4" = 2L,
    "SUGGESTIVE_H4" = 3L,
    "H3_DOMINANT" = 4L,
    "INCONCLUSIVE" = 5L
  )

  result <- unname(ranks[normalized])
  result[is.na(result)] <- 99L
  as.integer(result)
}

###############################################################################
# 4. Locate GENCODE v50 gene annotation
###############################################################################

find_gencode_annotation <- function(project_dir) {

  all_files <- list.files(
    project_dir,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
  )

  file_names <- basename(all_files)

  candidates <- all_files[
    grepl(
      "gencode",
      file_names,
      ignore.case = TRUE
    ) &
      grepl(
        "v50|version.?50",
        file_names,
        ignore.case = TRUE
      ) &
      grepl(
        "\\.(gtf|gtf\\.gz|tsv|tsv\\.gz|csv|csv\\.gz)$",
        file_names,
        ignore.case = TRUE
      )
  ]

  candidates <- unique(candidates)

  if (length(candidates) == 0L) {
    stop(
      paste(
        "No GENCODE v50 annotation file was found.",
        "",
        "Expected a filename containing both:",
        "  gencode",
        "  v50",
        "",
        "Supported formats:",
        "  .gtf",
        "  .gtf.gz",
        "  .tsv",
        "  .tsv.gz",
        "  .csv",
        "  .csv.gz",
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  candidate_info <- data.table(
    file = candidates,
    size_bytes = file.info(candidates)$size
  )

  fwrite(
    candidate_info,
    file.path(
      output_dir,
      "module12_6_gencode_candidate_files.tsv"
    ),
    sep = "\t"
  )

  # Prefer files that look like gene-level tables.
  candidate_info[
    ,
    score := 0L
  ]

  candidate_info[
    grepl("gene", basename(file), ignore.case = TRUE),
    score := score + 10L
  ]

  candidate_info[
    grepl("annotation", basename(file), ignore.case = TRUE),
    score := score + 5L
  ]

  candidate_info[
    grepl("\\.tsv|\\.csv", basename(file), ignore.case = TRUE),
    score := score + 3L
  ]

  candidate_info[
    order(-score, -size_bytes)
  ]$file[[1L]]
}

###############################################################################
# 5. Read GENCODE annotation
###############################################################################

extract_gtf_attribute <- function(attribute, key) {
  pattern <- paste0(
    ".*(?:^|;[[:space:]]*)",
    key,
    "[[:space:]]+\"([^\"]+)\".*"
  )

  result <- sub(pattern, "\\1", attribute)

  no_match <- !grepl(
    paste0(
      "(?:^|;[[:space:]]*)",
      key,
      "[[:space:]]+\""
    ),
    attribute
  )

  result[no_match] <- NA_character_
  result
}

read_gencode_annotation <- function(path) {

  message("GENCODE annotation selected:")
  message(path)

  if (grepl("\\.gtf(\\.gz)?$", path, ignore.case = TRUE)) {

    gtf <- fread(
      path,
      sep = "\t",
      header = FALSE,
      quote = "",
      fill = TRUE,
      comment.char = "#",
      showProgress = TRUE
    )

    if (ncol(gtf) < 9L) {
      stop(
        "The selected GTF file has fewer than 9 columns.",
        call. = FALSE
      )
    }

    setnames(
      gtf,
      names(gtf)[1:9],
      c(
        "chromosome",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute"
      )
    )

    genes <- gtf[
      feature == "gene"
    ]

    genes[
      ,
      gene_id := extract_gtf_attribute(attribute, "gene_id")
    ]

    genes[
      ,
      gene_symbol := extract_gtf_attribute(attribute, "gene_name")
    ]

    genes[
      ,
      gene_biotype := extract_gtf_attribute(attribute, "gene_type")
    ]

    genes <- genes[
      ,
      .(
        chromosome,
        gene_start = as.integer(start),
        gene_end = as.integer(end),
        strand,
        gene_id,
        gene_symbol,
        gene_biotype
      )
    ]

  } else {

    annotation <- fread(
      path,
      showProgress = TRUE
    )

    chr_col <- first_existing_column(
      annotation,
      c(
        "chromosome",
        "chr",
        "seqname",
        "seqnames",
        "gene_chromosome"
      )
    )

    start_col <- first_existing_column(
      annotation,
      c(
        "gene_start",
        "start",
        "start_position",
        "gene_start_position"
      )
    )

    end_col <- first_existing_column(
      annotation,
      c(
        "gene_end",
        "end",
        "end_position",
        "gene_end_position"
      )
    )

    id_col <- first_existing_column(
      annotation,
      c(
        "gene_id",
        "gene_id_stable",
        "ensembl_gene_id"
      )
    )

    symbol_col <- first_existing_column(
      annotation,
      c(
        "gene_symbol",
        "gene_name",
        "symbol"
      )
    )

    biotype_col <- first_existing_column(
      annotation,
      c(
        "gene_biotype",
        "gene_type",
        "biotype"
      )
    )

    required_mappings <- c(
      chromosome = chr_col,
      gene_start = start_col,
      gene_end = end_col,
      gene_id = id_col,
      gene_symbol = symbol_col,
      gene_biotype = biotype_col
    )

    if (anyNA(required_mappings)) {
      stop(
        paste(
          "The selected GENCODE table could not be mapped.",
          "",
          "Detected mappings:",
          paste(
            names(required_mappings),
            required_mappings,
            sep = " = ",
            collapse = "\n"
          ),
          "",
          "Available columns:",
          paste(names(annotation), collapse = ", "),
          sep = "\n"
        ),
        call. = FALSE
      )
    }

    strand_col <- first_existing_column(
      annotation,
      c("strand", "gene_strand")
    )

    genes <- annotation[
      ,
      .(
        chromosome = get(chr_col),
        gene_start = as.integer(get(start_col)),
        gene_end = as.integer(get(end_col)),
        strand = if (!is.na(strand_col)) {
          as.character(get(strand_col))
        } else {
          NA_character_
        },
        gene_id = as.character(get(id_col)),
        gene_symbol = as.character(get(symbol_col)),
        gene_biotype = as.character(get(biotype_col))
      )
    ]
  }

  genes[
    ,
    chromosome := standardize_chromosome(chromosome)
  ]

  genes[
    ,
    gene_id_stable := stable_ensembl_id(gene_id)
  ]

  genes <- genes[
    !is.na(chromosome) &
      !is.na(gene_start) &
      !is.na(gene_end) &
      !is.na(gene_id_stable)
  ]

  genes <- unique(
    genes,
    by = c(
      "chromosome",
      "gene_start",
      "gene_end",
      "gene_id_stable"
    )
  )

  genes
}

###############################################################################
# 6. Input validation
###############################################################################

message_section("MODULE 12.6 – INPUT VALIDATION")

assert_file_exists(
  clump_metadata_file,
  "Clump master metadata"
)

assert_file_exists(
  prioritization_file,
  "Module 12.5B prioritization table"
)

gencode_file <- find_gencode_annotation(project_dir)

cat("Clump metadata:\n", clump_metadata_file, "\n\n")
cat("Prioritization:\n", prioritization_file, "\n\n")
cat("GENCODE:\n", gencode_file, "\n")

###############################################################################
# 7. Read clump metadata
###############################################################################

message_section("READING CLUMP METADATA")

clumps <- fread(
  clump_metadata_file,
  showProgress = TRUE
)

assert_columns(
  clumps,
  c(
    "clump_id",
    "parent_locus",
    "chromosome",
    "index_variant",
    "index_position",
    "index_p"
  ),
  "Clump metadata"
)

clumps[
  ,
  chromosome := standardize_chromosome(chromosome)
]

clumps[
  ,
  index_position := as.integer(index_position)
]

clumps[
  ,
  index_p := as.numeric(index_p)
]

if (anyDuplicated(clumps$clump_id)) {
  stop(
    "Duplicated clump_id values were detected in clump metadata.",
    call. = FALSE
  )
}

cat("Total clumps:", nrow(clumps), "\n")
cat("Unique parent loci:", uniqueN(clumps$parent_locus), "\n")

###############################################################################
# 8. Read GENCODE genes
###############################################################################

message_section("READING GENCODE v50")

genes <- read_gencode_annotation(gencode_file)

protein_coding_genes <- genes[
  is_protein_coding(gene_biotype)
]

cat("All GENCODE genes:", nrow(genes), "\n")
cat("Protein-coding genes:", nrow(protein_coding_genes), "\n")

if (nrow(protein_coding_genes) == 0L) {
  stop(
    "No protein-coding genes were identified in the GENCODE annotation.",
    call. = FALSE
  )
}

###############################################################################
# 9. Conventional baseline assignment
###############################################################################

message_section("ASSIGNING CONVENTIONAL BASELINE GENES")

assign_baseline_gene <- function(clump_row, coding_genes) {

  chr <- clump_row$chromosome
  pos <- clump_row$index_position

  chr_genes <- coding_genes[
    chromosome == chr
  ]

  if (nrow(chr_genes) == 0L) {
    return(
      data.table(
        baseline_gene_id = NA_character_,
        baseline_gene_symbol = NA_character_,
        baseline_gene_biotype = NA_character_,
        baseline_assignment = "NO_CODING_GENE_ON_CHROMOSOME",
        baseline_distance_bp = NA_integer_,
        n_overlapping_coding_genes = 0L
      )
    )
  }

  overlapping <- chr_genes[
    gene_start <= pos &
      gene_end >= pos
  ]

  if (nrow(overlapping) > 0L) {

    overlapping[
      ,
      gene_length := gene_end - gene_start + 1L
    ]

    setorder(
      overlapping,
      gene_length,
      gene_id_stable
    )

    selected <- overlapping[1L]

    return(
      data.table(
        baseline_gene_id = selected$gene_id_stable,
        baseline_gene_symbol = selected$gene_symbol,
        baseline_gene_biotype = selected$gene_biotype,
        baseline_assignment = "OVERLAPPING_PROTEIN_CODING",
        baseline_distance_bp = 0L,
        n_overlapping_coding_genes = nrow(overlapping)
      )
    )
  }

  chr_genes[
    ,
    distance_bp := fifelse(
      pos < gene_start,
      gene_start - pos,
      pos - gene_end
    )
  ]

  setorder(
    chr_genes,
    distance_bp,
    gene_id_stable
  )

  selected <- chr_genes[1L]

  data.table(
    baseline_gene_id = selected$gene_id_stable,
    baseline_gene_symbol = selected$gene_symbol,
    baseline_gene_biotype = selected$gene_biotype,
    baseline_assignment = "NEAREST_PROTEIN_CODING",
    baseline_distance_bp = as.integer(selected$distance_bp),
    n_overlapping_coding_genes = 0L
  )
}

baseline_list <- vector(
  mode = "list",
  length = nrow(clumps)
)

for (i in seq_len(nrow(clumps))) {

  if (i %% 100L == 0L || i == 1L || i == nrow(clumps)) {
    message(
      "Baseline assignment: ",
      i,
      " / ",
      nrow(clumps)
    )
  }

  baseline_list[[i]] <- cbind(
    clumps[
      i,
      .(
        clump_id,
        parent_locus,
        chromosome,
        index_variant,
        index_position,
        index_p
      )
    ],
    assign_baseline_gene(
      clumps[i],
      protein_coding_genes
    )
  )
}

baseline <- rbindlist(
  baseline_list,
  use.names = TRUE,
  fill = TRUE
)

###############################################################################
# 10. Read Module 12.5B prioritization
###############################################################################

message_section("READING Module 12.5C gene-level prioritization")

priority <- fread(
  prioritization_file,
  showProgress = TRUE
)

cat("Priority rows:", nrow(priority), "\n")
cat("Priority columns:", ncol(priority), "\n")

clump_col <- first_existing_column(
  priority,
  c("clump_id", "clump")
)

gene_id_col <- first_existing_column(
  priority,
  c(
    "gene_id_stable",
    "gene_id",
    "ensembl_gene_id"
  )
)

gene_symbol_col <- first_existing_column(
  priority,
  c(
    "gene_symbol",
    "gene_name",
    "symbol"
  )
)

gene_biotype_col <- first_existing_column(
  priority,
  c(
    "gene_biotype",
    "gene_type",
    "biotype"
  )
)

candidate_class_col <- first_existing_column(
  priority,
  c(
    "candidate_class",
    "candidate_gene_class"
  )
)

priority_class_col <- first_existing_column(
  priority,
  c(
    "priority_class",
    "coloc_priority_class"
  )
)

priority_rank_col <- first_existing_column(
  priority,
  c(
    "priority_rank",
    "rank",
    "global_rank"
  )
)

candidate_score_col <- first_existing_column(
  priority,
  c(
    "candidate_score",
    "priority_score",
    "score"
  )
)

pph4_col <- first_existing_column(
  priority,
  c(
    "PP.H4",
    "PP.H4.abf",
    "pp_h4",
    "PP_H4"
  )
)

pph3_col <- first_existing_column(
  priority,
  c(
    "PP.H3",
    "PP.H3.abf",
    "pp_h3",
    "PP_H3"
  )
)

essential_mappings <- c(
  clump_id = clump_col,
  gene_id = gene_id_col,
  priority_class = priority_class_col,
  PP.H4 = pph4_col
)

if (anyNA(essential_mappings)) {
  stop(
    paste(
      "Could not map essential Module 12.5B columns.",
      "",
      "Detected mappings:",
      paste(
        names(essential_mappings),
        essential_mappings,
        sep = " = ",
        collapse = "\n"
      ),
      "",
      "Available columns:",
      paste(names(priority), collapse = ", "),
      sep = "\n"
    ),
    call. = FALSE
  )
}

priority_standard <- priority[
  ,
  .(
    clump_id = as.character(get(clump_col)),
    prioritized_gene_id = stable_ensembl_id(get(gene_id_col)),

    prioritized_gene_symbol =
      if (!is.na(gene_symbol_col)) {
        as.character(get(gene_symbol_col))
      } else {
        NA_character_
      },

    prioritized_gene_biotype =
      if (!is.na(gene_biotype_col)) {
        as.character(get(gene_biotype_col))
      } else {
        NA_character_
      },

    candidate_class =
      if (!is.na(candidate_class_col)) {
        as.character(get(candidate_class_col))
      } else {
        NA_character_
      },

    priority_class = as.character(get(priority_class_col)),

    source_priority_rank =
      if (!is.na(priority_rank_col)) {
        suppressWarnings(as.numeric(get(priority_rank_col)))
      } else {
        NA_real_
      },

    candidate_score =
      if (!is.na(candidate_score_col)) {
        suppressWarnings(as.numeric(get(candidate_score_col)))
      } else {
        NA_real_
      },

    PP.H4 = suppressWarnings(as.numeric(get(pph4_col))),

    PP.H3 =
      if (!is.na(pph3_col)) {
        suppressWarnings(as.numeric(get(pph3_col)))
      } else {
        NA_real_
      }
  )
]

priority_standard[
  ,
  priority_class_rank := priority_class_rank(priority_class)
]

priority_standard[
  ,
  prioritized_biotype_class := classify_candidate_biotype(
    prioritized_gene_biotype,
    candidate_class
  )
]

###############################################################################
# 11. Select best supported gene per clump
###############################################################################

best_gene <- priority_standard

cat(
  "Best genes loaded:",
  nrow(best_gene),
  "\n"
)

###############################################################################
# 12. Reclassification table
###############################################################################

message_section("BUILDING RECLASSIFICATION TABLE")

reclassification <- merge(
  baseline,
  best_gene,
  by = "clump_id",
  all.x = TRUE,
  sort = FALSE
)

reclassification[
  ,
  evaluable := !is.na(prioritized_gene_id)
]

reclassification[
  evaluable == TRUE,
  same_gene := baseline_gene_id == prioritized_gene_id
]

reclassification[
  evaluable == FALSE,
  same_gene := NA
]

reclassification[
  ,
  reclassified := fifelse(
    !evaluable,
    NA,
    !same_gene
  )
]

reclassification[
  evaluable == FALSE,
  reclassification_type := "NOT_EVALUABLE"
]

reclassification[
  evaluable == TRUE & same_gene == TRUE,
  reclassification_type := "UNCHANGED_PROTEIN_CODING"
]

reclassification[
  evaluable == TRUE &
    same_gene == FALSE &
    prioritized_biotype_class == "REGULATORY_RNA",
  reclassification_type := "CODING_TO_REGULATORY_RNA"
]

reclassification[
  evaluable == TRUE &
    same_gene == FALSE &
    prioritized_biotype_class == "PSEUDOGENE",
  reclassification_type := "CODING_TO_PSEUDOGENE"
]

reclassification[
  evaluable == TRUE &
    same_gene == FALSE &
    prioritized_biotype_class == "PROTEIN_CODING",
  reclassification_type := "CODING_TO_OTHER_PROTEIN_CODING"
]

reclassification[
  evaluable == TRUE &
    same_gene == FALSE &
    prioritized_biotype_class == "OTHER_OR_UNKNOWN",
  reclassification_type := "CODING_TO_OTHER_OR_UNKNOWN"
]

reclassification[
  ,
  interpretation_change := fifelse(
    reclassification_type %in% c(
      "CODING_TO_REGULATORY_RNA",
      "CODING_TO_PSEUDOGENE",
      "CODING_TO_OTHER_PROTEIN_CODING",
      "CODING_TO_OTHER_OR_UNKNOWN"
    ),
    TRUE,
    fifelse(
      reclassification_type == "UNCHANGED_PROTEIN_CODING",
      FALSE,
      NA
    )
  )
]

###############################################################################
# 13. Summary tables
###############################################################################

message_section("GENERATING SUMMARY TABLES")

overall_summary <- data.table(
  metric = c(
    "total_clumps",
    "unique_parent_loci",
    "evaluable_clumps",
    "non_evaluable_clumps",
    "unchanged_clumps",
    "reclassified_clumps"
  ),
  value = c(
    nrow(reclassification),
    uniqueN(reclassification$parent_locus),
    sum(reclassification$evaluable, na.rm = TRUE),
    sum(!reclassification$evaluable, na.rm = TRUE),
    sum(reclassification$interpretation_change == FALSE, na.rm = TRUE),
    sum(reclassification$interpretation_change == TRUE, na.rm = TRUE)
  )
)

n_evaluable <- sum(
  reclassification$evaluable,
  na.rm = TRUE
)

n_reclassified <- sum(
  reclassification$interpretation_change == TRUE,
  na.rm = TRUE
)

overall_summary <- rbind(
  overall_summary,
  data.table(
    metric = "reclassified_percent_of_evaluable",
    value = if (n_evaluable > 0L) {
      100 * n_reclassified / n_evaluable
    } else {
      NA_real_
    }
  )
)

type_summary <- reclassification[
  ,
  .N,
  by = reclassification_type
][
  order(-N)
]

type_summary <- reclassification[
  ,
  .N,
  by = reclassification_type
][
  order(-N)
]

type_summary[
  ,
  denominator := ifelse(
    reclassification_type == "NOT_EVALUABLE",
    nrow(reclassification),
    n_evaluable
  )
]

type_summary[
  ,
  percent := round(100 * N / denominator, 2)
]

type_summary[
  ,
  denominator := NULL
]

setnames(
  type_summary,
  "percent",
  "percent_of_reference"
)

priority_summary <- reclassification[
  evaluable == TRUE,
  .(
    n_clumps = .N,
    n_reclassified = sum(
      interpretation_change == TRUE,
      na.rm = TRUE
    ),
    percent_reclassified = 100 * mean(
      interpretation_change == TRUE,
      na.rm = TRUE
    )
  ),
  by = priority_class
][
  order(priority_class_rank(priority_class))
]

baseline_method_summary <- reclassification[
  ,
  .N,
  by = baseline_assignment
][
  order(-N)
]

reclassified_examples <- reclassification[
  interpretation_change == TRUE
][
  order(
    priority_class_rank(priority_class),
    -PP.H4,
    baseline_distance_bp
  )
]

###############################################################################
# 14. Parent-locus summary
###############################################################################

locus_summary <- reclassification[
  ,
  .(
    n_clumps = .N,
    n_evaluable_clumps = sum(evaluable, na.rm = TRUE),
    n_reclassified_clumps = sum(
      interpretation_change == TRUE,
      na.rm = TRUE
    ),
    any_reclassification = any(
      interpretation_change == TRUE,
      na.rm = TRUE
    ),
    reclassification_types = paste(
      sort(
        unique(
          reclassification_type[
            reclassification_type != "NOT_EVALUABLE"
          ]
        )
      ),
      collapse = ";"
    )
  ),
  by = parent_locus
]

locus_summary[
  n_evaluable_clumps == 0L,
  any_reclassification := NA
]

###############################################################################
# 15. QC tables
###############################################################################

unmatched_priority_clumps <- unique(
  priority_standard[
    !clump_id %in% clumps$clump_id,
    .(clump_id)
  ]
)

duplicate_best_gene_check <- best_gene[
  ,
  .N,
  by = clump_id
][
  N > 1L
]

###############################################################################
# 16. Write outputs
###############################################################################

message_section("WRITING OUTPUTS")

fwrite(
  reclassification,
  file.path(
    output_dir,
    "module12_6_clump_gene_reclassification.tsv.gz"
  ),
  sep = "\t"
)

fwrite(
  overall_summary,
  file.path(
    output_dir,
    "module12_6_reclassification_overall_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  type_summary,
  file.path(
    output_dir,
    "module12_6_reclassification_type_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  priority_summary,
  file.path(
    output_dir,
    "module12_6_reclassification_by_priority_class.tsv"
  ),
  sep = "\t"
)

fwrite(
  baseline_method_summary,
  file.path(
    output_dir,
    "module12_6_baseline_assignment_summary.tsv"
  ),
  sep = "\t"
)

fwrite(
  reclassified_examples,
  file.path(
    output_dir,
    "module12_6_reclassified_clumps.tsv"
  ),
  sep = "\t"
)

fwrite(
  locus_summary,
  file.path(
    output_dir,
    "module12_6_parent_locus_reclassification.tsv"
  ),
  sep = "\t"
)

fwrite(
  baseline,
  file.path(
    output_dir,
    "module12_6_conventional_gene_assignments.tsv.gz"
  ),
  sep = "\t"
)

fwrite(
  best_gene,
  file.path(
    output_dir,
    "module12_6_best_supported_gene_per_clump.tsv.gz"
  ),
  sep = "\t"
)

fwrite(
  unmatched_priority_clumps,
  file.path(
    output_dir,
    "module12_6_unmatched_priority_clumps.tsv"
  ),
  sep = "\t"
)

###############################################################################
# 17. Status and metadata
###############################################################################

status <- data.table(
  module = "Module 12.6",
  analysis = "Colocalization-driven Gene Reclassification",
  status = "SUCCESS",
  total_clumps = nrow(reclassification),
  evaluable_clumps = n_evaluable,
  reclassified_clumps = n_reclassified,
  reclassified_percent = if (n_evaluable > 0L) {
    100 * n_reclassified / n_evaluable
  } else {
    NA_real_
  },
  gencode_file = gencode_file,
  prioritization_file = prioritization_file,
  generated_at = format(
    Sys.time(),
    "%Y-%m-%d %H:%M:%S"
  )
)

fwrite(
  status,
  file.path(
    output_dir,
    "module12_6_status.tsv"
  ),
  sep = "\t"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    output_dir,
    "module12_6_sessionInfo.txt"
  )
)

###############################################################################
# 18. Console report
###############################################################################

message_section("MODULE 12.6 COMPLETED")

cat("Total clumps:       ", nrow(reclassification), "\n")
cat("Evaluable clumps:   ", n_evaluable, "\n")
cat("Reclassified clumps:", n_reclassified, "\n")

if (n_evaluable > 0L) {
  cat(
    "Reclassified:       ",
    sprintf(
      "%.2f%%",
      100 * n_reclassified / n_evaluable
    ),
    "\n"
  )
}

cat("\nReclassification types:\n")
print(type_summary)

cat("\nOutputs:\n")
cat(output_dir, "\n")