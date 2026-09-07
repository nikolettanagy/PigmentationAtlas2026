#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

cat("============================================================\n")
cat("PigmentationAtlas\n")
cat("MODULE 12.5B – GENCODE v50 ANNOTATION MERGE\n")
cat("============================================================\n\n")

start_time <- Sys.time()
project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# ------------------------------------------------------------------
# Command-line argument parser
# Optional:
#   --gencode path/to/gencode.v50.annotation.gtf.gz
#   --input   path/to/module12_5_candidate_gene_prioritization.tsv.gz
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(NULL)
  args[[idx + 1L]]
}

input_override <- get_arg("--input")
gencode_override <- get_arg("--gencode")

default_input <- file.path(
  project_dir,
  "results", "module12", "colocalization", "biological_annotation",
  "tables", "module12_5_candidate_gene_prioritization.tsv.gz"
)

input_file <- if (!is.null(input_override)) {
  normalizePath(input_override, winslash = "/", mustWork = FALSE)
} else {
  default_input
}

outdir <- file.path(
  project_dir,
  "results", "module12", "colocalization", "biological_annotation"
)
tables_dir <- file.path(outdir, "tables")
logs_dir <- file.path(outdir, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  logs_dir,
  paste0(
    "module12_5B_gencode_annotation_",
    format(start_time, "%Y%m%d_%H%M%S"),
    ".log"
  )
)

log_message <- function(...) {
  msg <- paste0(..., collapse = "")
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

fail <- function(msg) {
  log_message("FATAL: ", msg)
  cat("\nFINAL RESULT: FAIL\n")
  quit(save = "no", status = 1L)
}

if (!file.exists(input_file)) {
  fail(paste0("Missing input candidate catalogue: ", input_file))
}

# ------------------------------------------------------------------
# 1. Locate GENCODE v50 annotation source
# ------------------------------------------------------------------

candidate_roots <- unique(c(
  file.path(project_dir, "reference"),
  file.path(project_dir, "reference", "gencode"),
  file.path(project_dir, "reference", "GENCODE"),
  file.path(project_dir, "data"),
  file.path(project_dir, "results", "module04"),
  file.path(project_dir, "results", "module4"),
  file.path(project_dir, "results", "gene_annotation")
))

candidate_roots <- candidate_roots[dir.exists(candidate_roots)]

find_annotation_files <- function(roots) {
  if (length(roots) == 0L) return(character())
  files <- unlist(
    lapply(
      roots,
      function(x) {
        list.files(
          x,
          recursive = TRUE,
          full.names = TRUE,
          pattern = "\\.(gtf|gtf\\.gz|gff|gff3|gff\\.gz|gff3\\.gz|tsv|tsv\\.gz|txt|txt\\.gz)$",
          ignore.case = TRUE
        )
      }
    ),
    use.names = FALSE
  )
  unique(files)
}

if (!is.null(gencode_override)) {
  annotation_source <- normalizePath(
    gencode_override,
    winslash = "/",
    mustWork = FALSE
  )
  if (!file.exists(annotation_source)) {
    fail(paste0("Specified --gencode file does not exist: ", annotation_source))
  }
} else {
  all_candidates <- find_annotation_files(candidate_roots)

  if (length(all_candidates) == 0L) {
    fail(
      paste0(
        "No GENCODE annotation file found under reference/, data/, or Module 04 results. ",
        "Run again with --gencode <path>."
      )
    )
  }

  score_file <- function(x) {
    b <- tolower(basename(x))
    score <- 0L
    if (grepl("gencode", b, fixed = TRUE)) score <- score + 10L
    if (grepl("v50", b, fixed = TRUE)) score <- score + 10L
    if (grepl("annotation", b, fixed = TRUE)) score <- score + 4L
    if (grepl("\\.gtf(\\.gz)?$", b)) score <- score + 8L
    if (grepl("\\.gff3?(\\.gz)?$", b)) score <- score + 6L
    if (grepl("gene", b, fixed = TRUE)) score <- score + 3L
    score
  }

  scores <- vapply(all_candidates, score_file, integer(1L))
  ranked <- all_candidates[order(scores, decreasing = TRUE)]

  annotation_source <- ranked[[1L]]

  fwrite(
    data.table(
      candidate_file = ranked,
      discovery_score = scores[order(scores, decreasing = TRUE)]
    ),
    file.path(tables_dir, "module12_5B_annotation_source_candidates.tsv"),
    sep = "\t",
    quote = FALSE
  )
}

log_message("Candidate catalogue: ", input_file)
log_message("Annotation source: ", annotation_source)

# ------------------------------------------------------------------
# 2. Read candidate catalogue
# ------------------------------------------------------------------

dt <- fread(input_file, showProgress = FALSE)

if (!"gene_id_stable" %in% names(dt)) {
  if (!"gene_id" %in% names(dt)) {
    fail("Candidate catalogue contains neither gene_id_stable nor gene_id.")
  }
  dt[, gene_id_stable := sub("\\.[0-9]+$", "", as.character(gene_id))]
}

dt[, gene_id_stable := sub("\\.[0-9]+$", "", as.character(gene_id_stable))]

# Remove placeholder annotation columns before merge.
drop_cols <- intersect(
  c(
    "gene_symbol", "gene_name", "gene_biotype",
    "gene_chromosome", "gene_start", "gene_end", "gene_strand",
    "candidate_class"
  ),
  names(dt)
)

if (length(drop_cols) > 0L) {
  dt[, (drop_cols) := NULL]
}

# ------------------------------------------------------------------
# 3. Parse GTF/GFF or tabular gene annotation
# ------------------------------------------------------------------

parse_gtf_attributes <- function(x, key) {
  pattern <- paste0("(?:^|;[[:space:]]*)", key, "[[:space:]]+\"([^\"]+)\"")
  out <- sub(paste0(".*", pattern, ".*"), "\\1", x, perl = TRUE)
  miss <- !grepl(pattern, x, perl = TRUE)
  out[miss] <- NA_character_
  out
}

parse_gff_attributes <- function(x, keys) {
  out <- rep(NA_character_, length(x))
  for (key in keys) {
    pattern <- paste0("(?:^|;)", key, "=([^;]+)")
    hit <- grepl(pattern, x, perl = TRUE)
    tmp <- sub(paste0(".*", pattern, ".*"), "\\1", x, perl = TRUE)
    out[is.na(out) & hit] <- tmp[is.na(out) & hit]
  }
  out
}

read_gencode_gtf <- function(path) {
  log_message("Reading GTF/GFF annotation...")

  gtf <- fread(
    path,
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "#",
    fill = TRUE,
    showProgress = FALSE,
    col.names = c(
      "seqname", "source", "feature", "start", "end",
      "score", "strand", "frame", "attribute"
    )
  )

  gtf <- gtf[feature == "gene"]

  if (nrow(gtf) == 0L) {
    stop("No gene features found in GTF/GFF source.", call. = FALSE)
  }

  is_gff <- grepl("\\.gff3?(\\.gz)?$", path, ignore.case = TRUE)

  if (is_gff) {
    gene_id <- parse_gff_attributes(gtf$attribute, c("gene_id", "ID"))
    gene_name <- parse_gff_attributes(
      gtf$attribute,
      c("gene_name", "Name", "gene")
    )
    gene_type <- parse_gff_attributes(
      gtf$attribute,
      c("gene_type", "gene_biotype", "biotype")
    )
  } else {
    gene_id <- parse_gtf_attributes(gtf$attribute, "gene_id")
    gene_name <- parse_gtf_attributes(gtf$attribute, "gene_name")
    gene_type <- parse_gtf_attributes(gtf$attribute, "gene_type")

    missing_type <- is.na(gene_type)
    if (any(missing_type)) {
      gene_type2 <- parse_gtf_attributes(gtf$attribute, "gene_biotype")
      gene_type[missing_type] <- gene_type2[missing_type]
    }
  }

  data.table(
    gene_id_stable = sub("\\.[0-9]+$", "", gene_id),
    gene_id_gencode = gene_id,
    gene_symbol = gene_name,
    gene_name = gene_name,
    gene_biotype = gene_type,
    gene_chromosome = as.character(gtf$seqname),
    gene_start = as.integer(gtf$start),
    gene_end = as.integer(gtf$end),
    gene_strand = as.character(gtf$strand)
  )
}

choose_column <- function(nms, candidates) {
  hit <- candidates[candidates %chin% nms]
  if (length(hit) == 0L) return(NA_character_)
  hit[[1L]]
}

read_tabular_annotation <- function(path) {
  log_message("Reading tabular annotation...")

  x <- fread(path, showProgress = FALSE)
  nms <- names(x)

  id_col <- choose_column(
    nms,
    c("gene_id_stable", "gene_id", "ensembl_gene_id", "ensembl_id")
  )
  symbol_col <- choose_column(
    nms,
    c("gene_symbol", "gene_name", "symbol", "hgnc_symbol")
  )
  biotype_col <- choose_column(
    nms,
    c("gene_biotype", "gene_type", "biotype")
  )
  chr_col <- choose_column(
    nms,
    c("gene_chromosome", "chromosome", "chr", "seqname")
  )
  start_col <- choose_column(
    nms,
    c("gene_start", "start", "start_position")
  )
  end_col <- choose_column(
    nms,
    c("gene_end", "end", "end_position")
  )
  strand_col <- choose_column(
    nms,
    c("gene_strand", "strand")
  )

  if (is.na(id_col)) {
    stop(
      paste0(
        "Tabular annotation lacks a recognizable Ensembl gene ID column. ",
        "Columns: ", paste(nms, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  out <- data.table(
    gene_id_stable = sub("\\.[0-9]+$", "", as.character(x[[id_col]])),
    gene_id_gencode = as.character(x[[id_col]]),
    gene_symbol = if (!is.na(symbol_col)) as.character(x[[symbol_col]]) else NA_character_,
    gene_name = if (!is.na(symbol_col)) as.character(x[[symbol_col]]) else NA_character_,
    gene_biotype = if (!is.na(biotype_col)) as.character(x[[biotype_col]]) else NA_character_,
    gene_chromosome = if (!is.na(chr_col)) as.character(x[[chr_col]]) else NA_character_,
    gene_start = if (!is.na(start_col)) as.integer(x[[start_col]]) else NA_integer_,
    gene_end = if (!is.na(end_col)) as.integer(x[[end_col]]) else NA_integer_,
    gene_strand = if (!is.na(strand_col)) as.character(x[[strand_col]]) else NA_character_
  )

  out
}

annotation <- tryCatch(
  {
    if (grepl("\\.(gtf|gtf\\.gz|gff|gff3|gff\\.gz|gff3\\.gz)$",
              annotation_source, ignore.case = TRUE)) {
      read_gencode_gtf(annotation_source)
    } else {
      read_tabular_annotation(annotation_source)
    }
  },
  error = function(e) {
    fail(paste0("Annotation parsing failed: ", conditionMessage(e)))
  }
)

annotation <- annotation[
  !is.na(gene_id_stable) & gene_id_stable != ""
]

annotation <- unique(annotation, by = "gene_id_stable")

log_message("GENCODE gene records: ", nrow(annotation))

# ------------------------------------------------------------------
# 4. Merge annotation
# ------------------------------------------------------------------

annotated <- merge(
  dt,
  annotation,
  by = "gene_id_stable",
  all.x = TRUE,
  sort = FALSE
)

annotated[, candidate_class := fifelse(
  is.na(gene_biotype) | gene_biotype == "",
  "Unknown",
  fifelse(
    gene_biotype == "protein_coding",
    "Protein_coding",
    fifelse(
      grepl(
        "lncRNA|lincRNA|antisense|sense_intronic|sense_overlapping|processed_transcript|macro_lncRNA|bidirectional_promoter_lncRNA",
        gene_biotype,
        ignore.case = TRUE
      ),
      "Regulatory_RNA",
      fifelse(
        grepl("pseudogene", gene_biotype, ignore.case = TRUE),
        "Pseudogene",
        fifelse(
          grepl("miRNA|snoRNA|snRNA|rRNA|scaRNA|vault_RNA|ribozyme",
                gene_biotype, ignore.case = TRUE),
          "Small_or_structural_RNA",
          "Other"
        )
      )
    )
  )
)]

if ("candidate_score" %in% names(annotated)) {
  annotated[
    candidate_class == "Regulatory_RNA",
    candidate_score := round(candidate_score + 5, 3)
  ]
}

if (all(c("priority_rank", "candidate_score", "PP.H4", "global_rank") %in% names(annotated))) {
  setorder(
    annotated,
    priority_rank,
    -candidate_score,
    -PP.H4,
    global_rank
  )
}

# ------------------------------------------------------------------
# 5. Output tables
# ------------------------------------------------------------------

write_tsv <- function(x, filename, gzip = FALSE) {
  path <- file.path(tables_dir, filename)
  fwrite(
    x,
    path,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = if (gzip) "gzip" else "none"
  )
  cat("[1] \"", normalizePath(path, winslash = "/", mustWork = FALSE), "\"\n", sep = "")
  invisible(path)
}

full_file <- write_tsv(
  annotated,
  "module12_5B_candidate_gene_prioritization_GENCODEv50.tsv.gz",
  gzip = TRUE
)

strong <- annotated[priority_class == "STRONG_H4"]
strong_moderate <- annotated[
  priority_class %chin% c("STRONG_H4", "MODERATE_H4")
]
lncrna <- annotated[candidate_class == "Regulatory_RNA"]
protein_coding <- annotated[candidate_class == "Protein_coding"]
pseudogene <- annotated[candidate_class == "Pseudogene"]

write_tsv(strong, "module12_5B_strong_H4_candidates_GENCODEv50.tsv")
write_tsv(
  strong_moderate,
  "module12_5B_strong_and_moderate_candidates_GENCODEv50.tsv"
)
write_tsv(lncrna, "module12_5B_regulatory_RNA_candidates.tsv")
write_tsv(protein_coding, "module12_5B_protein_coding_candidates.tsv")
write_tsv(pseudogene, "module12_5B_pseudogene_candidates.tsv")

# Gene-level summary
safe_max <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

safe_first_character <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0L) return(NA_character_)
  as.character(x[[1L]])
}

safe_first_integer <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_integer_)
  as.integer(x[[1L]])
}


gene_summary <- annotated[
  ,
  .(
    gene_id_gencode = safe_first_character(gene_id_gencode),
    gene_symbol = safe_first_character(gene_symbol),
    gene_name = safe_first_character(gene_name),
    gene_biotype = safe_first_character(gene_biotype),
    candidate_class = first(candidate_class),
    gene_chromosome = safe_first_character(gene_chromosome),
    gene_start = safe_first_integer(gene_start),
    gene_end = safe_first_integer(gene_end),
    gene_strand = safe_first_character(gene_strand),
    n_coloc_units = .N,
    n_clumps = uniqueN(clump_id),
    n_tissues = uniqueN(tissue),
    n_phenotypes = uniqueN(phenotype_id),
    best_PP.H4 = safe_max(PP.H4),
    best_PP.H3 = safe_max(PP.H3),
    best_priority_rank = min(priority_rank, na.rm = TRUE),
    best_priority_class = priority_class[
      order(priority_rank, -PP.H4)
    ][1L],
    best_coloc_unit_id = coloc_unit_id[
      order(priority_rank, -PP.H4)
    ][1L]
  ),
  by = gene_id_stable
]

# Fix zero-length first() cases safely
for (col in c(
  "gene_id_gencode", "gene_symbol", "gene_name", "gene_biotype",
  "gene_chromosome", "gene_start", "gene_end", "gene_strand"
)) {
  if (!col %in% names(gene_summary)) next
}

setorder(gene_summary, best_priority_rank, -best_PP.H4)
gene_summary[, gene_rank := seq_len(.N)]
setcolorder(gene_summary, c("gene_rank", setdiff(names(gene_summary), "gene_rank")))

write_tsv(gene_summary, "module12_5B_gene_priority_summary_GENCODEv50.tsv")

# Biotype counts
biotype_counts <- annotated[
  ,
  .(
    n_units = .N,
    n_genes = uniqueN(gene_id_stable),
    n_clumps = uniqueN(clump_id),
    median_PP.H4 = median(PP.H4, na.rm = TRUE),
    max_PP.H4 = max(PP.H4, na.rm = TRUE)
  ),
  by = .(priority_class, candidate_class, gene_biotype)
][order(priority_rank <- match(
  priority_class,
  c("STRONG_H4", "MODERATE_H4", "SUGGESTIVE_H4", "H3_DOMINANT", "INCONCLUSIVE")
), candidate_class, gene_biotype)]

biotype_counts[, priority_rank := NULL]
write_tsv(biotype_counts, "module12_5B_biotype_counts.tsv")

# Strong-H4 unique-gene summary
strong_gene_summary <- gene_summary[
  gene_id_stable %chin% unique(strong$gene_id_stable)
]
write_tsv(
  strong_gene_summary,
  "module12_5B_strong_H4_unique_genes_GENCODEv50.tsv"
)

# ------------------------------------------------------------------
# 6. QC
# ------------------------------------------------------------------

n_annotated_units <- annotated[
  !is.na(gene_symbol) | !is.na(gene_biotype),
  .N
]
n_annotated_genes <- uniqueN(
  annotated[!is.na(gene_symbol) | !is.na(gene_biotype), gene_id_stable]
)
n_missing_genes <- uniqueN(
  annotated[is.na(gene_symbol) & is.na(gene_biotype), gene_id_stable]
)

missing_annotation <- unique(
  annotated[
    is.na(gene_symbol) & is.na(gene_biotype),
    .(gene_id_stable)
  ]
)

write_tsv(
  missing_annotation,
  "module12_5B_unmatched_Ensembl_gene_ids.tsv"
)

critical_failures <- 0L
warnings <- 0L

if (nrow(annotated) != nrow(dt)) critical_failures <- critical_failures + 1L
if (anyDuplicated(annotated$coloc_unit_id) > 0L) critical_failures <- critical_failures + 1L
if (n_missing_genes > 0L) warnings <- warnings + 1L
if (nrow(strong) != 25L) warnings <- warnings + 1L

final_status <- if (critical_failures > 0L) {
  "FAIL"
} else if (warnings > 0L) {
  "PASS_WITH_WARNINGS"
} else {
  "PASS"
}

end_time <- Sys.time()

status <- data.table(
  module = "12.5B",
  start_time = format(start_time, "%Y-%m-%d %H:%M:%S"),
  end_time = format(end_time, "%Y-%m-%d %H:%M:%S"),
  runtime_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),
  annotation_source = annotation_source,
  input_coloc_units = nrow(dt),
  output_coloc_units = nrow(annotated),
  unique_genes = uniqueN(annotated$gene_id_stable),
  annotated_units = n_annotated_units,
  annotated_genes = n_annotated_genes,
  unmatched_genes = n_missing_genes,
  strong_h4_units = nrow(strong),
  strong_h4_unique_genes = uniqueN(strong$gene_id_stable),
  regulatory_RNA_units = nrow(lncrna),
  regulatory_RNA_genes = uniqueN(lncrna$gene_id_stable),
  protein_coding_units = nrow(protein_coding),
  protein_coding_genes = uniqueN(protein_coding$gene_id_stable),
  pseudogene_units = nrow(pseudogene),
  pseudogene_genes = uniqueN(pseudogene$gene_id_stable),
  critical_failures = critical_failures,
  warnings = warnings,
  final_status = final_status
)

write_tsv(status, "module12_5B_status.tsv")

cat("\n============================================================\n")
cat("MODULE 12.5B – FINAL SUMMARY\n")
cat("============================================================\n\n")
print(status)

cat("\nBiotype distribution among Strong H4 units:\n")
print(
  strong[
    ,
    .(
      n_units = .N,
      n_genes = uniqueN(gene_id_stable)
    ),
    by = .(candidate_class, gene_biotype)
  ][order(-n_units)]
)

cat("\n============================================================\n")
cat("FINAL RESULT: ", final_status, "\n", sep = "")
cat("============================================================\n")

quit(
  save = "no",
  status = if (identical(final_status, "FAIL")) 1L else 0L
)
