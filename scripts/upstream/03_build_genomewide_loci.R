library(data.table)

# ============================================================
# PigmentRegAtlas
# Module 03: Build genome-wide loci from PLINK clumps
# ============================================================

clump_file <- paste0(
  "results/",
  "GCST90691600_clump_uniqueID_r2_0.1_kb1000.clumps"
)

output_clumps <- "results/GCST90691600_clump_intervals.tsv"
output_loci <- "results/GCST90691600_genomic_loci.tsv"
output_qc <- "results/GCST90691600_locus_QC.tsv"

padding_bp <- 250000L
merge_distance_bp <- 250000L

cat("Reading PLINK clumps:\n", clump_file, "\n")

clumps <- fread(
  clump_file,
  sep = "\t",
  header = TRUE,
  fill = TRUE,
  quote = ""
)

cat("Clumps read:", nrow(clumps), "\n")
cat("Columns detected:\n")
print(names(clumps))

required_columns <- c(
  "#CHROM",
  "POS",
  "ID",
  "P",
  "TOTAL",
  "SP2"
)

missing_columns <- setdiff(
  required_columns,
  names(clumps)
)

if (length(missing_columns) > 0L) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

setnames(
  clumps,
  old = c("#CHROM", "POS", "ID", "P"),
  new = c("chr", "lead_pos", "lead_id", "lead_p")
)

# ------------------------------------------------------------
# Extract genomic position from IDs such as:
# 16:89246170:A:T
# ------------------------------------------------------------

extract_position <- function(ids) {

  if (
    length(ids) == 0L ||
    all(is.na(ids)) ||
    all(ids == "") ||
    all(ids == ".")
  ) {
    return(integer())
  }

  ids <- ids[
    !is.na(ids) &
      ids != "" &
      ids != "."
  ]

  if (length(ids) == 0L) {
    return(integer())
  }

  fields <- strsplit(
    ids,
    ":",
    fixed = TRUE
  )

  positions <- suppressWarnings(
    as.integer(
      vapply(
        fields,
        function(x) {
          if (length(x) >= 2L) {
            x[2L]
          } else {
            NA_character_
          }
        },
        character(1)
      )
    )
  )

  positions[
    !is.na(positions)
  ]
}

# ------------------------------------------------------------
# Determine the interval represented by each PLINK clump
# ------------------------------------------------------------

clump_intervals_list <- vector(
  mode = "list",
  length = nrow(clumps)
)

for (i in seq_len(nrow(clumps))) {

  lead_position <- as.integer(
    clumps$lead_pos[i]
  )

  sp2_text <- clumps$SP2[i]

  member_ids <- character()

  if (
    !is.na(sp2_text) &&
      sp2_text != "" &&
      sp2_text != "."
  ) {

    member_ids <- trimws(
      unlist(
        strsplit(
          sp2_text,
          ",",
          fixed = TRUE
        )
      )
    )
  }

  member_positions <- extract_position(
    member_ids
  )

  all_positions <- c(
    lead_position,
    member_positions
  )

  all_positions <- all_positions[
    !is.na(all_positions)
  ]

  if (length(all_positions) == 0L) {
    warning(
      "No valid positions found for clump row ",
      i,
      ". Skipping."
    )
    next
  }

  raw_start <- min(all_positions)
  raw_end <- max(all_positions)

  clump_intervals_list[[i]] <- data.table(
    clump_id = sprintf(
      "CLUMP_%04d",
      i
    ),
    chr = as.character(
      clumps$chr[i]
    ),
    lead_id = as.character(
      clumps$lead_id[i]
    ),
    lead_pos = lead_position,
    lead_p = as.numeric(
      clumps$lead_p[i]
    ),
    plink_total = as.integer(
      clumps$TOTAL[i]
    ),
    n_ids_parsed = length(
      all_positions
    ),
    raw_start = raw_start,
    raw_end = raw_end,
    padded_start = max(
      1L,
      raw_start - padding_bp
    ),
    padded_end = raw_end + padding_bp
  )
}

clump_intervals <- rbindlist(
  clump_intervals_list,
  fill = TRUE
)

if (nrow(clump_intervals) == 0L) {
  stop("No valid clump intervals were created.")
}

clump_intervals[, chr_num := suppressWarnings(
  as.integer(chr)
)]

if (any(is.na(clump_intervals$chr_num))) {
  warning(
    "Some chromosome values could not be converted to integers."
  )
}

setorder(
  clump_intervals,
  chr_num,
  padded_start,
  padded_end,
  lead_p
)

clump_intervals[, chr_num := NULL]

fwrite(
  clump_intervals,
  output_clumps,
  sep = "\t"
)

cat(
  "Clump interval table saved to:\n",
  output_clumps,
  "\n"
)

# ------------------------------------------------------------
# Merge overlapping or nearby clump intervals
# ------------------------------------------------------------

merged_loci <- list()
locus_counter <- 0L

chromosomes <- unique(
  clump_intervals$chr
)

for (chromosome in chromosomes) {

  x <- clump_intervals[
    chr == chromosome
  ]

  if (nrow(x) == 0L) {
    next
  }

  current_start <- x$padded_start[1]
  current_end <- x$padded_end[1]
  member_rows <- 1L

  save_locus <- function(
    rows,
    start_pos,
    end_pos
  ) {

    locus_counter <<- locus_counter + 1L

    members <- x[rows]

    lead_row_index <- which.min(
      members$lead_p
    )

    lead_row <- members[
      lead_row_index
    ]

    merged_loci[[locus_counter]] <<- data.table(
      locus_id = sprintf(
        "LOCUS_%04d_chr%s",
        locus_counter,
        chromosome
      ),
      chr = chromosome,
      start = as.integer(
        start_pos
      ),
      end = as.integer(
        end_pos
      ),
      locus_size_bp = as.integer(
        end_pos - start_pos + 1L
      ),
      lead_id = lead_row$lead_id,
      lead_pos = lead_row$lead_pos,
      lead_p = lead_row$lead_p,
      n_clumps = nrow(
        members
      ),
      n_clump_variants = sum(
        members$plink_total,
        na.rm = TRUE
      ),
      member_clumps = paste(
        members$clump_id,
        collapse = ";"
      ),
      member_leads = paste(
        members$lead_id,
        collapse = ";"
      )
    )
  }

  if (nrow(x) > 1L) {

    for (i in 2:nrow(x)) {

      next_start <- x$padded_start[i]
      next_end <- x$padded_end[i]

      if (
        next_start <=
          current_end + merge_distance_bp
      ) {

        current_end <- max(
          current_end,
          next_end
        )

        member_rows <- c(
          member_rows,
          i
        )

      } else {

        save_locus(
          member_rows,
          current_start,
          current_end
        )

        current_start <- next_start
        current_end <- next_end
        member_rows <- i
      }
    }
  }

  save_locus(
    member_rows,
    current_start,
    current_end
  )
}

loci <- rbindlist(
  merged_loci,
  fill = TRUE
)

if (nrow(loci) == 0L) {
  stop("No genomic loci were created.")
}

loci[, chr_num := suppressWarnings(
  as.integer(chr)
)]

setorder(
  loci,
  chr_num,
  start,
  end
)

loci[, chr_num := NULL]

# Renumber after genomic sorting
loci[, locus_id := sprintf(
  "LOCUS_%04d_chr%s",
  .I,
  chr
)]

fwrite(
  loci,
  output_loci,
  sep = "\t"
)

# ------------------------------------------------------------
# QC summary
# ------------------------------------------------------------

qc <- data.table(
  parameter = c(
    "input_clumps",
    "final_genomic_loci",
    "padding_bp",
    "merge_distance_bp",
    "median_clumps_per_locus",
    "maximum_clumps_per_locus",
    "median_locus_size_bp",
    "maximum_locus_size_bp"
  ),
  value = c(
    nrow(clump_intervals),
    nrow(loci),
    padding_bp,
    merge_distance_bp,
    median(
      loci$n_clumps,
      na.rm = TRUE
    ),
    max(
      loci$n_clumps,
      na.rm = TRUE
    ),
    median(
      loci$locus_size_bp,
      na.rm = TRUE
    ),
    max(
      loci$locus_size_bp,
      na.rm = TRUE
    )
  )
)

fwrite(
  qc,
  output_qc,
  sep = "\t"
)

cat("\nLocus QC summary:\n")
print(qc)

cat(
  "\nFinal locus table saved to:\n",
  output_loci,
  "\n"
)

cat(
  "Locus QC saved to:\n",
  output_qc,
  "\n"
)

cat("\nTop 10 loci by significance:\n")

top_n <- min(
  10L,
  nrow(loci)
)

print(
  loci[
    order(lead_p)
  ][1:top_n]
)