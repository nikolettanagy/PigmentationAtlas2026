suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# PigmentRegAtlas
# Module 03 v3: Lead-position-based genomic locus construction
#
# Core principle:
#   - Sort PLINK clumps by chromosome and lead position
#   - Start a new locus when the gap between consecutive
#     clump lead variants exceeds a configurable threshold
#   - Define primary locus boundaries from lead positions
#   - Retain full raw clump support as separate QC information
# ============================================================

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

trait_id <- "GCST90691600"

lead_gap_threshold_bp <- 250000L
locus_padding_bp <- 50000L

input_file <- file.path(
  "results",
  paste0(trait_id, "_clump_intervals.tsv")
)

output_clumps <- file.path(
  "results",
  paste0(trait_id, "_clump_intervals_v3.tsv")
)

output_loci <- file.path(
  "results",
  paste0(trait_id, "_genomic_loci_v3.tsv")
)

output_membership <- file.path(
  "results",
  paste0(trait_id, "_locus_clump_membership_v3.tsv")
)

output_qc <- file.path(
  "results",
  paste0(trait_id, "_locus_QC_v3.tsv")
)

output_largest <- file.path(
  "results",
  paste0(trait_id, "_largest_loci_v3.tsv")
)

output_histogram <- file.path(
  "results",
  paste0(trait_id, "_locus_size_histogram_v3.png")
)

output_gap_histogram <- file.path(
  "results",
  paste0(trait_id, "_lead_gap_histogram_v3.png")
)

output_tpcn2 <- file.path(
  "results",
  paste0(trait_id, "_TPCN2_region_loci_v3.tsv")
)

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

collapse_unique <- function(x) {

  x <- unique(
    x[
      !is.na(x) &
        x != ""
    ]
  )

  if (length(x) == 0L) {
    return(NA_character_)
  }

  paste(
    x,
    collapse = ";"
  )
}

chromosome_order <- function(chr) {

  chr <- gsub(
    "^chr",
    "",
    as.character(chr),
    ignore.case = TRUE
  )

  result <- suppressWarnings(
    as.integer(chr)
  )

  result[chr == "X"] <- 23L
  result[chr == "Y"] <- 24L
  result[chr %in% c("M", "MT")] <- 25L

  result
}

safe_quantile <- function(x, probability) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  as.numeric(
    quantile(
      x,
      probability,
      na.rm = TRUE,
      names = FALSE
    )
  )
}

safe_median <- function(x) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  median(x)
}

safe_mean <- function(x) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  mean(x)
}

safe_max <- function(x) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  max(x)
}

safe_min <- function(x) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  min(x)
}

# ------------------------------------------------------------
# Input validation
# ------------------------------------------------------------

if (!file.exists(input_file)) {

  stop(
    "Input file not found: ",
    input_file
  )
}

cat("Reading:\n")
cat(input_file, "\n\n")

clumps <- fread(
  input_file,
  na.strings = c(
    "",
    "NA",
    "NaN"
  )
)

required_columns <- c(
  "clump_id",
  "chr",
  "lead_id",
  "lead_pos",
  "lead_p",
  "plink_total",
  "n_ids_parsed",
  "raw_start",
  "raw_end"
)

missing_columns <- setdiff(
  required_columns,
  names(clumps)
)

if (length(missing_columns) > 0L) {

  stop(
    "Missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

n_input_rows <- nrow(clumps)

cat(
  "Input clump rows:",
  n_input_rows,
  "\n"
)

# ------------------------------------------------------------
# Standardize columns
# ------------------------------------------------------------

clumps[, chr := gsub(
  "^chr",
  "",
  as.character(chr),
  ignore.case = TRUE
)]

numeric_columns <- c(
  "lead_pos",
  "lead_p",
  "plink_total",
  "n_ids_parsed",
  "raw_start",
  "raw_end"
)

for (column_name in numeric_columns) {

  clumps[
    ,
    (column_name) := as.numeric(
      get(column_name)
    )
  ]
}

# ------------------------------------------------------------
# Remove invalid records
# ------------------------------------------------------------

invalid_record <- (
  is.na(clumps$clump_id) |
    clumps$clump_id == "" |
    is.na(clumps$chr) |
    clumps$chr == "" |
    is.na(clumps$lead_id) |
    clumps$lead_id == "" |
    is.na(clumps$lead_pos) |
    is.na(clumps$lead_p) |
    is.na(clumps$raw_start) |
    is.na(clumps$raw_end) |
    clumps$lead_pos < 1 |
    clumps$raw_start < 1 |
    clumps$raw_end < clumps$raw_start
)

n_invalid_records <- sum(
  invalid_record
)

if (n_invalid_records > 0L) {

  warning(
    n_invalid_records,
    " invalid clump records will be removed."
  )

  clumps <- clumps[
    !invalid_record
  ]
}

# ------------------------------------------------------------
# Resolve duplicated clump IDs
# ------------------------------------------------------------

n_duplicate_clump_ids <- clumps[
  ,
  sum(
    duplicated(clump_id)
  )
]

if (n_duplicate_clump_ids > 0L) {

  warning(
    n_duplicate_clump_ids,
    " duplicated clump records detected. ",
    "The strongest row per clump will be retained."
  )
}

setorder(
  clumps,
  clump_id,
  lead_p
)

clumps <- clumps[
  ,
  .SD[1L],
  by = clump_id
]

# ------------------------------------------------------------
# Chromosome ordering
# ------------------------------------------------------------

clumps[
  ,
  chr_order := chromosome_order(chr)
]

if (anyNA(clumps$chr_order)) {

  warning(
    "Some chromosome names could not be assigned ",
    "a standard chromosome order."
  )
}

setorder(
  clumps,
  chr_order,
  lead_pos,
  lead_p,
  clump_id
)

# ------------------------------------------------------------
# Calculate lead SNP gaps
# ------------------------------------------------------------

clumps[
  ,
  previous_lead_pos := shift(
    lead_pos
  ),
  by = chr
]

clumps[
  ,
  lead_gap_bp := lead_pos -
    previous_lead_pos
]

clumps[
  ,
  starts_new_locus := fifelse(
    is.na(previous_lead_pos) |
      lead_gap_bp > lead_gap_threshold_bp,
    1L,
    0L
  )
]

clumps[
  ,
  chromosome_locus_group := cumsum(
    starts_new_locus
  ),
  by = chr
]

# ------------------------------------------------------------
# Create locus summaries
# ------------------------------------------------------------

loci <- clumps[
  ,
  {

    best_index <- which.min(
      lead_p
    )

    lead_minimum <- min(
      lead_pos,
      na.rm = TRUE
    )

    lead_maximum <- max(
      lead_pos,
      na.rm = TRUE
    )

    primary_start <- max(
      1,
      lead_minimum - locus_padding_bp
    )

    primary_end <- lead_maximum +
      locus_padding_bp

    support_start <- min(
      raw_start,
      na.rm = TRUE
    )

    support_end <- max(
      raw_end,
      na.rm = TRUE
    )

    list(
      start = as.integer(
        primary_start
      ),

      end = as.integer(
        primary_end
      ),

      locus_size_bp = as.integer(
        primary_end -
          primary_start +
          1
      ),

      minimum_lead_pos = as.integer(
        lead_minimum
      ),

      maximum_lead_pos = as.integer(
        lead_maximum
      ),

      lead_span_bp = as.integer(
        lead_maximum -
          lead_minimum
      ),

      raw_support_start = as.integer(
        support_start
      ),

      raw_support_end = as.integer(
        support_end
      ),

      raw_support_span_bp = as.integer(
        support_end -
          support_start +
          1
      ),

      lead_id = lead_id[
        best_index
      ],

      lead_pos = as.integer(
        lead_pos[
          best_index
        ]
      ),

      lead_p = lead_p[
        best_index
      ],

      n_clumps = uniqueN(
        clump_id
      ),

      sum_plink_total = sum(
        plink_total,
        na.rm = TRUE
      ),

      sum_ids_parsed = sum(
        n_ids_parsed,
        na.rm = TRUE
      ),

      minimum_internal_lead_gap_bp = {

        gaps <- lead_gap_bp[
          !is.na(lead_gap_bp)
        ]

        if (length(gaps) == 0L) {
          NA_integer_
        } else {
          as.integer(
            min(gaps)
          )
        }
      },

      maximum_internal_lead_gap_bp = {

        gaps <- lead_gap_bp[
          !is.na(lead_gap_bp)
        ]

        if (length(gaps) == 0L) {
          NA_integer_
        } else {
          as.integer(
            max(gaps)
          )
        }
      },

      member_clumps = collapse_unique(
        clump_id
      ),

      member_leads = collapse_unique(
        lead_id
      )
    )

  },
  by = .(
    chr,
    chromosome_locus_group
  )
]

loci[
  ,
  chr_order := chromosome_order(chr)
]

setorder(
  loci,
  chr_order,
  start,
  end,
  lead_p
)

loci[
  ,
  locus_number := seq_len(.N)
]

loci[
  ,
  locus_id := sprintf(
    "LOCUS_%04d_chr%s",
    locus_number,
    chr
  )
]

setcolorder(
  loci,
  c(
    "locus_id",
    "chr",
    "start",
    "end",
    "locus_size_bp",
    "minimum_lead_pos",
    "maximum_lead_pos",
    "lead_span_bp",
    "raw_support_start",
    "raw_support_end",
    "raw_support_span_bp",
    "lead_id",
    "lead_pos",
    "lead_p",
    "n_clumps",
    "sum_plink_total",
    "sum_ids_parsed",
    "minimum_internal_lead_gap_bp",
    "maximum_internal_lead_gap_bp",
    "member_clumps",
    "member_leads",
    "chromosome_locus_group",
    "locus_number",
    "chr_order"
  )
)

# ------------------------------------------------------------
# Assign locus IDs back to individual clumps
# ------------------------------------------------------------

membership_lookup <- loci[
  ,
  .(
    chr,
    chromosome_locus_group,
    locus_id,
    locus_start = start,
    locus_end = end,
    locus_size_bp,
    minimum_lead_pos,
    maximum_lead_pos,
    raw_support_start,
    raw_support_end
  )
]

membership <- merge(
  clumps,
  membership_lookup,
  by = c(
    "chr",
    "chromosome_locus_group"
  ),
  all.x = TRUE,
  sort = FALSE
)

setorder(
  membership,
  chr_order,
  locus_start,
  lead_pos,
  lead_p
)

membership[
  ,
  distance_from_locus_lead_min := lead_pos -
    minimum_lead_pos
]

membership[
  ,
  distance_from_locus_lead_max := maximum_lead_pos -
    lead_pos
]

# ------------------------------------------------------------
# Save clump-level v3 table
# ------------------------------------------------------------

clumps_to_save <- copy(
  clumps
)

clumps_to_save[
  ,
  chr_order := NULL
]

fwrite(
  clumps_to_save,
  output_clumps,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# Save final locus table
# ------------------------------------------------------------

loci_to_save <- copy(
  loci
)

loci_to_save[
  ,
  c(
    "chromosome_locus_group",
    "locus_number",
    "chr_order"
  ) := NULL
]

fwrite(
  loci_to_save,
  output_loci,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# Save membership table
# ------------------------------------------------------------

membership_to_save <- copy(
  membership
)

membership_to_save[
  ,
  c(
    "chr_order",
    "previous_lead_pos",
    "starts_new_locus"
  ) := NULL
]

setcolorder(
  membership_to_save,
  c(
    "locus_id",
    "locus_start",
    "locus_end",
    "locus_size_bp",
    "chr",
    "clump_id",
    "lead_id",
    "lead_pos",
    "lead_p",
    "previous_lead_pos",
    "lead_gap_bp",
    "plink_total",
    "n_ids_parsed",
    "raw_start",
    "raw_end",
    "minimum_lead_pos",
    "maximum_lead_pos",
    "raw_support_start",
    "raw_support_end",
    "distance_from_locus_lead_min",
    "distance_from_locus_lead_max",
    "chromosome_locus_group"
  ),
  skip_absent = TRUE
)

fwrite(
  membership_to_save,
  output_membership,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# QC statistics
# ------------------------------------------------------------

locus_sizes <- loci$locus_size_bp
lead_spans <- loci$lead_span_bp
support_spans <- loci$raw_support_span_bp
clumps_per_locus <- loci$n_clumps

all_lead_gaps <- clumps[
  !is.na(lead_gap_bp),
  lead_gap_bp
]

largest_index <- which.max(
  loci$locus_size_bp
)

qc <- data.table(
  parameter = c(
    "trait_id",
    "algorithm",
    "lead_gap_threshold_bp",
    "locus_padding_bp",
    "number_of_input_rows",
    "number_of_invalid_rows_removed",
    "number_of_valid_unique_clumps",
    "number_of_final_loci",
    "single_clump_loci",
    "multi_clump_loci",
    "minimum_locus_size_bp",
    "median_locus_size_bp",
    "mean_locus_size_bp",
    "locus_size_75th_percentile_bp",
    "locus_size_90th_percentile_bp",
    "locus_size_95th_percentile_bp",
    "locus_size_99th_percentile_bp",
    "maximum_locus_size_bp",
    "median_lead_span_bp",
    "maximum_lead_span_bp",
    "median_raw_support_span_bp",
    "maximum_raw_support_span_bp",
    "median_clumps_per_locus",
    "mean_clumps_per_locus",
    "maximum_clumps_per_locus",
    "median_consecutive_lead_gap_bp",
    "lead_gaps_over_100kb",
    "lead_gaps_over_250kb",
    "lead_gaps_over_500kb",
    "loci_larger_than_500kb",
    "loci_larger_than_1Mb",
    "loci_larger_than_2Mb",
    "largest_locus_id",
    "largest_locus_chr",
    "largest_locus_start",
    "largest_locus_end",
    "largest_locus_lead_id"
  ),

  value = c(
    trait_id,
    "consecutive_lead_gap_clustering",
    lead_gap_threshold_bp,
    locus_padding_bp,
    n_input_rows,
    n_invalid_records,
    nrow(clumps),
    nrow(loci),
    sum(
      loci$n_clumps == 1L
    ),
    sum(
      loci$n_clumps > 1L
    ),
    safe_min(
      locus_sizes
    ),
    safe_median(
      locus_sizes
    ),
    round(
      safe_mean(
        locus_sizes
      ),
      2
    ),
    safe_quantile(
      locus_sizes,
      0.75
    ),
    safe_quantile(
      locus_sizes,
      0.90
    ),
    safe_quantile(
      locus_sizes,
      0.95
    ),
    safe_quantile(
      locus_sizes,
      0.99
    ),
    safe_max(
      locus_sizes
    ),
    safe_median(
      lead_spans
    ),
    safe_max(
      lead_spans
    ),
    safe_median(
      support_spans
    ),
    safe_max(
      support_spans
    ),
    safe_median(
      clumps_per_locus
    ),
    round(
      safe_mean(
        clumps_per_locus
      ),
      2
    ),
    safe_max(
      clumps_per_locus
    ),
    safe_median(
      all_lead_gaps
    ),
    sum(
      all_lead_gaps > 100000,
      na.rm = TRUE
    ),
    sum(
      all_lead_gaps > 250000,
      na.rm = TRUE
    ),
    sum(
      all_lead_gaps > 500000,
      na.rm = TRUE
    ),
    sum(
      locus_sizes > 500000
    ),
    sum(
      locus_sizes > 1000000
    ),
    sum(
      locus_sizes > 2000000
    ),
    loci$locus_id[
      largest_index
    ],
    loci$chr[
      largest_index
    ],
    loci$start[
      largest_index
    ],
    loci$end[
      largest_index
    ],
    loci$lead_id[
      largest_index
    ]
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
# Largest loci
# ------------------------------------------------------------

largest_loci <- copy(
  loci_to_save
)

setorder(
  largest_loci,
  -locus_size_bp
)

largest_loci <- largest_loci[
  1:min(
    25L,
    .N
  )
]

fwrite(
  largest_loci,
  output_largest,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# TPCN2-region inspection
#
# Broad inspection window:
# chr11:67–70 Mb
# ------------------------------------------------------------

tpcn2_region <- loci_to_save[
  chr == "11" &
    end >= 67000000 &
    start <= 70000000
]

fwrite(
  tpcn2_region,
  output_tpcn2,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# Locus-size histogram
# ------------------------------------------------------------

png(
  filename = output_histogram,
  width = 1800,
  height = 1200,
  res = 180
)

hist(
  locus_sizes / 1000,
  breaks = 40,
  main = paste0(
    trait_id,
    ": lead-gap-defined locus sizes"
  ),
  xlab = "Primary locus size (kb)",
  ylab = "Number of loci"
)

abline(
  v = median(
    locus_sizes / 1000,
    na.rm = TRUE
  ),
  lty = 2,
  lwd = 2
)

legend(
  "topright",
  legend = paste0(
    "Median = ",
    round(
      median(
        locus_sizes / 1000,
        na.rm = TRUE
      ),
      1
    ),
    " kb"
  ),
  lty = 2,
  lwd = 2,
  bty = "n"
)

dev.off()

# ------------------------------------------------------------
# Lead-gap histogram
# ------------------------------------------------------------

png(
  filename = output_gap_histogram,
  width = 1800,
  height = 1200,
  res = 180
)

hist(
  all_lead_gaps / 1000,
  breaks = 60,
  main = paste0(
    trait_id,
    ": consecutive clump lead gaps"
  ),
  xlab = "Distance between consecutive lead variants (kb)",
  ylab = "Number of gaps"
)

abline(
  v = lead_gap_threshold_bp / 1000,
  lty = 2,
  lwd = 2
)

legend(
  "topright",
  legend = paste0(
    "Locus split threshold = ",
    lead_gap_threshold_bp / 1000,
    " kb"
  ),
  lty = 2,
  lwd = 2,
  bty = "n"
)

dev.off()

# ------------------------------------------------------------
# Console report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("PigmentRegAtlas Module 03 v3 completed\n")
cat("============================================================\n")

cat("\nConfiguration:\n")
cat(
  "Lead-gap threshold:",
  lead_gap_threshold_bp,
  "bp\n"
)

cat(
  "Primary locus padding:",
  locus_padding_bp,
  "bp per side\n"
)

cat("\nQC summary:\n")
print(qc)

cat("\nTPCN2-region loci:\n")

if (nrow(tpcn2_region) == 0L) {

  cat(
    "No loci intersect chr11:67,000,000-70,000,000.\n"
  )

} else {

  print(
    tpcn2_region[
      ,
      .(
        locus_id,
        chr,
        start,
        end,
        locus_size_bp,
        minimum_lead_pos,
        maximum_lead_pos,
        lead_span_bp,
        lead_id,
        lead_pos,
        lead_p,
        n_clumps,
        raw_support_start,
        raw_support_end,
        raw_support_span_bp
      )
    ]
  )
}

cat("\nLargest ten loci:\n")

print(
  largest_loci[
    1:min(
      10L,
      .N
    ),
    .(
      locus_id,
      chr,
      start,
      end,
      locus_size_bp,
      lead_id,
      lead_p,
      n_clumps
    )
  ]
)

cat("\nOutput files:\n")
cat(output_clumps, "\n")
cat(output_loci, "\n")
cat(output_membership, "\n")
cat(output_qc, "\n")
cat(output_largest, "\n")
cat(output_histogram, "\n")
cat(output_gap_histogram, "\n")
cat(output_tpcn2, "\n")