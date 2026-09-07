############################################################
# PigmentationAtlas
#
# Module: 07A
# Name: Build clump master metadata
#
# Purpose:
# Creates one master metadata table for all independent
# PLINK clumps. This table will control downstream analyses,
# including GWAS extraction, LD calculation, SuSiE
# fine-mapping, colocalization and figure generation.
#
# Inputs:
#   results/GCST90691600_clump_uniqueID_r2_0.1_kb1000.clumps
#   results/GCST90691600_genomic_loci_v3.tsv
#
# Output:
#   results/clump_master_metadata.tsv
#
# Project:
# PigmentationAtlas
############################################################

############################################################
## Start
############################################################

start_time <- Sys.time()

options(
    stringsAsFactors = FALSE,
    scipen = 999
)

suppressPackageStartupMessages(
    library(data.table)
)

cat("\n")
cat("============================================\n")
cat("PigmentationAtlas Module 07A\n")
cat("Build clump master metadata\n")
cat("============================================\n")

############################################################
## File paths
############################################################

clump_file <- file.path(
    "results",
    "GCST90691600_clump_uniqueID_r2_0.1_kb1000.clumps"
)

locus_file <- file.path(
    "results",
    "GCST90691600_genomic_loci_v3.tsv"
)

output_file <- file.path(
    "results",
    "clump_master_metadata.tsv"
)

############################################################
## Parameters
############################################################

padding_bp <- 50000L

############################################################
## Check input files
############################################################

required_files <- c(
    clump_file,
    locus_file
)

missing_files <- required_files[
    !file.exists(required_files)
]

if (length(missing_files) > 0L)
{
    stop(
        paste0(
            "Missing input file(s):\n",
            paste(
                missing_files,
                collapse = "\n"
            )
        )
    )
}

if (!dir.exists("results"))
{
    dir.create(
        "results",
        recursive = TRUE
    )
}

############################################################
## Read input tables
############################################################

cat("\nReading PLINK clumps...\n")

clumps <- fread(
    clump_file
)

cat(
    "PLINK clumps loaded:",
    nrow(clumps),
    "\n"
)

cat("\nReading genomic loci...\n")

loci <- fread(
    locus_file
)

cat(
    "Genomic loci loaded:",
    nrow(loci),
    "\n"
)

############################################################
## Validate input columns
############################################################

required_clump_columns <- c(
    "#CHROM",
    "POS",
    "ID",
    "P",
    "SP2"
)

missing_clump_columns <- setdiff(
    required_clump_columns,
    names(clumps)
)

if (length(missing_clump_columns) > 0L)
{
    stop(
        paste0(
            "Missing columns in clump file: ",
            paste(
                missing_clump_columns,
                collapse = ", "
            )
        )
    )
}

required_locus_columns <- c(
    "locus_id",
    "chr",
    "start",
    "end"
)

missing_locus_columns <- setdiff(
    required_locus_columns,
    names(loci)
)

if (length(missing_locus_columns) > 0L)
{
    stop(
        paste0(
            "Missing columns in locus file: ",
            paste(
                missing_locus_columns,
                collapse = ", "
            )
        )
    )
}

############################################################
## Standardize input column types
############################################################

clumps[, `#CHROM` := as.integer(`#CHROM`)]
clumps[, POS := as.integer(POS)]
clumps[, P := as.numeric(P)]
clumps[, ID := as.character(ID)]
clumps[, SP2 := as.character(SP2)]

loci[, chr := as.integer(chr)]
loci[, start := as.integer(start)]
loci[, end := as.integer(end)]
loci[, locus_id := as.character(locus_id)]

############################################################
## Helper function: parse PLINK SP2 field
############################################################

parse_sp2 <- function(index_id, sp2_string)
{
    index_id <- trimws(
        as.character(index_id)
    )

    variants <- index_id

    invalid_sp2_values <- c(
        "",
        ".",
        "NONE",
        "NA"
    )

    if (
        !is.na(sp2_string) &&
        !toupper(trimws(sp2_string)) %in%
            invalid_sp2_values
    )
    {
        sp2_variants <- trimws(
            unlist(
                strsplit(
                    as.character(sp2_string),
                    ",",
                    fixed = TRUE
                )
            )
        )

        valid_sp2 <- (
            !is.na(sp2_variants) &
            nzchar(sp2_variants) &
            !toupper(sp2_variants) %in%
                invalid_sp2_values &
            grepl(
                "^[0-9]+:[0-9]+:",
                sp2_variants
            )
        )

        sp2_variants <- sp2_variants[
            valid_sp2
        ]

        variants <- c(
            variants,
            sp2_variants
        )
    }

    unique(variants)
}

############################################################
## Helper function: extract chromosome and position
############################################################

extract_variant_positions <- function(variant_ids)
{
    variant_ids <- as.character(
        variant_ids
    )

    valid_format <- grepl(
        "^[0-9]+:[0-9]+:",
        variant_ids
    )

    if (any(!valid_format))
    {
        stop(
            paste0(
                "Invalid variant ID format: ",
                paste(
                    variant_ids[!valid_format],
                    collapse = ", "
                )
            )
        )
    }

    parts <- tstrsplit(
        variant_ids,
        ":",
        fixed = TRUE
    )

    result <- data.table(
        variant_id = variant_ids,
        chromosome = as.integer(parts[[1]]),
        position = as.integer(parts[[2]])
    )

    if (
        anyNA(result$chromosome) ||
        anyNA(result$position)
    )
    {
        stop(
            paste0(
                "Failed to parse chromosome or position for: ",
                paste(
                    variant_ids,
                    collapse = ", "
                )
            )
        )
    }

    result
}

############################################################
## Build master metadata
############################################################

cat("\nBuilding clump metadata...\n")

metadata_list <- vector(
    mode = "list",
    length = nrow(clumps)
)

for (i in seq_len(nrow(clumps)))
{
    clump_row <- clumps[i]

    index_variant <- as.character(
        clump_row$ID
    )

    expected_chr <- as.integer(
        clump_row[["#CHROM"]]
    )

    index_position <- as.integer(
        clump_row$POS
    )

    index_p <- as.numeric(
        clump_row$P
    )

    variant_ids <- parse_sp2(
        index_id = index_variant,
        sp2_string = clump_row$SP2
    )

    variant_positions <- extract_variant_positions(
        variant_ids
    )

    if (
        any(
            variant_positions$chromosome !=
                expected_chr
        )
    )
    {
        stop(
            paste0(
                "Chromosome mismatch in clump ",
                index_variant
            )
        )
    }

    raw_start <- min(
        variant_positions$position
    )

    raw_end <- max(
        variant_positions$position
    )

    region_start <- max(
        1L,
        raw_start - padding_bp
    )

    region_end <- raw_end + padding_bp

    parent <- loci[
        chr == expected_chr &
        start <= index_position &
        index_position <= end
    ]

    if (nrow(parent) != 1L)
    {
        stop(
            paste0(
                "Parent locus error for ",
                index_variant,
                ". Number of matching loci: ",
                nrow(parent)
            )
        )
    }

    metadata_list[[i]] <- data.table(
        clump_id = sprintf(
            "CLUMP_%04d",
            i
        ),
        parent_locus = parent$locus_id[1],
        chromosome = expected_chr,
        index_variant = index_variant,
        index_position = index_position,
        index_p = index_p,
        raw_start = raw_start,
        raw_end = raw_end,
        raw_span_bp = raw_end - raw_start + 1L,
        region_start = region_start,
        region_end = region_end,
        region_span_bp =
            region_end - region_start + 1L,
        n_variants = length(variant_ids),
        n_sp2_variants =
            length(variant_ids) - 1L,
        padding_bp = padding_bp,
        region_size_mb =
            round(
                (
                    region_end -
                    region_start +
                    1L
                ) / 1000000,
                6
            ),
        analysis_status = "pending",
        susie_status = "pending",
        coloc_status = "pending",
        figure_status = "pending"
    )

    if (
        i %% 100L == 0L ||
        i == nrow(clumps)
    )
    {
        cat(
            "Processed",
            i,
            "of",
            nrow(clumps),
            "clumps\n"
        )
    }
}

############################################################
## Combine rows
############################################################

clump_metadata <- rbindlist(
    metadata_list,
    use.names = TRUE,
    fill = FALSE
)

############################################################
## Sort by genomic position
############################################################

setorder(
    clump_metadata,
    chromosome,
    index_position
)

############################################################
## Reassign clump IDs after genomic sorting
############################################################

clump_metadata[
    ,
    clump_id := sprintf(
        "CLUMP_%04d",
        seq_len(.N)
    )
]

############################################################
## Quality control
############################################################

cat("\nRunning quality control...\n")

if (nrow(clump_metadata) != nrow(clumps))
{
    stop(
        "QC failure: output row count does not match clump count."
    )
}

if (anyDuplicated(clump_metadata$clump_id))
{
    stop(
        "QC failure: duplicated clump IDs."
    )
}

if (anyDuplicated(clump_metadata$index_variant))
{
    stop(
        "QC failure: duplicated index variants."
    )
}

if (anyNA(clump_metadata$parent_locus))
{
    stop(
        "QC failure: missing parent locus."
    )
}

if (anyNA(clump_metadata$raw_start))
{
    stop(
        "QC failure: missing raw_start."
    )
}

if (anyNA(clump_metadata$raw_end))
{
    stop(
        "QC failure: missing raw_end."
    )
}

if (anyNA(clump_metadata$region_start))
{
    stop(
        "QC failure: missing region_start."
    )
}

if (anyNA(clump_metadata$region_end))
{
    stop(
        "QC failure: missing region_end."
    )
}

if (anyNA(clump_metadata$region_span_bp))
{
    stop(
        "QC failure: missing region_span_bp."
    )
}

if (anyNA(clump_metadata$n_variants))
{
    stop(
        "QC failure: missing n_variants."
    )
}

if (
    any(
        clump_metadata$raw_start >
            clump_metadata$index_position
    ) ||
    any(
        clump_metadata$raw_end <
            clump_metadata$index_position
    )
)
{
    stop(
        "QC failure: index variant is outside raw clump interval."
    )
}

if (
    any(
        clump_metadata$region_start >
            clump_metadata$raw_start
    ) ||
    any(
        clump_metadata$region_end <
            clump_metadata$raw_end
    )
)
{
    stop(
        "QC failure: padded region does not contain raw clump interval."
    )
}

single_snp_clumps <- clump_metadata[
    n_variants == 1L
]

if (
    nrow(single_snp_clumps) > 0L &&
    any(
        single_snp_clumps$region_span_bp !=
            100001L
    )
)
{
    stop(
        "QC failure: invalid region size for single-SNP clumps."
    )
}

############################################################
## Save output
############################################################

fwrite(
    clump_metadata,
    output_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

############################################################
## Final report
############################################################

end_time <- Sys.time()

elapsed_seconds <- as.numeric(
    difftime(
        end_time,
        start_time,
        units = "secs"
    )
)

cat("\n")
cat("============================================\n")
cat("Module 07A completed successfully\n")
cat("============================================\n")
cat(
    "PLINK clumps          :",
    nrow(clump_metadata),
    "\n"
)
cat(
    "Genomic loci          :",
    uniqueN(clump_metadata$parent_locus),
    "\n"
)
cat(
    "Chromosomes           :",
    uniqueN(clump_metadata$chromosome),
    "\n"
)
cat(
    "Single-SNP clumps     :",
    sum(clump_metadata$n_variants == 1L),
    "\n"
)
cat(
    "Multi-variant clumps  :",
    sum(clump_metadata$n_variants > 1L),
    "\n"
)
cat(
    "Minimum region        :",
    min(clump_metadata$region_span_bp),
    "bp\n"
)
cat(
    "Median region         :",
    median(clump_metadata$region_span_bp),
    "bp\n"
)
cat(
    "Mean region           :",
    round(
        mean(clump_metadata$region_span_bp),
        1
    ),
    "bp\n"
)
cat(
    "Largest region        :",
    max(clump_metadata$region_span_bp),
    "bp\n"
)
cat(
    "Missing region spans  :",
    sum(
        is.na(
            clump_metadata$region_span_bp
        )
    ),
    "\n"
)
cat(
    "Output file           :",
    output_file,
    "\n"
)
cat(
    "Finished at           :",
    format(
        end_time,
        "%Y-%m-%d %H:%M:%S"
    ),
    "\n"
)
cat(
    "Elapsed time          :",
    round(
        elapsed_seconds,
        1
    ),
    "seconds\n"
)
cat("============================================\n")