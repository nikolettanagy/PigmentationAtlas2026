############################################################
# PigmentationAtlas
#
# Module: 08
# Name: Extract clump-specific GWAS regions
#
# Purpose:
# Extracts the harmonized GWAS variants falling inside each
# clump-specific genomic interval defined by Module 07A.
#
# Input:
#   results/clump_master_metadata.tsv
#   results/GCST90691600_harmonized.tsv.gz
#
# Output for each clump:
#   results/clumps/CLUMP_xxxx/gwas/gwas_region.tsv.gz
#
# Global output:
#   results/clump_gwas_extraction_manifest.tsv
#
# Notes:
# - The harmonized GWAS is read once.
# - Each clump interval is then extracted using data.table.
# - Overlapping clumps may contain overlapping GWAS variants.
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
cat("PigmentationAtlas Module 08\n")
cat("Extract clump-specific GWAS regions\n")
cat("============================================\n")

############################################################
## File paths
############################################################

metadata_file <- file.path(
    "results",
    "clump_master_metadata.tsv"
)

gwas_file <- file.path(
    "results",
    "GCST90691600_harmonized.tsv.gz"
)

workspace_root <- file.path(
    "results",
    "clumps"
)

manifest_file <- file.path(
    "results",
    "clump_gwas_extraction_manifest.tsv"
)

############################################################
## Parameters
############################################################

overwrite_existing <- TRUE

progress_interval <- 100L

############################################################
## Helper functions
############################################################

format_seconds <- function(seconds)
{
    if (seconds < 60)
    {
        return(
            paste0(
                round(seconds, 1),
                " seconds"
            )
        )
    }

    minutes <- floor(seconds / 60)
    remaining_seconds <- round(
        seconds - minutes * 60,
        1
    )

    paste0(
        minutes,
        " min ",
        remaining_seconds,
        " sec"
    )
}

find_column <- function(
    column_names,
    candidates,
    field_description
)
{
    match_index <- match(
        tolower(candidates),
        tolower(column_names)
    )

    match_index <- match_index[
        !is.na(match_index)
    ]

    if (length(match_index) == 0L)
    {
        stop(
            paste0(
                "Could not identify the ",
                field_description,
                " column. Tried: ",
                paste(
                    candidates,
                    collapse = ", "
                )
            )
        )
    }

    column_names[
        match_index[1]
    ]
}

normalise_chromosome <- function(x)
{
    x <- as.character(x)

    x <- sub(
        "^chr",
        "",
        x,
        ignore.case = TRUE
    )

    x
}
normalise_variant_id <- function(x)
{
    x <- as.character(x)

    x <- trimws(x)

    x <- sub(
        "^chr",
        "",
        x,
        ignore.case = TRUE
    )

    x <- sub(
        "_b3[78]$",
        "",
        x,
        ignore.case = TRUE
    )

    x <- gsub(
        "_",
        ":",
        x,
        fixed = TRUE
    )

    toupper(x)
}

############################################################
## Validate input files and directories
############################################################

if (!file.exists(metadata_file))
{
    stop(
        paste0(
            "Metadata file not found: ",
            metadata_file
        )
    )
}

if (!file.exists(gwas_file))
{
    stop(
        paste0(
            "Harmonized GWAS file not found: ",
            gwas_file
        )
    )
}

if (!dir.exists(workspace_root))
{
    stop(
        paste0(
            "Workspace root not found: ",
            workspace_root,
            "\nRun Module 07B first."
        )
    )
}

############################################################
## Read clump metadata
############################################################

cat("\nReading clump metadata...\n")

clump_metadata <- fread(
    metadata_file
)

cat(
    "Clumps loaded:",
    nrow(clump_metadata),
    "\n"
)

required_metadata_columns <- c(
    "clump_id",
    "chromosome",
    "index_variant",
    "index_position",
    "region_start",
    "region_end"
)

missing_metadata_columns <- setdiff(
    required_metadata_columns,
    names(clump_metadata)
)

if (length(missing_metadata_columns) > 0L)
{
    stop(
        paste0(
            "Missing metadata columns: ",
            paste(
                missing_metadata_columns,
                collapse = ", "
            )
        )
    )
}

if (nrow(clump_metadata) == 0L)
{
    stop(
        "The clump metadata table is empty."
    )
}

if (anyDuplicated(clump_metadata$clump_id))
{
    stop(
        "Duplicated clump IDs detected."
    )
}

if (
    anyNA(clump_metadata$chromosome) ||
    anyNA(clump_metadata$region_start) ||
    anyNA(clump_metadata$region_end)
)
{
    stop(
        "Missing chromosome or region coordinates detected."
    )
}

if (
    any(
        clump_metadata$region_start >
            clump_metadata$region_end
    )
)
{
    stop(
        "Invalid clump intervals detected."
    )
}

clump_metadata[
    ,
    chromosome_normalised :=
        normalise_chromosome(chromosome)
]

############################################################
## Check workspaces
############################################################

workspace_paths <- file.path(
    workspace_root,
    clump_metadata$clump_id
)

missing_workspaces <- clump_metadata$clump_id[
    !dir.exists(workspace_paths)
]

if (length(missing_workspaces) > 0L)
{
    stop(
        paste0(
            "Missing clump workspaces: ",
            paste(
                missing_workspaces,
                collapse = ", "
            )
        )
    )
}

gwas_directories <- file.path(
    workspace_paths,
    "gwas"
)

missing_gwas_directories <- clump_metadata$clump_id[
    !dir.exists(gwas_directories)
]

if (length(missing_gwas_directories) > 0L)
{
    stop(
        paste0(
            "Missing gwas directories: ",
            paste(
                missing_gwas_directories,
                collapse = ", "
            )
        )
    )
}

############################################################
## Read harmonized GWAS
############################################################

cat("\nReading harmonized GWAS...\n")
cat("Input:", gwas_file, "\n")

gwas_read_start <- Sys.time()

gwas <- fread(
    gwas_file,
    showProgress = TRUE
)

gwas_read_end <- Sys.time()

cat(
    "GWAS variants loaded:",
    nrow(gwas),
    "\n"
)

cat(
    "GWAS columns:",
    ncol(gwas),
    "\n"
)

cat(
    "Read time:",
    format_seconds(
        as.numeric(
            difftime(
                gwas_read_end,
                gwas_read_start,
                units = "secs"
            )
        )
    ),
    "\n"
)

if (nrow(gwas) == 0L)
{
    stop(
        "The harmonized GWAS file contains no variants."
    )
}

############################################################
## Identify essential GWAS columns
############################################################

chromosome_column <- find_column(
    names(gwas),
    c(
        "chromosome",
        "chr",
        "chrom"
    ),
    "chromosome"
)

position_column <- find_column(
    names(gwas),
    c(
        "base_pair_location",
        "position",
        "pos",
        "bp"
    ),
    "base-pair position"
)

variant_id_column <- find_column(
    names(gwas),
    c(
        "variant_id",
        "variant",
        "id"
    ),
    "variant identifier"
)

p_value_column <- find_column(
    names(gwas),
    c(
        "p_value",
        "pvalue",
        "p",
        "pval"
    ),
    "P-value"
)

rsid_candidates <- c(
    "rsid",
    "rs_id",
    "snp"
)

rsid_match <- match(
    tolower(rsid_candidates),
    tolower(names(gwas))
)

rsid_match <- rsid_match[
    !is.na(rsid_match)
]

if (length(rsid_match) > 0L)
{
    rsid_column <- names(gwas)[
        rsid_match[1]
    ]
} else {
    rsid_column <- NA_character_
}

cat("\nDetected GWAS columns:\n")
cat(
    "Chromosome :",
    chromosome_column,
    "\n"
)
cat(
    "Position   :",
    position_column,
    "\n"
)
cat(
    "Variant ID :",
    variant_id_column,
    "\n"
)
cat(
    "P-value    :",
    p_value_column,
    "\n"
)
cat(
    "rsID       :",
    ifelse(
        is.na(rsid_column),
        "not detected",
        rsid_column
    ),
    "\n"
)

############################################################
## Normalise GWAS chromosome and coordinate columns
############################################################

gwas[
    ,
    chromosome_module08 :=
        normalise_chromosome(
            get(chromosome_column)
        )
]

gwas[
    ,
    position_module08 :=
        suppressWarnings(
            as.integer(
                get(position_column)
            )
        )
]

if (anyNA(gwas$position_module08))
{
    n_missing_positions <- sum(
        is.na(gwas$position_module08)
    )

    stop(
        paste0(
            "The GWAS contains ",
            n_missing_positions,
            " missing or invalid genomic positions."
        )
    )
}

############################################################
## GWAS-level QC
############################################################

available_chromosomes <- unique(
    gwas$chromosome_module08
)

required_chromosomes <- unique(
    clump_metadata$chromosome_normalised
)

missing_chromosomes <- setdiff(
    required_chromosomes,
    available_chromosomes
)

if (length(missing_chromosomes) > 0L)
{
    stop(
        paste0(
            "The following clump chromosomes are absent from the GWAS: ",
            paste(
                missing_chromosomes,
                collapse = ", "
            )
        )
    )
}

invalid_p_values <- (
    is.na(gwas[[p_value_column]]) |
    gwas[[p_value_column]] < 0 |
    gwas[[p_value_column]] > 1
)

n_invalid_p_values <- sum(
    invalid_p_values
)

if (n_invalid_p_values > 0L)
{
    stop(
        paste0(
            "Invalid GWAS P-values detected: ",
            n_invalid_p_values
        )
    )
}

############################################################
## Index GWAS table
############################################################

cat("\nIndexing GWAS table by chromosome and position...\n")

setkeyv(
    gwas,
    c(
        "chromosome_module08",
        "position_module08"
    )
)

############################################################
## Prepare extraction manifest
############################################################

manifest_list <- vector(
    mode = "list",
    length = nrow(clump_metadata)
)

############################################################
## Extract clump regions
############################################################

cat("\nExtracting GWAS regions...\n")

extraction_start <- Sys.time()

for (i in seq_len(nrow(clump_metadata)))
{
    metadata_row <- clump_metadata[i]

    clump_id <- metadata_row$clump_id

    chromosome_value <-
        metadata_row$chromosome_normalised

    region_start <-
        as.integer(metadata_row$region_start)

    region_end <-
        as.integer(metadata_row$region_end)

    index_position <-
        as.integer(metadata_row$index_position)

    output_file <- file.path(
        workspace_root,
        clump_id,
        "gwas",
        "gwas_region.tsv.gz"
    )

    output_preexisting <- file.exists(
        output_file
    )

    extraction_status <- "pending"
    output_written <- FALSE

    region_gwas <- gwas[
        chromosome_module08 ==
            chromosome_value &
        position_module08 >=
            region_start &
        position_module08 <=
            region_end
    ]

    n_region_variants <- nrow(
        region_gwas
    )

    index_position_present <- any(
        region_gwas$position_module08 ==
            index_position
    )

   normalised_index_variant <- normalise_variant_id(
    metadata_row$index_variant
    )

   normalised_region_variants <- normalise_variant_id(
    region_gwas[[variant_id_column]]
    )

   index_variant_present <- any(
    normalised_region_variants ==
        normalised_index_variant
    )

    minimum_position <- if (
        n_region_variants > 0L
    )
    {
        min(
            region_gwas$position_module08
        )
    } else {
        NA_integer_
    }

    maximum_position <- if (
        n_region_variants > 0L
    )
    {
        max(
            region_gwas$position_module08
        )
    } else {
        NA_integer_
    }

    minimum_p_value <- if (
        n_region_variants > 0L
    )
    {
        min(
            region_gwas[[p_value_column]],
            na.rm = TRUE
        )
    } else {
        NA_real_
    }

    if (n_region_variants == 0L)
    {
        extraction_status <- "empty_region"
    } else if (!index_position_present)
    {
        extraction_status <-
            "index_position_missing"
    } else if (!index_variant_present)
    {
        extraction_status <-
            "index_variant_id_missing"
    } else {
        extraction_status <- "ready"
    }

    if (
        n_region_variants > 0L &&
        (
            overwrite_existing ||
            !output_preexisting
        )
    )
    {
        output_columns <- setdiff(
            names(region_gwas),
            c(
                "chromosome_module08",
                "position_module08"
            )
        )

        fwrite(
            region_gwas[
                ,
                ..output_columns
            ],
            output_file,
            sep = "\t",
            quote = FALSE,
            na = "NA",
            compress = "gzip"
        )

        output_written <- TRUE
    }

    output_exists_after <- file.exists(
        output_file
    )

    output_size_bytes <- if (
        output_exists_after
    )
    {
        file.info(
            output_file
        )$size
    } else {
        NA_real_
    }

    manifest_list[[i]] <- data.table(
        clump_id = clump_id,
        parent_locus =
            metadata_row$parent_locus,
        chromosome =
            metadata_row$chromosome,
        index_variant =
            metadata_row$index_variant,
        index_position =
            index_position,
        region_start =
            region_start,
        region_end =
            region_end,
        region_span_bp =
            metadata_row$region_span_bp,
        n_gwas_variants =
            n_region_variants,
        minimum_gwas_position =
            minimum_position,
        maximum_gwas_position =
            maximum_position,
        minimum_p_value =
            minimum_p_value,
        index_position_present =
            index_position_present,
        index_variant_present =
            index_variant_present,
        output_preexisting =
            output_preexisting,
        output_written =
            output_written,
        output_exists =
            output_exists_after,
        output_size_bytes =
            output_size_bytes,
        extraction_status =
            extraction_status,
        output_file =
            normalizePath(
                output_file,
                winslash = "/",
                mustWork = FALSE
            )
    )

    if (
        i %% progress_interval == 0L ||
        i == nrow(clump_metadata)
    )
    {
        elapsed_so_far <- as.numeric(
            difftime(
                Sys.time(),
                extraction_start,
                units = "secs"
            )
        )

        cat(
            "Processed",
            i,
            "of",
            nrow(clump_metadata),
            "clumps | elapsed:",
            format_seconds(
                elapsed_so_far
            ),
            "\n"
        )
    }
}

############################################################
## Combine extraction manifest
############################################################

extraction_manifest <- rbindlist(
    manifest_list,
    use.names = TRUE,
    fill = FALSE
)

setorder(
    extraction_manifest,
    chromosome,
    index_position
)

############################################################
## Global extraction QC
############################################################

cat("\nRunning extraction quality control...\n")

expected_clumps <- nrow(
    clump_metadata
)

n_output_files <- sum(
    extraction_manifest$output_exists
)

n_ready <- sum(
    extraction_manifest$extraction_status ==
        "ready"
)

n_empty_regions <- sum(
    extraction_manifest$extraction_status ==
        "empty_region"
)

n_missing_index_positions <- sum(
    extraction_manifest$extraction_status ==
        "index_position_missing"
)

n_missing_index_variant_ids <- sum(
    extraction_manifest$extraction_status ==
        "index_variant_id_missing"
)

n_zero_variant_outputs <- sum(
    extraction_manifest$n_gwas_variants == 0L
)

if (n_empty_regions > 0L)
{
    warning(
        paste0(
            n_empty_regions,
            " clump regions contained no GWAS variants."
        )
    )
}

if (n_missing_index_positions > 0L)
{
    warning(
        paste0(
            n_missing_index_positions,
            " clumps did not contain their index position."
        )
    )
}

if (n_missing_index_variant_ids > 0L)
{
    warning(
        paste0(
            n_missing_index_variant_ids,
            " clumps contained the index position but not the exact index variant ID."
        )
    )
}

############################################################
## Validate written region files
############################################################

cat("Validating written region files...\n")

files_to_validate <- extraction_manifest[
    output_exists == TRUE
]

file_validation <- logical(
    nrow(files_to_validate)
)

if (nrow(files_to_validate) > 0L)
{
    for (i in seq_len(nrow(files_to_validate)))
    {
        validation_header <- names(
            fread(
                files_to_validate$output_file[i],
                nrows = 0L
            )
        )

        file_validation[i] <- all(
            c(
                chromosome_column,
                position_column,
                variant_id_column,
                p_value_column
            ) %in%
                validation_header
        )
    }
}

n_invalid_output_files <- sum(
    !file_validation
)

if (n_invalid_output_files > 0L)
{
    stop(
        paste0(
            "QC failure: ",
            n_invalid_output_files,
            " output files have invalid headers."
        )
    )
}

############################################################
## Save manifest
############################################################

fwrite(
    extraction_manifest,
    manifest_file,
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

region_variant_summary <- summary(
    extraction_manifest$n_gwas_variants
)

cat("\n")
cat("============================================\n")
cat("Module 08 completed\n")
cat("============================================\n")
cat(
    "Expected clumps              :",
    expected_clumps,
    "\n"
)
cat(
    "GWAS variants loaded         :",
    nrow(gwas),
    "\n"
)
cat(
    "Output files                 :",
    n_output_files,
    "\n"
)
cat(
    "Ready clumps                 :",
    n_ready,
    "\n"
)
cat(
    "Empty regions                :",
    n_empty_regions,
    "\n"
)
cat(
    "Missing index positions      :",
    n_missing_index_positions,
    "\n"
)
cat(
    "Missing exact index IDs      :",
    n_missing_index_variant_ids,
    "\n"
)
cat(
    "Invalid output files         :",
    n_invalid_output_files,
    "\n"
)
cat(
    "Minimum variants per clump   :",
    min(
        extraction_manifest$n_gwas_variants
    ),
    "\n"
)
cat(
    "Median variants per clump    :",
    median(
        extraction_manifest$n_gwas_variants
    ),
    "\n"
)
cat(
    "Mean variants per clump      :",
    round(
        mean(
            extraction_manifest$n_gwas_variants
        ),
        1
    ),
    "\n"
)
cat(
    "Maximum variants per clump   :",
    max(
        extraction_manifest$n_gwas_variants
    ),
    "\n"
)
cat(
    "Manifest file                :",
    manifest_file,
    "\n"
)
cat(
    "Finished at                  :",
    format(
        end_time,
        "%Y-%m-%d %H:%M:%S"
    ),
    "\n"
)
cat(
    "Elapsed time                 :",
    format_seconds(
        elapsed_seconds
    ),
    "\n"
)
cat("============================================\n")

############################################################
## Strict failure criteria
############################################################

if (n_output_files != expected_clumps)
{
    stop(
        paste0(
            "Module 08 QC failure: expected ",
            expected_clumps,
            " output files, but found ",
            n_output_files,
            ". Review ",
            manifest_file,
            "."
        )
    )
}

if (n_zero_variant_outputs > 0L)
{
    stop(
        paste0(
            "Module 08 QC failure: ",
            n_zero_variant_outputs,
            " regions contained zero GWAS variants. Review ",
            manifest_file,
            "."
        )
    )
}

cat("\nModule 08 completed successfully.\n")