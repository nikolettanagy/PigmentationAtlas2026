############################################################
# PigmentationAtlas
#
# Module: 08B
# Name: Validate clump-specific GWAS regions
#
# Purpose:
# Validates the GWAS region files created by Module 08
# without re-reading the full harmonized GWAS file.
#
# Input:
#   results/clump_master_metadata.tsv
#   results/clump_gwas_extraction_manifest.tsv
#   results/clumps/CLUMP_xxxx/gwas/gwas_region.tsv.gz
#
# Output:
#   results/clump_gwas_validation_manifest.tsv
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
cat("PigmentationAtlas Module 08B\n")
cat("Validate clump-specific GWAS regions\n")
cat("============================================\n")

############################################################
## File paths
############################################################

metadata_file <- file.path(
    "results",
    "clump_master_metadata.tsv"
)

extraction_manifest_file <- file.path(
    "results",
    "clump_gwas_extraction_manifest.tsv"
)

workspace_root <- file.path(
    "results",
    "clumps"
)

validation_manifest_file <- file.path(
    "results",
    "clump_gwas_validation_manifest.tsv"
)

############################################################
## Parameters
############################################################

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

    x <- trimws(x)

    x <- sub(
        "^chr",
        "",
        x,
        ignore.case = TRUE
    )

    toupper(x)
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
        "([_:])b3[78]$",
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
## Validate required inputs
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

if (!file.exists(extraction_manifest_file))
{
    stop(
        paste0(
            "Extraction manifest not found: ",
            extraction_manifest_file,
            "\nRun Module 08 first."
        )
    )
}

if (!dir.exists(workspace_root))
{
    stop(
        paste0(
            "Workspace root not found: ",
            workspace_root
        )
    )
}

############################################################
## Read metadata and Module 08 manifest
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
        "Duplicated clump IDs detected in metadata."
    )
}

cat("Reading Module 08 extraction manifest...\n")

extraction_manifest <- fread(
    extraction_manifest_file
)

if (
    !"clump_id" %in%
        names(extraction_manifest)
)
{
    stop(
        "The Module 08 extraction manifest has no clump_id column."
    )
}

if (anyDuplicated(extraction_manifest$clump_id))
{
    stop(
        "Duplicated clump IDs detected in the Module 08 extraction manifest."
    )
}

missing_from_extraction_manifest <- setdiff(
    clump_metadata$clump_id,
    extraction_manifest$clump_id
)

if (length(missing_from_extraction_manifest) > 0L)
{
    stop(
        paste0(
            "Clumps missing from the Module 08 extraction manifest: ",
            paste(
                missing_from_extraction_manifest,
                collapse = ", "
            )
        )
    )
}

############################################################
## Prepare validation
############################################################

clump_metadata[
    ,
    chromosome_normalised :=
        normalise_chromosome(chromosome)
]

validation_list <- vector(
    mode = "list",
    length = nrow(clump_metadata)
)

cat("\nValidating existing GWAS region files...\n")

validation_start <- Sys.time()

############################################################
## Validate each clump region file
############################################################

for (i in seq_len(nrow(clump_metadata)))
{
    metadata_row <- clump_metadata[i]

    clump_id <- metadata_row$clump_id

    region_file <- file.path(
        workspace_root,
        clump_id,
        "gwas",
        "gwas_region.tsv.gz"
    )

    file_exists <- file.exists(
        region_file
    )

    file_size_bytes <- if (file_exists)
    {
        file.info(
            region_file
        )$size
    } else {
        NA_real_
    }

    file_readable <- FALSE
    header_valid <- FALSE
    region_nonempty <- FALSE
    chromosome_consistent <- FALSE
    positions_within_interval <- FALSE
    index_position_present <- FALSE
    index_variant_present <- FALSE
    n_gwas_variants <- 0L
    minimum_gwas_position <- NA_integer_
    maximum_gwas_position <- NA_integer_
    detected_chromosome_column <- NA_character_
    detected_position_column <- NA_character_
    detected_variant_id_column <- NA_character_
    validation_status <- "pending"
    validation_message <- NA_character_

    if (!file_exists)
    {
        validation_status <- "missing_file"
        validation_message <- "GWAS region file does not exist."
    } else if (
        is.na(file_size_bytes) ||
        file_size_bytes <= 0
    )
    {
        validation_status <- "empty_file"
        validation_message <- "GWAS region file has zero size."
    } else {
        region_gwas <- tryCatch(
            fread(
                region_file,
                showProgress = FALSE
            ),
            error = function(e)
            {
                validation_message <<-
                    conditionMessage(e)

                NULL
            }
        )

        if (is.null(region_gwas))
        {
            validation_status <- "unreadable_file"
        } else {
            file_readable <- TRUE

            detected_chromosome_column <- tryCatch(
                find_column(
                    names(region_gwas),
                    c(
                        "chromosome",
                        "chr",
                        "chrom"
                    ),
                    "chromosome"
                ),
                error = function(e)
                {
                    NA_character_
                }
            )

            detected_position_column <- tryCatch(
                find_column(
                    names(region_gwas),
                    c(
                        "base_pair_location",
                        "position",
                        "pos",
                        "bp"
                    ),
                    "base-pair position"
                ),
                error = function(e)
                {
                    NA_character_
                }
            )

            detected_variant_id_column <- tryCatch(
                find_column(
                    names(region_gwas),
                    c(
                        "variant_id",
                        "variant",
                        "id"
                    ),
                    "variant identifier"
                ),
                error = function(e)
                {
                    NA_character_
                }
            )

            header_valid <- all(
                !is.na(
                    c(
                        detected_chromosome_column,
                        detected_position_column,
                        detected_variant_id_column
                    )
                )
            )

            if (!header_valid)
            {
                validation_status <- "invalid_header"
                validation_message <-
                    "Required chromosome, position, or variant ID column is missing."
            } else {
                n_gwas_variants <- nrow(
                    region_gwas
                )

                region_nonempty <-
                    n_gwas_variants > 0L

                if (!region_nonempty)
                {
                    validation_status <- "empty_region"
                    validation_message <-
                        "GWAS region file contains no variants."
                } else {
                    region_chromosome <- normalise_chromosome(
                        region_gwas[[detected_chromosome_column]]
                    )

                    region_position <- suppressWarnings(
                        as.integer(
                            region_gwas[[detected_position_column]]
                        )
                    )

                    chromosome_consistent <- all(
                        !is.na(region_chromosome) &
                        region_chromosome ==
                            metadata_row$chromosome_normalised
                    )

                    positions_within_interval <- all(
                        !is.na(region_position) &
                        region_position >=
                            as.integer(metadata_row$region_start) &
                        region_position <=
                            as.integer(metadata_row$region_end)
                    )

                    minimum_gwas_position <- min(
                        region_position,
                        na.rm = TRUE
                    )

                    maximum_gwas_position <- max(
                        region_position,
                        na.rm = TRUE
                    )

                    index_position_present <- any(
                        region_position ==
                            as.integer(metadata_row$index_position),
                        na.rm = TRUE
                    )

                    normalised_region_variants <-
                        normalise_variant_id(
                            region_gwas[[detected_variant_id_column]]
                        )

                    normalised_index_variant <-
                        normalise_variant_id(
                            metadata_row$index_variant
                        )

                    index_variant_present <- any(
                        normalised_region_variants ==
                            normalised_index_variant,
                        na.rm = TRUE
                    )

                    if (!chromosome_consistent)
                    {
                        validation_status <-
                            "chromosome_mismatch"

                        validation_message <-
                            "One or more variants have a chromosome inconsistent with metadata."
                    } else if (!positions_within_interval)
                    {
                        validation_status <-
                            "position_outside_interval"

                        validation_message <-
                            "One or more variants fall outside the metadata interval."
                    } else if (!index_position_present)
                    {
                        validation_status <-
                            "index_position_missing"

                        validation_message <-
                            "The clump index position is absent from the region file."
                    } else if (!index_variant_present)
                    {
                        validation_status <-
                            "index_variant_missing"

                        validation_message <-
                            "The normalized clump index variant ID is absent from the region file."
                    } else {
                        validation_status <- "ready"
                        validation_message <- "PASS"
                    }
                }
            }
        }
    }

    validation_list[[i]] <- data.table(
        clump_id = clump_id,
        parent_locus =
            metadata_row$parent_locus,
        chromosome =
            metadata_row$chromosome,
        index_variant =
            metadata_row$index_variant,
        index_position =
            as.integer(metadata_row$index_position),
        region_start =
            as.integer(metadata_row$region_start),
        region_end =
            as.integer(metadata_row$region_end),
        file_exists =
            file_exists,
        file_size_bytes =
            file_size_bytes,
        file_readable =
            file_readable,
        header_valid =
            header_valid,
        region_nonempty =
            region_nonempty,
        chromosome_consistent =
            chromosome_consistent,
        positions_within_interval =
            positions_within_interval,
        n_gwas_variants =
            n_gwas_variants,
        minimum_gwas_position =
            minimum_gwas_position,
        maximum_gwas_position =
            maximum_gwas_position,
        index_position_present =
            index_position_present,
        index_variant_present =
            index_variant_present,
        validation_status =
            validation_status,
        validation_message =
            validation_message,
        region_file =
            normalizePath(
                region_file,
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
                validation_start,
                units = "secs"
            )
        )

        cat(
            "Validated",
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
## Combine validation manifest
############################################################

validation_manifest <- rbindlist(
    validation_list,
    use.names = TRUE,
    fill = FALSE
)

setorder(
    validation_manifest,
    chromosome,
    index_position
)

############################################################
## Save validation manifest
############################################################

fwrite(
    validation_manifest,
    validation_manifest_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

############################################################
## Final QC summary
############################################################

expected_clumps <- nrow(
    clump_metadata
)

n_files_present <- sum(
    validation_manifest$file_exists
)

n_ready <- sum(
    validation_manifest$validation_status ==
        "ready"
)

n_missing_files <- sum(
    validation_manifest$validation_status ==
        "missing_file"
)

n_empty_files <- sum(
    validation_manifest$validation_status ==
        "empty_file"
)

n_unreadable_files <- sum(
    validation_manifest$validation_status ==
        "unreadable_file"
)

n_invalid_headers <- sum(
    validation_manifest$validation_status ==
        "invalid_header"
)

n_empty_regions <- sum(
    validation_manifest$validation_status ==
        "empty_region"
)

n_chromosome_mismatches <- sum(
    validation_manifest$validation_status ==
        "chromosome_mismatch"
)

n_position_interval_failures <- sum(
    validation_manifest$validation_status ==
        "position_outside_interval"
)

n_missing_index_positions <- sum(
    validation_manifest$validation_status ==
        "index_position_missing"
)

n_missing_index_variants <- sum(
    validation_manifest$validation_status ==
        "index_variant_missing"
)

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
cat("Module 08B completed\n")
cat("============================================\n")
cat(
    "Expected clumps              :",
    expected_clumps,
    "\n"
)
cat(
    "Files present                :",
    n_files_present,
    "\n"
)
cat(
    "Ready clumps                 :",
    n_ready,
    "\n"
)
cat(
    "Missing files                :",
    n_missing_files,
    "\n"
)
cat(
    "Empty files                  :",
    n_empty_files,
    "\n"
)
cat(
    "Unreadable files             :",
    n_unreadable_files,
    "\n"
)
cat(
    "Invalid headers              :",
    n_invalid_headers,
    "\n"
)
cat(
    "Empty regions                :",
    n_empty_regions,
    "\n"
)
cat(
    "Chromosome mismatches        :",
    n_chromosome_mismatches,
    "\n"
)
cat(
    "Positions outside interval   :",
    n_position_interval_failures,
    "\n"
)
cat(
    "Missing index positions      :",
    n_missing_index_positions,
    "\n"
)
cat(
    "Missing index variants       :",
    n_missing_index_variants,
    "\n"
)
cat(
    "Validation manifest          :",
    validation_manifest_file,
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

if (n_ready != expected_clumps)
{
    failed_statuses <- validation_manifest[
        validation_status != "ready",
        .N,
        by = validation_status
    ]

    stop(
        paste0(
            "Module 08B QC failure: ",
            expected_clumps - n_ready,
            " clumps failed validation. Review ",
            validation_manifest_file,
            ".\nFailure summary:\n",
            paste(
                paste0(
                    failed_statuses$validation_status,
                    ": ",
                    failed_statuses$N
                ),
                collapse = "\n"
            )
        )
    )
}

cat("\nModule 08B completed successfully.\n")
