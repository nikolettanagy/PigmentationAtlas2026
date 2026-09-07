############################################################
# PigmentationAtlas
#
# Module: 09
# Name: Calculate clump-specific LD matrices
#
# Purpose:
# Uses the 1000 Genomes GRCh38 EUR PLINK2 reference panel
# to calculate signed, unphased LD correlation matrices for
# GWAS variants present in each validated clump region.
#
# Input:
#   results/clump_gwas_validation_manifest.tsv
#   results/clumps/CLUMP_xxxx/gwas/gwas_region.tsv.gz
#   reference/1000G_GRCh38/1000G_GRCh38_EUR.{pgen,pvar,psam}
#
# Output for each clump:
#   results/clumps/CLUMP_xxxx/ld/gwas_rsids.txt
#   results/clumps/CLUMP_xxxx/ld/ld_matrix.unphased.vcor1.bin
#   results/clumps/CLUMP_xxxx/ld/ld_matrix.unphased.vcor1.vars
#   results/clumps/CLUMP_xxxx/ld/ld_variant_map.tsv
#   results/clumps/CLUMP_xxxx/logs/module09_plink.log
#
# Global output:
#   results/clump_ld_manifest.tsv
#
# Notes:
# - The initial configuration is PILOT mode (first 5 clumps).
# - Change run_mode to "full" only after the pilot passes.
# - Existing valid LD outputs are skipped when overwrite_existing
#   is FALSE.
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
cat("PigmentationAtlas Module 09\n")
cat("Calculate clump-specific LD matrices\n")
cat("============================================\n")

############################################################
## File paths
############################################################

validation_manifest_file <- file.path(
    "results",
    "clump_gwas_validation_manifest.tsv"
)

workspace_root <- file.path(
    "results",
    "clumps"
)

reference_prefix <- file.path(
    "reference",
    "1000G_GRCh38",
    "1000G_GRCh38_EUR_uniqueID"
)

plink2_executable <- file.path(
    "tools",
    "plink2",
    "plink2.exe"
)

ld_manifest_file <- file.path(
    "results",
    "clump_ld_manifest.tsv"
)

############################################################
## Parameters
############################################################

# Use "pilot" for the first test and "full" after the pilot passes.
run_mode <- "full"

pilot_n_clumps <- 5L

overwrite_existing <- FALSE

progress_interval <- 1L

# Reference-panel filtering.
minimum_reference_maf <- 0.01

# Dense matrix safety limit.
# A bin4 square matrix requires approximately 4 * n^2 bytes.
maximum_ld_variants <- 10000L

# PLINK resource settings.
plink_threads <- 4L
plink_memory_mb <- 8000L

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

normalise_chromosome <- function(x)
{
    x <- as.character(x)
    x <- trimws(x)
    x <- sub("^chr", "", x, ignore.case = TRUE)
    toupper(x)
}

normalise_variant_id <- function(x)
{
    x <- as.character(x)
    x <- trimws(x)
    x <- sub("^chr", "", x, ignore.case = TRUE)
    x <- sub("([_:])b3[78]$", "", x, ignore.case = TRUE)
    x <- gsub("_", ":", x, fixed = TRUE)
    toupper(x)
}

is_missing_rsid <- function(x)
{
    x <- trimws(as.character(x))

    is.na(x) |
        x == "" |
        x == "." |
        toupper(x) == "NA" |
        !grepl("^rs[0-9]+$", x, ignore.case = TRUE)
}

find_column <- function(
    column_names,
    candidates,
    field_description,
    required = TRUE
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
        if (required)
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

        return(NA_character_)
    }

    column_names[
        match_index[1]
    ]
}

read_variant_ids <- function(filename)
{
    if (!file.exists(filename))
    {
        return(character())
    }

    lines <- readLines(
        filename,
        warn = FALSE
    )

    trimws(lines[
        nzchar(trimws(lines))
    ])
}

locate_ld_outputs <- function(output_prefix)
{
    output_directory <- dirname(output_prefix)
    output_basename <- basename(output_prefix)

    matrix_candidates <- list.files(
        output_directory,
        pattern = paste0(
            "^",
            output_basename,
            "\\.unphased\\.vcor1\\.bin$"
        ),
        full.names = TRUE
    )

    vars_candidates <- list.files(
        output_directory,
        pattern = paste0(
            "^",
            output_basename,
            "\\.unphased\\.vcor1\\.bin\\.vars$"
        ),
        full.names = TRUE
    )

    list(
        matrix_file = if (length(matrix_candidates) == 1L)
        {
            matrix_candidates[1]
        } else {
            NA_character_
        },
        vars_file = if (length(vars_candidates) == 1L)
        {
            vars_candidates[1]
        } else {
            NA_character_
        }
    )
}

############################################################
## Validate software and input files
############################################################

reference_pvar_file <- if (
    file.exists(
        paste0(
            reference_prefix,
            ".pvar"
        )
    )
)
{
    paste0(
        reference_prefix,
        ".pvar"
    )
} else if (
    file.exists(
        paste0(
            reference_prefix,
            ".pvar.zst"
        )
    )
)
{
    paste0(
        reference_prefix,
        ".pvar.zst"
    )
} else {
    NA_character_
}

required_reference_files <- c(
    paste0(reference_prefix, ".pgen"),
    paste0(reference_prefix, ".pvar.zst"),
    paste0(reference_prefix, ".psam")
)

missing_reference_files <- required_reference_files[
    !file.exists(required_reference_files)
]

if (length(missing_reference_files) > 0L)
{
    stop(
        paste0(
            "Missing PLINK reference files:\n",
            paste(
                missing_reference_files,
                collapse = "\n"
            )
        )
    )
}

if (!file.exists(plink2_executable))
{
    stop(
        paste0(
            "PLINK2 executable not found: ",
            plink2_executable
        )
    )
}

if (!file.exists(validation_manifest_file))
{
    stop(
        paste0(
            "Validation manifest not found: ",
            validation_manifest_file,
            "\nRun Module 08B first."
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

plink_version <- tryCatch(
    system2(
        plink2_executable,
        args = "--version",
        stdout = TRUE,
        stderr = TRUE
    ),
    error = function(e)
    {
        conditionMessage(e)
    }
)

cat("\nPLINK2:\n")
cat(
    paste(
        plink_version,
        collapse = "\n"
    ),
    "\n"
)

############################################################
## Read validated clumps
############################################################

cat("\nReading Module 08B validation manifest...\n")

validation_manifest <- fread(
    validation_manifest_file
)

required_validation_columns <- c(
    "clump_id",
    "parent_locus",
    "chromosome",
    "index_variant",
    "index_position",
    "region_start",
    "region_end",
    "validation_status",
    "region_file"
)

missing_validation_columns <- setdiff(
    required_validation_columns,
    names(validation_manifest)
)

if (length(missing_validation_columns) > 0L)
{
    stop(
        paste0(
            "Missing validation-manifest columns: ",
            paste(
                missing_validation_columns,
                collapse = ", "
            )
        )
    )
}

validated_clumps <- validation_manifest[
    validation_status == "ready"
]

if (nrow(validated_clumps) == 0L)
{
    stop(
        "No validated clumps are available for LD calculation."
    )
}

setorder(
    validated_clumps,
    chromosome,
    index_position
)

if (run_mode == "full")
{
    n_to_run <- min(
        pilot_n_clumps,
        nrow(validated_clumps)
    )

    validated_clumps <- validated_clumps[
        seq_len(n_to_run)
    ]
} else if (run_mode != "full")
{
    stop(
        "run_mode must be either 'pilot' or 'full'."
    )
}

cat(
    "Run mode:",
    run_mode,
    "\n"
)

cat(
    "Clumps selected:",
    nrow(validated_clumps),
    "\n"
)

cat(
    "Reference MAF threshold:",
    minimum_reference_maf,
    "\n"
)

cat(
    "Dense LD safety limit:",
    maximum_ld_variants,
    "variants\n"
)

############################################################
## Process clumps
############################################################

manifest_list <- vector(
    mode = "list",
    length = nrow(validated_clumps)
)

processing_start <- Sys.time()

cat("\nCalculating LD matrices...\n")

for (i in seq_len(nrow(validated_clumps)))
{
    clump_start <- Sys.time()

    metadata_row <- validated_clumps[i]

    clump_id <- metadata_row$clump_id

    chromosome_value <- normalise_chromosome(
        metadata_row$chromosome
    )

    region_file <- metadata_row$region_file

    ld_directory <- file.path(
        workspace_root,
        clump_id,
        "ld"
    )

    logs_directory <- file.path(
        workspace_root,
        clump_id,
        "logs"
    )

    dir.create(
        ld_directory,
        recursive = TRUE,
        showWarnings = FALSE
    )

    dir.create(
        logs_directory,
        recursive = TRUE,
        showWarnings = FALSE
    )

    extract_file <- file.path(
        ld_directory,
        "gwas_rsids.txt"
    )

    output_prefix <- file.path(
        ld_directory,
        "ld_matrix"
    )

    plink_log_copy <- file.path(
        logs_directory,
        "module09_plink.log"
    )

    variant_map_file <- file.path(
        ld_directory,
        "ld_variant_map.tsv"
    )

    status <- "pending"
    message <- NA_character_
    plink_exit_status <- NA_integer_
    n_region_variants <- 0L
    n_valid_rsids <- 0L
    n_unique_rsids <- 0L
    n_ld_variants <- 0L
    estimated_matrix_bytes <- NA_real_
    matrix_file <- NA_character_
    vars_file <- NA_character_
    matrix_exists <- FALSE
    vars_exists <- FALSE
    variant_map_exists <- FALSE
    output_reused <- FALSE

    existing_outputs <- locate_ld_outputs(
        output_prefix
    )

    existing_matrix_valid <- (
        !is.na(existing_outputs$matrix_file) &&
        file.exists(existing_outputs$matrix_file) &&
        file.info(existing_outputs$matrix_file)$size > 0
    )

    existing_vars_valid <- (
        !is.na(existing_outputs$vars_file) &&
        file.exists(existing_outputs$vars_file) &&
        length(
            read_variant_ids(
                existing_outputs$vars_file
            )
        ) >= 2L
    )

    if (
        !overwrite_existing &&
        existing_matrix_valid &&
        existing_vars_valid
    )
    {
        matrix_file <- existing_outputs$matrix_file
        vars_file <- existing_outputs$vars_file

        ld_variant_ids <- read_variant_ids(
            vars_file
        )

        n_ld_variants <- length(
            ld_variant_ids
        )

        expected_matrix_bytes <- 4 *
            as.numeric(n_ld_variants)^2

        observed_matrix_bytes <- file.info(
            matrix_file
        )$size

        matrix_size_valid <- (
            observed_matrix_bytes ==
            expected_matrix_bytes
        )

        if (matrix_size_valid)
        {
            status <- "ready"
            message <- "Existing LD outputs reused."
            output_reused <- TRUE
            matrix_exists <- TRUE
            vars_exists <- TRUE
            estimated_matrix_bytes <-
                expected_matrix_bytes
        }
    }

    if (status != "ready")
    {
        region_gwas <- tryCatch(
            fread(
                region_file,
                showProgress = FALSE
            ),
            error = function(e)
            {
                message <<- conditionMessage(e)
                NULL
            }
        )

        if (is.null(region_gwas))
        {
            status <- "region_read_error"
        } else {
            n_region_variants <- nrow(
                region_gwas
            )

            rsid_column <- find_column(
                names(region_gwas),
                c(
                    "rsid",
                    "rs_id",
                    "snp"
                ),
                "rsID",
                required = FALSE
            )

            variant_id_column <- find_column(
                names(region_gwas),
                c(
                    "variant_id",
                    "variant",
                    "id"
                ),
                "variant identifier"
            )

            position_column <- find_column(
                names(region_gwas),
                c(
                    "base_pair_location",
                    "position",
                    "pos",
                    "bp"
                ),
                "base-pair position"
            )

            effect_allele_column <- find_column(
                names(region_gwas),
                c(
                    "effect_allele",
                    "a1",
                    "effect"
                ),
                "effect allele",
                required = FALSE
            )

            other_allele_column <- find_column(
                names(region_gwas),
                c(
                    "other_allele",
                    "a2",
                    "non_effect_allele"
                ),
                "other allele",
                required = FALSE
            )

            if (is.na(rsid_column))
            {
                status <- "rsid_column_missing"
                message <-
                    "No rsID column was detected in the regional GWAS file."
            } else {
                valid_rsid_mask <- !is_missing_rsid(
                    region_gwas[[rsid_column]]
                )

                n_valid_rsids <- sum(
                    valid_rsid_mask
                )

            reference_variant_ids <- normalise_variant_id(
                region_gwas[[variant_id_column]]
                )

            reference_variant_ids <- reference_variant_ids[
                !is.na(reference_variant_ids) &
                reference_variant_ids != ""
                ]

            reference_variant_ids <- sort(
                unique(reference_variant_ids)
                )

            n_unique_rsids <- length(
                reference_variant_ids
            )

                if (n_unique_rsids < 2L)
                {
                    status <- "insufficient_gwas_rsids"
                    message <-
                        "Fewer than two valid unique GWAS rsIDs were available."
                } else {
                    writeLines(
                        reference_variant_ids,
                        extract_file,
                        useBytes = TRUE
                    )

                    old_outputs <- c(
                        paste0(
                            output_prefix,
                            ".unphased.vcor1.bin"
                        ),
                        paste0(
                            output_prefix,
                            ".unphased.vcor1.vars"
                        ),
                        paste0(
                            output_prefix,
                            ".log"
                        ),
                        paste0(
                            output_prefix,
                            ".nosex"
                        )
                    )

                    if (overwrite_existing)
                    {
                        unlink(
                            old_outputs,
                            force = TRUE
                        )
                    }

                    plink_arguments <- c(
                        "--pfile",
                        reference_prefix,
                        if (
                            grepl(
                                "\\.pvar\\.zst$",
                                reference_pvar_file
                            )
                        )
                        {
                            "vzs"
                        } else {
                            character()
                        },
                        "--chr",
                        chromosome_value,
                        "--from-bp",
                        as.character(
                            as.integer(metadata_row$region_start)
                        ),
                        "--to-bp",
                        as.character(
                            as.integer(metadata_row$region_end)
                        ),
                        "--extract",
                        extract_file,
                        "--force-intersect",
                        "--maf",
                        format(
                            minimum_reference_maf,
                            scientific = FALSE,
                            trim = TRUE
                        ),
                        "--max-alleles",
                        "2",
                        "--rm-dup",
                        "force-first",
                        "--r-unphased",
                        "square",
                        "bin4",
                        "ref-based",
                        "--threads",
                        as.character(plink_threads),
                        "--memory",
                        as.character(plink_memory_mb),
                        "--out",
                        output_prefix
                    )

                    plink_output <- tryCatch(
                        system2(
                            plink2_executable,
                            args = plink_arguments,
                            stdout = TRUE,
                            stderr = TRUE
                        ),
                        error = function(e)
                        {
                            message <<-
                                conditionMessage(e)

                            structure(
                                character(),
                                status = 1L
                            )
                        }
                    )

                    plink_exit_status <- attr(
                        plink_output,
                        "status"
                    )

                    if (is.null(plink_exit_status))
                    {
                        plink_exit_status <- 0L
                    }

                    writeLines(
                        plink_output,
                        plink_log_copy,
                        useBytes = TRUE
                    )

                    produced_outputs <- locate_ld_outputs(
                        output_prefix
                    )

                    matrix_file <-
                        produced_outputs$matrix_file

                    vars_file <-
                        produced_outputs$vars_file

                    matrix_exists <- (
                        !is.na(matrix_file) &&
                        file.exists(matrix_file)
                    )

                    vars_exists <- (
                        !is.na(vars_file) &&
                        file.exists(vars_file)
                    )

                    if (
                        plink_exit_status != 0L ||
                        !matrix_exists ||
                        !vars_exists
                    )
                    {
                        status <- "plink_failure"

                        if (is.na(message))
                        {
                            message <- paste0(
                                "PLINK failed or did not create both LD outputs. Exit status: ",
                                plink_exit_status
                            )
                        }
                    } else {
                        ld_variant_ids <- read_variant_ids(
                            vars_file
                        )

                        n_ld_variants <- length(
                            ld_variant_ids
                        )

                        estimated_matrix_bytes <- 4 *
                            as.numeric(n_ld_variants)^2

                        if (n_ld_variants < 2L)
                        {
                            status <- "insufficient_ld_variants"
                            message <-
                                "Fewer than two reference-matched variants remained."
                        } else if (
                            n_ld_variants >
                            maximum_ld_variants
                        )
                        {
                            status <- "ld_variant_limit_exceeded"

                            message <- paste0(
                                "LD matrix contains ",
                                n_ld_variants,
                                " variants, exceeding the safety limit of ",
                                maximum_ld_variants,
                                "."
                            )
                        } else {
                            observed_matrix_bytes <- file.info(
                                matrix_file
                            )$size

                            expected_matrix_bytes <-
                                estimated_matrix_bytes

                            if (
                                observed_matrix_bytes !=
                                expected_matrix_bytes
                            )
                            {
                                status <- "matrix_size_mismatch"

                                message <- paste0(
                                    "Observed binary matrix size was ",
                                    observed_matrix_bytes,
                                    " bytes; expected ",
                                    expected_matrix_bytes,
                                    " bytes."
                                )
                            } else {
                                region_variant_map <- data.table(
                                    rsid =
                                        as.character(
                                            region_gwas[[rsid_column]]
                                        ),
                                    variant_id =
                                        as.character(
                                            region_gwas[[variant_id_column]]
                                        ),
                                    position =
                                        as.integer(
                                            region_gwas[[position_column]]
                                        )
                                )

                                if (
                                    !is.na(
                                        effect_allele_column
                                    )
                                )
                                {
                                    region_variant_map[
                                        ,
                                        effect_allele :=
                                            as.character(
                                                region_gwas[[effect_allele_column]]
                                            )
                                    ]
                                }

                                if (
                                    !is.na(
                                        other_allele_column
                                    )
                                )
                                {
                                    region_variant_map[
                                        ,
                                        reference_variant_id :=
                                            as.character(
                                                region_gwas[[other_allele_column]]
                                            )
                                    ]
                                }

                                region_variant_map[
                                ,
                                ld_order := match(
                                    reference_variant_id,
                                    ld_variant_ids
                                    )
                                ]

                                region_variant_map <- region_variant_map[
                                    !is.na(ld_order)
                                ]

                                setorder(
                                    region_variant_map,
                                    ld_order
                                )

                                fwrite(
                                    region_variant_map,
                                    variant_map_file,
                                    sep = "\t",
                                    quote = FALSE,
                                    na = "NA"
                                )

                                variant_map_exists <- file.exists(
                                    variant_map_file
                                )

                                status <- "ready"
                                message <- "PASS"
                            }
                        }
                    }
                }
            }
        }
    }

    if (
        status == "ready" &&
        !variant_map_exists &&
        file.exists(variant_map_file)
    )
    {
        variant_map_exists <- TRUE
    }

    clump_end <- Sys.time()

    clump_elapsed_seconds <- as.numeric(
        difftime(
            clump_end,
            clump_start,
            units = "secs"
        )
    )

    manifest_list[[i]] <- data.table(
        clump_id =
            clump_id,
        parent_locus =
            metadata_row$parent_locus,
        chromosome =
            metadata_row$chromosome,
        index_variant =
            metadata_row$index_variant,
        index_position =
            metadata_row$index_position,
        region_start =
            metadata_row$region_start,
        region_end =
            metadata_row$region_end,
        n_region_variants =
            n_region_variants,
        n_valid_rsids =
            n_valid_rsids,
        n_unique_rsids =
            n_unique_rsids,
        n_ld_variants =
            n_ld_variants,
        reference_maf_threshold =
            minimum_reference_maf,
        estimated_matrix_bytes =
            estimated_matrix_bytes,
        matrix_exists =
            matrix_exists,
        vars_exists =
            vars_exists,
        variant_map_exists =
            variant_map_exists,
        output_reused =
            output_reused,
        plink_exit_status =
            plink_exit_status,
        ld_status =
            status,
        ld_message =
            message,
        elapsed_seconds =
            clump_elapsed_seconds,
        extract_file =
            normalizePath(
                extract_file,
                winslash = "/",
                mustWork = FALSE
            ),
        matrix_file =
            if (!is.na(matrix_file))
            {
                normalizePath(
                    matrix_file,
                    winslash = "/",
                    mustWork = FALSE
                )
            } else {
                NA_character_
            },
        vars_file =
            if (!is.na(vars_file))
            {
                normalizePath(
                    vars_file,
                    winslash = "/",
                    mustWork = FALSE
                )
            } else {
                NA_character_
            },
        variant_map_file =
            normalizePath(
                variant_map_file,
                winslash = "/",
                mustWork = FALSE
            ),
        plink_log =
            normalizePath(
                plink_log_copy,
                winslash = "/",
                mustWork = FALSE
            )
    )

    if (
        i %% progress_interval == 0L ||
        i == nrow(validated_clumps)
    )
    {
        elapsed_so_far <- as.numeric(
            difftime(
                Sys.time(),
                processing_start,
                units = "secs"
            )
        )

        cat(
            "Processed",
            i,
            "of",
            nrow(validated_clumps),
            "|",
            clump_id,
            "| status:",
            status,
            "| LD variants:",
            n_ld_variants,
            "| elapsed:",
            format_seconds(
                elapsed_so_far
            ),
            "\n"
        )
    }
}

############################################################
## Save manifest
############################################################

ld_manifest <- rbindlist(
    manifest_list,
    use.names = TRUE,
    fill = FALSE
)

fwrite(
    ld_manifest,
    ld_manifest_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

############################################################
## Final summary
############################################################

n_selected <- nrow(
    validated_clumps
)

n_ready <- sum(
    ld_manifest$ld_status == "ready"
)

n_failed <- n_selected - n_ready

n_reused <- sum(
    ld_manifest$output_reused
)

total_ld_variants <- sum(
    ld_manifest$n_ld_variants,
    na.rm = TRUE
)

total_matrix_bytes <- sum(
    ld_manifest$estimated_matrix_bytes,
    na.rm = TRUE
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
cat("Module 09 completed\n")
cat("============================================\n")
cat(
    "Run mode                     :",
    run_mode,
    "\n"
)
cat(
    "Selected clumps              :",
    n_selected,
    "\n"
)
cat(
    "Ready clumps                 :",
    n_ready,
    "\n"
)
cat(
    "Failed clumps                :",
    n_failed,
    "\n"
)
cat(
    "Reused outputs               :",
    n_reused,
    "\n"
)
cat(
    "Total LD variants            :",
    total_ld_variants,
    "\n"
)
cat(
    "Estimated matrix storage     :",
    round(
        total_matrix_bytes / 1024^3,
        3
    ),
    "GiB\n"
)
cat(
    "LD manifest                  :",
    ld_manifest_file,
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

if (n_failed > 0L)
{
    failure_summary <- ld_manifest[
        ld_status != "ready",
        .N,
        by = ld_status
    ]

    cat("\nFailure summary:\n")
    print(
        failure_summary
    )

    stop(
        paste0(
            "Module 09 completed with ",
            n_failed,
            " failed clumps. Review ",
            ld_manifest_file,
            " and the clump-specific PLINK logs."
        )
    )
}

cat("\nModule 09 completed successfully.\n")
