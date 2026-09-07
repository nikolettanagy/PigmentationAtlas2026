############################################################
# PigmentationAtlas
#
# Module: 07B
# Name: Create clump workspaces
#
# Purpose:
# Creates a separate analysis workspace for every PLINK
# clump listed in the master metadata table.
#
# Each workspace contains:
#   metadata.tsv
#   gwas/
#   ld/
#   susie/
#   coloc/
#   figures/
#   logs/
#   tmp/
#
# Input:
#   results/clump_master_metadata.tsv
#
# Outputs:
#   results/clumps/CLUMP_xxxx/
#   results/clump_workspace_manifest.tsv
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
cat("PigmentationAtlas Module 07B\n")
cat("Create clump workspaces\n")
cat("============================================\n")

############################################################
## File paths
############################################################

metadata_file <- file.path(
    "results",
    "clump_master_metadata.tsv"
)

workspace_root <- file.path(
    "results",
    "clumps"
)

manifest_file <- file.path(
    "results",
    "clump_workspace_manifest.tsv"
)

############################################################
## Parameters
############################################################

subdirectories <- c(
    "gwas",
    "ld",
    "susie",
    "coloc",
    "figures",
    "logs",
    "tmp"
)

overwrite_metadata <- TRUE

############################################################
## Check input
############################################################

if (!file.exists(metadata_file))
{
    stop(
        paste0(
            "Input file not found: ",
            metadata_file
        )
    )
}

if (!dir.exists("results"))
{
    stop(
        "The results directory does not exist."
    )
}

############################################################
## Read master metadata
############################################################

cat("\nReading master metadata...\n")

clump_metadata <- fread(
    metadata_file
)

cat(
    "Clumps loaded:",
    nrow(clump_metadata),
    "\n"
)

############################################################
## Validate required columns
############################################################

required_columns <- c(
    "clump_id",
    "parent_locus",
    "chromosome",
    "index_variant",
    "index_position",
    "index_p",
    "raw_start",
    "raw_end",
    "raw_span_bp",
    "region_start",
    "region_end",
    "region_span_bp",
    "n_variants",
    "n_sp2_variants",
    "padding_bp",
    "region_size_mb",
    "analysis_status",
    "susie_status",
    "coloc_status",
    "figure_status"
)

missing_columns <- setdiff(
    required_columns,
    names(clump_metadata)
)

if (length(missing_columns) > 0L)
{
    stop(
        paste0(
            "Missing required metadata columns: ",
            paste(
                missing_columns,
                collapse = ", "
            )
        )
    )
}

############################################################
## Basic metadata QC
############################################################

if (nrow(clump_metadata) == 0L)
{
    stop(
        "The master metadata table contains no rows."
    )
}

if (anyNA(clump_metadata$clump_id))
{
    stop(
        "Missing clump_id values detected."
    )
}

if (anyDuplicated(clump_metadata$clump_id))
{
    duplicated_ids <- unique(
        clump_metadata[
            duplicated(clump_id) |
            duplicated(clump_id, fromLast = TRUE),
            clump_id
        ]
    )

    stop(
        paste0(
            "Duplicated clump IDs detected: ",
            paste(
                duplicated_ids,
                collapse = ", "
            )
        )
    )
}

valid_clump_ids <- grepl(
    "^CLUMP_[0-9]{4,}$",
    clump_metadata$clump_id
)

if (any(!valid_clump_ids))
{
    stop(
        paste0(
            "Invalid clump ID format: ",
            paste(
                clump_metadata$clump_id[
                    !valid_clump_ids
                ],
                collapse = ", "
            )
        )
    )
}

if (
    anyNA(clump_metadata$region_start) ||
    anyNA(clump_metadata$region_end)
)
{
    stop(
        "Missing region coordinates detected."
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
        "Invalid region coordinates detected."
    )
}

############################################################
## Create workspace root
############################################################

if (!dir.exists(workspace_root))
{
    dir.create(
        workspace_root,
        recursive = TRUE,
        showWarnings = FALSE
    )

    cat(
        "\nCreated workspace root:",
        workspace_root,
        "\n"
    )
} else {
    cat(
        "\nWorkspace root already exists:",
        workspace_root,
        "\n"
    )
}

############################################################
## Prepare manifest
############################################################

manifest_list <- vector(
    mode = "list",
    length = nrow(clump_metadata)
)

############################################################
## Create workspaces
############################################################

cat("\nCreating clump workspaces...\n")

for (i in seq_len(nrow(clump_metadata)))
{
    metadata_row <- clump_metadata[i]

    clump_id <- metadata_row$clump_id

    clump_directory <- file.path(
        workspace_root,
        clump_id
    )

    clump_directory_existed <- dir.exists(
        clump_directory
    )

    if (!clump_directory_existed)
    {
        dir.create(
            clump_directory,
            recursive = TRUE,
            showWarnings = FALSE
        )
    }

    subdirectory_paths <- file.path(
        clump_directory,
        subdirectories
    )

    for (directory_path in subdirectory_paths)
    {
        if (!dir.exists(directory_path))
        {
            dir.create(
                directory_path,
                recursive = TRUE,
                showWarnings = FALSE
            )
        }
    }

    clump_metadata_file <- file.path(
        clump_directory,
        "metadata.tsv"
    )

    metadata_file_existed <- file.exists(
        clump_metadata_file
    )

    if (
        overwrite_metadata ||
        !metadata_file_existed
    )
    {
        fwrite(
            metadata_row,
            clump_metadata_file,
            sep = "\t",
            quote = FALSE,
            na = "NA"
        )
    }

    expected_paths <- c(
        clump_directory,
        subdirectory_paths
    )

    all_directories_present <- all(
        dir.exists(expected_paths)
    )

    metadata_present <- file.exists(
        clump_metadata_file
    )

    manifest_list[[i]] <- data.table(
        clump_id = clump_id,
        parent_locus =
            metadata_row$parent_locus,
        chromosome =
            metadata_row$chromosome,
        index_variant =
            metadata_row$index_variant,
        index_position =
            metadata_row$index_position,
        workspace_path =
            normalizePath(
                clump_directory,
                winslash = "/",
                mustWork = FALSE
            ),
        workspace_preexisting =
            clump_directory_existed,
        metadata_preexisting =
            metadata_file_existed,
        metadata_written =
            overwrite_metadata ||
            !metadata_file_existed,
        metadata_present =
            metadata_present,
        subdirectories_present =
            all_directories_present,
        workspace_status =
            if (
                metadata_present &&
                all_directories_present
            )
            {
                "ready"
            } else {
                "incomplete"
            }
    )

    if (
        i %% 100L == 0L ||
        i == nrow(clump_metadata)
    )
    {
        cat(
            "Processed",
            i,
            "of",
            nrow(clump_metadata),
            "workspaces\n"
        )
    }
}

############################################################
## Combine manifest
############################################################

workspace_manifest <- rbindlist(
    manifest_list,
    use.names = TRUE,
    fill = FALSE
)

setorder(
    workspace_manifest,
    chromosome,
    index_position
)

############################################################
## Detailed QC
############################################################

cat("\nRunning workspace quality control...\n")

expected_clumps <- nrow(
    clump_metadata
)

actual_clump_directories <- list.dirs(
    workspace_root,
    full.names = FALSE,
    recursive = FALSE
)

actual_clump_directories <- actual_clump_directories[
    grepl(
        "^CLUMP_[0-9]{4,}$",
        actual_clump_directories
    )
]

n_clump_directories <- length(
    actual_clump_directories
)

metadata_paths <- file.path(
    workspace_root,
    clump_metadata$clump_id,
    "metadata.tsv"
)

n_metadata_files <- sum(
    file.exists(metadata_paths)
)

subdirectory_counts <- vapply(
    subdirectories,
    FUN = function(subdirectory_name)
    {
        paths <- file.path(
            workspace_root,
            clump_metadata$clump_id,
            subdirectory_name
        )

        sum(
            dir.exists(paths)
        )
    },
    FUN.VALUE = integer(1)
)

missing_workspace_ids <- setdiff(
    clump_metadata$clump_id,
    actual_clump_directories
)

unexpected_workspace_ids <- setdiff(
    actual_clump_directories,
    clump_metadata$clump_id
)

incomplete_manifest_rows <- workspace_manifest[
    workspace_status != "ready"
]

############################################################
## Stop on QC failure
############################################################

if (n_clump_directories < expected_clumps)
{
    stop(
        paste0(
            "QC failure: expected ",
            expected_clumps,
            " clump directories, but found ",
            n_clump_directories,
            "."
        )
    )
}

if (length(missing_workspace_ids) > 0L)
{
    stop(
        paste0(
            "QC failure: missing workspace directories: ",
            paste(
                missing_workspace_ids,
                collapse = ", "
            )
        )
    )
}

if (n_metadata_files != expected_clumps)
{
    stop(
        paste0(
            "QC failure: expected ",
            expected_clumps,
            " metadata files, but found ",
            n_metadata_files,
            "."
        )
    )
}

if (any(subdirectory_counts != expected_clumps))
{
    failed_subdirectories <- names(
        subdirectory_counts[
            subdirectory_counts !=
                expected_clumps
        ]
    )

    stop(
        paste0(
            "QC failure: incomplete subdirectories: ",
            paste(
                failed_subdirectories,
                collapse = ", "
            )
        )
    )
}

if (nrow(incomplete_manifest_rows) > 0L)
{
    stop(
        paste0(
            "QC failure: incomplete workspaces detected: ",
            paste(
                incomplete_manifest_rows$clump_id,
                collapse = ", "
            )
        )
    )
}

############################################################
## Verify individual metadata files
############################################################

cat("Checking individual metadata files...\n")

metadata_validation <- vector(
    mode = "logical",
    length = expected_clumps
)

for (i in seq_len(expected_clumps))
{
    expected_row <- clump_metadata[i]

    individual_metadata_file <- file.path(
        workspace_root,
        expected_row$clump_id,
        "metadata.tsv"
    )

    individual_metadata <- fread(
        individual_metadata_file
    )

    metadata_validation[i] <- (
        nrow(individual_metadata) == 1L &&
        individual_metadata$clump_id[1] ==
            expected_row$clump_id &&
        individual_metadata$index_variant[1] ==
            expected_row$index_variant &&
        individual_metadata$region_start[1] ==
            expected_row$region_start &&
        individual_metadata$region_end[1] ==
            expected_row$region_end
    )
}

if (any(!metadata_validation))
{
    invalid_ids <- clump_metadata$clump_id[
        !metadata_validation
    ]

    stop(
        paste0(
            "QC failure: invalid individual metadata files: ",
            paste(
                invalid_ids,
                collapse = ", "
            )
        )
    )
}

############################################################
## Save manifest
############################################################

fwrite(
    workspace_manifest,
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

cat("\n")
cat("============================================\n")
cat("Module 07B completed successfully\n")
cat("============================================\n")
cat(
    "Expected clumps        :",
    expected_clumps,
    "\n"
)
cat(
    "Clump directories      :",
    n_clump_directories,
    "\n"
)
cat(
    "metadata.tsv files     :",
    n_metadata_files,
    "\n"
)

for (subdirectory_name in subdirectories)
{
    cat(
        sprintf(
            "%-22s: %d\n",
            paste0(
                subdirectory_name,
                " directories"
            ),
            subdirectory_counts[
                subdirectory_name
            ]
        )
    )
}

cat(
    "Preexisting workspaces :",
    sum(
        workspace_manifest$
            workspace_preexisting
    ),
    "\n"
)
cat(
    "New workspaces         :",
    sum(
        !workspace_manifest$
            workspace_preexisting
    ),
    "\n"
)
cat(
    "Ready workspaces       :",
    sum(
        workspace_manifest$
            workspace_status ==
            "ready"
    ),
    "\n"
)
cat(
    "Incomplete workspaces  :",
    sum(
        workspace_manifest$
            workspace_status !=
            "ready"
    ),
    "\n"
)
cat(
    "Unexpected directories :",
    length(
        unexpected_workspace_ids
    ),
    "\n"
)
cat(
    "Workspace root         :",
    workspace_root,
    "\n"
)
cat(
    "Manifest file          :",
    manifest_file,
    "\n"
)
cat(
    "Finished at            :",
    format(
        end_time,
        "%Y-%m-%d %H:%M:%S"
    ),
    "\n"
)
cat(
    "Elapsed time           :",
    round(
        elapsed_seconds,
        1
    ),
    "seconds\n"
)
cat("============================================\n")