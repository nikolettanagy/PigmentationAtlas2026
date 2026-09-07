############################################################
# PigmentationAtlas
# Module 11.2
# GTEx v10 apaQTL regional extraction
#
# Strategy:
#   - Read each tissue × chromosome Parquet file once
#   - Extract variant position from variant_id
#   - Cut all publication loci on that chromosome
#   - Parse alleles and APA phenotype metadata only after
#     regional extraction
#   - Write one compressed TSV per locus × tissue
#
# Input:
#   Publication locus manifest
#   GTEx v10 apaQTL chromosome-level Parquet files
#
# Outputs:
#   Regional apaQTL files
#   Extraction manifest
#   QC summaries
#   Processing log
############################################################

suppressPackageStartupMessages({
    library(data.table)
    library(arrow)
})

cat("\n")
cat("============================================================\n")
cat("MODULE 11.2: GTEx v10 apaQTL REGION EXTRACTION\n")
cat("============================================================\n")

############################################################
# 1. Configuration
############################################################

project_dir <- "C:/Users/User/Desktop/PigmentationAtlas"
setwd(project_dir)

apaqtl_dir <- file.path(
    project_dir,
    "data",
    "GTEx_v10_apaQTL"
)

output_root <- file.path(
    project_dir,
    "results",
    "module11",
    "eqtl_regions"
)

summary_dir <- file.path(
    project_dir,
    "results",
    "module11",
    "eqtl_regions",
    "summary"
)

log_dir <- file.path(
    project_dir,
    "logs",
    "module11"
)

dir.create(
    output_root,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    summary_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    log_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

log_file <- file.path(
    log_dir,
    "module11_2_extract_apaqtl_regions.log"
)

overwrite_existing <- FALSE

# Only autosomes are currently present in the publication manifest.
allowed_chromosomes <- as.character(1:22)

tissues <- c(
    "Skin_Sun_Exposed_Lower_leg",
    "Skin_Not_Sun_Exposed_Suprapubic"
)

apaqtl_columns <- c(
    "phenotype_id",
    "variant_id",
    "tss_distance",
    "af",
    "ma_samples",
    "ma_count",
    "pval_nominal",
    "slope",
    "slope_se"
)

############################################################
# 2. Logging
############################################################

sink(
    log_file,
    append = FALSE,
    split = TRUE
)

on.exit(
    {
        while (sink.number() > 0L) {
            sink()
        }
    },
    add = TRUE
)

cat("Project directory: ", project_dir, "\n", sep = "")
cat("apaQTL directory:   ", apaqtl_dir, "\n", sep = "")
cat("Output directory:   ", output_root, "\n", sep = "")
cat("Overwrite existing: ", overwrite_existing, "\n", sep = "")

############################################################
# 3. Helper functions
############################################################

normalize_chromosome <- function(x) {

    x <- as.character(x)
    x <- trimws(x)
    x <- sub("^chr", "", x, ignore.case = TRUE)

    x
}

sanitize_filename <- function(x) {

    x <- as.character(x)
    x <- gsub("[^A-Za-z0-9_.-]", "_", x)
    x <- gsub("_+", "_", x)
    x <- sub("^_+", "", x)
    x <- sub("_+$", "", x)

    x
}

format_bytes <- function(x) {

    if (is.na(x)) {
        return(NA_character_)
    }

    units <- c("B", "KB", "MB", "GB", "TB")
    index <- 1L
    value <- as.numeric(x)

    while (value >= 1024 && index < length(units)) {
        value <- value / 1024
        index <- index + 1L
    }

    paste0(
        format(
            round(value, 2),
            nsmall = 2,
            trim = TRUE
        ),
        " ",
        units[index]
    )
}

apaqtl_file_path <- function(tissue, chromosome) {

    file.path(
        apaqtl_dir,
        paste0(
            tissue,
            ".v10.cis_apaqtl.allpairs.chr",
            chromosome,
            ".parquet"
        )
    )
}

############################################################
# 4. Find the publication locus manifest
############################################################

find_locus_manifest <- function() {

    preferred_candidates <- c(
        file.path(
            project_dir,
            "results",
            "module11",
            "module11_locus_manifest.tsv"
        ),
        file.path(
            project_dir,
            "results",
            "module11",
            "publication_locus_manifest.tsv"
        ),
        file.path(
            project_dir,
            "results",
            "module10",
            "publication_manifest.tsv"
        ),
        file.path(
            project_dir,
            "results",
            "publication_manifest.tsv"
        )
    )

    existing_preferred <- preferred_candidates[
        file.exists(preferred_candidates)
    ]

    search_candidates <- list.files(
        file.path(project_dir, "results"),
        pattern = "(publication|locus).*manifest.*\\.tsv$",
        full.names = TRUE,
        recursive = TRUE,
        ignore.case = TRUE
    )

    candidates <- unique(
        c(
            existing_preferred,
            search_candidates
        )
    )

    if (length(candidates) == 0L) {
        stop(
            paste0(
                "No candidate publication manifest was found under:\n",
                file.path(project_dir, "results")
            )
        )
    }

    required_coordinate_columns <- c(
        "chromosome",
        "region_start",
        "region_end"
    )

    for (candidate in candidates) {

        header <- tryCatch(
            fread(
                candidate,
                nrows = 0L,
                showProgress = FALSE
            ),
            error = function(e) NULL
        )

        if (is.null(header)) {
            next
        }

        has_coordinates <- all(
            required_coordinate_columns %in% names(header)
        )

        has_identifier <- any(
            c(
                "locus_id",
                "clump_id"
            ) %in% names(header)
        )

        if (has_coordinates && has_identifier) {
            return(candidate)
        }
    }

    stop(
        paste0(
            "Candidate manifest files were found, but none contained ",
            "chromosome, region_start, region_end and locus_id/clump_id."
        )
    )
}

locus_manifest_file <- find_locus_manifest()

cat("\nLocus manifest selected:\n")
cat(locus_manifest_file, "\n")

############################################################
# 5. Read and standardize locus manifest
############################################################

loci <- fread(
    locus_manifest_file,
    showProgress = FALSE
)

if (!"locus_id" %in% names(loci)) {

    if ("clump_id" %in% names(loci)) {

        loci[
            ,
            locus_id := as.character(clump_id)
        ]

    } else {

        stop(
            "Neither locus_id nor clump_id is present in the manifest."
        )
    }
}

required_locus_columns <- c(
    "locus_id",
    "chromosome",
    "region_start",
    "region_end"
)

missing_locus_columns <- setdiff(
    required_locus_columns,
    names(loci)
)

if (length(missing_locus_columns) > 0L) {

    stop(
        "Missing required locus columns: ",
        paste(
            missing_locus_columns,
            collapse = ", "
        )
    )
}

loci[
    ,
    chromosome := normalize_chromosome(chromosome)
]

loci[
    ,
    `:=`(
        locus_id = as.character(locus_id),
        region_start = suppressWarnings(
            as.integer(region_start)
        ),
        region_end = suppressWarnings(
            as.integer(region_end)
        )
    )
]

loci <- loci[
    chromosome %in% allowed_chromosomes
]

invalid_loci <- loci[
    is.na(locus_id) |
    locus_id == "" |
    is.na(region_start) |
    is.na(region_end) |
    region_start < 1L |
    region_end < region_start
]

if (nrow(invalid_loci) > 0L) {

    invalid_file <- file.path(
        summary_dir,
        "module11_2_invalid_loci.tsv"
    )

    fwrite(
        invalid_loci,
        invalid_file,
        sep = "\t",
        na = "NA"
    )

    stop(
        nrow(invalid_loci),
        " invalid loci were detected. See: ",
        invalid_file
    )
}

duplicate_ids <- loci[
    duplicated(locus_id) |
    duplicated(locus_id, fromLast = TRUE)
]

if (nrow(duplicate_ids) > 0L) {

    duplicate_file <- file.path(
        summary_dir,
        "module11_2_duplicate_locus_ids.tsv"
    )

    fwrite(
        duplicate_ids,
        duplicate_file,
        sep = "\t",
        na = "NA"
    )

    stop(
        "Duplicate locus identifiers detected. See: ",
        duplicate_file
    )
}

setorder(
    loci,
    chromosome,
    region_start,
    region_end
)

cat("\nPublication loci retained: ", nrow(loci), "\n", sep = "")
cat(
    "Chromosomes represented: ",
    paste(
        sort(
            unique(
                as.integer(loci$chromosome)
            )
        ),
        collapse = ", "
    ),
    "\n",
    sep = ""
)

############################################################
# 6. Parse chromosome-level variant positions
############################################################

add_variant_position <- function(chr_dt, expected_chr) {

    # Expected format:
    # chr11_69049973_T_G_b38

    chr_dt[
        ,
        variant_position := suppressWarnings(
            as.integer(
                sub(
                    "^chr[^_]+_([0-9]+)_.*$",
                    "\\1",
                    variant_id
                )
            )
        )
    ]

    chr_dt[
        ,
        parsed_chromosome := sub(
            "^chr([^_]+)_.*$",
            "\\1",
            variant_id
        )
    ]

    invalid <- chr_dt[
        is.na(variant_position) |
        parsed_chromosome != expected_chr
    ]

    list(
        data = chr_dt,
        invalid = invalid
    )
}

############################################################
# 7. Parse regional variant identifiers
############################################################

parse_regional_variants <- function(region_dt) {

    if (nrow(region_dt) == 0L) {
        return(region_dt)
    }

    parts <- tstrsplit(
        region_dt$variant_id,
        "_",
        fixed = TRUE,
        keep = 1:5,
        type.convert = FALSE
    )

    if (length(parts) != 5L) {

        stop(
            "Unexpected variant_id structure in regional data."
        )
    }

    region_dt[
        ,
        `:=`(
            variant_chromosome = sub(
                "^chr",
                "",
                parts[[1L]]
            ),
            variant_position_parsed = suppressWarnings(
                as.integer(parts[[2L]])
            ),
            ref = toupper(parts[[3L]]),
            alt = toupper(parts[[4L]]),
            genome_build = parts[[5L]]
        )
    ]

    inconsistent_position <- region_dt[
        is.na(variant_position_parsed) |
        variant_position_parsed != variant_position
    ]

    if (nrow(inconsistent_position) > 0L) {

        stop(
            nrow(inconsistent_position),
            " inconsistent variant positions were detected."
        )
    }

    region_dt[
        ,
        variant_position_parsed := NULL
    ]

    # GTEx/QTLtools-style variant_id:
    # slope is interpreted relative to the ALT allele.
    region_dt[
        ,
        `:=`(
            effect_allele = alt,
            other_allele = ref,
            effect_allele_frequency = as.numeric(af)
        )
    ]

    region_dt
}

############################################################
# 8. Parse APA phenotype identifiers
############################################################

parse_apa_phenotypes <- function(region_dt) {

    if (nrow(region_dt) == 0L) {
        return(region_dt)
    }

    # Example:
    # ENSG00000177963.15_chr11:209407-209811

    region_dt[
        ,
        gene_id := sub(
            "_chr[^:]+:[0-9]+-[0-9]+$",
            "",
            phenotype_id
        )
    ]

    region_dt[
        ,
        gene_id_unversioned := sub(
            "\\.[0-9]+$",
            "",
            gene_id
        )
    ]

    region_dt[
        ,
        apa_interval := sub(
            "^.*_(chr[^:]+:[0-9]+-[0-9]+)$",
            "\\1",
            phenotype_id
        )
    ]

    region_dt[
        ,
        apa_chromosome := sub(
            "^chr([^:]+):.*$",
            "\\1",
            apa_interval
        )
    ]

    region_dt[
        ,
        apa_start := suppressWarnings(
            as.integer(
                sub(
                    "^chr[^:]+:([0-9]+)-[0-9]+$",
                    "\\1",
                    apa_interval
                )
            )
        )
    ]

    region_dt[
        ,
        apa_end := suppressWarnings(
            as.integer(
                sub(
                    "^chr[^:]+:[0-9]+-([0-9]+)$",
                    "\\1",
                    apa_interval
                )
            )
        )
    ]

    region_dt
}

############################################################
# 9. Existing output summary
############################################################

summarize_existing_output <- function(output_file) {

    existing <- tryCatch(
        fread(
            output_file,
            select = c(
                "variant_id",
                "phenotype_id",
                "gene_id_unversioned"
            ),
            showProgress = FALSE
        ),
        error = function(e) NULL
    )

    if (is.null(existing)) {
        return(NULL)
    }

    data.table(
        n_rows = nrow(existing),
        n_variants = uniqueN(existing$variant_id),
        n_phenotypes = uniqueN(existing$phenotype_id),
        n_genes = uniqueN(existing$gene_id_unversioned)
    )
}

############################################################
# 10. Prepare result containers
############################################################

n_expected_combinations <- (
    nrow(loci) *
    length(tissues)
)

extraction_results <- vector(
    mode = "list",
    length = n_expected_combinations
)

chromosome_qc_results <- vector(
    mode = "list",
    length = length(tissues) *
        uniqueN(loci$chromosome)
)

result_index <- 0L
chromosome_qc_index <- 0L

############################################################
# 11. Tissue × chromosome processing
############################################################

for (tissue in tissues) {

    cat("\n")
    cat("############################################################\n")
    cat("TISSUE: ", tissue, "\n", sep = "")
    cat("############################################################\n")

    tissue_output_dir <- file.path(
        output_root,
        tissue
    )

    dir.create(
        tissue_output_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    tissue_chromosomes <- sort(
        unique(
            as.integer(loci$chromosome)
        )
    )

    for (chr_integer in tissue_chromosomes) {

        chr_value <- as.character(chr_integer)

        chr_loci <- loci[
            chromosome == chr_value
        ]

        parquet_file <- apaqtl_file_path(
            tissue,
            chr_value
        )

        chromosome_qc_index <- chromosome_qc_index + 1L

        cat("\n")
        cat(
            "------------------------------------------------------------\n"
        )
        cat(
            tissue,
            " | chr",
            chr_value,
            " | loci: ",
            nrow(chr_loci),
            "\n",
            sep = ""
        )

        if (!file.exists(parquet_file)) {

            cat("Parquet file missing: ", parquet_file, "\n", sep = "")

            chromosome_qc_results[[chromosome_qc_index]] <- data.table(
                tissue = tissue,
                chromosome = chr_value,
                parquet_file = parquet_file,
                file_size_bytes = NA_real_,
                file_size_display = NA_character_,
                read_status = "MISSING_PARQUET",
                read_seconds = NA_real_,
                n_rows = 0L,
                n_variants = 0L,
                n_phenotypes = 0L,
                n_invalid_variant_ids = 0L,
                memory_mb = NA_real_,
                error = NA_character_
            )

            for (j in seq_len(nrow(chr_loci))) {

                locus <- chr_loci[j]
                result_index <- result_index + 1L

                extraction_results[[result_index]] <- data.table(
                    locus_id = locus$locus_id,
                    tissue = tissue,
                    chromosome = chr_value,
                    region_start = locus$region_start,
                    region_end = locus$region_end,
                    extraction_status = "MISSING_PARQUET",
                    n_rows = 0L,
                    n_variants = 0L,
                    n_phenotypes = 0L,
                    n_genes = 0L,
                    minimum_p_value = NA_real_,
                    output_file = NA_character_,
                    source_parquet = parquet_file,
                    error = NA_character_
                )
            }

            next
        }

        file_info <- file.info(parquet_file)

        read_result <- tryCatch(
            {
                read_time <- system.time({

                    chr_tbl <- arrow::read_parquet(
                        parquet_file,
                        col_select = apaqtl_columns,
                        as_data_frame = TRUE
                    )

                })

                chr_dt <- as.data.table(chr_tbl)
                rm(chr_tbl)

                position_result <- add_variant_position(
                    chr_dt,
                    expected_chr = chr_value
                )

                chr_dt <- position_result$data
                invalid_variant_ids <- position_result$invalid

                list(
                    success = TRUE,
                    data = chr_dt,
                    invalid = invalid_variant_ids,
                    elapsed = unname(
                        read_time[["elapsed"]]
                    ),
                    error = NA_character_
                )
            },
            error = function(e) {

                list(
                    success = FALSE,
                    data = NULL,
                    invalid = NULL,
                    elapsed = NA_real_,
                    error = conditionMessage(e)
                )
            }
        )

        if (!read_result$success) {

            cat(
                "Chromosome read failed: ",
                read_result$error,
                "\n",
                sep = ""
            )

            chromosome_qc_results[[chromosome_qc_index]] <- data.table(
                tissue = tissue,
                chromosome = chr_value,
                parquet_file = parquet_file,
                file_size_bytes = file_info$size,
                file_size_display = format_bytes(file_info$size),
                read_status = "READ_ERROR",
                read_seconds = read_result$elapsed,
                n_rows = 0L,
                n_variants = 0L,
                n_phenotypes = 0L,
                n_invalid_variant_ids = 0L,
                memory_mb = NA_real_,
                error = read_result$error
            )

            for (j in seq_len(nrow(chr_loci))) {

                locus <- chr_loci[j]
                result_index <- result_index + 1L

                extraction_results[[result_index]] <- data.table(
                    locus_id = locus$locus_id,
                    tissue = tissue,
                    chromosome = chr_value,
                    region_start = locus$region_start,
                    region_end = locus$region_end,
                    extraction_status = "CHROMOSOME_READ_ERROR",
                    n_rows = 0L,
                    n_variants = 0L,
                    n_phenotypes = 0L,
                    n_genes = 0L,
                    minimum_p_value = NA_real_,
                    output_file = NA_character_,
                    source_parquet = parquet_file,
                    error = read_result$error
                )
            }

            next
        }

        chr_dt <- read_result$data
        invalid_variant_ids <- read_result$invalid

        chromosome_memory_mb <- round(
            as.numeric(object.size(chr_dt)) / 1024^2,
            2
        )

        cat(
            "Rows loaded:      ",
            format(nrow(chr_dt), big.mark = ","),
            "\n",
            sep = ""
        )

        cat(
            "Unique variants:  ",
            format(
                uniqueN(chr_dt$variant_id),
                big.mark = ","
            ),
            "\n",
            sep = ""
        )

        cat(
            "APA phenotypes:   ",
            format(
                uniqueN(chr_dt$phenotype_id),
                big.mark = ","
            ),
            "\n",
            sep = ""
        )

        cat(
            "Memory:           ",
            chromosome_memory_mb,
            " MB\n",
            sep = ""
        )

        cat(
            "Read time:        ",
            round(read_result$elapsed, 2),
            " sec\n",
            sep = ""
        )

        n_invalid_variant_ids <- nrow(invalid_variant_ids)

        if (n_invalid_variant_ids > 0L) {

            invalid_output <- file.path(
                summary_dir,
                paste0(
                    "invalid_variant_ids_",
                    tissue,
                    "_chr",
                    chr_value,
                    ".tsv.gz"
                )
            )

            fwrite(
                invalid_variant_ids,
                invalid_output,
                sep = "\t",
                na = "NA",
                compress = "gzip"
            )

            cat(
                "Invalid variant IDs: ",
                n_invalid_variant_ids,
                "\n",
                sep = ""
            )

            chr_dt <- chr_dt[
                !is.na(variant_position) &
                parsed_chromosome == chr_value
            ]
        }

        chromosome_qc_results[[chromosome_qc_index]] <- data.table(
            tissue = tissue,
            chromosome = chr_value,
            parquet_file = normalizePath(
                parquet_file,
                winslash = "/",
                mustWork = FALSE
            ),
            file_size_bytes = file_info$size,
            file_size_display = format_bytes(file_info$size),
            read_status = "READY",
            read_seconds = read_result$elapsed,
            n_rows = nrow(chr_dt),
            n_variants = uniqueN(chr_dt$variant_id),
            n_phenotypes = uniqueN(chr_dt$phenotype_id),
            n_invalid_variant_ids = n_invalid_variant_ids,
            memory_mb = chromosome_memory_mb,
            error = NA_character_
        )

        # parsed_chromosome is not needed after QC.
        chr_dt[
            ,
            parsed_chromosome := NULL
        ]

        setkey(
            chr_dt,
            variant_position
        )

        ######################################################
        # Extract every locus on this chromosome
        ######################################################

        for (j in seq_len(nrow(chr_loci))) {

            locus <- chr_loci[j]
            result_index <- result_index + 1L

            locus_id <- as.character(locus$locus_id)
            region_start <- as.integer(locus$region_start)
            region_end <- as.integer(locus$region_end)

            output_stub <- sanitize_filename(
                paste0(
                    locus_id,
                    "_chr",
                    chr_value,
                    "_",
                    region_start,
                    "_",
                    region_end
                )
            )

            output_file <- file.path(
                tissue_output_dir,
                paste0(
                    output_stub,
                    ".apaqtl.tsv.gz"
                )
            )

            cat(
                "  ",
                locus_id,
                " | ",
                region_start,
                "-",
                region_end,
                sep = ""
            )

            if (
                file.exists(output_file) &&
                !overwrite_existing
            ) {

                existing_summary <- summarize_existing_output(
                    output_file
                )

                if (!is.null(existing_summary)) {

                    cat(
                        " | EXISTING | ",
                        format(
                            existing_summary$n_rows,
                            big.mark = ","
                        ),
                        " rows\n",
                        sep = ""
                    )

                    extraction_results[[result_index]] <- data.table(
                        locus_id = locus_id,
                        tissue = tissue,
                        chromosome = chr_value,
                        region_start = region_start,
                        region_end = region_end,
                        extraction_status = "EXISTING",
                        n_rows = existing_summary$n_rows,
                        n_variants = existing_summary$n_variants,
                        n_phenotypes = existing_summary$n_phenotypes,
                        n_genes = existing_summary$n_genes,
                        minimum_p_value = NA_real_,
                        output_file = normalizePath(
                            output_file,
                            winslash = "/",
                            mustWork = FALSE
                        ),
                        source_parquet = normalizePath(
                            parquet_file,
                            winslash = "/",
                            mustWork = FALSE
                        ),
                        error = NA_character_
                    )

                    next
                }
            }

            extraction <- tryCatch(
                {
                    # data.table indexed range extraction
                    region_dt <- chr_dt[
                        variant_position >= region_start &
                        variant_position <= region_end
                    ]

                    # Explicit copy prevents regional modifications from
                    # affecting the chromosome-level cache.
                    region_dt <- copy(region_dt)

                    if (nrow(region_dt) == 0L) {

                        list(
                            success = TRUE,
                            data = region_dt,
                            error = NA_character_
                        )

                    } else {

                        region_dt <- parse_regional_variants(
                            region_dt
                        )

                        region_dt <- parse_apa_phenotypes(
                            region_dt
                        )

                        region_dt[
                            ,
                            `:=`(
                                locus_id = locus_id,
                                tissue = tissue,
                                chromosome = chr_value,
                                region_start = region_start,
                                region_end = region_end
                            )
                        ]

                        preferred_columns <- c(
                            "locus_id",
                            "tissue",
                            "chromosome",
                            "region_start",
                            "region_end",
                            "phenotype_id",
                            "gene_id",
                            "gene_id_unversioned",
                            "apa_interval",
                            "apa_chromosome",
                            "apa_start",
                            "apa_end",
                            "variant_id",
                            "variant_chromosome",
                            "variant_position",
                            "ref",
                            "alt",
                            "effect_allele",
                            "other_allele",
                            "effect_allele_frequency",
                            "tss_distance",
                            "ma_samples",
                            "ma_count",
                            "pval_nominal",
                            "slope",
                            "slope_se",
                            "genome_build"
                        )

                        setcolorder(
                            region_dt,
                            preferred_columns
                        )

                        setorder(
                            region_dt,
                            phenotype_id,
                            variant_position,
                            ref,
                            alt
                        )

                        fwrite(
                            region_dt,
                            output_file,
                            sep = "\t",
                            na = "NA",
                            quote = FALSE,
                            compress = "gzip"
                        )

                        list(
                            success = TRUE,
                            data = region_dt,
                            error = NA_character_
                        )
                    }
                },
                error = function(e) {

                    list(
                        success = FALSE,
                        data = NULL,
                        error = conditionMessage(e)
                    )
                }
            )

            if (!extraction$success) {

                cat(
                    " | ERROR: ",
                    extraction$error,
                    "\n",
                    sep = ""
                )

                extraction_results[[result_index]] <- data.table(
                    locus_id = locus_id,
                    tissue = tissue,
                    chromosome = chr_value,
                    region_start = region_start,
                    region_end = region_end,
                    extraction_status = "ERROR",
                    n_rows = 0L,
                    n_variants = 0L,
                    n_phenotypes = 0L,
                    n_genes = 0L,
                    minimum_p_value = NA_real_,
                    output_file = NA_character_,
                    source_parquet = normalizePath(
                        parquet_file,
                        winslash = "/",
                        mustWork = FALSE
                    ),
                    error = extraction$error
                )

                next
            }

            region_dt <- extraction$data

            if (nrow(region_dt) == 0L) {

                cat(" | NO_APAQTL_ROWS\n")

                extraction_results[[result_index]] <- data.table(
                    locus_id = locus_id,
                    tissue = tissue,
                    chromosome = chr_value,
                    region_start = region_start,
                    region_end = region_end,
                    extraction_status = "NO_APAQTL_ROWS",
                    n_rows = 0L,
                    n_variants = 0L,
                    n_phenotypes = 0L,
                    n_genes = 0L,
                    minimum_p_value = NA_real_,
                    output_file = NA_character_,
                    source_parquet = normalizePath(
                        parquet_file,
                        winslash = "/",
                        mustWork = FALSE
                    ),
                    error = NA_character_
                )

                rm(region_dt)
                next
            }

            region_summary <- region_dt[
                ,
                .(
                    n_rows = .N,
                    n_variants = uniqueN(variant_id),
                    n_phenotypes = uniqueN(phenotype_id),
                    n_genes = uniqueN(gene_id_unversioned),
                    minimum_p_value = suppressWarnings(
                        min(
                            pval_nominal,
                            na.rm = TRUE
                        )
                    )
                )
            ]

            if (
                !is.finite(
                    region_summary$minimum_p_value
                )
            ) {
                region_summary$minimum_p_value <- NA_real_
            }

            cat(
                " | EXTRACTED | ",
                format(
                    region_summary$n_rows,
                    big.mark = ","
                ),
                " rows | ",
                region_summary$n_variants,
                " variants | ",
                region_summary$n_phenotypes,
                " phenotypes\n",
                sep = ""
            )

            extraction_results[[result_index]] <- data.table(
                locus_id = locus_id,
                tissue = tissue,
                chromosome = chr_value,
                region_start = region_start,
                region_end = region_end,
                extraction_status = "EXTRACTED",
                n_rows = region_summary$n_rows,
                n_variants = region_summary$n_variants,
                n_phenotypes = region_summary$n_phenotypes,
                n_genes = region_summary$n_genes,
                minimum_p_value = region_summary$minimum_p_value,
                output_file = normalizePath(
                    output_file,
                    winslash = "/",
                    mustWork = FALSE
                ),
                source_parquet = normalizePath(
                    parquet_file,
                    winslash = "/",
                    mustWork = FALSE
                ),
                error = NA_character_
            )

            rm(
                region_dt,
                extraction,
                region_summary
            )
        }

        ######################################################
        # Release chromosome cache
        ######################################################

        rm(
            chr_dt,
            read_result,
            invalid_variant_ids
        )

        invisible(gc())

        cat(
            "Chromosome chr",
            chr_value,
            " completed and memory released.\n",
            sep = ""
        )
    }
}

############################################################
# 12. Combine extraction results
############################################################

extraction_results <- extraction_results[
    seq_len(result_index)
]

extraction_manifest <- rbindlist(
    extraction_results,
    fill = TRUE
)

chromosome_qc_results <- chromosome_qc_results[
    seq_len(chromosome_qc_index)
]

chromosome_qc <- rbindlist(
    chromosome_qc_results,
    fill = TRUE
)

setorder(
    extraction_manifest,
    tissue,
    chromosome,
    region_start,
    region_end
)

setorder(
    chromosome_qc,
    tissue,
    chromosome
)

############################################################
# 13. Summary tables
############################################################

status_summary <- extraction_manifest[
    ,
    .(
        n_locus_tissue_combinations = .N,
        total_rows = sum(
            n_rows,
            na.rm = TRUE
        ),
        total_reported_variants = sum(
            n_variants,
            na.rm = TRUE
        ),
        total_reported_phenotypes = sum(
            n_phenotypes,
            na.rm = TRUE
        )
    ),
    by = .(
        tissue,
        extraction_status
    )
]

setorder(
    status_summary,
    tissue,
    extraction_status
)

tissue_summary <- extraction_manifest[
    ,
    .(
        n_loci = .N,
        n_extracted = sum(
            extraction_status %in% c(
                "EXTRACTED",
                "EXISTING"
            )
        ),
        n_empty = sum(
            extraction_status == "NO_APAQTL_ROWS"
        ),
        n_errors = sum(
            extraction_status %in% c(
                "ERROR",
                "CHROMOSOME_READ_ERROR",
                "MISSING_PARQUET"
            )
        ),
        total_rows = sum(
            n_rows,
            na.rm = TRUE
        ),
        median_variants_per_locus = median(
            n_variants[
                extraction_status %in% c(
                    "EXTRACTED",
                    "EXISTING"
                )
            ],
            na.rm = TRUE
        ),
        median_phenotypes_per_locus = median(
            n_phenotypes[
                extraction_status %in% c(
                    "EXTRACTED",
                    "EXISTING"
                )
            ],
            na.rm = TRUE
        )
    ),
    by = tissue
]

############################################################
# 14. Write manifests and QC outputs
############################################################

manifest_output <- file.path(
    summary_dir,
    "module11_2_apaqtl_region_manifest.tsv"
)

status_output <- file.path(
    summary_dir,
    "module11_2_extraction_status_summary.tsv"
)

tissue_output <- file.path(
    summary_dir,
    "module11_2_tissue_summary.tsv"
)

chromosome_qc_output <- file.path(
    summary_dir,
    "module11_2_chromosome_qc.tsv"
)

fwrite(
    extraction_manifest,
    manifest_output,
    sep = "\t",
    na = "NA"
)

fwrite(
    status_summary,
    status_output,
    sep = "\t",
    na = "NA"
)

fwrite(
    tissue_summary,
    tissue_output,
    sep = "\t",
    na = "NA"
)

fwrite(
    chromosome_qc,
    chromosome_qc_output,
    sep = "\t",
    na = "NA"
)

############################################################
# 15. Final report
############################################################

cat("\n")
cat("============================================================\n")
cat("MODULE 11.2 COMPLETE\n")
cat("============================================================\n")

cat(
    "Publication loci:              ",
    nrow(loci),
    "\n",
    sep = ""
)

cat(
    "Expected locus × tissue pairs: ",
    n_expected_combinations,
    "\n",
    sep = ""
)

cat(
    "Manifest rows produced:        ",
    nrow(extraction_manifest),
    "\n",
    sep = ""
)

cat(
    "Extracted:                     ",
    extraction_manifest[
        extraction_status == "EXTRACTED",
        .N
    ],
    "\n",
    sep = ""
)

cat(
    "Existing:                      ",
    extraction_manifest[
        extraction_status == "EXISTING",
        .N
    ],
    "\n",
    sep = ""
)

cat(
    "Empty:                         ",
    extraction_manifest[
        extraction_status == "NO_APAQTL_ROWS",
        .N
    ],
    "\n",
    sep = ""
)

cat(
    "Errors:                        ",
    extraction_manifest[
        extraction_status %in% c(
            "ERROR",
            "CHROMOSOME_READ_ERROR",
            "MISSING_PARQUET"
        ),
        .N
    ],
    "\n",
    sep = ""
)

cat("\nStatus summary:\n")
print(status_summary)

cat("\nTissue summary:\n")
print(tissue_summary)

cat("\nOutputs:\n")
cat("Region manifest: ", manifest_output, "\n", sep = "")
cat("Status summary:  ", status_output, "\n", sep = "")
cat("Tissue summary:  ", tissue_output, "\n", sep = "")
cat("Chromosome QC:   ", chromosome_qc_output, "\n", sep = "")
cat("Regional files:  ", output_root, "\n", sep = "")
cat("Log:             ", log_file, "\n", sep = "")

cat("============================================================\n")
