############################################################
# PigmentationAtlas
# Module 11.0
# Colocalization input audit and eQTL inventory
############################################################

suppressPackageStartupMessages({
    library(data.table)
})

cat("\n")
cat("============================================================\n")
cat("MODULE 11.0: COLOCALIZATION INPUT AUDIT\n")
cat("============================================================\n")

project_dir <- "C:/Users/User/Desktop/PigmentationAtlas"
setwd(project_dir)

publication_manifest_file <-
    "results/susie_manifest_publication.tsv"

production_manifest_file <-
    "results/susie_manifest_production.tsv"

eqtl_root <- "data/eqtl"

output_dir <- "results/module11"
log_dir <- "logs/module11"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
    log_dir,
    "module11_0_input_audit.log"
)

sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

############################################################
# Helper functions
############################################################

first_existing_column <- function(dt, candidates) {

    hit <- candidates[candidates %in% names(dt)]

    if (length(hit) == 0L) {
        return(NA_character_)
    }

    hit[1L]
}

safe_min <- function(x) {

    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]

    if (length(x) == 0L) {
        return(NA_real_)
    }

    min(x)
}

safe_max <- function(x) {

    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]

    if (length(x) == 0L) {
        return(NA_real_)
    }

    max(x)
}

safe_read_header <- function(path) {

    result <- tryCatch(
        {
            dt <- fread(
                path,
                nrows = 5L,
                showProgress = FALSE
            )

            list(
                readable = TRUE,
                n_columns = ncol(dt),
                columns = paste(names(dt), collapse = ";"),
                preview_rows = nrow(dt),
                error = NA_character_
            )
        },
        error = function(e) {
            list(
                readable = FALSE,
                n_columns = NA_integer_,
                columns = NA_character_,
                preview_rows = NA_integer_,
                error = conditionMessage(e)
            )
        }
    )

    result
}

############################################################
# 1. Publication manifest
############################################################

if (!file.exists(publication_manifest_file)) {
    stop(
        "Publication manifest not found: ",
        publication_manifest_file
    )
}

publication <- fread(publication_manifest_file)

cat("\nPublication manifest loaded\n")
cat("Rows:    ", nrow(publication), "\n", sep = "")
cat("Columns: ", ncol(publication), "\n", sep = "")

cat("\nPublication manifest columns:\n")
print(names(publication))

if ("publication_status" %in% names(publication)) {

    cat("\nPublication status distribution:\n")
    print(publication[, .N, by = publication_status][order(publication_status)])
}

############################################################
# 2. Identify essential columns
############################################################

clump_col <- first_existing_column(
    publication,
    c(
        "clump_id",
        "locus_id",
        "region_id"
    )
)

chr_col <- first_existing_column(
    publication,
    c(
        "chromosome",
        "chr",
        "CHR"
    )
)

start_col <- first_existing_column(
    publication,
    c(
        "region_start",
        "locus_start",
        "start",
        "START"
    )
)

end_col <- first_existing_column(
    publication,
    c(
        "region_end",
        "locus_end",
        "end",
        "END"
    )
)

lead_col <- first_existing_column(
    publication,
    c(
        "gwas_lead_rsid",
        "lead_rsid",
        "lead_snp",
        "gwas_lead",
        "top_gwas_rsid"
    )
)

top_pip_col <- first_existing_column(
    publication,
    c(
        "top_pip_rsid",
        "top_pip_snp",
        "susie_top_rsid"
    )
)

cat("\nDetected essential columns:\n")

detected_columns <- data.table(
    field = c(
        "locus identifier",
        "chromosome",
        "region start",
        "region end",
        "GWAS lead",
        "top PIP SNP"
    ),
    column = c(
        clump_col,
        chr_col,
        start_col,
        end_col,
        lead_col,
        top_pip_col
    )
)

print(detected_columns)

############################################################
# 3. Check required fields
############################################################

required_fields <- data.table(
    field = c(
        "locus_id",
        "chromosome",
        "region_start",
        "region_end"
    ),
    detected_column = c(
        clump_col,
        chr_col,
        start_col,
        end_col
    )
)

required_fields[
    ,
    available := !is.na(detected_column)
]

cat("\nRequired field audit:\n")
print(required_fields)

if (any(!required_fields$available)) {

    cat("\nWARNING:\n")
    cat(
        "Some essential locus columns could not be identified.\n",
        "The manifest will still be audited, but Module 11.1 will require\n",
        "explicit column mapping.\n",
        sep = ""
    )
}

############################################################
# 4. Build standardized locus table
############################################################

locus_table <- data.table(
    source_row = seq_len(nrow(publication))
)

if (!is.na(clump_col)) {
    locus_table[, locus_id := as.character(publication[[clump_col]])]
} else {
    locus_table[, locus_id := sprintf("LOCUS_%04d", source_row)]
}

if (!is.na(chr_col)) {
    locus_table[, chromosome := publication[[chr_col]]]
} else {
    locus_table[, chromosome := NA_integer_]
}

if (!is.na(start_col)) {
    locus_table[, region_start := publication[[start_col]]]
} else {
    locus_table[, region_start := NA_integer_]
}

if (!is.na(end_col)) {
    locus_table[, region_end := publication[[end_col]]]
} else {
    locus_table[, region_end := NA_integer_]
}

if (!is.na(lead_col)) {
    locus_table[, gwas_lead_rsid := as.character(publication[[lead_col]])]
} else {
    locus_table[, gwas_lead_rsid := NA_character_]
}

if (!is.na(top_pip_col)) {
    locus_table[, top_pip_rsid := as.character(publication[[top_pip_col]])]
} else {
    locus_table[, top_pip_rsid := NA_character_]
}

locus_table[
    ,
    chromosome := suppressWarnings(
        as.integer(
            gsub("^chr", "", as.character(chromosome), ignore.case = TRUE)
        )
    )
]

locus_table[
    ,
    `:=`(
        region_start = suppressWarnings(as.integer(region_start)),
        region_end = suppressWarnings(as.integer(region_end))
    )
]

locus_table[
    ,
    region_width := region_end - region_start + 1L
]

locus_table[
    ,
    coordinate_complete :=
        !is.na(chromosome) &
        !is.na(region_start) &
        !is.na(region_end) &
        region_start <= region_end
]

cat("\nStandardized locus table:\n")
cat("Loci:                 ", nrow(locus_table), "\n", sep = "")
cat(
    "Complete coordinates:  ",
    sum(locus_table$coordinate_complete),
    "\n",
    sep = ""
)
cat(
    "Incomplete coordinates:",
    sum(!locus_table$coordinate_complete),
    "\n",
    sep = ""
)

cat("\nLoci by chromosome:\n")
print(
    locus_table[
        ,
        .N,
        by = chromosome
    ][order(chromosome)]
)

cat("\nRegion-width summary:\n")
print(summary(locus_table$region_width))

############################################################
# 5. Check duplicate locus IDs and coordinates
############################################################

duplicate_locus_ids <- locus_table[
    !is.na(locus_id),
    .N,
    by = locus_id
][N > 1L]

duplicate_coordinates <- locus_table[
    coordinate_complete == TRUE,
    .N,
    by = .(
        chromosome,
        region_start,
        region_end
    )
][N > 1L]

cat("\nDuplicate locus identifiers: ")
cat(nrow(duplicate_locus_ids), "\n")

cat("Duplicate coordinate intervals: ")
cat(nrow(duplicate_coordinates), "\n")

############################################################
# 6. eQTL file inventory
############################################################

if (!dir.exists(eqtl_root)) {
    dir.create(eqtl_root, recursive = TRUE)
}

eqtl_files <- list.files(
    eqtl_root,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
)

supported_pattern <-
    "\\.(tsv|txt|csv|gz|bgz|parquet|rds|rdata|xlsx)$"

eqtl_files <- eqtl_files[
    grepl(
        supported_pattern,
        eqtl_files,
        ignore.case = TRUE
    )
]

cat("\neQTL root directory:\n")
cat(normalizePath(eqtl_root, winslash = "/", mustWork = FALSE), "\n")

cat("\nCandidate eQTL files found: ")
cat(length(eqtl_files), "\n")

if (length(eqtl_files) == 0L) {

    eqtl_inventory <- data.table(
        file_path = character(),
        file_name = character(),
        extension = character(),
        size_bytes = numeric(),
        readable_with_fread = logical(),
        n_columns_preview = integer(),
        columns = character(),
        read_error = character()
    )

    cat("\nNo eQTL files were found under data/eqtl.\n")

} else {

    file_info <- file.info(eqtl_files)

    eqtl_inventory <- data.table(
        file_path = normalizePath(
            eqtl_files,
            winslash = "/",
            mustWork = FALSE
        ),
        file_name = basename(eqtl_files),
        extension = tolower(tools::file_ext(eqtl_files)),
        size_bytes = file_info$size
    )

    text_like <- eqtl_inventory$extension %in%
        c("tsv", "txt", "csv", "gz", "bgz")

    eqtl_inventory[
        ,
        `:=`(
            readable_with_fread = NA,
            n_columns_preview = NA_integer_,
            columns = NA_character_,
            read_error = NA_character_
        )
    ]

    for (i in which(text_like)) {

        audit <- safe_read_header(eqtl_inventory$file_path[i])

        eqtl_inventory$readable_with_fread[i] <- audit$readable
        eqtl_inventory$n_columns_preview[i] <- audit$n_columns
        eqtl_inventory$columns[i] <- audit$columns
        eqtl_inventory$read_error[i] <- audit$error
    }

    cat("\neQTL file extensions:\n")
    print(
        eqtl_inventory[
            ,
            .N,
            by = extension
        ][order(-N)]
    )

    cat("\nLargest candidate eQTL files:\n")
    print(
        eqtl_inventory[
            order(-size_bytes),
            .(
                file_name,
                extension,
                size_bytes,
                readable_with_fread,
                n_columns_preview
            )
        ][1:min(.N, 20L)]
    )
}

############################################################
# 7. Search for known eQTL column names
############################################################

column_patterns <- list(
    variant = c(
        "variant_id",
        "variant",
        "snp",
        "rsid",
        "rs_id"
    ),
    chromosome = c(
        "chromosome",
        "chr"
    ),
    position = c(
        "position",
        "pos",
        "base_pair_location"
    ),
    effect_allele = c(
        "effect_allele",
        "alt",
        "a1"
    ),
    other_allele = c(
        "other_allele",
        "ref",
        "a2"
    ),
    beta = c(
        "beta",
        "slope",
        "nes",
        "effect"
    ),
    standard_error = c(
        "standard_error",
        "se",
        "slope_se"
    ),
    p_value = c(
        "p_value",
        "pval_nominal",
        "pval",
        "p"
    ),
    maf = c(
        "maf",
        "minor_allele_frequency"
    ),
    sample_size = c(
        "n",
        "sample_size"
    ),
    gene = c(
        "gene_id",
        "gene",
        "molecular_trait_id",
        "phenotype_id"
    )
)

if (nrow(eqtl_inventory) > 0L) {

    eqtl_inventory[
        ,
        detected_fields := vapply(
            columns,
            function(column_string) {

                if (is.na(column_string)) {
                    return(NA_character_)
                }

                x <- tolower(
                    unlist(
                        strsplit(column_string, ";", fixed = TRUE)
                    )
                )

                fields <- names(column_patterns)[
                    vapply(
                        column_patterns,
                        function(patterns) {
                            any(tolower(patterns) %in% x)
                        },
                        logical(1L)
                    )
                ]

                paste(fields, collapse = ";")
            },
            character(1L)
        )
    ]
}

############################################################
# 8. Write outputs
############################################################

fwrite(
    locus_table,
    file.path(
        output_dir,
        "module11_locus_manifest.tsv"
    ),
    sep = "\t",
    na = "NA"
)

fwrite(
    required_fields,
    file.path(
        output_dir,
        "module11_required_field_audit.tsv"
    ),
    sep = "\t",
    na = "NA"
)

fwrite(
    duplicate_locus_ids,
    file.path(
        output_dir,
        "module11_duplicate_locus_ids.tsv"
    ),
    sep = "\t",
    na = "NA"
)

fwrite(
    duplicate_coordinates,
    file.path(
        output_dir,
        "module11_duplicate_coordinates.tsv"
    ),
    sep = "\t",
    na = "NA"
)

fwrite(
    eqtl_inventory,
    file.path(
        output_dir,
        "module11_eqtl_file_inventory.tsv"
    ),
    sep = "\t",
    na = "NA"
)

############################################################
# 9. Final summary
############################################################

summary_table <- data.table(
    metric = c(
        "publication_loci",
        "loci_with_complete_coordinates",
        "loci_with_incomplete_coordinates",
        "duplicate_locus_ids",
        "duplicate_coordinate_intervals",
        "candidate_eqtl_files"
    ),
    value = c(
        nrow(locus_table),
        sum(locus_table$coordinate_complete),
        sum(!locus_table$coordinate_complete),
        nrow(duplicate_locus_ids),
        nrow(duplicate_coordinates),
        nrow(eqtl_inventory)
    )
)

fwrite(
    summary_table,
    file.path(
        output_dir,
        "module11_0_summary.tsv"
    ),
    sep = "\t"
)

cat("\n============================================================\n")
cat("MODULE 11.0 COMPLETE\n")
cat("============================================================\n")

print(summary_table)

cat("\nOutputs:\n")
cat(
    "Locus manifest:       ",
    file.path(output_dir, "module11_locus_manifest.tsv"),
    "\n",
    sep = ""
)
cat(
    "Required-field audit: ",
    file.path(output_dir, "module11_required_field_audit.tsv"),
    "\n",
    sep = ""
)
cat(
    "eQTL inventory:       ",
    file.path(output_dir, "module11_eqtl_file_inventory.tsv"),
    "\n",
    sep = ""
)
cat(
    "Summary:              ",
    file.path(output_dir, "module11_0_summary.tsv"),
    "\n",
    sep = ""
)
cat(
    "Log:                  ",
    log_file,
    "\n",
    sep = ""
)

cat("============================================================\n")