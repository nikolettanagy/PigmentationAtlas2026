############################################################
# PigmentationAtlas
# Module 11.1
# GTEx v10 apaQTL Parquet registry
############################################################

suppressPackageStartupMessages({
    library(data.table)
    library(arrow)
})

cat("\n")
cat("============================================================\n")
cat("MODULE 11.1: GTEx v10 apaQTL PARQUET REGISTRY\n")
cat("============================================================\n")

project_dir <- "C:/Users/User/Desktop/PigmentationAtlas"
setwd(project_dir)

apaqtl_dir <- "data/GTEx_v10_apaQTL"
output_dir <- "results/module11/eqtl_registry"
log_dir <- "logs/module11"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
    log_dir,
    "module11_1_apaqtl_parquet_registry.log"
)

sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

############################################################
# Discover files
############################################################

files <- list.files(
    apaqtl_dir,
    pattern = "\\.parquet$",
    full.names = TRUE,
    recursive = FALSE
)

cat("\nParquet files found: ", length(files), "\n", sep = "")

if (length(files) == 0L) {
    stop("No Parquet files found in: ", apaqtl_dir)
}

############################################################
# Parse filename metadata
############################################################

parse_tissue <- function(filename) {

    if (grepl(
        "^Skin_Sun_Exposed_Lower_leg",
        filename
    )) {
        return("Skin_Sun_Exposed_Lower_leg")
    }

    if (grepl(
        "^Skin_Not_Sun_Exposed_Suprapubic",
        filename
    )) {
        return("Skin_Not_Sun_Exposed_Suprapubic")
    }

    NA_character_
}

parse_chromosome <- function(filename) {

    match <- regmatches(
        filename,
        regexpr(
            "chr([0-9]+|X)\\.parquet$",
            filename
        )
    )

    if (length(match) == 0L || !nzchar(match)) {
        return(NA_character_)
    }

    sub(
        "^chr",
        "",
        sub("\\.parquet$", "", match)
    )
}

############################################################
# Inspect each file
############################################################

registry_list <- vector(
    "list",
    length(files)
)

for (i in seq_along(files)) {

    path <- files[i]
    filename <- basename(path)
    info <- file.info(path)

    cat(
        "[",
        i,
        "/",
        length(files),
        "] ",
        filename,
        "\n",
        sep = ""
    )

    result <- tryCatch(
        {
            ds <- open_dataset(
                path,
                format = "parquet"
            )

            schema_names <- names(ds)

            preview <- ds |>
                head(5) |>
                collect()

            list(
                readable = TRUE,
                columns = schema_names,
                preview_rows = nrow(preview),
                error = NA_character_
            )
        },
        error = function(e) {
            list(
                readable = FALSE,
                columns = character(),
                preview_rows = NA_integer_,
                error = conditionMessage(e)
            )
        }
    )

    registry_list[[i]] <- data.table(
        file_id = sprintf("APAQTL_%03d", i),
        file_name = filename,
        file_path = normalizePath(
            path,
            winslash = "/",
            mustWork = FALSE
        ),
        tissue = parse_tissue(filename),
        chromosome = parse_chromosome(filename),
        genome_build = "GRCh38",
        source = "GTEx",
        release = "v10",
        molecular_trait = "alternative_polyadenylation",
        size_bytes = info$size,
        size_mb = round(
            info$size / 1024^2,
            2
        ),
        readable = result$readable,
        n_columns = length(result$columns),
        columns = paste(
            result$columns,
            collapse = ";"
        ),
        preview_rows = result$preview_rows,
        read_error = result$error
    )
}

registry <- rbindlist(
    registry_list,
    fill = TRUE
)

############################################################
# Detect important fields
############################################################

detect_column <- function(column_string, candidates) {

    if (
        is.na(column_string) ||
        !nzchar(column_string)
    ) {
        return(NA_character_)
    }

    columns <- strsplit(
        column_string,
        ";",
        fixed = TRUE
    )[[1L]]

    lower_columns <- tolower(columns)

    hit <- match(
        tolower(candidates),
        lower_columns
    )

    hit <- hit[!is.na(hit)]

    if (length(hit) == 0L) {
        return(NA_character_)
    }

    columns[hit[1L]]
}

registry[
    ,
    variant_id_col := vapply(
        columns,
        detect_column,
        character(1L),
        candidates = c(
            "variant_id",
            "variant",
            "snp_id"
        )
    )
]

registry[
    ,
    trait_id_col := vapply(
        columns,
        detect_column,
        character(1L),
        candidates = c(
            "molecular_trait_id",
            "phenotype_id",
            "gene_id",
            "apa_id",
            "phenotype"
        )
    )
]

registry[
    ,
    chromosome_col := vapply(
        columns,
        detect_column,
        character(1L),
        candidates = c(
            "chromosome",
            "chr",
            "chrom"
        )
    )
]

registry[
    ,
    position_col := vapply(
        columns,
        detect_column,
        character(1L),
        candidates = c(
            "position",
            "pos",
            "base_pair_location"
        )
    )
]

registry[
    ,
    beta_col := vapply(
        columns,
        detect_column,
        character(1L),
        candidates = c(
            "beta",
            "slope",
            "effect_size",
            "nes"
        )
    )
]

registry[
    ,
    se_col := vapply(
        columns,
        detect_column,
        character(1L),
        candidates = c(
            "standard_error",
            "se",
            "slope_se",
            "beta_se"
        )
    )
]

registry[
    ,
    p_value_col := vapply(
        columns,
        detect_column,
        character(1L),
        candidates = c(
            "p_value",
            "pval_nominal",
            "pval",
            "p"
        )
    )
]

registry[
    ,
    maf_col := vapply(
        columns,
        detect_column,
        character(1L),
        candidates = c(
            "maf",
            "minor_allele_frequency",
            "effect_allele_frequency",
            "eaf"
        )
    )
]

registry[
    ,
    core_fields_complete :=
        !is.na(variant_id_col) &
        !is.na(trait_id_col) &
        !is.na(beta_col) &
        !is.na(p_value_col)
]

registry[
    ,
    registry_status := fifelse(
        !readable,
        "UNREADABLE",
        fifelse(
            core_fields_complete,
            "READY",
            "MISSING_CORE_FIELDS"
        )
    )
]

############################################################
# Completeness checks
############################################################

expected_chromosomes <- c(
    as.character(1:22),
    "X"
)

expected_grid <- CJ(
    tissue = c(
        "Skin_Sun_Exposed_Lower_leg",
        "Skin_Not_Sun_Exposed_Suprapubic"
    ),
    chromosome = expected_chromosomes
)

coverage <- merge(
    expected_grid,
    registry[
        ,
        .(
            tissue,
            chromosome,
            file_name,
            readable,
            registry_status
        )
    ],
    by = c(
        "tissue",
        "chromosome"
    ),
    all.x = TRUE
)

coverage[
    ,
    file_present := !is.na(file_name)
]

############################################################
# Summaries
############################################################

cat("\nFiles by tissue:\n")
print(
    registry[
        ,
        .N,
        by = tissue
    ]
)

cat("\nRegistry status:\n")
print(
    registry[
        ,
        .N,
        by = registry_status
    ]
)

cat("\nChromosome coverage:\n")
print(
    coverage[
        ,
        .(
            expected = .N,
            present = sum(file_present),
            missing = sum(!file_present)
        ),
        by = tissue
    ]
)

cat("\nDetected schemas:\n")
print(
    unique(
        registry[
            ,
            .(
                n_columns,
                columns
            )
        ]
    )
)

############################################################
# Write outputs
############################################################

fwrite(
    registry,
    file.path(
        output_dir,
        "module11_apaqtl_parquet_registry.tsv"
    ),
    sep = "\t",
    na = "NA"
)

fwrite(
    coverage,
    file.path(
        output_dir,
        "module11_apaqtl_chromosome_coverage.tsv"
    ),
    sep = "\t",
    na = "NA"
)

schema_summary <- registry[
    ,
    .N,
    by = .(
        n_columns,
        columns,
        variant_id_col,
        trait_id_col,
        beta_col,
        se_col,
        p_value_col,
        maf_col
    )
]

fwrite(
    schema_summary,
    file.path(
        output_dir,
        "module11_apaqtl_schema_summary.tsv"
    ),
    sep = "\t",
    na = "NA"
)

############################################################
# Final report
############################################################

cat("\n============================================================\n")
cat("MODULE 11.1 COMPLETE\n")
cat("============================================================\n")

cat("Parquet files: ", nrow(registry), "\n", sep = "")
cat(
    "Readable:      ",
    registry[readable == TRUE, .N],
    "\n",
    sep = ""
)
cat(
    "Ready:         ",
    registry[registry_status == "READY", .N],
    "\n",
    sep = ""
)

cat("\nOutputs:\n")
cat(
    "Registry: ",
    file.path(
        output_dir,
        "module11_apaqtl_parquet_registry.tsv"
    ),
    "\n",
    sep = ""
)
cat(
    "Coverage: ",
    file.path(
        output_dir,
        "module11_apaqtl_chromosome_coverage.tsv"
    ),
    "\n",
    sep = ""
)
cat(
    "Schema:   ",
    file.path(
        output_dir,
        "module11_apaqtl_schema_summary.tsv"
    ),
    "\n",
    sep = ""
)
cat(
    "Log:      ",
    log_file,
    "\n",
    sep = ""
)

cat("============================================================\n")