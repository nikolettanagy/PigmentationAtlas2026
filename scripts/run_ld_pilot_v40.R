#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE)
i <- match("--workspace", args)
if (is.na(i) || i == length(args)) stop("Usage: --workspace <clean workspace>", call.=FALSE)
root <- normalizePath(args[[i+1L]], winslash="/", mustWork=TRUE)
setwd(root)
Sys.setenv(PIGMENTATIONATLAS_ROOT=root, PIGMENTATIONATLAS_LD_MODE="pilot")
script <- "scripts/upstream/Code S14_v40.R"
if (!file.exists(script)) stop("Missing ", script, call.=FALSE)
manifest <- "results/clump_gwas_validation_manifest.tsv"
if (!file.exists(manifest)) stop("S13 manifest missing", call.=FALSE)
x <- fread(manifest)
if (nrow(x)!=1325L || !all(x$validation_status=="ready")) stop("S13 validation incomplete", call.=FALSE)
if (file.exists("results/clump_ld_manifest.tsv")) stop("LD manifest already exists; review it before proceeding", call.=FALSE)
for (p in c("reference/1000G_GRCh38/1000G_GRCh38_EUR_uniqueID.pgen", "reference/1000G_GRCh38/1000G_GRCh38_EUR_uniqueID.pvar.zst", "reference/1000G_GRCh38/1000G_GRCh38_EUR_uniqueID.psam", "tools/plink2/plink2.exe")) if (!file.exists(p)) stop("Missing input: ",p,call.=FALSE)
logdir <- "logs/ld_pilot_v40"
dir.create(logdir, recursive=TRUE,showWarnings=FALSE)
rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) stop("Rscript.exe missing",call.=FALSE)
code <- suppressWarnings(system2(rscript, shQuote(script), stdout=file.path(logdir,"stdout.log"), stderr=file.path(logdir,"stderr.log"), wait=TRUE))
p <- "results/clump_ld_manifest.tsv"
if (!file.exists(p)) stop("LD run produced no manifest; see ",logdir,call.=FALSE)
m <- fread(p)
if (as.integer(code)!=0L || nrow(m)!=5L) stop("Pilot failed or row count mismatch. Statuses: ",paste(m$ld_status,collapse=", "),"; see ",logdir,call.=FALSE)
for (j in seq_len(nrow(m))) {
  row <- m[j]
  if (row$ld_status!="ready") next
  ids <- readLines(row$vars_file, warn=FALSE)
  bytes <- file.info(row$matrix_file)$size
  map <- fread(row$variant_map_file)
  if (length(ids)<2L || is.na(bytes) || bytes != 4*length(ids)^2 || nrow(map)!=length(ids) || anyDuplicated(map$ld_order) || !identical(sort(map$ld_order), seq_along(ids)) || !identical(as.character(map$reference_variant_id[order(map$ld_order)]),ids)) stop("LD/map QC failed: ",row$clump_id,call.=FALSE)
}
print(m[, .N, by=ld_status][order(ld_status)])
if (!all(m$ld_status=="ready")) stop("Some pilot clumps failed; inspect ",p," and ",logdir,call.=FALSE)
cat("PASS LD pilot: 5 ready clumps; matrix and variant order checked.\n")
