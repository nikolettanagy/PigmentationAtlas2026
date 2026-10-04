#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly=TRUE)
i <- match("--workspace",a)
if (is.na(i) || i==length(a)) stop("Usage: --workspace <clean workspace>",call.=FALSE)
root <- normalizePath(a[[i+1L]],winslash="/",mustWork=TRUE)
setwd(root)
Sys.setenv(PIGMENTATIONATLAS_ROOT=root)
script <- "scripts/upstream/Code S20.R"
loci_path <- "results/module11/module11_locus_manifest.tsv"
registry_path <- "results/module11/eqtl_registry_v49/module11_apaqtl_parquet_registry.tsv"
manifest_path <- "results/module11/eqtl_regions/summary/module11_2_apaqtl_region_manifest.tsv"
if (!all(file.exists(c(script,loci_path,registry_path)))) stop("Missing S20 script, locus manifest, or v49 registry",call.=FALSE)
loci <- fread(loci_path)
registry <- fread(registry_path)
if (nrow(loci)!=265L || anyDuplicated(loci$locus_id) || nrow(registry)!=46L || !all(registry$registry_status=="READY")) stop("Input counts/status invalid",call.=FALSE)
if (file.exists(manifest_path)) stop("Region manifest exists; review before rerun",call.=FALSE)
logdir <- "logs/apaqtl_regions_v50"
dir.create(logdir,recursive=TRUE,showWarnings=FALSE)
rscript <- file.path(R.home("bin"),"Rscript.exe")
code <- suppressWarnings(system2(rscript,shQuote(script),stdout=file.path(logdir,"stdout.log"),stderr=file.path(logdir,"stderr.log"),wait=TRUE))
if (as.integer(code)!=0L || !file.exists(manifest_path)) stop("Module 11.2 failed; see ",logdir,call.=FALSE)
m <- fread(manifest_path)
if (nrow(m)!=530L || anyDuplicated(paste(m$locus_id,m$tissue)) || !setequal(m$locus_id,loci$locus_id) || any(!m$extraction_status %in% c("EXTRACTED","EXISTING","NO_APAQTL_ROWS"))) {
  print(m[, .N, by=extraction_status][order(extraction_status)])
  stop("Regional extraction QC failed; see ",manifest_path,call.=FALSE)
}
print(m[, .N, by=extraction_status][order(extraction_status)])
cat("PASS Module 11.2: 530 locus-tissue combinations audited; no extraction errors.\n")
