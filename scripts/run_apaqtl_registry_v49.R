#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly=TRUE)
i <- match("--workspace",a)
if (is.na(i) || i==length(a)) stop("Usage: --workspace <clean workspace>",call.=FALSE)
root <- normalizePath(a[[i+1L]],winslash="/",mustWork=TRUE)
setwd(root)
Sys.setenv(PIGMENTATIONATLAS_ROOT=root)
script <- "scripts/upstream/Code S19_v49.R"
input <- "data/GTEx_v10_apaQTL"
output <- "results/module11/eqtl_registry_v49/module11_apaqtl_parquet_registry.tsv"
if (!file.exists(script) || !dir.exists(input)) stop("Missing S19 script or apaQTL inputs",call.=FALSE)
if (!requireNamespace("arrow",quietly=TRUE) || !requireNamespace("dplyr",quietly=TRUE)) stop("arrow and dplyr packages are required",call.=FALSE)
if (length(list.files(input,pattern="[.]parquet$"))!=46L) stop("Expected 46 input Parquet files",call.=FALSE)
if (file.exists(output)) stop("v49 registry exists; refusing overwrite",call.=FALSE)
logdir <- "logs/apaqtl_registry_v49"
dir.create(logdir,recursive=TRUE,showWarnings=FALSE)
rscript <- file.path(R.home("bin"),"Rscript.exe")
code <- suppressWarnings(system2(rscript,shQuote(script),stdout=file.path(logdir,"stdout.log"),stderr=file.path(logdir,"stderr.log"),wait=TRUE))
if (as.integer(code)!=0L || !file.exists(output)) stop("Registry script failed; see ",logdir,call.=FALSE)
r <- fread(output)
if (nrow(r)!=46L || sum(r$readable==TRUE)!=46L || sum(r$registry_status=="READY")!=46L || anyDuplicated(r$file_name)) {
  print(r[, .N, by=registry_status][order(registry_status)])
  stop("Registry QC failed; see ",output,call.=FALSE)
}
cat("PASS apaQTL registry: 46/46 readable and READY.\n")
