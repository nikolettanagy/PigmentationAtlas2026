#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly=TRUE)
i <- match("--workspace",a)
if (is.na(i) || i==length(a)) stop("Usage: --workspace <clean workspace>",call.=FALSE)
root <- normalizePath(a[[i+1L]],winslash="/",mustWork=TRUE)
setwd(root)
Sys.setenv(PIGMENTATIONATLAS_ROOT=root,PIGMENTATIONATLAS_SUSIE_MODE="full")
old_path <- "results/susie_manifest.tsv"
new_path <- "results/susie_manifest_psd_v46.tsv"
script <- "scripts/upstream/Code S15_psd_v46.R"
if (!file.exists(script) || !file.exists(old_path)) stop("Missing v46 script or original full manifest",call.=FALSE)
old <- fread(old_path)
if (nrow(old)!=1325L || sum(old$status=="failed")!=430L || sum(old$status=="ready")!=714L || sum(old$status=="warning_nonconverged")!=181L) stop("Unexpected old manifest status counts",call.=FALSE)
psd_ids <- old[status=="failed" & grepl("minimum eigenvalue",reason)]$clump_id
few_ids <- old[status=="failed" & grepl("Too few LD variants: 2 < 10",reason)]$clump_id
if (length(psd_ids)!=429L || length(few_ids)!=1L) stop("Unexpected failure classes",call.=FALSE)
psd_values <- as.numeric(sub(".*minimum eigenvalue = ","",old[clump_id %in% psd_ids]$reason))
if (anyNA(psd_values) || min(psd_values)< -1e-5 || max(psd_values)>= -1e-6) stop("PSD values outside justified correction range",call.=FALSE)
if (file.exists(new_path)) stop("v46 manifest exists; inspect before rerun",call.=FALSE)
logdir <- "logs/susie_psd_v46"
dir.create(logdir,recursive=TRUE,showWarnings=FALSE)
rscript <- file.path(R.home("bin"),"Rscript.exe")
cat("Starting v46: 429 near-PSD retries, 895 existing fits reused, one two-variant exclusion.\n")
code <- suppressWarnings(system2(rscript,shQuote(script),stdout=file.path(logdir,"stdout.log"),stderr=file.path(logdir,"stderr.log"),wait=TRUE))
if (as.integer(code)!=0L || !file.exists(new_path)) stop("v46 run failed; see ",logdir,call.=FALSE)
new <- fread(new_path)
if (nrow(new)!=1325L || anyDuplicated(new$clump_id)) stop("v46 manifest incomplete",call.=FALSE)
print(new[, .N, by=status][order(status)])
remaining <- new[status=="failed"]$clump_id
if (!identical(remaining,few_ids) || any(new[clump_id %in% psd_ids]$status=="failed")) stop("Unexpected remaining failures; inspect ",new_path,call.=FALSE)
if (any(grepl("(^|;)NA($|;)",new$qc_flags),na.rm=TRUE)) stop("Malformed QC flags remain",call.=FALSE)
cat("PASS v46: 429 PSD failures resolved; one two-variant clump remains excluded. Original manifest preserved.\n")
