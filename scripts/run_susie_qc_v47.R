#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly=TRUE)
i <- match("--workspace",a)
if (is.na(i) || i==length(a)) stop("Usage: --workspace <clean workspace>",call.=FALSE)
root <- normalizePath(a[[i+1L]],winslash="/",mustWork=TRUE)
setwd(root)
Sys.setenv(PIGMENTATIONATLAS_ROOT=root)
input <- "results/susie_manifest_psd_v46.tsv"
script <- "scripts/upstream/Code S17_v47.R"
production <- "results/susie_manifest_production_v47.tsv"
publication <- "results/susie_manifest_publication_v47.tsv"
if (!file.exists(script) || !file.exists(input)) stop("Missing QC script or v46 manifest",call.=FALSE)
if (file.exists(production) || file.exists(publication)) stop("v47 QC outputs already exist",call.=FALSE)
m <- fread(input)
if (nrow(m)!=1325L || anyDuplicated(m$clump_id) || sum(m$status=="failed")!=1L || sum(m$status=="warning_nonconverged")!=337L || sum(m$status=="ready")!=987L) stop("Unexpected SuSiE status counts",call.=FALSE)
rscript <- file.path(R.home("bin"),"Rscript.exe")
logdir <- "logs/susie_qc_v47"
dir.create(logdir,recursive=TRUE,showWarnings=FALSE)
code <- suppressWarnings(system2(rscript,shQuote(script),stdout=file.path(logdir,"stdout.log"),stderr=file.path(logdir,"stderr.log"),wait=TRUE))
if (as.integer(code)!=0L || !file.exists(production) || !file.exists(publication)) stop("Module 10.3 failed; see ",logdir,call.=FALSE)
p <- fread(production)
u <- fread(publication)
if (nrow(p)!=1325L || anyDuplicated(p$clump_id) || nrow(u)!=sum(p$include_in_publication) || anyDuplicated(u$clump_id)) stop("QC manifest row audit failed",call.=FALSE)
if (any(p[status=="failed"]$include_in_publication) || any(p[status=="warning_nonconverged"]$include_in_publication) || any(u$converged!=TRUE) || any(u$reconstruction_status!="ready")) stop("Invalid publication inclusion",call.=FALSE)
print(p[, .N, by=publication_status][order(publication_status)])
cat("PASS production QC: 1325 rows; publication rows=",nrow(u),"; nonconverged/failed excluded.\n",sep="")
