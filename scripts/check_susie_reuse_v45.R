#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly=TRUE)
i <- match("--workspace",a)
if (is.na(i) || i==length(a)) stop("Usage: --workspace <clean workspace>",call.=FALSE)
root <- normalizePath(a[[i+1L]],winslash="/",mustWork=TRUE)
setwd(root)
Sys.setenv(PIGMENTATIONATLAS_ROOT=root,PIGMENTATIONATLAS_SUSIE_MODE="pilot")
script <- "scripts/upstream/Code S15_v45.R"
output <- "results/susie_manifest_pilot_v45.tsv"
if (!file.exists(script)) stop("Missing script: ",script,call.=FALSE)
if (file.exists(output)) stop("Pilot reuse audit output exists",call.=FALSE)
rscript <- file.path(R.home("bin"),"Rscript.exe")
logdir <- "logs/susie_reuse_v45"
dir.create(logdir,recursive=TRUE,showWarnings=FALSE)
code <- suppressWarnings(system2(rscript,shQuote(script),stdout=file.path(logdir,"stdout.log"),stderr=file.path(logdir,"stderr.log"),wait=TRUE))
if (as.integer(code)!=0L || !file.exists(output)) stop("Reuse audit failed; see ",logdir,call.=FALSE)
m <- fread(output)
if (nrow(m)!=5L || sum(m$status=="ready")!=3L || sum(m$status=="warning_nonconverged")!=2L ||
    !setequal(m[status=="warning_nonconverged"]$clump_id,c("CLUMP_0004","CLUMP_0005")) ||
    any(grepl("(^|;)NA($|;)",m$qc_flags),na.rm=TRUE) ||
    any(!grepl("nonconverged",m[status=="warning_nonconverged"]$qc_flags))) stop("Reuse status/QC failed",call.=FALSE)
print(m[,.(clump_id,status,converged,n_iterations,qc_flags)])
cat("PASS SuSiE reuse QC: 3 ready, 2 nonconverged retained and labeled.\n")
