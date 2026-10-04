#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly=TRUE)
i <- match("--workspace",a)
if (is.na(i) || i==length(a)) stop("Usage: --workspace <clean workspace>",call.=FALSE)
root <- normalizePath(a[[i+1L]],winslash="/",mustWork=TRUE)
setwd(root)
Sys.setenv(PIGMENTATIONATLAS_ROOT=root,PIGMENTATIONATLAS_SUSIE_MODE="pilot")
input <- "results/clump_ld_manifest_recovered_v42.tsv"
output <- "results/susie_manifest_pilot_v43.tsv"
script <- "scripts/upstream/Code S15_v43.R"
if (!file.exists(script) || !file.exists(input)) stop("Missing SuSiE script or recovered LD manifest",call.=FALSE)
ld <- fread(input)
if (nrow(ld)!=1325L || any(ld$ld_status!="ready") || anyDuplicated(ld$clump_id)) stop("LD manifest is not 1325/1325 ready",call.=FALSE)
if (!requireNamespace("susieR",quietly=TRUE)) stop("susieR missing",call.=FALSE)
if (file.exists(output)) stop("Pilot manifest exists; review before rerun",call.=FALSE)
rscript <- file.path(R.home("bin"),"Rscript.exe")
logdir <- "logs/susie_pilot_v43"
dir.create(logdir,recursive=TRUE,showWarnings=FALSE)
code <- suppressWarnings(system2(rscript,shQuote(script),stdout=file.path(logdir,"stdout.log"),stderr=file.path(logdir,"stderr.log"),wait=TRUE))
if (as.integer(code)!=0L || !file.exists(output)) stop("S15 failed; see ",logdir,call.=FALSE)
m <- fread(output)
if (nrow(m)!=5L) stop("Pilot must have exactly five rows",call.=FALSE)
print(m[, .N, by=status][order(status)])
if (any(!m$status %in% c("ready","warning_nonconverged"))) stop("SuSiE pilot has failed clumps; inspect ",output,call.=FALSE)
cat("PASS SuSiE pilot: five clumps processed; review QC and diagnostic files before full run.\n")
