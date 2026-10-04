#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly=TRUE)
i <- match("--workspace",a)
if (is.na(i) || i==length(a)) stop("Usage: --workspace <clean workspace>",call.=FALSE)
root <- normalizePath(a[[i+1L]],winslash="/",mustWork=TRUE)
setwd(root)
Sys.setenv(PIGMENTATIONATLAS_ROOT=root,PIGMENTATIONATLAS_SUSIE_MODE="pilot")
script <- "scripts/upstream/Code S15_retry_v44.R"
pilot_path <- "results/susie_manifest_pilot_v43.tsv"
retry_path <- "results/susie_manifest_retry_v44.tsv"
backup <- "results/susie_pilot_v43_backup"
if (!file.exists(script) || !file.exists(pilot_path)) stop("Missing retry script or pilot manifest",call.=FALSE)
if (file.exists(retry_path) || dir.exists(backup)) stop("Retry or backup already exists; refusing overwrite",call.=FALSE)
pilot <- fread(pilot_path)
targets <- c("CLUMP_0004","CLUMP_0005")
if (nrow(pilot)!=5L || !setequal(pilot[status=="warning_nonconverged"]$clump_id,targets) || !all(pilot[clump_id %in% targets]$n_iterations==500L)) stop("Pilot status changed; review before retry",call.=FALSE)
for (id in targets) {
  source_dir <- file.path("results","clumps",id,"susie")
  files <- file.path(source_dir,c("susie_fit.rds","susie_pip.tsv","susie_credible_sets.tsv","susie_summary.tsv","susie_diagnostics.tsv"))
  if (!all(file.exists(files))) stop("Missing pilot outputs: ",id,call.=FALSE)
}
for (id in targets) {
  source_dir <- file.path("results","clumps",id,"susie")
  files <- file.path(source_dir,c("susie_fit.rds","susie_pip.tsv","susie_credible_sets.tsv","susie_summary.tsv","susie_diagnostics.tsv"))
  destination <- file.path(backup,id)
  dir.create(destination,recursive=TRUE,showWarnings=FALSE)
  if (!all(file.copy(files,destination,overwrite=FALSE))) stop("Backup failed: ",id,call.=FALSE)
  cat("BACKED UP",id,"\n")
}
logdir <- "logs/susie_retry_v44"
dir.create(logdir,recursive=TRUE,showWarnings=FALSE)
rscript <- file.path(R.home("bin"),"Rscript.exe")
code <- suppressWarnings(system2(rscript,shQuote(script),stdout=file.path(logdir,"stdout.log"),stderr=file.path(logdir,"stderr.log"),wait=TRUE))
if (as.integer(code)!=0L || !file.exists(retry_path)) stop("Targeted SuSiE retry failed; see ",logdir,call.=FALSE)
retry <- fread(retry_path)
if (nrow(retry)!=2L || !setequal(retry$clump_id,targets)) stop("Retry manifest mismatch",call.=FALSE)
print(retry[,.(clump_id,status,converged,n_iterations,n_credible_sets,max_pip,qc_flags)])
cat("Retry output: ",file.path(root,retry_path),"\n",sep="")
if (any(retry$status=="failed")) stop("Targeted retry has failed clumps; review output",call.=FALSE)
