#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly=TRUE)
i <- match("--workspace", a)
if (is.na(i) || i==length(a)) stop("Usage: Rscript run_apaqtl_catalogue_v52.R --workspace PATH",call.=FALSE)
root <- normalizePath(a[[i+1L]],winslash="/",mustWork=TRUE)
setwd(root)
Sys.setenv(PIGMENTATIONATLAS_ROOT=root)
base <- "results/module11/harmonization/catalogue"
pre <- file.path(base,"module11_3B_2_test_status.tsv")
queue <- file.path(base,"module11_3B_work_queue.tsv")
status <- file.path(base,"module11_3B_build_status.tsv")
final <- file.path(base,"module11_3B_harmonized_apaQTL.tsv.gz")
source_script <- "scripts/upstream/Code S24.R"
reader <- "scripts/module11_3B_2_read_and_qc_apaqtl.R"
for (p in c(pre,queue,source_script,reader)) if (!file.exists(p)) stop("Missing prerequisite: ",p,call.=FALSE)
x <- fread(pre)
if (nrow(x)!=1L || !"final_status" %in% names(x) || x$final_status[[1L]]!="PASS") stop("S23 did not PASS",call.=FALSE)
if (nrow(fread(queue))!=516L) stop("Work queue must have 516 rows",call.=FALSE)
check <- function() {
  if (!file.exists(status) || !file.exists(final) || file.info(final)$size<=0) return(FALSE)
  z <- fread(status)
  fields <- c("final_status","expected_files","processed_files","failed_files","valid_chunks","final_catalogue_exists")
  nrow(z)==1L && all(fields %in% names(z)) && z$final_status[[1L]]=="PASS" &&
    z$expected_files[[1L]]==516L && z$processed_files[[1L]]==516L &&
    z$failed_files[[1L]]==0L && z$valid_chunks[[1L]]==516L &&
    isTRUE(as.logical(z$final_catalogue_exists[[1L]]))
}
if (check()) {cat("PASS S24 existing catalogue verified: 516/516 files.\n");quit(status=0L)}
logs <- "logs/apaqtl_catalogue_v52"
dir.create(logs,recursive=TRUE,showWarnings=FALSE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type=="windows") "Rscript.exe" else "Rscript")
if (!file.exists(rscript)) stop("Rscript missing: ",rscript,call.=FALSE)
cat("Starting S24; valid existing chunks can be reused. Logs: ",logs,"\n",sep="")
code <- suppressWarnings(system2(rscript,shQuote(source_script),stdout=file.path(logs,"S24.stdout.log"),stderr=file.path(logs,"S24.stderr.log"),wait=TRUE))
if (as.integer(code)!=0L || !check()) {
  if (file.exists(status)) print(fread(status))
  stop("S24 incomplete or failed (exit ",code,"); inspect ",logs," and ",status," before rerunning",call.=FALSE)
}
cat("PASS S24: 516/516 files; harmonized catalogue: ",final,"\n",sep="")
