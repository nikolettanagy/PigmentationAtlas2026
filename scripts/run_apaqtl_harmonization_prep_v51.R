#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly=TRUE)
opt <- function(flag, default=NULL) {i<-match(flag,a); if(is.na(i)) return(default); if(i==length(a)) stop("Missing ",flag,call.=FALSE); a[[i+1L]]}
root <- normalizePath(opt("--workspace"),winslash="/",mustWork=TRUE)
start <- opt("--start-at","S21")
stages <- c("S21","S22","S23")
if (!start %in% stages) stop("--start-at must be S21, S22 or S23",call.=FALSE)
setwd(root)
Sys.setenv(PIGMENTATIONATLAS_ROOT=root)
base <- "results/module11/harmonization"
paths <- c(S21=file.path(base,"region_mapping_v4","module11_3A_v4_status.tsv"),S22=file.path(base,"catalogue","module11_3B_initialization_status.tsv"),S23=file.path(base,"catalogue","module11_3B_2_test_status.tsv"))
check <- function(stage) {
  p <- paths[[stage]]
  if (!file.exists(p)) return(FALSE)
  x <- fread(p)
  if (nrow(x)!=1L || !"final_status" %in% names(x) || x$final_status[1L]!="PASS") return(FALSE)
  if (stage=="S21") {q<-file.path(base,"region_mapping_v4","module11_3A_v4_direct_clump_region_map.tsv");return(file.exists(q) && nrow(fread(q))==530L)}
  if (stage=="S22") {q<-file.path(base,"catalogue","module11_3B_work_queue.tsv");return(file.exists(q) && nrow(fread(q))==516L)}
  TRUE
}
if (!file.exists("results/module11/eqtl_regions/summary/module11_2_apaqtl_region_manifest.tsv")) stop("Module 11.2 output missing",call.=FALSE)
logdir <- "logs/apaqtl_prep_v51"
dir.create(logdir,recursive=TRUE,showWarnings=FALSE)
rscript <- file.path(R.home("bin"),"Rscript.exe")
first <- match(start,stages)
if (first>1L) for (stage in stages[seq_len(first-1L)]) {
  if (!check(stage)) stop("Cannot resume; previous stage invalid: ",stage,call.=FALSE)
  cat("VERIFIED",stage,"\n")
}
for (stage in stages[first:length(stages)]) {
  if (file.exists(paths[[stage]])) stop("Existing status; review before rerun: ",paths[[stage]],call.=FALSE)
  script <- file.path("scripts","upstream",paste0("Code ",stage,".R"))
  if (!file.exists(script)) stop("Missing script: ",script,call.=FALSE)
  code <- suppressWarnings(system2(rscript,shQuote(script),stdout=file.path(logdir,paste0(stage,".stdout.log")),stderr=file.path(logdir,paste0(stage,".stderr.log")),wait=TRUE))
  if (as.integer(code)!=0L || !check(stage)) stop(stage," failed or audit mismatch; see ",logdir,call.=FALSE)
  cat("PASS",stage,paths[[stage]],"\n")
}
cat("PASS apaQTL harmonization preparation: 530 mapped pairs; 516-file work queue; single-file QC.\n")
