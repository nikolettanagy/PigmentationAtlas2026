#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly=TRUE)
i <- match("--workspace",a)
if (is.na(i) || i==length(a)) stop("Usage: --workspace <clean workspace>",call.=FALSE)
root <- normalizePath(a[[i+1L]],winslash="/",mustWork=TRUE)
setwd(root)
Sys.setenv(PIGMENTATIONATLAS_ROOT=root)
script <- "scripts/upstream/Code S18_v48.R"
pub <- "results/susie_manifest_publication_v47.tsv"
out <- "results/module11/module11_locus_manifest.tsv"
if (!file.exists(script) || !file.exists(pub)) stop("Missing S18 script or v47 publication manifest",call.=FALSE)
if (file.exists(out)) stop("Module 11 locus output exists; refusing overwrite",call.=FALSE)
expected <- fread(pub)
if (nrow(expected)!=265L || anyDuplicated(expected$clump_id)) stop("Unexpected publication manifest",call.=FALSE)
rscript <- file.path(R.home("bin"),"Rscript.exe")
logdir <- "logs/module11_audit_v48"
dir.create(logdir,recursive=TRUE,showWarnings=FALSE)
code <- suppressWarnings(system2(rscript,shQuote(script),stdout=file.path(logdir,"stdout.log"),stderr=file.path(logdir,"stderr.log"),wait=TRUE))
if (as.integer(code)!=0L || !file.exists(out)) stop("Module 11.0 failed; see ",logdir,call.=FALSE)
loci <- fread(out)
if (nrow(loci)!=265L || anyDuplicated(loci$locus_id) || !setequal(loci$locus_id,expected$clump_id)) stop("Locus output does not match publication clumps",call.=FALSE)
cat("PASS Module 11.0: 265 publication loci audited.\n")
