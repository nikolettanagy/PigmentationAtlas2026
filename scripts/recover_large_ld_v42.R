#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE)
i <- match("--workspace", args)
if (is.na(i) || i==length(args)) stop("Usage: --workspace <clean workspace>",call.=FALSE)
root <- normalizePath(args[[i+1L]],winslash="/",mustWork=TRUE)
setwd(root)
input <- "results/clump_ld_manifest.tsv"
output <- "results/clump_ld_manifest_recovered_v42.tsv"
audit_path <- "results/clump_ld_recovery_audit_v42.tsv"
if (file.exists(output) || file.exists(audit_path)) stop("Recovery outputs already exist; refusing overwrite",call.=FALSE)
if (!file.exists(input)) stop("Missing full LD manifest",call.=FALSE)
m <- fread(input)
if (nrow(m)!=1325L || anyDuplicated(m$clump_id)) stop("Expected 1325 unique clumps",call.=FALSE)
failed <- m[ld_status!="ready"]
if (nrow(failed)!=2L || !setequal(failed$clump_id,c("CLUMP_0266","CLUMP_0267")) || !all(failed$ld_status=="ld_variant_limit_exceeded")) stop("Unexpected failures; recovery is not applicable",call.=FALSE)
normalize_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^chr","",x,ignore.case=TRUE)
  x <- sub("([_:])b3[78]$","",x,ignore.case=TRUE)
  toupper(gsub("_",":",x,fixed=TRUE))
}
prepared <- vector("list",nrow(failed))
audit <- vector("list",nrow(failed))
for (j in seq_len(nrow(failed))) {
  row <- failed[j]
  matrix_path <- row$matrix_file
  vars_path <- row$vars_file
  map_path <- row$variant_map_file
  if (anyNA(c(matrix_path,vars_path,map_path)) || !all(file.exists(c(matrix_path,vars_path)))) stop("Missing LD outputs: ",row$clump_id,call.=FALSE)
  if (file.exists(map_path)) stop("Existing map needs manual review: ",map_path,call.=FALSE)
  ids <- readLines(vars_path,warn=FALSE)
  n <- length(ids)
  if (n!=as.integer(row$n_ld_variants) || n<=10000L || n>20000L || anyDuplicated(ids)) stop("Invalid .vars count/IDs: ",row$clump_id,call.=FALSE)
  bytes <- file.info(matrix_path)$size
  if (is.na(bytes) || bytes!=4*as.double(n)^2) stop("Binary matrix byte count mismatch: ",row$clump_id,call.=FALSE)
  region <- file.path("results","clumps",row$clump_id,"gwas","gwas_region.tsv.gz")
  if (!file.exists(region)) stop("Missing GWAS region: ",region,call.=FALSE)
  g <- fread(region)
  required <- c("rsid","variant_id","base_pair_location","effect_allele","other_allele")
  if (!all(required %in% names(g))) stop("Missing GWAS columns: ",row$clump_id,call.=FALSE)
  ids_norm <- normalize_id(ids)
  g_norm <- normalize_id(g$variant_id)
  if (anyDuplicated(ids_norm) || anyDuplicated(g_norm)) stop("Duplicate normalized IDs: ",row$clump_id,call.=FALSE)
  ix <- match(ids_norm,g_norm)
  if (anyNA(ix)) stop("Unmatched PLINK IDs: ",row$clump_id,"; count=",sum(is.na(ix)),call.=FALSE)
  selected <- g[ix]
  map <- data.table(rsid=as.character(selected$rsid),variant_id=as.character(selected$variant_id),position=as.integer(selected$base_pair_location),effect_allele=as.character(selected$effect_allele),reference_variant_id=ids,other_allele=as.character(selected$other_allele),ld_order=seq_len(n))
  if (!identical(normalize_id(map$variant_id),ids_norm) || anyNA(map$position)) stop("Map order/position invalid: ",row$clump_id,call.=FALSE)
  prepared[[j]] <- list(path=map_path, map=map)
  audit[[j]] <- data.table(clump_id=row$clump_id,n_ld_variants=n,matrix_bytes=bytes,gwas_rows=nrow(g),matched_ids=nrow(map),source_manifest=input)
  cat("VERIFIED",row$clump_id,"variants=",n,"matrix_bytes=",bytes,"\n")
}
# All checks above complete before any result is written.
for (j in seq_along(prepared)) {
  p <- prepared[[j]]$path
  dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE)
  tmp <- paste0(p,".v42.tmp")
  if (file.exists(tmp)) stop("Temporary output already exists: ",tmp,call.=FALSE)
  fwrite(prepared[[j]]$map,tmp,sep="\t",quote=FALSE,na="NA")
  if (!file.rename(tmp,p)) stop("Cannot install map: ",p,call.=FALSE)
}
for (id in failed$clump_id) {
  k <- match(id,m$clump_id)
  m[k,`:=`(ld_status="ready",ld_message="v42 recovered map from existing PLINK LD matrix",variant_map_exists=TRUE,matrix_exists=TRUE,vars_exists=TRUE)]
}
if (nrow(m)!=1325L || any(m$ld_status!="ready") || !all(file.exists(m$variant_map_file))) stop("Recovered manifest audit failed",call.=FALSE)
fwrite(rbindlist(audit),audit_path,sep="\t")
fwrite(m,output,sep="\t",quote=FALSE,na="NA")
cat("PASS recovered LD manifest: 1325 ready; maps for CLUMP_0266 and CLUMP_0267; original manifest preserved.\n")
