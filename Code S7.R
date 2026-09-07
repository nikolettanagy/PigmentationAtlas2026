#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

args <- commandArgs(trailingOnly = TRUE)
candidate_file <- if (length(args) >= 1) args[1] else "results/candidate_genes/candidate_gene_prioritization_v2.tsv"
manifest_file  <- if (length(args) >= 2) args[2] else "config/expression_datasets.tsv"
out_root       <- if (length(args) >= 3) args[3] else "results/expression"

db_dir <- file.path(out_root, "database")
qc_dir <- file.path(out_root, "QC")
log_dir <- file.path(out_root, "logs")
dir.create(db_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(qc_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(log_dir, recursive=TRUE, showWarnings=FALSE)
log_file <- file.path(log_dir, paste0("05A_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
logmsg <- function(...) {x <- paste0(...); cat(x, "\n"); cat(x, "\n", file=log_file, append=TRUE)}
fail <- function(...) {x <- paste0(...); logmsg("ERROR: ", x); stop(x, call.=FALSE)}
strip_ver <- function(x) sub("\\.[0-9]+$", "", as.character(x))
first_col <- function(nms, aliases, required=TRUE, label="column") {
  hit <- aliases[aliases %in% nms]
  if (length(hit)) return(hit[1])
  if (required) fail("Missing ", label, ". Tried: ", paste(aliases, collapse=", "))
  NA_character_
}
num <- function(x) suppressWarnings(as.numeric(x))

if (!file.exists(candidate_file)) fail("Candidate file not found: ", candidate_file)
if (!file.exists(manifest_file)) fail("Manifest not found: ", manifest_file)

cand <- fread(candidate_file)
man <- fread(manifest_file)

gene_col <- first_col(names(cand), c("gene_id","ensembl_gene_id","gene_id_clean"), TRUE, "candidate gene ID")
locus_col <- first_col(names(cand), c("locus_id","locus"), TRUE, "locus ID")
symbol_col <- first_col(names(cand), c("gene_name","gene_symbol","symbol"), FALSE)
biotype_col <- first_col(names(cand), c("gene_type","gene_biotype","biotype"), FALSE)

cand_lookup <- unique(cand[, .(
  locus_id=as.character(get(locus_col)),
  gene_id=strip_ver(get(gene_col)),
  gene_symbol=if (!is.na(symbol_col)) as.character(get(symbol_col)) else NA_character_,
  gene_biotype=if (!is.na(biotype_col)) as.character(get(biotype_col)) else NA_character_
)])[!is.na(gene_id) & gene_id != ""]

req <- c("dataset_id","database","tissue","tissue_category","trait_relevance","file_path","enabled")
miss <- setdiff(req, names(man)); if (length(miss)) fail("Manifest missing: ", paste(miss, collapse=", "))
man[, enabled := tolower(as.character(enabled)) %in% c("true","1","yes","y")]
man <- man[enabled == TRUE]
if (!nrow(man)) fail("No enabled datasets in manifest.")

logmsg("Candidate locus-gene pairs: ", nrow(cand_lookup))
logmsg("Unique candidate genes: ", uniqueN(cand_lookup$gene_id))
logmsg("Enabled datasets: ", nrow(man))

idx_list <- list(); qc_list <- list()
for (i in seq_len(nrow(man))) {
  m <- man[i]; id <- as.character(m$dataset_id); src <- as.character(m$file_path)
  logmsg("[",i,"/",nrow(man),"] ", id, " -> ", src)
  if (!file.exists(src)) {
    logmsg("WARNING: file not found; skipped")
    qc_list[[length(qc_list)+1]] <- data.table(dataset_id=id,status="file_not_found",source_file=src)
    next
  }
    file_format <- if ("file_format" %in% names(m)) {
    tolower(as.character(m$file_format))
  } else {
    tolower(tools::file_ext(src))
  }

  if (file_format %in% c("parquet", "pq")) {
    ds <- open_dataset(src, format = "parquet")
    hdr <- ds$schema$names
  } else {
    hdr <- names(fread(src, nrows = 0))
  }
  choose_manifest <- function(field, aliases, required=TRUE) {
    if (field %in% names(m) && !is.na(m[[field]]) && nzchar(as.character(m[[field]])) && as.character(m[[field]]) %in% hdr) as.character(m[[field]])
    else first_col(hdr, aliases, required, paste0(id," ",field))
  }
  gc <- choose_manifest("gene_id_column", c("gene_id","phenotype_id","molecular_trait_id","gene"))
  vc <- choose_manifest("variant_id_column", c("variant_id","variant","snp_id","SNP","rsid"))
  pc <- choose_manifest("p_column", c("pval_nominal","p_value","pvalue","p","pval"))
  rc <- first_col(hdr, c("rsid","rs_id","rs_id_dbSNP151_GRCh38p7","snp"), FALSE)
  bc <- first_col(hdr, c("slope","beta","effect_size","nes","NES"), FALSE)
  sc <- first_col(hdr, c("slope_se","se","standard_error","beta_se"), FALSE)
  qc <- first_col(hdr, c("qval","q_value","fdr","FDR"), FALSE)
  mc <- first_col(hdr, c("maf","MAF","minor_allele_frequency"), FALSE)
  ac <- first_col(hdr, c("eaf","af","allele_frequency","effect_allele_frequency"), FALSE)
  sel <- unique(na.omit(c(gc,vc,pc,rc,bc,sc,qc,mc,ac)))
   if (file_format %in% c("parquet", "pq")) {
    x <- as.data.table(
      read_parquet(
        src,
        col_select = sel,
        as_data_frame = TRUE
      )
    )
  } else {
    x <- fread(src, select = sel, showProgress = TRUE)
  }
  n_source <- nrow(x)
  x[, gene_id := strip_ver(get(gc))]
  x <- x[gene_id %chin% cand_lookup$gene_id]
  if (!nrow(x)) {
    qc_list[[length(qc_list)+1]] <- data.table(dataset_id=id,status="no_candidate_records",source_file=src,source_rows=n_source,candidate_rows=0)
    next
  }
  z <- data.table(
    dataset_id=id,
    database=as.character(m$database),
    database_version=if ("database_version" %in% names(m)) as.character(m$database_version) else NA_character_,
    tissue=as.character(m$tissue),
    tissue_category=as.character(m$tissue_category),
    trait_relevance=as.character(m$trait_relevance),
    genome_build=if ("genome_build" %in% names(m)) as.character(m$genome_build) else "GRCh38",
    gene_id=x$gene_id,
    variant_id=as.character(x[[vc]]),
    rsid=if (!is.na(rc)) as.character(x[[rc]]) else NA_character_,
    p_value=num(x[[pc]]),
    beta=if (!is.na(bc)) num(x[[bc]]) else NA_real_,
    standard_error=if (!is.na(sc)) num(x[[sc]]) else NA_real_,
    q_value=if (!is.na(qc)) num(x[[qc]]) else NA_real_,
    maf=if (!is.na(mc)) num(x[[mc]]) else NA_real_,
    allele_frequency=if (!is.na(ac)) num(x[[ac]]) else NA_real_
  )[!is.na(p_value)]
  z[, minus_log10_p := -log10(pmax(p_value, .Machine$double.xmin))]
  setorder(z, gene_id, p_value)
  out <- file.path(db_dir, paste0(id, ".candidate_eqtl.tsv.gz"))
  fwrite(z, out, sep="\t")
  idx_list[[length(idx_list)+1]] <- data.table(
    dataset_id=id,database=as.character(m$database),
    database_version=if ("database_version" %in% names(m)) as.character(m$database_version) else NA_character_,
    tissue=as.character(m$tissue),tissue_category=as.character(m$tissue_category),trait_relevance=as.character(m$trait_relevance),
    genome_build=if ("genome_build" %in% names(m)) as.character(m$genome_build) else "GRCh38",
    source_file=src,standardized_file=out,n_records=nrow(z),n_genes=uniqueN(z$gene_id)
  )
  qc_list[[length(qc_list)+1]] <- data.table(dataset_id=id,status="ok",source_file=src,source_rows=n_source,candidate_rows=nrow(z),candidate_genes_found=uniqueN(z$gene_id),output_file=out)
  logmsg("Retained rows: ",nrow(z),"; genes: ",uniqueN(z$gene_id))
  rm(x,z); invisible(gc())
}

idx <- if (length(idx_list)) rbindlist(idx_list, fill=TRUE) else data.table()
qct <- rbindlist(qc_list, fill=TRUE)
fwrite(idx, file.path(db_dir,"expression_database_index.tsv"), sep="\t")
fwrite(cand_lookup, file.path(db_dir,"candidate_gene_lookup.tsv"), sep="\t")
fwrite(qct, file.path(qc_dir,"05A_expression_database_QC.tsv"), sep="\t")
if (!nrow(idx)) fail("No expression dataset was successfully built.")
logmsg("Module 05A completed successfully.")
