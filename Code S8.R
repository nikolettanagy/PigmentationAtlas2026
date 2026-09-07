#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
candidate_file <- if (length(args) >= 1) args[1] else "results/candidate_genes/candidate_gene_prioritization_v2.tsv"
out_root <- if (length(args) >= 2) args[2] else "results/expression"
thr_file <- if (length(args) >= 3) args[3] else "config/expression_thresholds.tsv"
db_dir <- file.path(out_root,"database"); qc_dir <- file.path(out_root,"QC"); log_dir <- file.path(out_root,"logs")
dir.create(qc_dir,recursive=TRUE,showWarnings=FALSE); dir.create(log_dir,recursive=TRUE,showWarnings=FALSE)
log_file <- file.path(log_dir,paste0("05B_",format(Sys.time(),"%Y%m%d_%H%M%S"),".log"))
logmsg <- function(...) {x<-paste0(...);cat(x,"\n");cat(x,"\n",file=log_file,append=TRUE)}
fail <- function(...) {x<-paste0(...);logmsg("ERROR: ",x);stop(x,call.=FALSE)}
strip_ver <- function(x) sub("\\.[0-9]+$","",as.character(x))
first_col <- function(nms,a,req=TRUE,label="column") {h<-a[a%in%nms];if(length(h))return(h[1]);if(req)fail("Missing ",label,". Tried: ",paste(a,collapse=", "));NA_character_}

idx_file <- file.path(db_dir,"expression_database_index.tsv")
if(!file.exists(candidate_file)) fail("Candidate file not found: ",candidate_file)
if(!file.exists(idx_file)) fail("Database index not found: ",idx_file)
if(!file.exists(thr_file)) fail("Threshold file not found: ",thr_file)

cand <- fread(candidate_file); idx <- fread(idx_file); thr <- fread(thr_file)
val <- function(name,default){x<-thr[parameter==name,value];if(length(x)&&!is.na(suppressWarnings(as.numeric(x[1]))))as.numeric(x[1])else default}
nominal_p <- val("nominal_p",0.05); strong_p <- val("strong_p",5e-8); q_cut <- val("q_value",0.05); min_strong <- as.integer(val("minimum_strong_variants",2))

gc <- first_col(names(cand),c("gene_id","ensembl_gene_id","gene_id_clean"),TRUE,"gene ID")
lc <- first_col(names(cand),c("locus_id","locus"),TRUE,"locus ID")
sc <- first_col(names(cand),c("gene_name","gene_symbol","symbol"),FALSE)
bc <- first_col(names(cand),c("gene_type","gene_biotype","biotype"),FALSE)
cc <- first_col(names(cand),c("candidate_source","source"),FALSE)
rc <- first_col(names(cand),c("candidate_rank","rank","structural_rank"),FALSE)

lookup <- unique(cand[,.(
  locus_id=as.character(get(lc)), gene_id=strip_ver(get(gc)),
  gene_symbol=if(!is.na(sc))as.character(get(sc))else NA_character_,
  gene_biotype=if(!is.na(bc))as.character(get(bc))else NA_character_,
  candidate_source=if(!is.na(cc))as.character(get(cc))else NA_character_,
  structural_rank=if(!is.na(rc))suppressWarnings(as.numeric(get(rc)))else NA_real_
)])

summary_list <- list(); variant_list <- list()
for(i in seq_len(nrow(idx))){
  m<-idx[i]; f<-as.character(m$standardized_file); id<-as.character(m$dataset_id)
  logmsg("[",i,"/",nrow(idx),"] ",id)
  if(!file.exists(f)){logmsg("WARNING missing: ",f);next}
  e <- fread(f)
  e[, gene_id := strip_ver(gene_id)]
  e[, rsid := as.character(rsid)]
  hit<-merge(lookup,e,by="gene_id",allow.cartesian=TRUE)
  if(nrow(hit)) variant_list[[length(variant_list)+1]]<-hit
  obs<-if(nrow(hit)) hit[,{
    j<-which.min(p_value); nnom<-sum(p_value<nominal_p,na.rm=TRUE); nstrong<-sum(p_value<strong_p,na.rm=TRUE); nq<-sum(!is.na(q_value)&q_value<=q_cut,na.rm=TRUE)
    .(eqtl_records=.N,unique_variants=uniqueN(variant_id),nominal_variants=nnom,strong_p_variants=nstrong,q_significant_variants=nq,
      lead_variant=variant_id[j],lead_rsid=rsid[j],lead_p=p_value[j],lead_beta=beta[j],lead_standard_error=standard_error[j],lead_q_value=q_value[j],
      eqtl_detected=nnom>0,strong_statistical_signal=(nstrong>=min_strong||nq>=1))
  },by=.(locus_id,gene_id,gene_symbol,gene_biotype,candidate_source,structural_rank,dataset_id,database,database_version,tissue,tissue_category,trait_relevance,genome_build)] else data.table()
  dataset_meta <- unique(
  m[, .(
    dataset_id,
    database,
    database_version,
    tissue,
    tissue_category,
    trait_relevance,
    genome_build
  )]
)

grid <- cbind(
  lookup[rep(seq_len(.N), each = nrow(dataset_meta))],
  dataset_meta[rep(seq_len(.N), times = nrow(lookup))]
)
  z<-merge(grid,obs,by=c("locus_id","gene_id","gene_symbol","gene_biotype","candidate_source","structural_rank","dataset_id","database","database_version","tissue","tissue_category","trait_relevance","genome_build"),all.x=TRUE)
  for(k in c("eqtl_records","unique_variants","nominal_variants","strong_p_variants","q_significant_variants")) z[is.na(get(k)),(k):=0L]
  z[is.na(eqtl_detected),eqtl_detected:=FALSE];z[is.na(strong_statistical_signal),strong_statistical_signal:=FALSE]
  z[,expression_dataset_status:=fifelse(eqtl_records>0,"records_available","no_records_for_gene")]
  z[,expression_support_class:=fcase(strong_statistical_signal,"Strong",eqtl_detected,"Nominal",eqtl_records>0,"Detected_but_not_nominal",default="None")]
  summary_list[[length(summary_list)+1]]<-z
  rm(e,hit,obs,grid,z);invisible(gc())
}

ds<-rbindlist(summary_list,fill=TRUE); if(!nrow(ds))fail("No summaries generated.")
variants<-if(length(variant_list))rbindlist(variant_list,fill=TRUE)else data.table()

gs<-ds[,{
  available<-.SD[eqtl_records>0]; best<-if(nrow(available))available[which.min(lead_p)]else .SD[1]
  nd<-sum(eqtl_records>0); nn<-sum(eqtl_detected); ns<-sum(strong_statistical_signal); ndb <- uniqueN(database[strong_statistical_signal == TRUE])
nt  <- uniqueN(tissue[strong_statistical_signal == TRUE])
  cls<-if(ndb>=2||nt>=2)"Replicated"else if(ns>=1)"Strong"else if(nn>=1)"Nominal"else if(nd>=1)"Detected_but_not_nominal"else "None"
  .(candidate_source=candidate_source[1],structural_rank=structural_rank[1],datasets_tested=.N,datasets_with_records=nd,datasets_with_nominal_eqtl=nn,datasets_with_strong_signal=ns,
    databases_with_strong_signal=ndb,tissues_with_strong_signal=nt,best_dataset=if(nrow(available))best$dataset_id[1]else NA_character_,best_database=if(nrow(available))best$database[1]else NA_character_,
    best_tissue=if(nrow(available))best$tissue[1]else NA_character_,best_trait_relevance=if(nrow(available))best$trait_relevance[1]else NA_character_,best_variant=if(nrow(available))best$lead_variant[1]else NA_character_,
    best_rsid = if (nrow(available)) as.character(best$lead_rsid[1]) else NA_character_,best_p=if(nrow(available))best$lead_p[1]else NA_real_,best_beta=if(nrow(available))best$lead_beta[1]else NA_real_,expression_support_class=cls)
},by=.(locus_id,gene_id,gene_symbol,gene_biotype)]
ord<-c(None=0L,Detected_but_not_nominal=1L,Nominal=2L,Strong=3L,Replicated=4L);gs[,expression_support_order:=unname(ord[expression_support_class])]
gs[, expression_rank_within_locus := frankv(
  .SD,
  cols = c("expression_support_order", "best_p"),
  order = c(-1L, 1L),
  ties.method = "min",
  na.last = "keep"
), by = locus_id,
.SDcols = c("expression_support_order", "best_p")]
setorder(gs,locus_id,-expression_support_order,best_p,structural_rank)

fwrite(ds,file.path(out_root,"expression_evidence.tsv"),sep="\t")
if(nrow(variants))fwrite(variants,file.path(out_root,"expression_evidence_variants.tsv.gz"),sep="\t")
fwrite(gs,file.path(out_root,"expression_summary.tsv"),sep="\t")
fwrite(gs[,.(locus_id,gene_id,gene_symbol,gene_biotype,expression_support_class,expression_support_order,expression_rank_within_locus,datasets_with_records,datasets_with_nominal_eqtl,datasets_with_strong_signal,best_database,best_tissue,best_trait_relevance,best_variant,best_rsid,best_p,best_beta)],file.path(out_root,"expression_support.tsv"),sep="\t")
qc<-rbindlist(list(
 data.table(metric="candidate_locus_gene_pairs",value=nrow(lookup)),
 data.table(metric="unique_candidate_genes",value=uniqueN(lookup$gene_id)),
 data.table(metric="datasets_tested",value=nrow(idx)),
 data.table(metric="genes_with_records",value=uniqueN(ds[eqtl_records>0,gene_id])),
 data.table(metric="genes_with_nominal_eqtl",value=uniqueN(ds[eqtl_detected==TRUE,gene_id])),
 data.table(metric="genes_with_strong_signal",value=uniqueN(ds[strong_statistical_signal==TRUE,gene_id])),
 data.table(metric="locus_gene_pairs_replicated",value=nrow(gs[expression_support_class=="Replicated"])),
 data.table(metric="locus_gene_pairs_strong",value=nrow(gs[expression_support_class=="Strong"])),
 data.table(metric="locus_gene_pairs_nominal",value=nrow(gs[expression_support_class=="Nominal"])),
 data.table(metric="locus_gene_pairs_none",value=nrow(gs[expression_support_class=="None"]))
),fill=TRUE)
fwrite(qc,file.path(qc_dir,"05B_expression_evidence_QC.tsv"),sep="\t")
logmsg("Module 05B completed successfully.")
