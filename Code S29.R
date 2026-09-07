#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
options(stringsAsFactors=FALSE, scipen=999, warn=1)
start_time <- Sys.time()
project_dir <- normalizePath(getwd(), winslash='/', mustWork=TRUE)
base_dir <- file.path(project_dir,'results','module12','colocalization')
in_summary <- file.path(base_dir,'catalogue','summary','module12_3_coloc_catalogue_summary.tsv.gz')
in_qc <- file.path(base_dir,'catalogue','summary','module12_3_coloc_catalogue_qc.tsv.gz')
out_dir <- file.path(base_dir,'prioritization')
tables_dir <- file.path(out_dir,'tables'); qc_dir <- file.path(out_dir,'qc'); logs_dir <- file.path(out_dir,'logs')
dir.create(tables_dir,recursive=TRUE,showWarnings=FALSE); dir.create(qc_dir,recursive=TRUE,showWarnings=FALSE); dir.create(logs_dir,recursive=TRUE,showWarnings=FALSE)
log_file <- file.path(logs_dir,paste0('module12_4_',format(Sys.time(),'%Y%m%d_%H%M%S'),'.log'))
log_msg <- function(...) {x<-paste0(...,collapse=''); y<-paste0('[',format(Sys.time(),'%Y-%m-%d %H:%M:%S'),'] ',x); cat(y,'\n'); cat(y,'\n',file=log_file,append=TRUE)}
stop_mod <- function(x){log_msg('FATAL: ',x); stop(x,call.=FALSE)}
write_tab <- function(x,n,gz=FALSE){p<-file.path(tables_dir,n); fwrite(x,p,sep='\t',quote=FALSE,na='NA',compress=if(gz)'gzip' else 'none'); p}
snum <- function(x) suppressWarnings(as.numeric(as.character(x)))
smax <- function(x){x<-snum(x); if(all(is.na(x))) NA_real_ else max(x,na.rm=TRUE)}
smed <- function(x){x<-snum(x); if(all(is.na(x))) NA_real_ else median(x,na.rm=TRUE)}
cat('============================================================\nPigmentationAtlas\nMODULE 12.4 – COLOCALIZATION PRIORITIZATION\n============================================================\n\n')
if(!file.exists(in_summary)) stop_mod(paste0('Missing input: ',in_summary))
dt <- fread(in_summary,showProgress=FALSE)
req <- c('coloc_unit_id','clump_id','tissue','phenotype_id','gene_id','nsnps','PP.H0','PP.H1','PP.H2','PP.H3','PP.H4','PP.H4_over_H3_H4','posterior_sum','top_coloc_variant','top_coloc_variant_PP.H4')
miss <- setdiff(req,names(dt)); if(length(miss)) stop_mod(paste('Missing columns:',paste(miss,collapse=', ')))
if(anyDuplicated(dt$coloc_unit_id)){log_msg('WARNING: duplicate coloc_unit_id rows removed.'); dt<-unique(dt,by='coloc_unit_id')}
for(n in intersect(c('nsnps','PP.H0','PP.H1','PP.H2','PP.H3','PP.H4','PP.H4_over_H3_H4','posterior_sum','top_coloc_variant_PP.H4'),names(dt))) set(dt,j=n,value=snum(dt[[n]]))
strong <- 0.80; moderate <- 0.50; suggestive <- 0.20; h3dom <- 0.80; cond_thr <- 0.80; min_snps <- 50L
dt[, posterior_sum_deviation:=abs(posterior_sum-1)]
dt[, posterior_sum_valid:=!is.na(posterior_sum_deviation)&posterior_sum_deviation<=1e-6]
dt[, priority_class:=fcase(PP.H4>=strong,'STRONG_H4',PP.H4>=moderate,'MODERATE_H4',PP.H4>=suggestive,'SUGGESTIVE_H4',PP.H3>=h3dom,'H3_DOMINANT',default='INCONCLUSIVE')]
dt[, priority_rank:=fcase(priority_class=='STRONG_H4',1L,priority_class=='MODERATE_H4',2L,priority_class=='SUGGESTIVE_H4',3L,priority_class=='H3_DOMINANT',4L,default=5L)]
dt[, conditional_h4_support:=!is.na(PP.H4_over_H3_H4)&PP.H4_over_H3_H4>=cond_thr]
dt[, adequate_variant_count:=!is.na(nsnps)&nsnps>=min_snps]
dt[, high_confidence_h4:=PP.H4>=strong&conditional_h4_support&adequate_variant_count&posterior_sum_valid]
dt[, unconditional_top_variant_shared_posterior:=PP.H4*top_coloc_variant_PP.H4]
dt[, h4_minus_h3:=PP.H4-PP.H3]
dt[, evidence_label:=fcase(high_confidence_h4,'HIGH_CONFIDENCE_SHARED_SIGNAL',priority_class=='STRONG_H4','STRONG_SHARED_SIGNAL',priority_class=='MODERATE_H4'&conditional_h4_support,'MODERATE_SHARED_SIGNAL_WITH_CONDITIONAL_SUPPORT',priority_class=='MODERATE_H4','MODERATE_SHARED_SIGNAL',priority_class=='SUGGESTIVE_H4'&conditional_h4_support,'SUGGESTIVE_SHARED_SIGNAL_WITH_CONDITIONAL_SUPPORT',priority_class=='SUGGESTIVE_H4','SUGGESTIVE_SHARED_SIGNAL',priority_class=='H3_DOMINANT','DISTINCT_CAUSAL_SIGNALS',default='INCONCLUSIVE')]
setorderv(dt,c('priority_rank','PP.H4','PP.H4_over_H3_H4','top_coloc_variant_PP.H4','nsnps','coloc_unit_id'),c(1,-1,-1,-1,-1,1),na.last=TRUE)
dt[,global_rank:=.I]
ranked_file <- write_tab(dt,'module12_4_ranked_coloc_catalogue.tsv.gz',TRUE)
write_tab(dt[priority_class=='STRONG_H4'],'module12_4_strong_H4_hits.tsv')
write_tab(dt[priority_class=='MODERATE_H4'],'module12_4_moderate_H4_hits.tsv')
write_tab(dt[priority_class=='SUGGESTIVE_H4'],'module12_4_suggestive_H4_hits.tsv')
write_tab(dt[priority_class=='H3_DOMINANT'],'module12_4_H3_dominant_hits.tsv')
write_tab(dt[high_confidence_h4 == TRUE],'module12_4_high_confidence_H4_hits.tsv')
write_tab(dt[conditional_h4_support == TRUE & PP.H4 < moderate],'module12_4_high_conditional_ratio_low_absolute_H4.tsv')
top_by <- function(x,grp){setorderv(x,c(grp,'priority_rank','PP.H4','PP.H4_over_H3_H4','nsnps'),c(rep(1,length(grp)),1,-1,-1,-1),na.last=TRUE); x[, .SD[1L], by=grp]}
write_tab(top_by(copy(dt),'clump_id'),'module12_4_top_coloc_per_clump.tsv')
write_tab(top_by(copy(dt),'gene_id'),'module12_4_top_coloc_per_gene.tsv')
write_tab(top_by(copy(dt),'tissue'),'module12_4_top_coloc_per_tissue.tsv')
write_tab(top_by(copy(dt),'phenotype_id'),'module12_4_top_coloc_per_phenotype.tsv')
write_tab(top_by(copy(dt),c('clump_id','tissue')),'module12_4_top_coloc_per_clump_tissue.tsv')
summary_by <- function(x,grp) x[,.(n_coloc_units=.N,n_strong_h4=sum(priority_class=='STRONG_H4'),n_moderate_h4=sum(priority_class=='MODERATE_H4'),n_suggestive_h4=sum(priority_class=='SUGGESTIVE_H4'),n_h3_dominant=sum(priority_class=='H3_DOMINANT'),n_inconclusive=sum(priority_class=='INCONCLUSIVE'),n_high_confidence_h4=sum(high_confidence_h4),max_PP.H4=smax(PP.H4),median_PP.H4=smed(PP.H4),max_conditional_H4=smax(PP.H4_over_H3_H4),median_nsnps=smed(nsnps)),by=grp]
gene_sum<-summary_by(dt,'gene_id'); setorder(gene_sum,-n_high_confidence_h4,-n_strong_h4,-n_moderate_h4,-max_PP.H4,gene_id); write_tab(gene_sum,'module12_4_gene_summary.tsv')
tissue_sum<-summary_by(dt,'tissue'); setorder(tissue_sum,-n_high_confidence_h4,-n_strong_h4,-n_moderate_h4,-max_PP.H4,tissue); write_tab(tissue_sum,'module12_4_tissue_summary.tsv')
clump_sum<-summary_by(dt,'clump_id'); setorder(clump_sum,-n_high_confidence_h4,-n_strong_h4,-n_moderate_h4,-max_PP.H4,clump_id); write_tab(clump_sum,'module12_4_clump_summary.tsv')
phen_sum<-summary_by(dt,c('phenotype_id','gene_id')); setorder(phen_sum,-n_high_confidence_h4,-n_strong_h4,-n_moderate_h4,-max_PP.H4,phenotype_id); write_tab(phen_sum,'module12_4_phenotype_summary.tsv')
counts<-dt[,.(n_units=.N,percent_units=100*.N/nrow(dt),n_unique_clumps=uniqueN(clump_id),n_unique_genes=uniqueN(gene_id),n_unique_phenotypes=uniqueN(phenotype_id),n_unique_tissues=uniqueN(tissue),median_PP.H4=smed(PP.H4),max_PP.H4=smax(PP.H4),median_nsnps=smed(nsnps)),by=.(priority_rank,priority_class)][order(priority_rank)]
write_tab(counts,'module12_4_priority_class_counts.tsv')
overview<-data.table(metric=c('Total eligible colocalization units','Unique PLINK clumps','Unique genes','Unique APA phenotypes','Unique tissues','Strong H4 units (PP.H4 >= 0.80)','Moderate H4 units (0.50 <= PP.H4 < 0.80)','Suggestive H4 units (0.20 <= PP.H4 < 0.50)','H3-dominant units (PP.H3 >= 0.80)','Inconclusive units','High-confidence H4 units','Units with conditional H4 support >= 0.80','Median variants per coloc unit','Minimum variants per coloc unit','Maximum variants per coloc unit','Posterior-sum QC failures'),value=as.character(c(nrow(dt),uniqueN(dt$clump_id),uniqueN(dt$gene_id),uniqueN(dt$phenotype_id),uniqueN(dt$tissue),dt[priority_class=='STRONG_H4',.N],dt[priority_class=='MODERATE_H4',.N],dt[priority_class=='SUGGESTIVE_H4',.N],dt[priority_class=='H3_DOMINANT',.N],dt[priority_class=='INCONCLUSIVE',.N],dt[high_confidence_h4 == TRUE, .N],dt[conditional_h4_support == TRUE, .N],smed(dt$nsnps),min(dt$nsnps,na.rm=TRUE),max(dt$nsnps,na.rm=TRUE),dt[posterior_sum_valid == FALSE, .N])))
write_tab(overview,'module12_4_manuscript_overview.tsv')
qc<-data.table(metric=c('catalogue_rows','unique_coloc_unit_ids','duplicate_coloc_unit_ids','missing_PP.H4','missing_PP.H3','missing_conditional_H4','missing_nsnps','nsnps_below_high_confidence_minimum','invalid_posterior_sum','missing_top_coloc_variant','missing_top_variant_posterior'),value=c(nrow(dt),uniqueN(dt$coloc_unit_id),nrow(dt)-uniqueN(dt$coloc_unit_id),sum(is.na(dt$PP.H4)),sum(is.na(dt$PP.H3)),sum(is.na(dt$PP.H4_over_H3_H4)),sum(is.na(dt$nsnps)),sum(!is.na(dt$nsnps)&dt$nsnps<min_snps),sum(!dt$posterior_sum_valid,na.rm=TRUE),sum(is.na(dt$top_coloc_variant)|!nzchar(dt$top_coloc_variant)),sum(is.na(dt$top_coloc_variant_PP.H4))))
fwrite(qc,file.path(qc_dir,'module12_4_qc_summary.tsv'),sep='\t',quote=FALSE,na='NA')
critical<-qc[metric%in%c('duplicate_coloc_unit_ids','missing_PP.H4','missing_PP.H3','invalid_posterior_sum'),sum(value)]
warns<-qc[metric%in%c('missing_conditional_H4','missing_nsnps','nsnps_below_high_confidence_minimum','missing_top_coloc_variant','missing_top_variant_posterior'),sum(value)]
final_status<-if(critical>0)'FAIL' else if(warns>0)'PASS_WITH_WARNINGS' else 'PASS'
status<-data.table(module='12.4',start_time=format(start_time,'%Y-%m-%d %H:%M:%S'),end_time=format(Sys.time(),'%Y-%m-%d %H:%M:%S'),runtime_seconds=as.numeric(difftime(Sys.time(),start_time,units='secs')),input_coloc_units=nrow(dt),unique_clumps=uniqueN(dt$clump_id),unique_genes=uniqueN(dt$gene_id),unique_phenotypes=uniqueN(dt$phenotype_id),unique_tissues=uniqueN(dt$tissue),strong_h4_units=dt[priority_class=='STRONG_H4',.N],moderate_h4_units=dt[priority_class=='MODERATE_H4',.N],suggestive_h4_units=dt[priority_class=='SUGGESTIVE_H4',.N],h3_dominant_units=dt[priority_class=='H3_DOMINANT',.N],high_confidence_h4_units=dt[high_confidence_h4 == TRUE, .N],critical_failures=critical,warnings=warns,final_status=final_status)
fwrite(status,file.path(out_dir,'module12_4_status.tsv'),sep='\t',quote=FALSE,na='NA')
cat('\n============================================================\nMODULE 12.4 – FINAL SUMMARY\n============================================================\n\n'); print(status); cat('\nPriority-class counts:\n'); print(counts); cat('\nRanked catalogue:\n  ',ranked_file,'\n',sep=''); cat('\n============================================================\nFINAL RESULT: ',final_status,'\n============================================================\n',sep='')
quit(save='no',status=if(final_status=='FAIL')1L else 0L)
