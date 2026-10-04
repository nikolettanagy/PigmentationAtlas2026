#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a<-commandArgs(trailingOnly=TRUE);i<-match('--workspace',a);if(is.na(i)||i==length(a))stop('Usage: Rscript audit_release_figures_v70.R --workspace PATH',call.=FALSE)
root<-normalizePath(a[[i+1L]],winslash='/',mustWork=TRUE)
read<-function(...) {p<-file.path(root,...);if(!file.exists(p))stop('Missing: ',p,call.=FALSE);fread(p,showProgress=FALSE)}
prior<-read('results','module12','colocalization','prioritization','module12_4_status.tsv')
reclass<-read('results','module12','reclassification','module12_6_status.tsv')
classes<-read('results','module12','colocalization','prioritization','tables','module12_4_priority_class_counts.tsv')
stopifnot(nrow(prior)==1L,prior$final_status[[1L]]=='PASS',nrow(reclass)==1L,reclass$status[[1L]]=='SUCCESS')
expected<-as.integer(prior$input_coloc_units[[1L]])
rc<-as.integer(reclass$evaluable_clumps[[1L]]);changed<-as.integer(reclass$reclassified_clumps[[1L]])
assert<-function(label,observed,required) data.table(check=label,observed=as.character(observed),expected=as.character(required),status=if(length(observed)==1L&&!is.na(observed)&&isTRUE(all.equal(as.numeric(observed),as.numeric(required))))'PASS' else 'FAIL')
checks<-list(assert('clean catalogue units',expected,9814),assert('evaluable clumps',rc,130),assert('reclassified clumps',changed,111))
get_class<-function(label){x<-classes[priority_class==label,n_units];if(length(x)!=1L)stop('Missing class ',label,call.=FALSE);as.integer(x)}
labels<-c('STRONG_H4','MODERATE_H4','SUGGESTIVE_H4','H3_DOMINANT','INCONCLUSIVE')
for(j in seq_along(labels))checks[[length(checks)+1L]]<-assert(paste0('priority ',labels[[j]]),get_class(labels[[j]]),c(30,54,239,424,9067)[[j]])
f2<-read('figures','main','Figure_2_colocalization_category_counts.tsv')
s2<-read('results','Supplementary_Figure_S2_evidence_class_counts.tsv')
s3a<-read('results','module13','Supplementary_Figure_S3A_candidate_class_statistics.tsv')
s3b<-read('results','module13','Supplementary_Figure_S3B_h4_supported_statistics.tsv')
s4<-read('results','Supplementary_Figure_S4_tissue_summary.tsv')
s5<-read('results','figure_source_data','Supplementary_Figure_S5_extended_H3_H4_source_data.tsv.gz')
s6<-read('results','figure_source_data','Supplementary_Figure_S6_reclassification_source_data.tsv.gz')
checks<-c(checks,list(assert('Figure 2 units',sum(f2$n_units),expected),assert('S2 units',sum(s2$n),expected),assert('S3A units',sum(s3a$count),expected),assert('S3B H4 units',sum(s3b$count),sum(classes[priority_class%in%labels[1:3],n_units])),assert('S4 tissue units',sum(s4$analyses),expected),assert('S5 source rows',nrow(s5),expected),assert('S6 clumps',nrow(s6),as.integer(reclass$total_clumps[[1L]])),assert('S6 evaluable',sum(as.logical(s6$evaluable),na.rm=TRUE),rc),assert('S6 reclassified',sum(as.logical(s6$interpretation_change),na.rm=TRUE),changed)))
fig3<-read('figures','Figure3_Gene_Reclassification_status.tsv');fig4<-read('figures','figure_data','Figure4_status.tsv');fig5<-read('figures','figure_data','Figure5_status.tsv')
checks<-c(checks,list(assert('Figure 3 evaluable',fig3$evaluable_clumps[[1L]],rc),assert('Figure 3 reclassified',fig3$reclassified_clumps[[1L]],changed),assert('Figure 4 evaluable',fig4$evaluable_clumps[[1L]],rc),assert('Figure 4 reclassified',fig4$reclassified_clumps[[1L]],changed),assert('Figure 5 evaluable',fig5$evaluable_clumps[[1L]],rc),assert('Figure 5 reclassified',fig5$reclassified_clumps[[1L]],changed)))
map<-data.table(label=c(paste0('Figure',1:5),paste0('S',1:6)),base=c('figures/Figure1_PigmentationAtlas_Workflow','figures/main/Figure_2_PigmentationAtlas_colocalization_summary','figures/Figure3_Gene_Reclassification','figures/Figure4_Genomewide_Reclassification_Landscape','figures/Figure5_Conceptual_Biological_Summary','figures/Supplementary_Figure_S1_pipeline_audit','figures/Supplementary_Figure_S2_posterior_probabilities','figures/supplementary/Supplementary_Figure_S3_candidate_gene_biotype_landscape','figures/supplementary/Supplementary_Figure_S4_tissue_specific_colocalization','figures/supplementary/Supplementary_Figure_S5_extended_H3_H4_posterior_analyses','figures/supplementary/Supplementary_Figure_S6_genomewide_reclassification_analyses'))
manifest<-map[,.(extension=c('pdf','png','tiff')),by=.(label,base)]
manifest[,path:=file.path(root,paste0(base,'.',extension))]
manifest[,bytes:=file.info(path)$size]
manifest[,status:=ifelse(!is.na(bytes)&bytes>0,'PASS','FAIL')]
audit<-rbindlist(checks)
out<-file.path(root,'results','release_figure_audit_v70');dir.create(out,recursive=TRUE,showWarnings=FALSE)
fwrite(audit,file.path(out,'checks.tsv'),sep='\t');fwrite(manifest,file.path(out,'figure_manifest.tsv'),sep='\t')
cat('Cross-table checks:\n');print(audit)
cat('Figure files: ',sum(manifest$status=='PASS'),'/',nrow(manifest),'\n',sep='')
pass<-all(audit$status=='PASS')&&all(manifest$status=='PASS')
cat('Result: ',if(pass)'PASS' else 'FAIL','\nOutput: ',out,'\n',sep='')
if(!pass)quit(save='no',status=1L)
