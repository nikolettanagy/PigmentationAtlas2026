#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a<-commandArgs(trailingOnly=TRUE);opt<-function(k){i<-match(k,a);if(is.na(i)||i==length(a))stop('Missing ',k,call.=FALSE);a[[i+1L]]}
clean<-normalizePath(opt('--workspace'),winslash='/',mustWork=TRUE);old<-normalizePath(opt('--old'),winslash='/',mustWork=TRUE)
out<-file.path(clean,'results','module12','reconciliation_v59');if(!dir.exists(out))stop('Run v59 first',call.=FALSE)
x<-fread(file.path(old,'module12_3_coloc_catalogue_summary.tsv.gz'),showProgress=FALSE)
y<-fread(file.path(clean,'results','module12','colocalization','catalogue','summary','module12_3_coloc_catalogue_summary.tsv.gz'),showProgress=FALSE)
keys<-c('tissue','phenotype_id','gene_id');pp<-paste0('PP.H',0:4)
if(!all(c(keys,'clump_id',pp)%in%names(x))||!all(c(keys,'clump_id',pp)%in%names(y)))stop('Missing comparison columns',call.=FALSE)
for(k in keys){set(x,j=k,value=as.character(x[[k]]));set(y,j=k,value=as.character(y[[k]]))}
x[,pair_count:=.N,by=keys];y[,pair_count:=.N,by=keys]
ux<-x[pair_count==1L,c(keys,'clump_id','coloc_unit_id',pp),with=FALSE]
uy<-y[pair_count==1L,c(keys,'clump_id','coloc_unit_id',pp),with=FALSE]
setnames(ux,c('clump_id','coloc_unit_id',pp),paste0('old_',c('clump_id','coloc_unit_id',pp)))
setnames(uy,c('clump_id','coloc_unit_id',pp),paste0('new_',c('clump_id','coloc_unit_id',pp)))
z<-merge(ux,uy,by=keys,all=FALSE,sort=FALSE)
for(p in pp)set(z,j=paste0('delta_',p),value=as.numeric(z[[paste0('new_',p)]])-as.numeric(z[[paste0('old_',p)]]))
z[,max_abs_delta:=do.call(pmax,c(lapply(.SD,abs),na.rm=FALSE)),.SDcols=paste0('delta_',pp)]
z[,clump_relabelled:=old_clump_id!=new_clump_id]
z[,posterior_changed:=is.na(max_abs_delta)|max_abs_delta>1e-6]
fwrite(z,file.path(out,'unique_biological_key_matches.tsv.gz'),sep='\t',na='NA',compress='gzip')
map<-z[,.(matched_units=.N,posterior_changed=sum(posterior_changed)),by=.(old_clump_id,new_clump_id)][order(-matched_units)]
fwrite(map,file.path(out,'clump_relabel_mapping.tsv'),sep='\t',na='NA')
metrics<-data.table(metric=c('old_rows','new_rows','old_ambiguous_without_clump','new_ambiguous_without_clump','unique_biological_key_matches','matched_same_clump_id','matched_relabelled_clump_id','relabelled_with_same_posterior','matched_posterior_changed','old_unmatched_unique_key','new_unmatched_unique_key'),value=c(nrow(x),nrow(y),x[pair_count>1L,.N],y[pair_count>1L,.N],nrow(z),z[clump_relabelled==FALSE,.N],z[clump_relabelled==TRUE,.N],z[clump_relabelled==TRUE&posterior_changed==FALSE,.N],z[posterior_changed==TRUE,.N],nrow(ux)-nrow(z),nrow(uy)-nrow(z)))
fwrite(metrics,file.path(out,'clump_label_audit.tsv'),sep='\t');cat('Clump-label audit:\n');print(metrics);cat('Largest relabel mappings:\n');print(head(map[old_clump_id!=new_clump_id],20L))
