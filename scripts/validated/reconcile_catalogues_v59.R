#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a<-commandArgs(trailingOnly=TRUE)
opt<-function(k){i<-match(k,a);if(is.na(i)||i==length(a))stop('Missing ',k,call.=FALSE);a[[i+1L]]}
clean<-normalizePath(opt('--workspace'),winslash='/',mustWork=TRUE)
old<-normalizePath(opt('--old'),winslash='/',mustWork=TRUE)
oldroot<-normalizePath(opt('--old-project'),winslash='/',mustWork=TRUE)
out<-file.path(clean,'results','module12','reconciliation_v59')
if(dir.exists(out))stop('Output already exists: ',out,call.=FALSE)
dir.create(out,recursive=TRUE)
read<-function(p){if(!file.exists(p))stop('Missing: ',p,call.=FALSE);fread(p,showProgress=FALSE)}
oldfile<-file.path(old,'module12_3_coloc_catalogue_summary.tsv.gz')
newfile<-file.path(clean,'results','module12','colocalization','catalogue','summary','module12_3_coloc_catalogue_summary.tsv.gz')
oldqueue<-file.path(oldroot,'results','module12_rebuild_staging','colocalization','work_queue','module12_1_coloc_work_queue_eligible.tsv.gz')
newqueue<-file.path(clean,'results','module12','colocalization','work_queue','module12_1_coloc_work_queue_eligible.tsv.gz')
oldx<-read(oldfile);newx<-read(newfile);oq<-read(oldqueue);nq<-read(newqueue)
keys<-c('clump_id','tissue','phenotype_id','gene_id')
for(n in c('old catalogue','new catalogue','old queue','new queue')){
 x<-switch(n,'old catalogue'=oldx,'new catalogue'=newx,'old queue'=oq,'new queue'=nq)
 miss<-setdiff(c(keys,'coloc_unit_id'),names(x));if(length(miss))stop(n,' missing: ',paste(miss,collapse=', '),call.=FALSE)
 for(k in keys)set(x,j=k,value=as.character(x[[k]]))
 if(anyNA(x[,..keys])||anyDuplicated(x,by=keys))stop(n,' has missing or duplicate natural keys',call.=FALSE)
}
if(nrow(oldx)!=9471L||nrow(oq)!=9471L||nrow(newx)!=9814L||nrow(nq)!=9814L)stop('Unexpected row count',call.=FALSE)
keystring<-function(x)do.call(paste,c(x[,..keys],sep='|'))
if(!setequal(keystring(oldx),keystring(oq))||!setequal(keystring(newx),keystring(nq)))stop('Catalogue/queue key mismatch',call.=FALSE)
pp<-paste0('PP.H',0:4)
fields<-c('coloc_unit_id',intersect(c(pp,'PP.H4_over_H3_H4','posterior_sum','nsnps','interpretation'),intersect(names(oldx),names(newx))))
oldsmall<-oldx[,c(keys,fields),with=FALSE];newsmall<-newx[,c(keys,fields),with=FALSE]
setnames(oldsmall,fields,paste0('old_',fields));setnames(newsmall,fields,paste0('new_',fields))
m<-merge(oldsmall,newsmall,by=keys,all=TRUE,sort=FALSE)
m[,membership:=fifelse(is.na(old_coloc_unit_id),'NEW_ONLY',fifelse(is.na(new_coloc_unit_id),'OLD_ONLY','SHARED'))]
shared<-m[membership=='SHARED']
for(p in pp){o<-paste0('old_',p);n<-paste0('new_',p);if(!all(c(o,n)%in%names(shared)))stop('Missing posterior: ',p,call.=FALSE);set(shared,j=paste0('delta_',p),value=as.numeric(shared[[n]])-as.numeric(shared[[o]]))}
deltas<-paste0('delta_',pp)
shared[,max_abs_posterior_delta:=do.call(pmax,c(lapply(.SD,abs),na.rm=FALSE)),.SDcols=deltas]
shared[,posterior_changed:=is.na(max_abs_posterior_delta)|max_abs_posterior_delta>1e-6]
if(all(c('old_interpretation','new_interpretation')%in%names(shared)))shared[,interpretation_changed:=is.na(old_interpretation)!=is.na(new_interpretation)|(!is.na(old_interpretation)&!is.na(new_interpretation)&old_interpretation!=new_interpretation)]
write<-function(x,n)fwrite(x,file.path(out,n),sep='\t',quote=FALSE,na='NA')
write(m[membership=='OLD_ONLY'],'old_only.tsv');write(m[membership=='NEW_ONLY'],'new_only.tsv')
write(shared[posterior_changed==TRUE],'shared_posterior_changed.tsv')
write(shared,'shared_comparison.tsv.gz')
byclump<-m[,.(old_units=sum(membership!='NEW_ONLY'),new_units=sum(membership!='OLD_ONLY'),shared_units=sum(membership=='SHARED')),by=.(clump_id,tissue)][order(clump_id,tissue)]
byclump[,net_change:=new_units-old_units];write(byclump,'clump_tissue_changes.tsv')
byclass<-if('interpretation_changed'%in%names(shared))shared[, .N,by=.(old_interpretation,new_interpretation)][order(-N)] else data.table(note='No common interpretation column')
write(byclass,'shared_interpretation_transitions.tsv')
metrics<-data.table(metric=c('old_rows','new_rows','shared_keys','old_only','new_only','net_new_minus_old','shared_posterior_changed','shared_interpretation_changed','old_ids_reused_for_different_key'),value=c(nrow(oldx),nrow(newx),nrow(shared),m[membership=='OLD_ONLY',.N],m[membership=='NEW_ONLY',.N],nrow(newx)-nrow(oldx),shared[posterior_changed==TRUE,.N],if('interpretation_changed'%in%names(shared))shared[interpretation_changed==TRUE,.N] else NA_integer_,sum(shared$old_coloc_unit_id!=shared$new_coloc_unit_id)))
write(metrics,'audit.tsv');cat('Reconciliation:\n');print(metrics);cat('Output: ',out,'\n',sep='')
