#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
a<-commandArgs(trailingOnly=TRUE);get<-function(k){i<-match(k,a);if(is.na(i)||i==length(a))stop('Missing ',k,call.=FALSE);a[[i+1L]]}
clean<-normalizePath(get('--workspace'),winslash='/',mustWork=TRUE)
old<-normalizePath(get('--old-project'),winslash='/',mustWork=TRUE)
out<-file.path(clean,'results','module12','reconciliation_v59')
if(!dir.exists(out))stop('Run v59 first',call.=FALSE)
read<-function(p){if(!file.exists(p))return(NULL);fread(p,showProgress=FALSE)}
oldbase<-file.path(old,'results','module12_rebuild_staging','colocalization','work_queue')
newbase<-file.path(clean,'results','module12','colocalization','work_queue')
keys<-c('clump_id','tissue','phenotype_id','gene_id')
normalize<-function(x,label){if(is.null(x))return(NULL);if(!all(keys%in%names(x)))stop(label,' lacks keys',call.=FALSE);for(k in keys)set(x,j=k,value=as.character(x[[k]]));x}
oe<-normalize(read(file.path(oldbase,'module12_1_coloc_work_queue_eligible.tsv.gz')),'old eligible')
ne<-normalize(read(file.path(newbase,'module12_1_coloc_work_queue_eligible.tsv.gz')),'new eligible')
of<-normalize(read(file.path(oldbase,'module12_1_coloc_work_queue.tsv.gz')),'old full')
nf<-normalize(read(file.path(newbase,'module12_1_coloc_work_queue.tsv.gz')),'new full')
if(is.null(oe)||is.null(ne)||is.null(nf))stop('Missing eligible or new full queue',call.=FALSE)
if(nrow(oe)!=9471L||nrow(ne)!=9814L||nrow(nf)!=9915L)stop('Unexpected queue count',call.=FALSE)
key<-function(x)do.call(paste,c(x[,..keys],sep='|'))
ok<-key(oe);nk<-key(ne);fk<-key(nf);ofk<-if(!is.null(of))key(of) else character()
oldonly<-oe[!ok%in%nk];newonly<-ne[!nk%in%ok]
oldonly[,other_queue_presence:=ifelse(key(oldonly)%in%fk,'PRESENT_IN_NEW_FULL_QUEUE','ABSENT_FROM_NEW_FULL_QUEUE')]
newonly[,other_queue_presence:=if(!is.null(of))ifelse(key(newonly)%in%ofk,'PRESENT_IN_OLD_FULL_QUEUE','ABSENT_FROM_OLD_FULL_QUEUE') else 'OLD_FULL_QUEUE_UNAVAILABLE']
if('exclusion_reason'%in%names(nf))oldonly[,other_exclusion_reason:=nf$exclusion_reason[match(key(oldonly),fk)]]
if(!is.null(of)&&'exclusion_reason'%in%names(of))newonly[,other_exclusion_reason:=of$exclusion_reason[match(key(newonly),ofk)]]
fwrite(oldonly,file.path(out,'old_only_membership.tsv'),sep='\t',na='NA');fwrite(newonly,file.path(out,'new_only_membership.tsv'),sep='\t',na='NA')
counts<-rbindlist(list(oldonly[,.(N=.N),by=other_queue_presence][,side:='OLD_ONLY'],newonly[,.(N=.N),by=other_queue_presence][,side:='NEW_ONLY']))[order(side,other_queue_presence)]
fwrite(counts,file.path(out,'membership_explanation.tsv'),sep='\t');cat('Queue membership:\n');print(counts)
cat('Old full queue rows: ',if(is.null(of))'UNAVAILABLE' else nrow(of),'\nNew full queue rows: ',nrow(nf),'\n',sep='')
cat('Top clump/tissue differences:\n')
affected<-rbindlist(list(oldonly[,.(side='OLD_ONLY',clump_id,tissue)],newonly[,.(side='NEW_ONLY',clump_id,tissue)]))[, .N,by=.(side,clump_id,tissue)][order(-N)]
fwrite(affected,file.path(out,'membership_by_clump_tissue.tsv'),sep='\t');print(head(affected,20L))
