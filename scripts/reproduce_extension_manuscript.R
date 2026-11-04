# Format external-baseline, TF and marker display summaries from complete extension CSVs.
# No models, thresholds, or original expression measurements are changed.
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),'bootstrap.R'))
args<-commandArgs(TRUE)
ids<-c('table_external_baseline','table_TF_sensitivity','table_marker_baseline')
if(length(args)!=1L||!args%in%c('all','--list',ids))stop('Provide one table ID, all, or --list')
if(args=='--list'){cat(ids,sep='\n');quit(status=0)}
input_root<-Sys.getenv('SALT_EXTENSION_ROOT');output_root<-Sys.getenv('SALT_OUTPUT_ROOT')
if(!nzchar(input_root)||!nzchar(output_root))stop('Set SALT_EXTENSION_ROOT and SALT_OUTPUT_ROOT')
input_root<-normalizePath(input_root,mustWork=TRUE)
dir.create(output_root,recursive=TRUE,showWarnings=FALSE)
output_root<-normalizePath(output_root,mustWork=TRUE)
inside<-function(a,b)identical(a,b)||startsWith(a,paste0(b,'/'))
if(inside(input_root,project_root)||inside(output_root,project_root)||inside(input_root,output_root)||inside(output_root,input_root))
 stop('Keep code, source tables, and formatted tables in separate directories')
names_in<-c('external_signature_tests','TF_target_enrichment','marker_tests')
paths<-setNames(file.path(input_root,names_in,paste0(names_in,'.csv')),names_in)
if(any(!file.exists(paths))||any(!file.exists(file.path(dirname(paths),'RUN_COMPLETE'))))stop('Incomplete extension sources')
hash<-function(p)digest::digest(file=p,algo='sha256')
before<-vapply(paths,hash,'');code_paths<-source_files(project_root)
code_hashes<-setNames(vapply(code_paths,hash,''),substring(code_paths,nchar(project_root)+2L))
d<-lapply(paths,read.csv,check.names=FALSE,stringsAsFactors=FALSE)
b<-d$external_signature_tests
b<-b[b$contrast=='baseline_genotype',]
b<-b[order(b$dataset,b$reference_tissue,as.numeric(sub('.*hour','',b$reference))),]
key<-function(x)paste(x$dataset,x$contrast,if('TF'%in%names(x))x$TF else x$set,sep='::')
t<-d$TF_target_enrichment;p<-t[t$mode=='primary'&t$q_global<.05,]
s<-t[t$mode=='sensitivity',];k<-match(key(p),key(s));stopifnot(!anyNA(k),!anyDuplicated(key(s)))
p$q_global_correlation_005<-s$q_global[k];p$q_BY_correlation_005<-s$q_BY[k]
m<-d$marker_tests;a<-m[m$mode=='primary'&m$dataset=='Array2'&m$contrast=='baseline_genotype'&m$q_global<.05,]
s<-m[m$mode=='without_focal_candidates',];k<-match(key(a),key(s));stopifnot(!anyNA(k),!anyDuplicated(key(s)))
a$q_global_without_focal<-s$q_global[k];a$q_BY_without_focal<-s$q_BY[k]
a<-a[order(a$q_global,a$set),]
values<-list(table_external_baseline=b,table_TF_sensitivity=p,table_marker_baseline=a)
revision<-system2('git',c('-C',shQuote(project_root),'rev-parse','HEAD'),stdout=TRUE)
for(id in if(args=='all')ids else args){
 folder<-file.path(output_root,id);dir.create(folder,showWarnings=FALSE)
 marker<-file.path(folder,'RUN_COMPLETE');if(file.exists(marker))unlink(marker)
 path<-file.path(folder,paste0(id,'.csv'));write_table(values[[id]],path)
 stopifnot(identical(before,vapply(paths,hash,'')),identical(unname(code_hashes),unname(vapply(code_paths,hash,''))))
 provenance<-list(artifact_id=id,git_revision=revision,generated_utc=format(Sys.time(),tz='UTC',usetz=TRUE),
  source_sha256=before,source_paths=paths,code_sha256=code_hashes,output_sha256=setNames(hash(path),basename(path)),
  selection=switch(id,table_external_baseline='Every eligible baseline genotype salt-signature test, including nonsignificant tests',
   table_TF_sensitivity='Every primary global-BH TF finding, with the matched correlation-0.05 sensitivity',
   table_marker_baseline='Every primary global-BH Array2 baseline marker finding, with matched focal-gene exclusion sensitivity'),
  command=commandArgs(),session=capture.output(sessionInfo()))
 dput(provenance,file=file.path(folder,'provenance.R'));writeLines('Complete; see provenance.R',marker)
 message('Completed ',id)
}
