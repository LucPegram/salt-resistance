script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),'bootstrap.R'))
out<-Sys.getenv('SALT_OUTPUT_ROOT')
registry<-c('external_score_inference_sensitivity','marker_score_inference_sensitivity','table_basal_score_inference','figure_external_score_inference')
args<-commandArgs(TRUE);if(!length(args))args<-'all'
if(length(args)!=1L||!args%in%c('all','--list',registry))stop('Provide one output ID, all or --list')
if(args=='--list'){cat(registry,sep='\n');quit(status=0)}
if(!nzchar(out))stop('Set SALT_OUTPUT_ROOT')
root<-normalizePath(Sys.getenv('SALT_EXTENSION_ROOT'),mustWork=TRUE)
inside<-function(a,b)identical(a,b)||startsWith(a,paste0(b,'/'))
dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out,mustWork=TRUE)
if(any(vapply(c(root,project_root),function(p)inside(out,p)||inside(p,out),TRUE)))stop('Separate code, source and output roots')
ids<-c('project_samples','external_signature_scores','external_signature_tests','marker_scores','marker_tests')
paths<-setNames(file.path(root,ids,paste0(ids,'.csv')),ids)
stopifnot(all(file.exists(paths)),all(file.exists(file.path(dirname(paths),'RUN_COMPLETE'))))
hash<-function(p)digest::digest(file=p,algo='sha256')
before<-vapply(paths,hash,'');code<-source_files(project_root)
code_before<-setNames(vapply(code,hash,''),substring(code,nchar(project_root)+2L))
d<-lapply(paths,read.csv,check.names=FALSE,stringsAsFactors=FALSE)
a<-bmc_unmoderated_scores(d$external_signature_scores,d$external_signature_tests,d$project_samples,'reference')
b<-bmc_unmoderated_scores(d$marker_scores,d$marker_tests,d$project_samples,'set')
stopifnot(nrow(a)==141L,all(table(b$mode)==781L))
basal<-a[a$contrast=='baseline_genotype',]
basal<-basal[order(basal$dataset,sub('^AtGen_([^_]+)_.*','\\1',basal$feature),as.numeric(sub('.*hour','',basal$feature))),]
values<-list(external_score_inference_sensitivity=a,marker_score_inference_sensitivity=b,
 table_basal_score_inference=basal)
revision<-system2('git',c('-C',shQuote(project_root),'rev-parse','HEAD'),stdout=TRUE)
for(id in if(args=='all')registry else args){
 folder<-file.path(out,id);dir.create(folder,showWarnings=FALSE)
 marker<-file.path(folder,'RUN_COMPLETE');if(file.exists(marker))unlink(marker)
 is_figure<-id=='figure_external_score_inference'
 path<-file.path(folder,paste0(id,if(is_figure)'.pdf' else '.csv'))
 if(is_figure)bmc_score_inference_figure(d$external_signature_scores,a,d$project_samples,path) else write_table(values[[id]],path)
 outputs<-if(is_figure)c(paste0(sub('\\.pdf$','',path),'.',c('pdf','svg','png','tiff')),sub('\\.pdf$','_caption.txt',path)) else path
 stopifnot(identical(before,vapply(paths,hash,'')),identical(unname(code_before),unname(vapply(code,hash,''))))
 dput(list(artifact_id=id,git_revision=revision,input_sha256=before,code_sha256=code_before,
  output_sha256=setNames(vapply(outputs,hash,''),basename(outputs)),command=commandArgs(),session=capture.output(sessionInfo()),
  specification='Ordinary residual-variance t sensitivity, preserving complete eligible families and primary score covariance. Not an exact randomization test or independent validation.'),file=file.path(folder,'provenance.R'))
 writeLines('Complete; see provenance.R',marker);message('Completed ',id)
}
