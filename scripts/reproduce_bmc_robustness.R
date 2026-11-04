script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),'bootstrap.R'))
args<-commandArgs(TRUE)
registry<-c(setNames(rep('table',6),c('dose_omission_metrics','dose_omission_candidates','dose_omission_joint_candidates',
 'dose_omission_discoveries','dose_omission_patterns','dose_omission_gene_ranges')),figure_dose_omission='figure')
if(length(args)!=1L||!args%in%c('all','--list',names(registry)))stop('Provide one output ID, all or --list')
if(args=='--list'){cat(paste(names(registry),registry,sep=' : '),sep='\n');quit(status=0)}
root<-normalizePath(Sys.getenv('SALT_DATA_ROOT'),mustWork=TRUE)
raw_root<-normalizePath(Sys.getenv('SALT_RAW_ROOT'),mustWork=TRUE)
out<-Sys.getenv('SALT_OUTPUT_ROOT');cache<-Sys.getenv('SALT_QC_CACHE')
if(!nzchar(out)||!nzchar(cache))stop('Set SALT_OUTPUT_ROOT and SALT_QC_CACHE')
inside<-function(a,b)identical(a,b)||startsWith(a,paste0(b,'/'))
for(p in c(out,cache)){
 p<-normalizePath(p,mustWork=FALSE)
 if(any(vapply(c(root,raw_root,project_root),function(q)inside(p,q)||inside(q,p),TRUE)))stop('Separate code, inputs, output and cache roots')
 dir.create(p,recursive=TRUE,showWarnings=FALSE)
}
paths<-c(normalization_paths(raw_root),setNames(file.path(root,vapply(input_specs[c('microarray_root','microarray_shoot')],`[[`,'','path')),c('microarray_root','microarray_shoot')),
 normalization=file.path(root,'NORMALIZATION_PROVENANCE.R'))
hash<-function(p)digest::digest(file=p,algo='sha256')
before<-vapply(paths,hash,'')
code_paths<-source_files(project_root)
code_before<-setNames(vapply(code_paths,hash,''),substring(code_paths,nchar(project_root)+2L))
model_files<-file.path(project_root,c('R/bmc_robustness.R','R/normalization.R','R/inputs.R','R/revision_models.R','R/integration_models.R','renv.lock'))
fingerprint<-digest::digest(list(before,vapply(model_files,hash,''),R.version.string,as.character(packageVersion('limma'))),algo='sha256')
cachefile<-file.path(cache,paste0(fingerprint,'.rds'))
if(file.exists(cachefile)){
 saved<-readRDS(cachefile);stopifnot(identical(saved$fingerprint,fingerprint),identical(saved$sha256,digest::digest(saved$value,algo='sha256')))
 z<-saved$value;rm(saved);message('Verified pair-omission cache')
}else{
 inputs<-setNames(lapply(c('microarray_root','microarray_shoot'),read_input,data_root=root,specs=input_specs),c('microarray_root','microarray_shoot'))
 source<-normalization_sources(raw_root)
 z<-bmc_dose_omissions(inputs,source)
 saveRDS(list(fingerprint=fingerprint,value=z,sha256=digest::digest(z,algo='sha256')),cachefile)
}
revision<-system2('git',c('-C',shQuote(project_root),'rev-parse','HEAD'),stdout=TRUE)
for(id in if(args=='all')names(registry) else args){
 folder<-file.path(out,id);dir.create(folder,showWarnings=FALSE)
 marker<-file.path(folder,'RUN_COMPLETE');if(file.exists(marker))unlink(marker)
 path<-file.path(folder,paste0(id,if(registry[[id]]=='table')'.csv' else '.pdf'))
 if(registry[[id]]=='table')write_table(z[[id]],path) else bmc_omission_figure(z,path)
 outputs<-if(registry[[id]]=='table')path else c(paste0(sub('\\.pdf$','',path),'.',c('pdf','svg','png','tiff')),sub('\\.pdf$','_caption.txt',path))
 stopifnot(all(file.exists(outputs)),identical(before,vapply(paths,hash,'')),
  identical(unname(code_before),unname(vapply(code_paths,hash,''))))
 provenance<-list(artifact_id=id,git_revision=revision,code_sha256=code_before,input_sha256=before,
  output_sha256=setNames(vapply(outputs,hash,''),basename(outputs)),fingerprint=fingerprint,
  specification='All twelve harvested pair omissions, frozen and refitted quantile normalization; no primary exclusion; fixed 27089-locus/15-contrast family and separate three-joint-test family',
  normalization_reconstruction_max_abs_error=z$reconstruction_max_abs_error,
  command=commandArgs(),session=capture.output(sessionInfo()))
 dput(provenance,file=file.path(folder,'provenance.R'));writeLines('Complete; see provenance.R',marker)
 message('Completed ',id)
}
