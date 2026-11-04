script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),'bootstrap.R'))
source(file.path(project_root,'config','extensions.R'))
if(length(commandArgs(TRUE)))stop('No arguments; set external roots')
root<-normalizePath(Sys.getenv('SALT_DATA_ROOT'),mustWork=TRUE)
out<-normalizePath(Sys.getenv('SALT_OUTPUT_ROOT'),mustWork=FALSE)
if(!nzchar(out))stop('Set SALT_OUTPUT_ROOT')
public<-extension_paths();dapfile<-public[['AtRegNet.zip']]
inside<-function(a,b)identical(a,b)||startsWith(a,paste0(b,'/'))
if(any(vapply(c(root,dirname(dapfile),project_root),function(x)inside(out,x)||inside(x,out),TRUE)))stop('Keep code, inputs, references and outputs separate')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ids<-c('microarray_root','microarray_shoot','microarray_genotype','proteomics_root','proteomics_shoot','saved_proteomics_root','saved_proteomics_shoot')
paths<-c(setNames(file.path(root,vapply(input_specs[ids],`[[`,'','path')),ids),DAP=dapfile)
hash<-function(p)digest::digest(file=p,algo='sha256')
before<-vapply(paths,hash,'');code<-source_files(project_root)
code_before<-setNames(vapply(code,hash,''),substring(code,nchar(project_root)+2L))
inputs<-setNames(lapply(ids,read_input,data_root=root,specs=input_specs),ids)
dose<-extension_dose(inputs);objects<-extension_objects(inputs,dose)
dap<-extension_read_dap(dapfile);result<-bmc_tf_correlation(dap,objects)
stopifnot(nrow(result$TF_estimated_correlation_tests)==5935L)
revision<-system2('git',c('-C',shQuote(project_root),'rev-parse','HEAD'),stdout=TRUE)
for(id in names(result)){
 folder<-file.path(out,id);dir.create(folder,showWarnings=FALSE)
 marker<-file.path(folder,'RUN_COMPLETE');if(file.exists(marker))unlink(marker)
 path<-file.path(folder,paste0(id,'.csv'));write_table(result[[id]],path)
 stopifnot(identical(before,vapply(paths,hash,'')),identical(unname(code_before),unname(vapply(code,hash,''))))
 dput(list(artifact_id=id,git_revision=revision,input_sha256=before,code_sha256=code_before,
  output_sha256=setNames(hash(path),basename(path)),command=commandArgs(),
  specification='Same moderated feature statistics and DAP-seq universes as primary tests. Estimate target-set average residual correlation after paired-design whitening; truncate negative estimates to zero. Full separate BH/BY family. Sensitivity, not inferred TF activity.',
  session=capture.output(sessionInfo())),file=file.path(folder,'provenance.R'))
 writeLines('Complete; see provenance.R',marker);message('Completed ',id)
}
writeLines('Complete; each artifact has provenance.R',file.path(out,'RUN_COMPLETE'))
