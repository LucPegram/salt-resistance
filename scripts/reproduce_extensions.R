script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),"bootstrap.R"))
source(file.path(project_root,"config","extensions.R"))
args<-commandArgs(TRUE);registry<-extension_registry()
if(length(args)!=1L)stop("Provide one artifact ID, 'all', or '--list'")
if(args=="--list"){cat(paste(names(registry),registry,sep=" : "),sep="\n");quit(status=0)}
ids<-if(args=="all")names(registry) else args
if(any(!ids%in%names(registry)))stop("Unknown extension artifact")
data_root<-normalizePath(Sys.getenv("SALT_DATA_ROOT"),mustWork=TRUE)
if(!file.exists(file.path(data_root,"NORMALIZATION_PROVENANCE.R")))stop("Use the primary quantile-normalized input root")
paths<-extension_paths();if(any(!file.exists(paths)))stop("Fetch all extension references first")
output_root<-Sys.getenv("SALT_OUTPUT_ROOT");if(!nzchar(output_root))stop("Set SALT_OUTPUT_ROOT")
dir.create(output_root,recursive=TRUE,showWarnings=FALSE);output_root<-normalizePath(output_root,mustWork=TRUE)
inside<-function(a,b)identical(a,b)||startsWith(a,paste0(b,"/"))
if(inside(output_root,project_root)||inside(output_root,data_root)||inside(data_root,output_root)||
   inside(output_root,dirname(paths[[1]]))||inside(dirname(paths[[1]]),output_root))stop("Keep output, code, data, and references separate")
input_ids<-c("microarray_root","microarray_shoot","microarray_genotype","proteomics_root","proteomics_shoot",
  "saved_proteomics_root","saved_proteomics_shoot","gene_annotation")
input_paths<-setNames(vapply(input_ids,function(n)file.path(data_root,input_specs[[n]]$path),""),input_ids)
hash<-function(p)digest::digest(file=p,algo="sha256")
code_paths<-source_files(project_root);code_before<-setNames(vapply(code_paths,hash,""),substring(code_paths,nchar(project_root)+2L))
input_hashes<-vapply(input_paths,hash,"");reference_hashes<-vapply(paths,hash,"")
fingerprint<-digest::digest(list(code_before,input_hashes,reference_hashes,extension_settings(),
  R.version.string,vapply(c("limma","statmod","readxl"),function(p)as.character(packageVersion(p)),"")),algo="sha256")
cachefile<-file.path(output_root,"derived_analysis_cache.rds")
cached<-if(file.exists(cachefile))readRDS(cachefile) else NULL
if(!is.null(cached)&&identical(cached$fingerprint,fingerprint)&&
   identical(cached$sha256,digest::digest(cached$value,algo="sha256"))){z<-cached$value;rm(cached);message("Verified extension cache")
}else{
  inputs<-setNames(lapply(input_ids,read_input,data_root=data_root,specs=input_specs),input_ids)
  set.seed(extension_settings()$seed);z<-extension_compute(inputs,paths)
  saveRDS(list(fingerprint=fingerprint,value=z,sha256=digest::digest(z,algo="sha256")),cachefile,compress=FALSE)
}
tables<-extension_tables(z)
revision<-system2("git",c("-C",shQuote(project_root),"rev-parse","HEAD"),stdout=TRUE)
for(id in ids){
  folder<-file.path(output_root,id);dir.create(folder,showWarnings=FALSE)
  marker<-file.path(folder,"RUN_COMPLETE");if(file.exists(marker))unlink(marker)
  type<-registry[[id]];path<-file.path(folder,paste0(id,if(type=="table")".csv" else ".pdf"))
  if(type=="table")write_table(tables[[id]],path) else switch(type,
    dose=extension_dose_figure(z,path),candidates=extension_candidate_figure(z,path),
    benchmark=extension_benchmark_figure(z,path),concordance=extension_benchmark_figure(z,path,TRUE),
    TF=extension_tf_figure(z,path),markers=extension_marker_figure(z,path),report=extension_report(z,path))
  files<-if(type%in%c("table","report"))path else paste0(sub("\\.pdf$","",path),".",c("pdf","svg","png","tiff"))
  if(!type%in%c("table","report"))files<-c(files,sub("\\.pdf$","_caption.txt",path))
  stopifnot(all(file.exists(files)),all(file.info(files)$size>0))
  if(!identical(input_hashes,vapply(input_paths,hash,""))||!identical(reference_hashes,vapply(paths,hash,""))||
     !identical(unname(code_before),unname(vapply(code_paths,hash,""))))stop("Input/reference/code changed during run")
  provenance<-list(artifact_id=id,generated_utc=format(Sys.time(),tz="UTC",usetz=TRUE),git_revision=revision,
    code_sha256=code_before,input_sha256=input_hashes,input_relative_paths=vapply(input_ids,function(n)input_specs[[n]]$path,""),
    reference_sha256=reference_hashes,reference_urls=extension_resources(),settings=extension_settings(),
    normalization=dget(file.path(data_root,"NORMALIZATION_PROVENANCE.R")),output_sha256=setNames(vapply(files,hash,""),basename(files)),
    command=commandArgs(),session=capture.output(sessionInfo()),fingerprint=fingerprint)
  dput(provenance,file=file.path(folder,"provenance.R"));writeLines("Complete; see provenance.R",marker)
  message("Completed ",id)
}
