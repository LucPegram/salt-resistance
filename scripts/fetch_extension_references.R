script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),"bootstrap.R"))
source(file.path(project_root,"config","extensions.R"))
root<-Sys.getenv("SALT_REFERENCE_ROOT")
if(!nzchar(root))stop("Set SALT_REFERENCE_ROOT to an external reference directory")
dir.create(root,recursive=TRUE,showWarnings=FALSE)
paths<-extension_paths();urls<-extension_resources();options(timeout=600)
manifest<-list()
for(n in names(urls)){
  p<-paths[[n]]
  if(!file.exists(p)){
    tmp<-paste0(p,".partial")
    tryCatch({download.file(urls[[n]],tmp,mode="wb",method="libcurl");file.rename(tmp,p)},
      error=function(e){unlink(tmp);stop(e)})
  }
  if(file.info(p)$size<1000)stop("Incomplete reference: ",n)
  if(grepl("(xlsx|zip)$",n)){
    if(!length(unzip(p,list=TRUE)$Name))stop("Invalid ZIP resource: ",n)
  }else if(grepl("gz$",n))extension_read_matrix(p)
  else if(!any(startsWith(readLines(p,n=20,warn=FALSE),"^PLATFORM")))stop("Invalid GEO platform file")
  manifest[[n]]<-data.frame(resource=n,url=urls[[n]],bytes=file.info(p)$size,
    sha256=digest::digest(file=p,algo="sha256"),recorded_utc=format(Sys.time(),tz="UTC",usetz=TRUE))
  message("Verified ",n)
}
write.csv(do.call(rbind,manifest),file.path(root,"download_manifest.csv"),row.names=FALSE)
