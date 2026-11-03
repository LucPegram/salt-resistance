# Inventory source evidence without executing historical scripts.
script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),"bootstrap.R"))
root<-normalizePath(Sys.getenv("SALT_DATA_ROOT"),mustWork=TRUE)
out<-Sys.getenv("SALT_OUTPUT_ROOT")
if(!nzchar(out))stop("Set SALT_OUTPUT_ROOT outside the source directory")
out<-normalizePath(out,mustWork=FALSE)
if(startsWith(out,paste0(root,"/")) || identical(root,out))stop("Outputs must be separate")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
files<-sort(list.files(root,recursive=TRUE,full.names=TRUE,all.files=TRUE,include.dirs=FALSE))
relative<-substring(files,nchar(root)+2L)
keep<-!grepl("(^|/)(\\.Rproj.user|__MACOSX|\\.DS_Store)(/|$)|(^|/)~\\$",relative)
files<-files[keep];relative<-relative[keep]
role<-ifelse(grepl("\\.[Rr]$",relative),"supplied source code; execution unverified",
 ifelse(grepl("\\.(pdf|png|jpg)$",relative,ignore.case=TRUE),"supplied document or visualization",
 ifelse(grepl("results|summary|R_output|R_analysis|Heatmap|HeatMap|DEG",relative),"supplied analysis output or selected subset",
 "measurement, annotation or other source; inspect before use")))
d<-data.frame(path=relative,bytes=file.info(files)$size,
  sha256=vapply(files,function(p)digest::digest(file=p,algo="sha256"),""),role=role)
write.csv(d,file.path(out,"source_inventory.csv"),row.names=FALSE)
scripts<-files[grepl("\\.[Rr]$",relative)]
report<-lapply(scripts,function(p){
 text<-readLines(p,warn=FALSE)
 packages<-unique(unlist(regmatches(text,gregexpr("(?<=(library|require)\\()[A-Za-z0-9.]+",text,perl=TRUE))))
 parsed<-tryCatch({parse(p);"parses"},error=function(e)paste("parse error:",conditionMessage(e)))
 data.frame(path=substring(p,nchar(root)+2L),parse_status=parsed,
   referenced_packages=paste(sort(packages),collapse=";"),
   has_absolute_windows_path=any(grepl("[A-Za-z]:/",text)),
   execution_status="Not executed; presence is not proof of completed analysis")
})
write.csv(do.call(rbind,report),file.path(out,"legacy_analysis_inventory.csv"),row.names=FALSE)
cat("Inventoried",nrow(d),"source files and",length(scripts),"R scripts.\n")
