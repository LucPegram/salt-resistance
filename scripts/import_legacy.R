# Explicit maintenance command: copy R source only; never copy measurements.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
source(file.path(dirname(script),"bootstrap.R"))
root<-normalizePath(Sys.getenv("SALT_DATA_ROOT"),mustWork=TRUE)
locations<-c(microarray="Microarray_for_Stephano/Luc_Initial_Analysis/R_scripts",
             proteomics="Data_from_Luc_2/R_Scripts_Proteomics")
records<-list()
for(group in names(locations)){
  directory<-file.path(root,locations[group])
  files<-sort(list.files(directory,pattern="\\.[Rr]$",recursive=TRUE,full.names=TRUE))
  for(p in files){
    relative<-substring(p,nchar(directory)+2L)
    clean<-gsub("[^A-Za-z0-9_./-]","_",relative)
    destination<-file.path(project_root,"legacy",group,clean)
    if(file.exists(destination))stop("Refusing to replace archived source: ",destination)
    dir.create(dirname(destination),recursive=TRUE,showWarnings=FALSE)
    if(!file.copy(p,destination))stop("Failed copy")
    original<-digest::digest(file=p,algo="sha256")
    stopifnot(identical(original,digest::digest(file=destination,algo="sha256")))
    records[[length(records)+1L]]<-c(source=file.path(locations[group],relative),
      archived=substring(destination,nchar(project_root)+2L),sha256=original)
  }
}
header<-c("# Archived source scripts","",
  "These are unmodified historical R scripts from the external source archives.",
  "Filenames have been made portable; file bytes and original author comments are preserved.",
  "They are not run by the current workflow. Hard-coded paths, obsolete model choices,",
  "and incomplete exports remain as historical evidence, not endorsed analysis defaults.",
  "Execution cannot be inferred from the presence of a script. Supplied results are",
  "inventoried separately in local run outputs. Do not execute these scripts over original inputs.","",
  "| Original source path | Archived code path | SHA256 |","|---|---|---|")
lines<-vapply(records,function(r)paste0("| `",r[1],"` | `",r[2],"` | `",r[3],"` |"),"")
writeLines(c(header,lines),file.path(project_root,"legacy","README.md"))
cat("Archived",length(records),"R scripts; measurements were not copied.\n")
