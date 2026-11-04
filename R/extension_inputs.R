# Independent resources are stored in SALT_REFERENCE_ROOT, outside the repository.
extension_paths <- function() {
  p<-Sys.getenv("SALT_REFERENCE_ROOT")
  if(!nzchar(p))stop("Set SALT_REFERENCE_ROOT outside the code and experimental inputs")
  p<-normalizePath(p,mustWork=TRUE)
  if(identical(p,project_root)||startsWith(p,paste0(project_root,"/")))stop("References must be outside code")
  setNames(file.path(p,names(extension_resources())),names(extension_resources()))
}

extension_read_matrix <- function(path) {
  lines<-readLines(gzfile(path),warn=FALSE)
  a<-match("!series_matrix_table_begin",lines);b<-match("!series_matrix_table_end",lines)
  if(is.na(a)||is.na(b))stop("No GEO matrix: ",path)
  d<-read.delim(text=paste(lines[(a+1L):(b-1L)],collapse="\n"),check.names=FALSE)
  y<-as.matrix(d[-1]);rownames(y)<-as.character(d[[1]])
  if(anyDuplicated(rownames(y))||any(!is.finite(y)))stop("Invalid external expression matrix")
  field<-function(key){
    l<-lines[startsWith(lines,paste0(key,"\t"))]
    if(length(l)!=1L)stop("Expected one GEO metadata field: ",key)
    as.character(read.delim(text=l,header=FALSE,check.names=FALSE,quote='"')[-1])
  }
  list(y=y,title=field("!Sample_title"),source=field("!Sample_source_name_ch1"),
    processing=field("!Sample_data_processing"),sample=colnames(y))
}

extension_read_platform <- function(path,field) {
  con<-file(path,"rt");on.exit(close(con))
  repeat{line<-readLines(con,n=1,warn=FALSE);if(!length(line))stop("No platform table")
    if(line=="!platform_table_begin")break}
  header<-strsplit(readLines(con,n=1),"\t",fixed=TRUE)[[1]]
  keep<-match(c("ID",field),header);if(anyNA(keep))stop("Invalid platform annotation")
  classes<-rep("NULL",length(header));classes[keep]<-"character"
  d<-read.delim(con,header=FALSE,col.names=header,colClasses=classes,quote="",fill=TRUE,
    comment.char="!",check.names=FALSE,stringsAsFactors=FALSE)
  names(d)<-c("probe","annotation")
  hits<-regmatches(toupper(d$annotation),gregexpr("AT[1-5CM]G[0-9]{5}",toupper(d$annotation)))
  hits<-lapply(hits,unique)
  d$locus<-vapply(hits,function(v)if(length(v)==1L)v else NA_character_,"")
  d$mapping_status<-ifelse(lengths(hits)==1L,"unique locus",ifelse(lengths(hits)==0,"no AGI","ambiguous loci"))
  d
}

extension_map_matrix <- function(y,mapping) {
  if(anyDuplicated(mapping$probe))stop("Nonunique platform probe IDs")
  locus<-mapping$locus[match(rownames(y),mapping$probe)]
  ok<-!is.na(locus);n<-table(locus[ok])
  z<-rowsum(y[ok,,drop=FALSE],locus[ok],reorder=FALSE)
  z<-z/as.numeric(n[rownames(z)])
  stopifnot(all(is.finite(z)),!anyDuplicated(rownames(z)))
  list(y=z,audit=data.frame(probe=rownames(y),locus=locus,
    mapping_status=mapping$mapping_status[match(rownames(y),mapping$probe)]))
}

extension_read_dap <- function(path,min_edges=10000L) {
  # AGRIS includes malformed historical CSV rows. DAP-seq records used here must
  # have exactly 13 fields and explicit original-paper provenance; never shift columns.
  con<-unz(path,"AtRegNet.csv",open="r");on.exit(close(con));header<-readLines(con,n=1)
  stopifnot(startsWith(header,"TFName,TFLocus,TFFamily,TargetName,TargetLocus"))
  rows<-list();seen<-excluded<-0L
  repeat{
    x<-readLines(con,n=20000,warn=FALSE);if(!length(x))break
    x<-x[grepl("27203113",x,fixed=TRUE)];seen<-seen+length(x)
    a<-strsplit(x,",",fixed=TRUE);valid<-lengths(a)==13L;excluded<-excluded+sum(!valid)
    a<-a[valid];if(!length(a))next
    m<-matrix(trimws(unlist(a,use.names=FALSE)),ncol=13,byrow=TRUE)
    ok<-grepl("^AT[1-5CM]G[0-9]{5}$",toupper(m[,2])) &
      grepl("^AT[1-5CM]G[0-9]{5}$",toupper(m[,5])) &
      tolower(m[,7])=="direct" & m[,13]=="27203113" & grepl("DAP",m[,12],ignore.case=TRUE)
    excluded<-excluded+sum(!ok)
    rows[[length(rows)+1L]]<-data.frame(TF=toupper(m[ok,2]),TF_name=m[ok,1],
      TF_family=m[ok,3],target=toupper(m[ok,5]),stringsAsFactors=FALSE)
  }
  d<-unique(do.call(rbind,rows));if(nrow(d)<min_edges)stop("Incomplete DAP reference")
  list(edges=d,audit=data.frame(source="AGRIS; O'Malley et al. 2016 DAP-seq only",
    candidate_rows=seen,excluded_rows=excluded,unique_edges=nrow(d),TFs=length(unique(d$TF))))
}

extension_read_markers <- function(curated,clusters) {
  rows<-list();sheets<-c("Root Cell Type Markers","Mesophyll Markers","Vasculature Markers",
    "SAM Markers","Other Cell Type Markers")
  # Include any epidermal sheet under its actual published name.
  sheets<-unique(c(sheets,grep("Epiderm",readxl::excel_sheets(curated),value=TRUE)))
  for(s in sheets){
    d<-as.data.frame(readxl::read_excel(curated,sheet=s,skip=2,.name_repair="unique"))
    assert_columns(d,c("Gene ID","Cell type","Reference"),s)
    locus<-toupper(trimws(d[["Gene ID"]]));type<-trimws(d[["Cell type"]])
    valid<-grepl("^AT[1-5CM]G[0-9]{5}$",locus)&!is.na(type)&nzchar(type)
    if(s=="Other Cell Type Markers")valid<-valid&grepl("Dividing|G1/S",type)
    category<-if(s=="Root Cell Type Markers")"root_identity" else if(s=="Other Cell Type Markers")"proliferation" else "shoot_identity"
    rows[[length(rows)+1L]]<-data.frame(set=paste(s,type[valid],sep=" / "),locus=locus[valid],
      category=category,source="Lee et al. 2025 Table S1",original_reference=d$Reference[valid])
    if(s=="Root Cell Type Markers")for(stage in c("Meristem","Elongation","Maturation")){
      k<-valid&startsWith(type,paste0(stage,"_"))
      rows[[length(rows)+1L]]<-data.frame(set=paste("Root developmental stage",stage,sep=" / "),
        locus=locus[k],category="root_development",source="Lee et al. 2025 Table S1; Shahan et al. 2022 markers",
        original_reference=d$Reference[k])
    }
  }
  d<-as.data.frame(readxl::read_excel(clusters,sheet="Rosette_21d",skip=2))
  a<-as.data.frame(readxl::read_excel(clusters,sheet="Suggested Cell Type Annotations",skip=2))
  a<-a[!is.na(a$Organ)&a$Organ=="Rosette 21d",]
  if(!nrow(a)||anyDuplicated(a[["Major Cluster"]]))stop("Missing or ambiguous rosette cluster annotations")
  assert_columns(d,c("Gene","avg_log2FC","pct.1","pct.2","p_val_adj","cluster"),"rosette markers")
  label<-a[["Cell type suggested_v3"]][match(d$cluster,a[["Major Cluster"]])]
  s<-extension_settings();valid<-!is.na(label)&label!="Unknown"&d$p_val_adj<s$atlas_max_adj_P &
    d$avg_log2FC>=s$atlas_min_log2FC & d$pct.1>d$pct.2 & grepl("^AT[1-5CM]G[0-9]{5}$",d$Gene)
  d<-d[valid,];label<-label[valid]
  for(cl in sort(unique(d$cluster))){
    k<-which(d$cluster==cl);k<-head(k[order(-d$avg_log2FC[k],d$Gene[k])],s$atlas_top_per_cluster)
    rows[[length(rows)+1L]]<-data.frame(set=paste("Rosette21d",label[k],sep=" / "),locus=d$Gene[k],
      category=ifelse(grepl("[Dd]ividing",label[k]),"proliferation","shoot_identity"),
      source="Lee et al. 2025 Table S2; top 50 eligible markers per annotated cluster",
      original_reference="10.1038/s41477-025-02072-z")
  }
  out<-unique(do.call(rbind,rows));out[order(out$set,out$locus),]
}
