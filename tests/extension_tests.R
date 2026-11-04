script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),"..","scripts","bootstrap.R"))
source(file.path(project_root,"config","extensions.R"))
set.seed(918)
# Contrast arithmetic and the twelve paired harvest units must survive expansion.
ids<-paste0("AT1G",sprintf("%05d",seq_len(120)),".1")
m<-expand.grid(replicate=1:3,dose=c(0,50,75,100))
pool_noise<-matrix(rnorm(120*12,sd=.1),120)
base<-matrix(rnorm(120*12,10,.15),120)+pool_noise
root<-sweep(base,2,m$dose*.01,"+");shoot<-sweep(base,2,-m$dose*.005,"+")
input<-function(y,tissue){
  ss<-paste(ifelse(m$dose==0,"PN",paste(m$dose,"NaCl")),tolower(tissue),paste0("B",m$replicate))
  list(values=2^y,ids=ids,mapped=rep(TRUE,length(ids)),samples=ss,spec=list(tissue=tissue))
}
z<-extension_dose(list(microarray_root=input(root,"Root"),microarray_shoot=input(shoot,"Shoot")))
stopifnot(length(unique(z$metadata$pool))==12L,all(table(z$metadata$pool)==2L),ncol(z$contrasts)==15L)
a<-z$table[z$table$contrast=="Difference_75v50",]
manual<-(rowMeans(shoot[,m$dose==75])-rowMeans(shoot[,m$dose==50]))-
  (rowMeans(root[,m$dose==75])-rowMeans(root[,m$dose==50]))
stopifnot(max(abs(a$log2FC[match(sub("\\.1$","",ids),a$feature)]-manual))<1e-10)
stopifnot(identical(z$table$q_global,p.adjust(z$table$P,"BH")))
# Ambiguous mappings are excluded, never assigned to the first AGI.
mapping<-data.frame(probe=c("a","b","c"),locus=c("AT1G00010","AT1G00010",NA),mapping_status=c("unique","unique","ambiguous"))
y<-matrix(1:6,3,dimnames=list(c("a","b","c"),c("s1","s2")))
v<-extension_map_matrix(y,mapping)
stopifnot(nrow(v$y)==1L,all(v$y[1,]==c(1.5,4.5)),is.na(v$audit$locus[3]))
# Fractional rank scores have one value per biological sample, with fixed measured membership.
s<-extension_rank_scores(y,list(low="a",high="c",absent="x"),1L)
stopifnot(identical(dim(s),c(2L,2L)),all(s["high",]==1),all(s["low",]==1/3),!"absent"%in%rownames(s))
# Selection of plotted cells must reject duplicates and preserve untested cells.
d<-data.frame(set=c("one","two"),dataset="Array2",contrast="salt_WT")
cols<-data.frame(dataset="Array2",contrast=c("salt_WT","interaction"),key=paste("Array2",c("salt_WT","interaction"),sep="::"))
mm<-extension_matrix(d,"set",c("one","two"),cols,c(2,-3))
stopifnot(identical(unname(mm[,1]),c(2,-3)),all(is.na(mm[,2])))
stopifnot(inherits(try(extension_matrix(rbind(d,d[1,]),"set",c("one","two"),cols,1:3),silent=TRUE),"try-error"))
# Streaming an explicitly opened ZIP connection must reach EOF, ignore malformed
# rows, strip AGI case/whitespace, and accept only the original DAP paper.
tmp<-tempfile();dir.create(tmp);f<-file.path(tmp,"AtRegNet.csv")
header<-"TFName,TFLocus,TFFamily,TargetName,TargetLocus,IsTF,TargetType,Confirmation,Category,InteractionType,Reference,Note,PubMedID"
writeLines(c(header,"TF, at1g00010 ,bZIP,G, AT2G00020 ,No, Direct ,Unconfirmed,NA,Unknown,NA, DAP-Seq , 27203113",
  "TF,AT1G00010,bZIP,G,AT2G00030,No,Direct,Confirmed,NA,Unknown,NA,ChIP,27203113",
  "bad,row,27203113"),f)
old<-getwd();setwd(tmp);utils::zip("fixture.zip","AtRegNet.csv",flags="-q");setwd(old)
dap<-extension_read_dap(file.path(tmp,"fixture.zip"),1L)
stopifnot(nrow(dap$edges)==1L,dap$edges$TF=="AT1G00010",dap$edges$target=="AT2G00020",
  dap$audit$candidate_rows==3L,dap$audit$excluded_rows==2L)
unlink(tmp,recursive=TRUE)
cat("Extension contrast, pairing, mapping, rank-score, heatmap, and reference-stream tests passed.\n")
