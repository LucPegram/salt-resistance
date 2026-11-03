script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),"..","scripts","bootstrap.R"))
# Synthetic values verify correction BEFORE selection, sign and distinct loci.
d<-expand.grid(locus=paste0("AT1G",sprintf("%05d",1:50)),tissue=c("Root","Shoot"),stringsAsFactors=FALSE)
d$log2FC<-rep(c(seq(1,25),-seq(1,25)),2);d$log2FC[51:100]<-d$log2FC[51:100]/2
d$P<-seq_len(100)/1e6;d$eligible<-TRUE
a<-toprank_select(d)
stopifnot(max(abs(a$q_BH-p.adjust(d$P,"BH")))<1e-15,
  sum(a$selected_top20)==40L,sum(a$plotted_top10)==20L,
  !anyDuplicated(a$locus[a$selected_top20]),
  all(a$tissue[a$selected_top20]=="Root"),
  all(a$direction[a$selected_top20&a$log2FC>0]=="Up"),
  all(a$direction[a$selected_top20&a$log2FC<0]=="Down"))
for(direction in c("Up","Down")){
  z<-a[a$selected_top20&a$direction==direction,];z<-z[order(z$rank),]
  stopifnot(all(diff(abs(z$log2FC))<=0),identical(z$rank,1:20))
}
d$eligible[1]<-FALSE
b<-toprank_select(d)
stopifnot(is.na(b$q_BH[1]),!b$selected_top20[1],
  max(abs(b$q_BH[-1]-p.adjust(d$P[-1],"BH")))<1e-15)
annot<-list(data=data.frame('Gene ID'=c("AT1G00001.1","AT1G00001.2","AT1G00002.1"),
  'Gene info'=c(": TEST1, ALIAS1 | test description | chr1:1",": TEST1, ALIAS1 | test description | chr1:1",
    "| uncharacterized protein; FUNCTIONS IN: unknown | chr1:2"),check.names=FALSE))
m<-toprank_annotation(annot)
stopifnot(nrow(m)==2L,m$Gene_name[m$locus=="AT1G00001"]=="TEST1",
  m$Aliases[m$locus=="AT1G00001"]=="TEST1; ALIAS1",
  m$Gene_name[m$locus=="AT1G00002"]=="uncharacterized protein [AT1G00002]")
# Independent linear-model means reproduce the paired dose contrast estimates.
set.seed(921)
make<-function(tissue){
 samples<-unlist(lapply(c("PN","50 NaCl","75 NaCl","100 NaCl"),
   function(s)paste(s,tolower(tissue),paste0("B",1:3))))
 values<-2^matrix(rnorm(60*12,8,1),nrow=60)
 list(values=values,ids=paste0("AT1G",sprintf("%05d",1:60),".1"),mapped=rep(TRUE,60),
   samples=samples,spec=list(tissue=tissue,assay="microarray"))
}
x<-list(microarray_root=make("Root"),microarray_shoot=make("Shoot"))
z<-toprank_dose(x)
for(tissue in c("Root","Shoot")){
 v<-z$y[1,];m<-z$metadata
 expected<-mean(v[m$tissue==tissue&m$dose==100])-mean(v[m$tissue==tissue&m$dose==0])
 actual<-z$table$log2FC[z$table$locus==rownames(z$y)[1]&z$table$tissue==tissue]
 stopifnot(abs(expected-actual)<1e-12)
}
cat("Ranking multiplicity, distinct-locus selection, annotation and dose contrast checks passed.\n")
