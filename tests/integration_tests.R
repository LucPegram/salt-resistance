script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),"..","scripts","bootstrap.R"))

# Unequal duplicate-spot counts must not overweight one transcript at a locus.
x<-list(ids=c("AT1G00010.1","AT1G00010.1","AT1G00010.2","AT1G00020.1"),
  mapped=rep(TRUE,4),values=2^matrix(c(0,2,10,3,2,4,12,4),ncol=2))
z<-integration_loci(x)
stopifnot(z$y["AT1G00010",1]==5.5,z$y["AT1G00010",2]==7.5,
  z$transcripts["AT1G00010"]==2L,z$transcripts["AT1G00020"]==1L)

# Neither one significant assay nor opposite directions establish conjunction.
p<-integration_conjunction(c(.001,.001,.001,.5),c(1,1,-1,1),
  c(.02,.5,.02,.001),c(1,1,-1,-1))
stopifnot(max(abs(p-c(.02,.5,.02,1)))<1e-12)
stopifnot(integration_conjunction(1e-12,1,1e-12,-1)==1)

# Across-family correction must include every contrast and remain method-specific.
d<-expand.grid(method=c("primary","single_transcript"),comparison="No0_genotype",
  tissue=c("Root","Shoot"),contrast=c("baseline_genotype","interaction"),gene=1:5,
  stringsAsFactors=FALSE)
d$joint_P<-seq_len(nrow(d))/(nrow(d)+1)
a<-integration_adjust(d)
for(method in unique(d$method)){
  k<-d$method==method
  stopifnot(identical(a$joint_q_global[k],p.adjust(d$joint_P[k],"BH")),
    identical(a$joint_q_BY[k],p.adjust(d$joint_P[k],"BY")))
}
for(co in unique(d$contrast))for(t in unique(d$tissue))for(method in unique(d$method)){
  k<-d$contrast==co&d$tissue==t&d$method==method
  stopifnot(identical(a$joint_q_family[k],p.adjust(d$joint_P[k],"BH")))
}

# Mapping tests use invented accessions/AGIs and measurements only.
meta<-data.frame(genotype=rep(c("No-0","NIT2"),each=6),salt_mM=rep(rep(c(0,75),each=3),2))
d<-data.frame(Protein.Group=c("X1","X2","X3","X4;X5","X6","X7"),Genes=letters[1:6],
  N.Sequences=rep(2,6),N.Proteotypic.Sequences=c(1,0,1,1,1,1))
v<-matrix(10,nrow=6,ncol=12);v[6,1:2]<-NA_real_
x<-list(data=d,values=v,metadata=meta,samples=paste0("s",1:12),spec=list(assay="proteomics",tissue="Root"))
saved<-list(data=data.frame(Protein.Group=d$Protein.Group,AGI=c("AT1G00010","AT1G00020","AT1G00030","AT1G00040","AT1G00030","AT1G00060")))
m<-integration_protein_mapping(x,saved,paste0("AT1G000",c("10","20","30","40","60")))
stopifnot(identical(which(m$eligible),1L),!m$eligible[2],all(m$duplicate_locus[c(3,5)]),
  !m$single_accession[4],!m$primary_observed[6])
cat("Integration aggregation, conjunction, multiplicity and mapping tests passed.\n")

# Ontology aliases must not lose ancestry or count a gene multiple times.
fake<-list(data=data.frame(locus=c("AT1G00010","AT1G00010","AT1G00020"),
  GO=c("GO:0015992","GO:0015991","GO:0009651")),spec=list())
resolved<-integration_canonical_associations(fake)
stopifnot(nrow(resolved$data)==2L,
  identical(resolved$data$GO,c("GO:1902600","GO:0009651")),
  all(integration_settings()$pathway_GO%in%AnnotationDbi::keys(GO.db::GO.db,keytype="GOID")))
cat("Secondary GO identifier resolution and de-duplication passed.\n")
