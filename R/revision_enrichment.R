# Reanalyse the supplied association resource with an experiment-specific
# eligible universe. Annotation completeness/version remains a limitation.
read_revision_associations <- function(path,spec) {
  d<-read.delim(path,header=FALSE,stringsAsFactors=FALSE)
  if(ncol(d)!=4L || any(!grepl("^AT[1-5CM]G[0-9]+$",d[[3]])) ||
     any(!grepl("^GO:[0-9]{7}$",d[[4]])))stop("Invalid external gene-GO associations")
  list(data=unique(data.frame(locus=d[[3]],GO=d[[4]])),spec=spec)
}

revision_go_sets <- function(association) revision_cached("GO_sets",association,{
  if(!requireNamespace("GO.db",quietly=TRUE))stop("Restore GO.db")
  a<-association$data
  terms<-unique(a$GO)
  # Propagate through all ontology ancestors using the pinned GO.db release.
  ancestors<-list(BP=as.list(GO.db::GOBPANCESTOR),MF=as.list(GO.db::GOMFANCESTOR),
    CC=as.list(GO.db::GOCCANCESTOR))
  expanded<-lapply(split(a$locus,a$GO),unique)
  for(go in terms){
    parents<-unique(unlist(lapply(ancestors,function(x)x[[go]])))
    for(p in parents[!is.na(parents)])expanded[[p]]<-unique(c(expanded[[p]],expanded[[go]]))
  }
  expanded
})

build_revision_enrichment <- function(inputs,settings,destination) {
  sets<-revision_go_sets(inputs$go_associations)
  rows<-list()
  for(t in c("root","shoot")){
    z<-revision_protein(inputs[[paste0("proteomics_",t)]])
    saved<-inputs[[paste0("saved_proteomics_",t)]]$data
    d<-z$annotations
    # Only unambiguous single-accession groups with one unique legacy AGI are
    # eligible for this gene-level analysis; no first accession substitution.
    idx<-match(d$Protein.Group,saved$Protein.Group)
    loci<-saved$AGI[idx]
    valid<-!grepl(";",d$Protein.Group,fixed=TRUE) &
      !is.na(d$N.Proteotypic.Sequences) & d$N.Proteotypic.Sequences>=1 & !is.na(loci) &
      grepl("^AT[1-5CM]G[0-9]+$",loci) & z$selection$observed_min2
    repeated<-duplicated(loci)|duplicated(loci,fromLast=TRUE)
    valid<-valid & !repeated
    mapping<-setNames(loci[valid],d$Protein.Group[valid])
    annotated<-unique(unlist(sets));universe<-intersect(unname(mapping),annotated)
    eligible_sets<-lapply(sets,intersect,y=universe)
    eligible_sets<-eligible_sets[lengths(eligible_sets)>=5 & lengths(eligible_sets)<=.8*length(universe)]
    for(co in c("baseline_genotype","saline_genotype","interaction"))for(direction in c("higher","lower")){
      tab<-z$table[z$table$method=="observed_min2" & z$table$contrast==co,,drop=FALSE]
      hit<-!is.na(tab$q)&tab$q<.05 & if(direction=="higher")tab$log2FC>0 else tab$log2FC<0
      query<-intersect(unname(mapping[tab$feature[hit]]),universe)
      for(go in names(eligible_sets)){
        members<-eligible_sets[[go]];k<-sum(query%in%members)
        p<-phyper(k-1,length(members),length(universe)-length(members),length(query),lower.tail=FALSE)
        rows[[length(rows)+1L]]<-data.frame(tissue=t,contrast=co,direction=direction,
          GO=go,query_hits=k,query_total=length(query),universe_hits=length(members),
          universe_total=length(universe),P=p,
          genes=paste(sort(intersect(query,members)),collapse=";"))
      }
    }
  }
  out<-do.call(rbind,rows)
  out$q_family<-ave(out$P,interaction(out$tissue,out$contrast,out$direction,drop=TRUE),FUN=function(x)p.adjust(x,"BH"))
  out$q_global<-p.adjust(out$P,"BH")
  term<-AnnotationDbi::select(GO.db::GO.db,keys=unique(out$GO),keytype="GOID",columns=c("TERM","ONTOLOGY"))
  index<-match(out$GO,term$GOID)
  out$term<-term$TERM[index];out$ontology<-term$ONTOLOGY[index]
  write_table(out,destination)
}
