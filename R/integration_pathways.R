integration_pathway_inputs <- function(z,comparison,tissue,sets) {
  mapping<-z$mapping[[tissue]]
  mapping<-mapping[mapping$eligible&mapping$complete_observed,,drop=FALSE]
  array<-z$arrays[[if(comparison=="No0_genotype")"genotype" else paste0("Col0_",tissue)]]
  mapping<-mapping[mapping$locus%in%rownames(array$y),]
  p<-z$proteins[[tissue]]
  py<-p$y[match(mapping$protein_group,rownames(p$y)),,drop=FALSE]
  rownames(py)<-mapping$locus
  ay<-array$y[mapping$locus,,drop=FALSE]
  stopifnot(identical(rownames(ay),rownames(py)),all(is.finite(py)),all(is.finite(ay)))
  pm<-p$metadata;pd<-model.matrix(~0+cell,pm);colnames(pd)<-levels(pm$cell)
  list(RNA=ay,protein=py,array=array,protein_design=pd,protein_contrasts=factorial_contrasts(pd),
    sets=lapply(sets,intersect,y=mapping$locus))
}

integration_pathways <- function(inputs) integration_cached("focused_pathways",inputs,{
  z<-integration_results(inputs);settings<-integration_settings()
  allsets<-revision_go_sets(inputs$go_associations)
  sets<-setNames(lapply(settings$pathway_GO,function(go)allsets[[go]]),settings$pathway_GO)
  labels<-AnnotationDbi::select(GO.db::GO.db,keys=names(sets),columns=c("TERM","ONTOLOGY"),keytype="GOID")
  rows<-memberships<-list()
  for(tissue in c("Root","Shoot"))for(comparison in c("No0_genotype","Col0_context")){
    a<-integration_pathway_inputs(z,comparison,tissue,sets)
    n<-lengths(a$sets);eligible<-n>=settings$pathway_min_members & n<=settings$pathway_max_fraction*nrow(a$RNA)
    indices<-lapply(a$sets[eligible],function(m)match(m,rownames(a$RNA)))
    memberships[[length(memberships)+1L]]<-data.frame(comparison=comparison,tissue=tissue,
      set="eligible_universe",locus=rownames(a$RNA))
    for(go in names(a$sets))if(length(a$sets[[go]]))memberships[[length(memberships)+1L]]<-
      data.frame(comparison=comparison,tissue=tissue,set=go,locus=a$sets[[go]])
    contrasts<-if(comparison=="No0_genotype")settings$contrasts else "salt_WT"
    for(co in contrasts)for(mode in c("estimated_correlation","fixed_0_01")){
      correlation<-if(mode=="estimated_correlation")NA_real_ else .01
      if(length(indices)){
        ra<-limma::camera(a$RNA,index=indices,design=a$array$design,contrast=a$array$contrasts[,co],
          inter.gene.cor=correlation,allow.neg.cor=FALSE,trend.var=TRUE,sort=FALSE)
        pr<-limma::camera(a$protein,index=indices,design=a$protein_design,contrast=a$protein_contrasts[,co],
          inter.gene.cor=correlation,allow.neg.cor=FALSE,trend.var=TRUE,sort=FALSE)
      }
      for(go in names(sets)){
        ok<-eligible[[go]];r<-if(ok)ra[go,] else NULL;p<-if(ok)pr[go,] else NULL
        rr<-if(ok)ifelse(r$Direction=="Up",1,-1) else NA_real_
        pp<-if(ok)ifelse(p$Direction=="Up",1,-1) else NA_real_
        rows[[length(rows)+1L]]<-data.frame(method=mode,comparison=comparison,tissue=tissue,contrast=co,
          GO=go,term=labels$TERM[match(go,labels$GOID)],members=n[[go]],universe=nrow(a$RNA),
          RNA_direction=if(ok)r$Direction else NA_character_,protein_direction=if(ok)p$Direction else NA_character_,
          RNA_P=if(ok)r$PValue else NA_real_,protein_P=if(ok)p$PValue else NA_real_,
          RNA_correlation=if(ok && mode=="estimated_correlation")r$Correlation else correlation,
          protein_correlation=if(ok && mode=="estimated_correlation")p$Correlation else correlation,
          joint_P=if(ok)integration_conjunction(r$PValue,rr,p$PValue,pp) else NA_real_,
          status=if(ok)"tested" else "fewer than five measured members or more than 80 percent of universe")
      }
    }
  }
  d<-integration_adjust(do.call(rbind,rows))
  list(table=d,memberships=do.call(rbind,memberships))
})

build_integration_pathways <- function(inputs,settings,destination)
  write_table(integration_pathways(inputs)$table,destination)
build_integration_memberships <- function(inputs,settings,destination)
  write_table(integration_pathways(inputs)$memberships,destination)
