build_revision_primary_summary <- function(inputs,settings,destination) {
  rows<-list()
  z<-revision_discoveries(revision_genotype(inputs$microarray_genotype)$table,"contrast")
  z$tissue<-"Combined";z$method<-"limma_genotype";rows[[1]]<-z
  for(x in inputs[c("proteomics_root","proteomics_shoot")]){
    d<-revision_protein(x)$table;d<-d[d$method=="observed_min2",]
    rows[[length(rows)+1L]]<-revision_discoveries(d,c("contrast","tissue","method"))
  }
  columns<-c("tissue","method","contrast","tested","q_lt_0_05","positive","negative","global_q_lt_0_05")
  d<-do.call(rbind,lapply(rows,function(x)x[columns]))
  d<-d[d$contrast%in%c("baseline_genotype","salt_WT","salt_OE","saline_genotype","interaction"),]
  write_table(d,destination)
}

build_revision_primary_candidates <- function(inputs,settings,destination) {
  rows<-list()
  for(x in inputs){z<-revision_protein(x)
    for(g in c("NIT2","NHX7","AHA2"))for(co in
      c("baseline_genotype","saline_genotype","salt_WT","interaction")){
      a<-z$table[z$table$method=="observed_min2"&z$table$Genes==g&z$table$contrast==co,]
      gene<-if(g=="NHX7")"SOS1" else g
      rows[[length(rows)+1L]]<-data.frame(Gene=gene,tissue=x$spec$tissue,contrast=co,
        log2FC=if(nrow(a))a$log2FC else NA_real_,
        lower95=if(nrow(a))a$lower95 else NA_real_,upper95=if(nrow(a))a$upper95 else NA_real_,
        q=if(nrow(a))a$q else NA_real_,q_all_contrasts=if(nrow(a))a$q_all_contrasts else NA_real_,
        status=if(nrow(a))"primary observed-data model" else "ineligible: fewer than two observations in at least one cell",
        identity=if(g=="AHA2")"AHA2-labelled group; zero reported proteotypic sequences" else "reported gene label")
    }
  }
  write_table(do.call(rbind,rows),destination)
}

build_revision_similarity_summary <- function(inputs,settings,destination) {
  rows<-list()
  a<-revision_genotype(inputs$microarray_genotype)
  sets<-list(Genotype_microarray=list(y=a$y,m=a$metadata))
  for(x in inputs[c("proteomics_root","proteomics_shoot")]){
    a<-revision_protein(x)
    sets[[paste0(x$spec$tissue,"_proteome_complete_cases")]]<-list(
      y=a$y[a$selection$complete_case,,drop=FALSE],m=a$metadata)
  }
  for(name in names(sets)){
    z<-prestress_splits(sets[[name]]$y,sets[[name]]$m,name)
    a<-z[z$comparison=="shared_controls",];b<-z[z$comparison=="disjoint_controls",]
    rows[[length(rows)+1L]]<-data.frame(dataset=name,features=a$features,
      shared_Pearson=a$Pearson,disjoint_Pearson_min=min(b$Pearson),disjoint_Pearson_max=max(b$Pearson),
      splits=nrow(b))
  }
  write_table(do.call(rbind,rows),destination)
}

build_revision_sensitivity_summary <- function(inputs,settings,destination) {
  d<-do.call(rbind,lapply(inputs,function(x)revision_protein(x)$table))
  d<-d[d$contrast%in%c("baseline_genotype","interaction"),]
  write_table(revision_discoveries(d,c("contrast","tissue","method")),destination)
}
