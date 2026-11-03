integration_settings <- function() list(
  contrasts=c("baseline_genotype","salt_WT","salt_OE","saline_genotype","interaction"),
  methods=c("primary","complete_case","single_transcript"),
  primary_protein_method="observed_min2", locus_summary="mean of unique transcript log2 signals",
  joint_threshold=.05, pathway_min_members=5L, pathway_max_fraction=.8,
  pathway_GO=c("GO:0009651","GO:0006814","GO:0006813","GO:1902600","GO:0009733",
    "GO:0000302","GO:0006979","GO:0042744","GO:0016209"),
  candidate_loci=c(NIT2="AT3G44300",SOS1="AT2G01980",AHA2="AT4G30190",
    HKT1="AT4G10310",HAK5="AT4G13420",KUP4="AT4G23640"))

integration_cached <- function(key,inputs,expression) {
  files<-file.path(project_root,c("R/integration_models.R","R/integration_pathways.R"))
  fingerprint<-digest::digest(vapply(files,function(p)digest::digest(file=p,algo="sha256"),""),algo="sha256")
  revision_cached(paste("integration",key,fingerprint,sep=":"),inputs,expression)
}

integration_loci <- function(x) {
  y<-collapse_microarray(x)
  locus<-sub("\\.[0-9]+$","",rownames(y))
  n<-table(locus)
  value<-rowsum(y,locus,reorder=FALSE)
  value<-value/as.numeric(n[rownames(value)])
  stopifnot(all(is.finite(value)),!anyDuplicated(rownames(value)))
  list(y=value,transcripts=setNames(as.integer(n),names(n)))
}

integration_array <- function(x) integration_cached(paste0("array_",x$spec$tissue),x,{
  loci<-integration_loci(x);y<-loci$y
  if(identical(x$spec$tissue,"Combined")){
    m<-factorial_metadata(x)
    design<-model.matrix(~0+cell,m);colnames(design)<-levels(m$cell)
    cm<-factorial_contrasts(design)
  }else{
    m<-data.frame(sample=x$samples,
      dose=factor(ifelse(startsWith(x$samples,"PN"),"0",sub(" .*","",x$samples)),
        levels=c("0","50","75","100")))
    design<-model.matrix(~0+dose,m);colnames(design)<-paste0("dose",levels(m$dose))
    cm<-limma::makeContrasts(salt_WT=dose75-dose0,levels=design)
  }
  fit<-moderated_contrasts(y,design,cm)
  list(y=y,metadata=m,design=design,contrasts=cm,table=fit$table,transcripts=loci$transcripts)
})

integration_protein_mapping <- function(x,saved,array_loci) {
  d<-x$data;s<-saved$data
  assert_columns(s,c("Protein.Group","AGI"),"integration protein annotation")
  if(anyDuplicated(s$Protein.Group))stop("Duplicate annotation protein-group keys")
  loci<-trimws(s$AGI[match(d$Protein.Group,s$Protein.Group)])
  retained<-proteomics_keep(x,analysis_settings)
  m<-factorial_metadata(x);counts<-protein_cell_counts(x$values,m)
  valid_locus<-!is.na(loci)&grepl("^AT[1-5CM]G[0-9]+$",loci)
  duplicated_locus<-rep(FALSE,nrow(d))
  k<-which(retained&valid_locus)
  duplicated_locus[k]<-duplicated(loci[k])|duplicated(loci[k],fromLast=TRUE)
  single_accession<-!grepl(";",d$Protein.Group,fixed=TRUE)
  unique_support<-!is.na(d$N.Proteotypic.Sequences)&d$N.Proteotypic.Sequences>=1
  primary<-retained&apply(counts,1,min)>=2
  complete<-retained&rowSums(counts)==12
  shared<-valid_locus&loci%in%array_loci
  eligible<-primary&single_accession&unique_support&valid_locus&!duplicated_locus&shared
  status<-vapply(seq_len(nrow(d)),function(i){
    reasons<-c(if(!retained[i])"fails original retention",if(!primary[i])"fewer than two observed pools in a cell or insufficient sequence support",
      if(!single_accession[i])"multiple accessions",if(!unique_support[i])"no reported proteotypic sequence",
      if(!valid_locus[i])"missing or invalid AGI",if(duplicated_locus[i])"duplicate locus among retained groups",
      if(!shared[i])"no corresponding array locus")
    if(length(reasons))paste(reasons,collapse="; ") else "eligible"
  },"")
  data.frame(tissue=x$spec$tissue,protein_group=d$Protein.Group,Gene_label=d$Genes,locus=loci,
    retained=retained,single_accession=single_accession,proteotypic_sequences=d$N.Proteotypic.Sequences,
    valid_AGI=valid_locus,duplicate_locus=duplicated_locus,array_covered=shared,
    primary_observed=primary,complete_observed=complete,eligible=eligible,status=status,
    stringsAsFactors=FALSE)
}

# Conservative directional conjunction; neither assay may dominate the other.
integration_conjunction <- function(p1,b1,p2,b2) {
  positive1<-ifelse(b1>=0,p1/2,1-p1/2);positive2<-ifelse(b2>=0,p2/2,1-p2/2)
  negative1<-ifelse(b1<=0,p1/2,1-p1/2);negative2<-ifelse(b2<=0,p2/2,1-p2/2)
  pmin(1,2*pmin(pmax(positive1,positive2),pmax(negative1,negative2)))
}

integration_adjust <- function(d,p="joint_P",sensitivity="method") {
  d$joint_q_family<-ave(d[[p]],interaction(d[[sensitivity]],d$comparison,d$tissue,d$contrast,drop=TRUE),
    FUN=function(x)p.adjust(x,"BH"))
  d$joint_q_global<-ave(d[[p]],d[[sensitivity]],FUN=function(x)p.adjust(x,"BH"))
  d$joint_q_BY<-ave(d[[p]],d[[sensitivity]],FUN=function(x)p.adjust(x,"BY"))
  d
}

integration_join <- function(array,protein,mapping,contrast,tissue,comparison,method) {
  eligible<-mapping$eligible
  if(method=="complete_case")eligible<-eligible&mapping$complete_observed
  if(method=="single_transcript")eligible<-eligible&array$transcripts[mapping$locus]==1L
  eligible[is.na(eligible)]<-FALSE
  map<-mapping[eligible,c("protein_group","locus"),drop=FALSE]
  pmethod<-if(method=="complete_case")"complete_case" else "observed_min2"
  p<-protein$table[protein$table$method==pmethod&protein$table$contrast==contrast,,drop=FALSE]
  p<-p[match(map$protein_group,p$feature),,drop=FALSE]
  a<-array$table[array$table$contrast==contrast,,drop=FALSE]
  a<-a[match(map$locus,a$feature),,drop=FALSE]
  stopifnot(!anyNA(a$feature),!anyNA(p$feature))
  out<-data.frame(method=method,comparison=comparison,tissue=tissue,contrast=contrast,
    locus=map$locus,protein_group=map$protein_group,
    transcript_ID_count=unname(array$transcripts[map$locus]),
    RNA_log2FC=a$log2FC,RNA_lower95=a$lower95,RNA_upper95=a$upper95,RNA_P=a$P,RNA_q=a$q,
    protein_log2FC=p$log2FC,protein_lower95=p$lower95,protein_upper95=p$upper95,protein_P=p$P,protein_q=p$q,
    stringsAsFactors=FALSE)
  out$same_direction<-sign(out$RNA_log2FC)==sign(out$protein_log2FC)
  out$joint_P<-integration_conjunction(out$RNA_P,out$RNA_log2FC,out$protein_P,out$protein_log2FC)
  out
}

integration_results <- function(inputs) {
  base<-c("microarray_genotype","microarray_root","microarray_shoot","proteomics_root",
    "proteomics_shoot","saved_proteomics_root","saved_proteomics_shoot")
  integration_cached("all_gene_pairs",inputs[base],{
    arrays<-list(genotype=integration_array(inputs$microarray_genotype),
      Col0_Root=integration_array(inputs$microarray_root),Col0_Shoot=integration_array(inputs$microarray_shoot))
    proteins<-maps<-list();rows<-list()
    for(t in c("Root","Shoot")){
      id<-paste0("proteomics_",tolower(t));p<-revision_protein(inputs[[id]])
      mapping<-integration_protein_mapping(inputs[[id]],inputs[[paste0("saved_",id)]],rownames(arrays$genotype$y))
      proteins[[t]]<-p;maps[[t]]<-mapping
      for(method in integration_settings()$methods){
        for(co in integration_settings()$contrasts)
          rows[[length(rows)+1L]]<-integration_join(arrays$genotype,p,mapping,co,t,"No0_genotype",method)
        rows[[length(rows)+1L]]<-integration_join(arrays[[paste0("Col0_",t)]],p,mapping,
          "salt_WT",t,"Col0_context",method)
      }
    }
    d<-integration_adjust(do.call(rbind,rows))
    list(table=d,arrays=arrays,proteins=proteins,mapping=maps)
  })
}

build_integration_genes <- function(inputs,settings,destination)
  write_table(integration_results(inputs)$table,destination)

build_integration_mapping <- function(inputs,settings,destination)
  write_table(do.call(rbind,integration_results(inputs)$mapping),destination)

integration_summary <- function(d) {
  indices<-split(seq_len(nrow(d)),interaction(d$method,d$comparison,d$tissue,d$contrast,drop=TRUE))
  do.call(rbind,lapply(indices,function(idx){
    z<-d[idx,];s<-z[1,c("method","comparison","tissue","contrast")]
    s$genes<-nrow(z);s$Pearson<-cor(z$RNA_log2FC,z$protein_log2FC)
    s$Spearman<-cor(z$RNA_log2FC,z$protein_log2FC,method="spearman")
    without_transgene<-z$locus!=unname(integration_settings()$candidate_loci["NIT2"])
    s$Pearson_without_NIT2<-cor(z$RNA_log2FC[without_transgene],z$protein_log2FC[without_transgene])
    s$Spearman_without_NIT2<-cor(z$RNA_log2FC[without_transgene],z$protein_log2FC[without_transgene],method="spearman")
    s$direction_agreement<-mean(z$same_direction)
    s$both_assay_q<-sum(z$RNA_q<.05&z$protein_q<.05)
    s$opposite_both_assay_q<-sum(z$RNA_q<.05&z$protein_q<.05&!z$same_direction)
    s$joint_family_BH<-sum(z$joint_q_family<.05)
    s$joint_global_BH<-sum(z$joint_q_global<.05)
    s$joint_global_BY<-sum(z$joint_q_BY<.05)
    s
  }))
}

build_integration_summary <- function(inputs,settings,destination)
  write_table(integration_summary(integration_results(inputs)$table),destination)

build_integration_candidates <- function(inputs,settings,destination) {
  z<-integration_results(inputs);d<-z$table[z$table$method=="primary",]
  settings<-integration_settings();rows<-list()
  families<-unique(d[c("comparison","tissue","contrast")])
  for(i in seq_len(nrow(families)))for(g in names(settings$candidate_loci)){
    f<-families[i,];locus<-unname(settings$candidate_loci[g])
    hit<-d[d$comparison==f$comparison&d$tissue==f$tissue&d$contrast==f$contrast&d$locus==locus,]
    mapping<-z$mapping[[f$tissue]]
    notes<-unique(mapping$status[!is.na(mapping$locus)&mapping$locus==locus])
    a<-z$arrays[[if(f$comparison=="No0_genotype")"genotype" else paste0("Col0_",f$tissue)]]$table
    a<-a[a$feature==locus&a$contrast==f$contrast,]
    rows[[length(rows)+1L]]<-data.frame(f,Gene=g,locus=locus,
      RNA_log2FC=if(nrow(a))a$log2FC else NA_real_,RNA_q=if(nrow(a))a$q else NA_real_,
      protein_log2FC=if(nrow(hit))hit$protein_log2FC else NA_real_,protein_q=if(nrow(hit))hit$protein_q else NA_real_,
      joint_q_global=if(nrow(hit))hit$joint_q_global else NA_real_,
      status=if(nrow(hit))"eligible" else if(length(notes))paste(notes,collapse=" | ") else "not mapped to a supplied protein group")
  }
  write_table(do.call(rbind,rows),destination)
}
