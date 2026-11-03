# Post-review diagnostics; these do not redefine the primary inference universe.
reviewer_loo_cor <- function(x,y) {
  stopifnot(length(x)==length(y),length(x)>3L,all(is.finite(c(x,y))))
  n<-length(x);dx<-x-mean(x);dy<-y-mean(y);scale<-n/(n-1)
  numerator<-sum(dx*dy)-scale*dx*dy
  denominator<-sqrt(pmax(0,sum(dx^2)-scale*dx^2)*pmax(0,sum(dy^2)-scale*dy^2))
  ifelse(denominator>0,numerator/denominator,NA_real_)
}

build_reviewer_influence <- function(inputs,settings,destination) {
  all<-integration_results(inputs)$table
  d<-all[all$method=="primary",]
  # Explicit post hoc identity set: the three loci highlighted by primary joint BH.
  highlighted<-sort(unique(d$locus[d$joint_q_global<.05]))
  indices<-split(seq_len(nrow(d)),interaction(d$comparison,d$tissue,d$contrast,drop=TRUE))
  rows<-lapply(indices,function(idx){
    z<-d[idx,];keep<-!z$locus%in%highlighted
    loo<-reviewer_loo_cor(z$RNA_log2FC,z$protein_log2FC)
    r<-cor(z$RNA_log2FC,z$protein_log2FC);influence<-abs(loo-r)
    j<-which.max(influence)
    data.frame(z[1,c("comparison","tissue","contrast")],genes=nrow(z),
      Pearson=r,Spearman=cor(z$RNA_log2FC,z$protein_log2FC,method="spearman"),
      omitted_loci=paste(highlighted,collapse=";"),genes_after_omission=sum(keep),
      Pearson_without_highlighted=cor(z$RNA_log2FC[keep],z$protein_log2FC[keep]),
      Spearman_without_highlighted=cor(z$RNA_log2FC[keep],z$protein_log2FC[keep],method="spearman"),
      leave_one_out_Pearson_min=min(loo,na.rm=TRUE),leave_one_out_Pearson_max=max(loo,na.rm=TRUE),
      most_influential_locus=z$locus[j],largest_absolute_r_change=influence[j],
      diagnostic="Post-review influence assessment; no hypothesis test or biological replicate exclusion")
  })
  write_table(do.call(rbind,rows),destination)
}

build_reviewer_joint_sensitivity <- function(inputs,settings,destination) {
  z<-integration_results(inputs);d<-z$table
  hit<-d[d$method=="primary"&d$joint_q_global<.05,]
  rows<-lapply(seq_len(nrow(hit)),function(i){
    h<-hit[i,];same<-d$comparison==h$comparison&d$tissue==h$tissue&d$contrast==h$contrast&d$locus==h$locus
    sensitivity<-d[same,]
    mapping<-z$mapping[[h$tissue]]
    m<-mapping[match(h$protein_group,mapping$protein_group),]
    value<-function(method,column){x<-sensitivity[sensitivity$method==method,column];if(length(x))x else NA_real_}
    data.frame(h[c("comparison","tissue","contrast","locus","protein_group","RNA_log2FC","protein_log2FC")],
      primary_global_BH=h$joint_q_global,primary_global_BY=h$joint_q_BY,
      complete_case_global_BH=value("complete_case","joint_q_global"),
      single_transcript_global_BH=value("single_transcript","joint_q_global"),
      all_observed=m$complete_observed,proteotypic_sequences=m$proteotypic_sequences,
      transcript_ID_count=h$transcript_ID_count,
      all_BH_sensitivities_below_0_05=all(sensitivity$joint_q_global<.05)&&nrow(sensitivity)==3L)
  })
  write_table(do.call(rbind,rows),destination)
}

build_reviewer_measurement_audit <- function(inputs,settings,destination) {
  z<-integration_results(inputs);h<-z$table[z$table$method=="primary"&z$table$joint_q_global<.05,]
  rows<-list()
  append_cells<-function(y,m,locus,assay,tissue,protein_group=NA_character_,support=NA_integer_,transcripts=NA_integer_){
    for(cell in levels(m$cell)){
      v<-as.numeric(y[m$cell==cell]);observed<-v[is.finite(v)]
      rows[[length(rows)+1L]]<<-data.frame(assay=assay,tissue=tissue,locus=locus,protein_group=protein_group,
        cell=cell,expected_pools=length(v),observed_pools=length(observed),
        minimum_log2=min(observed),mean_log2=mean(observed),maximum_log2=max(observed),
        proteotypic_sequences=support,transcript_ID_count=transcripts,
        scope="Supplied processed values and reported support; peptide or probe specificity not independently validated")
    }
  }
  a<-z$arrays$genotype
  for(locus in unique(h$locus))append_cells(a$y[locus,],a$metadata,locus,"RNA","Combined",
    transcripts=unname(a$transcripts[locus]))
  identities<-unique(h[c("tissue","locus","protein_group")])
  for(i in seq_len(nrow(identities))){
    id<-identities[i,];p<-z$proteins[[id$tissue]];m<-z$mapping[[id$tissue]]
    mapping<-m[match(id$protein_group,m$protein_group),]
    stopifnot(mapping$eligible,mapping$complete_observed)
    append_cells(p$y[id$protein_group,],p$metadata,id$locus,"Protein",id$tissue,id$protein_group,
      support=mapping$proteotypic_sequences)
  }
  # Independently recalculate displayed basal/saline differences from cell means.
  audit<-do.call(rbind,rows)
  means<-function(assay,tissue,locus,cell)audit$mean_log2[audit$assay==assay&audit$tissue==tissue&
    audit$locus==locus&audit$cell==cell]
  for(i in seq_len(nrow(h))){
    cells<-if(h$contrast[i]=="baseline_genotype")c("OE0","WT0") else if(h$contrast[i]=="saline_genotype")c("OES","WTS") else stop("Unexpected highlighted contrast")
    r<-means("RNA","Combined",h$locus[i],cells[1])-means("RNA","Combined",h$locus[i],cells[2])
    p<-means("Protein",h$tissue[i],h$locus[i],cells[1])-means("Protein",h$tissue[i],h$locus[i],cells[2])
    stopifnot(length(r)==1L,length(p)==1L,abs(r-h$RNA_log2FC[i])<1e-10,abs(p-h$protein_log2FC[i])<1e-10)
  }
  write_table(audit,destination)
}

build_reviewer_candidate_dose <- function(inputs,settings,destination) {
  z<-revision_dose(inputs)
  loci<-c(NIT2="AT3G44300",SOS1="AT2G01980",AHA2="AT4G30190",KUP4="AT4G23640")
  rows<-lapply(seq_along(loci),function(i){
    id<-revision_candidate_transcript(z$y,loci[i]);a<-z$joint[z$joint$feature==id,]
    data.frame(Gene=names(loci)[i],locus=unname(loci[i]),a,
      correction="BH across all 29391 transcript IDs separately for each joint test")
  })
  write_table(do.call(rbind,rows),destination)
}
