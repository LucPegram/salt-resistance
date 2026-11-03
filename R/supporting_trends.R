# Post hoc prioritization of significant, concordant transcript/protein trends.
# Preserve all original model estimates, testing universes and adjusted P values.
.supporting_cache <- new.env(parent=emptyenv())
supporting_settings <- function() list(assay_BH=.05,profile_Pearson_min=0,
  top_per_tissue=10L,profile_genes=3L,profile_cells=c('WT0','WTS','OE0','OES'),
  ranking='ascending maximum assay q; descending profile r; descending smaller absolute effect',
  scope='post hoc candidate screen; original joint global BH/BY retained')

supporting_cells <- function(y,m) {
  cells<-supporting_settings()$profile_cells
  stopifnot(all(cells%in%as.character(m$cell)))
  out<-vapply(cells,function(cell)rowMeans(y[,m$cell==cell,drop=FALSE],na.rm=TRUE),numeric(nrow(y)))
  rownames(out)<-rownames(y);out
}

supporting_profile <- function(rna,protein) {
  stopifnot(identical(dim(rna),dim(protein)),ncol(rna)==4L)
  data.frame(profile_Pearson=vapply(seq_len(nrow(rna)),function(i){
    if(!all(is.finite(c(rna[i,],protein[i,])))||sd(rna[i,])==0||sd(protein[i,])==0)NA_real_ else cor(rna[i,],protein[i,])},0.),
    profile_Spearman=vapply(seq_len(nrow(rna)),function(i){
      if(!all(is.finite(c(rna[i,],protein[i,])))||sd(rna[i,])==0||sd(protein[i,])==0)NA_real_ else cor(rna[i,],protein[i,],method='spearman')},0.))
}

supporting_flag <- function(d) {
  finite<-is.finite(d$RNA_q)&is.finite(d$protein_q)&is.finite(d$RNA_log2FC)&is.finite(d$protein_log2FC)
  d$assaywise_supported<-finite&d$RNA_q<supporting_settings()$assay_BH&
    d$protein_q<supporting_settings()$assay_BH&d$RNA_log2FC!=0&d$protein_log2FC!=0&sign(d$RNA_log2FC)==sign(d$protein_log2FC)
  d$profile_positive<-is.finite(d$profile_Pearson)&d$profile_Pearson>supporting_settings()$profile_Pearson_min
  d$selected_support<-d$assaywise_supported&d$profile_positive&d$comparison=='No0_genotype'
  d$balanced_assay_q<-pmax(d$RNA_q,d$protein_q)
  d$smaller_abs_effect<-pmin(abs(d$RNA_log2FC),abs(d$protein_log2FC))
  d$direction<-ifelse(d$RNA_log2FC>0,'Up','Down')
  d$evidence_level<-ifelse(d$joint_q_global<.05,'joint global BH',
    ifelse(d$joint_q_family<.05,'joint family BH','assay-wise BH only'))
  d
}

supporting_order <- function(d) {
  d[order(d$balanced_assay_q,-d$profile_Pearson,-d$smaller_abs_effect,d$locus,d$tissue,d$contrast),,drop=FALSE]
}

supporting_dual <- function(d) {
  rows<-list()
  for(idx in split(seq_len(nrow(d)),paste(d$tissue,d$locus))){
    z<-d[idx,,drop=FALSE]
    salts<-which(z$contrast%in%c('salt_WT','salt_OE'))
    genotypes<-which(z$contrast%in%c('baseline_genotype','saline_genotype'))
    if(!length(salts)||!length(genotypes))next
    pairs<-expand.grid(salt=salts,genotype=genotypes)
    pairs$worst_q<-pmax(z$balanced_assay_q[pairs$salt],z$balanced_assay_q[pairs$genotype])
    pairs<-pairs[order(pairs$worst_q,z$contrast[pairs$salt],z$contrast[pairs$genotype]),]
    salt<-z[pairs$salt[1],,drop=FALSE];genotype<-z[pairs$genotype[1],,drop=FALSE]
    row<-salt[,intersect(c('Gene_name','Aliases','Description','Name_source','locus','protein_group','tissue',
      'proteotypic_sequences','transcript_ID_count','profile_Pearson','profile_Spearman'),names(salt)),drop=FALSE]
    row$salt_contrast<-salt$contrast;row$genotype_contrast<-genotype$contrast
    fields<-intersect(c('RNA_log2FC','RNA_lower95','RNA_upper95','RNA_P','RNA_q',
      'protein_log2FC','protein_lower95','protein_upper95','protein_P','protein_q',
      'joint_q_family','joint_q_global','joint_q_BY','complete_case_supported',
      'complete_case_joint_global_BH','single_transcript_supported','single_transcript_joint_global_BH'),names(salt))
    for(field in fields){row[[paste0('salt_',field)]]<-salt[[field]];row[[paste0('genotype_',field)]]<-genotype[[field]]}
    row$worst_assay_q<-pairs$worst_q[1]
    row$both_contrasts_joint_global_BH<-salt$joint_q_global<.05&genotype$joint_q_global<.05
    rows[[length(rows)+1L]]<-row
  }
  if(!length(rows))return(data.frame())
  d<-do.call(rbind,rows);d<-d[order(d$worst_assay_q,-d$profile_Pearson,d$locus,d$tissue),]
  d$rank<-seq_len(nrow(d));d
}

supporting_three_datasets <- function(d) {
  a<-d[d$selected_support&d$contrast=='salt_WT',,drop=FALSE]
  b<-d[d$comparison=='Col0_context'&d$assaywise_supported,,drop=FALSE]
  key<-function(x)paste(x$tissue,x$locus,x$protein_group)
  j<-match(key(a),key(b));keep<-!is.na(j)
  a<-a[keep,,drop=FALSE];b<-b[j[keep],,drop=FALSE]
  agree<-sign(a$RNA_log2FC)==sign(b$RNA_log2FC)
  a<-a[agree,,drop=FALSE];b<-b[agree,,drop=FALSE]
  stopifnot(identical(a$protein_group,b$protein_group),all(abs(a$protein_log2FC-b$protein_log2FC)<1e-10))
  out<-a[,c('Gene_name','Aliases','Description','Name_source','locus','protein_group','tissue',
    'proteotypic_sequences','profile_Pearson','profile_Spearman'),drop=FALSE]
  out$Array1_log2FC<-b$RNA_log2FC;out$Array1_q<-b$RNA_q
  out$Array2_log2FC<-a$RNA_log2FC;out$Array2_q<-a$RNA_q
  out$protein_log2FC<-a$protein_log2FC;out$protein_q<-a$protein_q
  out$Array2_joint_q_global<-a$joint_q_global;out$Array1_joint_q_global<-b$joint_q_global
  out$worst_assay_q<-pmax(out$Array1_q,out$Array2_q,out$protein_q)
  out<-out[order(out$worst_assay_q,-out$profile_Pearson,out$locus,out$tissue),,drop=FALSE]
  out$rank<-seq_len(nrow(out));out
}

supporting_results <- function(inputs) {
  key<-digest::digest(list(inputs,supporting_settings()),algo='sha256')
  if(exists(key,.supporting_cache,inherits=FALSE))return(get(key,.supporting_cache,inherits=FALSE))
  z<-integration_results(inputs);d<-z$table[z$table$method=='primary',,drop=FALSE]
  profiles<-list();array_means<-supporting_cells(z$arrays$genotype$y,z$arrays$genotype$metadata)
  for(t in c('Root','Shoot')){
    map<-z$mapping[[t]];map<-map[map$eligible,]
    protein_means<-supporting_cells(z$proteins[[t]]$y,z$proteins[[t]]$metadata)
    metrics<-supporting_profile(array_means[map$locus,,drop=FALSE],protein_means[map$protein_group,,drop=FALSE])
    profiles[[t]]<-cbind(map[,c('locus','protein_group','tissue','Gene_label','proteotypic_sequences'),drop=FALSE],metrics)
  }
  profile<-do.call(rbind,profiles)
  joined<-match(paste(d$tissue,d$locus,d$protein_group),paste(profile$tissue,profile$locus,profile$protein_group))
  stopifnot(!anyNA(joined))
  for(field in c('profile_Pearson','profile_Spearman','Gene_label','proteotypic_sequences'))d[[field]]<-profile[[field]][joined]
  d<-toprank_add_names(d,toprank_annotation(inputs$gene_annotation))
  # Prefer a meaningful supplied protein gene symbol when the array has none.
  fallback<-!grepl('symbol$',d$Name_source)&!is.na(d$Gene_label)&nzchar(d$Gene_label)&
    !grepl('^(AT[1-5CM]G[0-9]+|[FTM][A-Z0-9]+[._][0-9]+)$',d$Gene_label,ignore.case=TRUE)
  d$Gene_name[fallback]<-d$Gene_label[fallback]
  d$Name_source[fallback]<-'supplied protein gene label; no array symbol supplied'
  primary_key<-paste(d$comparison,d$tissue,d$contrast,d$locus)
  for(method in c('complete_case','single_transcript')){
    s<-z$table[z$table$method==method,];j<-match(primary_key,paste(s$comparison,s$tissue,s$contrast,s$locus))
    d[[paste0(method,'_supported')]]<-s$RNA_q[j]<.05&s$protein_q[j]<.05&s$same_direction[j]
    d[[paste0(method,'_joint_global_BH')]]<-s$joint_q_global[j]
  }
  d<-supporting_flag(d)
  ranked<-supporting_order(d[d$selected_support,,drop=FALSE])
  ranked$rank_within_comparison<-ave(seq_len(nrow(ranked)),interaction(ranked$tissue,ranked$contrast,drop=TRUE),FUN=seq_along)
  dual<-supporting_dual(ranked)
  triple<-supporting_three_datasets(d)
  summary<-do.call(rbind,lapply(split(seq_len(nrow(d)),interaction(d$comparison,d$tissue,d$contrast,drop=TRUE)),function(idx){
    a<-d[idx,,drop=FALSE];out<-a[1,c('comparison','tissue','contrast')]
    out$eligible<-nrow(a);out$both_assay_significant_concordant<-sum(a$assaywise_supported)
    out$positive_profile_supported<-if(a$comparison[1]=='No0_genotype')sum(a$selected_support) else NA_integer_
    out$joint_family_BH<-sum(a$joint_q_family<.05);out$joint_global_BH<-sum(a$joint_q_global<.05)
    out$joint_global_BY<-sum(a$joint_q_BY<.05);out
  }))
  result<-list(audit=d,ranked=ranked,dual=dual,triple=triple,summary=summary,model=z)
  assign(key,result,.supporting_cache);result
}

supporting_signal <- function(result) {
  selected<-head(result$dual,supporting_settings()$profile_genes);rows<-list()
  for(i in seq_len(nrow(selected))){
    a<-selected[i,,drop=FALSE]
    for(assay in c('Array 2','Protein')){
      z<-if(assay=='Array 2')result$model$arrays$genotype else result$model$proteins[[a$tissue]]
      id<-if(assay=='Array 2')a$locus else a$protein_group
      rows[[length(rows)+1L]]<-data.frame(Gene_name=a$Gene_name,locus=a$locus,protein_group=a$protein_group,
        protein_tissue=a$tissue,assay=assay,sample=z$metadata$sample,cell=as.character(z$metadata$cell),
        genotype=z$metadata$genotype,salt_mM=z$metadata$salt_mM,
        log2_signal=as.numeric(z$y[id,]),observed=is.finite(z$y[id,]),stringsAsFactors=FALSE)
    }
  }
  if(length(rows))do.call(rbind,rows) else data.frame()
}

supporting_label <- function(x) unname(c(baseline_genotype='Basal genotype',saline_genotype='Genotype in salt',
  salt_WT='No-0 salt',salt_OE='NIT2 salt',interaction='Genotype x salt')[x])

supporting_page <- function(title,subtitle='') {
  par(mar=c(0,0,0,0));plot.new();plot.window(xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i')
  text(.035,.953,title,adj=0,font=2,cex=1.45,col=publication_colors['ink'])
  if(nzchar(subtitle))text(.035,.913,subtitle,adj=0,cex=.84,col=publication_colors['muted'])
}

supporting_text <- function(text,y,cex=.95,width=140) {
  lines<-strwrap(text,width=width)
  for(i in seq_along(lines))graphics::text(.035,y-(i-1)*.027,lines[i],adj=0,cex=cex,col=publication_colors['ink'])
  invisible(y-length(lines)*.027)
}

supporting_print_table <- function(d,fields,headers,x,y=.81,step=.031,cex=.78) {
  rect(.035,y-.019,.965,y+.019,col='#E8EEF2',border=NA)
  for(j in seq_along(x))text(x[j],y,headers[j],adj=if(j==1)0 else .5,cex=cex,font=2)
  for(i in seq_len(nrow(d))){
    yy<-y-i*step
    if(i%%2==0)rect(.035,yy-step*.48,.965,yy+step*.48,col='#F5F7F9',border=NA)
    for(j in seq_along(fields)){
      value<-d[[fields[j]]][i]
      if(is.numeric(value))value<-if(is.na(value))'NA' else formatC(value,digits=3,format='g')
      value<-as.character(value)
      if(j==1&&nchar(value)>37){
        prefix<-if(fields[j]=='Gene_name'&&'locus'%in%names(d))paste0(d$locus[i],' | ') else ''
        value<-paste0(prefix,substr(value,1,34-nchar(prefix)),'...')
      }
      text(x[j],yy,value,adj=if(j==1)0 else .5,cex=cex)
    }
  }
  invisible(y-(nrow(d)+1)*step)
}

supporting_report <- function(result,destination) {
  z<-result;d<-z$ranked
  grDevices::cairo_pdf(destination,width=11.7,height=8.3,family='sans',pointsize=10)
  on.exit(grDevices::dev.off(),add=TRUE)
  supporting_page('Genes with supporting transcript and protein trends','Candidate prioritization from Array 2 and tissue-resolved proteomics; Array 1 provides additional salt context')
  salt<-d[d$contrast%in%c('salt_WT','salt_OE'),];genotype<-d[d$contrast%in%c('baseline_genotype','saline_genotype'),]
  supporting_text(paste(nrow(d),'supported gene-by-tissue-by-comparison results represent',length(unique(d$locus)),
    'loci:',nrow(salt),'salt responses and',nrow(genotype),'genotype differences.',nrow(z$dual),'gene/tissue candidates below support both categories in the same protein tissue.'),.85)
  supporting_text('Screen: same-direction effects, original within-assay BH q < 0.05 in both assays, and positive Pearson correlation of the four condition means. Correlation uses No-0 control, No-0 salt, NIT2 control and NIT2 salt; it is descriptive, with no correlation P value and no sample pairing between assays.',.73)
  dual<-z$dual
  if(nrow(dual)){
    dual$display<-paste(ifelse(nchar(dual$Gene_name)>20,paste0(substr(dual$Gene_name,1,17),'...'),dual$Gene_name),dual$locus,sep=' | ')
    supporting_print_table(dual,c('display','tissue','salt_RNA_log2FC','salt_protein_log2FC','genotype_RNA_log2FC','genotype_protein_log2FC','worst_assay_q','profile_Pearson'),
      c('Gene / AGI','Tissue','Salt RNA','Salt protein','Geno RNA','Geno protein','Worst q','Profile r'),
      c(.045,.39,.49,.59,.69,.79,.87,.94),y=.56,step=.047,cex=.78)
  }
  supporting_text('Values are log2 differences. Each selected gene has significant, concordant changes in both assays in one salt and one genotype comparison. The selected contrasts for each row are listed below; full effects, pointwise intervals and q values are in supporting_dual.csv.',.32)
  supporting_text(paste(vapply(seq_len(nrow(dual)),function(i)paste0(dual$locus[i],': ',supporting_label(dual$salt_contrast[i]),' + ',supporting_label(dual$genotype_contrast[i])),''),collapse='; '),.225,cex=.78,width=175)
  supporting_text('Rank by the largest assay q (smaller is stronger), then profile correlation and the smaller absolute effect. These are post hoc supporting candidates, not a new joint-FDR discovery set. Original directional-conjunction family/global BH and global BY results remain in every audit; the stricter global-joint results are shown on page 5.',.155)
  supporting_text('Array 2 uses combined tissues at 0/100 mM; proteins use separate tissues at 0/75 mM. Agreement is useful evidence for prioritization. It does not require equal signal scales or a resolved explanation for other discordant genes.',.055,cex=.85)

  supporting_page('Top supporting salt responses','Ten per protein tissue, ordered by the weaker of the two assay q values; all are No-0 salt versus control')
  for(t in c('Root','Shoot')){
    a<-head(salt[salt$tissue==t,],supporting_settings()$top_per_tissue)
    text(.035,if(t=='Root').85 else .46,paste(t,'proteins |',sum(salt$tissue==t),'supported loci'),adj=0,font=2,cex=1.03)
    supporting_print_table(a,c('Gene_name','locus','RNA_log2FC','protein_log2FC','RNA_q','protein_q','profile_Pearson','joint_q_global'),
      c('Gene name','AGI','RNA log2','Protein log2','RNA q','Protein q','Profile r','Global q'),
      c(.045,.39,.49,.59,.69,.79,.87,.94),y=if(t=='Root').81 else .42,step=.029,cex=.76)
  }
  supporting_text(paste(sum(salt$contrast=='salt_OE'),'NIT2 salt responses pass the screen;',sum(z$audit$comparison=='No0_genotype' & z$audit$assaywise_supported & !z$audit$profile_positive),'same-direction, dual-significant result(s) fail positive profile correlation and remain in the audit.',sum(salt$joint_q_global<.05),'selected salt results pass original joint global BH.'),.075,cex=.82)

  supporting_page('Top supporting genotype differences','Ten distinct loci per genotype comparison; the strongest eligible tissue represents each displayed gene')
  for(co in c('baseline_genotype','saline_genotype')){
    a<-genotype[genotype$contrast==co,];a<-a[!duplicated(a$locus),];a<-head(a,10)
    text(.035,if(co=='baseline_genotype').85 else .46,supporting_label(co),adj=0,font=2,cex=1.03)
    supporting_print_table(a,c('Gene_name','tissue','RNA_log2FC','protein_log2FC','RNA_q','protein_q','profile_Pearson','joint_q_global'),
      c('Gene name','Tissue','RNA log2','Protein log2','RNA q','Protein q','Profile r','Global q'),
      c(.045,.39,.49,.59,.69,.79,.87,.94),y=if(co=='baseline_genotype').81 else .42,step=.029,cex=.76)
  }
  supporting_text(paste('NIT2 minus No-0 is evaluated separately at control and at salt. These contrasts differ from a genotype-by-salt interaction;',sum(d$contrast=='interaction'),'interactions meet this screen.'),.075,cex=.84)

  supporting_page('Salt support across all three datasets','Array 1: Col-0, 75 minus 0 mM; Array 2: No-0, 100 minus 0 mM; protein: No-0, 75 minus 0 mM')
  supporting_print_table(z$triple,c('Gene_name','tissue','Array1_log2FC','Array2_log2FC','protein_log2FC','Array1_q','Array2_q','protein_q'),
    c('Gene name','Tissue','Array 1 log2','Array 2 log2','Protein log2','Array 1 q','Array 2 q','Protein q'),
    c(.045,.38,.48,.58,.68,.78,.86,.94),y=.81,step=.045,cex=.79)
  supporting_text('All three effects have the same sign and each assay q is below 0.05. The same protein salt estimate is shared by the two array comparisons; it is not counted as two independent protein validations. Array 2/protein four-condition profile correlation is positive. Array 1 remains a separate background and dose-series context.',.30)
  supporting_text(paste(length(intersect(z$dual$locus,z$triple$locus)),'genes in the salt-plus-genotype priority list also meet this three-dataset salt screen. This list identifies additional candidates with broad agreement across datasets. Full names and all q values are preserved in supporting_three_datasets.csv.'),.17)

  supporting_page('Strongest original joint evidence','Directional-conjunction BH below 0.05 across the full original gene-by-family testing universe')
  a<-supporting_order(z$audit[z$audit$joint_q_global<.05&z$audit$comparison=='No0_genotype',])
  a$comparison_label<-supporting_label(a$contrast)
  supporting_print_table(a,c('Gene_name','tissue','comparison_label','RNA_log2FC','protein_log2FC','joint_q_global','joint_q_BY'),
    c('Gene name','Tissue','Comparison','RNA log2','Protein log2','Global BH','Global BY'),
    c(.045,.36,.49,.64,.75,.86,.94),y=.81,step=.05,cex=.82)
  supporting_text(paste(nrow(a),'gene-by-comparison results represent',length(unique(a$locus)),'loci.',sum(a$complete_case_joint_global_BH<.05,na.rm=TRUE),'pass the complete-case global BH sensitivity;',sum(a$single_transcript_joint_global_BH<.05,na.rm=TRUE),'pass the single-transcript global BH sensitivity;',sum(a$joint_q_BY<.05),'pass the original dependence-robust global BY sensitivity. These retain the original complete testing families.'),.43)
  supporting_text('The broader assay-wise candidate screen answers a different prioritization question: which significant measured changes support each other across assays? Its lists retain the original assay q values, joint q values, identity support and sensitivity outcomes, without recalculating significance after selection or treating molecular features or four means as biological replication.',.29)
  supporting_text('Names use the supplied array annotation, or a meaningful supplied protein gene label when the array has no symbol. Functional descriptions and AGI identifiers remain where no symbol is supplied. Protein assignments require a single accession, unique AGI, reported proteotypic support and at least two observed pools in each genotype/treatment cell. Missing values are not filled.',.14)
}

supporting_figure <- function(result,destination) {
  selected<-head(result$dual,supporting_settings()$profile_genes);signal<-supporting_signal(result)
  if(!nrow(selected)){
    publication_figure(destination,9,5,function(){plot.new();text(.5,.6,'No gene meets both salt and genotype support criteria',cex=1.1);text(.5,.4,'No nonsignificant candidate has been substituted.',cex=.9)});return(invisible(destination))
  }
  n<-nrow(selected)
  publication_figure(destination,9,9.1,function(){
    par(mfrow=c(3,n),mar=c(4.5,4,3.4,.8),oma=c(1.5,0,1.9,0))
    for(assay in c('Array 2','Protein','Correlation'))for(i in seq_len(n)){
      a<-selected[i,];s<-signal[signal$locus==a$locus&signal$protein_tissue==a$tissue,]
      label<-sub(' \\[AT[1-5CM]G[0-9]+\\]$','',a$Gene_name)
      if(nchar(label)>29)label<-paste0(substr(label,1,26),'...')
      if(assay!='Correlation'){
        s<-s[s$assay==assay,];means<-vapply(supporting_settings()$profile_cells,function(cell)mean(s$log2_signal[s$cell==cell],na.rm=TRUE),0.)
        counts<-vapply(supporting_settings()$profile_cells,function(cell)sum(s$observed[s$cell==cell]),0L)
        salt<-if(assay=='Array 2')100 else 75
        publication_axes(c(.55,4.45),publication_range(s$log2_signal),ylab=if(assay=='Array 2')'Log2 RNA signal' else 'Log2 protein abundance',at=1:4,
          labels=paste0(c('No-0\n0','No-0\n', 'NIT2\n0','NIT2\n'),c('',salt,'',salt),'\nn=',counts))
        for(j in 1:4){
          v<-s[s$cell==supporting_settings()$profile_cells[j],]
          color<-publication_colors[if(j<=2)'No-0' else 'NIT2']
          segments(j-.23,means[j],j+.23,means[j],lwd=1.5,col=color)
          points(j+seq(-.13,.13,length.out=nrow(v)),v$log2_signal,pch=if(j%%2==1)21 else 24,
            bg=color,col='white',cex=1.1,lwd=.5)
        }
        publication_title(LETTERS[(if(assay=='Array 2')0 else n)+i],label,paste(a$locus,'|',if(assay=='Array 2')'Array 2 combined' else paste(a$tissue,'protein')))
      }else{
        means<-lapply(c('Array 2','Protein'),function(type)vapply(supporting_settings()$profile_cells,function(cell)mean(s$log2_signal[s$cell==cell&s$assay==type],na.rm=TRUE),0.))
        publication_axes(publication_range(means[[1]],.23),publication_range(means[[2]],.25),
          xlab='RNA condition mean (log2)',ylab='Protein condition mean (log2)')
        points(means[[1]],means[[2]],pch=c(21,24,21,24),bg=publication_colors[c('No-0','No-0','NIT2','NIT2')],col='white',cex=1.3)
        publication_pool_labels(means[[1]],means[[2]],c('WT0','WTS','OE0','OES'))
        publication_title(LETTERS[2*n+i],paste0('Four means | r = ',sprintf('%.2f',a$profile_Pearson)),
          paste0('Spearman = ',sprintf('%.2f',a$profile_Spearman),' | descriptive; n = 4 means'))
      }
    }
    publication_key(FALSE,'Circles: control; triangles: salt. Raw pool points and condition means; no imputation.\nWT: No-0; OE: NIT2; S: salt. Four means are not four paired biological replicates.')
  })
}

build_supporting_audit <- function(inputs,settings,destination)write_table(supporting_results(inputs)$audit,destination)
build_supporting_ranked <- function(inputs,settings,destination)write_table(supporting_results(inputs)$ranked,destination)
build_supporting_dual <- function(inputs,settings,destination)write_table(supporting_results(inputs)$dual,destination)
build_supporting_three <- function(inputs,settings,destination)write_table(supporting_results(inputs)$triple,destination)
build_supporting_summary <- function(inputs,settings,destination)write_table(supporting_results(inputs)$summary,destination)
build_supporting_signal <- function(inputs,settings,destination)write_table(supporting_signal(supporting_results(inputs)),destination)
build_supporting_report <- function(inputs,settings,destination)supporting_report(supporting_results(inputs),destination)
build_supporting_profiles <- function(inputs,settings,destination)supporting_figure(supporting_results(inputs),destination)

supporting_registry <- function() {
  base<-c('microarray_genotype','microarray_root','microarray_shoot','proteomics_root','proteomics_shoot',
    'saved_proteomics_root','saved_proteomics_shoot','gene_annotation')
  entry<-function(builder,filename,figure=FALSE){
    out<-list(inputs=base,builder=builder,filename=filename,status='post hoc supporting-trend candidate prioritization')
    if(figure)out$formats<-c('pdf','svg','tiff','png');out
  }
  list(table_supporting_audit=entry('build_supporting_audit','supporting_audit.csv'),
    table_supporting_ranked=entry('build_supporting_ranked','supporting_ranked.csv'),
    table_supporting_dual=entry('build_supporting_dual','supporting_dual.csv'),
    table_supporting_three=entry('build_supporting_three','supporting_three_datasets.csv'),
    table_supporting_summary=entry('build_supporting_summary','supporting_summary.csv'),
    table_supporting_signal=entry('build_supporting_signal','supporting_profile_signals.csv'),
    report_supporting_trends=entry('build_supporting_report','supporting_trends_report.pdf'),
    figure_supporting_profiles=entry('build_supporting_profiles','supporting_profiles.pdf',TRUE))
}
