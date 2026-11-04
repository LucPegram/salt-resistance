# Dataset-specific descriptive rankings. No experimental values are embedded.
.toprank_cache <- new.env(parent=emptyenv())

toprank_contrast_label <- function(dataset) {
  if(dataset=="array1")return(paste0(analysis_settings$toprank$array1_endpoint_mM," mM minus control; Col-0"))
  salt<-if(dataset=="array2")100L else 75L
  switch(analysis_settings$toprank$factorial_contrast,
    salt_WT=paste0("No-0: ",salt," mM minus control"),
    salt_OE=paste0("NIT2: ",salt," mM minus control"),
    baseline_genotype="NIT2 minus No-0 at 0 mM",
    saline_genotype=paste0("NIT2 minus No-0 at ",salt," mM"),
    interaction="Difference between NIT2 and No-0 salt responses",
    average_genotype="NIT2 minus No-0, averaged across treatments",
    average_salt="Salt minus control, averaged across genotypes",stop("Unknown ranking contrast"))
}

read_toprank_annotation <- function(path,spec) {
  d<-read.csv(path,check.names=FALSE,stringsAsFactors=FALSE)
  assert_columns(d,c("Gene ID","Gene info"),"supplied gene annotation")
  list(data=d,spec=spec)
}

toprank_annotation <- function(x) {
  d<-x$data;d$locus<-sub("\\.[0-9]+$","",trimws(d[["Gene ID"]]))
  d<-d[grepl("^AT[1-5CM]G[0-9]+$",d$locus),]
  rows<-lapply(split(seq_len(nrow(d)),d$locus),function(i){
    info<-unique(trimws(d[["Gene info"]][i]));info<-info[!is.na(info)]
    aliases<-unique(unlist(lapply(info,function(s){
      prefix<-trimws(strsplit(s,"|",fixed=TRUE)[[1]][1])
      if(!startsWith(prefix,":"))return(character())
      trimws(strsplit(sub("^: *","",prefix),",",fixed=TRUE)[[1]])
    })))
    aliases<-aliases[nzchar(aliases)]
    descriptions<-unique(vapply(info,function(s){
      parts<-trimws(strsplit(s,"|",fixed=TRUE)[[1]])
      parts<-parts[nzchar(parts)&!startsWith(parts,":")&!startsWith(parts,"chr")]
      if(length(parts))trimws(strsplit(parts[1],";",fixed=TRUE)[[1]][1]) else "Uncharacterized gene"
    },""))
    locus<-d$locus[i[1]]
    data.frame(locus=locus,Gene_name=if(length(aliases))aliases[1] else paste0(descriptions[1]," [",locus,"]"),
      Aliases=paste(aliases,collapse="; "),Description=paste(descriptions,collapse=" | "),
      Name_source=if(length(aliases))"supplied Gene info symbol" else "supplied functional name; no symbol supplied",
      stringsAsFactors=FALSE)
  })
  do.call(rbind,rows)
}

toprank_add_names <- function(d,annotation) {
  idx<-match(d$locus,annotation$locus)
  for(col in c("Gene_name","Aliases","Description","Name_source"))d[[col]]<-annotation[[col]][idx]
  missing<-is.na(d$Gene_name)|!nzchar(d$Gene_name)
  d$Gene_name[missing]<-paste0("Uncharacterized gene [",d$locus[missing],"]")
  d$Name_source[missing]<-"no supplied symbol or description; AGI retained"
  d
}

toprank_dose <- function(inputs) {
  # Locus-level counterpart of the established categorical paired dose model.
  xs<-inputs[c("microarray_root","microarray_shoot")]
  ys<-lapply(xs,function(x)integration_loci(x)$y)
  stopifnot(identical(rownames(ys[[1]]),rownames(ys[[2]])))
  y<-do.call(cbind,ys)
  m<-do.call(rbind,lapply(xs,function(x)data.frame(sample=x$samples,tissue=x$spec$tissue,
    dose=as.integer(ifelse(startsWith(x$samples,"PN"),"0",sub(" .*","",x$samples))),
    replicate=as.integer(sub(".*B","",x$samples)))))
  m$pool<-factor(paste(m$dose,m$replicate,sep="_"))
  m$cell<-factor(paste0(m$tissue,m$dose),levels=unlist(lapply(c("Root","Shoot"),
    function(t)paste0(t,c(0,50,75,100)))))
  design<-model.matrix(~0+cell,m);colnames(design)<-levels(m$cell)
  rho<-limma::duplicateCorrelation(y,design,block=m$pool)$consensus.correlation
  endpoint<-analysis_settings$toprank$array1_endpoint_mM
  cm<-limma::makeContrasts(contrasts=c(paste0("Root",endpoint,"-Root0"),
    paste0("Shoot",endpoint,"-Shoot0")),levels=design)
  colnames(cm)<-c("Root","Shoot")
  fit<-moderated_contrasts(y,design,cm,block=m$pool,correlation=rho)$table
  fit$locus<-fit$feature;fit$tissue<-fit$contrast
  fit$contrast<-paste0(endpoint," mM minus control; Col-0")
  fit$eligible<-TRUE;fit$eligibility<-"locus-level paired dose model"
  list(table=fit,y=y,metadata=m,pair_correlation=rho)
}

toprank_protein <- function(inputs) {
  results<-list();objects<-list()
  for(tissue in c("Root","Shoot")){
    id<-paste0("proteomics_",tolower(tissue));x<-inputs[[id]]
    z<-revision_protein(x);objects[[tissue]]<-z
    saved<-inputs[[paste0("saved_",id)]]$data
    assert_columns(saved,c("Protein.Group","AGI"),"protein gene annotation")
    if(anyDuplicated(saved$Protein.Group))stop("Duplicate protein annotation keys")
    a<-z$annotations;key<-match(a$Protein.Group,saved$Protein.Group)
    locus<-trimws(saved$AGI[key]);valid<-!is.na(locus)&grepl("^AT[1-5CM]G[0-9]+$",locus)
    duplicate<-duplicated(locus)|duplicated(locus,fromLast=TRUE)
    single<-!grepl(";",a$Protein.Group,fixed=TRUE)
    support<-!is.na(a$N.Proteotypic.Sequences)&a$N.Proteotypic.Sequences>=1
    eligible<-valid&!duplicate&single&support&z$selection$observed_min2
    status<-vapply(seq_len(nrow(a)),function(i){
      why<-c(if(!valid[i])"invalid/unavailable AGI",if(duplicate[i])"duplicate locus",
        if(!single[i])"multiple accessions",if(!support[i])"no reported proteotypic support",
        if(!z$selection$observed_min2[i])"fewer than two observed pools in a cell")
      if(length(why))paste(why,collapse="; ") else "eligible uniquely assigned protein group"
    },"")
    # Use the primary-model statistics. Other retained groups remain in the audit.
    tab<-z$table[z$table$method=="observed_min2"&z$table$contrast==analysis_settings$toprank$factorial_contrast,]
    out<-data.frame(feature=a$Protein.Group,locus=locus,tissue=tissue,contrast=analysis_settings$toprank$factorial_contrast,
      eligible=eligible,eligibility=status,proteotypic_sequences=a$N.Proteotypic.Sequences,
      reported_gene_label=a$Genes,stringsAsFactors=FALSE)
    for(col in c("log2FC","SE","lower95","upper95","P"))out[[col]]<-tab[[col]][match(out$feature,tab$feature)]
    results[[tissue]]<-out
  }
  list(table=do.call(rbind,results),objects=objects)
}

toprank_select <- function(d,n=20L,cutoff=.05) {
  # Correction precedes tissue selection and distinct-locus ranking.
  d$q_BH<-NA_real_;k<-which(d$eligible&is.finite(d$P))
  d$q_BH[k]<-p.adjust(d$P[k],"BH")
  d$direction<-ifelse(d$log2FC>0,"Up",ifelse(d$log2FC<0,"Down","Unchanged"))
  ok<-d$eligible&is.finite(d$q_BH)&d$q_BH<cutoff&is.finite(d$log2FC)
  candidate<-which(ok)
  ordered<-candidate[order(-abs(d$log2FC[candidate]),d$q_BH[candidate],d$locus[candidate],d$tissue[candidate])]
  chosen<-ordered[!duplicated(d$locus[ordered])]
  d$representative<-seq_len(nrow(d))%in%chosen
  d$rank<-NA_integer_;d$selected_top20<-FALSE;d$plotted_top10<-FALSE
  for(direction in c("Up","Down")){
    indices<-chosen[d$direction[chosen]==direction]
    if(length(indices)<n)stop("Fewer than ",n," eligible significant ",direction," loci; do not fill with nonsignificant rows")
    d$rank[indices]<-seq_along(indices)
    d$selected_top20[head(indices,n)]<-TRUE;d$plotted_top10[head(indices,10)]<-TRUE
  }
  d$opposite_sign_other_tissue<-FALSE
  for(indices in split(seq_len(nrow(d)),d$locus)){
    for(i in indices)d$opposite_sign_other_tissue[i]<-any(d$tissue[indices]!=d$tissue[i]&
      is.finite(d$log2FC[indices])&sign(d$log2FC[indices])!=sign(d$log2FC[i]),na.rm=TRUE)
  }
  d
}

toprank_result <- function(dataset,inputs) {
  key<-digest::digest(list(dataset,inputs,analysis_settings$toprank),algo="sha256")
  if(exists(key,.toprank_cache,inherits=FALSE))return(get(key,.toprank_cache,inherits=FALSE))
  annotation<-toprank_annotation(inputs$gene_annotation)
  if(dataset=="array1"){
    z<-toprank_dose(inputs);d<-z$table
  }else if(dataset=="array2"){
    z<-integration_array(inputs$microarray_genotype)
    d<-z$table[z$table$contrast==analysis_settings$toprank$factorial_contrast,]
    d$locus<-d$feature;d$tissue<-"Combined";d$eligible<-TRUE;d$eligibility<-"locus-level genotype array model"
  }else if(dataset=="proteomics"){
    z<-toprank_protein(inputs);d<-z$table
  }else stop("Unknown top-ranking dataset")
  d<-toprank_add_names(d,annotation)
  d<-toprank_select(d,analysis_settings$toprank$table_per_direction,analysis_settings$toprank$BH_cutoff)
  selected<-d[d$selected_top20,];selected<-selected[order(match(selected$direction,c("Up","Down")),selected$rank),]
  selected$fold_change<-2^selected$log2FC
  selected$dataset<-dataset
  selected$contrast_label<-toprank_contrast_label(dataset)
  selected$BH_scope<-if(dataset=="array2")"all tested loci for the selected contrast" else
    "all eligible locus-by-tissue comparisons for the selected contrast"
  cols<-c("dataset","direction","rank","Gene_name","Aliases","Description","locus","feature","tissue",
    "contrast","contrast_label","log2FC","fold_change","SE","lower95","upper95","P","q_BH","BH_scope",
    "Name_source","opposite_sign_other_tissue","eligibility")
  if(dataset=="proteomics")cols<-c(cols,"proteotypic_sequences","reported_gene_label")
  result<-list(table=selected[,cols],audit=d,model=z,dataset=dataset)
  assign(key,result,.toprank_cache);result
}

toprank_signal <- function(result) {
  selected<-result$table[result$table$rank<=analysis_settings$toprank$figure_per_direction,]
  rows<-lapply(seq_len(nrow(selected)),function(i){
    a<-selected[i,];z<-result$model
    if(result$dataset=="proteomics"){
      z<-z$objects[[a$tissue]];m<-z$metadata;values<-z$y[a$feature,]
      tissue<-rep(a$tissue,length(values));genotype<-m$genotype;salt<-m$salt_mM
    }else if(result$dataset=="array1"){
      m<-z$metadata;values<-z$y[a$locus,];tissue<-m$tissue
      genotype<-rep("Col-0",length(values));salt<-m$dose
    }else{
      m<-z$metadata;values<-z$y[a$locus,];tissue<-rep("Combined",nrow(m));genotype<-m$genotype;salt<-m$salt_mM
    }
    data.frame(direction=a$direction,rank=a$rank,Gene_name=a$Gene_name,locus=a$locus,feature=a$feature,
      selected_tissue=a$tissue,sample=m$sample,tissue=tissue,genotype=genotype,salt_mM=salt,
      log2_signal=as.numeric(values),observed=is.finite(values),stringsAsFactors=FALSE)
  })
  do.call(rbind,rows)
}

toprank_table_pdf <- function(result,destination) {
  grDevices::cairo_pdf(destination,width=11.7,height=8.3,family="sans",pointsize=10)
  on.exit(grDevices::dev.off(),add=TRUE)
  labels<-c(array1="Array 1 | Col-0 dose response",array2="Array 2 | genotype experiment",proteomics="Proteomics")
  for(direction in c("Up","Down")){
    d<-result$table[result$table$direction==direction,]
    par(mar=c(0,0,0,0));plot.new();plot.window(xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i")
    text(.035,.95,paste(labels[result$dataset],"|",direction,"regulated"),adj=0,cex=1.5,font=2,col="#253342")
    contrast<-toprank_contrast_label(result$dataset)
    text(.035,.913,paste("Ranking contrast:",contrast,"| BH q < 0.05; greatest absolute log2 difference per distinct locus"),adj=0,cex=.86,col="#5D6975")
    xx<-c(.04,.10,.54,.68,.80,.91)
    rect(.035,.847,.965,.886,col="#E8EEF2",border=NA)
    for(j in seq_along(xx))text(xx[j],.866,c("Rank","Gene name / AGI","Tissue","Log2 change","95% interval","BH q")[j],adj=if(j==2)0 else .5,font=2,cex=.95)
    for(i in seq_len(nrow(d))){
      y<-.824-(i-1)*.035
      if(i%%2==0)rect(.035,y-.017,.965,y+.017,col="#F5F7F9",border=NA)
      name<-d$Gene_name[i];if(nchar(name)>58)name<-paste0(substr(name,1,55),"...")
      text(xx[1],y,d$rank[i],cex=.83)
      text(xx[2],y+.006,name,adj=0,cex=.84)
      text(xx[2],y-.0085,d$locus[i],adj=0,cex=.62,col="#687985")
      text(xx[3],y,d$tissue[i],cex=.82)
      text(xx[4],y,sprintf("%+.2f",d$log2FC[i]),cex=.88,font=2)
      text(xx[5],y,sprintf("%.2f to %.2f",d$lower95[i],d$upper95[i]),cex=.8)
      text(xx[6],y,formatC(d$q_BH[i],format="g",digits=3),cex=.82)
    }
    text(.035,.075,"Intervals are pointwise. Ranking and selection are descriptive; separate tissues are not pooled measurements.",adj=0,cex=.82,col="#5D6975")
    text(.035,.048,"Gene symbols follow supplied annotation; functional names retain AGI where a symbol is unavailable. Full names and statistics: companion CSV.",adj=0,cex=.77,col="#5D6975")
  }
}

toprank_figure <- function(result,destination) {
  d<-result$table[result$table$rank<=analysis_settings$toprank$figure_per_direction,]
  signals<-toprank_signal(result);dataset<-result$dataset
  publication_figure(destination,14,15.4,function(){
    layout(cbind(matrix(1:10,nrow=5,byrow=TRUE),matrix(11:20,nrow=5,byrow=TRUE)))
    par(mar=c(3.5,3.7,2.8,.9),oma=c(.4,0,5.1,0),cex=1)
    for(i in seq_len(nrow(d))){
      a<-d[i,];s<-signals[signals$locus==a$locus,]
      groups<-if(dataset=="array1")s$tissue else s$genotype
      lev<-if(dataset=="array1")c("Root","Shoot") else c("No-0","NIT2")
      x<-sort(unique(s$salt_mM));gap<-max(x)*.025
      yr<-publication_range(s$log2_signal,.17)
      publication_axes(c(min(x)-gap*3,max(x)+gap*3),yr,"NaCl (mM)",
        if(dataset=="proteomics")"Log2 abundance" else "Log2 normalized signal",at=x)
      for(j in seq_along(lev)){
        k<-groups==lev[j]&s$observed;v<-s[k,];color<-publication_colors[lev[j]]
        means<-tapply(v$log2_signal,v$salt_mM,mean)
        lines(as.numeric(names(means)),means,col=color,lwd=1.7,lty=j)
        offsets<-ave(seq_len(nrow(v)),v$salt_mM,FUN=function(z)seq(-.35,.35,length.out=length(z)))
        points(v$salt_mM+(j-1.5)*gap+offsets*gap,v$log2_signal,pch=c(21,24)[j],
          bg=color,col="white",lwd=.5,cex=1.04)
      }
      label<-sub(" \\[AT[1-5CM]G[0-9]+\\]$","",a$Gene_name)
      title_lines<-strwrap(label,width=36)
      if(length(title_lines)>2L)title_lines<-c(title_lines[1],paste0(title_lines[2],"..."))
      title_lines[1]<-paste0(substr(a$direction,1,1),sprintf("%02d",a$rank),"  ",title_lines[1])
      for(j in seq_along(title_lines))mtext(title_lines[j],side=3,line=.9+(length(title_lines)-j)*.9,
        adj=0,font=2,cex=.87)
      publication_caption(sprintf("Panel %s%02d (%s): %s; %s; log2 change %+.2f; BH q = %.3g.",
        substr(a$direction,1,1),a$rank,a$Gene_name,a$locus,
        if(dataset!="array2")paste("selected tissue",a$tissue) else "combined tissues",a$log2FC,a$q_BH))
    }
    title<-c(array1="Array 1 | Col-0 salt dose response",array2="Array 2 | combined-tissue genotype experiment",proteomics="Proteomics | tissue-specific protein abundance")[dataset]
    mtext(title,side=3,outer=TRUE,line=3.4,font=2,cex=1.32,col="#253342")
    mtext("Top 10 up-regulated",side=3,outer=TRUE,line=1.75,adj=.23,font=2,cex=1.05,col="#007C5B")
    mtext("Top 10 down-regulated",side=3,outer=TRUE,line=1.75,adj=.77,font=2,cex=1.05,col="#B2473E")
    key<-if(dataset=="array1")"Blue circles / solid: roots    Orange triangles / dashed: shoots" else
      "Blue circles / solid: No-0    Purple triangles / dashed: NIT2 overexpressor"
    mtext(key,side=3,outer=TRUE,line=.55,cex=.87,col="#5D6975")
    contrast<-toprank_contrast_label(dataset)
    publication_caption(paste("Selection:",contrast,"| BH q < 0.05; one representative tissue per locus where applicable"))
    publication_caption("Points: independently grown plant pools. Lines: condition means, not paired trajectories. No missing protein value is imputed.")
    publication_caption("Gene selection and ranking use the same data shown; panels are descriptive and have independent signal scales. Full names: companion table.")
  })
}

build_toprank <- function(inputs,settings,destination,dataset,kind) {
  result<-toprank_result(dataset,inputs)
  switch(kind,table=write_table(result$table,destination),audit=write_table(result$audit,destination),
    signal=write_table(toprank_signal(result),destination),print=toprank_table_pdf(result,destination),
    figure=toprank_figure(result,destination),stop("Unknown ranking output kind"))
}

toprank_registry <- function() {
  specs<-list(array1=c("microarray_root","microarray_shoot","gene_annotation"),
    array2=c("microarray_genotype","gene_annotation"),
    proteomics=c("proteomics_root","proteomics_shoot","saved_proteomics_root","saved_proteomics_shoot","gene_annotation"))
  out<-list()
  for(dataset in names(specs))for(kind in c("table","audit","signal","print","figure")){
    id<-paste("toprank",dataset,kind,sep="_")
    out[[id]]<-list(inputs=specs[[dataset]],builder=paste0("build_",id),
      filename=paste0(dataset,"_",kind,if(kind%in%c("figure","print"))".pdf" else ".csv"),
      status="post hoc dataset ranking; descriptive selected trajectories")
    if(kind=="figure")out[[id]]$formats<-c("pdf","svg","tiff","png")
  }
  out
}

# Explicit wrappers keep every registry builder independently addressable.
build_toprank_array1_table<-function(inputs,settings,destination)build_toprank(inputs,settings,destination,"array1","table")
build_toprank_array1_audit<-function(inputs,settings,destination)build_toprank(inputs,settings,destination,"array1","audit")
build_toprank_array1_signal<-function(inputs,settings,destination)build_toprank(inputs,settings,destination,"array1","signal")
build_toprank_array1_print<-function(inputs,settings,destination)build_toprank(inputs,settings,destination,"array1","print")
build_toprank_array1_figure<-function(inputs,settings,destination)build_toprank(inputs,settings,destination,"array1","figure")
build_toprank_array2_table<-function(inputs,settings,destination)build_toprank(inputs,settings,destination,"array2","table")
build_toprank_array2_audit<-function(inputs,settings,destination)build_toprank(inputs,settings,destination,"array2","audit")
build_toprank_array2_signal<-function(inputs,settings,destination)build_toprank(inputs,settings,destination,"array2","signal")
build_toprank_array2_print<-function(inputs,settings,destination)build_toprank(inputs,settings,destination,"array2","print")
build_toprank_array2_figure<-function(inputs,settings,destination)build_toprank(inputs,settings,destination,"array2","figure")
build_toprank_proteomics_table<-function(inputs,settings,destination)build_toprank(inputs,settings,destination,"proteomics","table")
build_toprank_proteomics_audit<-function(inputs,settings,destination)build_toprank(inputs,settings,destination,"proteomics","audit")
build_toprank_proteomics_signal<-function(inputs,settings,destination)build_toprank(inputs,settings,destination,"proteomics","signal")
build_toprank_proteomics_print<-function(inputs,settings,destination)build_toprank(inputs,settings,destination,"proteomics","print")
build_toprank_proteomics_figure<-function(inputs,settings,destination)build_toprank(inputs,settings,destination,"proteomics","figure")
