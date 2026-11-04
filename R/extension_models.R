extension_fit <- function(y,metadata,design,contrasts,block=NULL,correlation=NULL) {
  z<-moderated_contrasts(y,design,contrasts,block,correlation)
  z$y<-y;z$metadata<-metadata;z$design<-design;z$contrasts<-contrasts
  z$block<-block;z$correlation<-correlation;z
}

extension_dose <- function(inputs) {
  xs<-inputs[c("microarray_root","microarray_shoot")]
  yy<-lapply(xs,function(x)integration_loci(x)$y)
  stopifnot(identical(rownames(yy[[1]]),rownames(yy[[2]])))
  y<-do.call(cbind,yy)
  m<-do.call(rbind,lapply(xs,function(x)data.frame(sample=x$samples,tissue=x$spec$tissue,
    dose=as.integer(ifelse(startsWith(x$samples,"PN"),"0",sub(" .*","",x$samples))),
    replicate=as.integer(sub(".*B","",x$samples)))))
  m$pool<-factor(paste(m$dose,m$replicate,sep="_"))
  m$cell<-factor(paste0(m$tissue,m$dose),levels=unlist(lapply(c("Root","Shoot"),
    function(t)paste0(t,c(0,50,75,100)))))
  design<-model.matrix(~0+cell,m);colnames(design)<-levels(m$cell)
  rho<-limma::duplicateCorrelation(y,design,block=m$pool)$consensus.correlation
  pairs<-list(`50v0`=c(50,0),`75v0`=c(75,0),`100v0`=c(100,0),
    `75v50`=c(75,50),`100v75`=c(100,75))
  cm<-matrix(0,ncol(design),15,dimnames=list(colnames(design),
    unlist(lapply(c("Root","Shoot","Difference"),function(t)paste(t,names(pairs),sep="_")))))
  for(k in names(pairs)){
    doses<-pairs[[k]]
    for(t in c("Root","Shoot"))cm[paste0(t,doses),paste(t,k,sep="_")]<-c(1,-1)
    cm[,paste0("Difference_",k)]<-cm[,paste0("Shoot_",k)]-cm[,paste0("Root_",k)]
  }
  out<-extension_fit(y,m,design,cm,m$pool,rho)
  out$table$q_global<-p.adjust(out$table$P,"BH")
  out$joint<-do.call(rbind,lapply(c("Root","Shoot","Difference"),function(t){
    z<-limma::topTable(out$contrast_fit,coef=paste(t,names(pairs)[1:3],sep="_"),number=Inf,sort.by="none")
    data.frame(locus=rownames(z),test=t,F=z$F,P=z$P.Value,q=z$adj.P.Val)
  }))
  out$joint$q_global<-p.adjust(out$joint$P,"BH")
  means<-vapply(levels(m$cell),function(g)rowMeans(y[,m$cell==g,drop=FALSE]),numeric(nrow(y)))
  rownames(means)<-rownames(y);out$means<-means
  patterns<-lapply(c("Root","Shoot"),function(t){
    mu<-means[,paste0(t,c(0,50,75,100)),drop=FALSE]
    increments<-mu[,2:4]-mu[,1:3]
    pattern<-ifelse(apply(increments,1,min)>=0,"increasing means",
      ifelse(apply(increments,1,max)<=0,"decreasing means","non-monotonic means"))
    joint<-out$joint[out$joint$test==t,]
    qq<-vapply(paste(t,names(pairs)[1:3],sep="_"),function(co){
      z<-out$table[out$table$contrast==co,];z$q_global[match(rownames(y),z$feature)]
    },numeric(nrow(y)))
    first<-rep(NA_real_,nrow(y))
    for(i in 3:1)first[qq[,i]<.05]<-c(50,75,100)[i]
    data.frame(locus=rownames(y),tissue=t,mean_pattern=pattern,first_tested_dose_q_global_lt_05=first,
      endpoint_log2FC=mu[,4]-mu[,1],joint_P=joint$P[match(rownames(y),joint$locus)],
      joint_q_global=joint$q_global[match(rownames(y),joint$locus)],
      interpretation="Mean shape is descriptive; first significant tested dose is not a physiological threshold")
  })
  out$patterns<-do.call(rbind,patterns);out
}

extension_external <- function(paths) {
  p198<-extension_read_platform(paths[["GPL198.soft"]],"AGI")
  p17133<-extension_read_platform(paths[["GPL17133.soft"]],"gene_assignment")
  control<-extension_read_matrix(paths[["GSE5620_series_matrix.txt.gz"]])
  salt<-extension_read_matrix(paths[["GSE5623_series_matrix.txt.gz"]])
  stopifnot(identical(rownames(control$y),rownames(salt$y)),all(control$y>0),all(salt$y>0))
  titles<-c(control$title,salt$title)
  m<-data.frame(sample=c(control$sample,salt$sample),title=titles,
    tissue=ifelse(grepl("-Roots-",titles),"Root","Shoot"),
    treatment=ifelse(grepl("_Control-",titles),"control","salt"),
    hours=as.numeric(sub(".*-([0-9.]+)h_.*","\\1",titles)),
    replicate=as.integer(sub(".*_Rep","",titles)),source="GSE5620/GSE5623")
  times<-c(.5,1,3,6,12,24)
  y<-cbind(control$y,salt$y);out<-list();audit<-list();samples<-list()
  for(t in c("Root","Shoot")){
    k<-m$tissue==t&m$hours%in%times;meta<-m[k,]
    # Re-normalize deposited positive MAS5 values jointly within tissue, before
    # unique-AGI collapse. No project arrays enter this transformation.
    raw<-limma::normalizeBetweenArrays(log2(y[,k,drop=FALSE]),method="quantile")
    mapped<-extension_map_matrix(raw,p198)
    meta$cell<-factor(paste0(meta$treatment,"_",meta$hours))
    design<-model.matrix(~0+cell,meta);colnames(design)<-levels(meta$cell)
    cm<-matrix(0,ncol(design),length(times),dimnames=list(colnames(design),paste0("hour",times)))
    for(i in seq_along(times))cm[c(paste0("salt_",times[i]),paste0("control_",times[i])),i]<-c(1,-1)
    z<-extension_fit(mapped$y,meta,design,cm)
    z$table$reference<-paste0("AtGen_",t,"_",z$table$contrast)
    z$table$tissue<-t;z$table$source<-"GSE5620/GSE5623";z$table$salt_mM<-150
    z$table$hours<-as.numeric(sub("hour","",z$table$contrast))
    out[[t]]<-z$table;aa<-mapped$audit;aa$source<-paste0("AtGen_",t);audit[[t]]<-aa
    samples[[t]]<-meta[,names(m)]
  }
  hydro<-extension_read_matrix(paths[["GSE223207_series_matrix.txt.gz"]])
  stopifnot(min(hydro$y)>0,max(hydro$y)<30,
    all(grepl("GCRMA",hydro$processing,fixed=TRUE)))
  # All deposited GCRMA log2 observations are retained, including non-NaCl
  # treatments, as distinct categorical cells for residual variance estimation.
  meta<-data.frame(sample=hydro$sample,title=hydro$title,tissue="Root",
    treatment=hydro$source,hours=.25,replicate=NA_integer_,source="GSE223207")
  meta$cell<-factor(hydro$source)
  design<-model.matrix(~0+cell,meta);colnames(design)<-levels(meta$cell)
  cm<-matrix(0,ncol(design),4,dimnames=list(colnames(design),paste0("dose",c(25,50,75,100))))
  for(dose in c(25,50,75,100))cm[c(paste0("NaCl ",dose,"mM"),"Control hydroponic solution"),paste0("dose",dose)]<-c(1,-1)
  mapped<-extension_map_matrix(hydro$y,p17133)
  z<-extension_fit(mapped$y,meta,design,cm)
  z$table$reference<-paste0("Hydro15min_Root_",z$table$contrast);z$table$tissue<-"Root"
  z$table$source<-"GSE223207";z$table$salt_mM<-as.numeric(sub("dose","",z$table$contrast));z$table$hours<-.25
  out$hydro<-z$table;aa<-mapped$audit;aa$source<-"GSE223207";audit$hydro<-aa
  samples$hydro<-meta[,names(m)]
  list(table=do.call(rbind,out),mapping=do.call(rbind,audit),samples=do.call(rbind,samples))
}

extension_objects <- function(inputs,dose) {
  a<-integration_array(inputs$microarray_genotype)
  aa<-extension_fit(a$y,a$metadata,a$design,a$contrasts)
  selected<-integration_settings()$contrasts
  out<-list(Array1=within(dose,{dataset<-"Array1";tissue<-"Root/Shoot";
      assay<-"RNA";selected_contrasts<-grep("^(Root|Shoot)_(50v0|75v0|100v0)$",colnames(contrasts),value=TRUE)}),
    Array2=within(aa,{dataset<-"Array2";tissue<-"Combined";assay<-"RNA";selected_contrasts<-selected}))
  for(t in c("Root","Shoot")){
    id<-paste0("proteomics_",tolower(t));x<-inputs[[id]];p<-revision_protein(x)
    mapping<-integration_protein_mapping(x,inputs[[paste0("saved_",id)]],rownames(a$y))
    k<-which(mapping$eligible & mapping$complete_observed)
    y<-p$y[match(mapping$protein_group[k],rownames(p$y)),,drop=FALSE];rownames(y)<-mapping$locus[k]
    stopifnot(all(is.finite(y)),!anyDuplicated(rownames(y)))
    design<-model.matrix(~0+cell,p$metadata);colnames(design)<-levels(p$metadata$cell)
    z<-extension_fit(y,p$metadata,design,factorial_contrasts(design))
    z$dataset<-paste0("Protein_",t);z$tissue<-t;z$assay<-"protein"
    z$selected_contrasts<-selected;out[[z$dataset]]<-z
  }
  out
}

extension_rank_scores <- function(y,sets,min_members=5L) {
  # Fractional expression ranks within each biological sample. These are
  # descriptive expression signatures, not estimated cell percentages.
  stopifnot(all(is.finite(y)),!anyDuplicated(rownames(y)))
  ranks<-apply(y,2,rank,ties.method="average")/nrow(y)
  rownames(ranks)<-rownames(y)
  indices<-lapply(sets,function(g)match(intersect(g,rownames(y)),rownames(y)))
  indices<-indices[lengths(indices)>=min_members]
  if(!length(indices))return(matrix(numeric(),nrow=0,ncol=ncol(y),dimnames=list(NULL,colnames(y))))
  scores<-t(vapply(indices,function(i)colMeans(ranks[i,,drop=FALSE]),numeric(ncol(y))))
  colnames(scores)<-colnames(y);scores
}

extension_score_fit <- function(scores,object) {
  if(nrow(scores)<3L)stop("Too few eligible signatures for moderated analysis")
  correlation<-if(is.null(object$block))NULL else
    limma::duplicateCorrelation(scores,object$design,block=object$block)$consensus.correlation
  z<-moderated_contrasts(scores,object$design,object$contrasts,
    block=object$block,correlation=correlation)$table
  z<-z[z$contrast%in%object$selected_contrasts,]
  names(z)[names(z)=="log2FC"]<-"score_difference"
  z$dataset<-object$dataset;z$assay<-object$assay;z$tissue<-object$tissue
  z$pair_correlation<-if(is.null(correlation))NA_real_ else correlation
  z
}

extension_benchmark <- function(external,objects) {
  s<-extension_settings();members<-list();coverage<-list();results<-list();correlations<-list();scores_all<-list()
  refs<-split(external$table,external$table$reference)
  sets<-list()
  for(n in names(refs)){
    d<-refs[[n]];hit<-d$q<s$salt_signature_q & abs(d$log2FC)>=s$salt_signature_abs_log2FC
    for(direction in c("Up","Down")){
      k<-hit & if(direction=="Up")d$log2FC>0 else d$log2FC<0
      sets[[paste(n,direction,sep="__")]]<-d$feature[k]
      if(any(k))members[[length(members)+1L]]<-data.frame(reference=n,direction=direction,locus=d$feature[k],
        external_log2FC=d$log2FC[k],external_q=d$q[k],source=d$source[k],tissue=d$tissue[k])
    }
  }
  for(o in objects){
    halves<-extension_rank_scores(o$y,sets,s$signature_min_members)
    usable<-names(refs)[vapply(names(refs),function(n)all(paste(n,c("Up","Down"),sep="__")%in%rownames(halves)),TRUE)]
    for(n in names(refs)){
      d<-refs[[n]];allowed<-o$tissue%in%c("Root/Shoot","Combined",d$tissue[1])
      coverage[[length(coverage)+1L]]<-data.frame(dataset=o$dataset,reference=n,reference_tissue=d$tissue[1],
        Up_total=length(sets[[paste0(n,"__Up")]]),Down_total=length(sets[[paste0(n,"__Down")]]),
        Up_measured=sum(sets[[paste0(n,"__Up")]]%in%rownames(o$y)),
        Down_measured=sum(sets[[paste0(n,"__Down")]]%in%rownames(o$y)),
        eligible=n%in%usable && allowed)
    }
    usable<-usable[vapply(usable,function(n)o$tissue%in%c("Root/Shoot","Combined",refs[[n]]$tissue[1]),TRUE)]
    if(length(usable)>=3L){
      scores<-halves[paste0(usable,"__Up"),,drop=FALSE]-halves[paste0(usable,"__Down"),,drop=FALSE]
      rownames(scores)<-usable
      z<-extension_score_fit(scores,o);z$reference<-z$feature
      rt<-vapply(refs,function(d)d$tissue[1],"");z$reference_tissue<-unname(rt[z$reference])
      if(o$dataset=="Array1")z<-z[startsWith(z$contrast,z$reference_tissue),]
      results[[o$dataset]]<-z
      scores_all[[o$dataset]]<-data.frame(dataset=o$dataset,reference=rep(rownames(scores),ncol(scores)),
        sample=rep(colnames(scores),each=nrow(scores)),score=as.vector(scores))
    }
    for(co in o$selected_contrasts)for(n in names(refs)){
      d<-refs[[n]];tissue<-if(o$dataset=="Array1")sub("_.*","",co) else o$tissue
      if(!tissue%in%c("Combined",d$tissue[1]))next
      a<-o$table[o$table$contrast==co,];k<-match(a$feature,d$feature);ok<-!is.na(k)
      if(sum(ok)<100)next
      x<-a$log2FC[ok];y<-d$log2FC[k[ok]]
      correlations[[length(correlations)+1L]]<-data.frame(dataset=o$dataset,contrast=co,reference=n,
        reference_tissue=d$tissue[1],shared_loci=length(x),Pearson=cor(x,y),
        Spearman=cor(x,y,method="spearman"),direction_agreement=mean(sign(x)==sign(y)),
        interpretation="Descriptive across loci; loci are not biological replicates; no correlation P-value")
    }
  }
  d<-do.call(rbind,results);d$q_global<-p.adjust(d$P,"BH")
  list(table=d,coverage=do.call(rbind,coverage),members=do.call(rbind,members),
    correlations=do.call(rbind,correlations),scores=do.call(rbind,scores_all))
}

extension_markers <- function(markers,objects) {
  s<-extension_settings();sets<-lapply(split(markers$locus,markers$set),unique)
  annotation<-unique(markers[c("set","category")]);stopifnot(!anyDuplicated(annotation$set))
  rows<-scores_all<-coverage<-list()
  for(o in objects)for(mode in c("primary","without_focal_candidates")){
    ss<-sets
    if(mode!="primary")ss<-lapply(ss,setdiff,y=c(unname(integration_settings()$candidate_loci),
      "AT3G48750","AT5G44720","AT1G18390"))
    scores<-extension_rank_scores(o$y,ss,s$marker_min_members)
    category<-annotation$category[match(names(ss),annotation$set)]
    applicable<-if(o$tissue%in%c("Root/Shoot","Combined"))rep(TRUE,length(ss)) else
      if(o$tissue=="Root")category!="shoot_identity" else !category%in%c("root_development","root_identity")
    cov<-data.frame(dataset=o$dataset,mode=mode,set=names(ss),reference_members=lengths(ss),
      measured_members=vapply(ss,function(g)sum(g%in%rownames(o$y)),0L),
      eligible=names(ss)%in%rownames(scores)&applicable,anatomically_applicable=applicable)
    coverage[[paste(o$dataset,mode)]]<-cov
    z<-extension_score_fit(scores,o);names(z)[names(z)=="feature"]<-"set"
    category<-annotation$category[match(z$set,annotation$set)]
    applicable<-if(o$dataset=="Array1")ifelse(startsWith(z$contrast,"Root_"),category!="shoot_identity",
      !category%in%c("root_development","root_identity")) else cov$anatomically_applicable[match(z$set,cov$set)]
    z<-z[applicable,]
    z$mode<-mode;rows[[paste(o$dataset,mode)]]<-z
    scores_all[[paste(o$dataset,mode)]]<-data.frame(dataset=o$dataset,mode=mode,set=rep(rownames(scores),ncol(scores)),
      sample=rep(colnames(scores),each=nrow(scores)),score=as.vector(scores))
  }
  d<-do.call(rbind,rows);d$q_global<-ave(d$P,d$mode,FUN=function(p)p.adjust(p,"BH"))
  d$category<-annotation$category[match(d$set,annotation$set)]
  list(table=d,coverage=do.call(rbind,coverage),scores=do.call(rbind,scores_all))
}

extension_tf <- function(dap,objects) {
  s<-extension_settings();sets<-lapply(split(dap$edges$target,dap$edges$TF),unique)
  labels<-unique(dap$edges[c("TF","TF_name","TF_family")]);labels<-labels[!duplicated(labels$TF),]
  rows<-coverage<-list()
  for(o in objects){
    measured<-lapply(sets,intersect,y=rownames(o$y))
    valid<-lengths(measured)>=s$TF_min_members & lengths(measured)<=s$TF_max_universe_fraction*nrow(o$y)
    coverage[[o$dataset]]<-data.frame(dataset=o$dataset,TF=names(sets),
      reference_targets=lengths(sets),measured_targets=lengths(measured),universe=nrow(o$y),eligible=valid)
    indices<-lapply(measured[valid],function(g)match(g,rownames(o$y)))
    for(co in o$selected_contrasts){
      a<-o$table[o$table$contrast==co,];a<-a[match(rownames(o$y),a$feature),]
      stopifnot(!anyNA(a$feature),all(is.finite(a$log2FC/a$SE)))
      for(mode in names(s$camera_correlations)){
        z<-limma::cameraPR(a$log2FC/a$SE,index=indices,use.ranks=TRUE,
          inter.gene.cor=s$camera_correlations[[mode]],sort=FALSE)
        rows[[length(rows)+1L]]<-data.frame(dataset=o$dataset,contrast=co,mode=mode,
          TF=rownames(z),NGenes=z$NGenes,Direction=z$Direction,P=z$PValue,
          inter_gene_correlation=s$camera_correlations[[mode]],assay=o$assay,tissue=o$tissue)
      }
    }
  }
  d<-do.call(rbind,rows)
  d$q_global<-ave(d$P,d$mode,FUN=function(p)p.adjust(p,"BH"))
  d$q_family<-ave(d$P,interaction(d$mode,d$dataset,d$contrast,drop=TRUE),FUN=function(p)p.adjust(p,"BH"))
  d$TF_name<-labels$TF_name[match(d$TF,labels$TF)];d$TF_family<-labels$TF_family[match(d$TF,labels$TF)]
  list(table=d,coverage=do.call(rbind,coverage))
}
