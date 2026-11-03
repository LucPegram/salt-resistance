# New fits use external measurement matrices. Legacy result workbooks are never
# treated as inferential inputs. No experimental values are embedded here.
.revision_cache <- new.env(parent = emptyenv())

revision_cached <- function(key, inputs, expression) {
  model_files<-file.path(project_root,c("R/revision_models.R","R/inputs.R",
    "R/revision_enrichment.R","config/analysis.R","renv.lock"))
  fingerprint <- digest::digest(list(key,inputs, analysis_settings,
    vapply(model_files,function(p)digest::digest(file=p,algo="sha256"),""),
    vapply(c("limma","statmod","proDA","GO.db"),function(p)as.character(packageVersion(p)),"")), algo = "sha256")
  key <- paste(key, fingerprint, sep = ":")
  directory<-Sys.getenv("SALT_MODEL_CACHE")
  if(nzchar(directory)){
    directory<-normalizePath(directory,mustWork=FALSE)
    roots<-c(project_root,Sys.getenv("SALT_DATA_ROOT"))
    for(root in roots[nzchar(roots)]){
      root<-normalizePath(root,mustWork=TRUE)
      if(identical(directory,root)||startsWith(directory,paste0(root,"/")))
        stop("Keep model caches outside code and original inputs")
    }
  }
  cachefile<-if(nzchar(directory))file.path(directory,paste0(fingerprint,".rds")) else ""
  if(!exists(key,.revision_cache,inherits=FALSE)){
    if(nzchar(cachefile)&&file.exists(cachefile)){
      saved<-readRDS(cachefile)
      if(!identical(saved$fingerprint,fingerprint)||
         !identical(saved$sha256,digest::digest(saved$value,algo="sha256")))stop("Invalid model cache")
      value<-saved$value
    }else{
      value<-force(expression)
      if(nzchar(cachefile)){
        dir.create(directory,recursive=TRUE,showWarnings=FALSE)
        saveRDS(list(fingerprint=fingerprint,value=value,
          sha256=digest::digest(value,algo="sha256")),cachefile)
      }
    }
    assign(key,value,.revision_cache)
  }
  get(key, .revision_cache, inherits = FALSE)
}

factorial_metadata <- function(x) {
  if (x$spec$assay == "microarray") {
    s <- x$samples
    genotype <- ifelse(startsWith(s, "OE"), "NIT2", "No-0")
    salt <- ifelse(grepl(" PN ", s), 0L, 100L)
  } else {
    s <- x$samples; genotype <- x$metadata$genotype; salt <- x$metadata$salt_mM
  }
  data.frame(sample = s, genotype = genotype, salt_mM = salt,
    cell = factor(paste0(ifelse(genotype == "No-0", "WT", "OE"),
      ifelse(salt == 0, "0", "S")), levels = c("WT0", "WTS", "OE0", "OES")),
    stringsAsFactors = FALSE)
}

factorial_contrasts <- function(design) {
  limma::makeContrasts(
    baseline_genotype = OE0-WT0, saline_genotype = OES-WTS,
    salt_WT = WTS-WT0, salt_OE = OES-OE0,
    interaction = (OES-OE0)-(WTS-WT0),
    average_genotype = (OE0+OES-WT0-WTS)/2,
    average_salt = (WTS+OES-WT0-OE0)/2, levels = design)
}

moderated_contrasts <- function(y, design, contrasts, block = NULL, correlation = NULL) {
  stopifnot(ncol(y) == nrow(design), nrow(contrasts) == ncol(design))
  fit <- if (is.null(block)) limma::lmFit(y, design) else
    limma::lmFit(y, design, block = block, correlation = correlation)
  cf <- limma::eBayes(limma::contrasts.fit(fit, contrasts), trend = TRUE, robust = TRUE)
  tables <- lapply(seq_len(ncol(contrasts)), function(j) {
    se <- cf$stdev.unscaled[,j] * sqrt(cf$s2.post)
    critical <- stats::qt(.975, cf$df.total)
    data.frame(feature = rownames(y), contrast = colnames(contrasts)[j],
      log2FC = cf$coefficients[,j], SE = se,
      lower95 = cf$coefficients[,j] - critical*se,
      upper95 = cf$coefficients[,j] + critical*se,
      P = cf$p.value[,j], q = p.adjust(cf$p.value[,j], "BH"),
      df = cf$df.total, stringsAsFactors = FALSE)
  })
  table <- do.call(rbind, tables)
  # BH across features within each contrast, then Bonferroni across the seven
  # planned contrast families. This conservative second criterion is nested.
  table$q_all_contrasts <- pmin(1,ncol(contrasts)*table$q)
  list(table = table, fit = fit, contrast_fit = cf)
}

revision_genotype <- function(x) revision_cached("genotype", x, {
  y <- collapse_microarray(x)
  stopifnot(all(is.finite(y)))
  m <- factorial_metadata(x)
  design <- model.matrix(~0+cell, m); colnames(design) <- levels(m$cell)
  out <- moderated_contrasts(y, design, factorial_contrasts(design))
  out$y <- y; out$metadata <- m; out
})

revision_dose <- function(inputs) revision_cached("dose", inputs, {
  xs <- inputs[c("microarray_root", "microarray_shoot")]
  matrices <- lapply(xs, collapse_microarray)
  stopifnot(identical(rownames(matrices[[1]]), rownames(matrices[[2]])))
  y <- do.call(cbind, matrices)
  m <- do.call(rbind, lapply(xs, function(x) data.frame(
    sample = x$samples, tissue = x$spec$tissue,
    dose = as.integer(ifelse(startsWith(x$samples,"PN"), "0", sub(" .*", "",x$samples))),
    replicate = as.integer(sub(".*B", "",x$samples)))))
  m$pool <- factor(paste(m$dose,m$replicate,sep="_"))
  m$cell <- factor(paste0(m$tissue,m$dose),
    levels = unlist(lapply(c("Root","Shoot"),function(t)paste0(t,c(0,50,75,100)))))
  design <- model.matrix(~0+cell,m);colnames(design) <- levels(m$cell)
  rho <- limma::duplicateCorrelation(y, design, block=m$pool)$consensus.correlation
  cm <- limma::makeContrasts(
    root50=Root50-Root0, root75=Root75-Root0, root100=Root100-Root0,
    shoot50=Shoot50-Shoot0, shoot75=Shoot75-Shoot0, shoot100=Shoot100-Shoot0,
    tissue50=(Shoot50-Shoot0)-(Root50-Root0),
    tissue75=(Shoot75-Shoot0)-(Root75-Root0),
    tissue100=(Shoot100-Shoot0)-(Root100-Root0),levels=design)
  out <- moderated_contrasts(y,design,cm,block=m$pool,correlation=rho)
  ftests <- lapply(list(root=1:3,shoot=4:6,tissue_by_dose=7:9),function(j){
    z<-limma::topTable(out$contrast_fit,coef=j,number=Inf,sort.by="none")
    data.frame(feature=rownames(z),F=z$F,P=z$P.Value,q=z$adj.P.Val)
  })
  names(ftests)<-c("root","shoot","tissue_by_dose")
  out$joint<-do.call(rbind,lapply(names(ftests),function(n){z<-ftests[[n]];z$test<-n;z}))
  out$y<-y;out$metadata<-m;out$pair_correlation<-rho;out
})

protein_cell_counts <- function(y, m) vapply(levels(m$cell),function(g)
  rowSums(!is.na(y[,m$cell==g,drop=FALSE])), numeric(nrow(y)))

revision_protein <- function(x) revision_cached(paste0("protein_",x$spec$tissue),x,{
  keep <- proteomics_keep(x,analysis_settings)
  d <- x$data[keep,,drop=FALSE]
  y <- log2(x$values[keep,,drop=FALSE]);rownames(y)<-d$Protein.Group
  m <- factorial_metadata(x)
  counts <- protein_cell_counts(y,m)
  design <- model.matrix(~0+cell,m);colnames(design)<-levels(m$cell)
  cm <- factorial_contrasts(design)
  selection <- list(observed_min2=apply(counts,1,min)>=2,
    observed_available=apply(counts,1,min)>=1 & rowSums(counts)>=8,
    complete_case=rowSums(counts)==12)
  results<-lapply(names(selection),function(method){
    out<-moderated_contrasts(y[selection[[method]],,drop=FALSE],design,cm)$table
    out$method<-method;out
  })
  # Normalization sensitivity only; supplied measurements remain unchanged.
  yn<-sweep(y,2,apply(y,2,median,na.rm=TRUE),"-")
  out<-moderated_contrasts(yn[selection$observed_min2,,drop=FALSE],design,cm)$table
  out$method<-"observed_min2_median_centered";results[[length(results)+1L]]<-out
  if(!requireNamespace("proDA",quietly=TRUE))stop("Restore proDA before revision analyses")
  warnings<-character()
  pd<-withCallingHandlers(proDA::proDA(y,design=design,verbose=FALSE),
    warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  pdresults<-lapply(seq_len(ncol(cm)),function(j){
    contrast<-paste(paste0("(",cm[,j],")*",rownames(cm)),collapse=" + ")
    z<-proDA::test_diff(pd,contrast=contrast,sort_by=NULL)
    critical<-qt(.975,z$df)
    data.frame(feature=z$name,contrast=colnames(cm)[j],log2FC=z$diff,SE=z$se,
      lower95=z$diff-critical*z$se,upper95=z$diff+critical*z$se,P=z$pval,
      q=z$adj_pval,df=z$df,method="proDA_dropout",stringsAsFactors=FALSE)
  })
  pdtable<-do.call(rbind,pdresults)
  pdtable$q_all_contrasts<-pmin(1,ncol(cm)*pdtable$q)
  # Same field order as limma tables.
  results[[length(results)+1L]]<-pdtable[,names(results[[1]])]
  table<-do.call(rbind,results)
  table$tissue<-x$spec$tissue
  table$Genes<-d$Genes[match(table$feature,d$Protein.Group)]
  list(y=y,metadata=m,counts=counts,annotations=d,table=table,
    proDA_warnings=unique(warnings),proDA_convergence=proDA::convergence(pd),
    selection=selection)
})

# Disjoint WT-control subsets avoid shared-denominator correlation. Their
# similarity range is descriptive; molecular rows are not biological replicates.
prestress_splits <- function(y,m,label) {
  wt0<-which(m$cell=="WT0");oe0<-which(m$cell=="OE0");wts<-which(m$cell=="WTS")
  summarize<-function(a,b,kind,split){
    ok<-is.finite(a)&is.finite(b)
    a<-a[ok];b<-b[ok]
    data.frame(dataset=label,comparison=kind,split=split,features=length(a),
      Pearson=cor(a,b),Spearman=cor(a,b,method="spearman"),
      cosine=sum(a*b)/sqrt(sum(a*a)*sum(b*b)),
      direction_agreement=mean(sign(a)==sign(b)),
      median_abs_baseline=median(abs(a)),median_abs_WT_salt=median(abs(b)))
  }
  fullB<-rowMeans(y[,oe0,drop=FALSE])-rowMeans(y[,wt0,drop=FALSE])
  fullS<-rowMeans(y[,wts,drop=FALSE])-rowMeans(y[,wt0,drop=FALSE])
  rows<-list(summarize(fullB,fullS,"shared_controls","all"))
  # All nonempty proper subsets, not a data-selected favorable split.
  for(k in seq_len(length(wt0)-1L))for(a in combn(wt0,k,simplify=FALSE)){
    b<-setdiff(wt0,a)
    A<-rowMeans(y[,oe0,drop=FALSE])-rowMeans(y[,a,drop=FALSE])
    B<-rowMeans(y[,wts,drop=FALSE])-rowMeans(y[,b,drop=FALSE])
    rows[[length(rows)+1L]]<-summarize(A,B,"disjoint_controls",
      paste(paste(a,collapse=","),paste(b,collapse=","),sep=" vs "))
  }
  do.call(rbind,rows)
}

build_revision_microarray <- function(inputs,settings,destination)
  write_table(revision_genotype(inputs$microarray_genotype)$table,destination)
build_revision_dose <- function(inputs,settings,destination)
  write_table(revision_dose(inputs)$joint,destination)
build_revision_protein <- function(inputs,settings,destination)
  write_table(do.call(rbind,lapply(inputs, function(x)revision_protein(x)$table)),destination)

revision_discoveries <- function(d, fields) {
  groups<-split(seq_len(nrow(d)),interaction(d[fields],drop=TRUE,lex.order=TRUE))
  do.call(rbind,lapply(groups,function(i){
    z<-d[i,,drop=FALSE];v<-z[1,fields,drop=FALSE]
    v$tested<-sum(!is.na(z$P));v$q_lt_0_05<-sum(z$q<.05,na.rm=TRUE)
    v$positive<-sum(z$q<.05 & z$log2FC>0,na.rm=TRUE)
    v$negative<-sum(z$q<.05 & z$log2FC<0,na.rm=TRUE)
    v$global_q_lt_0_05<-sum(z$q_all_contrasts<.05,na.rm=TRUE);v
  }))
}

build_revision_summary <- function(inputs,settings,destination) {
  rows<-list()
  if("microarray_genotype"%in%names(inputs)){
    z<-revision_discoveries(revision_genotype(inputs$microarray_genotype)$table,"contrast")
    z$tissue<-"Combined";z$method<-"limma_genotype";rows[[1]]<-z
  }
  for(x in Filter(function(x)x$spec$assay=="proteomics",inputs))
    rows[[length(rows)+1L]]<-revision_discoveries(revision_protein(x)$table,c("contrast","tissue","method"))
  columns<-c("tissue","method","contrast","tested","q_lt_0_05","positive","negative","global_q_lt_0_05")
  write_table(do.call(rbind,lapply(rows,function(x)x[columns])),destination)
}

build_revision_candidates <- function(inputs,settings,destination) {
  rows<-lapply(inputs,function(x){
    out<-revision_protein(x)
    hit<-out$table$Genes%in%c("NIT2","AHA2","NHX7")
    z<-out$table[hit,,drop=FALSE]
    z$Gene<-ifelse(z$Genes=="NHX7","SOS1",z$Genes)
    z
  })
  write_table(do.call(rbind,rows),destination)
}

build_revision_missingness <- function(inputs,settings,destination) {
  rows<-list()
  for(x in inputs){z<-revision_protein(x)
    for(j in seq_len(ncol(z$counts)))rows[[length(rows)+1L]]<-data.frame(
      tissue=x$spec$tissue,cell=colnames(z$counts)[j],retained=nrow(z$y),
      missing_fraction=mean(is.na(z$y[,z$metadata$cell==colnames(z$counts)[j],drop=FALSE])),
      completely_missing=sum(z$counts[,j]==0),single_observation=sum(z$counts[,j]==1),
      primary_eligible=sum(z$selection$observed_min2),complete=sum(z$selection$complete_case),
      proDA_warnings=paste(z$proDA_warnings,collapse=" | "))
  }
  write_table(do.call(rbind,rows),destination)
}

build_revision_identity <- function(inputs,settings,destination) {
  rows<-lapply(inputs,function(x){
    d<-x$data[x$data$Genes%in%c("NIT2","AHA2","NHX7"),,drop=FALSE]
    out<-d[,c("Protein.Group","Genes","N.Sequences","N.Proteotypic.Sequences")]
    out$tissue<-x$spec$tissue;out$accessions_in_group<-lengths(strsplit(out$Protein.Group,";",fixed=TRUE))
    out$complete_observed<-rowSums(!is.na(as.matrix(d[x$samples])))
    counts<-protein_cell_counts(as.matrix(d[x$samples]),factorial_metadata(x))
    out<-cbind(out,as.data.frame(counts))
    out
  })
  write_table(do.call(rbind,rows),destination)
}

build_revision_model_qc <- function(inputs,settings,destination) {
  d<-revision_dose(inputs)
  rows<-list(data.frame(dataset="paired_dose",diagnostic="consensus_root_shoot_correlation",
    value=as.character(d$pair_correlation)))
  for(x in Filter(function(x)x$spec$assay=="proteomics",inputs)){
    z<-revision_protein(x)
    rows[[length(rows)+1L]]<-data.frame(dataset=x$spec$tissue,
      diagnostic="proDA_convergence",value=paste(capture.output(str(z$proDA_convergence)),collapse=" "))
    rows[[length(rows)+1L]]<-data.frame(dataset=x$spec$tissue,
      diagnostic="proDA_warnings",value=paste(z$proDA_warnings,collapse=" | "))
  }
  write_table(do.call(rbind,rows),destination)
}

build_revision_prestress <- function(inputs,settings,destination) {
  rows<-list()
  if("microarray_genotype"%in%names(inputs)){
    z<-revision_genotype(inputs$microarray_genotype)
    rows[[1]]<-prestress_splits(z$y,z$metadata,"Genotype microarray")
  }
  for(x in Filter(function(x)x$spec$assay=="proteomics",inputs)){
    z<-revision_protein(x)
    # Use fully observed groups so both control splits contain observations.
    y<-z$y[z$selection$complete_case,,drop=FALSE]
    rows[[length(rows)+1L]]<-prestress_splits(y,z$metadata,paste(x$spec$tissue,"proteome complete cases"))
  }
  write_table(do.call(rbind,rows),destination)
}

revision_qc <- function(x) {
  y<-if(x$spec$assay=="microarray")collapse_microarray(x) else {
    k<-proteomics_keep(x,analysis_settings);log2(x$values[k,,drop=FALSE])}
  complete<-rowSums(is.finite(y))==ncol(y)
  pca<-prcomp(t(y[complete,,drop=FALSE]),center=TRUE,scale.=FALSE)
  fraction<-pca$sdev^2/sum(pca$sdev^2)
  data.frame(sample=x$samples,tissue=x$spec$tissue,assay=x$spec$assay,
    median=apply(y,2,median,na.rm=TRUE),Q25=apply(y,2,quantile,.25,na.rm=TRUE),
    Q75=apply(y,2,quantile,.75,na.rm=TRUE),missing_fraction=colMeans(is.na(y)),
    PCA_features=sum(complete),PC1=pca$x[,1],PC2=pca$x[,2],
    PC1_fraction=fraction[1],PC2_fraction=fraction[2],stringsAsFactors=FALSE)
}
build_revision_qc <- function(inputs,settings,destination)
  write_table(do.call(rbind,lapply(inputs,revision_qc)),destination)

build_revision_legacy_sensitivity <- function(inputs,settings,destination) {
  rows<-list()
  for(t in c("root","shoot")){
    new<-revision_protein(inputs[[paste0("proteomics_",t)]])$table
    old<-inputs[[paste0("saved_proteomics_",t)]]$data
    for(co in c(average_genotype="p_FDR_genotype",average_salt="p_FDR_treatment",interaction="p_FDR_interaction")){
      contrast<-names(c(average_genotype="p_FDR_genotype",average_salt="p_FDR_treatment",interaction="p_FDR_interaction"))[
        match(co,c("p_FDR_genotype","p_FDR_treatment","p_FDR_interaction"))]
      oldsig<-old$Protein.Group[checked_q(old,co)<.05]
      oldsig<-oldsig[!is.na(oldsig)]
      for(method in unique(new$method)){
        z<-new[new$method==method & new$contrast==contrast,,drop=FALSE]
        keys<-sub(";.*","",z$feature)
        if(anyDuplicated(keys))stop("Ambiguous legacy group match")
        nsig<-keys[!is.na(z$q)&z$q<.05]
        rows[[length(rows)+1L]]<-data.frame(tissue=t,contrast=contrast,method=method,
          old_significant=length(oldsig),old_significant_testable=sum(oldsig%in%keys),
          new_significant=length(nsig),overlap=sum(nsig%in%oldsig))
      }
    }
  }
  write_table(do.call(rbind,rows),destination)
}
