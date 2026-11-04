# Residual-correlation sensitivity for the already defined DAP-seq target sets.
# The feature statistics and target universes are held fixed. Correlation is
# estimated from residual coordinates after whitening paired harvests, and is
# not selected to improve enrichment significance.
bmc_residual_coordinates <- function(y,design,block=NULL,correlation=NULL) {
 stopifnot(all(is.finite(y)),ncol(y)==nrow(design))
 v<-diag(ncol(y))
 if(!is.null(block)){
  stopifnot(length(block)==ncol(y),length(correlation)==1L,is.finite(correlation))
  v[outer(as.character(block),as.character(block),`==`)&row(v)!=col(v)]<-correlation
 }
 a<-solve(chol(v))
 x<-crossprod(a,design);q<-qr.Q(qr(x),complete=TRUE)
 rank<-qr(x)$rank;df<-nrow(x)-rank
 if(df<2L)stop('Insufficient residual degrees of freedom')
 r<-(y%*%a)%*%q[,seq.int(rank+1L,nrow(q)),drop=FALSE]
 ss<-rowSums(r*r)
 if(any(!is.finite(ss)|ss<=0))stop('Zero or invalid residual variance')
 list(coordinates=r/sqrt(ss),df=df,covariance=v,whitener=a)
}

bmc_mean_residual_correlation <- function(coordinates,index) {
 n<-length(index);if(n<2L)stop('At least two targets are required')
 # Equivalent to the off-diagonal mean of the residual correlation matrix,
 # without forming a potentially large gene-by-gene matrix.
 (sum(colSums(coordinates[index,,drop=FALSE])^2)-n)/(n*(n-1))
}

bmc_tf_correlation <- function(dap,objects) {
 s<-extension_settings();sets<-lapply(split(dap$edges$target,dap$edges$TF),unique)
 labels<-unique(dap$edges[c('TF','TF_name','TF_family')]);labels<-labels[!duplicated(labels$TF),]
 rows<-audit<-list()
 for(o in objects){
  measured<-lapply(sets,intersect,y=rownames(o$y))
  valid<-lengths(measured)>=s$TF_min_members & lengths(measured)<=s$TF_max_universe_fraction*nrow(o$y)
  indices<-lapply(measured[valid],function(g)match(g,rownames(o$y)))
  r<-bmc_residual_coordinates(o$y,o$design,o$block,o$correlation)
  rho<-vapply(indices,function(i)bmc_mean_residual_correlation(r$coordinates,i),0.0)
  used<-pmin(1-1e-8,pmax(0,rho))
  audit[[o$dataset]]<-data.frame(dataset=o$dataset,TF=names(indices),targets=lengths(indices),
   universe=nrow(o$y),residual_df=r$df,raw_residual_correlation=rho,used_residual_correlation=used,
   harvested_pair_correlation=if(is.null(o$correlation))NA_real_ else o$correlation)
  for(co in o$selected_contrasts){
   a<-o$table[o$table$contrast==co,];a<-a[match(rownames(o$y),a$feature),]
   stopifnot(!anyNA(a$feature),all(is.finite(a$log2FC/a$SE)))
   z<-limma::cameraPR(a$log2FC/a$SE,index=indices,use.ranks=TRUE,inter.gene.cor=used,sort=FALSE)
   rows[[length(rows)+1L]]<-data.frame(dataset=o$dataset,contrast=co,TF=rownames(z),
    NGenes=z$NGenes,Direction=z$Direction,P=z$PValue,
    residual_correlation=rho[rownames(z)],used_correlation=used[rownames(z)],
    residual_df=r$df,assay=o$assay,tissue=o$tissue)
  }
 }
 d<-do.call(rbind,rows);d$q_global<-p.adjust(d$P,'BH');d$q_BY<-p.adjust(d$P,'BY')
 d$TF_name<-labels$TF_name[match(d$TF,labels$TF)];d$TF_family<-labels$TF_family[match(d$TF,labels$TF)]
 list(TF_residual_correlation=do.call(rbind,audit),TF_estimated_correlation_tests=d)
}
