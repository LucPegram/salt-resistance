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
  used<-setNames(pmin(1-1e-8,pmax(0,rho)),names(rho))
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
 d<-do.call(rbind,rows)
 stopifnot(all(is.finite(d$used_correlation)),all(is.finite(d$residual_correlation)))
 d$q_global<-p.adjust(d$P,'BH');d$q_BY<-p.adjust(d$P,'BY')
 d$TF_name<-labels$TF_name[match(d$TF,labels$TF)];d$TF_family<-labels$TF_family[match(d$TF,labels$TF)]
 list(TF_residual_correlation=do.call(rbind,audit),TF_estimated_correlation_tests=d)
}

# Show every setting for the union of findings across the complete original
# families. Never recalculate a multiple-testing correction on selected rows.
bmc_tf_comparison <- function(fixed,estimated) {
 key<-function(d)paste(d$dataset,d$contrast,d$TF,sep='::')
 stopifnot(setequal(unique(fixed$mode),c('primary','sensitivity')),
  !anyDuplicated(paste(key(fixed),fixed$mode)),!anyDuplicated(key(estimated)),
  setequal(key(fixed[fixed$mode=='primary',]),key(estimated)),
  setequal(key(fixed[fixed$mode=='sensitivity',]),key(estimated)))
 selected<-unique(c(key(fixed)[fixed$q_global<.05],key(estimated)[estimated$q_global<.05]))
 columns<-c('TF','TF_name','dataset','contrast')
 a<-data.frame(fixed[columns],mode=ifelse(fixed$mode=='primary','Fixed 0.01','Fixed 0.05'),
  targets=fixed$NGenes,correlation=fixed$inter_gene_correlation,
  BH=fixed$q_global,BY=fixed$q_BY,Direction=fixed$Direction)
 b<-data.frame(estimated[columns],mode='Estimated residual',targets=estimated$NGenes,
  correlation=estimated$used_correlation,BH=estimated$q_global,BY=estimated$q_BY,
  Direction=estimated$Direction)
 d<-rbind(a[key(fixed)%in%selected,],b[key(estimated)%in%selected,])
 if(!nrow(d))return(d)
 stopifnot(!anyDuplicated(paste(key(d),d$mode)),all(table(key(d))==3L))
 best<-tapply(d$BH,key(d),min);keys<-names(best)[order(best,names(best))]
 d[order(match(key(d),keys),match(d$mode,c('Fixed 0.01','Estimated residual','Fixed 0.05'))),]
}

bmc_tf_comparison_figure <- function(d,path) {
 stopifnot(nrow(d)>0,all(is.finite(d$BH)),all(d$BH>=0&d$BH<=1))
 key<-paste(d$dataset,d$contrast,d$TF,sep='::')
 keys<-unique(key);n<-length(keys)
 labels<-vapply(keys,function(k){
  x<-d[match(k,key),];name<-if(x$TF_name==x$TF)x$TF else paste(x$TF_name,x$TF,sep=' / ')
  if(length(unique(paste(d$dataset,d$contrast)))>1)paste(name,x$dataset,x$contrast,sep=' | ') else name
 },'')
 modes<-c('Fixed 0.01','Estimated residual','Fixed 0.05')
 colors<-c('#0072B2','#8B4D9C','#D55E00');shapes<-c(16,17,15);offset<-c(.19,0,-.19)
 publication_figure(path,9.5,max(4.8,.45*n+1.7),function(){
  layout(matrix(c(1,2),1,2),widths=c(2.25,1))
  par(mar=c(4.4,14.8,4.1,1.2),cex=1)
  score<- -log10(pmax(d$BH,1e-300));xmax<-max(2,max(score)*1.12)
  plot(NA,xlim=c(-.025*xmax,xmax),ylim=c(.45,n+1.05),axes=FALSE,xlab='-log10 global BH q',ylab='',xaxs='i')
  abline(h=seq_len(n),col='#E6EBEF',lwd=.7);abline(v=-log10(.05),col='#7B8794',lty=2)
  axis(1);axis(2,at=n:1,labels=labels,las=1,tick=FALSE,cex.axis=.86)
  for(j in seq_along(modes)){
   z<-d[d$mode==modes[j],];ii<-match(paste(z$dataset,z$contrast,z$TF,sep='::'),keys)
   points(-log10(pmax(z$BH,1e-300)),n-ii+1+offset[j],pch=shapes[j],col=colors[j],cex=.95)
  }
  title(main='A  Target enrichment',adj=0,cex.main=1.1)
  legend('top',inset=c(0,.015),legend=modes,pch=shapes,col=colors,bty='n',horiz=TRUE,cex=.78)
  par(mar=c(4.4,.8,4.1,1),cex=1)
  z<-d[d$mode=='Estimated residual',];ii<-match(paste(z$dataset,z$contrast,z$TF,sep='::'),keys)
  xmax<-max(.055,max(z$correlation)*1.12)
  plot(NA,xlim=c(0,xmax),ylim=c(.45,n+1.05),axes=FALSE,xlab='Residual target correlation',ylab='',xaxs='i')
  abline(h=seq_len(n),col='#E6EBEF',lwd=.7);abline(v=c(.01,.05),col='#7B8794',lty=2)
  axis(1);points(z$correlation,n-ii+1,pch=17,col=colors[2],cex=1)
  title(main='B  Estimated correlation',adj=0,cex.main=1.1)
  publication_caption('Every TF-by-assay-by-contrast finding with global BH q < 0.05 under any of the three correlation settings is displayed under all three settings. Selection follows the complete families; no correction is recalculated on these rows.')
  publication_caption('Panel A retains each setting\'s original full-family global BH q values. The dashed vertical line is q = 0.05. Symbols are vertically offset to distinguish coincident values. Global BY values and target counts are in the companion table.')
  publication_caption('Panel B shows the nonnegative estimated average residual target correlation after whitening the fitted harvested-pair covariance; dashed lines mark the two fixed assumptions, 0.01 and 0.05. The settings are sensitivities on the same data, not independent confirmations; no setting establishes TF activity.')
 })
}
