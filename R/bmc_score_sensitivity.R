# Ordinary residual-variance sensitivity for biological-pool signature scores.
# It preserves each primary design, paired-score covariance and complete test
# family, and changes only empirical-Bayes variance borrowing.
bmc_score_design <- function(metadata) {
 m<-metadata
 if(unique(m$dataset)=='Array1'){
  m$cell<-factor(paste0(m$tissue,m$salt_mM),levels=unlist(lapply(c('Root','Shoot'),function(t)paste0(t,c(0,50,75,100)))))
  x<-model.matrix(~0+cell,m);colnames(x)<-levels(m$cell)
  cm<-matrix(0,8,6,dimnames=list(colnames(x),unlist(lapply(c('Root','Shoot'),function(t)paste0(t,'_',c('50v0','75v0','100v0'))))))
  for(t in c('Root','Shoot'))for(d in c(50,75,100))cm[paste0(t,c(d,0)),paste0(t,'_',d,'v0')]<-c(1,-1)
  block<-m$pool
 }else{
  m$cell<-factor(paste0(ifelse(m$genotype=='No-0','WT','OE'),ifelse(m$salt_mM==0,'0','S')),levels=c('WT0','WTS','OE0','OES'))
  x<-model.matrix(~0+cell,m);colnames(x)<-levels(m$cell)
  cm<-factorial_contrasts(x);block<-NULL
 }
 list(design=x,contrasts=cm,block=block)
}

bmc_unmoderated_scores <- function(scores,primary,metadata,feature_column) {
 stopifnot(all(c('feature','contrast','score_difference','SE','P','q_global','q_BY')%in%names(primary)))
 result<-list()
 modes<-if('mode'%in%names(scores))unique(scores$mode) else 'primary'
 if(!'mode'%in%names(scores))scores$mode<-'primary'
 if(!'mode'%in%names(primary))primary$mode<-'primary'
 for(ds in unique(scores$dataset))for(mode in modes){
  s<-scores[scores$dataset==ds & scores$mode==mode,]
  p<-primary[primary$dataset==ds & primary$mode==mode,]
  m<-metadata[metadata$dataset==ds,]
  if(!nrow(p))next
  features<-unique(s[[feature_column]])
  stopifnot(nrow(s)==length(features)*nrow(m),
   identical(s$sample,rep(m$sample,each=length(features))),
   identical(s[[feature_column]],rep(features,nrow(m))))
  y<-matrix(s$score,nrow=length(features),dimnames=list(features,m$sample))
  des<-bmc_score_design(m)
  rho<-if(is.null(des$block))NULL else unique(p$pair_correlation)
  if(!is.null(rho))stopifnot(length(rho)==1L,is.finite(rho))
  fit<-if(is.null(des$block))limma::lmFit(y,des$design) else limma::lmFit(y,des$design,block=des$block,correlation=rho)
  cf<-limma::contrasts.fit(fit,des$contrasts)
  rows<-lapply(seq_len(ncol(des$contrasts)),function(j){
   se<-cf$stdev.unscaled[,j]*cf$sigma;mu<-cf$coefficients[,j]
   critical<-qt(.975,cf$df.residual)
   data.frame(dataset=ds,mode=mode,feature=rownames(y),contrast=colnames(des$contrasts)[j],
    score_difference=mu,ordinary_SE=se,ordinary_lower95=mu-critical*se,ordinary_upper95=mu+critical*se,
    ordinary_P=2*pt(-abs(mu/se),df=cf$df.residual),ordinary_df=cf$df.residual)
  })
  a<-do.call(rbind,rows);key<-function(z)paste(z$feature,z$contrast,sep='::')
  a<-a[match(key(p),key(a)),]
  stopifnot(nrow(a)==nrow(p),nrow(a)>0L,!anyNA(a$feature),max(abs(a$score_difference-p$score_difference))<1e-10,
   all(is.finite(a$ordinary_P)))
  a$primary_SE<-p$SE;a$primary_lower95<-p$lower95;a$primary_upper95<-p$upper95
  a$primary_P<-p$P;a$primary_global_BH<-p$q_global;a$primary_global_BY<-p$q_BY
  result[[length(result)+1L]]<-a
 }
 d<-do.call(rbind,result)
 d$ordinary_global_BH<-ave(d$ordinary_P,d$mode,FUN=function(p)p.adjust(p,'BH'))
 d$ordinary_global_BY<-ave(d$ordinary_P,d$mode,FUN=function(p)p.adjust(p,'BY'))
 d
}
