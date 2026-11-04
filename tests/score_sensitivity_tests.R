script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),'..','scripts','bootstrap.R'))
set.seed(583)
m<-data.frame(dataset='test',sample=paste0('s',1:8),genotype=rep(c('No-0','NIT2'),each=4),salt_mM=rep(rep(c(0,100),each=2),2))
des<-bmc_score_design(m);y<-matrix(rnorm(80),10,8,dimnames=list(paste0('g',1:10),m$sample))
p<-moderated_contrasts(y,des$design,des$contrasts)$table
names(p)[names(p)=='log2FC']<-'score_difference'
p$dataset<-'test';p$q_global<-p.adjust(p$P,'BH');p$q_BY<-p.adjust(p$P,'BY')
s<-data.frame(dataset='test',reference=rep(rownames(y),ncol(y)),sample=rep(m$sample,each=nrow(y)),score=as.vector(y))
ans<-bmc_unmoderated_scores(s,p,m,'reference')
fit<-lm(y[1,]~0+des$design)
estimate<-coef(fit)[3]-coef(fit)[1]
variance<-sum(vcov(fit)[c(1,3),c(1,3)]*matrix(c(1,-1,-1,1),2))
a<-ans[ans$feature=='g1'&ans$contrast=='baseline_genotype',]
stopifnot(nrow(ans)==70L,nrow(a)==1L,abs(a$score_difference-estimate)<1e-12,abs(a$ordinary_SE-sqrt(variance))<1e-12,
 abs(a$ordinary_P-2*pt(-abs(estimate/sqrt(variance)),4))<1e-12,all(ans$ordinary_df==4))
cat('Ordinary score effects, variance and P values agree with independent lm calculation\n')
