script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),'..','scripts','bootstrap.R'))
set.seed(418)
x<-model.matrix(~factor(rep(1:4,each=3)))
y<-matrix(rnorm(60*12),60,12)
z<-bmc_residual_coordinates(y,x)
i<-c(1,3,5,6,8)
reference<-limma::interGeneCorrelation(y[i,],x)$correlation
stopifnot(abs(bmc_mean_residual_correlation(z$coordinates,i)-reference)<1e-12)
b<-rep(1:6,2);rho<-.35
z<-bmc_residual_coordinates(y,x,b,rho)
a<-z$whitener
stopifnot(max(abs(crossprod(a,z$covariance%*%a)-diag(12)))<1e-12)
reference<-limma::interGeneCorrelation((y%*%a)[i,],crossprod(a,x))$correlation
stopifnot(abs(bmc_mean_residual_correlation(z$coordinates,i)-reference)<1e-12)
source(file.path(project_root,'config','extensions.R'))
rownames(y)<-paste0('g',1:60)
dap<-list(edges=data.frame(TF=rep(c('tf1','tf2'),each=10),target=paste0('g',1:20),
 TF_name=rep(c('One','Two'),each=10),TF_family='test'))
o<-list(y=y,design=x,block=b,correlation=rho,dataset='synthetic',assay='RNA',tissue='paired',
 selected_contrasts='test',table=data.frame(feature=rownames(y),contrast='test',log2FC=rnorm(60),SE=1))
ans<-bmc_tf_correlation(dap,list(o))
stopifnot(all(is.finite(ans$TF_estimated_correlation_tests$used_correlation)),
 identical(ans$TF_estimated_correlation_tests$used_correlation,
  ans$TF_residual_correlation$used_residual_correlation))
cat('Residual correlation agrees with limma; paired covariance whitening verified\n')

# A finding selected under one setting must retain all other settings, even when
# their q values are nonsignificant. Values are copied from full-family results.
fixed<-data.frame(TF=rep(c('A','B','C'),2),TF_name=rep(c('Alpha','Beta','Gamma'),2),
 dataset='test',contrast='salt',mode=rep(c('primary','sensitivity'),each=3),
 NGenes=12,inter_gene_correlation=rep(c(.01,.05),each=3),
 q_global=c(.01,.8,.9,1,1,1),q_BY=c(.04,1,1,1,1,1),Direction='Up')
estimated<-data.frame(TF=c('C','A','B'),TF_name=c('Gamma','Alpha','Beta'),
 dataset='test',contrast='salt',NGenes=12,used_correlation=c(.01,.03,.005),
 q_global=c(.9,.4,.02),q_BY=c(1,1,.08),Direction='Up')
comparison<-bmc_tf_comparison(fixed,estimated)
stopifnot(nrow(comparison)==6L,setequal(comparison$TF,c('A','B')),
 all(table(comparison$TF)==3L),
 comparison$BH[comparison$TF=='A'&comparison$mode=='Estimated residual']==.4,
 comparison$BH[comparison$TF=='B'&comparison$mode=='Fixed 0.01']==.8)
cat('All-setting display retains the union of full-family findings and unchanged q values\n')
