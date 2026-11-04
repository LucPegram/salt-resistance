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
cat('Residual correlation agrees with limma; paired covariance whitening verified\n')
