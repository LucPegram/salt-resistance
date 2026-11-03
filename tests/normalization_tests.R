script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),'..','scripts','bootstrap.R'))
# Synthetic continuous signals: verify quantile normalization against rank means,
# exclusion of control spots, and independence of root/shoot fitting groups.
set.seed(902)
n<-200L;regular<-seq_len(n)<=196L
samples<-unlist(lapply(c('PN','50N','75N','100N'),function(g)
 unlist(lapply(c('root','shoot'),function(t)paste(g,t,1:3)))))
x<-matrix(exp(rnorm(n*24,5,1)),nrow=n,dimnames=list(NULL,samples))
z<-list(values=x,samples=samples,probe=c(rep('duplicate',2),paste0('p',3:n)),
 control=c(rep(0L,196),1L,1L,-1L,-1L))
a<-normalization_variants(list(objects=list(dose=z)))
root<-a$matrices$Array1_Root;shoot<-a$matrices$Array1_Shoot
raw<-log2(x[regular,root$columns]);target<-rowMeans(apply(raw,2,sort))
expected<-apply(raw,2,function(v)target[rank(v,ties.method='first')])
stopifnot(max(abs(expected-root$methods$quantile[regular,]))<1e-12,
 nrow(root$methods$quantile)==n,all(is.finite(root$methods$quantile)),
 identical(root$columns,grep('root',samples)))
z$values[!regular,]<-z$values[!regular,]*1000
z$values[,grep('shoot',samples)]<-z$values[,grep('shoot',samples)]^1.3
b<-normalization_variants(list(objects=list(dose=z)))
stopifnot(identical(root$methods$quantile[regular,],b$matrices$Array1_Root$methods$quantile[regular,]),
 !isTRUE(all.equal(shoot$methods$quantile[regular,],b$matrices$Array1_Shoot$methods$quantile[regular,])))
# Ties share their averaged normalized ranks; no arbitrary spot ordering.
tied<-matrix(c(1,1,3,4,2,4,6,8),ncol=2)
qn<-limma::normalizeBetweenArrays(tied,method='quantile')
stopifnot(qn[1,1]==qn[2,1])
cat('Independent rank-mean, control exclusion, separate-tissue fitting and tie checks passed.\n')
