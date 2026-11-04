script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),'..','scripts','bootstrap.R'))
set.seed(20261104)
m<-expand.grid(replicate=1:3,dose=c(0,50,75,100),tissue=c('Root','Shoot'),stringsAsFactors=FALSE)
m$sample<-paste(m$tissue,m$dose,m$replicate)
y<-matrix(rnorm(200*24,sd=.05),200,24,dimnames=list(paste0('L',1:200),m$sample))
y<-sweep(y,2,ifelse(m$tissue=='Root',m$dose/100,m$dose/50),'+')
# Check against a separate cell-mean calculation after deleting a complete pair.
keep<-!(m$dose==75&m$replicate==2)
z<-bmc_fit_dose(y[,keep],m[keep,])
stopifnot(nrow(z$metadata)==22,nlevels(z$metadata$pool)==11,all(table(z$metadata$pool)==2))
for(co in c('Root_75v0','Shoot_75v0','Difference_75v0')){
 a<-z$table[z$table$contrast==co,]
 mu<-function(t,d)rowMeans(y[,keep & m$tissue==t & m$dose==d,drop=FALSE])
 expected<-switch(co,Root_75v0=mu('Root',75)-mu('Root',0),Shoot_75v0=mu('Shoot',75)-mu('Shoot',0),
  Difference_75v0=mu('Shoot',75)-mu('Shoot',0)-mu('Root',75)+mu('Root',0))
 stopifnot(max(abs(a$log2FC-expected))<1e-10,all(a$lower95<a$upper95))
}
# Repeated spots are technical measurements, not added biological replicates.
a<-matrix(c(1,3,9,5,7,11),3,2);mappings<-c('AT1G00010.1','AT1G00010.1','AT1G00020.1')
b<-bmc_collapse_logspots(a,mappings)
stopifnot(identical(as.numeric(b['AT1G00010',]),c(2,6)),ncol(b)==2)
cat('Paired omission, independent contrast arithmetic and spot-collapse checks passed\n')
