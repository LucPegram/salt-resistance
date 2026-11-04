# Prespecified response to the array-QC concern: all twelve harvested pairs,
# both frozen and leave-pair-out refitted quantile normalization. No sample is
# removed from the primary analysis. Loci and hypothesis families are fixed.
bmc_dose_metadata <- function(samples,tissue) data.frame(sample=samples,tissue=tissue,
 dose=as.integer(ifelse(startsWith(samples,'PN'),'0',sub(' .*','',samples))),
 replicate=as.integer(sub('.*B','',samples)),stringsAsFactors=FALSE)

bmc_collapse_logspots <- function(y,mapping) {
 ok<-!is.na(mapping)&nzchar(mapping)&mapping!='#N/A'
 ids<-mapping[ok];n<-table(ids);a<-rowsum(y[ok,,drop=FALSE],ids,reorder=FALSE)
 a<-a/as.numeric(n[rownames(a)])
 loci<-sub('\\.[0-9]+$','',rownames(a));n<-table(loci)
 b<-rowsum(a,loci,reorder=FALSE);b/as.numeric(n[rownames(b)])
}

bmc_fit_dose <- function(y,m) {
 m$pool<-factor(paste(m$dose,m$replicate,sep='_'))
 m$cell<-factor(paste0(m$tissue,m$dose),levels=unlist(lapply(c('Root','Shoot'),function(t)paste0(t,c(0,50,75,100)))))
 stopifnot(all(table(m$pool)==2L),all(table(m$cell)>=2L),ncol(y)==nrow(m))
 design<-model.matrix(~0+cell,m);colnames(design)<-levels(m$cell)
 rho<-limma::duplicateCorrelation(y,design,block=m$pool)$consensus.correlation
 pairs<-list(`50v0`=c(50,0),`75v0`=c(75,0),`100v0`=c(100,0),`75v50`=c(75,50),`100v75`=c(100,75))
 cm<-matrix(0,8,15,dimnames=list(colnames(design),unlist(lapply(c('Root','Shoot','Difference'),function(t)paste(t,names(pairs),sep='_')))))
 for(k in names(pairs)){
  for(t in c('Root','Shoot'))cm[paste0(t,pairs[[k]]),paste(t,k,sep='_')]<-c(1,-1)
  cm[,paste0('Difference_',k)]<-cm[,paste0('Shoot_',k)]-cm[,paste0('Root_',k)]
 }
 z<-moderated_contrasts(y,design,cm,m$pool,rho)
 z$table$q_global<-p.adjust(z$table$P,'BH')
 joint<-do.call(rbind,lapply(c('Root','Shoot','Difference'),function(t){
  a<-limma::topTable(z$contrast_fit,coef=paste(t,names(pairs)[1:3],sep='_'),number=Inf,sort.by='none')
  data.frame(locus=rownames(a),test=t,F=a$F,P=a$P.Value)
 }))
 joint$q_global<-p.adjust(joint$P,'BH')
 means<-vapply(levels(m$cell),function(cell)rowMeans(y[,m$cell==cell,drop=FALSE]),numeric(nrow(y)))
 list(table=z$table,joint=joint,means=means,pair_correlation=rho,metadata=m)
}

bmc_dose_omissions <- function(inputs,raw_source) {
 inputs<-inputs[c('microarray_root','microarray_shoot')]
 primary_y<-do.call(cbind,lapply(inputs,function(x)integration_loci(x)$y))
 m<-do.call(rbind,lapply(inputs,function(x)bmc_dose_metadata(x$samples,x$spec$tissue)))
 primary<-bmc_fit_dose(primary_y,m);base<-primary$table
 raw<-raw_source$objects$dose;regular<-raw$control==0
 target<-sub('([0-9]+)N ','\\1 NaCl ',raw$samples)
 target<-sub(' ([123])$',' B\\1',target)
 j<-match(m$sample,target);stopifnot(!anyNA(j),!anyDuplicated(target))
 raw_log<-log2(raw$values[regular,j,drop=FALSE]);mapping<-raw$data[[2]][regular]
 normalized<-lapply(c('Root','Shoot'),function(t){
  k<-m$tissue==t;bmc_collapse_logspots(limma::normalizeBetweenArrays(raw_log[,k,drop=FALSE],method='quantile'),mapping)
 })
 reconstructed<-do.call(cbind,normalized)
 reconstructed<-reconstructed[rownames(primary_y),,drop=FALSE]
 stopifnot(max(abs(primary_y-reconstructed))<1e-10)
 focal<-c(NIT2='AT3G44300',SOS1='AT2G01980',AHA2='AT4G30190',KUP4='AT4G23640',
  HKT1='AT4G10310',HAK5='AT4G13420',`CDKA;1`='AT3G48750',AT5G44720='AT5G44720',`LRK10L-1.2`='AT1G18390')
 candidates<-joint_candidates<-discovery<-patterns<-metrics<-list()
 collect<-function(z,mode,pool){
  a<-z$table[z$table$feature%in%focal,];a$Gene<-names(focal)[match(a$feature,focal)]
  a$normalization<-mode;a$omitted_pool<-pool;candidates[[length(candidates)+1L]]<<-a
  a<-z$joint[z$joint$locus%in%focal,];a$Gene<-names(focal)[match(a$locus,focal)]
  a$normalization<-mode;a$omitted_pool<-pool;joint_candidates[[length(joint_candidates)+1L]]<<-a
  for(t in c('Root','Shoot','Difference')){
   a<-z$joint[z$joint$test==t,]
   discovery[[length(discovery)+1L]]<<-data.frame(normalization=mode,omitted_pool=pool,test=t,
    loci=nrow(a),global_BH_lt_05=sum(a$q_global<.05),pair_correlation=z$pair_correlation,samples=nrow(z$metadata))
   if(t!='Difference'){
    mu<-z$means[,paste0(t,c(0,50,75,100)),drop=FALSE]
    increments<-mu[,2:4]-mu[,1:3]
    type<-ifelse(apply(increments,1,min)>=0,'increasing means',ifelse(apply(increments,1,max)<=0,'decreasing means','non-monotonic means'))
    q<-a$q_global[match(rownames(mu),a$locus)]
    counts<-table(factor(type[q<.05],levels=c('increasing means','decreasing means','non-monotonic means')))
    patterns[[length(patterns)+1L]]<<-data.frame(normalization=mode,omitted_pool=pool,tissue=t,mean_pattern=names(counts),loci=as.integer(counts))
   }
  }
 }
 collect(primary,'primary','none')
 pools<-unique(paste(m$dose,m$replicate,sep='_'))
 gene_summaries<-list()
 for(mode in c('frozen','refitted')){
  accumulator<-data.frame(feature=base$feature,contrast=base$contrast,normalization=mode,
   primary_log2FC=base$log2FC,primary_q_global=base$q_global,
   min_log2FC=base$log2FC,max_log2FC=base$log2FC,max_abs_effect_change=0,
   min_q_global=base$q_global,max_q_global=base$q_global,sign_changes=0L,largest_change_pool='none')
  for(pool in pools){
   keep<-paste(m$dose,m$replicate,sep='_')!=pool
   if(mode=='frozen')y<-primary_y[,keep,drop=FALSE] else {
    mats<-lapply(c('Root','Shoot'),function(t){
     k<-keep&m$tissue==t
     bmc_collapse_logspots(limma::normalizeBetweenArrays(raw_log[,k,drop=FALSE],method='quantile'),mapping)
    })
    y<-do.call(cbind,mats);y<-y[rownames(primary_y),,drop=FALSE]
   }
   z<-bmc_fit_dose(y,m[keep,]);a<-z$table
   stopifnot(identical(a$feature,base$feature),identical(a$contrast,base$contrast))
   collect(z,mode,pool)
   delta<-a$log2FC-base$log2FC;change<-abs(delta)>accumulator$max_abs_effect_change
   accumulator$min_log2FC<-pmin(accumulator$min_log2FC,a$log2FC)
   accumulator$max_log2FC<-pmax(accumulator$max_log2FC,a$log2FC)
   accumulator$max_abs_effect_change<-pmax(accumulator$max_abs_effect_change,abs(delta))
   accumulator$largest_change_pool[change]<-pool
   accumulator$min_q_global<-pmin(accumulator$min_q_global,a$q_global)
   accumulator$max_q_global<-pmax(accumulator$max_q_global,a$q_global)
   accumulator$sign_changes<-accumulator$sign_changes+as.integer(sign(a$log2FC)!=sign(base$log2FC))
   for(co in unique(base$contrast)){
    k<-base$contrast==co;informative<-k&abs(base$log2FC)>=.5
    metrics[[length(metrics)+1L]]<-data.frame(normalization=mode,omitted_pool=pool,contrast=co,
     flagged_pair=pool=='75_2',loci=sum(k),Pearson=cor(base$log2FC[k],a$log2FC[k]),
     Spearman=cor(base$log2FC[k],a$log2FC[k],method='spearman'),
     median_absolute_effect_change=median(abs(delta[k])),p95_absolute_effect_change=unname(quantile(abs(delta[k]),.95)),
     maximum_absolute_effect_change=max(abs(delta[k])),primary_abs_log2FC_ge_05=sum(informative),
     sign_agreement_among_effects_ge_05=if(any(informative))mean(sign(base$log2FC[informative])==sign(a$log2FC[informative])) else NA_real_,
     global_BH_lt_05=sum(a$q_global[k]<.05))
   }
   message('Completed pair omission: ',mode,' ',pool)
  }
  gene_summaries[[mode]]<-accumulator
 }
 list(dose_omission_metrics=do.call(rbind,metrics),dose_omission_candidates=do.call(rbind,candidates),
  dose_omission_joint_candidates=do.call(rbind,joint_candidates),dose_omission_discoveries=do.call(rbind,discovery),
  dose_omission_patterns=do.call(rbind,patterns),dose_omission_gene_ranges=do.call(rbind,gene_summaries),
  reconstruction_max_abs_error=max(abs(primary_y-reconstructed)),flagged_pair='75_2')
}

bmc_omission_figure <- function(z,path) {
 d<-z$dose_omission_candidates
 publication_figure(path,9.5,7.0,function(){
  par(mfrow=c(2,2),mar=c(5.3,4,2.8,.7),oma=c(.3,0,1.8,0))
  genes<-c('NIT2','SOS1','AHA2','KUP4')
  pools<-unique(d$omitted_pool[d$normalization=='frozen'])
  for(i in seq_along(genes)){
   a<-d[d$Gene==genes[i]&d$contrast=='Root_75v0',]
   base<-a$log2FC[a$normalization=='primary']
   publication_axes(c(.5,12.5),publication_range(c(a$lower95,a$upper95)),
    ylab='Root 75 versus 0 mM (log2)',at=1:12,labels=rep('',12))
   axis(1,at=1:12,labels=gsub('_',' B',pools),las=2,cex.axis=.7,tick=FALSE)
   abline(h=base,lty=2,col='#253342')
   for(mode in c('frozen','refitted')){
    b<-a[a$normalization==mode,];b<-b[match(pools,b$omitted_pool),]
    xx<-1:12+if(mode=='frozen')-.12 else .12
    color<-if(mode=='frozen')'#0072B2' else '#D55E00'
    segments(xx,b$lower95,xx,b$upper95,col=adjustcolor(color,.6),lwd=.9)
    points(xx,b$log2FC,pch=if(mode=='frozen')21 else 24,bg=color,col='white',cex=1.05)
   }
   abline(v=match('75_2',pools),col='#8B4D9C',lty=3)
   publication_title(LETTERS[i],genes[i])
  }
  mtext('Frozen normalization',3,outer=TRUE,line=.3,adj=.3,col='#0072B2',font=2,cex=.9)
  mtext('Refitted normalization',3,outer=TRUE,line=.3,adj=.75,col='#D55E00',font=2,cex=.9)
  publication_caption('Each x-axis label identifies one omitted harvested root-shoot pair; 0 B1 denotes the control pair labelled B1, not cross-treatment batch matching. All twelve omissions are shown. Points and bars are effect estimates and pointwise moderated 95% intervals from 22 remaining tissue profiles; overlapping analyses are not independent replicates. The horizontal dashed line is the unchanged full-data estimate. The vertical purple line identifies the 75 mM B2 pair containing the flagged root array. No primary sample was excluded; both normalization strategies use the same fixed gene universe and testing families.')
 })
}
