# Reconstruct between-array normalization from supplied pre-normalization signals.
normalization_settings <- function() list(primary='quantile',source='Agilent gProcessedSignal',
  groups=c('Array1_Root','Array1_Shoot','Array2_Combined'),controls='ControlType == 0 for fitting',
  q75_probability=.75,q75_type=7L,q75_target=1000,cyclicloess_iterations=3L,
  collapse='normalize spots first; preserve existing exact-transcript/locus summaries',
  background='retain vendor processed signals; no second background correction',
  filtering='no new sample or expressed-gene exclusion')

normalization_paths <- function(root) {
  p<-c(dose='Microarray_for_Stephano/Hydroponics array/Hydro all gProcessedSignals.xlsx',
    genotype='Microarray_for_Stephano/Nit2 OX array/Nit2 OX All gProcessedSignals.xlsx')
  directories<-file.path(root,c('Microarray_for_Stephano/Hydroponics array/Text files',
    'Microarray_for_Stephano/Nit2 OX array/Text files'))
  files<-unlist(lapply(directories,function(d)list.files(d,pattern='\\.txt$',full.names=TRUE)))
  files<-files[!grepl('background24',basename(files),fixed=TRUE)]
  names(files)<-substring(files,nchar(root)+2L)
  c(setNames(file.path(root,p),p),files)
}

normalization_feature_file <- function(path) {
  first<-readLines(path,n=20,warn=FALSE);i<-grep('^FEATURES\t',first)
  if(length(i)!=1L)return(NULL)
  d<-read.delim(path,skip=i-1L,check.names=FALSE,stringsAsFactors=FALSE)
  assert_columns(d,c('FeatureNum','ControlType','ProbeName','gProcessedSignal'),path)
  d
}

normalization_sources <- function(root) {
  files<-normalization_paths(root);fe_files<-files[grepl('/Text files/',names(files),fixed=TRUE)]
  reference_path<-fe_files[grepl('BG PN 1 -',basename(fe_files),fixed=TRUE)][1]
  reference<-normalization_feature_file(reference_path)
  canonical<-function(x)gsub(',','.',trimws(x),fixed=TRUE)
  objects<-list()
  for(dataset in c('dose','genotype')){
    path<-files[if(dataset=='dose')1 else 2]
    d<-as.data.frame(readxl::read_excel(path,sheet='gProcessedSignals',col_types='text',.name_repair='minimal'))
    d<-d[-1,nzchar(names(d)),drop=FALSE]
    stopifnot(nrow(d)==nrow(reference),identical(canonical(d[[1]]),canonical(reference$ProbeName)))
    samples<-names(d)[-c(1,2)];values<-vapply(d[-c(1,2)],as_numeric_checked,numeric(nrow(d)),context=path)
    stopifnot(all(is.finite(values)),all(values>0))
    objects[[dataset]]<-list(data=d,values=values,samples=samples,probe=canonical(reference$ProbeName),
      control=reference$ControlType,feature=reference$FeatureNum,reference_path=reference_path)
  }
  qc<-list()
  for(path in fe_files){
    name<-sub(' - .*','',basename(path));dataset<-if(grepl('^B[123] ',name))'dose' else 'genotype'
    sample<-if(dataset=='dose'){
      a<-strsplit(name,' ',fixed=TRUE)[[1]];paste(a[2],a[3],sub('B','',a[1],fixed=TRUE))
    }else name
    col<-match(sample,objects[[dataset]]$samples);stopifnot(!is.na(col))
    full<-normalization_feature_file(path)
    row<-data.frame(dataset=dataset,sample=sample,full_FE=!is.null(full),rows=NA_integer_,
      max_processed_relative_difference=NA_real_,saturated=NA_integer_,feature_nonuniform=NA_integer_,
      well_above_background=NA_integer_,stringsAsFactors=FALSE)
    if(!is.null(full)){
      stopifnot(identical(canonical(full$ProbeName),objects[[dataset]]$probe),
        identical(full$FeatureNum,reference$FeatureNum),identical(full$ControlType,reference$ControlType))
      delta<-max(abs(full$gProcessedSignal-objects[[dataset]]$values[,col])/pmax(full$gProcessedSignal,1))
      if(delta>1e-10)stop('Workbook/FE values do not agree: ',path)
      row$rows<-nrow(full);row$max_processed_relative_difference<-delta
      for(pair in list(c('saturated','gIsSaturated'),c('feature_nonuniform','gIsFeatNonUnifOL'),c('well_above_background','gIsWellAboveBG')))
        if(pair[2]%in%names(full))row[[pair[1]]]<-sum(full[[pair[2]]]>0,na.rm=TRUE)
    }else{
      reduced<-read.delim(path,header=FALSE,skip=1,dec=',',fill=TRUE,stringsAsFactors=FALSE)
      numeric<-suppressWarnings(as.numeric(sub(',','.',as.character(reduced[[1]]),fixed=TRUE)))
      stopifnot(length(numeric)==nrow(reference),all(is.finite(numeric)))
      delta<-max(abs(numeric-objects[[dataset]]$values[,col])/pmax(numeric,1))
      if(delta>1e-10)stop('Reduced export/workbook signals do not agree: ',path)
      row$rows<-length(numeric);row$max_processed_relative_difference<-delta
    }
    qc[[length(qc)+1L]]<-row
  }
  list(objects=objects,feature_qc=do.call(rbind,qc),input_paths=files)
}

normalization_duplicate_mad <- function(y,probe) {
  groups<-split(seq_len(nrow(y)),probe);groups<-groups[lengths(groups)>1L]
  two<-groups[lengths(groups)==2L]
  result<-numeric()
  if(length(two)){
    i<-vapply(two,`[`,0L,1L);j<-vapply(two,`[`,0L,2L)
    result<-apply(abs(y[i,,drop=FALSE]-y[j,,drop=FALSE])*.7413,1,median)
  }
  more<-groups[lengths(groups)>2L]
  if(length(more))result<-c(result,vapply(more,function(k)median(apply(y[k,,drop=FALSE],2,mad)),0.))
  median(result)
}

normalization_variants <- function(source) {
  rows<-list();matrices<-list()
  for(dataset in names(source$objects)){
    z<-source$objects[[dataset]];regular<-z$control==0
    for(tissue in if(dataset=='dose')c('Root','Shoot') else 'Combined'){
      columns<-if(dataset=='dose')grep(tolower(tissue),z$samples) else seq_along(z$samples)
      raw<-log2(z$values[,columns,drop=FALSE]);groups<-sub(' [123]$','',colnames(raw))
      factors<-log2(apply(z$values[regular,columns,drop=FALSE],2,quantile,.75,type=normalization_settings()$q75_type))
      methods<-list(processed=raw,q75=sweep(raw,2,factors,'-')+log2(normalization_settings()$q75_target))
      for(method in c('quantile','cyclicloess')){
        y<-raw;y[regular,]<-limma::normalizeBetweenArrays(raw[regular,,drop=FALSE],method=method)
        for(i in seq_along(columns))y[!regular,i]<-approx(raw[regular,i],y[regular,i],xout=raw[!regular,i],rule=2,ties=mean)$y
        methods[[method]]<-y
      }
      low<-apply(raw[regular,,drop=FALSE],2,min);high<-apply(raw[regular,,drop=FALSE],2,max)
      inrange<-rowSums(sweep(raw,2,low,'>=')&sweep(raw,2,high,'<='))==ncol(raw)
      for(method in names(methods)){
        y<-methods[[method]];signal<-y[regular,,drop=FALSE];residual<-signal
        for(g in unique(groups)){k<-which(groups==g);residual[,k]<-sweep(signal[,k,drop=FALSE],1,rowMeans(signal[,k,drop=FALSE]),'-')}
        rms<-sqrt(rowMeans(residual^2));reference<-rowMeans(signal)
        control_rms<-function(k){
          a<-y[k,,drop=FALSE]
          for(g in unique(groups)){j<-which(groups==g);a[,j]<-sweep(a[,j,drop=FALSE],1,rowMeans(a[,j,drop=FALSE]),'-')}
          median(sqrt(rowMeans(a^2)))
        }
        positive<-z$control==1&!z$probe%in%c('GE_BrightCorner','DarkCorner')&inrange
        negative<-z$control==-1&inrange
        bias<-vapply(seq_len(ncol(signal)),function(j){
          peers<-which(groups==groups[j]&seq_len(ncol(signal))!=j)
          ref<-rowMeans(signal[,peers,drop=FALSE]);a<-(signal[,j]+ref)/2;m<-signal[,j]-ref
          breaks<-unique(quantile(a,seq(0,1,length.out=11)))
          median_differences<-tapply(m,cut(a,breaks,include.lowest=TRUE),median)
          sqrt(mean(median_differences^2))
        },0.)
        rows[[length(rows)+1L]]<-data.frame(dataset=dataset,tissue=tissue,method=method,regular_probes=sum(regular),
          median_within_group_RMS=median(rms),upper_signal_RMS=median(rms[reference>quantile(reference,.75)]),
          median_duplicate_spot_MAD=normalization_duplicate_mad(signal,z$probe[regular]),
          positive_control_within_group_RMS=control_rms(positive),negative_control_within_group_RMS=control_rms(negative),
          positive_inrange_features=sum(positive),negative_inrange_features=sum(negative),median_intensity_bias_RMS=median(bias))
      }
      group<-paste(if(dataset=='dose')'Array1' else 'Array2',tissue,sep='_')
      matrices[[group]]<-list(methods=methods,regular=regular,columns=columns,dataset=dataset,tissue=tissue,
        groups=groups,source=z)
    }
  }
  list(metrics=do.call(rbind,rows),matrices=matrices)
}

normalization_comparison <- function(result,source_root) {
  rows<-list()
  for(group in names(result$matrices)){
    z<-result$matrices[[group]];id<-if(z$dataset=='genotype')'microarray_genotype' else paste0('microarray_',tolower(z$tissue))
    old<-read_input(id,source_root,input_specs)
    target<-if(z$dataset=='genotype')sub('100N','100 NaCl',colnames(z$methods$processed),fixed=TRUE) else {
      a<-sub('([0-9]+)N ','\\1 NaCl ',colnames(z$methods$processed));sub(' ([123])$',' B\\1',a)
    }
    j<-match(old$samples,target);stopifnot(!anyNA(j))
    mapped<-!is.na(z$source$data[[2]])&nzchar(z$source$data[[2]])&z$source$data[[2]]!='#N/A'
    stopifnot(identical(z$source$data[[2]][mapped],old$ids[old$mapped]))
    supplied<-log2(old$values[old$mapped,,drop=FALSE]);raw<-z$methods$processed[mapped,j,drop=FALSE]
    ratio<-raw-supplied
    regenerated<-z$methods$quantile[mapped,j,drop=FALSE]
    rows[[group]]<-data.frame(group=group,sample=old$samples,
      supplied_log2_scale=apply(ratio,2,median),supplied_scale_residual_MAD=apply(ratio,2,mad),
      supplied_scale_residual_range=apply(ratio,2,function(x)diff(range(x))),
      quantile_vs_supplied_median_abs_log2=apply(abs(regenerated-supplied),2,median),
      quantile_vs_supplied_max_abs_log2=apply(abs(regenerated-supplied),2,max))
  }
  do.call(rbind,rows)
}

normalization_diagnostics_figure <- function(result,path) {
  publication_figure(path,12,8.7,function(){
    par(mfrow=c(3,3),mar=c(4.2,4.1,2.8,.8),oma=c(1.2,0,1,0))
    for(group in names(result$matrices))for(method in c('processed','q75','quantile')){
      z<-result$matrices[[group]];y<-z$methods[[method]][z$regular,,drop=FALSE]
      boxplot(y,outline=FALSE,axes=FALSE,col=grDevices::adjustcolor('#0072B2',.22),border='#0072B2',
        xlab='',ylab='Log2 processed signal')
      axis(2,cex.axis=.75);axis(1,at=seq_len(ncol(y)),labels=seq_len(ncol(y)),cex.axis=.7)
      box(bty='l');title(paste(group,method),cex.main=.95)
      mtext('Sample order in normalization_samples.csv',side=1,line=2.5,cex=.66)
    }
    mtext('Equal distributions are imposed by quantile normalization; method choice also uses independent control and replicate diagnostics.',side=1,outer=TRUE,cex=.78)
  })
}

normalization_prepare_root <- function(source_root,destination,result) {
  source_root<-normalizePath(source_root,mustWork=TRUE);destination<-normalizePath(destination,mustWork=FALSE)
  inside<-function(p,parent)identical(p,parent)||startsWith(p,paste0(parent,'/'))
  if(inside(destination,source_root)||inside(source_root,destination)||inside(destination,project_root))stop('Keep normalized inputs separate from code and original data')
  if(dir.exists(destination)&&length(list.files(destination,all.files=TRUE,no..=TRUE))&&
     !file.exists(file.path(destination,'NORMALIZATION_PROVENANCE.R')))stop('Destination is not an owned normalized-input root')
  dir.create(destination,recursive=TRUE,showWarnings=FALSE)
  original_files<-list.files(source_root,recursive=TRUE,all.files=TRUE,full.names=TRUE,include.dirs=FALSE)
  for(p in original_files){
    rel<-substring(p,nchar(source_root)+2L);to<-file.path(destination,rel)
    dir.create(dirname(to),recursive=TRUE,showWarnings=FALSE)
    if(!file.exists(to))if(!file.symlink(p,to))stop('Cannot link unchanged source: ',rel)
  }
  normalized<-character()
  for(group in names(result$matrices)){
    z<-result$matrices[[group]];id<-if(z$dataset=='genotype')'microarray_genotype' else paste0('microarray_',tolower(z$tissue))
    old<-read_input(id,source_root,input_specs)
    target<-if(z$dataset=='genotype')sub('100N','100 NaCl',colnames(z$methods$processed),fixed=TRUE) else {
      a<-sub('([0-9]+)N ','\\1 NaCl ',colnames(z$methods$processed));sub(' ([123])$',' B\\1',a)
    }
    index<-match(old$samples,target);stopifnot(!anyNA(index))
    d<-data.frame(mapping=z$source$data[[2]][z$regular],check.names=FALSE)
    names(d)<-'CORRECTED MAPPING';d$`CORRECTED MAPPING`[is.na(d$`CORRECTED MAPPING`)|!nzchar(d$`CORRECTED MAPPING`)]<-'#N/A'
    y<-2^z$methods$quantile[z$regular,index,drop=FALSE];colnames(y)<-old$samples
    d<-cbind(d,as.data.frame(y,check.names=FALSE));path<-file.path(destination,input_specs[[id]]$path)
    if(file.exists(path))unlink(path)
    write_table(d,path);normalized[id]<-path
  }
  normalized
}
