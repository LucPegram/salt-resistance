extension_columns <- function(objects) {
  do.call(rbind,lapply(objects,function(o)data.frame(dataset=o$dataset,contrast=o$selected_contrasts,
    key=paste(o$dataset,o$selected_contrasts,sep="::"),stringsAsFactors=FALSE)))
}

extension_column_label <- function(dataset,contrast) {
  if(dataset=="Array1")return(paste("A1",substr(contrast,1,1),sub(".*_([0-9]+)v0","\\1",contrast)))
  ds<-switch(dataset,Array2="A2",Protein_Root="PR",Protein_Shoot="PS")
  co<-c(baseline_genotype="basal",salt_WT="WT salt",salt_OE="OE salt",saline_genotype="salt geno",interaction="G x S")
  paste(ds,co[contrast])
}

extension_heatmap <- function(path,mat,title,subtitle,columns,stars=NULL,height=7.4,
  labels=rownames(mat),legend="Signed -log10(global BH q); color capped at 6",limit=6) {
  palette<-colorRampPalette(c("#0072B2","#F8FAFC","#C44E52"))(201)
  publication_figure(path,10.5,height,function(){
    par(mar=c(6.7,17,3.8,1.2),oma=c(1.4,0,0,0))
    nr<-nrow(mat);nc<-ncol(mat)
    plot(NA,xlim=c(.5,nc+.5),ylim=c(.5,nr+.5),axes=FALSE,xlab="",ylab="",xaxs="i",yaxs="i")
    for(i in seq_len(nr))for(j in seq_len(nc)){
      v<-mat[i,j];color<-if(is.finite(v))palette[round((max(-limit,min(limit,v))+limit)/(2*limit)*200)+1] else "#E3E8ED"
      rect(j-.5,nr-i+.5,j+.5,nr-i+1.5,col=color,border="white",lwd=.35)
      if(!is.null(stars)&&isTRUE(stars[i,j]))text(j,nr-i+1,"*",cex=.92,col=if(abs(v)>3)"white" else "#253342")
    }
    axis(2,at=nr:1,labels=labels,las=1,tick=FALSE,cex.axis=.72,line=-.5)
    axis(1,at=seq_len(nc),labels=mapply(extension_column_label,columns$dataset,columns$contrast),
      las=2,tick=FALSE,cex.axis=.74,line=-.4)
    boundaries<-which(head(columns$dataset,-1)!=tail(columns$dataset,-1))+.5
    abline(v=boundaries,col="#596775",lwd=1.4)
    title(main=title,sub="",adj=0,cex.main=1.1,line=2.3)
    mtext(subtitle,3,line=.9,adj=0,cex=.78,col=publication_colors["muted"])
    mtext(legend,1,outer=TRUE,line=.2,cex=.75)
    mtext("A1: Col-0 dose series | A2: No-0/NIT2 array | PR/PS: root/shoot protein; grey: untested",1,outer=TRUE,line=1.1,cex=.7)
  })
}

extension_matrix <- function(d,rowfield,rows,columns,value) {
  m<-matrix(NA_real_,length(rows),nrow(columns),dimnames=list(rows,columns$key))
  ri<-match(d[[rowfield]],rows);ci<-match(paste(d$dataset,d$contrast,sep="::"),columns$key)
  ok<-!is.na(ri)&!is.na(ci)
  if(anyDuplicated(paste(ri[ok],ci[ok])))stop("Duplicate heatmap cells")
  m[cbind(ri[ok],ci[ok])]<-value[ok];m
}

extension_benchmark_figure <- function(z,path,correlation=FALSE) {
  d<-if(correlation)z$benchmark$correlations else z$benchmark$table
  cols<-extension_columns(z$objects);refs<-sort(unique(d$reference))
  if(correlation){v<-d$Spearman;limit<-.5;legend<-"Spearman correlation of gene-level effects; color capped at +/-0.5; no gene-wise correlation P-values";stars<-NULL}
  else{v<-sign(d$score_difference)*pmin(6,-log10(pmax(d$q_global,1e-300)));limit<-6
    legend<-"Signed -log10(global BH q); color capped at 6; * global q < 0.05; positive = stronger external salt signature"
    stars<-extension_matrix(d,"reference",refs,cols,d$q_global<.05)}
  mat<-extension_matrix(d,"reference",refs,cols,v)
  labels<-gsub("AtGen_","AtGenExpress ",refs,fixed=TRUE)
  labels<-gsub("Hydro15min_Root_dose","Hydroponic root 15 min / ",labels,fixed=TRUE)
  labels<-gsub("_hour"," / ",labels,fixed=TRUE);labels<-gsub("_"," ",labels,fixed=TRUE)
  labels<-ifelse(startsWith(refs,"AtGen"),paste0(labels," h"),paste0(labels," mM"))
  extension_heatmap(path,mat,if(correlation)"External salt-effect concordance" else "Independent salt-response signatures",
    "References: AtGenExpress 150 mM, 0.5-24 h; GSE223207 root, 25-100 mM, 15 min",cols,stars,
    labels=labels,height=7.1,limit=limit,legend=legend)
}

extension_tf_figure <- function(z,path) {
  d<-z$tf$table[z$tf$table$mode=="primary",];cols<-extension_columns(z$objects)
  rank<-aggregate(q_global~TF,d,min);rank<-rank[order(rank$q_global,rank$TF),]
  selected<-head(rank$TF,20)
  labels<-unique(d[c("TF","TF_name")]);labels<-labels$TF_name[match(selected,labels$TF)]
  mat<-extension_matrix(d,"TF",selected,cols,ifelse(d$Direction=="Up",1,-1)*pmin(6,-log10(pmax(d$q_global,1e-300))))
  stars<-extension_matrix(d,"TF",selected,cols,d$q_global<.05)
  extension_heatmap(path,mat,"Transcription-factor target enrichment",
    "20 TFs with the smallest global q; DAP-seq targets; competitive rank-based CAMERA, correlation 0.01",cols,stars,
    labels=paste(labels,selected,sep=" / "),height=8.2,
    legend="Signed -log10(global BH q); * q < 0.05; positive = targets shifted upward; enrichment does not establish TF activity")
}

extension_marker_figure <- function(z,path) {
  d<-z$marker_result$table[z$marker_result$table$mode=="primary",];cols<-extension_columns(z$objects)
  preferred<-unique(d$set[d$category%in%c("root_development","proliferation")])
  rank<-aggregate(q_global~set,d,min);rank<-rank[order(rank$q_global,rank$set),]
  selected<-unique(c(sort(preferred),head(setdiff(rank$set,preferred),20)))
  mat<-extension_matrix(d,"set",selected,cols,sign(d$score_difference)*pmin(6,-log10(pmax(d$q_global,1e-300))))
  stars<-extension_matrix(d,"set",selected,cols,d$q_global<.05)
  labels<-sub("Root Cell Type Markers / ","Root / ",selected,fixed=TRUE)
  labels<-sub("Root developmental stage / ","Root stage / ",labels,fixed=TRUE)
  labels<-sub("Other Cell Type Markers / ","Cycle / ",labels,fixed=TRUE)
  extension_heatmap(path,mat,"Developmental and tissue-marker signatures",
    "All eligible stage/cycle sets plus 20 strongest other marker sets; independent published atlas",cols,stars,
    labels=labels,height=max(8.5,length(selected)*.24+2.2),
    legend="Signed -log10(global BH q); * q < 0.05; scores reflect marker expression, not cell proportions or inferred cell loss")
}

extension_dose_figure <- function(z,path) {
  p<-z$dose$patterns
  publication_figure(path,9.5,7.3,function(){
    par(mfrow=c(2,2),mar=c(4,6.5,3.2,1))
    for(t in c("Root","Shoot")){
      d<-p[p$tissue==t&p$joint_q_global<.05,]
      types<-c("increasing means","decreasing means","non-monotonic means")
      counts<-table(factor(d$mean_pattern,levels=types))
      par(mar=c(4,10,3.2,1))
      bp<-barplot(counts,horiz=TRUE,las=1,col=c("#C44E52","#0072B2","#8C95A4"),border=NA,
        xlim=c(0,max(counts)*1.2),xlab="Loci with dose-joint global q < 0.05",cex.names=.83)
      text(as.numeric(counts)+max(counts)*.035,bp,counts,adj=0,cex=.85)
      title(paste(t,"dose-response mean shapes"),adj=0,cex.main=1.05)
      mtext("Shapes describe observed means; no monotonicity test is implied",3,line=.25,adj=0,cex=.68)
      par(mar=c(4,6.5,3.2,1))
      d<-head(d[order(d$joint_q_global,-abs(d$endpoint_log2FC),d$locus),],20)
      mat<-z$dose$means[d$locus,paste0(t,c(0,50,75,100)),drop=FALSE]
      mat<-mat-mat[,1];limit<-max(abs(mat));pal<-colorRampPalette(c("#0072B2","#FAFBFC","#C44E52"))(201)
      image(1:4,1:nrow(mat),t(mat[nrow(mat):1,,drop=FALSE]),col=pal,zlim=c(-limit,limit),axes=FALSE,xlab="NaCl (mM)",ylab="")
      axis(1,at=1:4,labels=c(0,50,75,100),tick=FALSE)
      names<-z$annotation$Gene_name[match(d$locus,z$annotation$locus)]
      names[is.na(names)|nchar(names)>18]<-d$locus[is.na(names)|nchar(names)>18]
      axis(2,at=nrow(mat):1,labels=names,las=1,tick=FALSE,cex.axis=.66)
      title(paste(t,"20 strongest dose-joint loci"),adj=0,cex.main=.97)
      mtext(sprintf("Mean log2 signal minus control; blue < 0 < red; range +/- %.1f",limit),3,line=.25,adj=0,cex=.68)
    }
  })
}

extension_candidate_figure <- function(z,path) {
  loci<-c(NIT2="AT3G44300",SOS1="AT2G01980",AHA2="AT4G30190",`KUP4 / POT3`="AT4G23640",
    `CDKA;1`="AT3G48750",`LRK10L-1.2`="AT1G18390")
  d<-z$dose$table
  publication_figure(path,9.6,6.6,function(){
    par(mfrow=c(2,3),mar=c(4,4,3,.7),oma=c(1,0,1.5,0))
    for(i in seq_along(loci)){
      dd<-d[d$feature==loci[i]&grepl("^(Root|Shoot)_(50v0|75v50|100v75)$",d$contrast),]
      publication_axes(c(.65,3.35),publication_range(c(dd$lower95,dd$upper95)),"Dose increment (mM)","Change in log2 signal",
        at=1:3,labels=c("0 to 50","50 to 75","75 to 100"));abline(h=0,lty=3,col="#8C95A4")
      for(t in c("Root","Shoot")){
        a<-dd[match(paste(t,c("50v0","75v50","100v75"),sep="_"),dd$contrast),]
        x<-1:3+if(t=="Root")-.08 else .08;color<-publication_colors[t]
        segments(x,a$lower95,x,a$upper95,col=color,lwd=1.3)
        points(x,a$log2FC,pch=if(t=="Root")21 else 24,bg=color,col="white",cex=1.25)
      }
      publication_title(LETTERS[i],names(loci)[i],paste0(loci[i]," | pointwise 95% intervals"))
    }
    publication_key(footnote="Intervals describe each increment; direct tissue-difference tests and global q-values are in the companion tables.")
  })
}
