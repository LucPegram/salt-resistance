# Publication styling only; model fits and prespecified feature selection are unchanged.
publication_colors <- c(Root="#0072B2", Shoot="#D55E00", `No-0`="#0072B2",
  NIT2="#8B4D9C", ink="#253342", muted="#5D6975", grid="#E6EBEF")

# Explanatory prose belongs to the figure caption, never inside the artwork.
.publication_caption <- new.env(parent=emptyenv())
publication_caption <- function(text) {
  if(length(text)) .publication_caption$lines <- unique(c(.publication_caption$lines,
    as.character(text[nzchar(text)])))
  invisible(NULL)
}

publication_figure <- function(path,width,height,draw) {
  stem<-sub("\\.pdf$","",path)
  .publication_caption$lines <- character()
  for(format in c("pdf","svg","tiff","png")) {
    target<-paste0(stem,".",format)
    switch(format,
      pdf=grDevices::cairo_pdf(target,width=width,height=height,family="sans",pointsize=10,bg="white"),
      svg=grDevices::svg(target,width=width,height=height,family="sans",pointsize=10,bg="white"),
      tiff=grDevices::tiff(target,width=width,height=height,units="in",res=600,
        compression="lzw",type="cairo",pointsize=10,bg="white"),
      png=grDevices::png(target,width=width,height=height,units="in",res=200,
        type="cairo",pointsize=10,bg="white"))
    tryCatch({
      par(family="sans",las=1,fg=publication_colors["ink"],col.axis=publication_colors["muted"],
        col.lab=publication_colors["ink"],col.main=publication_colors["ink"],
        cex.axis=.85,cex.lab=.9,mgp=c(2.25,.6,0),tck=-.018,bty="l",lend="round")
      draw()
    },finally=grDevices::dev.off())
  }
  writeLines(.publication_caption$lines,paste0(stem,"_caption.txt"),useBytes=TRUE)
}

publication_range <- function(y,padding=.12) {
  r<-range(y,finite=TRUE);r+c(-1,1)*max(diff(r),.2)*padding
}

publication_axes <- function(xlim,ylim,xlab="",ylab="",at=NULL,labels=at) {
  # mfrow otherwise silently shrinks axis text in dense multi-panel figures.
  par(cex=1)
  plot(NA,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,axes=FALSE,xaxs="i",yaxs="i")
  yt<-pretty(ylim,n=4);yt<-yt[yt>=ylim[1]&yt<=ylim[2]]
  abline(h=yt,col=publication_colors["grid"],lwd=.7)
  axis(2,at=yt,lwd=0,lwd.ticks=.6)
  if(is.null(at))axis(1,lwd=0,lwd.ticks=.6) else if(any(grepl("\n",labels,fixed=TRUE))){
    axis(1,at=at,labels=FALSE,lwd=0,lwd.ticks=.6)
    parts<-strsplit(labels,"\n",fixed=TRUE)
    for(i in seq_along(at))for(j in seq_along(parts[[i]]))
      mtext(parts[[i]][j],side=1,at=at[i],line=.6+(j-1)*.9,cex=.75,col=publication_colors["muted"])
  }else axis(1,at=at,labels=labels,lwd=0,lwd.ticks=.6)
  box(bty="l",col="#A9B4BE",lwd=.65)
}

publication_title <- function(letter,title,subtitle="") {
  mtext(paste(letter,title,sep="   "),side=3,line=1.1,adj=0,font=2,cex=.98)
  if(nzchar(subtitle))publication_caption(paste0("Panel ",letter," (",title,"): ",subtitle,"."))
}

publication_key <- function(tissue=TRUE,footnote="") {
  labels<-if(tissue)c("Root","Shoot") else c("No-0","NIT2")
  text<-if(tissue)labels else c("No-0","NIT2 overexpressor")
  for(i in 1:2)mtext(text[i],side=3,outer=TRUE,line=.45,adj=c(.32,.68)[i],
    col=publication_colors[labels[i]],cex=.9,font=2)
  if(nzchar(footnote))publication_caption(gsub("\n"," ",footnote,fixed=TRUE))
}

# Small deterministic label placement avoids points and earlier label rectangles.
publication_pool_labels <- function(x,y,labels) {
  usr<-par("usr");sx<-diff(usr[1:2]);sy<-diff(usr[3:4]);prior<-list()
  directions<-rbind(c(0,1),c(0,-1),c(1,0),c(-1,0),c(1,1),c(-1,1),c(1,-1),c(-1,-1))
  for(i in seq_along(x)) {
    w<-strwidth(labels[i],cex=.7)/2+sx*.008;h<-strheight(labels[i],cex=.7)/2+sy*.008
    choices<-cbind(x[i]+directions[,1]*(w+sx*.018),y[i]+directions[,2]*(h+sy*.025))
    score<-apply(choices,1,function(p){
      b<-c(p[1]-w,p[1]+w,p[2]-h,p[2]+h)
      hits<-sum(x>=b[1]&x<=b[2]&y>=b[3]&y<=b[4])
      overlap<-sum(vapply(prior,function(a)b[1]<a[2]&&b[2]>a[1]&&b[3]<a[4]&&b[4]>a[3],logical(1)))
      outside<-b[1]<usr[1]||b[2]>usr[2]||b[3]<usr[3]||b[4]>usr[4]
      100*(hits+overlap)+1000*outside
    })
    p<-choices[which.min(score),];text(p[1],p[2],labels[i],cex=.7)
    prior[[i]]<-c(p[1]-w,p[1]+w,p[2]-h,p[2]+h)
  }
}

# All combinations share one PCA within each assay/tissue; do not refit by cell.
publication_pca_data <- function(inputs) {
  rows<-lapply(inputs,function(x){
    d<-revision_qc(x)
    factorial<-x$spec$assay=="proteomics"||identical(x$spec$tissue,"Combined")
    if(factorial){
      m<-factorial_metadata(x);d$genotype<-m$genotype;d$salt_mM<-m$salt_mM
      d$dataset<-if(x$spec$assay=="proteomics")"Proteomics" else "Array 2"
    }else{
      d$genotype<-"Col-0"
      d$salt_mM<-as.integer(ifelse(startsWith(x$samples,"PN"),"0",sub(" .*","",x$samples)))
      d$dataset<-"Array 1"
    }
    d$replicate<-sub(".*[ _]","",x$samples)
    d$group<-paste(d$genotype,d$salt_mM,"mM")
    d
  })
  do.call(rbind,rows)
}

build_manuscript_pca_coordinates <- function(inputs,settings,destination)
  write_table(publication_pca_data(inputs),destination)

build_revision_qc_figure <- function(inputs,settings,destination) {
  d<-publication_pca_data(inputs)
  dose_colors<-c("#677487","#0072B2","#009E73","#D55E00")
  publication_figure(destination,7.2,7.4,function(){
    par(mfrow=c(3,2),mar=c(3.6,4,2.7,.8),oma=c(.5,0,.1,0))
    for(i in seq_along(inputs)){
      x<-inputs[[i]];z<-d[d$dataset==if(x$spec$assay=="proteomics")"Proteomics" else
        if(identical(x$spec$tissue,"Combined"))"Array 2" else "Array 1",]
      z<-z[z$tissue==x$spec$tissue,]
      if(z$dataset[1]=="Array 1"){
        k<-match(z$salt_mM,c(0,50,75,100));colors<-dose_colors[k];shapes<-c(21,22,23,24)[k]
      }else{
        colors<-publication_colors[z$genotype];shapes<-ifelse(z$salt_mM==0,21,24)
      }
      publication_axes(publication_range(z$PC1,.22),publication_range(z$PC2,.22),
        sprintf("PC1 (%.1f%%)",100*z$PC1_fraction[1]),sprintf("PC2 (%.1f%%)",100*z$PC2_fraction[1]))
      abline(h=0,v=0,col="#CCD5DE",lty=3,lwd=.65)
      points(z$PC1,z$PC2,pch=shapes,bg=colors,col="white",cex=1.4,lwd=.6)
      publication_pool_labels(z$PC1,z$PC2,z$replicate)
      title<-paste(z$dataset[1],if(x$spec$tissue=="Combined")"combined tissue" else tolower(x$spec$tissue))
      units<-if(x$spec$assay=="proteomics")"protein groups" else "transcript IDs"
      publication_title(LETTERS[i],title,
        paste0(format(z$PCA_features[1],big.mark=",")," ",units," | ",nrow(z)," pools"))
    }
    plot.new();par(usr=c(0,1,0,1))
    key<-function(x,y,label,shape,color){
      points(x,y,pch=shape,bg=color,col="white",cex=1.25)
      text(x+.035,y,label,adj=0,cex=.76)
    }
    text(0,.96,"Treatment and genotype keys",adj=0,font=2,cex=.94)
    text(0,.84,"Array 1 | Col-0",adj=0,font=2,cex=.84)
    for(j in 1:4)key(c(.04,.04,.55,.55)[j],c(.73,.62,.73,.62)[j],
      paste(c(0,50,75,100)[j],"mM"),c(21,22,23,24)[j],dose_colors[j])
    text(0,.48,"Array 2 and proteomics",adj=0,font=2,cex=.84)
    labels<-c("No-0 control","No-0 salt","NIT2 control","NIT2 salt")
    for(j in 1:4)key(c(.04,.04,.55,.55)[j],c(.37,.26,.37,.26)[j],labels[j],
      c(21,24,21,24)[j],publication_colors[c("No-0","No-0","NIT2","NIT2")][j])
    publication_caption("Salt: Array 2, 100 mM; proteomics, 75 mM. Labels identify pool IDs. Features are centered without variance scaling; no imputation.")

  })
}

revision_candidate_transcript <- function(y,locus) {
  hits<-grep(paste0("^",locus,"(\\.|$)"),rownames(y),value=TRUE)
  if(!length(hits))stop("Prespecified locus absent: ",locus)
  if(paste0(locus,".1")%in%hits)paste0(locus,".1") else sort(hits)[1]
}

build_revision_dose_figure <- function(inputs,settings,destination) {
  out<-revision_dose(inputs)
  loci<-c(NIT2="AT3G44300",SOS1="AT2G01980",AHA2="AT4G30190",KUP4="AT4G23640")
  publication_figure(destination,7.2,5.7,function(){
    par(mfrow=c(2,2),mar=c(3.8,4.3,2.7,.8),oma=c(.5,0,1.8,0))
    for(i in seq_along(loci)){
      id<-revision_candidate_transcript(out$y,loci[i]);v<-out$y[id,];m<-out$metadata
      publication_axes(c(-6,106),publication_range(v),"NaCl (mM)","Log2 normalized signal",at=c(0,50,75,100))
      for(t in c("Root","Shoot")){
        k<-m$tissue==t;color<-publication_colors[t];means<-tapply(v[k],m$dose[k],mean)
        lines(as.numeric(names(means)),means,col=color,lwd=1.7,lty=if(t=="Root")1 else 2)
        offset<-if(t=="Root")-1.4 else 1.4
        points(m$dose[k]+offset+(m$replicate[k]-2)*.7,v[k],bg=color,col="white",lwd=.4,
          pch=if(t=="Root")21 else 24,cex=1.05)
      }
      joint<-out$joint[out$joint$feature==id&out$joint$test=="tissue_by_dose",]
      stopifnot(nrow(joint)==1L)
      publication_title(LETTERS[i],if(names(loci)[i]=="KUP4")"KUP4 / POT3" else names(loci)[i],
        paste0(id," | tissue x dose q = ",formatC(joint$q,digits=2,format="g")))
    }
    publication_key()
  })
}

build_revision_protein_figure <- function(inputs,settings,destination) {
  panels<-lapply(inputs,revision_protein)
  publication_figure(destination,7.2,6.1,function(){
    par(mfrow=c(2,3),mar=c(4.5,3.8,2.9,.65),oma=c(1.9,0,1.8,0))
    for(i in seq_along(panels))for(j in 1:3){
      z<-panels[[i]];g<-c("NIT2","NHX7","AHA2")[j];k<-which(z$annotations$Genes==g)
      if(length(k)!=1L)stop("Ambiguous candidate protein")
      y<-z$y[k,];m<-z$metadata
      counts<-vapply(levels(m$cell),function(cell)sum(is.finite(y[m$cell==cell])),integer(1))
      publication_axes(c(.55,4.45),publication_range(y),ylab="Log2 abundance",at=1:4,
        labels=paste0(c("No-0\n0","No-0\n75","NIT2\n0","NIT2\n75"),"\nn=",counts))
      for(cell in 1:4){
        v<-y[m$cell==levels(m$cell)[cell]];genotype<-if(cell<=2)"No-0" else "NIT2"
        observed<-v[is.finite(v)]
        # With one observed pool, show only that pool rather than a degenerate box.
        if(length(observed)>=2L)boxplot(list(observed),at=cell,add=TRUE,axes=FALSE,
          boxwex=.5,outline=FALSE,range=1.5,
          col=grDevices::adjustcolor(publication_colors[genotype],alpha.f=.16),
          border=publication_colors[genotype],lwd=.9,medlwd=1.4,
          medcol=publication_colors["ink"],whisklty=1,staplewex=.5)
        points(cell+c(-.14,0,.14),v,pch=if(cell%%2==1)21 else 24,
          bg=publication_colors[genotype],col="white",lwd=.5,cex=1.15)
      }
      label<-if(g=="NHX7")"SOS1 / NHX7" else if(g=="AHA2")"AHA2 group*" else g
      publication_title(LETTERS[(i-1)*3+j],paste(label,inputs[[i]]$spec$tissue),
        paste0(z$annotations$Protein.Group[k],"   |   ",sum(is.finite(y)),"/12 observed"))
    }
    publication_key(FALSE,paste0("Circles: control; triangles: salt; boxes: median / hinges; whiskers: 1.5 x IQR; n=1: point only\n",
      "NaCl (mM) and observed n below genotype. *AHA2-labelled group; isoform specificity is unresolved"))
  })
}

build_revision_contrast_figure <- function(inputs,settings,destination) {
  d<-do.call(rbind,lapply(inputs,function(x)revision_protein(x)$table))
  d<-d[d$method=="observed_min2"&d$Genes%in%c("NIT2","NHX7","AHA2"),]
  contrasts<-c("baseline_genotype","salt_WT","saline_genotype","interaction")
  titles<-c("Basal genotype","Salt response in No-0","Genotype under salt","Genotype by salt")
  subtitles<-c("NIT2 - No-0 at 0 mM","75 mM - control","NIT2 - No-0 at 75 mM","Difference between salt responses")
  publication_figure(destination,7.2,5.8,function(){
    par(mfrow=c(2,2),mar=c(3.7,4.9,2.8,.8),oma=c(1.2,0,1.8,0))
    par(cex=1)
    for(i in seq_along(contrasts)){
      z<-d[d$contrast==contrasts[i],];z<-z[order(match(z$Genes,c("NIT2","NHX7","AHA2")),z$tissue),]
      labels<-paste(ifelse(z$Genes=="NHX7","SOS1",ifelse(z$Genes=="AHA2","AHA2*",z$Genes)),z$tissue)
      yy<-rev(seq_len(nrow(z)));xr<-publication_range(c(0,z$lower95,z$upper95),.1)
      plot(NA,xlim=c(xr[1],xr[2]+diff(xr)*.29),ylim=c(.5,nrow(z)+.65),
        axes=FALSE,xaxs="i",xlab="Log2 difference",ylab="")
      for(r in seq_along(yy))if(r%%2==1)rect(xr[1],yy[r]-.45,xr[2],yy[r]+.45,col="#F3F6F8",border=NA)
      abline(v=0,col="#A9B4BE",lty=2,lwd=.9)
      colors<-publication_colors[z$tissue]
      segments(z$lower95,yy,z$upper95,yy,col=colors,lwd=1.5)
      points(z$log2FC,yy,bg=colors,col="white",pch=ifelse(z$tissue=="Root",21,24),cex=1.2,lwd=.5)
      ticks<-pretty(xr,n=4);ticks<-ticks[ticks>=xr[1]&ticks<=xr[2]]
      axis(1,at=ticks,lwd=0,lwd.ticks=.6);axis(2,at=yy,labels=labels,tick=FALSE,cex.axis=.79)
      segments(xr[1],.5,xr[2],.5,col="#A9B4BE",lwd=.65)
      text(xr[2]+diff(xr)*.07,yy,formatC(z$q,digits=2,format="g"),adj=0,xpd=NA,cex=.76,col=publication_colors["muted"])
      text(xr[2]+diff(xr)*.07,nrow(z)+.65,"BH q",adj=c(0,0),xpd=NA,cex=.72,font=2,col=publication_colors["muted"])
      publication_title(LETTERS[i],titles[i],subtitles[i])
    }
    publication_key(footnote="Pointwise 95% intervals. Root NIT2 is ineligible; AHA2 specificity is unresolved.")
  })
}

build_revision_prestress_figure <- function(inputs,settings,destination) {
  a<-revision_genotype(inputs$microarray_genotype)
  sets<-list("Genotype microarray"=list(y=a$y,m=a$metadata))
  for(x in Filter(function(x)x$spec$assay=="proteomics",inputs)){
    z<-revision_protein(x);sets[[paste(x$spec$tissue,"proteome")]]<-list(
      y=z$y[z$selection$complete_case,,drop=FALSE],m=z$metadata)
  }
  similarities<-lapply(seq_along(sets),function(i)prestress_splits(sets[[i]]$y,sets[[i]]$m,names(sets)[i]))
  publication_figure(destination,7.2,5.7,function(){
    layout(rbind(c(1,1,1),2:4),heights=c(1,1.65))
    par(oma=c(.6,0,.6,0),cex=1)
    par(mar=c(3.2,8.8,2.7,1))
    plot(NA,xlim=c(-1,1),ylim=c(.5,3.5),xlab="Pearson correlation",ylab="",yaxt="n",axes=FALSE,xaxs="i")
    abline(v=seq(-1,1,by=.25),col=publication_colors["grid"],lwd=.6);abline(v=0,col="#A9B4BE",lty=2)
    axis(1,at=seq(-1,1,by=.25),lwd=0,lwd.ticks=.6)
    axis(2,at=3:1,labels=names(sets),tick=FALSE,cex.axis=.86)
    for(i in seq_along(sets)){
      z<-similarities[[i]];v<-z$Pearson[z$comparison=="disjoint_controls"];r<-z$Pearson[z$comparison=="shared_controls"]
      segments(min(v),4-i,max(v),4-i,col="#687985",lwd=1.3)
      points(v,4-i+seq(-.12,.12,length.out=length(v)),pch=21,bg="white",col="#687985",cex=.95)
      points(r,4-i,pch=18,col=publication_colors["NIT2"],cex=1.45)
    }
    publication_title("A","Sensitivity to control sharing")
    legend("topright",c("Shared controls","Each disjoint split"),pch=c(18,21),
      col=c(publication_colors["NIT2"],"#687985"),pt.bg="white",bty="n",cex=.78,horiz=TRUE,
      inset=c(0,0))
    par(mar=c(3.8,3.8,2.6,.6))
    for(i in seq_along(sets)){
      y<-sets[[i]]$y;m<-sets[[i]]$m
      b<-rowMeans(y[,m$cell=="OE0",drop=FALSE])-rowMeans(y[,m$cell=="WT0",drop=FALSE])
      s<-rowMeans(y[,m$cell=="WTS",drop=FALSE])-rowMeans(y[,m$cell=="WT0",drop=FALSE])
      # All measured features and extremes remain visible on symmetric axes.
      lim<-max(abs(c(b,s)))*1.06
      publication_axes(c(-lim,lim),c(-lim,lim),"No-0 salt response (log2)",
        if(i==1)"Basal genotype difference (log2)" else "")
      abline(v=0,h=0,col="#A9B4BE",lty=2,lwd=.7)
      points(s,b,pch=16,cex=.27,col=adjustcolor(publication_colors["No-0"],alpha.f=.20))
      publication_title(LETTERS[i+1L],names(sets)[i],paste(format(nrow(y),big.mark=","),"measured features"))
    }

  })
}

build_publication_qpcr <- function(inputs,settings,destination) {
  d<-inputs$qpcr_summary$data;genes<-c("AHA2","SOS1","HKT1","HAK5")
  expected<-expand.grid(Treatment=settings$qpcr_treatment_labels,Genotype=settings$qpcr_genotype_labels,
    Tissue=settings$tissue_levels,stringsAsFactors=FALSE)
  key<-function(z)paste(z$Tissue,z$Genotype,z$Treatment)
  for(g in genes){z<-d[d$Gene==g,];if(nrow(z)!=nrow(expected)||!setequal(key(z),key(expected)))stop("Incomplete qPCR candidate panel")}
  publication_figure(destination,7.2,5.2,function(){
    par(mfrow=c(2,4),mar=c(3,3.6,2.7,.5),oma=c(2,.1,1.8,0))
    for(i in seq_along(settings$tissue_levels))for(j in seq_along(genes)){
      z<-d[d$Gene==genes[j]&d$Tissue==settings$tissue_levels[i],];v<-z[["Average Normalized expression"]];se<-z$STE
      limits<-range(c(0,v,v-se,v+se),na.rm=TRUE);limits[2]<-max(limits[2],.001)*1.15
      if(limits[1]<0)limits[1]<-limits[1]*1.1
      publication_axes(c(-.3,1.3),limits,"NaCl (mM)",if(j==1)"Normalized expression" else "",
        at=c(0,1),labels=c("0","75"))
      for(k in seq_len(nrow(z))){
        genotype<-if(z$Genotype[k]=="No-0")"No-0" else "NIT2";salt<-z$Treatment[k]!="Untreated"
        xx<-as.integer(salt)+if(genotype=="No-0")-.09 else .09;color<-publication_colors[genotype]
        if(is.finite(se[k])&&se[k]>0){
          segments(xx,v[k]-se[k],xx,v[k]+se[k],col=color,lwd=1)
          segments(xx-.035,v[k]-se[k],xx+.035,v[k]-se[k],col=color,lwd=1)
          segments(xx-.035,v[k]+se[k],xx+.035,v[k]+se[k],col=color,lwd=1)
        }
        points(xx,v[k],bg=color,col="white",lwd=.5,pch=if(salt)24 else 21,cex=1.2)
      }
      publication_title(LETTERS[(i-1)*4+j],genes[j],settings$tissue_levels[i])
    }
    publication_key(FALSE,"Group means +/- source STE; short bars may sit inside symbols. n and STE definition unknown.")
  })
}
