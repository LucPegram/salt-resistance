revision_pdf <- function(path,width=8,height=6) {
  grDevices::pdf(path,width=width,height=height,useDingbats=FALSE)
  graphics::par(family="sans",las=1)
}

build_revision_qc_figure <- function(inputs,settings,destination) {
  revision_pdf(destination,9,8);on.exit(dev.off(),add=TRUE)
  par(mfrow=c(3,2),mar=c(4,4,3,1))
  for(id in names(inputs)){
    x<-inputs[[id]];d<-revision_qc(x)
    if(x$spec$assay=="proteomics" || identical(x$spec$tissue,"Combined")){
      m<-factorial_metadata(x)
      color<-ifelse(m$genotype=="No-0","#4477AA","#CC6677")
      shape<-ifelse(m$salt_mM==0,1,16)
    }else{
      dose<-as.numeric(ifelse(startsWith(x$samples,"PN"),"0",sub(" .*","",x$samples)))
      color<-c("#777777","#4477AA","#228833","#CC6677")[match(dose,c(0,50,75,100))]
      shape<-rep(16,length(color))
    }
    plot(d$PC1,d$PC2,col=color,pch=shape,cex=1.1,
      xlab=sprintf("PC1 (%.1f%%)",100*d$PC1_fraction[1]),
      ylab=sprintf("PC2 (%.1f%%)",100*d$PC2_fraction[1]),
      main=paste(LETTERS[match(id,names(inputs))],x$spec$tissue,
        if(x$spec$assay=="proteomics")"proteome" else "microarray"))
    text(d$PC1,d$PC2,labels=sub(".*[ _]","",x$samples),pos=3,cex=.65)
    if(x$spec$assay=="proteomics" || identical(x$spec$tissue,"Combined"))
      legend("topright",c("No-0 control","No-0 salt","NIT2 control","NIT2 salt"),
        col=rep(c("#4477AA","#CC6677"),each=2),pch=rep(c(1,16),2),bty="n",cex=.65)
    else legend("topright",paste(c(0,50,75,100),"mM"),col=c("#777777","#4477AA","#228833","#CC6677"),pch=16,bty="n",cex=.65)
  }
  plot.new();text(.05,.9,"PCA uses fully observed features",adj=0,cex=.9)
  text(.05,.78,"No imputation and no sample exclusions",adj=0,cex=.9)
  text(.05,.66,"Labels identify biological pool replicates",adj=0,cex=.9)
}

revision_candidate_transcript <- function(y,locus) {
  hits<-grep(paste0("^",locus,"(\\.|$)"),rownames(y),value=TRUE)
  if(!length(hits))stop("Prespecified locus absent: ",locus)
  # Stable, documented selection, independent of observed effects or P values.
  if(paste0(locus,".1")%in%hits)paste0(locus,".1") else sort(hits)[1]
}

build_revision_dose_figure <- function(inputs,settings,destination) {
  out<-revision_dose(inputs)
  loci<-c(NIT2="AT3G44300",SOS1="AT2G01980",AHA2="AT4G30190",HKT1="AT4G10310")
  revision_pdf(destination,8,7);on.exit(dev.off(),add=TRUE)
  par(mfrow=c(2,2),mar=c(4,4,3,1))
  for(g in names(loci)){
    id<-revision_candidate_transcript(out$y,loci[g]);v<-out$y[id,];m<-out$metadata
    plot(NA,xlim=c(-5,105),ylim=range(v),xlab="NaCl (mM)",
      ylab="Log2 supplied normalized signal",main=paste(LETTERS[match(g,names(loci))],g,id),cex.main=.85)
    for(t in c("Root","Shoot")){
      k<-m$tissue==t;col<-if(t=="Root")"#4477AA" else "#CC6677"
      offset<-if(t=="Root")-1.5 else 1.5
      points(m$dose[k]+offset+(m$replicate[k]-2)*.6,v[k],col=col,pch=16,cex=.8)
      means<-tapply(v[k],m$dose[k],mean);lines(as.numeric(names(means))+offset,means,col=col,lwd=1.5)
    }
    legend("topleft",c("Root","Shoot"),col=c("#4477AA","#CC6677"),pch=16,bty="n",cex=.7)
  }
}

build_revision_protein_figure <- function(inputs,settings,destination) {
  revision_pdf(destination,9,6.5);on.exit(dev.off(),add=TRUE)
  par(mfrow=c(2,3),mar=c(4.7,4.2,3,1))
  for(x in inputs){z<-revision_protein(x)
    for(g in c("NIT2","NHX7","AHA2")){
      k<-which(z$annotations$Genes==g)
      if(length(k)!=1L)stop("Ambiguous candidate protein")
      y<-z$y[k,];m<-z$metadata
      plot(NA,xlim=c(.5,4.5),ylim=range(y,na.rm=TRUE),xaxt="n",
        xlab="",ylab="Log2 supplied abundance",
        main=paste(LETTERS[match(g,c("NIT2","NHX7","AHA2"))+if(x$spec$tissue=="Shoot")3 else 0],
          if(g=="NHX7")"SOS1" else if(g=="AHA2")"AHA2-labelled group" else g,x$spec$tissue),cex.main=.85)
      axis(1,at=1:4,labels=c("WT 0","WT 75","NIT2 0","NIT2 75"),cex.axis=.72)
      colors<-c("#4477AA","#4477AA","#CC6677","#CC6677")
      for(j in 1:4){v<-y[m$cell==levels(m$cell)[j]]
        points(j+c(-.11,0,.11),v,pch=16,col=colors[j],cex=.85)
        if(any(is.finite(v)))segments(j-.23,mean(v,na.rm=TRUE),j+.23,mean(v,na.rm=TRUE),lwd=1.8)
      }
      mtext(paste0(z$annotations$Protein.Group[k],"; ",sum(!is.na(y)),"/12 observed"),cex=.7)
    }
  }
}

build_revision_contrast_figure <- function(inputs,settings,destination) {
  d<-do.call(rbind,lapply(inputs,function(x)revision_protein(x)$table))
  d<-d[d$method=="observed_min2"&d$Genes%in%c("NIT2","NHX7","AHA2"),]
  revision_pdf(destination,8,7);on.exit(dev.off(),add=TRUE)
  par(mfrow=c(2,2),mar=c(4,5.2,3,1))
  for(co in c("baseline_genotype","salt_WT","saline_genotype","interaction")){
    z<-d[d$contrast==co,];z<-z[order(match(z$Genes,c("NIT2","NHX7","AHA2")),z$tissue),]
    labels<-paste(ifelse(z$Genes=="NHX7","SOS1",z$Genes),z$tissue)
    yy<-rev(seq_len(nrow(z)))
    plot(z$log2FC,yy,xlim=range(c(0,z$lower95,z$upper95)),ylim=c(.5,nrow(z)+.5),
      yaxt="n",ylab="",xlab="Log2 contrast with pointwise 95% interval",
      main=paste(LETTERS[match(co,c("baseline_genotype","salt_WT","saline_genotype","interaction"))],
        gsub("_"," ",co)),cex.main=.9,pch=16,col=ifelse(z$tissue=="Root","#4477AA","#CC6677"))
    axis(2,at=yy,labels=labels,las=1,cex.axis=.75);abline(v=0,lty=2,col="gray60")
    segments(z$lower95,yy,z$upper95,yy,col=ifelse(z$tissue=="Root","#4477AA","#CC6677"))
    text(z$log2FC,yy+.23,labels=paste0("q=",formatC(z$q,digits=2,format="g")),cex=.6)
  }
}

build_revision_prestress_figure <- function(inputs,settings,destination) {
  sets<-list()
  a<-revision_genotype(inputs$microarray_genotype)
  sets[["Genotype microarray"]]<-list(y=a$y,m=a$metadata)
  for(x in Filter(function(x)x$spec$assay=="proteomics",inputs)){
    z<-revision_protein(x);sets[[paste(x$spec$tissue,"proteome")]]<-list(
      y=z$y[z$selection$complete_case,,drop=FALSE],m=z$metadata)
  }
  revision_pdf(destination,10,4);on.exit(dev.off(),add=TRUE)
  par(mfrow=c(1,3),mar=c(4.5,4.2,3,1))
  for(n in names(sets)){
    y<-sets[[n]]$y;m<-sets[[n]]$m
    b<-rowMeans(y[,m$cell=="OE0",drop=FALSE])-rowMeans(y[,m$cell=="WT0",drop=FALSE])
    s<-rowMeans(y[,m$cell=="WTS",drop=FALSE])-rowMeans(y[,m$cell=="WT0",drop=FALSE])
    plot(s,b,pch=16,cex=.24,col=adjustcolor("#4477AA",alpha.f=.25),
      xlab="WT salt response (log2)",ylab="Baseline genotype difference (log2)",
      main=paste(LETTERS[match(n,names(sets))],n),cex.main=.85)
    abline(h=0,v=0,col="gray65",lty=2)
    z<-prestress_splits(y,m,n);v<-z$Pearson[z$comparison=="disjoint_controls"]
    legend("topleft",c(sprintf("Shared control r = %.2f",cor(b,s)),
      sprintf("Disjoint control r: %.2f to %.2f",min(v),max(v))),bty="n",cex=.65)
  }
}
