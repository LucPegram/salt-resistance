build_integration_concordance_figure <- function(inputs,settings,destination) {
  z<-integration_results(inputs);d<-z$table[z$table$method=="primary",]
  specs<-data.frame(comparison=c("No0_genotype","No0_genotype","Col0_context"),
    contrast=c("baseline_genotype","salt_WT","salt_WT"),
    title=c("Basal genotype","No-0 salt","Col-0 context"),
    subtitle=c("Combined RNA / tissue protein","RNA 100 mM / protein 75 mM","RNA Col-0 / protein No-0; 75 mM"))
  targets<-integration_settings()$candidate_loci
  publication_figure(destination,7.2,6.5,function(){
    par(mfrow=c(2,3),mar=c(3.8,4.0,3.4,.6),oma=c(2.0,0,.1,0))
    for(i in 1:2)for(j in 1:3){
      tissue<-c("Root","Shoot")[i]
      a<-d[d$tissue==tissue&d$comparison==specs$comparison[j]&d$contrast==specs$contrast[j],]
      column<-d[d$comparison==specs$comparison[j]&d$contrast==specs$contrast[j],]
      publication_axes(publication_range(c(0,column$RNA_log2FC),.15),
        publication_range(c(0,column$protein_log2FC),.18),"RNA log2 difference","Protein log2 difference")
      abline(v=0,h=0,col="#A9B4BE",lty=2,lwd=.8)
      joint<-a$joint_q_global<.05
      points(a$RNA_log2FC[!joint],a$protein_log2FC[!joint],pch=16,cex=.45,col=adjustcolor("#687985",alpha.f=.3))
      points(a$RNA_log2FC[joint],a$protein_log2FC[joint],pch=21,bg="#009E73",col="white",lwd=.4,cex=.85)
      hit<-which(a$locus%in%targets)
      if(length(hit)){
        points(a$RNA_log2FC[hit],a$protein_log2FC[hit],pch=ifelse(joint[hit],23,21),
          bg="#E69F00",col=ifelse(joint[hit],publication_colors["ink"],"white"),cex=1.05,lwd=.6)
      }
      labelled<-which(a$locus%in%targets|joint)
      label<-names(targets)[match(a$locus[labelled],targets)]
      label[is.na(label)]<-a$locus[labelled][is.na(label)]
      publication_pool_labels(a$RNA_log2FC[labelled],a$protein_log2FC[labelled],label)
      publication_title(LETTERS[(i-1)*3+j],paste(specs$title[j],tissue),
        sprintf("%s genes; r = %.2f; rho = %.2f",format(nrow(a),big.mark=","),
          cor(a$RNA_log2FC,a$protein_log2FC),cor(a$RNA_log2FC,a$protein_log2FC,method="spearman")))
      publication_caption(paste0("Panels ",LETTERS[j]," and ",LETTERS[j+3],": ",specs$subtitle[j],"."))
    }
    publication_caption("Joint global BH < 0.05: green points or amber diamonds. Other candidates: amber circles.")
    publication_caption("Cross-experiment effects; correlations are descriptive, without gene-as-replicate significance tests.")
  })
}

build_integration_pathway_figure <- function(inputs,settings,destination) {
  all<-integration_pathways(inputs)$table
  d<-all[all$method=="estimated_correlation",]
  gos<-integration_settings()$pathway_GO
  labels<-vapply(gos,function(g)unique(d$term[d$GO==g])[1],"")
  # Short display names retain the full GO identifier/term in the result table.
  short<-c("Salt stress","Sodium transport","Potassium transport","Proton transport",
    "Auxin response","ROS response","Oxidative stress","H2O2 catabolism","Antioxidant activity")
  cols<-c(integration_settings()$contrasts,"Col0_context")
  xlabs<-c("Basal\ngeno.","No-0\nsalt","NIT2\nsalt","Salt\ngeno.","Geno.\nby salt","Col-0\n75 mM")
  publication_figure(destination,7.2,4.7,function(){
    layout(matrix(1:2,nrow=1),widths=c(1.5,1))
    par(oma=c(2.5,0,.1,0),mar=c(3.4,7.8,2.6,.5),cex=1)
    for(i in 1:2){
      tissue<-c("Root","Shoot")[i]
      if(i==2)par(mar=c(3.4,1.0,2.6,.6))
      plot(NA,xlim=c(.45,6.55),ylim=c(.45,9.55),axes=FALSE,xaxs="i",yaxs="i",xlab="",ylab="")
      for(row in 1:9)if(row%%2==1)rect(.45,10-row-.45,6.55,10-row+.45,col="#F3F6F8",border=NA)
      if(i==1)axis(2,at=9:1,labels=short,tick=FALSE,cex.axis=.78)
      for(j in 1:6){
        comparison<-if(j==6)"Col0_context" else "No0_genotype"
        co<-if(j==6)"salt_WT" else cols[j]
        for(k in 1:9){
          a<-d[d$tissue==tissue&d$comparison==comparison&d$contrast==co&d$GO==gos[k],]
          stopifnot(nrow(a)==1L)
          if(a$status!="tested"){points(j,10-k,pch=4,col="#B3BDC5",cex=.65);next}
          direction<-if(a$RNA_direction==a$protein_direction)a$RNA_direction else "Opposed"
          color<-c(Up="#D55E00",Down="#0072B2",Opposed="#A9B4BE")[direction]
          significant<-a$joint_q_global<.05
          size<-.85+min(4,-log10(max(a$joint_P,1e-12)))*.28
          points(j,10-k,pch=21,bg=if(significant)color else "white",col=color,cex=size,lwd=1.1)
        }
        parts<-strsplit(xlabs[j],"\n",fixed=TRUE)[[1]]
        for(k in seq_along(parts))mtext(parts[k],side=1,at=j,line=.4+(k-1)*.85,cex=.68)
      }
      publication_title(LETTERS[i],paste(tissue,"protein integration"),"RNA combined except Col-0 context")
    }
    publication_caption("Red: higher relative enrichment in both assays; blue: lower; grey: opposite directions.")
    publication_caption("Filled: joint global BH < 0.05; size: -log10 joint P (capped at 4); cross: too few measured members.")
  })
}
