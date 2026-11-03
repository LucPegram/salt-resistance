script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),'..','scripts','bootstrap.R'))
# Positive profile alone is not enough; both assays must be significant in the same direction.
d<-data.frame(comparison='No0_genotype',RNA_q=c(.01,.01,.01,.2,NA,.05),protein_q=c(.02,.02,.02,.001,.001,.001),
 RNA_log2FC=1,protein_log2FC=c(1,-1,1,1,1,1),profile_Pearson=c(.7,.7,-.2,.8,.9,.7),
 joint_q_global=.2,joint_q_family=.1)
f<-supporting_flag(d)
stopifnot(identical(which(f$selected_support),1L),f$assaywise_supported[3],!f$profile_positive[3],
 !f$assaywise_supported[6])
# Four condition means have equal weight; correlation is not a sample-level test.
x<-rbind(c(1,2,3,4),c(1,2,3,4),rep(2,4));y<-rbind(c(3,5,7,9),c(4,3,2,1),c(1,2,3,4))
p<-supporting_profile(x,y)
stopifnot(abs(p$profile_Pearson[1]-1)<1e-12,abs(p$profile_Pearson[2]+1)<1e-12,
 is.na(p$profile_Pearson[3]),!any(grepl('P_value|pval',names(p))))
m<-data.frame(cell=factor(rep(c('WT0','WTS','OE0','OES'),each=2),levels=c('WT0','WTS','OE0','OES')))
y<-matrix(1:16,nrow=2);rownames(y)<-c('a','b');means<-supporting_cells(y,m)
stopifnot(all(means[1,]==c(2,6,10,14)),identical(colnames(means),supporting_settings()$profile_cells))
# A salt hit in roots and a genotype hit in shoots must not become a dual-tissue claim.
d<-data.frame(Gene_name=c('A','A','B','B'),locus=c('AT1G00010','AT1G00010','AT1G00020','AT1G00020'),
 tissue=c('Root','Shoot','Root','Root'),protein_group=c('X','Y','Z','Z'),
 contrast=c('salt_WT','saline_genotype','salt_WT','saline_genotype'),
 balanced_assay_q=c(.01,.01,.02,.03),profile_Pearson=.8,profile_Spearman=.8,
 RNA_log2FC=1,protein_log2FC=2,RNA_q=.02,protein_q=.03,joint_q_global=.5)
a<-supporting_dual(d)
stopifnot(nrow(a)==1L,a$locus=='AT1G00020',a$tissue=='Root',a$worst_assay_q==.03,
 !a$both_contrasts_joint_global_BH)
# Significance values stay as calculated in their original universe after ranking.
d$smaller_abs_effect<-1;before<-d$RNA_q;a<-supporting_order(d)
stopifnot(setequal(a$RNA_q,before),identical(a$balanced_assay_q,sort(d$balanced_assay_q)))
cat('Supporting-trend criteria, missing/boundary handling, means/correlation and same-tissue dual support checks passed.\n')
