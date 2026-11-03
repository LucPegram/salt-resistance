# Compare two completed analysis releases without fitting or changing models.
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),'bootstrap.R'))
required<-c('SALT_BASELINE_REVISION_ROOT','SALT_BASELINE_SUPPORTING_ROOT','SALT_OUTPUT_ROOT','SALT_COMPARISON_ROOT')
roots<-Sys.getenv(required)
if(any(!nzchar(roots)))stop('Set ',paste(required,collapse=', '))
roots[1:3]<-vapply(roots[1:3],normalizePath,'',mustWork=TRUE)
roots[4]<-normalizePath(roots[4],mustWork=FALSE)
inside<-function(a,b)identical(a,b)||startsWith(a,paste0(b,'/'))
if(any(vapply(c(project_root,roots[1:3]),function(r)inside(roots[4],r)||inside(r,roots[4]),FALSE)))stop('Separate comparison output from code and analysis inputs')
dir.create(roots[4],recursive=TRUE,showWarnings=FALSE)
used<-character()
read<-function(root,id,file){
 path<-file.path(root,id,file)
 if(!file.exists(file.path(root,id,'RUN_COMPLETE')))stop('Incomplete artifact: ',id)
 used<<-unique(c(used,path));read.csv(path,check.names=FALSE,stringsAsFactors=FALSE)
}
a<-read(roots[1],'table_revision_genotype','genotype_contrasts.csv')
b<-read(roots[3],'table_revision_genotype','genotype_contrasts.csv')
compare<-function(a,b,label,group,key='feature',effect=TRUE){
 stopifnot(!anyDuplicated(a[[key]]),!anyDuplicated(b[[key]]),setequal(a[[key]],b[[key]]))
 b<-b[match(a[[key]],b[[key]]),];sa<-a$q<.05;sb<-b$q<.05
 union<-sum(sa|sb)
 data.frame(dataset=label,test=group,tested=nrow(a),supplied_scaling_significant=sum(sa),quantile_significant=sum(sb),
 retained=sum(sa&sb),lost=sum(sa&!sb),gained=sum(!sa&sb),Jaccard=if(union)sum(sa&sb)/union else NA_real_,
 effect_Pearson=if(effect)cor(a$log2FC,b$log2FC) else NA_real_,
 median_absolute_log2FC_change=if(effect)median(abs(a$log2FC-b$log2FC)) else NA_real_,
 significant_both_direction_reversals=if(effect)sum(sa&sb&sign(a$log2FC)!=sign(b$log2FC)) else NA_integer_)
}
rows<-lapply(sort(unique(a$contrast)),function(g)compare(a[a$contrast==g,],b[b$contrast==g,],'Array2',g))
a<-read(roots[1],'table_revision_dose','dose_joint_tests.csv');b<-read(roots[3],'table_revision_dose','dose_joint_tests.csv')
rows<-c(rows,lapply(sort(unique(a$test)),function(g)compare(a[a$test==g,],b[b$test==g,],'Array1',g,effect=FALSE)))
write_table(do.call(rbind,rows),file.path(roots[4],'normalization_impact.csv'))
a<-read(roots[2],'table_supporting_ranked','supporting_ranked.csv')
b<-read(roots[3],'table_supporting_ranked','supporting_ranked.csv')
key<-function(d)paste(d$comparison,d$tissue,d$contrast,d$locus,sep='|')
a$key<-key(a);b$key<-key(b)
fields<-c('key','Gene_name','locus','tissue','contrast','RNA_log2FC','RNA_q','protein_log2FC','protein_q','profile_Pearson')
d<-merge(a[fields],b[fields],by='key',all=TRUE,suffixes=c('_supplied_scaling','_quantile'))
d$status<-ifelse(is.na(d$locus_supplied_scaling),'gained',ifelse(is.na(d$locus_quantile),'lost','retained'))
write_table(d,file.path(roots[4],'supporting_candidate_changes.csv'))
used_hash<-setNames(vapply(used,function(p)digest::digest(file=p,algo='sha256'),''),used)
outputs<-list.files(roots[4],pattern='\\.csv$',full.names=TRUE)
code<-source_files(project_root)
dput(list(generated_UTC=format(Sys.time(),tz='UTC',usetz=TRUE),input_sha256=used_hash,
 git_revision=system2('git',c('-C',shQuote(project_root),'rev-parse','HEAD'),stdout=TRUE),
 code_sha256=setNames(vapply(code,function(p)digest::digest(file=p,algo='sha256'),''),substring(code,nchar(project_root)+2L)),
 output_sha256=setNames(vapply(outputs,function(p)digest::digest(file=p,algo='sha256'),''),basename(outputs)),
 interpretation='Sensitivity to between-array normalization; no method chosen using these discovery counts.'),
 file=file.path(roots[4],'provenance.R'))
writeLines(capture.output(sessionInfo()),file.path(roots[4],'sessionInfo.txt'))
writeLines('Complete; see provenance.R.',file.path(roots[4],'RUN_COMPLETE'))
