# Select display rows from complete reviewer-response outputs. No fitted model,
# significance cutoff or source adjustment is recomputed on the display subset.
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),'bootstrap.R'))
ids<-c('table_probe_audit','table_pair_influence','table_supporting_candidates','table_strict_probe_baseline','table_TF_correlations','table_TF_overlap','table_TF_target_members','figure_TF_correlation_sensitivity')
args<-commandArgs(TRUE)
if(length(args)!=1L||!args%in%c('all','--list',ids))stop('Provide one table ID, all or --list')
if(args=='--list'){cat(ids,sep='\n');quit(status=0)}
roots<-c(primary=Sys.getenv('SALT_PRIMARY_RESULTS'),probe=Sys.getenv('SALT_PROBE_RESULTS'),
 qc=Sys.getenv('SALT_PAIR_RESULTS'),tf=Sys.getenv('SALT_TF_RESULTS'),extension=Sys.getenv('SALT_EXTENSION_ROOT'))
out<-Sys.getenv('SALT_OUTPUT_ROOT')
if(any(!nzchar(c(roots,out))))stop('Set all result roots and output root')
roots<-vapply(roots,normalizePath,'',mustWork=TRUE)
inside<-function(a,b)identical(a,b)||startsWith(a,paste0(b,'/'))
dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out,mustWork=TRUE)
if(any(vapply(c(roots,project_root),function(p)inside(out,p)||inside(p,out),TRUE)))stop('Separate source, code and output roots')
sources<-c(
 audit=file.path(roots['probe'],'probe_feature_audit/probe_feature_audit.csv'),
 reconstruction=file.path(roots['probe'],'probe_reconstruction_audit/probe_reconstruction_audit.csv'),
 metrics=file.path(roots['qc'],'dose_omission_metrics/dose_omission_metrics.csv'),
 dual=file.path(roots['primary'],'table_supporting_dual/supporting_dual.csv'),
 three=file.path(roots['primary'],'table_supporting_three/supporting_three_datasets.csv'),
 baseline=file.path(roots['probe'],'probe_sensitive_external_tests/probe_sensitive_external_tests.csv'),
 tf=file.path(roots['tf'],'TF_estimated_correlation_tests/TF_estimated_correlation_tests.csv'),
 fixed=file.path(roots['extension'],'TF_target_enrichment/TF_target_enrichment.csv'))
if(any(!file.exists(sources))||any(!file.exists(file.path(dirname(sources),'RUN_COMPLETE'))))stop('Incomplete source outputs')
sources<-c(sources,sequences=file.path(roots['probe'],'reconstructed_complete_probes.fasta'))
if(!file.exists(sources['sequences'])||!file.exists(file.path(roots['probe'],'RUN_COMPLETE')))stop('Incomplete reconstructed probe sequences')
public_root<-normalizePath(Sys.getenv('SALT_REFERENCE_ROOT'),mustWork=TRUE)
if(inside(out,public_root)||inside(public_root,out))stop('Separate public references and outputs')
sources<-c(sources,DAP=file.path(public_root,'AtRegNet.zip'),
 universe=file.path(roots['extension'],'dose_joint_tests/dose_joint_tests.csv'))
if(any(!file.exists(sources)))stop('Missing DAP reference or measured-universe table')
hash<-function(p)digest::digest(file=p,algo='sha256')
before<-vapply(sources,hash,'');code<-source_files(project_root)
code_before<-setNames(vapply(code,hash,''),substring(code,nchar(project_root)+2L))
d<-lapply(sources[!names(sources)%in%c('sequences','DAP')],read.csv,stringsAsFactors=FALSE,check.names=FALSE)
seqs<-as.character(Biostrings::readDNAStringSet(sources['sequences']))
a<-d$audit
probe<-rbind(data.frame(item='Complete consistent queries',count=d$reconstruction$value[d$reconstruction$quantity=='queries_with_complete_consistent_sequence'],unit='query identifiers'),
 data.frame(item=names(table(a$status)),count=as.integer(table(a$status)),unit='ordinary physical features'),
 data.frame(item='Strict eligible probe features',count=sum(a$strict_probe_sensitivity_eligible),unit='ordinary physical features'))
focal<-c(NIT2='AT3G44300',SOS1='AT2G01980',AHA2='AT4G30190',KUP4='AT4G23640',HKT1='AT4G10310',HAK5='AT4G13420',`CDKA;1`='AT3G48750',AT5G44720='AT5G44720',`LRK10L-1.2`='AT1G18390')
probe<-rbind(probe,do.call(rbind,lapply(seq_along(focal),function(i)data.frame(item=paste(names(focal)[i],focal[i],sep=' / '),
 count=sum(a$supplied_locus==focal[i]&a$strict_probe_sensitivity_eligible),unit='strict eligible physical features'))))
probe$distinct_query_ids<-probe$distinct_sequences<-NA_integer_
for(i in seq_along(focal)){
 queries<-unique(a$probe[a$supplied_locus==focal[i]&a$strict_probe_sensitivity_eligible])
 stopifnot(all(queries%in%names(seqs)))
 j<-match(paste(names(focal)[i],focal[i],sep=' / '),probe$item)
 probe$distinct_query_ids[j]<-length(queries)
 probe$distinct_sequences[j]<-length(unique(seqs[queries]))
}
co<-c('Root_75v0','Root_75v50','Root_100v75','Shoot_75v0','Difference_75v0','Difference_75v50','Difference_100v75')
q<-d$metrics[d$metrics$omitted_pool=='75_2',]
f<-q[q$normalization=='frozen',];r<-q[q$normalization=='refitted',]
f<-f[match(co,f$contrast),];r<-r[match(co,r$contrast),]
qc<-data.frame(contrast=co,frozen_Pearson=f$Pearson,refitted_Pearson=r$Pearson,
 refitted_median_absolute_effect_change=r$median_absolute_effect_change,
 refitted_p95_absolute_effect_change=r$p95_absolute_effect_change,
 refitted_sign_agreement=r$sign_agreement_among_effects_ge_05)
support<-rbind(data.frame(Gene=d$dual$Gene_name,locus=d$dual$locus,tissue=d$dual$tissue,screen='Salt and saline genotype',
 RNA_salt_log2FC=d$dual$salt_RNA_log2FC,protein_salt_log2FC=d$dual$salt_protein_log2FC,
 worst_assay_q=d$dual$worst_assay_q,joint_global_BH=d$dual$salt_joint_q_global),
 data.frame(Gene=d$three$Gene_name,locus=d$three$locus,tissue=d$three$tissue,screen='Three-dataset salt context',
 RNA_salt_log2FC=d$three$Array2_log2FC,protein_salt_log2FC=d$three$protein_log2FC,
 worst_assay_q=d$three$worst_assay_q,joint_global_BH=d$three$Array2_joint_q_global))
b<-d$baseline[d$baseline$contrast=='baseline_genotype',];b<-b[order(b$dataset,b$reference_tissue,as.numeric(sub('.*hour','',b$reference))),]
tf<-bmc_tf_comparison(d$fixed,d$tf)
stopifnot(all(tf$dataset=='Array1'),all(tf$contrast=='Shoot_100v0'))
dap<-extension_read_dap(sources['DAP']);universe<-unique(d$universe$locus)
selected<-unique(tf$TF);sets<-lapply(selected,function(id)intersect(dap$edges$target[dap$edges$TF==id],universe));names(sets)<-selected
expected<-tf$targets[match(selected,tf$TF)];stopifnot(identical(unname(lengths(sets)),as.integer(expected)))
overlap<-bmc_target_overlap(sets)
overlap$TF1_name<-tf$TF_name[match(overlap$TF1,tf$TF)];overlap$TF2_name<-tf$TF_name[match(overlap$TF2,tf$TF)]
members<-do.call(rbind,lapply(selected,function(id)data.frame(TF=id,TF_name=tf$TF_name[match(id,tf$TF)],locus=sets[[id]],dataset='Array1')))
values<-list(table_probe_audit=probe,table_pair_influence=qc,table_supporting_candidates=support,
 table_strict_probe_baseline=b,table_TF_correlations=tf,table_TF_overlap=overlap,table_TF_target_members=members)
revision<-system2('git',c('-C',shQuote(project_root),'rev-parse','HEAD'),stdout=TRUE)
for(id in if(args=='all')ids else args){
 folder<-file.path(out,id);dir.create(folder,showWarnings=FALSE)
 marker<-file.path(folder,'RUN_COMPLETE');if(file.exists(marker))unlink(marker)
 if(id=='figure_TF_correlation_sensitivity'){
  path<-file.path(folder,paste0(id,'.pdf'));bmc_tf_comparison_figure(tf,path,overlap)
  outputs<-list.files(folder,pattern='\\.(pdf|svg|png|tiff)$|_caption.txt$',full.names=TRUE)
 }else{
  path<-file.path(folder,paste0(id,'.csv'));write_table(values[[id]],path);outputs<-path
 }
 stopifnot(identical(before,vapply(sources,hash,'')),identical(unname(code_before),unname(vapply(code,hash,''))))
 dput(list(artifact_id=id,git_revision=revision,input_sha256=before,code_sha256=code_before,
  output_sha256=setNames(vapply(outputs,hash,''),basename(outputs)),command=commandArgs(),session=capture.output(sessionInfo()),
  specification='Display selection only. Complete source-family adjustments preserved; no correction recalculated on selected rows.'),file=file.path(folder,'provenance.R'))
 writeLines('Complete; see provenance.R',marker);message('Completed ',id)
}
