script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),'bootstrap.R'))
source(file.path(project_root,'config','extensions.R'))
if(length(commandArgs(TRUE)))stop('No arguments; set external roots')
root<-normalizePath(Sys.getenv('SALT_DATA_ROOT'),mustWork=TRUE)
raw<-normalizePath(Sys.getenv('SALT_RAW_ROOT'),mustWork=TRUE)
refs<-normalizePath(Sys.getenv('SALT_PROBE_REFERENCE_ROOT'),mustWork=TRUE)
out<-normalizePath(Sys.getenv('SALT_OUTPUT_ROOT'),mustWork=FALSE)
if(!nzchar(out))stop('Set SALT_OUTPUT_ROOT')
inside<-function(a,b)identical(a,b)||startsWith(a,paste0(b,'/'))
if(any(vapply(c(root,raw,refs,project_root),function(x)inside(out,x)||inside(x,out),TRUE)))stop('Keep inputs, references, code and output separate')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
workbook<-file.path(raw,'Microarray_for_Stephano/Probes renamed (by Stu).xlsx')
reference<-file.path(refs,'TAIR10_release63_CDS.fa.gz')
fe_paths<-normalization_paths(raw);fe<-fe_paths[grepl('BG PN 1 -',basename(fe_paths),fixed=TRUE)][1]
ids<-c('microarray_root','microarray_shoot','microarray_genotype','proteomics_root','proteomics_shoot','saved_proteomics_root','saved_proteomics_shoot')
input_paths<-setNames(file.path(root,vapply(input_specs[ids],`[[`,'','path')),ids)
public<-extension_paths()
used<-c(input_paths,probe_workbook=workbook,feature_layout=fe,CDS_reference=reference,public[c('GSE5620_series_matrix.txt.gz','GSE5623_series_matrix.txt.gz','GSE223207_series_matrix.txt.gz','GPL198.soft','GPL17133.soft')])
hash<-function(p)digest::digest(file=p,algo='sha256')
before<-vapply(used,hash,'');code_paths<-source_files(project_root)
code_before<-setNames(vapply(code_paths,hash,''),substring(code_paths,nchar(project_root)+2L))
inputs<-setNames(lapply(ids,read_input,data_root=root,specs=input_specs),ids)
e<-bmc_probe_evidence(workbook);message('Recovered complete query sequences')
m<-bmc_exact_cds_map(e$probes,reference);message('Completed independent exact CDS mapping')
f<-normalization_feature_file(fe);f<-f[f$ControlType==0,]
a<-bmc_probe_audit(e,m,f,inputs$microarray_root$ids)
external<-extension_external(public)
sensitivity<-bmc_probe_sensitivity(inputs,a,external);message('Completed strict probe sensitivity')
summary<-as.data.frame(table(a$status),stringsAsFactors=FALSE);names(summary)<-c('status','physical_features')
q<-e$probes;dna<-Biostrings::DNAStringSet(q$query_sequence);names(dna)<-q$probe
fasta<-file.path(out,'reconstructed_complete_probes.fasta');Biostrings::writeXStringSet(dna,fasta)
artifacts<-c(list(probe_reconstruction_audit=e$audit,probe_reference_audit=m$reference,probe_feature_audit=a,
 probe_status_summary=summary,probe_exact_CDS_hits=m$hits,probe_supplied_near_matches=e$near_edges),sensitivity)
revision<-system2('git',c('-C',shQuote(project_root),'rev-parse','HEAD'),stdout=TRUE)
for(id in names(artifacts)){
 folder<-file.path(out,id);dir.create(folder,showWarnings=FALSE)
 marker<-file.path(folder,'RUN_COMPLETE');if(file.exists(marker))unlink(marker)
 path<-file.path(folder,paste0(id,'.csv'));write_table(artifacts[[id]],path)
 stopifnot(identical(before,vapply(used,hash,'')),identical(unname(code_before),unname(vapply(code_paths,hash,''))))
 provenance<-list(artifact_id=id,git_revision=revision,code_sha256=code_before,input_sha256=before,
  output_sha256=setNames(hash(path),basename(path)),reconstructed_probes_fasta_sha256=hash(fasta),
  reference_URL='https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-63/plants/fasta/arabidopsis_thaliana/cds/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz',
  specification='Reconstruct only complete consistent ACGT query sequences; exact full-length matches on both strands against all reference CDS; strict sensitivity requires one exact locus agreeing with supplied mapping and no alternative supplied >=95%-identity/>=95%-query-coverage locus. No claim of experimentally proven hybridization specificity or isoform resolution.',
  command=commandArgs(),session=capture.output(sessionInfo()))
 dput(provenance,file=file.path(folder,'provenance.R'));writeLines('Complete; see provenance.R',marker)
}
writeLines('Complete; each artifact has provenance.R',file.path(out,'RUN_COMPLETE'))
message('Completed probe audit and sensitivity outputs')
