script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),'bootstrap.R'))
if(length(commandArgs(TRUE)))stop('No arguments; set SALT_PROBE_REFERENCE_ROOT')
root<-Sys.getenv('SALT_PROBE_REFERENCE_ROOT');if(!nzchar(root))stop('Set SALT_PROBE_REFERENCE_ROOT')
root<-normalizePath(root,mustWork=FALSE)
if(identical(root,project_root)||startsWith(root,paste0(project_root,'/')))stop('Reference data must stay outside code')
dir.create(root,recursive=TRUE,showWarnings=FALSE)
url<-'https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-63/plants/fasta/arabidopsis_thaliana/cds/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz'
path<-file.path(root,'TAIR10_release63_CDS.fa.gz')
expected<-'3580fd49dc079892359c1e5b7ffe76b7c2aa13718daa4daa6f397882ad41795d'
if(!file.exists(path))download.file(url,path,mode='wb',method='libcurl')
actual<-digest::digest(file=path,algo='sha256');if(actual!=expected)stop('Pinned reference checksum mismatch; inspect release changes before use')
write_table(data.frame(resource=basename(path),URL=url,release='Ensembl Plants 63; TAIR10 assembly',sha256=actual,
 verified_UTC=format(Sys.time(),tz='UTC',usetz=TRUE)),file.path(root,'download_manifest.csv'))
message('Pinned reference verified')
