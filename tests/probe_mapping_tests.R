script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
source(file.path(dirname(script),'..','scripts','bootstrap.R'))
# Independent synthetic targets include an off-target duplicate, a reverse
# complement match, and a sequence that would match only across a CDS boundary.
ref<-tempfile(fileext='.fa')
seqs<-Biostrings::DNAStringSet(c(AT1G00010.1='AACCGGTTACGT',AT1G00020.1='TTCCGGTTGGAA',AT1G00030.1='TATAGGGGTTAT'))
Biostrings::writeXStringSet(seqs,ref)
p<-data.frame(probe=c('unique','shared','reverse','boundary'),query_sequence=c('ACCGGT','CCGGTT','ATAACC','ACGTTT'))
z<-bmc_exact_cds_map(p,ref)
stopifnot(setequal(unique(z$hits$locus[z$hits$sequence=='CCGGTT']),c('AT1G00010','AT1G00020')),
 'AT1G00030'%in%z$hits$locus[z$hits$sequence=='ATAACC'],!any(z$hits$sequence=='ACGTTT'))
e<-list(probes=p,near_edges=data.frame(probe=c('unique','shared'),target_locus=c('AT1G00010','AT1G00020')))
f<-data.frame(ProbeName=p$probe,FeatureNum=1:4)
a<-bmc_probe_audit(e,z,f,c('AT1G00010.1','AT1G00010.1','AT1G00030.1','AT1G00010.1'))
stopifnot(a$exact_CDS_loci[2]==2L,a$strict_probe_sensitivity_eligible[3],!a$strict_probe_sensitivity_eligible[4])
unlink(ref);cat('Independent matching, reverse-strand, multi-locus and boundary checks passed\n')
