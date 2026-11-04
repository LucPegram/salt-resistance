# Recover complete query sequences from the supplied alignment workbook. Partial
# alignments are never filled with guessed sequence. Exact CDS matching is an
# annotation diagnostic, not experimental proof of hybridization specificity.
bmc_probe_evidence <- function(workbook) {
 d<-as.data.frame(readxl::read_excel(workbook,sheet='probes-all_v_TAIR10_cds_2010121',col_types='text',.name_repair='minimal'))
 num<-function(x)suppressWarnings(as.numeric(sub(',','.',x,fixed=TRUE)))
 a<-data.frame(probe=trimws(d[[1]]),rank=num(d[[2]]),target=trimws(d[[3]]),
  identity=num(d[[9]]),alignment_length=num(d[[10]]),query_length=num(d[[14]]),
  query_start=num(d[[15]]),query_end=num(d[[16]]),query_sequence=toupper(gsub('-','',d[[17]],fixed=TRUE)),stringsAsFactors=FALSE)
 a$target_locus<-sub('\\.[0-9]+$','',a$target)
 valid<-!is.na(a$probe)&nzchar(a$probe)&grepl('^AT[1-5CM]G[0-9]{5}(\\.[0-9]+)?$',a$target)
 raw_rows<-nrow(a);raw_queries<-length(unique(a$probe[!is.na(a$probe)&nzchar(a$probe)]))
 excluded_rows<-sum(!valid)
 a<-a[valid,]
 complete<-pmin(a$query_start,a$query_end)==1 & pmax(a$query_start,a$query_end)==a$query_length &
  nchar(a$query_sequence)==a$query_length & grepl('^[ACGT]+$',a$query_sequence)
 complete[is.na(complete)]<-FALSE
 seqs<-a[complete,c('probe','query_sequence','query_start','query_end')]
 rev<-seqs$query_start>seqs$query_end
 if(any(rev))seqs$query_sequence[rev]<-as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seqs$query_sequence[rev])))
 pairs<-unique(seqs[c('probe','query_sequence')]);counts<-table(pairs$probe)
 conflicts<-names(counts[counts>1L]);pairs<-pairs[!pairs$probe%in%conflicts,]
 stopifnot(!anyDuplicated(pairs$probe))
 coverage<-abs(a$query_end-a$query_start)+1
 close<-a$identity>=95 & coverage/a$query_length>=.95
 close[is.na(close)]<-FALSE
 edges<-unique(a[close,c('probe','target_locus')])
 list(probes=pairs,near_edges=edges,alignments=a,
  audit=data.frame(quantity=c('input_alignment_rows','input_distinct_queries','excluded_invalid_alignment_rows','alignment_rows','distinct_queries','queries_with_complete_consistent_sequence','queries_with_conflicting_complete_sequence','near_match_edges'),
   value=c(raw_rows,raw_queries,excluded_rows,nrow(a),length(unique(a$probe)),nrow(pairs),length(conflicts),nrow(edges))))
}

bmc_exact_cds_map <- function(probes,reference) {
 cds<-Biostrings::readDNAStringSet(reference)
 transcript<-sub(' .*','',names(cds));locus<-sub('\\.[0-9]+$','',transcript)
 stopifnot(all(grepl('^AT[1-5CM]G[0-9]{5}$',locus)),!anyDuplicated(transcript))
 ends<-cumsum(Biostrings::width(cds));begins<-c(1,head(ends,-1)+1)
 subject<-unlist(cds,use.names=FALSE)
 seqs<-unique(probes$query_sequence)
 hits<-list()
 for(len in sort(unique(nchar(seqs)))){
  s<-seqs[nchar(seqs)==len]
  for(strand in c('+','-')){
   dna<-Biostrings::DNAStringSet(s);if(strand=='-')dna<-Biostrings::reverseComplement(dna)
   mi<-Biostrings::matchPDict(Biostrings::PDict(dna),subject,max.mismatch=0)
   h<-as.data.frame(mi)
   if(!nrow(h))next
   ti<-findInterval(h$start,begins)
   ok<-h$end<=ends[ti]
   hits[[length(hits)+1L]]<-data.frame(sequence=s[h$group[ok]],transcript=transcript[ti[ok]],locus=locus[ti[ok]],strand=strand,
    start_in_CDS=h$start[ok]-begins[ti[ok]]+1L,stringsAsFactors=FALSE)
  }
 }
 h<-unique(do.call(rbind,hits))
 list(hits=h,reference=data.frame(transcripts=length(cds),loci=length(unique(locus)),coding_bases=sum(Biostrings::width(cds))))
}

bmc_probe_audit <- function(evidence,matches,feature,corrected_ids) {
 stopifnot(nrow(feature)==length(corrected_ids))
 ids<-as.character(feature$ProbeName);loci<-sub('\\.[0-9]+$','',corrected_ids)
 seqs<-evidence$probes$query_sequence[match(ids,evidence$probes$probe)]
 exact<-lapply(split(matches$hits$locus,matches$hits$sequence),unique)
 near<-lapply(split(evidence$near_edges$target_locus,evidence$near_edges$probe),unique)
 targets<-lapply(seqs,function(s)if(is.na(s))character() else exact[[s]])
 near_targets<-lapply(ids,function(p)near[[p]])
 n<-lengths(targets);nn<-lengths(near_targets)
 contains<-vapply(seq_along(loci),function(i)loci[i]%in%targets[[i]],TRUE)
 near_agree<-vapply(seq_along(loci),function(i)!length(setdiff(near_targets[[i]],loci[i])),TRUE)
 status<-ifelse(is.na(seqs),'complete sequence unavailable',ifelse(n==0,'no exact CDS target',
  ifelse(n>1,'multiple exact CDS loci',ifelse(!contains,'exact CDS locus differs from supplied mapping','one exact matching CDS locus'))))
 strict<-!is.na(seqs)&n==1L&contains&near_agree
 data.frame(feature=feature$FeatureNum,probe=ids,supplied_transcript=corrected_ids,supplied_locus=loci,
  sequence_reconstructed=!is.na(seqs),exact_CDS_loci=n,exact_matches_supplied_locus=contains,
  exact_target_loci=vapply(targets,paste,'',collapse=';'),reported_near_match_loci=nn,
  reported_near_target_loci=vapply(near_targets,paste,'',collapse=';'),
  status=status,strict_probe_sensitivity_eligible=strict,stringsAsFactors=FALSE)
}

bmc_probe_sensitivity <- function(inputs,audit,external) {
 masked<-inputs
 for(id in c('microarray_root','microarray_shoot','microarray_genotype')){
  stopifnot(identical(masked[[id]]$ids,audit$supplied_transcript),length(masked[[id]]$mapped)==nrow(audit))
  masked[[id]]$mapped<-masked[[id]]$mapped & audit$strict_probe_sensitivity_eligible
 }
 dose<-extension_dose(masked);objects<-extension_objects(masked,dose)
 benchmark<-extension_benchmark(external,objects)
 benchmark$table$q_BY<-p.adjust(benchmark$table$P,'BY')
 primary_geno<-integration_array(inputs$microarray_genotype)
 strict_geno<-objects$Array2
 common<-intersect(rownames(primary_geno$y),rownames(strict_geno$y))
 tables<-list(probe_sensitive_dose_contrasts=dose$table,probe_sensitive_dose_joint=dose$joint,
  probe_sensitive_genotype_contrasts=strict_geno$table,probe_sensitive_external_tests=benchmark$table,
  probe_sensitive_external_coverage=benchmark$coverage,
  probe_sensitive_universes=do.call(rbind,lapply(objects,function(o)data.frame(dataset=o$dataset,loci=nrow(o$y),samples=ncol(o$y)))))
 comparisons<-do.call(rbind,lapply(unique(strict_geno$table$contrast),function(co){
  a<-primary_geno$table[primary_geno$table$contrast==co,];b<-strict_geno$table[strict_geno$table$contrast==co,]
  a<-a[match(common,a$feature),];b<-b[match(common,b$feature),]
  data.frame(contrast=co,shared_loci=length(common),Pearson=cor(a$log2FC,b$log2FC),Spearman=cor(a$log2FC,b$log2FC,method='spearman'),
   median_absolute_effect_change=median(abs(a$log2FC-b$log2FC)),maximum_absolute_effect_change=max(abs(a$log2FC-b$log2FC)),
   primary_BH_lt_05=sum(a$q<.05),strict_BH_lt_05=sum(b$q<.05),
   interpretation='Different mapping universes; discovery counts are not a comparison of biological effect sizes')
 }))
 tables$probe_sensitive_genotype_agreement<-comparisons
 tables
}
