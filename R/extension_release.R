extension_compute <- function(inputs,paths) {
  message("Fitting paired root-shoot dose model")
  dose<-extension_dose(inputs)
  message("Preparing complete-observation gene/protein universes")
  objects<-extension_objects(inputs,dose)
  message("Fitting independent public salt datasets")
  external<-extension_external(paths)
  message("Reading independent atlas marker sets")
  markers<-extension_read_markers(paths[["atlas2025_markers.xlsx"]],paths[["atlas2025_clusters.xlsx"]])
  message("Reading original-paper DAP-seq targets from AGRIS")
  dap<-extension_read_dap(paths[["AtRegNet.zip"]])
  message("Testing independent salt signatures")
  benchmark<-extension_benchmark(external,objects)
  message("Testing developmental/tissue signatures and focal-gene exclusion sensitivity")
  marker_result<-extension_markers(markers,objects)
  message("Testing TF target sets and inter-gene correlation sensitivity")
  tf<-extension_tf(dap,objects)
  for(n in c("benchmark","marker_result")){
    z<-get(n);z$table$q_BY<-if("mode"%in%names(z$table))
      ave(z$table$P,z$table$mode,FUN=function(p)p.adjust(p,"BY")) else p.adjust(z$table$P,"BY")
    assign(n,z)
  }
  tf$table$q_BY<-ave(tf$table$P,tf$table$mode,FUN=function(p)p.adjust(p,"BY"))
  list(dose=dose,objects=objects,external=external,markers=markers,dap=dap,
    benchmark=benchmark,marker_result=marker_result,tf=tf,annotation=toprank_annotation(inputs$gene_annotation))
}

extension_universe_table <- function(z) {
  do.call(rbind,lapply(z$objects,function(o)data.frame(dataset=o$dataset,tissue=o$tissue,assay=o$assay,
    loci=nrow(o$y),samples=ncol(o$y),missing_values=sum(!is.finite(o$y)),
    source_cell=paste(colnames(o$y),collapse="; "),
    eligibility=if(o$assay=="protein")"Unique single-accession AGI; proteotypic support; all 12 pools observed; no imputation" else
      "Mean of unique transcript log2 signals per locus; primary quantile-normalized project inputs")))
}

extension_dose_summary <- function(z) {
  d<-z$dose$table
  a<-do.call(rbind,lapply(split(d,d$contrast),function(x)data.frame(contrast=x$contrast[1],tested=nrow(x),
    global_q_lt_05=sum(x$q_global<.05),higher=sum(x$q_global<.05&x$log2FC>0),lower=sum(x$q_global<.05&x$log2FC<0))))
  a
}

extension_sample_table <- function(z) {
  do.call(rbind,lapply(z$objects,function(o){
    m<-o$metadata
    data.frame(dataset=o$dataset,sample=colnames(o$y),
      tissue=if("tissue"%in%names(m))m$tissue else o$tissue,
      genotype=if("genotype"%in%names(m))m$genotype else "Col-0",
      salt_mM=if("dose"%in%names(m))m$dose else m$salt_mM,
      pool=if("pool"%in%names(m))as.character(m$pool) else colnames(o$y))
  }))
}

extension_tables <- function(z) {
  named<-function(d){d$locus<-d$feature;toprank_add_names(d,z$annotation)}
  means<-z$dose$means
  dose_means<-data.frame(locus=rep(rownames(means),ncol(means)),cell=rep(colnames(means),each=nrow(means)),
    mean_log2_signal=as.vector(means))
  list(
    assay_universes=extension_universe_table(z),project_samples=extension_sample_table(z),
    dose_contrasts=named(z$dose$table),dose_joint_tests=toprank_add_names(z$dose$joint,z$annotation),
    dose_patterns=toprank_add_names(z$dose$patterns,z$annotation),dose_cell_means=toprank_add_names(dose_means,z$annotation),
    dose_discoveries=extension_dose_summary(z),
    external_samples=z$external$samples,external_probe_mapping=z$external$mapping,
    external_salt_contrasts=named(z$external$table),external_signature_members=z$benchmark$members,
    external_signature_coverage=z$benchmark$coverage,external_signature_scores=z$benchmark$scores,
    external_signature_tests=z$benchmark$table,external_effect_concordance=z$benchmark$correlations,
    TF_reference_audit=z$dap$audit,TF_target_coverage=z$tf$coverage,TF_target_enrichment=z$tf$table,
    TF_focal_gene_targets=merge(z$dap$edges[z$dap$edges$target%in%c(unname(integration_settings()$candidate_loci),
      "AT3G48750","AT5G44720","AT1G18390"),],z$annotation,by.x="target",by.y="locus",all.x=TRUE),
    marker_members=z$markers,marker_coverage=z$marker_result$coverage,marker_scores=z$marker_result$scores,
    marker_tests=z$marker_result$table)
}

extension_report_pages <- function(z) {
  fmt<-function(x)formatC(x,digits=3,format="g")
  d<-z$dose$patterns;j<-z$dose$joint
  counts<-vapply(c("Root","Shoot","Difference"),function(t)sum(j$test==t&j$q_global<.05),0L)
  b<-z$benchmark$table;t<-z$tf$table;m<-z$marker_result$table
  basal<-b[b$contrast=="baseline_genotype",];basal<-basal[order(basal$q_global),]
  sig_basal<-basal[basal$q_global<.05,]
  topbasal<-if(nrow(sig_basal))vapply(seq_len(min(6,nrow(sig_basal))),function(i){x<-sig_basal[i,]
    nit2<-any(z$benchmark$members$reference==x$reference & z$benchmark$members$locus=="AT3G44300")
    paste(x$dataset,x$reference,": score change",fmt(x$score_difference),"; global BH q",fmt(x$q_global),
      "; BY q",fmt(x$q_BY),"; NIT2 in reference signature:",if(nit2)"yes" else "no")},"") else
    "No baseline-genotype salt-signature score passed global BH q < 0.05."
  tfp<-t[t$mode=="primary",];tfp<-tfp[order(tfp$q_global,tfp$TF),]
  top_tf<-vapply(seq_len(min(5,nrow(tfp))),function(i){x<-tfp[i,]
    ss<-t[t$mode=="sensitivity"&t$dataset==x$dataset&t$contrast==x$contrast&t$TF==x$TF,]
    paste(x$TF_name,x$dataset,x$contrast,x$Direction,": global q",fmt(x$q_global),
      "; q with correlation 0.05",fmt(ss$q_global))},"")
  mg<-m[m$mode=="primary"&m$dataset=="Array2"&m$contrast%in%c("baseline_genotype","interaction"),]
  mg<-mg[order(mg$q_global,mg$set),]
  marker_lines<-vapply(seq_len(min(7,nrow(mg))),function(i){x<-mg[i,]
    ss<-m[m$mode=="without_focal_candidates"&m$dataset==x$dataset&m$contrast==x$contrast&m$set==x$set,]
    paste(x$set,":",x$contrast,"change",fmt(x$score_difference),"; BH q",fmt(x$q_global),"; BY q",fmt(x$q_BY),
      "; focal-gene-excluded q",if(nrow(ss))fmt(ss$q_global) else "not eligible")},"")
  list(
    `Design and dose-response findings`=c(
      "Computational extension | Arabidopsis NIT2 and salt response | 4 November 2026",
      "These analyses extend the primary quantile-normalized release. The requested 10x dataset was excluded. External resources, output tables, and figures remain outside the code-only repository.",
      paste("Array 1:",format(nrow(z$objects$Array1$y),big.mark=","),"loci, roots and shoots, 0/50/75/100 mM, three independent harvested pools per dose; root-shoot pairing retained. Array 2: No-0/NIT2, 0/100 mM, two pools per cell. Protein signature/TF tests use unambiguous, fully observed mapped loci only:",
        format(nrow(z$objects$Protein_Root$y),big.mark=","),"root and",format(nrow(z$objects$Protein_Shoot$y),big.mark=","),"shoot loci, 12 pools each."),
      "Dose model: categorical tissue x dose, paired-pool consensus correlation, limma trend/robust moderation. Fifteen planned contrasts cover control comparisons, adjacent increments, and direct differences between tissue responses. BH is calculated globally over all locus-contrast tests; the three joint tests have a separate global family.",
      paste("Joint global BH q < 0.05:",paste(names(counts),counts,collapse="; "),"loci. These locus-level results are distinct from earlier transcript-ID summaries."),
      paste("Among dose-responsive root loci:",paste(names(table(d$mean_pattern[d$tissue=="Root"&d$joint_q_global<.05])),
        as.integer(table(d$mean_pattern[d$tissue=="Root"&d$joint_q_global<.05])),collapse="; ")),
      paste("Among dose-responsive shoot loci:",paste(names(table(d$mean_pattern[d$tissue=="Shoot"&d$joint_q_global<.05])),
        as.integer(table(d$mean_pattern[d$tissue=="Shoot"&d$joint_q_global<.05])),collapse="; ")),
      "Mean-shape labels are descriptive. A first significant tested concentration is not a physiological threshold. Separate quantile normalization by tissue prevents absolute transcript-abundance or total-RNA comparisons between organs. A root-only significant test is not evidence of a tissue difference; use the direct Difference contrasts."
    ),
    `Independent salt-response benchmarking`=c(
      "GSE5623 salt + GSE5620 matched-time controls: Col-0, 16-day plants, roots/shoots, 150 mM, 0.5/1/3/6/12/24 h, two arrays per cell. Deposited positive MAS5 values were log2 transformed and quantile normalized together within each tissue; uniquely mapped probes were collapsed to AGI loci.",
      "GSE223207: Col-0 hydroponic roots, age 16 days, 25/50/75/100 mM, 15 minutes, two arrays per cell. Deposited GCRMA log2 values were used without another log transform. All 16 deposited treatment groups remain distinct when estimating residual variance; only four NaCl-control contrasts define salt signatures.",
      "Each public signature requires external within-contrast BH q < 0.05 and |log2FC| >= log2(1.5). Up/down sets must each have at least ten measured members. A score is mean within-sample fractional rank of Up genes minus that of Down genes. Differential score tests use biological samples, paired root-shoot structure where applicable, and a global BH/BY family. No current-project outcome selected signature genes.",
      paste(nrow(b),"eligible signature-contrast tests;",sum(b$q_global<.05),"global BH discoveries;",sum(b$q_BY<.05),"global BY discoveries."),
      paste("Eligible thresholded reference signatures:",length(unique(b$reference[startsWith(b$reference,"AtGen")])),
        "AtGenExpress time/tissue contexts and",length(unique(b$reference[startsWith(b$reference,"Hydro")])),
        "15-minute hydroponic contexts. All hydroponic gene-effect correlations are retained, even when no signature meets the significance/effect/coverage criteria. No cutoff was relaxed to obtain signatures."),
      "Baseline NIT2 minus No-0 results with global BH q < 0.05:",topbasal,
      "Effect-vector correlations and direction agreement are descriptive across shared loci, with no gene-count-based significance tests. A positive baseline signature is compatible with salt-response resemblance; it cannot establish protection, stress memory, or equivalence to a salt-treated plant.",
      "These references capture early responses and differ in dose, age, medium, platform, and genotype from the current experiments. They are independent contextual benchmarks, not matched replication of the two-week response. Cackett et al. 2022 material was not used as an independent reference."
    ),
    `Transcription-factor target evidence`=c(
      paste("AGRIS download filtered to the original O'Malley et al. 2016 DAP-seq paper (PMID 27203113), explicit DAP provenance, direct target type, valid AGI pairs:",format(nrow(z$dap$edges),big.mark=","),"unique edges from",length(unique(z$dap$edges$TF)),"TFs.",z$dap$audit$excluded_rows,"invalid candidate rows were excluded and audited."),
      "TF target tests: competitive rank-based CAMERA on moderated gene t statistics; 10 members minimum, at most half the measured universe. Each inter-gene correlation assumption (0.01 primary; 0.05 sensitivity) has a global BH/BY family across all TFs, assays, and selected contrasts. Fixed correlation assumptions are not estimates of TF activation.",
      paste("TF global BH discoveries:",sum(t$mode=="primary"&t$q_global<.05),"at 0.01 correlation;",sum(t$mode=="sensitivity"&t$q_global<.05),"at 0.05 correlation."),top_tf,
      "The two correlation settings represent sensitivity assumptions, not two opportunities to choose a favorable finding. A factor with q < 0.05 only under the primary setting remains a correlation-sensitive candidate. Targets may overlap extensively among TF family members. A DAP-seq binding target is an in-vitro DNA-binding association, not proof of binding, activation, or repression in the current tissues.",
      "The full target-enrichment table contains every tested TF and contrast under both settings. The focal-target table records which reference TFs have binding associations with the transporters and cross-assay candidate loci; these individual associations are reference annotations, not new experimental validation."
    ),
    `Developmental and tissue-marker signatures`=c(
      "Developmental/cell identity sets were obtained from Lee et al. 2025 Tables S1-S2, including Shahan et al. 2022 root markers, three broad root stages, literature proliferation markers, and top 50 eligible markers per annotated 21-day rosette cluster. Tissue-specific tests are restricted to anatomically applicable signatures; combined-tissue Array 2 allows both organ signatures.",
      "Marker scores are mean within-sample fractional expression ranks (minimum five measured genes). They measure signature expression, not cell abundance. Primary global BH/BY and a separate sensitivity excluding focal candidates are supplied. No marker score was regressed out of the primary expression models.",
      paste("Array 2 baseline marker discoveries:",sum(m$mode=="primary"&m$dataset=="Array2"&m$contrast=="baseline_genotype"&m$q_global<.05),
        "under global BH;",sum(m$mode=="primary"&m$dataset=="Array2"&m$contrast=="baseline_genotype"&m$q_BY<.05),
        "under global BY. These changes motivate developmental alternatives but do not establish altered cell proportions."),
      "Strongest Array 2 basal/interaction marker tests (ranked regardless of significance):",marker_lines
    ),
    `Interpretation, reproducibility and sources`=c(
      "Interpretation: these analyses can identify organ-dependent salt programs, resemblance to independent early stress, candidate regulator targets, and developmental alternatives to the pre-stress hypothesis. Neither a marker-score shift nor DAP target enrichment resolves cell composition versus within-cell regulation. They do not establish ion flux, ROS concentration, or transporter activity.",
      "Protein signatures use complete observations for a stable per-sample gene universe; they are supplementary sensitivity analyses rather than replacements for the primary available-observation protein contrasts. Combined-tissue RNA cannot identify which organ produced a genotype effect. No cross-assay sample pairing or inferred batch blocking was introduced.",
      "Every CSV and figure can be reproduced individually with scripts/reproduce_extensions.R <artifact_id>, or together with 'all'. Run fetch_extension_references.R first. Set SALT_DATA_ROOT to the quantile input root, SALT_REFERENCE_ROOT to external references, SALT_MODEL_CACHE to the existing external cache, and SALT_OUTPUT_ROOT to a new external output directory. Output provenance records original-input, normalization, reference, code, and output hashes, settings, versions, and Git revision.",
      "Sources and primary methods:",
      "AtGenExpress: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5623 ; matched controls GSE5620. Kilian et al. 2007, PMID 17376166.",
      "Early hydroponic dose series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223207 ; associated publication PMID 37603421.",
      "TF binding: O'Malley et al. 2016, Cell, https://doi.org/10.1016/j.cell.2016.04.038 ; AGRIS https://agris-knowledgebase.org/Downloads/AtRegNet.zip . PlantRegMap downloads were inaccessible, so no PlantRegMap prediction file was substituted or claimed.",
      "Atlas: Lee et al. 2025, A single-cell, spatial transcriptomic atlas of the Arabidopsis life cycle. https://doi.org/10.1038/s41477-025-02072-z . Only published marker tables were used; no 10x spatial dataset was used.",
      "Statistics: limma User's Guide, https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf . Hypothesis families and score definitions are detailed in docs/computational_extensions.md."
    ))
}

extension_report <- function(z,path) {
  pages<-extension_report_pages(z)
  grDevices::cairo_pdf(path,width=8.3,height=11.7,pointsize=10,family="sans")
  on.exit(dev.off())
  for(i in seq_along(pages)){
    par(mar=c(.8,.8,.8,.8));plot.new();par(usr=c(0,1,0,1))
    text(0,.985,names(pages)[i],adj=c(0,1),font=2,cex=1.5,col="#253342")
    y<-.94
    for(p in pages[[i]]){
      lines<-strwrap(p,width=103)
      for(line in lines){text(0,y,line,adj=c(0,1),cex=.84,col="#253342");y<-y-.017}
      y<-y-.012
    }
    if(y<.035)stop("Report page overflows; revise pagination")
    text(0,.01,"NIT2 salt-resistance | computational extension | supplementary analyses",adj=0,cex=.7,col="#5D6975")
    text(1,.01,paste(i,"/",length(pages)),adj=1,cex=.7,col="#5D6975")
  }
}

extension_registry <- function() c(
  setNames(rep("table",23),c("assay_universes","project_samples","dose_contrasts","dose_joint_tests","dose_patterns",
    "dose_cell_means","dose_discoveries","external_samples","external_probe_mapping","external_salt_contrasts",
    "external_signature_members","external_signature_coverage","external_signature_scores","external_signature_tests",
    "external_effect_concordance","TF_reference_audit","TF_target_coverage","TF_target_enrichment","TF_focal_gene_targets",
    "marker_members","marker_coverage","marker_scores","marker_tests")),
  figure_dose_programs="dose",figure_dose_candidates="candidates",figure_external_signatures="benchmark",
  figure_external_concordance="concordance",figure_TF_targets="TF",figure_marker_signatures="markers",
  analysis_report="report")
