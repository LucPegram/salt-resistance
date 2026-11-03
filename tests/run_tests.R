script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(script), "..", "scripts", "bootstrap.R"))
expect_error <- function(expr) {
  error <- tryCatch({force(expr); FALSE}, error = function(e) TRUE)
  stopifnot(error)
}
# Artificial values only; no experimental input or copied research result.
x <- list(ids = c("AT1G00010.1", "AT1G00010.1", "AT1G00010.2", "#N/A"),
  mapped = c(TRUE, TRUE, TRUE, FALSE),
  values = matrix(c(1, 4, 16, 100, 4, 16, NA, 100), ncol = 2))
collapsed <- collapse_microarray(x)
stopifnot(nrow(collapsed) == 2, collapsed["AT1G00010.1", 1] == 1,
          collapsed["AT1G00010.1", 2] == 3,
          collapsed["AT1G00010.2", 1] == 4,
          is.na(collapsed["AT1G00010.2", 2]))
expect_error(as_numeric_checked(c("2", "bad"), "test"))
expect_error(as_numeric_checked(c("2", "Inf"), "test"))
stopifnot(identical(as_numeric_checked(c("2", "", "NA"), "test"), c(2, NA_real_, NA_real_)))
protein <- list(values = rbind(c(1, 2, NA, NA, NA, NA),
  c(1, NA, NA, 2, NA, NA), c(1, 2, 3, 4, 5, 6)),
  data = data.frame(N.Sequences = c(2, 2, 1)),
  metadata = data.frame(genotype = rep("No-0", 6), salt_mM = rep(c(0, 75), each = 3)))
stopifnot(identical(proteomics_keep(protein, analysis_settings), c(TRUE, FALSE, FALSE)))
expect_error(assert_columns(data.frame(a = 1), c("a", "b"), "test"))
stopifnot(!anyDuplicated(names(artifact_registry)))
for (entry in artifact_registry) {
  stopifnot(all(entry$inputs %in% names(input_specs)), exists(entry$builder, mode = "function"))
}
mapped <- map_modified_peptide("AST(UniMod:21)K", "QQASTKQQ")
stopifnot(mapped$start == 3L, mapped$phosphosites == 5L)
expect_error(map_modified_peptide("AST(UniMod:21)K", "ASTKQQASTK"))
set.seed(11)
artificial_abundance <- matrix(exp(rnorm(36)), nrow = 3)
p <- balanced_protein_anova(artificial_abundance)
metadata <- data.frame(genotype = rep(c("A", "B"), each = 6),
  treatment = rep(rep(c("control", "salt"), each = 3), 2))
expected <- stats::anova(stats::lm(log(artificial_abundance[1, ]) ~ genotype * treatment,
  data = metadata))[["Pr(>F)"]][1:3]
stopifnot(max(abs(p[1, ] - expected)) < 1e-12)
cat("Artificial-example tests passed.\n")

# Statistical checks use synthetic measurements and independent lm calculations.
set.seed(17)
m<-data.frame(cell=factor(rep(c("WT0","WTS","OE0","OES"),each=3),
  levels=c("WT0","WTS","OE0","OES")))
design<-model.matrix(~0+cell,m);colnames(design)<-levels(m$cell)
y<-matrix(rnorm(1200),nrow=100,dimnames=list(paste0("feature",1:100),NULL))
y[1,]<-y[1,]+rep(c(0,1,2,4),each=3)
cm<-factorial_contrasts(design)
answer<-moderated_contrasts(y,design,cm)$table
expected<-coef(lm(y[1,]~0+cell,m))
z<-answer[answer$feature=="feature1"&answer$contrast=="interaction",]
stopifnot(abs(z$log2FC-((expected[4]-expected[3])-(expected[2]-expected[1])))<1e-12,
  all(answer$lower95<=answer$log2FC),all(answer$upper95>=answer$log2FC),
  all(answer$q_all_contrasts>=answer$q),
  all(answer$q_all_contrasts==pmin(1,7*answer$q)))
splits<-prestress_splits(y,m,"synthetic")
stopifnot(sum(splits$comparison=="disjoint_controls")==6)
for(label in splits$split[splits$comparison=="disjoint_controls"]){
  parts<-strsplit(label," vs ",fixed=TRUE)[[1]]
  a<-as.integer(strsplit(parts[1],",",fixed=TRUE)[[1]])
  b<-as.integer(strsplit(parts[2],",",fixed=TRUE)[[1]])
  stopifnot(length(intersect(a,b))==0,setequal(c(a,b),which(m$cell=="WT0")))
}
cat("Synthetic contrast, interval, multiplicity and disjoint-control checks passed.\n")

