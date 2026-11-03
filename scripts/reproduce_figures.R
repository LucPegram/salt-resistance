# Rebuild four current main figures and the qPCR supplementary figure.
# The two integration figures are built by reproduce_integration.R.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
source(file.path(dirname(script), "bootstrap.R"))
if (length(commandArgs(TRUE))) stop("No arguments; set the external roots in the environment")
data_root <- Sys.getenv("SALT_DATA_ROOT")
output_root <- Sys.getenv("SALT_OUTPUT_ROOT")
if (!nzchar(data_root) || !nzchar(output_root)) stop("Set SALT_DATA_ROOT and SALT_OUTPUT_ROOT")
data_root <- normalizePath(data_root, mustWork=TRUE)
dir.create(output_root, recursive=TRUE, showWarnings=FALSE)
output_root <- normalizePath(output_root, mustWork=TRUE)
inside <- function(path,parent) identical(path,parent) || startsWith(path,paste0(parent,"/"))
if (inside(data_root,project_root) || inside(output_root,data_root) || inside(data_root,output_root) ||
    (inside(output_root,project_root) && !inside(output_root,file.path(project_root,"outputs"))))
  stop("Keep original data external and separate from generated outputs")
ids <- c("figure_revision_dose", "figure_revision_protein", "figure_revision_contrasts",
  "figure_revision_qc", "figure_manuscript_qpcr_candidates")
cache <- new.env(parent=emptyenv())
for (id in ids) run_artifact(id,project_root,data_root,output_root,cache)
