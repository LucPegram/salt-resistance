script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(script), "bootstrap.R"))
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("Provide one output ID, 'all', or '--list'")
if (args == "--list") {
  for (id in names(artifact_registry)) {
    cat(id, "[", artifact_registry[[id]]$status, "]\n", sep = " ")
  }
} else {
  for (p in c("readxl", "digest", "jsonlite")) {
    if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p, ". Restore renv.lock.")
  }
  data_root <- Sys.getenv("SALT_DATA_ROOT")
  if (!nzchar(data_root)) stop("Set SALT_DATA_ROOT to the external extracted-input directory")
  data_root <- normalizePath(data_root, mustWork = TRUE)
  inside <- function(path, parent) identical(path, parent) || startsWith(path, paste0(parent, "/"))
  if (inside(data_root, project_root)) stop("Keep experimental inputs outside this repository")
  output_root <- Sys.getenv("SALT_OUTPUT_ROOT", file.path(project_root, "outputs"))
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  output_root <- normalizePath(output_root, mustWork = TRUE)
  if (inside(output_root, data_root) || inside(data_root, output_root)) {
    stop("Input and output directories must not contain one another")
  }
  if (inside(output_root, project_root) &&
      !inside(output_root, file.path(project_root, "outputs"))) {
    stop("Only the ignored outputs/ directory may hold generated files inside the repository")
  }
  ids <- if (args == "all") names(artifact_registry) else args
  if (any(!ids %in% names(artifact_registry))) stop("Unregistered output ID: ", args)
  cache <- new.env(parent = emptyenv())
  for (id in ids) run_artifact(id, project_root, data_root, output_root, cache)
}

