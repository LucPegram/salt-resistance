source_files <- function(root) {
  paths <- list.files(root, recursive = TRUE, all.files = TRUE, full.names = TRUE,
                      include.dirs = FALSE)
  relative <- substring(paths, nchar(root) + 2L)
  keep <- grepl("^(R|config|scripts|tests)/.*\\.R$", relative) |
    grepl("^(docs/.*\\.md|legacy/.*\\.(R|md)|README\\.md|DESCRIPTION|renv\\.lock|\\.gitignore)$", relative)
  paths[keep & !relative %in% "config/local.R"]
}

run_artifact <- function(id, project_root, data_root, output_root, cache = new.env()) {
  entry <- artifact_registry[[id]]
  if (is.null(entry)) stop("Unregistered output: ", id)
  destination_dir <- file.path(output_root, id)
  dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
  marker <- file.path(destination_dir, "RUN_COMPLETE")
  if (file.exists(marker)) unlink(marker)
  input_paths <- vapply(entry$inputs, function(n) file.path(data_root, input_specs[[n]]$path), "")
  hash <- function(path) digest::digest(file = path, algo = "sha256")
  code_paths <- source_files(project_root)
  code_before <- setNames(vapply(code_paths, hash, ""),
                          substring(code_paths, nchar(project_root) + 2L))
  if (any(!file.exists(input_paths))) stop("Missing required inputs for ", id)
  before <- vapply(input_paths, hash, "")
  inputs <- setNames(lapply(entry$inputs, function(n) {
    if (!exists(n, cache, inherits = FALSE)) {
      assign(n, list(value = read_input(n, data_root, input_specs), hash = before[[n]]), cache)
    }
    saved <- get(n, cache, inherits = FALSE)
    if (!identical(saved$hash, before[[n]])) stop("Input changed during this run: ", n)
    saved$value
  }), entry$inputs)
  set.seed(analysis_settings$seed)
  destination <- file.path(destination_dir, entry$filename)
  get(entry$builder, mode = "function")(inputs, analysis_settings, destination)
  output_paths <- if (is.null(entry$formats)) destination else
    paste0(sub("\\.pdf$", "", destination), ".", entry$formats)
  if (any(!file.exists(output_paths)) || any(file.info(output_paths)$size == 0)) {
    stop("Output was not successfully written: ", destination)
  }
  after <- vapply(input_paths, hash, "")
  if (!identical(before, after)) stop("An input changed while producing ", id)
  if (!identical(code_before, setNames(vapply(code_paths, hash, ""), names(code_before)))) {
    stop("Code changed while producing ", id)
  }
  revision <- NA_character_
  if (nzchar(Sys.which("git")) && file.exists(file.path(project_root, ".git"))) {
    candidate <- tryCatch(suppressWarnings(system2("git",
      c("-C", shQuote(project_root), "rev-parse", "--verify", "HEAD"),
      stdout = TRUE, stderr = FALSE)), error = function(e) character())
    if (length(candidate) == 1L && is.null(attr(candidate, "status")) &&
        grepl("^[0-9a-f]{40}$", candidate)) revision <- candidate
  }
  provenance <- list(artifact_id = id, status = entry$status,
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    command = commandArgs(), git_revision = revision,
    input_relative_paths = setNames(vapply(entry$inputs, function(n) input_specs[[n]]$path, ""), entry$inputs),
    input_sha256 = before, settings = analysis_settings,
    integration_settings = if (grepl("integration", id, fixed = TRUE)) integration_settings() else NULL,
    code_sha256 = code_before,
    output_sha256 = setNames(vapply(output_paths, hash, ""), basename(output_paths)),
    R_version = R.version.string,
    package_versions = setNames(vapply(c("readxl", "digest", "jsonlite", "limma", "statmod", "proDA", "GO.db", "AnnotationDbi"),
      function(n) as.character(utils::packageVersion(n)), ""), c("readxl", "digest", "jsonlite", "limma", "statmod", "proDA", "GO.db", "AnnotationDbi")))
  dput(provenance, file = file.path(destination_dir, "provenance.R"))
  writeLines(capture.output(sessionInfo()), file.path(destination_dir, "sessionInfo.txt"))
  writeLines("Complete; see provenance.R for inputs, code, and output hashes.", marker)
  message("Completed ", id, ": ", destination)
  invisible(destination)
}

