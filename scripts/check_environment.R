script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(script), "bootstrap.R"))
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install renv and restore renv.lock")
lock <- jsonlite::fromJSON(file.path(project_root, "renv.lock"), simplifyVector = FALSE)
problems <- character()
if (as.character(getRversion()) != lock$R$Version) {
  problems <- c(problems, paste("R", lock$R$Version, "required; found", getRversion()))
}
for (p in names(lock$Packages)) {
  expected <- lock$Packages[[p]]$Version
  actual <- if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "missing"
  if (actual=="missing" || package_version(expected)!=package_version(actual))
    problems <- c(problems, paste(p, expected, "required; found", actual))
}
if (length(problems)) stop(paste(problems, collapse = "\n"))
cat("R and all locked package versions match.\n")

