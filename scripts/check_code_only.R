script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(script), "bootstrap.R"))
args <- commandArgs(TRUE)
if (length(args) > 1L) stop("Optional argument: a file containing repository-relative paths to check")
if (length(args)) {
  paths <- readLines(args, warn = FALSE)
} else if (file.exists(file.path(project_root, ".git"))) {
  paths <- system2("git", c("-C", shQuote(project_root), "ls-files"), stdout = TRUE)
} else {
  paths <- substring(list.files(project_root, recursive = TRUE, all.files = TRUE,
                                full.names = TRUE), nchar(project_root) + 2L)
}
allowed <- grepl("^(R|config|scripts|tests)/[A-Za-z0-9_./-]+\\.R$", paths) |
  grepl("^(README\\.md|docs/[A-Za-z0-9_./-]+\\.md|DESCRIPTION|renv\\.lock|\\.gitignore)$", paths)
forbidden <- grepl("(^|/)(data|source_data|outputs|results|figures|tables|\\.git|\\.Rproj.user)(/|$)", paths) |
  paths == "config/local.R" | grepl("(^|/)\\.\\.(/|$)", paths)
bad <- paths[!allowed | forbidden]
if (length(bad)) stop("Files prohibited in this code-only repository:\n", paste(bad, collapse = "\n"))
cat("Code-only check passed for", length(paths), "files.\n")

