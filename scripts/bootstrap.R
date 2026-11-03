script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("Use Rscript --vanilla scripts/<script>.R")
project_root <- normalizePath(file.path(dirname(sub("^--file=", "", script_argument)), ".."),
                              mustWork = TRUE)
Sys.setenv(TZ = "UTC")
if (nzchar(Sys.getenv("SALT_R_LIB"))) .libPaths(c(Sys.getenv("SALT_R_LIB"), .libPaths()))
invisible(Sys.setlocale("LC_COLLATE", "C"))
for (p in sort(list.files(file.path(project_root, "R"), pattern = "\\.R$", full.names = TRUE))) {
  source(p, local = .GlobalEnv)
}
source(file.path(project_root, "config", "inputs.R"))
source(file.path(project_root, "config", "analysis.R"))
source(file.path(project_root, "config", "artifacts.R"))

