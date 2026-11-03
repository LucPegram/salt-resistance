# Maintainer command: update the lock only after deliberately changing dependencies.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(script), "bootstrap.R"))
direct <- c("readxl", "digest", "jsonlite", "limma", "statmod", "proDA", "GO.db", "AnnotationDbi")
installed <- utils::installed.packages()
dependencies <- tools::package_dependencies(direct, db = installed,
  which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE)
packages <- sort(unique(c(direct, unlist(dependencies))))
packages <- packages[packages != "R" & packages %in% rownames(installed)]
packages <- packages[is.na(installed[packages, "Priority"])]
records <- setNames(lapply(packages, function(p) {
  description <- utils::packageDescription(p)
  repository<-description$Repository
  if(!is.null(repository)&&repository=="CRAN")
    return(list(Package=p,Version=description$Version,Source="Repository",Repository="CRAN"))
  if(!is.null(description$biocViews)||(!is.null(repository)&&grepl("Bioconductor",repository)))
    return(list(Package=p,Version=description$Version,Source="Bioconductor",Repository="Bioconductor"))
  stop("Package installation source needs explicit handling: ", p)
}), packages)
lock <- list(Bioconductor=list(Version="3.22"),R = list(Version = as.character(getRversion()),
  Repositories = list(list(Name = "CRAN", URL = "https://cloud.r-project.org"))), Packages = records)
jsonlite::write_json(lock, file.path(project_root, "renv.lock"), pretty = TRUE, auto_unbox = TRUE)
cat("Recorded installed dependency versions; restoration must be verified separately.\n")

