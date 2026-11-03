assert_columns <- function(x, required, context) {
  absent <- setdiff(required, names(x))
  if (length(absent)) stop(context, ": missing columns: ", paste(absent, collapse = ", "))
  if (anyDuplicated(names(x))) stop(context, ": duplicate column headings")
}

as_numeric_checked <- function(x, context) {
  text <- trimws(as.character(x))
  missing <- is.na(x) | text %in% c("", "NA", "NaN")
  value <- suppressWarnings(as.numeric(text))
  if (any(!missing & (is.na(value) | !is.finite(value)))) {
    stop(context, ": nonnumeric or nonfinite observations")
  }
  value[missing] <- NA_real_
  value
}

read_microarray <- function(path, spec) {
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                na.strings = c("", "NA"))
  assert_columns(d, "CORRECTED MAPPING", path)
  samples <- setdiff(names(d), "CORRECTED MAPPING")
  if (!length(samples)) stop("No expression columns: ", path)
  pattern <- if (spec$tissue == "Combined") {
    "^(BG|OE) (PN|100) (NaCl )?([12])$"
  } else {
    paste0("^(PN|50 NaCl|75 NaCl|100 NaCl) ", tolower(spec$tissue), " B[123]$")
  }
  if (any(!grepl(pattern, samples))) stop("Unexpected microarray sample headings: ", path)
  expected <- if (spec$tissue == "Combined") 8L else 12L
  if (length(samples) != expected) stop("Incomplete microarray sample design: ", path)
  for (s in samples) d[[s]] <- as_numeric_checked(d[[s]], paste(path, s))
  values <- as.matrix(d[samples])
  if (any(values <= 0, na.rm = TRUE)) stop("Microarray values must be positive before log2: ", path)
  ids <- trimws(d[["CORRECTED MAPPING"]])
  mapped <- !is.na(ids) & nzchar(ids) & ids != "#N/A"
  if (any(!grepl("^AT[1-5CM]G[0-9]+(\\.[0-9]+)?$", ids[mapped], ignore.case = TRUE))) {
    stop("Unrecognized mapped transcript identifiers: ", path)
  }
  list(data = d, values = values, ids = ids, mapped = mapped,
       samples = samples, spec = spec)
}

# Average repeated exact transcript IDs on the log scale. No annotation join.
# Transcript suffixes are retained: these are not automatically locus-level genes.
collapse_microarray <- function(x) {
  if (!any(x$mapped)) stop("No mapped microarray rows")
  values <- log2(x$values[x$mapped, , drop = FALSE])
  ids <- x$ids[x$mapped]
  present <- !is.na(values)
  values[!present] <- 0
  total <- rowsum(values, ids, reorder = FALSE)
  n <- rowsum(present * 1L, ids, reorder = FALSE)
  result <- total / n
  result[n == 0] <- NA_real_
  result
}

read_proteomics <- function(path, spec) {
  d <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                  na.strings = c("", "NA"), quote = "")
  samples <- unlist(lapply(c("No-0", "NIT2"), function(g) {
    unlist(lapply(c(0L, 75L), function(s) paste(g, paste0(s, "mM"), 1:3, sep = "_")))
  }))
  assert_columns(d, c("Protein.Group", "N.Sequences", "N.Proteotypic.Sequences", samples), path)
  quantitative <- grep("^(No-0|NIT2)_[0-9]+mM_[0-9]+$", names(d), value = TRUE)
  if (!setequal(samples, quantitative)) stop("Unexpected proteomics sample design: ", path)
  if (anyNA(d$Protein.Group) || any(!nzchar(d$Protein.Group)) || anyDuplicated(d$Protein.Group)) {
    stop("Protein-group identifiers must be nonempty and unique: ", path)
  }
  for (s in c("N.Sequences", "N.Proteotypic.Sequences", samples)) d[[s]] <- as_numeric_checked(d[[s]], paste(path, s))
  values <- as.matrix(d[samples])
  if (any(values <= 0, na.rm = TRUE)) stop("Nonpositive measured protein abundance: ", path)
  meta <- data.frame(sample = samples,
    genotype = sub("_.*", "", samples),
    salt_mM = as.integer(sub("^[^_]+_([0-9]+)mM_.*$", "\\1", samples)),
    replicate = as.integer(sub(".*_", "", samples)),
    tissue = spec$tissue, stringsAsFactors = FALSE)
  # Retain the complete protein group; do not choose its first accession as truth.
  list(data = d, values = values, samples = samples, metadata = meta, spec = spec)
}

proteomics_keep <- function(x, settings) {
  groups <- split(seq_len(ncol(x$values)),
                  paste(x$metadata$genotype, x$metadata$salt_mM, sep = "_"))
  eligible <- Reduce(`|`, lapply(groups, function(cols) {
    rowSums(!is.na(x$values[, cols, drop = FALSE])) >=
      settings$proteomics_min_observed_in_one_condition
  }))
  sequences <- x$data$N.Sequences
  eligible & !is.na(sequences) & sequences >= settings$proteomics_min_sequences
}

read_qpcr_summary <- function(path, spec) {
  if (!requireNamespace("readxl", quietly = TRUE)) stop("Restore the locked dependencies first")
  # Text import prevents automatic identifier/date or other type interpretation.
  d <- as.data.frame(readxl::read_excel(path, sheet = spec$sheet, col_types = "text",
                                      .name_repair = "minimal"), stringsAsFactors = FALSE)
  assert_columns(d, c("Gene", "Genotype", "Tissue", "Treatment",
                      "Average Normalized expression", "STE"), path)
  for (s in c("Gene", "Genotype", "Tissue", "Treatment")) {
    d[[s]] <- trimws(d[[s]])
    if (anyNA(d[[s]]) || any(!nzchar(d[[s]]))) stop("Missing qPCR group labels: ", s)
  }
  if (anyDuplicated(d[c("Gene", "Genotype", "Tissue", "Treatment")])) {
    stop("Duplicate qPCR group summaries")
  }
  for (s in c("Average Normalized expression", "STE")) {
    d[[s]] <- as_numeric_checked(d[[s]], s)
    if (any(d[[s]] < 0, na.rm = TRUE) ||
        (s == "Average Normalized expression" && anyNA(d[[s]]))) {
      stop("Invalid qPCR summary: ", s)
    }
  }
  list(data = d, spec = spec)
}

read_input <- function(id, data_root, specs) {
  spec <- specs[[id]]
  if (is.null(spec)) stop("Unknown input: ", id)
  path <- file.path(data_root, spec$path)
  if (!file.exists(path)) stop("Required external input not found: ", path)
  if (spec$assay == "go_associations") return(read_revision_associations(path, spec))
  if (startsWith(spec$assay, "legacy_")) return(read_legacy(path, spec))
  switch(spec$assay,
    microarray = read_microarray(path, spec),
    proteomics = read_proteomics(path, spec),
    qpcr_summary = read_qpcr_summary(path, spec),
    stop("Unknown assay: ", spec$assay))
}

