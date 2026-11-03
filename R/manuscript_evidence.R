# Reproduce summaries of supplied analysis outputs; these functions do not refit
# or validate the legacy models, nor convert their statistics into new contrasts.
legacy_result_specs <- function() {
  base <- "Microarray_for_Stephano/Luc_Initial_Analysis/"
  specs <- list()
  for (tissue in c("Root", "Shoot")) {
    prefix <- paste0(base, "Col-0_SaltHydroponics/", tissue, "/R_output/")
    specs[[paste0("saved_dose_linear_", tolower(tissue))]] <- list(
      path = paste0(prefix, "linear_trend_results_all.csv"), assay = "legacy_csv",
      tissue = tissue, model = "Linear dose trend")
    specs[[paste0("saved_dose_quadratic_", tolower(tissue))]] <- list(
      path = paste0(prefix, "poly2/poly2_results_all.csv"), assay = "legacy_csv",
      tissue = tissue, model = "Quadratic joint dose test")
    specs[[paste0("saved_proteomics_", tolower(tissue))]] <- list(
      path = paste0("Microarray_for_Stephano/Proteomics_Files/", tissue,
                    "Tissue_summary_NotImputed.xlsx"), assay = "legacy_xlsx",
      tissue = tissue, sheet = 1L)
  }
  prefix <- paste0(base, "NIT2vsNo-0_Treatments/R_analysis/")
  specs$saved_genotype_anova <- list(path = paste0(prefix, "Two-wayAnova/two_way_anova_results.csv"),
    assay = "legacy_csv", tissue = "Combined", model = "Two factor ANOVA")
  for (method in c("ordinary_t-test", "limma_eBayes_t-test")) {
    id <- if (method == "ordinary_t-test") "saved_genotype_ordinary" else "saved_genotype_moderated"
    specs[[id]] <- list(path = paste0(prefix, method, "/genotype_comparison_all.csv"),
      assay = "legacy_csv", tissue = "Combined", model = method)
  }
  for (comparison in c("OE NaCl v UT", "OE v WT NaCl", "OE v WT UT", "WT NaCl v UT")) {
    id <- paste0("saved_pairwise_", gsub(" ", "_", comparison))
    specs[[id]] <- list(path = paste0("Microarray_for_Stephano/",
      "Stu Analysis_Pairwise DEG Comparisons/DEG - ", comparison, ".xlsx"),
      assay = "legacy_pairwise", comparison = comparison)
  }
  prefix <- paste0("Data_from_Luc_2/Proteomics_AgriGO_Analysis/Root/Phospho/",
    "PhosphoPeptides_vs_DataSet/NA_4_or_less/")
  for (n in c("Response to Lipid", "Transport", "Transporter Activity_agi")) {
    specs[[paste0("saved_phospho_", gsub(" ", "_", n))]] <- list(
      path = paste0(prefix, n, ".csv"), assay = "legacy_csv", model = "Selected precursor subset")
  }
  specs
}

read_legacy <- function(path, spec) {
  if (spec$assay == "legacy_fasta") {
    lines <- readLines(path, warn = FALSE)
    if (sum(startsWith(lines, ">")) != 1L) stop("Expected one reference protein")
    sequence <- paste(lines[!startsWith(lines, ">")], collapse = "")
    if (!grepl("^[A-Z]+$", sequence)) stop("Invalid reference amino acid sequence")
    return(list(data = data.frame(accession = spec$accession, sequence = sequence), spec = spec))
  }
  if (spec$assay == "legacy_pairwise") {
    read_sheet <- function(sheet, skip = 0L, names = TRUE) {
      as.data.frame(readxl::read_excel(path, sheet = sheet, skip = skip,
        col_names = names, col_types = "text", .name_repair = "minimal"))
    }
    return(list(data = read_sheet(1L, names = FALSE), up = read_sheet("DEGs UP", 1L),
      down = read_sheet("DEGs DOWN", 1L), spec = spec))
  }
  d <- if (spec$assay == "legacy_xlsx") {
    as.data.frame(readxl::read_excel(path, sheet = spec$sheet, col_types = "text",
                                    .name_repair = "minimal"))
  } else if (spec$assay == "legacy_go") {
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "",
               comment.char = "", na.strings = c("", "NA"))
  } else read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                  na.strings = c("", "NA"))
  list(data = d, spec = spec)
}

checked_q <- function(x, column) {
  if (!column %in% names(x)) stop("Missing FDR column: ", column)
  q <- as_numeric_checked(x[[column]], column)
  if (any(q < 0 | q > 1, na.rm = TRUE)) stop("FDR outside [0,1]")
  q
}

build_existing_microarray_results <- function(inputs, settings, destination) {
  rows <- list()
  for (id in names(inputs)) {
    x <- inputs[[id]]; d <- x$data
    assert_columns(d, "AGI MAPPING", id)
    if (anyDuplicated(d[["AGI MAPPING"]])) stop("Duplicated result identifiers in ", id)
    terms <- if (id == "saved_genotype_anova") {
      c(Genotype = "adj.P_genotype", Salt = "adj.P_treatment", Interaction = "adj.P_interaction")
    } else if (grepl("quadratic", id)) {
      c(Joint = "adj.P_joint", Curvature = "adj.P_quadratic")
    } else c(Effect = "adj.P.Val")
    for (term in names(terms)) {
      q <- checked_q(d, terms[[term]])
      effect <- if ("slope_per_100" %in% names(d)) d$slope_per_100 else
        if ("log2FC" %in% names(d)) d$log2FC else rep(NA_real_, nrow(d))
      sig <- !is.na(q) & q < 0.05
      rows[[length(rows) + 1L]] <- data.frame(input = id, tissue = x$spec$tissue,
        model = x$spec$model, term = term, tested_ids = sum(!is.na(q)),
        FDR_lt_0_05 = sum(sig), positive = if (all(is.na(effect))) NA_integer_ else sum(sig & effect > 0),
        negative = if (all(is.na(effect))) NA_integer_ else sum(sig & effect < 0),
        source = x$spec$path)
    }
  }
  write_table(do.call(rbind, rows), destination)
}

build_manuscript_design <- function(inputs, settings, destination) {
  rows <- lapply(names(inputs), function(id) {
    x <- inputs[[id]]
    qpcr <- x$spec$assay == "qpcr_summary"
    combined <- identical(x$spec$tissue, "Combined")
    data.frame(dataset = id,
      genotype = if (qpcr || x$spec$assay == "proteomics" || combined) "No-0 and 35S::NIT2" else "Col-0",
      salt_mM = if (qpcr || x$spec$assay == "proteomics") "0 and 75" else if (combined) "0 and 100" else "0, 50, 75 and 100",
      tissue = if (qpcr) "Root and Shoot" else x$spec$tissue,
      biological_pools_per_cell = if (qpcr) NA_integer_ else if (combined) 2L else 3L,
      measurements = if (qpcr) "Group summaries; biological n unreported" else as.character(length(x$samples)))
  })
  write_table(do.call(rbind, rows), destination)
}

build_manuscript_qpcr_candidates <- function(inputs, settings, destination) {
  write_table(qpcr_summary_table(inputs$qpcr_summary$data, settings,
    c("AHA2", "SOS1", "HKT1", "HAK5")), destination)
}

build_manuscript_qpcr_full <- function(inputs, settings, destination) {
  d <- inputs$qpcr_summary$data
  write_table(qpcr_summary_table(d, settings, sort(unique(d$Gene))), destination)
}

qpcr_summary_table <- function(d, settings, genes) {
  rows <- list()
  for (gene in genes) for (tissue in c("Root", "Shoot")) {
    row <- data.frame(Gene = gene, Tissue = tissue)
    for (genotype in settings$qpcr_genotype_labels) for (treatment in settings$qpcr_treatment_labels) {
      z <- d[d$Gene == gene & d$Tissue == tissue & d$Genotype == genotype & d$Treatment == treatment, ]
      if (nrow(z) != 1L) stop("Missing or repeated candidate qPCR group")
      label <- paste(if (genotype == "No-0") "WT" else "NIT2",
        if (treatment == "Untreated") "0" else "75", sep = "_")
      row[[label]] <- paste0(format(signif(z[["Average Normalized expression"]], 3), trim = TRUE),
        " (", if (is.na(z$STE)) "not reported" else format(signif(z$STE, 3), trim = TRUE), ")")
    }
    rows[[length(rows) + 1L]] <- row
  }
  do.call(rbind, rows)
}

build_manuscript_qpcr_figure <- function(inputs, settings, destination) {
  settings$qpcr_genes <- c("AHA2", "SOS1", "HKT1", "HAK5")
  settings$qpcr_panel_columns <- 2L
  build_qpcr_figure(inputs, settings, destination)
}

build_existing_pairwise_results <- function(inputs, settings, destination) {
  rows <- lapply(inputs, function(x) {
    for (d in list(x$up, x$down)) assert_columns(d, c("ID_correct", "adj.P.Val"), x$spec$path)
    up <- x$up[!is.na(x$up$ID_correct), , drop = FALSE]
    down <- x$down[!is.na(x$down$ID_correct), , drop = FALSE]
    q <- c(checked_q(up, "adj.P.Val"), checked_q(down, "adj.P.Val"))
    if (anyNA(q) || any(q >= 0.05)) stop("Legacy DEG list contains entries without FDR < 0.05")
    total <- suppressWarnings(sum(as.numeric(x$data[-1, 2]), na.rm = TRUE))
    data.frame(comparison = x$spec$comparison, summary_entries = total,
      up_entries = nrow(up), down_entries = nrow(down),
      distinct_mapped_ids = length(unique(c(up$ID_correct, down$ID_correct))), source = x$spec$path)
  })
  write_table(do.call(rbind, rows), destination)
}

build_existing_proteomics_results <- function(inputs, settings, destination) {
  rows <- lapply(inputs, function(x) {
    d <- x$data
    qg <- checked_q(d, "p_FDR_genotype"); qt <- checked_q(d, "p_FDR_treatment")
    qi <- checked_q(d, "p_FDR_interaction")
    data.frame(tissue = x$spec$tissue, protein_groups = nrow(d),
      genotype_FDR_lt_0_05 = sum(qg < 0.05, na.rm = TRUE),
      salt_FDR_lt_0_05 = sum(qt < 0.05, na.rm = TRUE),
      interaction_FDR_lt_0_05 = sum(qi < 0.05, na.rm = TRUE),
      any_FDR_lt_0_05 = sum(qg < 0.05 | qt < 0.05 | qi < 0.05, na.rm = TRUE),
      genotype_and_salt = sum(qg < 0.05 & qt < 0.05, na.rm = TRUE),
      source = x$spec$path)
  })
  write_table(do.call(rbind, rows), destination)
}

# Ordinary balanced 2x2 ANOVA on fully observed groups only. This audit neither
# imputes missing values nor validates the legacy imputation-dependent tests.
balanced_protein_anova <- function(values) {
  stopifnot(ncol(values) == 12L, all(is.finite(values)), all(values > 0))
  y <- log(values)
  cells <- vapply(seq(1L, 12L, 3L), function(i) rowMeans(y[, i:(i + 2L), drop = FALSE]),
                  numeric(nrow(y)))
  cells <- matrix(cells, nrow = nrow(y), ncol = 4L)
  overall <- rowMeans(y)
  genotype <- cbind(rowMeans(y[, 1:6, drop = FALSE]), rowMeans(y[, 7:12, drop = FALSE]))
  treatment <- cbind(rowMeans(y[, c(1:3, 7:9), drop = FALSE]), rowMeans(y[, c(4:6, 10:12), drop = FALSE]))
  ss_genotype <- 6 * rowSums((genotype - overall)^2)
  ss_treatment <- 6 * rowSums((treatment - overall)^2)
  interaction <- cells - genotype[, c(1, 1, 2, 2), drop = FALSE] -
    treatment[, c(1, 2, 1, 2), drop = FALSE] + overall
  ss_interaction <- 3 * rowSums(interaction^2)
  mse <- rowSums((y - cells[, rep(1:4, each = 3L), drop = FALSE])^2) / 8
  if (any(mse <= 0)) stop("Zero residual variance in complete-case audit")
  cbind(p_genotype = stats::pf(ss_genotype / mse, 1, 8, lower.tail = FALSE),
    p_treatment = stats::pf(ss_treatment / mse, 1, 8, lower.tail = FALSE),
    p_interaction = stats::pf(ss_interaction / mse, 1, 8, lower.tail = FALSE))
}

build_completecase_audit <- function(inputs, settings, destination) {
  rows <- lapply(c("root", "shoot"), function(tissue) {
    raw <- inputs[[paste0("proteomics_", tissue)]]
    saved <- inputs[[paste0("saved_proteomics_", tissue)]]$data
    keep <- proteomics_keep(raw, settings)
    # First-accession keys are used here solely to audit the legacy export.
    keys <- sub(";.*", "", raw$data$Protein.Group[keep])
    if (anyDuplicated(keys)) stop("Ambiguous legacy accession keys in audit")
    index <- match(saved$Protein.Group, keys)
    if (anyNA(index)) stop("Legacy group absent from filtered input")
    values <- raw$values[keep, , drop = FALSE][index, , drop = FALSE]
    complete <- rowSums(!is.na(values)) == ncol(values)
    calculated <- balanced_protein_anova(values[complete, , drop = FALSE])
    observed <- as.matrix(data.frame(lapply(saved[complete, colnames(calculated)],
      as_numeric_checked, context = "legacy ANOVA P values")))
    data.frame(tissue = raw$spec$tissue, complete_groups_checked = sum(complete),
      max_absolute_P_difference = max(abs(calculated - observed)),
      assessed = "ordinary ANOVA calculations for fully observed groups only")
  })
  write_table(do.call(rbind, rows), destination)
}

build_existing_candidates <- function(inputs, settings, destination) {
  numeric <- c("log.No.0_salt/untreated", "log.NIT2_salt/untreated",
    "log.NIT2/No.0_untreated/untreated", "log.NIT2/No.0_salt/salt",
    "p_FDR_genotype", "p_FDR_treatment", "p_FDR_interaction")
  candidates <- c("NIT2", "AHA1", "AHA2", "SOS1", "NHX7", "HKT1", "HAK5", "AKT1",
    "GORK", "SKOR", "KUP1", "KUP4", "KUP7", "NHX1", "NHX3", "CHX17")
  rows <- lapply(inputs, function(x) {
    d <- x$data
    hit <- vapply(strsplit(ifelse(is.na(d$Genes), "", d$Genes), ";", fixed = TRUE),
      function(v) any(trimws(v) %in% candidates), logical(1))
    out <- d[hit, c("AGI", "Genes", "Protein.Group", numeric), drop = FALSE]
    for (n in numeric) out[[n]] <- as_numeric_checked(out[[n]], n)
    out$tissue <- x$spec$tissue
    out
  })
  write_table(do.call(rbind, rows), destination)
}

build_manuscript_candidates <- function(inputs, settings, destination) {
  rows <- list()
  for (gene in c("NIT2", "NHX7", "AHA2")) for (x in inputs) {
    d <- x$data[x$data$Genes == gene & !is.na(x$data$Genes), , drop = FALSE]
    if (nrow(d) != 1L) stop("Candidate assignment absent or ambiguous: ", gene)
    rows[[length(rows) + 1L]] <- data.frame(Protein = if (gene == "NHX7") "SOS1" else gene,
      Tissue = x$spec$tissue,
      Control_log2FC = as_numeric_checked(d[["log.NIT2/No.0_untreated/untreated"]], "control FC"),
      Salt_log2FC = as_numeric_checked(d[["log.NIT2/No.0_salt/salt"]], "salt FC"),
      Genotype_FDR = checked_q(d, "p_FDR_genotype"), Salt_FDR = checked_q(d, "p_FDR_treatment"),
      Interaction_FDR = checked_q(d, "p_FDR_interaction"))
  }
  write_table(do.call(rbind, rows), destination)
}

build_existing_go_inventory <- function(inputs, settings, destination) {
  rows <- lapply(inputs, function(x) {
    d <- x$data; q <- checked_q(d, "FDR")
    background <- if ("bgtotal" %in% names(d)) d$bgtotal else d$reftotal
    data.frame(source = x$spec$path, terms_reported = nrow(d), terms_FDR_lt_0_05 = sum(q < 0.05, na.rm = TRUE),
      query_total = paste(sort(unique(d$querytotal)), collapse = ";"),
      reference_total = paste(sort(unique(background)), collapse = ";"))
  })
  write_table(do.call(rbind, rows), destination)
}

build_existing_go_terms <- function(inputs, settings, destination) {
  target <- c("GO:0009651", "GO:0006468", "GO:0005215", "GO:0022857", "GO:0016491",
    "GO:0006952", "GO:0030003", "GO:0006813", "GO:0033993", "GO:0005886")
  rows <- lapply(inputs, function(x) {
    d <- x$data; d$FDR <- checked_q(d, "FDR")
    d <- d[d$GO_acc %in% target, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    data.frame(source = x$spec$path, GO = d$GO_acc, term = d$Term, query_hits = d$queryitem,
      query_total = d$querytotal, FDR = d$FDR)
  })
  write_table(do.call(rbind, rows), destination)
}

build_existing_phosphopeptides <- function(inputs, settings, destination) {
  peptide_inputs <- Filter(function(x) x$spec$assay == "legacy_csv", inputs)
  rows <- lapply(peptide_inputs, function(x) x$data)
  d <- do.call(rbind, rows)
  key <- paste(d$Protein.Group, d$Precursor.Id, sep = "|")
  d <- d[!duplicated(key), , drop = FALSE]
  samples <- grep("^(No\\.0|NIT2)_[0-9]+mM_[123]$", names(d), value = TRUE)
  if (length(samples) != 12L) stop("Expected twelve phosphopeptide sample columns")
  values <- as.matrix(data.frame(lapply(d[samples], as_numeric_checked, context = "phosphopeptide")))
  out <- d[c("AGI", "Genes", "Protein.Group", "Modified.Sequence", "Precursor.Id")]
  out$observed_samples <- rowSums(!is.na(values)); out$missing_samples <- rowSums(is.na(values))
  out$reference_phosphosite <- NA_character_
  references <- Filter(function(x) x$spec$assay == "legacy_fasta", inputs)
  for (ref in references) {
    hit <- which(out$Protein.Group == ref$data$accession)
    for (i in hit) {
      mapped <- map_modified_peptide(out$Modified.Sequence[i], ref$data$sequence)
      out$reference_phosphosite[i] <- paste(mapped$phosphosites, collapse = ";")
    }
  }
  write_table(out, destination)
}

map_modified_peptide <- function(modified, reference) {
  peptide <- gsub("\\([^)]*\\)", "", modified)
  positions <- gregexpr(peptide, reference, fixed = TRUE)[[1L]]
  if (length(positions) != 1L || positions[1L] < 1L) stop("Peptide mapping is absent or ambiguous")
  tokens <- regmatches(modified, gregexpr("[A-Z]|\\([^)]*\\)", modified, perl = TRUE))[[1L]]
  amino_acid <- 0L; sites <- integer()
  for (token in tokens) {
    if (grepl("^[A-Z]$", token)) amino_acid <- amino_acid + 1L
    else if (token == "(UniMod:21)") sites <- c(sites, positions[1L] + amino_acid - 1L)
  }
  list(start = positions[1L], phosphosites = sites)
}

build_existing_dose_figure <- function(inputs, settings, destination) {
  counts <- matrix(NA_real_, 2L, 2L, dimnames = list(c("Root", "Shoot"), c("Linear trend", "Quadratic joint test")))
  for (x in inputs) {
    q <- checked_q(x$data, if ("adj.P_joint" %in% names(x$data)) "adj.P_joint" else "adj.P.Val")
    column <- if ("adj.P_joint" %in% names(x$data)) 2L else 1L
    counts[x$spec$tissue, column] <- sum(q < 0.05, na.rm = TRUE)
  }
  grDevices::pdf(destination, width = 7, height = 5, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(5, 5, 3, 1))
  bars <- graphics::barplot(counts, beside = TRUE, col = c("#4477AA", "#CC6677"),
    ylim = c(0, max(counts) * 1.25), ylab = "Mapped transcript IDs with FDR < 0.05",
    main = "Existing dose response results")
  graphics::text(bars, counts, labels = counts, pos = 3, cex = 0.9)
  graphics::legend("topleft", legend = rownames(counts), fill = c("#4477AA", "#CC6677"), bty = "n")
}

build_existing_proteomics_figure <- function(inputs, settings, destination) {
  counts <- matrix(NA_real_, 2L, 3L, dimnames = list(c("Root", "Shoot"), c("Genotype", "Salt", "Interaction")))
  for (x in inputs) counts[x$spec$tissue, ] <- vapply(
    c("p_FDR_genotype", "p_FDR_treatment", "p_FDR_interaction"),
    function(n) sum(checked_q(x$data, n) < 0.05, na.rm = TRUE), 1L)
  grDevices::pdf(destination, width = 7, height = 5, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(5, 5, 3, 1))
  bars <- graphics::barplot(counts, beside = TRUE, col = c("#4477AA", "#CC6677"),
    ylim = c(0, max(counts) * 1.3), ylab = "Protein groups with FDR < 0.05",
    main = "Existing proteomics model terms")
  graphics::text(bars, counts, labels = counts, pos = 3, cex = 0.9)
  graphics::legend("topright", legend = rownames(counts), fill = c("#4477AA", "#CC6677"), bty = "n")
}
