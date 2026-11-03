write_table <- function(x, destination) {
  write.csv(x, destination, row.names = FALSE, na = "NA")
}

build_input_inventory <- function(inputs, settings, destination) {
  rows <- lapply(names(inputs), function(id) {
    x <- inputs[[id]]
    data.frame(input = id, assay = x$spec$assay, rows = nrow(x$data),
      distinct_feature_ids = switch(x$spec$assay,
        microarray = length(unique(x$ids[x$mapped])),
        proteomics = length(unique(x$data$Protein.Group)),
        qpcr_summary = length(unique(x$data$Gene))),
      sample_columns = if (x$spec$assay == "qpcr_summary") NA_integer_ else length(x$samples),
      unit = if (x$spec$assay == "qpcr_summary") "group summary; biological n unavailable" else
        "biological plant-pool measurement", stringsAsFactors = FALSE)
  })
  write_table(do.call(rbind, rows), destination)
}

build_mapping_audit <- function(inputs, settings, destination) {
  rows <- lapply(names(inputs), function(id) {
    x <- inputs[[id]]
    collapsed <- collapse_microarray(x)
    data.frame(input = id, probe_rows = length(x$ids), mapped_rows = sum(x$mapped),
      unmapped_rows = sum(!x$mapped), distinct_transcript_ids = nrow(collapsed),
      repeated_mapped_rows = sum(x$mapped) - nrow(collapsed),
      missing_expression_cells = sum(is.na(x$values)), stringsAsFactors = FALSE)
  })
  write_table(do.call(rbind, rows), destination)
}

build_proteomics_filtering <- function(inputs, settings, destination) {
  rows <- lapply(names(inputs), function(id) {
    x <- inputs[[id]]
    keep <- proteomics_keep(x, settings)
    data.frame(input = id, original_protein_groups = nrow(x$data),
      retained_protein_groups = sum(keep), excluded_protein_groups = sum(!keep),
      retained_missing_cells = sum(is.na(x$values[keep, , drop = FALSE])),
      ambiguous_retained_groups = sum(grepl(";", x$data$Protein.Group[keep], fixed = TRUE)),
      min_observed_in_one_condition = settings$proteomics_min_observed_in_one_condition,
      min_sequences = settings$proteomics_min_sequences, stringsAsFactors = FALSE)
  })
  write_table(do.call(rbind, rows), destination)
}

build_proteomics_missingness <- function(inputs, settings, destination) {
  grDevices::pdf(destination, width = 11, height = 7, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(mfrow = c(length(inputs), 1L), mar = c(7, 5, 3, 1))
  on.exit(graphics::par(old), add = TRUE, after = FALSE)
  for (x in inputs) {
    missing <- 100 * colMeans(is.na(x$values))
    graphics::barplot(missing, names.arg = x$samples, las = 2, cex.names = 0.8,
      col = ifelse(x$metadata$genotype == "No-0", "#4477AA", "#AA3377"),
      ylim = c(0, 100), ylab = "Missing measurements (%)",
      main = paste(x$spec$tissue, "protein groups: before filtering or imputation"))
  }
}

build_qpcr_table <- function(inputs, settings, destination) {
  write_table(inputs$qpcr_summary$data, destination)
}

build_qpcr_figure <- function(inputs, settings, destination) {
  d <- inputs$qpcr_summary$data
  expected <- expand.grid(Treatment = settings$qpcr_treatment_labels,
    Genotype = settings$qpcr_genotype_labels, Tissue = settings$tissue_levels,
    stringsAsFactors = FALSE)
  key <- function(z) paste(z$Tissue, z$Genotype, z$Treatment, sep = "|")
  expected_key <- key(expected)
  genes <- if (is.null(settings$qpcr_genes)) sort(unique(d$Gene)) else settings$qpcr_genes
  if (!all(genes %in% d$Gene)) stop("Selected qPCR gene absent from input")
  ordered <- lapply(genes, function(gene) {
    z <- d[d$Gene == gene, , drop = FALSE]
    if (!setequal(key(z), expected_key) || nrow(z) != nrow(expected)) {
      stop("Incomplete or unexpected qPCR groups for ", gene)
    }
    z[match(expected_key, key(z)), , drop = FALSE]
  })
  columns <- if (is.null(settings$qpcr_panel_columns)) 4L else settings$qpcr_panel_columns
  rows <- ceiling(length(genes) / columns)
  grDevices::pdf(destination, width = if (columns == 2L) 8 else 13,
                  height = max(4, rows * 3.5), useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(mfrow = c(rows, columns), mar = c(7, 4, 3, 1),
                        oma = c(0, 0, 3, 0))
  on.exit(graphics::par(old), add = TRUE, after = FALSE)
  for (i in seq_along(genes)) {
    z <- ordered[[i]]
    mean <- z[["Average Normalized expression"]]
    se <- z$STE
    labels <- paste(z$Tissue, ifelse(z$Genotype == "No-0", "WT", "NIT2"),
                    ifelse(z$Treatment == "Untreated", "0", "75"))
    limits <- range(c(0, mean, mean - se, mean + se), na.rm = TRUE)
    if (diff(limits) == 0) limits <- c(0, 1)
    limits[2] <- limits[2] * 1.1
    graphics::plot(seq_along(mean), mean, ylim = limits, xaxt = "n", xlab = "",
      ylab = "Reported normalized expression", main = paste0(LETTERS[i], "  ", genes[i]), pch = 19,
      col = ifelse(z$Treatment == "Untreated", "#4477AA", "#CC6677"))
    graphics::axis(1, at = seq_along(mean), labels = labels, las = 2, cex.axis = 0.65)
    has_error <- !is.na(se) & se > 0
    if (any(has_error)) {
      positions <- which(has_error)
      lower <- (mean - se)[has_error]
      upper <- (mean + se)[has_error]
      graphics::segments(positions, lower, positions, upper)
      graphics::segments(positions - 0.06, lower, positions + 0.06, lower)
      graphics::segments(positions - 0.06, upper, positions + 0.06, upper)
    }
    if (anyNA(se)) graphics::mtext("Some STE values not reported", side = 3, cex = 0.6)
  }
  graphics::mtext("qPCR group means; STE where reported; 0/75 mM NaCl",
                   outer = TRUE, line = 1, cex = 1)
  # Summary values cannot recover biological n or support replicate-level tests.
}

