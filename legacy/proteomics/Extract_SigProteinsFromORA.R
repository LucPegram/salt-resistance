# Extracts significant proteins from ORA output, returning data contained about them
# from the AnovaPvalues files
# Luc Pegram
# 2026/05/28

## Tissue
#===========================
tissue = "Shoot"
#===========================
tissue = "Root"
#===========================

## Type
#===========================
type = "Genotype"
#===========================
type = "Treatment"
#===========================
type = "Interaction"
#===========================


input_file_ORA <- paste0(
  "D:/Masters_PhD/Data/Proteomics/2026/ForAnalysis/R_Exports/ORA/",
  tissue,
  "/",
  type,
  "/",
  tissue,
  "_ORA_LLPS.csv")

input_file_anova <- paste0("D:/Masters_PhD/Data/Proteomics/2026/ForAnalysis/R_Exports/",
                     tissue,
                     "Tissue_AnovaPvalues.csv")

output_dir <- paste0(
  "D:/Masters_PhD/Data/Proteomics/2026/ForAnalysis/R_Exports/ORA/",
  tissue,
  "/",
  type,
  "/",
  "Data/")

output_file <- paste0(output_dir, tissue, "_LLPS_Data.csv")


# ==============================================================================
# LOAD AND PREPARE DATA
# ==============================================================================

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ora   <- read.csv(input_file_ORA, header = TRUE, stringsAsFactors = FALSE,
                  check.names = FALSE)

og_df <- read.csv(input_file_anova, header = TRUE, stringsAsFactors = FALSE,
                  check.names = FALSE)


# ── Extract all unique AGIs from ora$Sig_Proteins ─────────────────────────────
# Each cell contains AGIs separated by "; " — split, flatten, and deduplicate
sig_agi <- unique(strsplit(ora$Sig_Proteins[1], split = "; "))

## Remove any empty strings that may result from splitting
#sig_agi <- sig_agi[nzchar(sig_agi)]

# ── Filter og_df ──────────────────────────────────────────────────────────────
filtered_df <- og_df[og_df$AGI %in% sig_agi[[1]], ]

cat("Rows matched:", nrow(filtered_df), "out of", nrow(og_df), "\n")

# ── Write output ───────────────────────────────────────────────────────────────

write.csv(filtered_df, output_file, row.names = FALSE)

cat("Saved to:", output_file, "\n")
