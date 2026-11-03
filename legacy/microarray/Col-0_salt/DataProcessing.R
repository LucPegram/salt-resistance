# Data Processing, converting the data file into long format.
# Written by Claude AI and Luc Pegram.
# Input columns: "CORRECTED MAPPING", then "<treatment> root <batch>" columns
# Output columns: AGI MAPPING | treatment | batch | values
# Processing order:
#   1. log-transform all values (optional)
#   2. average rows sharing the same CORRECTED MAPPING ID (per column, NA ignored)
#   3. pivot to long format
# Averaging is done after the log transform, so duplicate probes are averaged on the
# log scale (the same approach as limma::avereps).
# 2026/10/28

library(readr)
library(dplyr)
library(tidyr)

rm(list=ls())

##---------------------------
#tissue = 'Root'

for (tissue in c("Shoot", "Root")) {

  
input_file <- paste0("D:/Masters_PhD/Data/Microarray Data_From Lara/My Analysis/Col-0_SaltHydroponics/",
                     tissue,
                     "/Col-0_SaltHydro_",
                     tissue,
                     ".csv")

output_file <- paste0("D:/Masters_PhD/Data/Microarray Data_From Lara/My Analysis/Col-0_SaltHydroponics/",
                      tissue,
                      "/Col-0_",
                      tissue,
                      "_long.csv")

log_transform    <- TRUE   # set FALSE if the values are already log-scale
log_base         <- 2      # 2 = log2 (consistent with linear_trend_analysis.R)
log_pseudocount  <- 0      # added before logging; use e.g. 1 if the data contain zeros

# ---- Read -------------------------------------------------------------------
# Name repair is off so the original headings are preserved; all data columns are
# read as numeric (any non-numeric entry is reported and stops the script)
wide <- read_csv(
  input_file,
  col_types      = cols(.default = col_double(), `CORRECTED MAPPING` = col_character()),
  name_repair    = "minimal",
  show_col_types = FALSE
)

nas <- which(wide$`CORRECTED MAPPING`=='#N/A')  #removes missing AGI mapping
wide <- wide[-nas,]

if (nrow(problems(wide)) > 0) {
  print(problems(wide))
  stop("Some values could not be read as numbers (see above)")
}

# Sanity check: every data column heading must look like "PN root B1" or "50 NaCl root B1"
data_cols <- setdiff(names(wide), "CORRECTED MAPPING")
pattern   <- "^(PN|[0-9]+).*(B[0-9]+)$"
bad <- data_cols[!grepl(pattern, data_cols)]
if (length(bad) > 0) {
  stop("Column headings not matching expected format: ", paste(bad, collapse = ", "))
}

# ---- 1. Log transformation --------------------------------------------------
if (log_transform) {
  vals <- as.matrix(wide[, data_cols])
  if (any(vals < 0, na.rm = TRUE)) {
    stop("Negative values found: cannot log-transform. Set log_transform <- FALSE ",
         "if the data are already log-scale or background-corrected.")
  }
  if (log_pseudocount == 0 && any(vals == 0, na.rm = TRUE)) {
    stop(sum(vals == 0, na.rm = TRUE), " zero values found: log(0) is -Inf. ",
         "Set log_pseudocount (e.g. 1) and rerun.")
  }
  wide <- wide %>%
    mutate(across(all_of(data_cols), ~ log(.x + log_pseudocount, base = log_base)))
  message("Values transformed: log", log_base, "(x", if (log_pseudocount > 0) paste0(" + ", log_pseudocount), ")")
}

# ---- 2. Average duplicated gene IDs -----------------------------------------
id_order <- unique(wide$`CORRECTED MAPPING`)       # keep original gene order
n_rows   <- nrow(wide)
n_unique <- length(id_order)
message(n_rows, " input rows; ", n_unique, " unique IDs; ",
        n_rows - n_unique, " duplicated-ID rows averaged into their first occurrence")

wide <- wide %>%
  group_by(`CORRECTED MAPPING`) %>%
  summarise(
    across(all_of(data_cols),
           ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(match(`CORRECTED MAPPING`, id_order))

# ---- 3. Wide -> long --------------------------------------------------------
long <- wide %>%
  rename(`AGI MAPPING` = `CORRECTED MAPPING`) %>%
  pivot_longer(
    cols           = all_of(data_cols),
    names_to       = c("treatment", "batch"),
    names_pattern  = pattern,
    values_to      = "values"
  ) %>%
  mutate(
    treatment = factor(treatment, levels = c("PN", "50", "75", "100")),
    batch     = factor(batch,     levels = c("B1", "B2", "B3"))
  ) %>%
  select(`AGI MAPPING`, treatment, batch, values)

# Checks: headings within expected levels; one row per AGI x treatment x batch
stopifnot(!anyNA(long$treatment), !anyNA(long$batch))
stopifnot("Duplicated AGI x treatment x batch combinations remain" =
            !anyDuplicated(long[, c("AGI MAPPING", "treatment", "batch")]))
stopifnot(nrow(long) == n_unique * length(data_cols))

#print(long, n = 10)

write_csv(long, output_file)
}
