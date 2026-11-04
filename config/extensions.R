# Public reference resources, not project measurements or derived gene lists.
extension_settings <- function() list(seed=20261104L,
  salt_signature_q=.05, salt_signature_abs_log2FC=log2(1.5),
  signature_min_members=10L, marker_min_members=5L,
  TF_min_members=10L, TF_max_universe_fraction=.5,
  camera_correlations=c(primary=.01,sensitivity=.05),
  atlas_top_per_cluster=50L, atlas_min_log2FC=.5, atlas_max_adj_P=.05)

extension_resources <- function() {
  geo <- "https://ftp.ncbi.nlm.nih.gov/geo/series/"
  publisher <- paste0("https://media.springernature.com/full/springer-static/esm/",
    "art%3A10.1038%2Fs41477-025-02072-z/MediaObjects/")
  c(
    GSE5620_series_matrix.txt.gz=paste0(geo,"GSE5nnn/GSE5620/matrix/GSE5620_series_matrix.txt.gz"),
    GSE5623_series_matrix.txt.gz=paste0(geo,"GSE5nnn/GSE5623/matrix/GSE5623_series_matrix.txt.gz"),
    GSE223207_series_matrix.txt.gz=paste0(geo,"GSE223nnn/GSE223207/matrix/GSE223207_series_matrix.txt.gz"),
    GPL198.soft="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL198&targ=self&form=text&view=full",
    GPL17133.soft="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL17133&targ=self&form=text&view=full",
    AtRegNet.zip="https://agris-knowledgebase.org/Downloads/AtRegNet.zip",
    atlas2025_markers.xlsx=paste0(publisher,"41477_2025_2072_MOESM3_ESM.xlsx"),
    atlas2025_clusters.xlsx=paste0(publisher,"41477_2025_2072_MOESM4_ESM.xlsx"))
}
