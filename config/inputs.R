# Paths are relative to the external SALT_DATA_ROOT, never to a user's computer.
input_specs <- list(
  microarray_root = list(
    path = paste0("Microarray_for_Stephano/Luc_Initial_Analysis/",
                  "Col-0_SaltHydroponics/Root/Col-0_SaltHydro_Root.csv"),
    assay = "microarray", tissue = "Root", ecotype = "Col-0"),
  microarray_shoot = list(
    path = paste0("Microarray_for_Stephano/Luc_Initial_Analysis/",
                  "Col-0_SaltHydroponics/Shoot/Col-0_SaltHydro_Shoot.csv"),
    assay = "microarray", tissue = "Shoot", ecotype = "Col-0"),
  microarray_genotype = list(
    path = paste0("Microarray_for_Stephano/Luc_Initial_Analysis/",
                  "NIT2vsNo-0_Treatments/",
                  "NIT2vsNo-0_Excel Normalized (31-01).csv"),
    assay = "microarray", tissue = "Combined", ecotype = "No-0"),
  proteomics_root = list(
    path = "Microarray_for_Stephano/Proteomics_Files/RootTissue_pg.tsv",
    assay = "proteomics", tissue = "Root", ecotype = "No-0"),
  proteomics_shoot = list(
    path = "Microarray_for_Stephano/Proteomics_Files/ShootTissue_pg.tsv",
    assay = "proteomics", tissue = "Shoot", ecotype = "No-0"),
  qpcr_summary = list(
    path = paste0("Microarray_for_Stephano/Paul's_qPCRs/",
                  "Transporters&channels_Average gene expression data_Paul thesis.xlsx"),
    assay = "qpcr_summary", sheet = "Average Normalized expression")
)
input_specs <- c(input_specs, legacy_result_specs())
input_specs$gene_annotation <- list(path="Microarray_for_Stephano/Stu_MappingDescription.csv",
  assay="gene_annotation")
input_specs$aha2_reference <- list(path = "references/P19456.fasta", assay = "legacy_fasta", accession = "P19456")
input_specs$go_associations <- list(
  path = "Data_from_Luc_2/Proteomics_AgriGO_Analysis/317_slimGO.tsv", assay = "go_associations")


# Exact legacy GO output paths; only file specifications, not results.
legacy_go_paths <- c(
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Root/All/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Root/Genotype_sig_pVal(0.05)/GenotypeFC_log2(1.5)/vs_Genome/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Root/Genotype_sig_pVal(0.05)/GenotypeFC_log2(1.5)/vs_ProteomicData/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Root/Genotype_sig_pVal(0.05)/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Root/Phospho/PhosphoPeptides_vs_DataSet/NA_4_or_less/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Root/Phospho/PhosphoPeptides_vs_DataSet/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Root/Phospho/PhosphoPeptides_vs_Genome/NA_4_or_less/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Root/Phospho/PhosphoPeptides_vs_Genome/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Root/sig_pVal(0.1)/FC_Genotype/1.5fold_UP&DOWN/vsDataSet/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Root/sig_pVal(0.1)/FC_Genotype/1.5fold_UP&DOWN/vsWholeGenome/old/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Root/sig_pVal(0.1)/FC_Genotype/1.5fold_UP&DOWN/vsWholeGenome/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Root/sig_pVal(0.1)/FC_Genotype/2fold_UP&DOWN/vsDataSet/old/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Root/sig_pVal(0.1)/FC_Genotype/2fold_UP&DOWN/vsWholeGenome/old/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Root/sig_pVal(0.1)/FC_SaltGentype/2fold_UP&DOWN/old/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Shoot/All/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Shoot/Genotype_sig_pVal(0.05)/GenotypeFC_log2(1.5)/vs_Genome/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Shoot/sig_pVal(0.1)/FC_Genotype/1.5fold_UP&DOWN/vsWholeGenome/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Shoot/sig_pVal(0.1)/FC_Genotype/2fold_UP&DOWN/vsDataSet/old/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Shoot/sig_pVal(0.1)/FC_Genotype/2fold_UP&DOWN/vsWholeGenome/old/results.txt",
  "Data_from_Luc_2/Proteomics_AgriGO_Analysis/Shoot/sig_pVal(0.1)/FC_SaltGentype/2fold_UP&DOWN/old/results.txt",
  "Microarray_for_Stephano/Luc_Initial_Analysis/AgriGO/NIT2vs_No-0/adj_pVal_0.01/results.txt",
  "Microarray_for_Stephano/Luc_Initial_Analysis/AgriGO/NIT2vs_No-0/adj_pVal_0.01/vsFullMicroarray/results.txt",
  "Microarray_for_Stephano/Luc_Initial_Analysis/AgriGO/Root/FDR_pVal_0.001/results.txt",
  "Microarray_for_Stephano/Luc_Initial_Analysis/AgriGO/Root/FDR_pVal_0.01/results.txt",
  "Microarray_for_Stephano/Luc_Initial_Analysis/AgriGO/Root/FDR_pVal_0.01/vs_microarray/results.txt",
  "Microarray_for_Stephano/Luc_Initial_Analysis/AgriGO/Shoot/FDR_pVal_0.01/results.txt",
  "Microarray_for_Stephano/Luc_Initial_Analysis/AgriGO/Shoot/FDR_pVal_0.01/vsMicroarray/results.txt"
)
for (i in seq_along(legacy_go_paths)) {
  input_specs[[sprintf("saved_go_%02d", i)]] <- list(path = legacy_go_paths[i], assay = "legacy_go")
}
