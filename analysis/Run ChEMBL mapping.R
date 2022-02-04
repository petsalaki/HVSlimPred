library(HVSlimPred)
chembl_data = read.csv("data/input/ChEMBL/filtred_chembl_map.csv.gz")
Drug_indication = read.csv("data/input/ChEMBL/Drug_indication_chembl_27.csv")
# Host viral ChEMBL mapping -----------------------------------------------
HV_final_hits = readRDS("data/output/HV_final_hits.rds")
HV_final_clinvar_path = readRDS("data/output/HV_Clinvar_and_PTM.rds")
HV_final_clinvar_path  = HV_final_clinvar_path$ClinVar_path$ClinVar_path

HV_H1_chembl = H1_to_ChEMBL(slim_hits = HV_final_hits, clinvar_df = HV_final_clinvar_path, host_viral = T, chemble_source = chembl_data, Drug_inds = Drug_indication, indirect_mapping_dist = 3)
saveRDS(HV_H1_chembl, "data/output/HV_H1_chembl.rds")

# Human only ChEMBL mapping -----------------------------------------------
human_final_hits = readRDS("data/output/human_final_hits.rds")
human_final_clinvar_path = readRDS("data/output/human_Clinvar_and_PTM.rds")
human_final_clinvar_path  = human_final_clinvar_path$ClinVar_path$ClinVar_path

human_H1_chembl = H1_to_ChEMBL(slim_hits = human_final_hits, clinvar_df = human_final_clinvar_path, host_viral = F, chemble_source = chembl_data, Drug_inds = Drug_indication, indirect_mapping_dist = 3)
saveRDS(human_H1_chembl, "data/output/human_H1_chembl.rds")


