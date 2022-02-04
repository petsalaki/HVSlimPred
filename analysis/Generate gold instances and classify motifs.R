library(HVSlimPred)

# Host viral gold and classification --------------------------------------
HV_clinvar_path = readRDS("data/output/HV_Clinvar_and_PTM.rds")
HV_PTM_info = HV_clinvar_path$PTM
HV_clinvar_path = HV_clinvar_path$ClinVar_path$ClinVar_path
All_pepsite = readRDS("data/output/HV_pepsite_final_with_AF.rds")
HV_final_hits = readRDS("data/output/HV_final_hits.rds")
HV_chembl = readRDS("data/output/HV_H1_chembl.rds")
HV_chembl = HV_chembl$H1_and_Clinvar
HV_chembl = HV_chembl[which(HV_chembl$distance <= 2),]
HV_HMM_iELM = readRDS("data/output/HV_HMM_iELM.rds")
pfam_doms = readRDS("data/input/pfam_doms_human.rds")
HV_Enriched_doms = readRDS("data/output/HV_enriched_doms_with_emppval.rds")
HV_all_pepsite = readRDS("data/output/HV_pepsite_final_with_AF.rds")
HV_input_prots = readRDS("data/input/HV_proteins_per_dataset_allhitsint.RDS")

HMM_insts = HMM_iELM_classify(new_hits_ints = HV_final_hits, HMM_iELM = HV_HMM_iELM,host_viral = T)
HMM_insts = HMM_insts[which(HMM_insts$passed_iELM == 1),]

filter_combs = prepare_required_data_OR(EM = HV_Enriched_doms, clin_path = HV_clinvar_path, iELM = HMM_insts,
                                        new_hits_ints = HV_final_hits, pepsite_res = HV_all_pepsite,
                                        min_domain = 5,adj.pval = 0.05,slim_input_prots = HV_input_prots, prot_only = F)

HV_gold_insts = generate_potential_candidates(pep_clin_EM = filter_combs$ABC, clinvar_path = HV_clinvar_path,
                                              PTM_info = HV_PTM_info,
                                              pfam_doms = pfam_doms, new_hits_ints = HV_final_hits, compari_cutoff = 0.66)

HV_classified_insts = get_classified_insts(new_hits_ints = HV_final_hits, clinvar_path = HV_clinvar_path,
                                           PTM_info = HV_PTM_info, pepsite_hits = HV_all_pepsite,
                                           Drug_hits = HV_chembl, HMM_iELM = HV_HMM_iELM, EM_list = HV_Enriched_doms)

saveRDS(HV_gold_insts, "data/output/HV_gold_insts.rds")
saveRDS(HV_classified_insts, "data/output/HV_classified_insts.rds")


# Human only gold and classification --------------------------------------
human_clinvar_path = readRDS("data/output/human_Clinvar_and_PTM.rds")
human_PTM_info = human_clinvar_path$PTM
human_clinvar_path = human_clinvar_path$ClinVar_path$ClinVar_path
All_pepsite = readRDS("data/output/All_human_only_pepsite_with_AF.rds")
human_final_hits = readRDS("data/output/human_final_hits.rds")
human_chembl = readRDS("data/output/human_H1_chembl.rds")
human_chembl = human_chembl$H1_and_Clinvar
human_chembl = human_chembl[which(human_chembl$distance <= 2),]
human_HMM_iELM = readRDS("data/output/human_HMM_iELM.rds")
pfam_doms = readRDS("data/input/pfam_doms_human.rds")
human_Enriched_doms = readRDS("data/output/human_only_enriched_doms.rds")
human_input_prots = readRDS("data/input/human_only_qslim_datasets.rds")

HMM_insts = HMM_iELM_classify(new_hits_ints = human_final_hits, HMM_iELM = human_HMM_iELM,host_viral = F)
HMM_insts = HMM_insts[which(HMM_insts$passed_iELM == 1),]

filter_combs = prepare_required_data_OR(EM = human_Enriched_doms, clin_path = human_clinvar_path, iELM = HMM_insts,
                                        new_hits_ints = human_final_hits, pepsite_res = All_pepsite,
                                        min_domain = 5,adj.pval = 0.05,slim_input_prots = human_input_prots, prot_only = F)

human_gold_insts = generate_potential_candidates(pep_clin_EM = filter_combs$ABC, clinvar_path = human_clinvar_path,
                                              PTM_info = human_PTM_info,
                                              pfam_doms = pfam_doms, new_hits_ints = human_final_hits, compari_cutoff = 0.66)

human_classified_insts = get_classified_insts(new_hits_ints = human_final_hits, clinvar_path = human_clinvar_path,
                                           PTM_info = human_PTM_info, pepsite_hits = All_pepsite,
                                           Drug_hits = human_chembl, HMM_iELM = human_HMM_iELM, EM_list = human_Enriched_doms, host_viral = F)

saveRDS(human_gold_insts, "data/output/human_gold_insts.rds")
saveRDS(human_classified_insts, "data/output/human_classified_insts.rds")
