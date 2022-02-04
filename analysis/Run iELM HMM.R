library(HVSlimPred)

#iELM HMMs was downloaded from http://elmint.embl.de/program_file/ and hmmsearch from HMMER3 toolkit was used to search each HMM profile against
# H1 proteins (see methods in manuscript for more details).

#H1 prots for both Host-viral and human-only approach can be retrieved from the final hits.
HV_final_hits = readRDS("data/output/HV_final_hits.rds")
human_final_hits = readRDS("data/output/human_final_hits.rds")

HV_H1_prots = HV_final_hits %>% tidyr::separate(Dataset, into = c("V1", "H1"), sep = "_") %>% dplyr::select(H1)
HV_H1_prots = unique(HV_H1_prots$H1)

human_H1_prots = unique(human_final_hits$Dataset)

#Using the above H1 prots, the bash script "analysis/H1_HMMER.sh" can be used to run HMMER3 against both pfam and iELM HMMs downloaded above.
# For the sake of reproducibility, the output of HMMER3 has been provided in "data/input/iELM".

ELM_renamed = read.delim("data/input/ELM/ELM_renamed.tsv")
HV_HMM_iELM = iELM_HMM(path_HMM_pfam_domtblouts = "data/input/iELM/HV/pfam_H1_domtblouts/", path_HMM_iELM_domtblouts = "data/input/iELM/HV/domain_all_H1_domtblouts/", ELM_renamed_classes = ELM_renamed)
human_HMM_iELM = iELM_HMM(path_HMM_pfam_domtblouts = "data/input/iELM/Human_only/pfam_H1_domtblouts/", path_HMM_iELM_domtblouts = "data/input/iELM/Human_only/domain_all_H1_domtblouts/", ELM_renamed_classes = ELM_renamed)

saveRDS(HV_HMM_iELM, "data/output/HV_HMM_iELM.rds")
saveRDS(human_HMM_iELM, "data/output/human_HMM_iELM.rds")
