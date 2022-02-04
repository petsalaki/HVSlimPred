library(HVSlimPred)

# Preprocess slimfinder/qslimfinder output by providing both occurences 'slimfinder.occ.csv' and main output 'slimfinder.csv'
# If no arguments are given, the main output will be downloaded from BioStudies
HV_main_out = read.csv("data/input/HV_allhitsperint.csv")
HV_occ = read.csv("data/input/HV_allhitsperint.occ.csv")

HV_cons_main_out = read.csv("data/input/HV_cons_allhitsperint.csv")
HV_cons_occ = read.csv("data/input/HV_cons_allhitsperint.occ.csv")

human_main_out = read.csv("data/input/human_allhitsint.csv")
human_occ = read.csv("data/input/human_allhitsint.occ.csv")

HV_all_hits = preprocess_slimfinder_out(main_output = HV_main_out, occ_out = HV_occ)
HV_cons_all_hits = preprocess_slimfinder_out(main_output = HV_cons_main_out, occ_out = HV_cons_occ)
human_all_hits = preprocess_slimfinder_out(main_output = human_main_out, occ_out = human_occ)

HV_cons_all_hits$all_hitsint = dplyr::select(HV_cons_all_hits$all_hitsint, -HomNum, -GlobID, -LocID)


#Prepare input for CompariMotif comparison of motif patterns to ELM classes
ELM_classes = Get_ELM_all(ELM_data_type = "classes")
HV_compari_motif_input = prepare_input_comparimotif(pred_hits = HV_all_hits[["all_hitsint"]], ELM_classes = ELM_classes)
human_compari_motif_input = prepare_input_comparimotif(pred_hits = human_all_hits[["all_hitsint"]], ELM_classes = ELM_classes)


output_path = "data/output"
write.table(HV_compari_motif_input[["hits_to_CompariMotif"]],
            file = file.path(output_path,"HV_allhits_slim_compari.txt"), row.names = F)
write.table(HV_compari_motif_input[["ELM_to_CompariMotif"]],
            file = file.path(output_path,"compari_ELM_ids.txt"), row.names = F)

write.table(human_compari_motif_input[["hits_to_CompariMotif"]],
            file = file.path(output_path,"human_allhits_slim_compari.txt"), row.names = F)

# Run the following command on a cluster or local shell to run CompariMotif comparison from the input above
#python /path/to/comparimotif/comparimotif_V3.py motifs=/data/output/HV_allhits_slim_compari.txt searchdb=/data/output/compari_ELM_ids.txt unmatched=T overlaps=F xgmml=F
#python /path/to/comparimotif/comparimotif_V3.py motifs=/data/output/human_allhits_slim_compari.txt searchdb=/data/output/compari_ELM_ids.txt unmatched=T overlaps=F xgmml=F


#Preprocess the comparimotif .tdt output and parse the final
HV_compari_output = readRDS("data/input/compari_motif_ELM_ids.RDS")
human_compari_output = read.table("data/input/human_allhits_slim_compari-compari_ELM_ids.compare.tdt", sep = "\t", header = T)

HV_final_hits = parse_qslim(main_output = HV_main_out, occ_out = HV_occ, compari_output = HV_compari_output, preprocess_compari = T, select_top_rank = F)
match_cols = c("Dataset", "Pattern", "uniprot", "Start_Pos", "End_Pos", "Match", "Cons")
x = HV_cons_all_hits$all_hitsint[,match_cols]
HV_final_hits = dplyr::left_join(HV_final_hits, x)

human_final_hits = parse_qslim(main_output = human_main_out, occ_out = human_occ, compari_output = human_compari_output, preprocess_compari = T, select_top_rank = F)
saveRDS(HV_final_hits, "data/output/HV_final_hits.rds")
saveRDS(human_final_hits, "data/output/human_final_hits.rds")

