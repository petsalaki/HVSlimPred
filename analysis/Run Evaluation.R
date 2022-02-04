library(HVSlimPred)
HV_enriched_doms = readRDS("data/output/HV_enriched_doms_with_emppval.rds")
human_enriched_doms = readRDS("data/output/human_only_enriched_doms.rds")

HV_final_hits = readRDS("data/output/HV_final_hits.rds")
human_final_hits = readRDS("data/output/human_final_hits.rds")

HV_all_pfam_doms = prepare_pfam_doms_enrich(HV_enriched_doms)
human_all_pfam_doms = prepare_pfam_doms_enrich(human_enriched_doms)

HV_qslim_datasets = readRDS("data/input/HV_proteins_per_dataset_allhitsint.RDS")
human_qslim_datasets = readRDS("data/input/human_only_qslim_datasets.rds")

pwm_prmdb_scan = read.csv("data/input/PRMDB/pwm_prmdb_scan.csv.gz")
all_seqs = read.csv("data/input/uniprot_sequences.csv.gz")
prm_pfam_final = read.csv("data/input/PRMDB/prm_pfam_final.csv.gz")


# Host-viral evaluation ---------------------------------------------------
HV_prot_level_eval_ELM = HV_prot_level_eval(slim_hits = HV_final_hits, bench_data_type = "ELM", all_pfam_dom = HV_all_pfam_doms, datasets_input = HV_qslim_datasets)
HV_prot_level_eval_PRMdb = HV_prot_level_eval(slim_hits = HV_final_hits,bench_data = pwm_prmdb_scan ,bench_data_type = "PRMdb", all_pfam_dom = HV_all_pfam_doms, datasets_input = HV_qslim_datasets)

HV_HH_level_eval_ELM = HV_motif_level_eval(slim_hits = HV_final_hits,bench_data_type = "ELM", all_pfam_dom = HV_all_pfam_doms, datasets_input = HV_qslim_datasets, prot_seqs = all_seqs)
HV_HH_level_eval_PRMdb = HV_motif_level_eval(slim_hits = HV_final_hits, bench_data = pwm_prmdb_scan, bench_data_type = "PRMdb", all_pfam_dom = HV_all_pfam_doms, datasets_input = HV_qslim_datasets, prot_seqs = all_seqs)

HV_prot_dom_level_eval_ELM = HV_prot_dom_int_eval(bench_data_type = "ELM", all_pfam_dom = HV_all_pfam_doms)
HV_prot_dom_level_eval_PRMdb = HV_prot_dom_int_eval(bench_data = prm_pfam_final,bench_data_type = "PRMdb", all_pfam_dom = HV_all_pfam_doms)


saveRDS(list("ELM" = list("Prot_level" = HV_prot_level_eval_ELM, "HH_motif" = HV_HH_level_eval_ELM, "Prot_dom" = HV_prot_dom_level_eval_ELM),
             "PRMdb" = list("Prot_level" = HV_prot_level_eval_PRMdb, "HH_motif" = HV_HH_level_eval_PRMdb, "Prot_dom" = HV_prot_dom_level_eval_PRMdb)),
        "data/output/HV_Complete_Evaluation.rds")

beep(3)

# Human only evaluation ---------------------------------------------------

human_prot_level_eval_ELM = HV_prot_level_eval(slim_hits = human_final_hits, bench_data_type = "ELM", all_pfam_dom = human_all_pfam_doms, datasets_input = human_qslim_datasets)
human_prot_level_eval_PRMdb = HV_prot_level_eval(slim_hits = human_final_hits,bench_data = pwm_prmdb_scan ,bench_data_type = "PRMdb", all_pfam_dom = human_all_pfam_doms, datasets_input = human_qslim_datasets)

human_HH_level_eval_ELM = HV_motif_level_eval(slim_hits = human_final_hits,bench_data_type = "ELM", all_pfam_dom = human_all_pfam_doms, datasets_input = human_qslim_datasets, prot_seqs = all_seqs)
human_HH_level_eval_PRMdb = HV_motif_level_eval(slim_hits = human_final_hits, bench_data = pwm_prmdb_scan, bench_data_type = "PRMdb", all_pfam_dom = human_all_pfam_doms, datasets_input = human_qslim_datasets, prot_seqs = all_seqs)

human_prot_dom_level_eval_ELM = HV_prot_dom_int_eval(bench_data_type = "ELM", all_pfam_dom = human_all_pfam_doms)
human_prot_dom_level_eval_PRMdb = HV_prot_dom_int_eval(bench_data = prm_pfam_final,bench_data_type = "PRMdb", all_pfam_dom = human_all_pfam_doms)


saveRDS(list("ELM" = list("Prot_level" = human_prot_level_eval_ELM, "HH_motif" = human_HH_level_eval_ELM, "Prot_dom" = human_prot_dom_level_eval_ELM),
             "PRMdb" = list("Prot_level" = human_prot_level_eval_PRMdb, "HH_motif" = human_HH_level_eval_PRMdb, "Prot_dom" = human_prot_dom_level_eval_PRMdb)),
        "data/output/human_Complete_Evaluation.rds")



# Shared inst Human only vs Host-viral ------------------------------------
host_viral_shared_prot = HV_final_hits[which(HV_final_hits$uniprot %in% human_final_hits$uniprot),]

human_only_common_insts = human_final_hits[which(human_final_hits$uniprot %in% HV_final_hits$uniprot),]

req_com_host_viral = list("new_hits" = host_viral_shared_prot, "new_all_pfam_doms" = HV_all_pfam_doms)
req_com_human = list("new_hits" = human_only_common_insts, "new_all_pfam_doms" = human_all_pfam_doms, "datasets_input" = human_qslim_datasets)
Evaluation_automate = function(required,host_viral = T){
  if (host_viral){
    HV_prot_level_eval_ELM = HV_prot_level_eval(slim_hits = required$new_hits, bench_data_type = "ELM", all_pfam_dom = HV_all_pfam_doms, datasets_input = HV_qslim_datasets)
    HV_prot_level_eval_PRMdb = HV_prot_level_eval(slim_hits = required$new_hits,bench_data = pwm_prmdb_scan ,bench_data_type = "PRMdb", all_pfam_dom = HV_all_pfam_doms, datasets_input = HV_qslim_datasets)

    HV_HH_level_eval_ELM = HV_motif_level_eval(slim_hits = required$new_hits,bench_data_type = "ELM", all_pfam_dom = HV_all_pfam_doms, datasets_input = HV_qslim_datasets, prot_seqs = all_seqs)
    HV_HH_level_eval_PRMdb = HV_motif_level_eval(slim_hits = required$new_hits, bench_data = pwm_prmdb_scan, bench_data_type = "PRMdb", all_pfam_dom = HV_all_pfam_doms, datasets_input = HV_qslim_datasets, prot_seqs = all_seqs)

    return(list("ELM" = list("Prot_level" = HV_prot_level_eval_ELM, "HH_motif" = HV_HH_level_eval_ELM),
                "PRMdb" = list("Prot_level" = HV_prot_level_eval_PRMdb, "HH_motif" = HV_HH_level_eval_PRMdb)))
  }
  else{
    human_prot_level_eval_ELM = HV_prot_level_eval(slim_hits = required$new_hits, bench_data_type = "ELM", all_pfam_dom = human_all_pfam_doms, datasets_input = human_qslim_datasets)
    human_prot_level_eval_PRMdb = HV_prot_level_eval(slim_hits = required$new_hits,bench_data = pwm_prmdb_scan ,bench_data_type = "PRMdb", all_pfam_dom = human_all_pfam_doms, datasets_input = human_qslim_datasets)

    human_HH_level_eval_ELM = HV_motif_level_eval(slim_hits = required$new_hits,bench_data_type = "ELM", all_pfam_dom = human_all_pfam_doms, datasets_input = human_qslim_datasets, prot_seqs = all_seqs)
    human_HH_level_eval_PRMdb = HV_motif_level_eval(slim_hits = required$new_hits, bench_data = pwm_prmdb_scan, bench_data_type = "PRMdb", all_pfam_dom = human_all_pfam_doms, datasets_input = human_qslim_datasets, prot_seqs = all_seqs)

    return(list("ELM" = list("Prot_level" = human_prot_level_eval_ELM, "HH_motif" = human_HH_level_eval_ELM),
                "PRMdb" = list("Prot_level" = human_prot_level_eval_PRMdb, "HH_motif" = human_HH_level_eval_PRMdb)))
  }
}

shared_host_viral_eval = Evaluation_automate(required = req_com_host_viral, host_viral = T)
saveRDS(shared_host_viral_eval,"data/output/shared_host_viral_eval.rds")
shared_human_eval = Evaluation_automate(required = req_com_human, host_viral = F)
saveRDS(shared_human_eval,"data/output/shared_human_eval.rds")

#Comparing human and host-viral with comparable number of motif instances

host_viral_shared_prot = HV_final_hits[which(HV_final_hits$uniprot %in% human_final_hits$uniprot),]
human_only_common_insts = human_final_hits[which(human_final_hits$uniprot %in% HV_final_hits$uniprot),]

x = count_motif_instances(human_only_common_insts)

sig_levels = seq(0.001,0.1,0.001)
sig_cutoffs = list()
for (i in 1:length(sig_levels)){
  host_viral_shared_prot = HV_final_hits[which(HV_final_hits$uniprot %in% human_final_hits$uniprot),]
  host_viral_shared_prot = host_viral_shared_prot[which(host_viral_shared_prot$Sig < sig_levels[i]),]
  human_only_common_insts = human_final_hits[which(human_final_hits$uniprot %in% host_viral_shared_prot$uniprot),]
  sig_cutoffs[[i]] = data.frame(sig_thresh = sig_levels[i], HV_insts_n = count_motif_instances(host_viral_shared_prot),
                                human_insts_n = count_motif_instances(human_only_common_insts))
}
sig_cutoffs = dplyr::bind_rows(sig_cutoffs)
sig_cutoffs$HV_human_diff = sig_cutoffs$HV_insts_n - sig_cutoffs$human_insts_n

optim_sig = sig_cutoffs$sig_thresh[which.min(abs(sig_cutoffs$HV_human_diff))]

host_viral_shared_prot = HV_final_hits[which(HV_final_hits$uniprot %in% human_final_hits$uniprot),]
host_viral_shared_prot = host_viral_shared_prot[which(host_viral_shared_prot$Sig < optim_sig),]
human_only_common_insts = human_final_hits[which(human_final_hits$uniprot %in% host_viral_shared_prot$uniprot),]

req_com_host_viral = list("new_hits" = host_viral_shared_prot, "new_all_pfam_doms" = HV_all_pfam_doms)
req_com_human = list("new_hits" = human_only_common_insts, "new_all_pfam_doms" = human_all_pfam_doms, "datasets_input" = human_qslim_datasets)


shared_host_viral_eval = Evaluation_automate(required = req_com_host_viral, host_viral = T)
saveRDS(shared_host_viral_eval,"data/output/shared_host_viral_eval_comparable_insts.rds")
shared_human_eval = Evaluation_automate(required = req_com_human, host_viral = F)
saveRDS(shared_human_eval,"data/output/shared_human_eval_comparable_insts.rds")
