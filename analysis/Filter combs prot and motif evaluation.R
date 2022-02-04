library(HVSlimPred)
HV_clinvar_path = readRDS("data/output/HV_Clinvar_and_PTM.rds")
HV_clinvar_path = HV_clinvar_path$ClinVar_path$ClinVar_path
HV_final_hits = readRDS("data/output/HV_final_hits.rds")
HV_HMM_iELM = readRDS("data/output/HV_HMM_iELM.rds")
pfam_doms = readRDS("data/input/pfam_doms_human.rds")
HV_Enriched_doms = readRDS("data/output/HV_enriched_doms_with_emppval.rds")
HV_all_pepsite = readRDS("data/output/HV_pepsite_final_with_AF.rds")
HV_input_prots = readRDS("data/input/HV_proteins_per_dataset_allhitsint.RDS")

HMM_insts = HMM_iELM_classify(new_hits_ints = HV_final_hits, HMM_iELM = HV_HMM_iELM,host_viral = T)
HMM_insts = HMM_insts[which(HMM_insts$passed_iELM == 1),]

pwm_prmdb_scan = read.csv("data/input/PRMDB/pwm_prmdb_scan.csv.gz")
pwm_ids = unique(pwm_prmdb_scan$Id)
domain_ids_pwm = list()
for (i in 1:length(pwm_ids)){
  if (nchar(pwm_ids[i]) == 8){
    domain_ids_pwm[[i]] = data.frame(Id = pwm_ids[i], domain_id = pwm_ids[i])
  }
  else{
    domain_ids_pwm[[i]] = data.frame(Id = pwm_ids[i],
                                     domain_id = paste(strsplit(pwm_ids[i], "_")[[1]][c(1,2)], collapse = "_"))
  }
}
domain_ids_pwm = dplyr::bind_rows(domain_ids_pwm)
pwm_prmdb_scan = dplyr::left_join(pwm_prmdb_scan, domain_ids_pwm)

PRM_master = read.csv("data/input/PRMDB/PRM_Master.csv")
PRM_master = PRM_master[which(PRM_master$Species == "Human"),]
PRM_master = PRM_master[,c(1:7)]
colnames(PRM_master)[1] = "domain_id"
pwm_prmdb_scan = pwm_prmdb_scan[which(pwm_prmdb_scan$domain_id %in% PRM_master$domain_id),]
names(PRM_master) = c("domain_id", "Domain_protein_name", "Domain_class", "Domain_class_number", "Domain_uniprot","Domain_start", "Domain_end")
PRMdb_data = dplyr::left_join(pwm_prmdb_scan, PRM_master)

all_seqs = read.csv("data/input/uniprot_sequences.csv.gz")


all_filters_data = list()
for (i in c(TRUE, FALSE)){
  if (i){
    all_filters_data[["conserved_only"]] = prepare_required_data_OR(EM = HV_Enriched_doms,
                                                                    clin_path = HV_clinvar_path,
                                                                    iELM = HMM_insts,
                                                                    new_hits_ints = HV_final_hits,
                                                                    pepsite_res = HV_all_pepsite,min_domain = 5,
                                                                    adj.pval = 0.05,slim_input_prots = HV_input_prots,cons = i, prot_only = F)
  }
  else{
    all_filters_data[["all_hits"]] = prepare_required_data_OR(EM = HV_Enriched_doms,
                                                                    clin_path = HV_clinvar_path,
                                                                    iELM = HMM_insts,
                                                                    new_hits_ints = HV_final_hits,
                                                                    pepsite_res = HV_all_pepsite,min_domain = 5,
                                                                    adj.pval = 0.05,slim_input_prots = HV_input_prots,cons = i, prot_only = F)
  }
}

saveRDS(all_filters_data, "data/output/HV_all_filters_data.rds")


run_evaluation_per_filter = function(EM_list, Clinvar_path, iELM_insts, new_hits_ints, pepsite_res,prot_seqs,
                                     min_dom = 5, adj_pval = 0.05,qslim_input = NULL,
                                     bench_type = c("ELM","PRMDB"), cons_only = F, eval_type = c("prot","motif")){

  if (cons_only){
    new_hits_ints = new_hits_ints[!is.na(new_hits_ints$Cons),]
  }

  if (eval_type == "prot"){
    if (bench_type == "ELM"){
      prot_comparison = All_prot_OR(EM_list = EM_list, Clinvar_path = Clinvar_path, iELM_insts = iELM_insts,
                                    new_hits_ints = new_hits_ints,pepsite_res = pepsite_res,bench_data = NULL,
                                    min_dom = 5,adj_pval = 0.05, qslim_input = qslim_input,
                                    comparison = T, bg_possible = T, prot_only = T)

      prot_no_comparison = All_prot_OR(EM_list = EM_list, Clinvar_path = Clinvar_path, iELM_insts = iELM_insts,
                                       new_hits_ints = new_hits_ints,pepsite_res = pepsite_res,bench_data = NULL,
                                       min_dom = 5,adj_pval = 0.05, qslim_input = qslim_input,
                                       comparison = F, bg_possible = T, prot_only = T)

      final = list("comparison_prot" = prot_comparison, "no_comparison_prot" = prot_no_comparison)
    }
    else{
      prot_comparison = All_prot_OR(EM_list = EM_list, Clinvar_path = Clinvar_path, iELM_insts = iELM_insts,
                                    new_hits_ints = new_hits_ints,pepsite_res = pepsite_res,bench_data = PRMdb_data,
                                    min_dom = 5,adj_pval = 0.05, qslim_input = qslim_input,
                                    comparison = T, bg_possible = T, prot_only = T)

      prot_no_comparison = All_prot_OR(EM_list = EM_list, Clinvar_path = Clinvar_path, iELM_insts = iELM_insts,
                                       new_hits_ints = new_hits_ints,pepsite_res = pepsite_res,bench_data = PRMdb_data,
                                       min_dom = 5,adj_pval = 0.05, qslim_input = qslim_input,
                                       comparison = F, bg_possible = T, prot_only = T)

      final = list("comparison_prot" = prot_comparison, "no_comparison_prot" = prot_no_comparison)
    }
  }
  else{
    if (bench_type == "ELM"){
      final = All_motif_inst_OR(EM_list = EM_list, Clinvar_path = Clinvar_path, iELM_insts = iELM_insts,
                          new_hits_ints = new_hits_ints,pepsite_res = pepsite_res,bench_data = NULL,
                          min_dom = 5,adj_pval = 0.05, datasets_input = qslim_input,
                          prot_only = F, prot_seqs = prot_seqs)
    }
    else{
      final = All_motif_inst_OR(EM_list = EM_list, Clinvar_path = Clinvar_path, iELM_insts = iELM_insts,
                          new_hits_ints = new_hits_ints,pepsite_res = pepsite_res,bench_data = PRMdb_data,
                          min_dom = 5,adj_pval = 0.05, datasets_input = qslim_input,
                          prot_only = F, prot_seqs = prot_seqs)
    }
  }
  return(final)
}

counter = 0
pb = txtProgressBar(min = 0, max = 8, style = 3)
OR_all_eval = list()
for (i in c("prot", "motif")){
  for (j in c("ELM","PRMDB")){
    for (k in c(TRUE, FALSE)){
      counter = counter + 1

      res = run_evaluation_per_filter(EM_list = HV_Enriched_doms,
                                Clinvar_path = HV_clinvar_path,
                                iELM_insts = HMM_insts,
                                new_hits_ints = HV_final_hits,
                                pepsite_res = HV_all_pepsite,prot_seqs = all_seqs,min_dom = 5,
                                adj_pval = 0.05, qslim_input = HV_input_prots,
                                bench_type = j, cons_only = k, eval_type = i)

      if (k){
        OR_all_eval[[i]][[j]][["conserved_only"]] = res
      }
      else{
        OR_all_eval[[i]][[j]][["all_hits"]] = res
      }
      setTxtProgressBar(pb, counter)
    }
  }
}

saveRDS(OR_all_eval, "data/output/HV_all_filters_eval_OR.rds")
