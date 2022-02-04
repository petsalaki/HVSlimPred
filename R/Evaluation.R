#' Protein-level evaluation of predicted hits
#'
#' Calculates performance metrics for evaluating predicted hits against ELM true positive instances on the protein-level.
#'
#' @details
#' For protein-level enrichment, we measured the enrichment of true-positives in our predicted dataset using a one-tailed fisher-exact test, where the odds ratio represents the magnitude of the enrichment. True-positives are the number of motif-carrying proteins present in both the predicted dataset and the ELM dataset regardless of whether the predicted protein has the right motif or found in the right location
#'
#' @author Bishoy Wadie, Evangelia Petsalaki
#' @return
#' A data frame containing all the relevant protein-level performance metrics for each domain enrichment filter in addition to the non-filtered qslim output.
#' @export
#' @importFrom utils read.delim read.table read.csv
#' @examples
#' prot_eval_metrics = HV_prot_level_eval()
HV_prot_level_eval = function(slim_hits,bench_data, bench_data_type = c("ELM", "PRMdb"), datasets_input = NULL, all_pfam_dom = NULL){
  all_hits = slim_hits
  # if (parse_occ & is.null(occ_out)){
  #   all_hits = parse_qslim(main_output = main_output,occ_out = occ_out, compari_output = compari_out, preprocess_compari = preprocess_compari)
  # }
  # else{
  #   all_hits = preprocess_slimfinder_out(main_output = main_output, occ_out = occ_out)
  #   all_hits = all_hits[[1]]
  # }
  if (is.null(datasets_input)){
    # TODO change input path
    qslim_input = download_data(path = "Evaluation_data/proteins_per_dataset_allhitsint.RDS", is.rds = T)
  }
  else{
    qslim_input = datasets_input
  }
  if (bench_data_type == "ELM" & missing(bench_data)){
    # FIXME fix getting ELM instances after changing the function GET_ELM_all
    ELM_all_instances = Get_ELM_all(ELM_data_type = "instances")
  }
  else{
    ELM_all_instances = bench_data
  }

  All_enrichment = list()
  for (i in 1:nrow(pfam_dom_enrich_info)){
    filtered = pfam_dom_prot_eval(pfam = pfam_dom_enrich_info$pfam_filter[i], ELM_data = ELM_all_instances,
                                  pred_qslim = all_hits, qslim_input_dat = qslim_input, all_pfam_dom_list = all_pfam_dom)
    All_enrichment[[i]] = filtered
  }

  All_enrichment = dplyr::bind_rows(All_enrichment)

  qslim_enrich = HV_ProtEval_metrics(is_pfam_dom = F, ELM_all = ELM_all_instances,
                                     pred_hits = all_hits, qslim_input = qslim_input)
  qslim_enrich$filter = "qslim"

  All_enrichment = rbind.data.frame(qslim_enrich, All_enrichment)

  All_enrichment$Accuracy = (All_enrichment$TP + All_enrichment$TN) / (All_enrichment$TP + All_enrichment$FN +
                                                                         All_enrichment$FP + All_enrichment$TN)

  for (i in 1:nrow(All_enrichment)){
    All_enrichment$MCC[i] = mltools::mcc(TP = All_enrichment$TP[i], FP = All_enrichment$FP[i],
                                       TN = All_enrichment$TN[i], FN = All_enrichment$FN[i])
  }
  return(All_enrichment)
}

fisher_and_F1 = function(numbers, pval_method = c("greater", "two.sided", "less")){
  ficher_mat = matrix(c(numbers[1], numbers[2], numbers[3], numbers[4]),
                      nrow = 2, ncol = 2, byrow = T)
  colnames(ficher_mat) = c("predicted", "not_predicted")
  rownames(ficher_mat) = c("ELM", "not_ELM")

  test = hypergea::hypergeom.test(ficher_mat, alternative = pval_method)

  results = data.frame(odds_ratio = as.numeric(test$estimate),p.val = test$p.value,TP = ficher_mat[1,1], FN = ficher_mat[1,2],
                       FP = ficher_mat[2,1],TN = ficher_mat[2,2])
  results$precision = results$TP / (results$TP + results$FP)
  results$recall = results$TP / (results$TP + results$FN)
  results$F1 = 2*(results$precision * results$recall) / (results$precision + results$recall)
  # TODO Add F0.5 calculation here as well
  results$log_OR = log2(results$odds_ratio)
  return(results)
}
HV_ProtEval_metrics = function(df, is_pfam_dom = T, ELM_all, pred_hits, qslim_input){
  ELM_tp = ELM_all[which(ELM_all$Logic == "true positive"),]
  if (!is_pfam_dom & missing(df)){
    filtered = pred_hits
  }
  else{
    filtered = pred_hits %>% dplyr::filter(uniprot %in% df$uniprot)
  }

  prots = qslim_input[which(names(qslim_input) %in% filtered$Dataset)]
  prots = unique(unlist(prots))
  possible_ELM_prots = dplyr::intersect(ELM_tp$uniprot, prots)
  possible_ELM = ELM_tp %>% dplyr::filter(uniprot %in% prots)
  TP = dplyr::intersect(filtered$uniprot, ELM_tp$uniprot)
  FP = unique(filtered$uniprot[which(filtered$uniprot %nin% TP)])
  FN = possible_ELM_prots[which(possible_ELM_prots %nin% TP)]
  TN = prots[which(prots %nin% c(ELM_tp$uniprot, filtered$uniprot))]
  numbers = c(length(TP), length(unique(FN)), length(FP), length(unique(TN)))
  return(fisher_and_F1(numbers, pval_method = "greater"))
}

clean_pfam_enrich = function(uniprot, domains){
  .Deprecated("pfam_dom_prot_eval")
  if (!is.na(domains)){
    dom_list = strsplit(domains, ";")[[1]]
    domain_table = data.frame(uniprot_id = uniprot, domain = dom_list)
    domain_table = domain_table %>% tidyr::separate(domain, into = c("PFAM","Int_freq", "bg_freq",
                                                              "enrichment_val", "pval", "adj_pval",
                                                              "emp_pval"), sep = ":")
    colnames(domain_table)[1] = "uniprot"
    return(domain_table)
  }
  else{
    domain_table = data.frame(matrix(nrow = 1, ncol = 8))
    colnames(domain_table) = c("uniprot", "PFAM","Int_freq", "bg_freq","enrichment_val",
                               "pval", "adj_pval","emp_pval")
    domain_table$uniprot = uniprot
    return(domain_table)
  }
}
clean_pfam_dom_data = function(df){
  .Deprecated("pfam_dom_prot_eval")
  pfam_clean = list()
  for (i in 1:nrow(df)){
    result = clean_pfam_enrich(df$uniprot[i], df$domains[i])
    result$total_domains_nr = df$total_domains_nr[i]
    pfam_clean[[i]] = result
  }
  pfam_clean = dplyr::bind_rows(pfam_clean)
  pfam_clean = pfam_clean[!is.na(pfam_clean$PFAM),]
  return(pfam_clean)
}
pfam_dom_prot_eval = function(pfam, ELM_data, pred_qslim, qslim_input_dat, all_pfam_dom_list){
  if (is.null(all_pfam_dom_list)){
    # TODO change input path
    all_pfam_dom_enrich = download_data(path = "Evaluation_data/All_pfam_dom_enrich.RDS", is.rds = T)
  }
  else{
    all_pfam_dom_enrich = all_pfam_dom_list
  }
  #file_name = pfam_dom_enrich_info[which(pfam_dom_enrich_info$pfam_filter == pfam),]$file_name
  cleaned = all_pfam_dom_enrich[[pfam]]

  enrichment = HV_ProtEval_metrics(cleaned, ELM_all = ELM_data, pred_hits = pred_qslim,
                                   qslim_input = qslim_input_dat)
  enrichment$filter = pfam
  return(enrichment)
}

#' Motif-level evaluation of predicted hits
#'
#' Calculates performance metrics for evaluating predicted hits against ELM true positive instances on the motif-level.
#'
#' @details
#' For motif-level enrichment, we simply cannot use a binary classification as we did for protein-level evaluation because in reality, predicted motifs are partially correct to some extent as they might contain true-positive residues in a given sequence stretch, and therefore we used a re-implemented version of the evaluation protocol proposed in [Prytuliak et al. 2017](https://academic.oup.com/nar/article/45/W1/W470/3782606) instead of binary classification, where we computed the common performance metrics (Recall, precision F1, etc .. ) both residue-wise and site-wise given that the motif-carrying proteins are also found in the ELM benchmarking dataset. So this analysis was not performed on proteins not reported in the ELM dataset.
#' @return
#' A data frame containing all the relevant motif-level performance metrics for each domain enrichment filter in addition to the non-filtered qslim output.
#' @export
#' @author Bishoy Wadie, Evangelia Petsalaki
#' @importFrom utils read.delim read.table read.csv
#' @references
#' Prytuliak, Roman, et al. "HH-MOTiF: de novo detection of short linear motifs in proteins by Hidden Markov Model comparisons." \emph{Nucleic acids research} 45.W1 (2017): W470-W477.
#' @examples
#' motif_eval_metrics = HV_motif_level_eval()
HV_motif_level_eval = function(slim_hits, bench_data, bench_data_type = c("ELM", "PRMdb"), datasets_input = NULL, all_pfam_dom = NULL, prot_seqs = NULL){
  all_hits = slim_hits
  # if (parse_occ & is.null(occ_out)){
  #   all_hits = parse_qslim(main_output = main_output,occ_out = occ_out, compari_output = compari_out, preprocess_compari = preprocess_compari)
  # }
  # else{
  #   all_hits = preprocess_slimfinder_out(main_output = main_output, occ_out = occ_out)
  #   all_hits = all_hits[[1]]
  # }
  if (is.null(datasets_input)){
    # TODO change input path
    qslim_input = download_data(path = "Evaluation_data/proteins_per_dataset_allhitsint.RDS", is.rds = T)
  }
  else{
    qslim_input = datasets_input
  }

  if (bench_data_type == "ELM" & missing(bench_data)){
    # FIXME fix getting ELM instances after changing the function GET_ELM_all
    ELM_all_instances = Get_ELM_all(ELM_data_type = "instances")
  }
  else{
    ELM_all_instances = bench_data
  }
  if (is.null(prot_seqs)){
    # TODO change input path
    uniprot_seqs = read.csv(download_data(path = "Evaluation_data/uniprot_sequences.csv.gz"), header = T)
  }
  else{
    uniprot_seqs = prot_seqs
  }
  # uniprot_seqs[which(stringr::str_detect(uniprot_seqs$uniprot, "-") == TRUE),]$uniprot =
  #   substr(uniprot_seqs[which(stringr::str_detect(uniprot_seqs$uniprot, "-") == TRUE),]$uniprot, 1, nchar(uniprot_seqs[which(stringr::str_detect(uniprot_seqs$uniprot, "-") == TRUE),]$uniprot) -2 )

  required = list("pred_hits" = all_hits, "qslim_input" = qslim_input, "ELM_all" = ELM_all_instances,
                  "uniprot_seqs"  = uniprot_seqs)
  all_input_prots = unique(unlist(qslim_input))

  ELM_possible = ELM_all_instances[which(ELM_all_instances$uniprot %in% all_input_prots),]
  ELM_possible = ELM_possible[which(ELM_possible$Logic == "true positive"),]
  qslim_instances = dplyr::select(all_hits, uniprot, Start_Pos, End_Pos)
  qslim_instances = qslim_instances[!duplicated(qslim_instances),]

  possible_motifs = unique(ELM_possible$Id)

  all_metrics_per_motif = list()
  weights_per_motif = list()
  for (i in 1:length(possible_motifs)){
    result = calc_metrics_per_motif(qslim_instances, possible_motifs[i], min_annot = 3, required_data = required)
    all_metrics_per_motif[[i]] = result[["metrics"]]
    weights_per_motif[[i]] = result[["weights"]]
  }
  all_metrics_per_motif = dplyr::bind_rows(all_metrics_per_motif)
  weights_per_motif = dplyr::bind_rows(weights_per_motif)
  qslim_F1 = calc_F1_score_HH(weights = weights_per_motif, metrics = all_metrics_per_motif)
  qslim_F1$filter = "qslim"

  All_F1_scores = list()
  for (i in 1:nrow(pfam_dom_enrich_info)){
    filtered = pfam_dom_motif_eval(pfam_dom_enrich_info$pfam_filter[i], prereqs = required, all_pfam_dom_list = all_pfam_dom)
    All_F1_scores[[i]] = filtered
  }
  All_F1_scores = dplyr::bind_rows(All_F1_scores)
  All_F1_scores = rbind.data.frame(qslim_F1, All_F1_scores)

  All_F1_scores$weight_site_resid = paste0(All_F1_scores$weight_type, " ", All_F1_scores$site_resid)
  return(list("All_F1_scores" = All_F1_scores, "all_metrics_per_motif" = all_metrics_per_motif,
              "weights_per_motif" = weights_per_motif))
}

calc_metrics_per_motif = function(df,motif, min_annot = 1, required_data, match_insts = F){
  ELM_all_instances = required_data[["ELM_all"]]
  qslim_input = required_data[["qslim_input"]]
  all_hits = required_data[["pred_hits"]]
  uniprot_seqs = required_data[["uniprot_seqs"]]

  all_input_prots = unique(unlist(qslim_input))
  ELM_possible = ELM_all_instances[which(ELM_all_instances$Logic == "true positive"),]
  ELM_possible = ELM_possible[which(ELM_possible$uniprot %in% all_input_prots),]
  qslim_instances = dplyr::select(all_hits, uniprot, Start_Pos, End_Pos)
  qslim_instances = qslim_instances[!duplicated(qslim_instances),]

  motif_filtered = ELM_possible[which(ELM_possible$Id == motif),]
  if (length(unique(motif_filtered$uniprot)) < min_annot){
    return(NULL)
  }
  metrics = rep(0,10)
  names(metrics) = c("TP_res", "PA_res","PNA_res","FP_res","FN_res","TN_res", "AP_site",
                     "PA_site","PNA_site","ANP_site")

  if (!match_insts){
    df_instances = qslim_instances[which(qslim_instances$uniprot %in% df$uniprot),]
  }
  else{
    df_instances = plyr::match_df(qslim_instances, df, on = c("uniprot", "Start_Pos", "End_Pos"))
  }
  qslim_filtered = df_instances[which(df_instances$uniprot %in% motif_filtered$uniprot),]
  motif_filtered = motif_filtered[which(motif_filtered$uniprot %in% qslim_filtered$uniprot),]

  if(nrow(qslim_filtered) == 0){
    return(NULL)
  }
  unique_prots = unique(qslim_filtered$uniprot)
  count_inst = list()
  for (i in 1:length(unique_prots)){
    filtered_prots = qslim_filtered[which(qslim_filtered$uniprot == unique_prots[i]),]
    per_prot = list()
    for (j in 1:nrow(filtered_prots)){
      per_prot[[j]] = c(filtered_prots$Start_Pos[j]:filtered_prots$End_Pos[j])
    }
    count_inst[[unique_prots[i]]] = unlist(per_prot)
  }

  ELM_prots = unique(motif_filtered$uniprot)
  ELM_inst = list()
  for (i in 1:length(ELM_prots)){
    filtered_prots = motif_filtered[which(motif_filtered$uniprot == ELM_prots[i]),]
    per_prot = list()
    for (j in 1:nrow(filtered_prots)){
      per_prot[[j]] = c(filtered_prots$start[j]:filtered_prots$end[j])
    }
    ELM_inst[[ELM_prots[i]]] = unlist(per_prot)
  }


  PA_prot = dplyr::intersect(motif_filtered$uniprot, qslim_filtered$uniprot)
  Annot_inst = ELM_inst[which(names(ELM_inst) %in% PA_prot)]
  for (i in 1:length(PA_prot)){
    res_annot = Annot_inst[which(names(Annot_inst) == PA_prot[i])][[1]]
    res_pred = count_inst[which(names(count_inst) == PA_prot[i])][[1]]
    metrics["TP_res"] = metrics["TP_res"] + length(unique(dplyr::intersect(res_annot, res_pred)))
    metrics["PA_res"] = metrics["PA_res"] + length(which(res_pred %in% res_annot))
    metrics["PNA_res"] = metrics["PNA_res"] + length(which(res_pred %nin% res_annot))
    metrics["FP_res"] = metrics["FP_res"] + length(unique(res_pred[which(res_pred %nin% res_annot)]))
    metrics["FN_res"] = metrics["FN_res"] + length(unique(res_annot[which(res_annot %nin% res_pred)]))
    metrics["TN_res"] = metrics["TN_res"] +
      (nchar(uniprot_seqs[which(uniprot_seqs$uniprot == PA_prot[i]),]$sequence) - length(unique(c(res_annot, res_pred))))
  }
  for (i in 1:length(PA_prot)){
    filtered_annot = motif_filtered[which(motif_filtered$uniprot == PA_prot[i]),]
    filtered_pred = qslim_filtered[which(qslim_filtered$uniprot == PA_prot[i]),]
    sites_annot = ELM_inst[which(names(ELM_inst) == PA_prot[i])]
    sites_pred = count_inst[which(names(count_inst) == PA_prot[i])]

    for (j in 1:nrow(filtered_pred)){
      resids_pred = c(filtered_pred$Start_Pos[j]:filtered_pred$End_Pos[j])
      if(length(dplyr::intersect(resids_pred,sites_annot[[1]])) >= 1){
        metrics["PA_site"] = metrics["PA_site"] + 1
      }
      else{
        metrics["PNA_site"] = metrics["PNA_site"] + 1
      }
    }
    for (k in 1:nrow(filtered_annot)){
      resids_annot = c(filtered_annot$start[k]:filtered_annot$end[k])
      if(length(dplyr::intersect(resids_annot,sites_pred[[1]])) >= 1){
        metrics["AP_site"] = metrics["AP_site"] + 1
      }
      else{
        metrics["ANP_site"] = metrics["ANP_site"] + 1
      }
    }
  }

  final_result = data.frame(Id = motif, Rc_res = NA, Rc_site = NA, Pr_res = NA,
                            Pr_site = NA, SP_res = NA, FPR_res = NA, PC_res = NA,
                            PC_site = NA)

  final_result$Rc_res = metrics["TP_res"] / (metrics["TP_res"] + metrics["FN_res"])
  final_result$Rc_site = metrics["AP_site"] / (metrics["AP_site"] + metrics["ANP_site"])
  final_result$Pr_res = metrics["PA_res"] / (metrics["PA_res"] + metrics["PNA_res"])
  final_result$Pr_site = metrics["PA_site"] / (metrics["PA_site"] + metrics["PNA_site"])
  final_result$SP_res = metrics["TN_res"] / (metrics["TN_res"] + metrics["FP_res"])
  final_result$FPR_res = metrics["FP_res"] / (metrics["TN_res"] + metrics["FP_res"])
  final_result$PC_res = metrics["PA_res"] / (metrics["PA_res"] + metrics["PNA_res"] + metrics["FN_res"])
  final_result$PC_site = metrics["PA_site"] / (metrics["PA_site"] + metrics["PNA_site"] + metrics["ANP_site"])

  weights = data.frame(Id = motif, No_proteins = NA, No_sites = NA, No_resids = NA)
  weights$No_proteins = length(unique(motif_filtered$uniprot))
  weights$No_sites = nrow(motif_filtered)
  counter = 0
  for (i in 1:length(ELM_inst)){
    counter = counter + length(unique(unlist(ELM_inst[[i]])))
  }
  weights$No_resids = counter

  return(list("metrics" = final_result, "weights" = weights))
}
weight_values  = function(metrics, weights, weight_type, site_resid){
  final_result = data.frame(Rc = NA, Pr = NA, PC = NA, SP = NA, weight_type = weight_type,
                            site_resid = site_resid)
  if(weight_type == "None"){
    if(site_resid == "residue_wise"){
      final_result$Rc = mean(metrics$Rc_res)
      final_result$Pr = mean(metrics$Pr_res)
      final_result$PC = mean(metrics$PC_res)
      final_result$SP = mean(metrics$SP_res)
    }
    else{
      final_result$Rc = mean(metrics$Rc_site)
      final_result$Pr = mean(metrics$Pr_site)
      final_result$PC = mean(metrics$PC_site)
    }
  }
  else if (weight_type == "prot"){
    if(site_resid == "residue_wise"){
      final_result$Rc = sum(weights$No_proteins * metrics$Rc_res) / sum(weights$No_proteins)
      final_result$Pr = sum(weights$No_proteins * metrics$Pr_res) / sum(weights$No_proteins)
      final_result$PC = sum(weights$No_proteins * metrics$PC_res) / sum(weights$No_proteins)
      final_result$SP = sum(weights$No_proteins * metrics$SP_res) / sum(weights$No_proteins)
    }
    else{
      final_result$Rc = sum(weights$No_proteins * metrics$Rc_site) / sum(weights$No_proteins)
      final_result$Pr = sum(weights$No_proteins * metrics$Pr_site) / sum(weights$No_proteins)
      final_result$PC = sum(weights$No_proteins * metrics$PC_site) / sum(weights$No_proteins)
    }
  }
  else if(weight_type == "site"){
    if(site_resid == "residue_wise"){
      final_result$Rc = sum(weights$No_sites * metrics$Rc_res) / sum(weights$No_sites)
      final_result$Pr = sum(weights$No_sites * metrics$Pr_res) / sum(weights$No_sites)
      final_result$PC = sum(weights$No_sites * metrics$PC_res) / sum(weights$No_sites)
      final_result$SP = sum(weights$No_sites * metrics$SP_res) / sum(weights$No_sites)
    }
    else{
      final_result$Rc = sum(weights$No_sites * metrics$Rc_site) / sum(weights$No_sites)
      final_result$Pr = sum(weights$No_sites * metrics$Pr_site) / sum(weights$No_sites)
      final_result$PC = sum(weights$No_sites * metrics$PC_site) / sum(weights$No_sites)
    }
  }
  else if (weight_type == "residue"){
    if(site_resid == "residue_wise"){
      final_result$Rc = sum(weights$No_resids * metrics$Rc_res) / sum(weights$No_resids)
      final_result$Pr = sum(weights$No_resids * metrics$Pr_res) / sum(weights$No_resids)
      final_result$PC = sum(weights$No_resids * metrics$PC_res) / sum(weights$No_resids)
      final_result$SP = sum(weights$No_resids * metrics$SP_res) / sum(weights$No_resids)
    }
    else{
      final_result$Rc = sum(weights$No_resids * metrics$Rc_site) / sum(weights$No_resids)
      final_result$Pr = sum(weights$No_resids * metrics$Pr_site) / sum(weights$No_resids)
      final_result$PC = sum(weights$No_resids * metrics$PC_site) / sum(weights$No_resids)
    }
  }
  return(final_result)
}
calc_F1_score_HH = function(metrics, weights){
  results = list()
  counter = 0
  for (i in c("None", "prot", "site", "residue")){
    for (j in c("site_wise", "residue_wise")){
      counter = counter + 1
      results[[counter]] = weight_values(weights = weights, metrics = metrics, weight_type = i, site_resid = j)
    }
  }
  results = dplyr::bind_rows(results)
  results$F1 = (2 * results$Pr * results$Rc) / (results$Pr + results$Rc)
  # TODO Add calculation of F0.5 score here as well
  results$BA = NA
  results[which(results$site_resid == "residue_wise"),]$BA = 0.5 * (results[which(results$site_resid == "residue_wise"),]$Rc + results[which(results$site_resid == "residue_wise"),]$SP)
  return(results)
}
pfam_dom_motif_eval = function(pfam, prereqs, all_pfam_dom_list){
  if (is.null(all_pfam_dom_list)){
    # TODO change input path
    all_pfam_dom_enrich = download_data(path = "Evaluation_data/All_pfam_dom_enrich.RDS", is.rds = T)
  }
  else{
    all_pfam_dom_enrich = all_pfam_dom_list
  }
  #file_name = pfam_dom_enrich_info[which(pfam_dom_enrich_info$pfam_filter == pfam),]$file_name
  cleaned = all_pfam_dom_enrich[[pfam]]

  ELM_all_instances = prereqs[["ELM_all"]]
  qslim_input = prereqs[["qslim_input"]]
  all_input_prots = unique(unlist(qslim_input))

  ELM_possible = ELM_all_instances[which(ELM_all_instances$Logic == "true positive"),]
  ELM_possible = ELM_possible[which(ELM_possible$uniprot %in% all_input_prots),]
  possible_motifs = unique(ELM_possible$Id)

  all_metrics_per_motif = list()
  weights_per_motif = list()
  for (i in 1:length(possible_motifs)){
    result = calc_metrics_per_motif(cleaned, possible_motifs[i], min_annot = 3, required_data = prereqs)
    all_metrics_per_motif[[i]] = result[["metrics"]]
    weights_per_motif[[i]] = result[["weights"]]
  }
  all_metrics_per_motif = dplyr::bind_rows(all_metrics_per_motif)
  weights_per_motif = dplyr::bind_rows(weights_per_motif)

  final_result = calc_F1_score_HH(weights = weights_per_motif, metrics = all_metrics_per_motif)
  final_result$filter = pfam
  return(final_result)
}

#' Evaluation of protein-domain interactions
#'
#' Calculates performance metrics for evaluating protein-domain interactions between the motif-carrying protein and the enriched domain against ELM interactions.
#'
#' @details
#' For evaluating protein-domain interactions we measured the enrichment of true-positive interactions between a given motif-carrying protein and its associated domains as reported in the ELM interaction dataset. As in the motif-level evaluation, this analysis was performed only on the motif-carrying proteins reported in the ELM interactions dataset, where true-positives represents the number of correctly associated domains for a given motif-carrying protein and then summed over all motif-carrying proteins in the predicted dataset.
#' @return
#' A data frame containing all the relevant protein-domain interactions' performance metrics for each domain enrichment filter.
#' @export
#' @author Bishoy Wadie, Evangelia Petsalaki
#' @importFrom utils read.delim read.table read.csv
#' @examples
#' ProtDom_int_eval_metrics = HV_prot_dom_int_eval()
HV_prot_dom_int_eval = function(bench_data, bench_data_type = c("ELM", "PRMdb"), all_pfam_dom = NULL){
  if (bench_data_type == "ELM" & missing(bench_data)){
    # TODO change input path
    ELM_ints = Get_ELM_all(ELM_data_type = "interactions")
  }
  else{
    ELM_ints = bench_data
  }
  required = list("ELM_ints" = ELM_ints)

  All_prot_domain_F1 = list()
  for (i in 1:nrow(pfam_dom_enrich_info)){
    All_prot_domain_F1[[i]] = pfam_prot_domain_F1(pfam_dom_enrich_info$pfam_filter[i], required_data = required, all_pfam_dom_list = all_pfam_dom)
  }
  All_prot_domain_F1 = dplyr::bind_rows(All_prot_domain_F1)
  All_prot_domain_F1$filter = pfam_dom_enrich_info$pfam_filter

  for (i in 1:nrow(All_prot_domain_F1)){
    All_prot_domain_F1$MCC[i] = mltools::mcc(TP = All_prot_domain_F1$TP[i], FP = All_prot_domain_F1$FP[i],
                                    TN = All_prot_domain_F1$TN[i], FN = All_prot_domain_F1$FN[i])
  }
  return(All_prot_domain_F1)
}

pfam_prot_domain_F1 = function(pfam, required_data, all_pfam_dom_list){
  if (is.null(all_pfam_dom_list)){
    # TODO change input path
    all_pfam_dom_enrich = download_data(path = "Evaluation_data/All_pfam_dom_enrich.RDS", is.rds = T)
  }
  else{
    all_pfam_dom_enrich = all_pfam_dom_list
  }
  #file_name = pfam_dom_enrich_info[which(pfam_dom_enrich_info$pfam_filter == pfam),]$file_name
  cleaned = all_pfam_dom_enrich[[pfam]]

  ELM_interactions = required_data[["ELM_ints"]]
  pfam_nr = all_pfam_dom_enrich[["PF_nr"]]

  filtered_pfam = dplyr::filter(cleaned, uniprot %in% ELM_interactions$interactorElm)
  filtered_pfam$PFAM_id = sub("\\--.*", "",filtered_pfam$PFAM )

  pfam_nr_filt = dplyr::filter(pfam_nr, uniprot %in% filtered_pfam$uniprot)
  pfam_nr_filt$PFAM_id = sub("\\--.*", "",pfam_nr_filt$PFAM)

  ELM_filtered = dplyr::filter(ELM_interactions, interactorElm %in% filtered_pfam$uniprot)

  metrics = rep(0,4)
  names(metrics) = c("TP", "FN", "FP", "TN")

  prots = unique(filtered_pfam$uniprot)
  for (i in 1:length(prots)){
    ELM_filtered_prot = dplyr::filter(ELM_filtered, interactorElm == prots[i])
    filtered_pfam_prot = dplyr::filter(filtered_pfam, uniprot == prots[i])
    pfam_nr_filt_prot = dplyr::filter(pfam_nr_filt, uniprot == prots[i])
    metrics["TP"] = metrics["TP"] + length(dplyr::intersect(ELM_filtered_prot$Domain, filtered_pfam_prot$PFAM_id))
    metrics["FP"] = metrics["FP"] + length(unique(filtered_pfam_prot[which(filtered_pfam_prot$PFAM_id %nin% ELM_filtered_prot$Domain),]$PFAM_id))
    metrics["FN"] = metrics["FN"] + length(unique(ELM_filtered_prot[which(ELM_filtered_prot$Domain %nin% filtered_pfam_prot$PFAM_id),]$Domain))
    metrics["TN"] = metrics["TN"] + length(unique(pfam_nr_filt_prot[which(pfam_nr_filt_prot$PFAM_id %nin% c(filtered_pfam_prot$PFAM_id, ELM_filtered_prot$Domain)),]$PFAM_id))
  }

  result = fisher_and_F1(metrics, pval_method = "greater")
  return(result)
}
