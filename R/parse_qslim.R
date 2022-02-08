#' Parse qslim output per host-viral interaction data
#'
#' Parses the raw combined qslim output per host-viral interaction dataset. Then it adds the corresponding ELM instances along with the CompariMotif similarity between the predicted regular expression motif patterns and that of ELM classes. It also includes the distance between the predicted motif and that of an ELM instance, only if the motif-carrying protein and ELM class match.
#'
#' @return
#' A data frame containing all the predicted hits annotated with corresponding ELM instances and CompariMotif similarity scores.
#'
#' @format
#'  - *Dataset*: Uniprot IDs of Host-viral interaction denoted as Viral-prot_Human-prot
#'  - *Pattern*: Motif's regular expression pattern
#'  - *uniprot*: Uniprot ID of motif-carrying protein
#'  - *Org*: Organism info about the motif-carrying protein
#'  - *Start_Pos*: First position of the predicted motif
#'  - *End_Pos*: Last position of the predicted motif
#'  - *Desc*: Uniprot description of the motif-carrying protein
#'  - *Id*: ELM Identifier
#'  - *Sim_rank*: Ranking of CompariMotif relationships. A value of 9 corresponds to Exact match, while a value of 1 corresponds to Complex Overlap (Check CompariMotif [Website](http://bioware.ucd.ie/~compass/biowareweb/Server_pages/help/slimfinder/comparimotif_help.html#relationships) for more details)
#'  - *Score*: CompariMotif heuristic similarity score, defined as the product of Matched positions and Normalized Information content (Check CompariMotif [paper](https://academic.oup.com/bioinformatics/article/24/10/1307/177233) for more details)
#'  - *Accession*: ELM Instance Accession
#'  - *start*: First position of the ELM motif
#'  - *end*: Last position of the ELM motif
#'  - *Logic*: ELM instance logic as reported by ELM database.
#'  - *motif_distance*: Distance between predicted and ELM motif on the same protein
#'  - *motif_overlap*: Overlap proportion between predicted and ELM motif on the same protein
#'
#' @author Bishoy Wadie, Evangelia Petsalaki
#'
#' @export
#' @importFrom utils read.csv
#' @examples
#' pred_hits = parse_qslim()
parse_qslim = function(main_output = NULL, occ_out = NULL, compari_output = NULL, ELM_input = NULL, preprocess_compari = T, select_top_rank = F, IC_cutoff = 3, remove_motif_variants = T){
  ELM_all_instances = Get_ELM_all(ELM_input = ELM_input, ELM_data_type = "instances")
  ELM_all_instances = ELM_all_instances %>% dplyr::select(Accession, Id, uniprot, start, end, Logic)
  ELM_all_instances = ELM_all_instances[!duplicated(ELM_all_instances),]

  if (remove_motif_variants){
    all_hitsint = remove_redundancy_clouds_variants(main_output = main_output, occ_out = occ_out, select_top_rank = select_top_rank, IC_cutoff = IC_cutoff)
  }
  else{
    all_hitsint = preprocess_slimfinder_out(main_output = main_output, occ_out = occ_out)
    all_hitsint = all_hitsint[["all_hitsint"]]
  }

  if (!preprocess_compari & !is.null(compari_output)){
    comparimotif_data = compari_output
  }
  else{
    comparimotif_data = preprocess_comparimotif(compari_tdt_output = compari_output)
  }
  colnames(comparimotif_data)[4] = "Pattern"


  comparimotif_data = comparimotif_data %>% dplyr::select(Dataset,uniprot, Pattern, Id, Sim_rank, MatchPos,NormIC,Score,new_score)
  comparimotif_data = comparimotif_data[!duplicated(comparimotif_data),]

  new_hits_ints = dplyr::left_join(all_hitsint, comparimotif_data, by = c("Dataset", "Pattern", "uniprot"))

  new_hits_ints = dplyr::left_join(new_hits_ints, ELM_all_instances, by = c("Id","uniprot"))
  new_hits_ints = qslim_ELM_motif_distance(new_hits_ints)
  return(new_hits_ints)
}
preprocess_slimfinder_out = function(main_output = NULL, occ_out = NULL){
  if (is.null(occ_out)){
    # FIXME change path of input data
    input_occ = download_data(path = "Evaluation_data/allhitsperint.occ.csv")
    all_hitsint = read.csv(input_data)
  }
  else{
    all_hitsint = occ_out
  }

  if (is.null(main_output)){
    # TODO add main qslimfinder.csv output to Biostudies and link it here
    input_main = download_data(path = "Evaluation_data/allhitsperint.occ.csv")
    metadata = read.csv(input_main)
  }
  else{
    metadata = main_output
  }

  metadata = metadata %>% dplyr::select(-RunID, -Masking, -Build, -Chance, -RunTime)
  all_hitsint = all_hitsint %>% dplyr::select(-RunID, -Prot_Len, -Variant, -MisMatch, -PepSeq, -PepDesign)

  if (diff(as.numeric(gregexpr("_", all_hitsint$Seq[1])[[1]]))[1] == 1){
    all_hitsint = all_hitsint %>% tidyr::separate(Seq, into = c("uniprot", "Prot_name"), sep = "__") %>%
      tidyr::separate(Prot_name, into = c("remove", "Org"), sep = "_") %>% dplyr::select(-remove)
  }
  else{
    all_hitsint = all_hitsint %>% tidyr::separate(Seq, into = c("Prot_name", "uniprot"), sep = "__")
    if (any(stringr::str_detect(all_hitsint$uniprot,"HUMAN"))){
      colnames(all_hitsint)[which(colnames(all_hitsint) == "uniprot")] = "Org"
      all_hitsint = all_hitsint %>% tidyr::separate(Prot_name, into = c("uniprot", "remove"), sep = "_") %>% dplyr::select(-remove)
    }
    else{
      all_hitsint = all_hitsint %>% tidyr::separate(Prot_name, into = c("remove", "Org"), sep = "_") %>% dplyr::select(-remove)
    }
  }
  all_hitsint = all_hitsint[!duplicated(all_hitsint),]
  #all_hitsint = all_hitsint %>% dplyr::select(-Sig)
  # all_hitsint[which(stringr::str_detect(all_hitsint$uniprot, "-") == TRUE),]$uniprot =
  #   substr(all_hitsint[which(stringr::str_detect(all_hitsint$uniprot, "-") == TRUE),]$uniprot, 1, nchar(all_hitsint[which(stringr::str_detect(all_hitsint$uniprot, "-") == TRUE),]$uniprot) -2 )

  return(list("all_hitsint" = all_hitsint, "metadata" = metadata))
}
qslim_ELM_motif_distance = function(df){
  df$motif_distance = NA
  df$motif_overlap = NA
  for (i in 1:nrow(df)){
    if (is.na(df$start[i])){
      next()
      df$motif_distance[i] = NA
      df$motif_overlap[i] = NA
    }
    if (abs(df$Start_Pos[i]) > abs(df$end[i])){
      df$motif_distance[i] = (abs(df$Start_Pos[i]) - abs(df$end[i]))
      df$motif_overlap[i] = 0
    }
    else if (abs(df$Start_Pos[i]) > abs(df$start[i])){
      if (abs(df$End_Pos[i]) <= abs(df$end[i])){
        df$motif_distance[i] = 0
        df$motif_overlap[i] = (abs(df$End_Pos[i]) - abs(df$Start_Pos[i])) /
          min((abs(df$end[i]) - abs(df$start[i])),(abs(df$End_Pos[i]) - abs(df$Start_Pos[i])))
      }
      else{
        df$motif_distance[i] = 0
        df$motif_overlap[i] = (abs(df$end[i]) - abs(df$Start_Pos[i])) /
          min((abs(df$end[i]) - abs(df$start[i])),(abs(df$End_Pos[i]) - abs(df$Start_Pos[i])))
      }
    }
    else if (abs(df$start[i]) > abs(df$End_Pos[i])){
      df$motif_distance[i] = (abs(df$start[i]) - abs(df$End_Pos[i]))
      df$motif_overlap[i] = 0
    }
    else{
      if (abs(df$end[i]) <= abs(df$End_Pos[i])){
        df$motif_distance[i] = 0
        df$motif_overlap[i] = (abs(df$end[i]) - abs(df$start[i])) /
          min((abs(df$end[i]) - abs(df$start[i])),(abs(df$End_Pos[i]) - abs(df$Start_Pos[i])))
      }
      else{
        df$motif_distance[i] = 0
        df$motif_overlap[i] = (abs(df$End_Pos[i]) - abs(df$start[i])) /
          min((abs(df$end[i]) - abs(df$start[i])),(abs(df$End_Pos[i]) - abs(df$Start_Pos[i])))
      }
    }
  }
  return(df)
}
motif_sim_to_score = function(Sim){
  results = c()
  for (i in 1:length(Sim)){
    rank = switch(Sim[i], "Exact Match" = 9, "Variant Match" = 8, "Degenerate Match" = 8,
                  "Complex Match" = 7, "Exact Parent" = 6, "Exact Subsequence" = 6,
                  "Degenerate Parent" = 5, "Degenerate Subsequence" = 5, "Variant Parent" = 5,
                  "Variant Subsequence" = 5, "Complex Parent" = 4, "Complex Subsequence" = 4,
                  "Exact Overlap" = 3, "Degenerate Overlap" = 2, "Variant Overlap" = 2,
                  "Complex Overlap" = 1)
    results[i] = rank
  }
  return(results)
}
preprocess_comparimotif = function(compari_tdt_output = NULL){
  if (is.null(compari_tdt_output)){
    # FIXME change input path
    comparimotif_Ids = download_data(path = "Evaluation_data/compari_motif_ELM_ids.RDS", is.rds = T)
  }
  else{
    comparimotif_all = compari_tdt_output %>% tidyr::separate(Name1, into = c("Name1", "motif_1_rm"), sep = "[|]")
    comparimotif_all = comparimotif_all %>% tidyr::separate(Name2, into = c("Name2", "motif_2_rm"), sep = "[|]")
    comparimotif_Ids = comparimotif_all %>% dplyr::select(-motif_1_rm, -motif_2_rm)
  }
  #According to qslim_finder paper will apply MatchIC >= 1.5, NormIC >=0.5, matchpos 2+
  comparimotif_Ids = comparimotif_Ids %>% dplyr::filter(MatchIC >= 1.5, NormIC >= 0.5, MatchPos >=2)
  comparimotif_Ids$Sim_rank = motif_sim_to_score(comparimotif_Ids$Sim1)

  comparimotif_Ids = comparimotif_Ids[,c(1:8,17,9:16)]

  all_motifs_regex = unique(c(comparimotif_Ids$Motif1, comparimotif_Ids$Motif2, comparimotif_Ids$Match))

  all_regex_lengths = list()
  for (i in 1:length(all_motifs_regex)){
    motif_memb = c(ifelse(all_motifs_regex[i] %in% comparimotif_Ids$Motif1, "Motif1",NA),
                   ifelse(all_motifs_regex[i] %in% comparimotif_Ids$Motif2, "Motif2",NA),
                   ifelse(all_motifs_regex[i] %in% comparimotif_Ids$Match, "Match",NA))
    motif_memb = motif_memb[!is.na(motif_memb)]
    all_regex_lengths[[i]] = data.frame(motif_regex = all_motifs_regex[i],
                                        motif_len = length_motif_regex(all_motifs_regex[i])[["non_wild"]],
                                        motif_col_memb = motif_memb)
  }
  all_regex_lengths = dplyr::bind_rows(all_regex_lengths)

  for (j in unique(all_regex_lengths$motif_col_memb)){
    x = all_regex_lengths[which(all_regex_lengths$motif_col_memb == j),]
    x = x[,-3]
    colnames(x) = c(j,paste0(j,"_len"))
    comparimotif_Ids = dplyr::left_join(comparimotif_Ids, x, by = j)
  }

  comparimotif_Ids = comparimotif_Ids %>% dplyr::mutate(max_motifs_nw_len=pmax(Motif1_len, Motif2_len))

  comparimotif_Ids$new_score = ifelse(comparimotif_Ids$Sim_rank > 6, comparimotif_Ids$MatchPos / comparimotif_Ids$Match_len,
                                      comparimotif_Ids$MatchPos / comparimotif_Ids$max_motifs_nw_len)
  comparimotif_Ids$new_score = comparimotif_Ids$new_score * comparimotif_Ids$NormIC

  qslim_elm = comparimotif_Ids %>% dplyr::group_by(Name1,Name2, Motif1) %>% dplyr::summarise(n = dplyr::n())
  qslim_elm = as.data.frame(qslim_elm)


  uniques = qslim_elm[which(qslim_elm$n == 1),]
  dupls = qslim_elm[which(qslim_elm$n > 1),]
  dupls = as.data.frame(dupls)
  result = list()
  for (i in 1:nrow(dupls)){
    filtered = dplyr::filter(comparimotif_Ids, Name1 == dupls$Name1[i], Name2 == dupls$Name2[i], Motif1 == dupls$Motif1[i])
    filtered = filtered[which.max(filtered$new_score),]
    result[[i]] = filtered
  }
  result = dplyr::bind_rows(result)
  unique_df = suppressMessages(plyr::match_df(comparimotif_Ids, uniques))
  Final_data = rbind.data.frame(result, unique_df)

  new_compari_results = Final_data %>% tidyr::separate(Name1, into = c("Dataset", "uniprot"), sep = ":")
  colnames(new_compari_results)[5] = "Id"
  new_compari_results = new_compari_results[,-c(1:2)]
  return(new_compari_results)
}
Get_ELM_all = function(ELM_input = NULL, ELM_data_type = c("instances","classes","interactions")){
  # FIXME Fix this function to get ELM from Biostudies or get the latest based on user input
  if (is.null(ELM_input)){
    input_path = switch(ELM_data_type, "instances" = "data/input/ELM/elm_all_instances.csv",
                                         "classes" = "data/input/ELM/elm_classes.csv",
                                         "interactions" = "data/input/ELM/elm_interactions.csv")
    elm_out = read.csv(input_path, header = T)
  }
  else{
    elm_out = ELM_input
  }

  if (ELM_data_type == "instances"){
    colnames(elm_out) = c("Accession","Class_type","Id","Entry_name","uniprot","other_uniprots","start","end","Pubmeds","Method","Logic","PDBs","Organism")
    elm_out = elm_out[which(elm_out$Organism == "Homo sapiens"),]
  }
  return(elm_out)
}
clean_up_patterns = function(new_hits_ints,dataset, motif_prot, remove_CLV_trg = T, sim_score_thres = 2.5){
  .Deprecated("remove_redundancy_clouds_variants")
  A = dplyr::filter(new_hits_ints, Dataset == dataset, uniprot == motif_prot)
  if (nrow(A) == 1){
    A$ELM_class = substr(A$Id, 1, 3)
    if (remove_CLV_trg){
      A = A[which(A$ELM_class %nin% c("CLV", "TRG")),]
    }
    A = A %>% dplyr::select(-ELM_class)
    return(A)
  }
  A = A %>% dplyr::select(Pattern, Start_Pos, End_Pos, Id, Score)
  A$ELM_class = substr(A$Id, 1, 3)
  if (remove_CLV_trg){
    A = A[which(A$ELM_class %nin% c("CLV", "TRG")),]
  }
  A = A[!duplicated(A),]
  A$resids = paste0(A$Start_Pos, ":", A$End_Pos)

  if (is.na(A$Id)){
    C = A
  }
  else{
    C = A[which(A$Score < sim_score_thres),]
  }
  A = A[which(A$Score >= sim_score_thres),]

  if (nrow(A) != 0){
    uni_Ids = unique(A$Id)

    # x = A %>% select(Pattern, Id)
    # x = x %>% group_by(Pattern) %>% dplyr::summarise(n = n())
    # most_common_pattern = x$Pattern[which.max(x$n)]

    final_list = list()
    counter = 0
    for (i in 1:length(uni_Ids)){
      filtered = A[which(A$Id == uni_Ids[i]),]
      clusts = cluster_resids_pairs(unique(filtered$resids))
      inner = list()
      for (j in clusts){
        counter = counter + 1
        overlapp_resids = unlist(strsplit(j,","))
        second = filtered[which(filtered$resid %in% overlapp_resids),]
        second = second[which(second$Score == max(second$Score)),]
        # if (nrow(second) > 1){
        #   if (most_common_pattern %in% second$Id){
        #     second = second[which(second$Id == most_common_pattern),]
        #   }
        #   else{
        #     second = second[which.max(second$Score),]
        #   }
        # }
        inner[[j]] = second
      }
      inner = dplyr::bind_rows(inner)
      inner = inner[!duplicated(inner),]
      inner$diff = inner$End_Pos - inner$Start_Pos
      inner = inner[which.max(inner$diff),]
      inner = inner %>% dplyr::select(-diff)
      final_list[[i]] = inner
    }
    final_list = dplyr::bind_rows(final_list)
    orig = dplyr::filter(new_hits_ints, Dataset == dataset, uniprot == motif_prot)
    orig_new = plyr::match_df(orig, final_list, on = c("Pattern", "Start_Pos", "End_Pos", "Id", "Score"))

    if (nrow(C) != 0){
      C = dplyr::filter(new_hits_ints, Dataset == dataset, uniprot == motif_prot)
      C = C[which(C$Score < sim_score_thres),]
      C$ELM_class = substr(C$Id, 1, 3)
      if (remove_CLV_trg){
        C = C[which(C$ELM_class %nin% c("CLV", "TRG")),]
      }
      C = C %>% dplyr::select(-ELM_class)
      final = rbind.data.frame(orig_new, C)
      final = final[!duplicated(final),]
      return(final)
    }
    else{
      return(orig_new)
    }
  }
  else{
    C = dplyr::filter(new_hits_ints, Dataset == dataset, uniprot == motif_prot)
    if (!is.na(C$Score)){
      C = C[which(C$Score < sim_score_thres),]
    }
    C$ELM_class = substr(C$Id, 1, 3)
    if (remove_CLV_trg){
      C = C[which(C$ELM_class %nin% c("CLV", "TRG")),]
    }
    C = C %>% dplyr::select(-ELM_class)
    return(C)
  }
}
remove_redundancy_clouds_variants = function(main_output = NULL, occ_out = NULL, select_top_rank = T, IC_cutoff = 3){
  # QUESTION If we choose the top ranked patterns, we would favor move occurrences over information content
  # and thus the comparimotif score would be lower when compared with
  all_hitsint = preprocess_slimfinder_out(main_output = main_output, occ_out = occ_out)

  if (is.null(all_hitsint[[2]])){
    message("No motif ranking found, please provide the main output of SLiMFinder")
    return(NULL)
  }
  all_data = suppressMessages(dplyr::left_join(all_hitsint[[1]], all_hitsint[[2]]))

  dtst_uniprot = all_data %>% dplyr::group_by(Dataset, uniprot) %>% dplyr::summarise(n=dplyr::n())

  unique_motifs = all_data[which(all_data$MotNum == 1),]
  variant_motifs = all_data[which(all_data$MotNum > 1),]

  x = plyr::match_df(dtst_uniprot, variant_motifs)
  x = as.data.frame(x)

  final_res = list()
  for (i in 1:nrow(x)){
    filtered = variant_motifs %>% dplyr::filter(Dataset == x$Dataset[i], uniprot == x$uniprot[i])
    if (length(unique(filtered$Rank)) > 1){
      if (select_top_rank | max(filtered$IC) <= IC_cutoff){
        filtered = filtered[which.min(filtered$Rank),]
      }
      else{
        high_IC_patterns = unique(filtered$Pattern[which(filtered$IC > IC_cutoff)])
        filtered = filtered[which(filtered$Pattern %in% high_IC_patterns),]
      }
    }
    final_res[[i]] = filtered
  }
  final_res = dplyr::bind_rows(final_res)
  all_res = rbind.data.frame(unique_motifs, final_res)

  output = plyr::match_df(all_hitsint[[1]], all_res)
  return(output)
}

combine_all_pfam_dom_enrich = function(PfamDomenrich_info = pfam_dom_enrich_info, pfam_dir_path = "Evaluation_data/pfam_dom_enrich"){
  All_cleaned = list()
  for (i in 1:nrow(PfamDomenrich_info)){
    data = read.table(download_data(path = file.path(pfam_dir_path, PfamDomenrich_info$file_name[i])), sep = "\t", na.strings = c("", NA))
    colnames(data) = c("uniprot", "domains", "total_domains_nr")
    cleaned = clean_pfam_dom_data(data)
    All_cleaned[[PfamDomenrich_info$pfam_filter[i]]] = cleaned
  }
  return(All_cleaned)
}
Motifcarr_to_pfam_dom_enrich = function(mot_car_prots, all_pfam_dom_list){
  prots = unique(mot_car_prots)
  if (is.null(all_pfam_dom_list)){
    # TODO change input path
    all_pfam_dom_enrich = download_data(path = "Evaluation_data/All_pfam_dom_enrich.RDS", is.rds = T)
  }
  else{
    all_pfam_dom_enrich = all_pfam_dom_list
  }
  filters = list()
  for (i in 1:length(all_pfam_dom_enrich)){
    cleaned = all_pfam_dom_enrich[[i]]
    if (length(intersect(cleaned$uniprot, prots)) != 0){
      result = data.frame(uniprot = intersect(cleaned$uniprot, prots), pfam_filter = names(all_pfam_dom_enrich)[i])
      filters[[i]] = result
    }
    else{
      filters[[i]] = NULL
    }
  }
  filters = dplyr::bind_rows(filters)

  prots = unique(filters$uniprot)
  pfam_filters_uniprot = data.frame(uniprot = prots, pfam_filter = NA, no_of_filters = NA)
  for (i in 1:nrow(pfam_filters_uniprot)){
    filtered = dplyr::filter(filters, uniprot == pfam_filters_uniprot$uniprot[i])
    result = paste(unique(filtered$pfam_filter), collapse = ";")
    pfam_filters_uniprot$pfam_filter[i] = result
    pfam_filters_uniprot$no_of_filters[i] = length(unique(filtered$pfam_filter))
  }
  return(pfam_filters_uniprot)
}

prepare_pfam_doms_enrich = function(EM_list){
  x = mapply(cbind, EM_list$Enriched_domains, "uniprot"=names(EM_list$Enriched_domains), SIMPLIFY=F)
  x = dplyr::bind_rows(x)

  new_all_pfam_doms = list()
  new_all_pfam_doms[["PF_nr"]] = x
  for (i in 2:nrow(pfam_dom_enrich_info)){
    result = dplyr::filter(x, TP >= pfam_dom_enrich_info$mindom[i],
                           OR > pfam_dom_enrich_info$OR[i],
                           adj.p.val < pfam_dom_enrich_info$adj[i])
    if ("emp_pval" %in% colnames(result)){
      result = result[which(result$emp_pval < 0.05),]
    }
    new_all_pfam_doms[[pfam_dom_enrich_info$pfam_filter[i]]] = result
  }
  return(new_all_pfam_doms)
}
cluster_resids_pairs = function(resids_vec, overlap_thres = 0.8){
  resid_pairs = gimme:::expand.grid.unique(resids_vec, resids_vec)
  resid_pairs = as.data.frame(resid_pairs)
  colnames(resid_pairs) = c("res1", "res2")
  resid_pairs$int = 0

  for (i in 1:nrow(resid_pairs)){
    x = c(as.numeric(strsplit(resid_pairs$res1[i], ":")[[1]][1]):as.numeric(strsplit(resid_pairs$res1[i], ":")[[1]][2]))
    y = c(as.numeric(strsplit(resid_pairs$res2[i], ":")[[1]][1]):as.numeric(strsplit(resid_pairs$res2[i], ":")[[1]][2]))
    resid_pairs$int[i] = length(intersect(x,y)) / min(length(x),length(y))
  }

  resids_pairs_reverse = resid_pairs
  colnames(resids_pairs_reverse) = c("res2","res1", "int")
  final_resid = rbind.data.frame(resid_pairs, resids_pairs_reverse)

  final_data = data.frame(seed = resids_vec, clust = NA)

  for (i in 1:nrow(final_data)){
    filtered = final_resid[which(final_resid$res1 == final_data$seed[i]),]
    filtered = filtered[which(filtered$int >= overlap_thres),]
    final_data$clust[i] = paste(unique(filtered$res2), collapse = ",")
  }
  final_data = final_data[!duplicated(final_data$clust),]
  return(final_data$clust)
}

generate_potential_candidates = function(new_hits_ints, pep_clin_EM, clinvar_path,
                                         PTM_info, pfam_doms, compari_cutoff = 0.66){
  pepsite_clinvar = plyr::match_df(clinvar_path, pep_clin_EM, on = c("uniprot", "Start_Pos", "End_Pos"))
  pepsite_clinvar = pepsite_clinvar[is.na(pepsite_clinvar$Accession),]

  ABC = dplyr::left_join(pep_clin_EM, pepsite_clinvar)

  PTM_info = PTM_info %>% dplyr::select(Dataset, uniprot, Start_Pos, End_Pos, PTM_res, PTM_type, PTM_comments)
  PTM_info = PTM_info[!duplicated(PTM_info),]
  PTM_info = PTM_info[!is.na(PTM_info$PTM_res),]
  ABC = dplyr::left_join(ABC, PTM_info)

  pfam_id_annot = pfam_doms %>% tidyr::separate(PFAM, into = c("pfam_id", "pfam_name"), sep = "--") %>% dplyr::select(pfam_id, pfam_name)
  pfam_id_annot = pfam_id_annot[!duplicated(pfam_id_annot),]


  ABC = dplyr::left_join(ABC, pfam_id_annot)

  # all_hitsint = dplyr::select(new_hits_ints, Dataset, Pattern, Sig, uniprot, Start_Pos, End_Pos)
  # ABC = dplyr::left_join(ABC, all_hitsint)

  pep_clin_inst = ABC %>% dplyr::select(uniprot, Start_Pos, End_Pos)
  pep_clin_inst = pep_clin_inst[!duplicated(pep_clin_inst),]

  compiled = list()
  for (i in 1:nrow(pep_clin_inst)){
    filtered = dplyr::filter(ABC, uniprot == pep_clin_inst$uniprot[i], Start_Pos == pep_clin_inst$Start_Pos[i],
                             End_Pos == pep_clin_inst$End_Pos[i])
    matched_orig = plyr::match_df(new_hits_ints, filtered, on = c("Dataset", "uniprot","Start_Pos", "End_Pos", "Pattern")) %>% dplyr::select(Dataset, uniprot,Start_Pos, End_Pos, Pattern,Id, new_score)
    filtered = dplyr::left_join(filtered, matched_orig)
    if (all(!is.na(filtered$Id))){
      filtered$Instance_type = ifelse(filtered$new_score > compari_cutoff, "New_instance", "Novel")
    }
    else{
      filtered$Instance_type = "Novel"
    }

    filtered = dplyr::select(filtered, Dataset, uniprot, GeneSymbol, ENTREZID, Desc,
                      Pattern, peptide_seq, Start_Pos, End_Pos, Sig, Type, Name,
                      ClinicalSignificance, RCVaccession, PhenotypeList,
                      OriginSimple, Chromosome, AA_change, fromAA,toAA, PTM_res, PTM_type, PTM_comments,
                      pfam_id, pfam_name, domain_protein, PDB_ID, Chain,
                      experimentalTechnique, resolution, ligandId, pdb_start, pdb_end,
                      uniprot_start, uniprot_end, pepsite_score, p.val, peptide_length,
                      no_of_binding, binding_resids, Id,Rank, new_score, Instance_type)

    if (all(filtered$Instance_type == "Novel")){
      filtered$Id = NA
      filtered$new_score = NA
      filtered = filtered[!duplicated(filtered),]
      compiled[[i]] = filtered
    }
    else{
      filtered = dplyr::filter(filtered, Instance_type == "New_instance")
      filtered = filtered[!duplicated(filtered),]
      compiled[[i]] = filtered
    }
  }
  compiled = dplyr::bind_rows(compiled)
  #compiled$Score = compiled$Score / 5
  compiled = compiled %>%
    dplyr::rename(motif_protein = uniprot, motif_prot_Desc = Desc,
                  qslim_pval = Sig, Variant_Type = Type, Variant_Name = Name,
                  Clinvar_accession = RCVaccession, Disease_name = PhenotypeList,
                  pepsite_p.val = p.val, ELM_Class = Id, CompariMotif_norm_Score = new_score)

  # if (is.null(all_hits_ints)){
  #   compiled = compiled %>% tidyr::separate(Dataset, into = c("V1", "H1"), sep = "_")
  #   colnames(compiled)[c(3,6,11,12,13,15,16,17,39,43,44)] = c("motif_protein", "motif_prot_Desc", "qslim_pval", "Variant_Type", "Variant_Name",
  #                                                             "dbSNP_id", "Clinvar_accession", "Disease_name", "pepsite_p.val",
  #                                                             "ELM_Class", "CompariMotif_Sim_Score")
  # }
  # else{
  #   colnames(compiled)[c(2,5,10,11,12,14,15,16,38,42,43)] = c("motif_protein", "motif_prot_Desc", "qslim_pval", "Variant_Type", "Variant_Name",
  #                                                             "dbSNP_id", "Clinvar_accession", "Disease_name", "pepsite_p.val",
  #                                                             "ELM_Class", "CompariMotif_Sim_Score")
  # }

  return(compiled)
}

get_classified_insts = function(new_hits_ints, clinvar_path, PTM_info, pepsite_hits, Drug_hits,HMM_iELM, EM_list, compari_cutoff = 0.5, mindom = 5, adjpval = 0.05, host_viral = T){
  qslim_instances = new_hits_ints %>% dplyr::select(uniprot, Start_Pos, End_Pos)
  qslim_instances = qslim_instances[!duplicated(qslim_instances),]

  classified_qslim_insts = list()
  pb = txtProgressBar(min = 0, max = nrow(qslim_instances), initial = 0, style = 3)
  for (i in 1:nrow(qslim_instances)){
    classified_qslim_insts[[i]] = Classify_each_instance(new_hits_ints = new_hits_ints, clinvar_path = clinvar_path,
                                                         PTM_info = PTM_info, pepsite_hits = pepsite_hits, Drug_hits = Drug_hits,
                                                         mot_prot = qslim_instances$uniprot[i], mot_start = qslim_instances$Start_Pos[i],
                                                         mot_end = qslim_instances$End_Pos[i])
    setTxtProgressBar(pb = pb, i)
  }
  classified_qslim_insts = dplyr::bind_rows(classified_qslim_insts)
  if (host_viral){
    classified_qslim_insts = classified_qslim_insts %>% tidyr::separate(Dataset, into = c("V1", "H1"), sep = "_")
  }
  passed_HMM_insts = HMM_iELM_classify(new_hits_ints = new_hits_ints, HMM_iELM = HMM_iELM, compariMotif_cutoff = compari_cutoff, host_viral = host_viral)
  classified_qslim_insts = dplyr::left_join(classified_qslim_insts, passed_HMM_insts)
  classified_qslim_insts[is.na(classified_qslim_insts$passed_iELM),]$passed_iELM = 0

  passed_EM = Classify_Enrich_dom_inst(EM_list = EM_list, new_hits_ints = new_hits_ints)
  passed_EM = passed_EM[which(passed_EM$ED_adj_pval == adjpval),]
  passed_EM = passed_EM[which(passed_EM$min_ED == mindom),]
  passed_EM$ED_M5_A0.05 = 1
  passed_EM = passed_EM[,c(1:3,6)]
  classified_qslim_insts = dplyr::left_join(classified_qslim_insts, passed_EM)
  classified_qslim_insts[is.na(classified_qslim_insts$ED_M5_A0.05),]$ED_M5_A0.05 = 0
  return(classified_qslim_insts)
}
Classify_each_instance = function(new_hits_ints, clinvar_path, PTM_info, pepsite_hits, Drug_hits, mot_prot, mot_start, mot_end){
  qslim_instances = new_hits_ints %>% dplyr::select(uniprot, Start_Pos, End_Pos)
  qslim_instances = qslim_instances[!duplicated(qslim_instances),]

  PTM_info = PTM_info[!is.na(PTM_info$PTM_res),]

  V1_H1_H2 = new_hits_ints %>% dplyr::filter(uniprot == mot_prot, Start_Pos == mot_start,
                                             End_Pos == mot_end) %>% dplyr::select(Dataset, uniprot)
  V1_H1_H2 = V1_H1_H2[!duplicated(V1_H1_H2),]

  insts = qslim_instances %>% dplyr::filter(uniprot == mot_prot, Start_Pos == mot_start,
                                            End_Pos == mot_end)
  insts = dplyr::left_join(insts, V1_H1_H2, by = "uniprot")
  insts = insts[,c(4,1:3)]

  insts_known_ELM = new_hits_ints[!is.na(new_hits_ints$Accession),] %>% dplyr::filter(motif_overlap != 0)

  has_clinvar = plyr::match_df(insts, clinvar_path, on = c("uniprot", "Start_Pos", "End_Pos"))
  has_PTM = plyr::match_df(insts, PTM_info, on = c("uniprot", "Start_Pos", "End_Pos"))
  known_ELM = plyr::match_df(insts, insts_known_ELM, on = c("uniprot", "Start_Pos", "End_Pos"))

  has_pepsite = pepsite_hits %>% dplyr::filter(motif_protein == mot_prot, start == mot_start,
                                                          end == mot_end)
  passed_chembl = plyr::match_df(Drug_hits, V1_H1_H2, on = c("Dataset", "uniprot"))

  insts$PepSite = ifelse(nrow(has_pepsite) != 0, 1, 0)
  insts$Clinvar_pathogenic = ifelse(nrow(has_clinvar) != 0, 1, 0)
  insts$has_PTM = ifelse(nrow(has_PTM) != 0, 1, 0)
  insts$Known_ELM = ifelse(nrow(known_ELM) != 0, 1, 0)
  insts$ChEMBL_repurp = ifelse(nrow(passed_chembl) != 0, 1, 0)


  return(insts)
}
HMM_iELM_classify = function(new_hits_ints, HMM_iELM, compariMotif_cutoff = 0.5, host_viral = T){
  if (host_viral){
    pass_new_hits_int = new_hits_ints %>% tidyr::separate(Dataset, into = c("V1", "H1"), sep = "_")
  }
  else{
    pass_new_hits_int = new_hits_ints
    colnames(pass_new_hits_int)[which(names(pass_new_hits_int) == "Dataset")] = "H1"
  }
  pass_new_hits_int = pass_new_hits_int[which(pass_new_hits_int$H1 %in% HMM_iELM$H1),]
  pass_new_hits_int = pass_new_hits_int[which(pass_new_hits_int$new_score >= compariMotif_cutoff),]


  qslim_instances = new_hits_ints %>% dplyr::select(Dataset, uniprot, Start_Pos, End_Pos)
  if (host_viral){
    qslim_instances = qslim_instances %>% tidyr::separate(Dataset, into = c("V1", "H1"), sep = "_")
  }
  else{
    colnames(qslim_instances)[which(names(qslim_instances) == "Dataset")] = "H1"
  }
  qslim_instances = qslim_instances[!duplicated(qslim_instances),]

  classified = plyr::match_df(qslim_instances, pass_new_hits_int)
  classified_insts = classified %>% dplyr::select(H1, uniprot, Start_Pos, End_Pos)
  classified_insts = classified_insts[!duplicated(classified_insts),]
  classified_insts$passed_iELM = 0
  H1s = unique(classified_insts$H1)
  H1_Ids = list()
  for (i in 1:length(H1s)){
    filtered = dplyr::filter(HMM_iELM, H1 == H1s[i])
    H1_Ids[[H1s[i]]] = unique(filtered$ELM_Class)
  }
  for (i in 1:nrow(classified_insts)){
    x = plyr::match_df(pass_new_hits_int, classified_insts[i,], on = c("H1", "uniprot", "Start_Pos", "End_Pos"))
    uni_Ids = unique(x$Id)
    for (k in uni_Ids){
      if (any(stringr::str_detect(k, H1_Ids[[classified_insts$H1[i]]]))){
        classified_insts$passed_iELM[i] = 1
        break()
      }
    }
  }
  return(classified_insts)
}
Classify_Enrich_dom_inst = function(EM_list, new_hits_ints){
  qslim_instances = dplyr::select(new_hits_ints, uniprot, Start_Pos, End_Pos)
  qslim_instances = qslim_instances[!duplicated(qslim_instances),]

  A = EM_list$Enriched_domains
  A = mapply(cbind, A, "uniprot"=names(A), SIMPLIFY=F)
  A = dplyr::bind_rows(A)

  combs = pfam_dom_enrich_info[-1,-1] %>% dplyr::select(-OR)
  combs = combs[!duplicated(combs),]

  final = list()
  for (i in 1:nrow(combs)){
    if ("emp_pval" %in% names(A)){
      filtA = A %>% dplyr::filter(TP >= combs$mindom[i], adj.p.val < combs$adj[i], emp_pval < 0.05)
    }
    else{
      filtA = A[which(A$TP >= combs$mindom[i] & A$adj.p.val < combs$adj[i]),]
    }
    filtA = filtA %>% tidyr::separate(PFAM, into = c("pfam_id", "pfam_name"), sep = "--")

    x = filtA %>% dplyr::select(pfam_id, uniprot)
    x = x[!duplicated(x),]
    filtA = dplyr::left_join(x,qslim_instances, by = "uniprot")
    filtA$min_ED = combs$mindom[i]
    filtA$ED_adj_pval = combs$adj[i]
    filtA = filtA[,-1]
    final[[i]] = filtA
  }
  final = dplyr::bind_rows(final)
  final = final[!duplicated(final),]
  return(final)
}

Classify_insts_by_score = function(new_hits_ints){
  known = new_hits_ints[!is.na(new_hits_ints$Accession),]
  known = known[which(known$motif_overlap != 0),]
  known$inst_type = "known"

  novel = new_hits_ints[which(new_hits_ints$new_score < 0.3),]
  new_insts = new_hits_ints[which(new_hits_ints$new_score > 0.66),]
  new_insts = mismatch_df(new_insts, known, on = c("uniprot", "Start_Pos", "End_Pos"))
  #new_insts = mismatch_df(new_insts, putative, on = c("uniprot", "Start_Pos", "End_Pos"))
  new_insts$inst_type = "new_insts"

  putative = new_hits_ints[which(new_hits_ints$new_score >= 0.3 & new_hits_ints$new_score <= 0.66),]
  putative = mismatch_df(putative, novel, on = c("uniprot", "Start_Pos", "End_Pos"))
  putative = mismatch_df(putative, new_insts, on = c("uniprot", "Start_Pos", "End_Pos"))
  putative = mismatch_df(putative, known, on = c("uniprot", "Start_Pos", "End_Pos"))
  putative$inst_type = "putative"

  novel = rbind.data.frame(novel, new_hits_ints[is.na(new_hits_ints$new_score),])
  novel = mismatch_df(novel, new_insts, on = c("uniprot", "Start_Pos", "End_Pos"))
  novel = mismatch_df(novel, known, on = c("uniprot", "Start_Pos", "End_Pos"))
  novel = mismatch_df(novel, putative, on = c("uniprot", "Start_Pos", "End_Pos"))
  novel$inst_type = "novel"

  final = rbind.data.frame(known, new_insts, putative, novel)
  final = final %>% dplyr::select(uniprot, Start_Pos, End_Pos, inst_type)
  final = final[!duplicated(final),]
  return(final)
}

compare_insts_uniprot_bed = function(new_hits_ints, path_to_bed_files_dir){
  bed_files = list.files(path_to_bed_files_dir, full.names = T)
  all_bed_parsed = list()
  for (i in 1:length(bed_files)){
    bed = read.table(bed_files[i], sep = "\t", quote = "")
    bed = bed[which(bed$V4 %in% new_hits_ints$uniprot),]
    bed = bed[,c(4,14)]
    bed = bed %>% tidyr::separate(V14, into = c("site", "comment"), sep = ";", extra = "merge")
    bed$bed_type = sub(".*UP000005640_9606_","",bed_files[i])
    bed$bed_type = sub("[.].*", "", bed$bed_type)
    all_bed_parsed[[i]] = bed
  }
  all_bed_parsed = dplyr::bind_rows(all_bed_parsed)
  single = all_bed_parsed[stringr::str_which(all_bed_parsed$site, "[,]", negate = T),]
  single = single %>% tidyr::separate(site, into = c("start", "end"), sep = "-")
  single$end[is.na(single$end)] = single$start[is.na(single$end)]

  multiple = all_bed_parsed[stringr::str_which(all_bed_parsed$site, "[,]", negate = F),]
  mult_resolved = list()
  for (i in 1:nrow(multiple)){
    filtered = multiple[i,]
    sites = unlist(strsplit(filtered$site,","))
    new = data.frame(V4 = filtered$V4, site = sites, comment = filtered$comment, bed_type = filtered$bed_type)
    new = new %>% tidyr::separate(site, into = c("start", "end"), sep = "-")
    new$end[is.na(new$end)] = new$start[is.na(new$end)]
    mult_resolved[[i]] = new
  }
  mult_resolved = dplyr::bind_rows(mult_resolved)

  all_bed_resolved = rbind.data.frame(single, mult_resolved)
  colnames(all_bed_resolved)[1] = "uniprot"
  all_bed_resolved$start = substr(all_bed_resolved$start, 2, nchar(all_bed_resolved$start))
  all_bed_resolved$end = substr(all_bed_resolved$end, 2, nchar(all_bed_resolved$end))
  all_bed_resolved$start = as.integer(all_bed_resolved$start)
  all_bed_resolved$end = as.integer(all_bed_resolved$end)

  qslim_instances = new_hits_ints %>% dplyr::select(uniprot, Start_Pos, End_Pos)
  qslim_instances = qslim_instances[!duplicated(qslim_instances),]

  qslim_instances = dplyr::left_join(qslim_instances, all_bed_resolved)

  final = qslim_ELM_motif_distance(qslim_instances)
  return(final)
}

prepare_input_comparimotif = function(pred_hits, ELM_classes){
  pred_hits = pred_hits[which(colnames(pred_hits) %in% c("Dataset", "uniprot", "Pattern")),]
  pred_hits$Name = paste0(pred_hits$Dataset, ":" ,pred_hits$uniprot)
  hits_to_compari = pred_hits %>% dplyr::select(Name, Pattern)

  ELM_to_compari = ELM_classes %>% dplyr::select(ELMIdentifier, Regex)
  colnames(ELM_to_compari) = c("Name", "Pattern")

  return(list("hits_to_CompariMotif" = hits_to_compari, "ELM_to_CompariMotif" = ELM_to_compari))
}
