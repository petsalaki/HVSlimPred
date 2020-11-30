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
parse_qslim = function(){
  input_data = download_data(path = "Evaluation_data/allhitsperint.occ.csv")
  comparimotif_data = preprocess_comparimotif()
  ELM_all_instances = Get_ELM_all()

  all_hitsint = read.csv(input_data)
  all_hitsint = all_hitsint[,-c(2,3,12,14)]
  all_hitsint$uniprot = NA
  all_hitsint$Org = NA
  for (i in 1:nrow(all_hitsint)){
    all_hitsint$uniprot[i] = strsplit(all_hitsint$Seq[i], split = "_")[[1]][1]
    all_hitsint$Org[i] = strsplit(all_hitsint$Seq[i], split = "_")[[1]][4]
  }
  all_hitsint = all_hitsint[,c(1:3,12,13,5:11)]
  all_hitsint = all_hitsint[!duplicated(all_hitsint),]

  colnames(comparimotif_data)[4] = "Pattern"
  new_hits_ints = dplyr::left_join(all_hitsint, comparimotif_data[,c(1,2,4,3,8,14)], by = c("Dataset", "Pattern", "uniprot"))
  new_hits_ints = new_hits_ints %>% dplyr::select(-Sig, -Prot_Len, -Variant, -Match)
  new_hits_ints[which(stringr::str_detect(new_hits_ints$uniprot, "-") == TRUE),]$uniprot =
    substr(new_hits_ints[which(stringr::str_detect(new_hits_ints$uniprot, "-") == TRUE),]$uniprot, 1, nchar(new_hits_ints[which(stringr::str_detect(new_hits_ints$uniprot, "-") == TRUE),]$uniprot) -2 )

  new_hits_ints = dplyr::left_join(new_hits_ints, ELM_all_instances[,c(1,3,5,7,8,11)], by = c("Id","uniprot"))
  new_hits_ints = qslim_ELM_motif_distance(new_hits_ints)
  new_hits_ints = new_hits_ints %>% dplyr::select(-PepDesign)
  return(new_hits_ints)
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
    if (abs(df$start[i]) > abs(df$End_Pos[i])){
      df$motif_distance[i] = (abs(df$start[i]) - abs(df$End_Pos[i]))
      df$motif_overlap[i] = 0
    }
    else if (abs(df$start[i]) > abs(df$Start_Pos[i])){
      df$motif_distance[i] = (abs(df$start[i]) - abs(df$Start_Pos[i]))
      df$motif_overlap[i] = (abs(df$End_Pos[i]) - abs(df$start[i])) / (abs(df$End_Pos[i]) - abs(df$Start_Pos[i]))
    }
    else if (abs(df$Start_Pos[i]) > abs(df$end[i])){
      df$motif_distance[i] = (abs(df$Start_Pos[i]) - abs(df$end[i]))
      df$motif_overlap[i] = 0
    }
    else{
      df$motif_distance[i] = (abs(df$Start_Pos[i]) - abs(df$start[i]))
      df$motif_overlap[i] = (abs(df$end[i]) - abs(df$Start_Pos[i])) / (abs(df$end[i]) - abs(df$start[i]))
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
preprocess_comparimotif = function(){
  comparimotif_Ids = download_data(path = "Evaluation_data/compari_motif_ELM_ids.RDS", is.rds = T)
  #According to qslim_finder paper will apply MatchIC >= 1.5, NormIC >=0.5, matchpos 2+
  comparimotif_Ids = comparimotif_Ids %>% dplyr::filter(MatchIC >= 1.5, NormIC >= 0.5, MatchPos >=2)
  comparimotif_Ids$Sim_rank = motif_sim_to_score(comparimotif_Ids$Sim1)

  comparimotif_Ids = comparimotif_Ids[,c(1:8,17,9:16)]

  qslim_elm = comparimotif_Ids %>% dplyr::group_by(Name1,Name2, Motif1) %>% dplyr::summarise(n = n())
  qslim_elm = as.data.frame(qslim_elm)

  uniques = qslim_elm[which(qslim_elm$n == 1),]
  dupls = qslim_elm[which(qslim_elm$n > 1),]
  dupls = as.data.frame(dupls)
  result = list()
  for (i in 1:nrow(dupls)){
    filtered = dplyr::filter(comparimotif_Ids, Name1 == dupls$Name1[i], Name2 == dupls$Name2[i], Motif1 == dupls$Motif1[i])
    filtered = filtered[which.max(filtered$Score),]
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
Get_ELM_all = function(){
  ELM_input = download_data(path = "Evaluation_data/elm_all_instances.csv")
  ELM_all_instances = read.csv(ELM_input)
  colnames(ELM_all_instances) = c("Accession","Class_type","Id","Entry_name","uniprot","other_uniprots","start","end","Pubmeds","Method","Logic","PDBs","Organism")
  ELM_all_instances[which(stringr::str_detect(ELM_all_instances$uniprot, "-") == TRUE),]$uniprot =
    substr(ELM_all_instances[which(stringr::str_detect(ELM_all_instances$uniprot, "-") == TRUE),]$uniprot, 1, nchar(ELM_all_instances[which(stringr::str_detect(ELM_all_instances$uniprot, "-") == TRUE),]$uniprot) -2 )

  return(ELM_all_instances)
}




