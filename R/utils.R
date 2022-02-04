#' @export
`%nin%` <- Negate("%in%")


mismatch_df = function (x, y, on = NULL){
  if (is.null(on)) {
    on <- dplyr::intersect(names(x), names(y))
    message("Matching on: ", paste(on, collapse = ", "))
  }
  keys <- plyr::join.keys(x, y, on)
  x[keys$x %nin% keys$y, , drop = FALSE]
}

#' @importFrom utils download.file
download_data = function(path, is.rds = F){
  # FIXME change ftp base to point to Biostudies
  ftp_base = "ftp://ftp.ebi.ac.uk/pub/contrib/petsalaki/Wadie_et_al"
  if (!is.rds){
    tmp = tempfile()
    data = download.file(file.path(ftp_base, path), destfile = tmp)
    return(tmp)
  }
  else{
    data = readRDS(url(file.path(ftp_base, path), "rb"))
    return(data)
  }
}

create_output_path = function(outdir_path, file_name){
  if (substr(outdir_path, nchar(outdir_path), nchar(outdir_path)) == "/"){
    output_file = paste0(outdir_path, file_name)
  }
  else{
    output_file = file.path(outdir_path, file_name)
  }
  return(output_file)
}

count_motif_instances = function(df, inst_cols){
  if (missing(inst_cols)){
    df = df %>% dplyr::select(uniprot, Start_Pos, End_Pos)
    df = df[!duplicated(df),]
    return(nrow(df))
  }
  else{
    df = df[,inst_cols]
    df = df[!duplicated(df),]
    return(nrow(df))
  }
}

count_peptide_domain_pepsite = function(new_hits_ints,human_human_int, EM_list, pfam_doms_list, pfam_pdb, reqs_raw, all_motif_domain_ints, host_viral = T){
  # TODO check this function if anything need to be added to account for alphafold results either seperately
  # or collectively.
  all_ints_pf_nr = list()
  x = unique(names(EM_list$pfam_doms_nr))
  for (i in 1:length(x)){
    ints = data.frame(motif_protein = x[i], domain_protein = get_human_interactors(human_human_intact = human_human_int, uniprot_id = x[i], return_fasta = F))
    all_ints_pf_nr[[i]] = ints
  }
  all_ints_pf_nr = dplyr::bind_rows(all_ints_pf_nr)

  H1_pred_motifs = new_hits_ints %>% dplyr::select(Dataset,uniprot, Match, Org, Start_Pos, End_Pos)

  H1_pred_motifs = H1_pred_motifs[!duplicated(H1_pred_motifs),]
  if (host_viral){
    H1_pred_motifs = H1_pred_motifs %>% tidyr::separate(Dataset, into = c("V1", "H1"), sep = "_")
  }
  else{
    colnames(H1_pred_motifs)[1] = "H1"
  }
  H1_pred_motifs = dplyr::left_join(H1_pred_motifs, pfam_doms_list, by = c("H1" = "uniprot"))
  H1_pred_motifs = H1_pred_motifs %>% tidyr::separate(PFAM, into = c("pfam_id", "pfam_name"), sep = "--")
  before_struct_H1 = H1_pred_motifs %>% dplyr::select(Match, pfam_id)
  before_struct_H1 = before_struct_H1[!duplicated(before_struct_H1),]

  x = pfam_pdb
  colnames(x)[which(names(x) == "uniprot")] = "H1"
  #H1_pred_motifs = plyr::match_df(H1_pred_motifs, x)
  H1_pred_motifs_AF = H1_pred_motifs[which(H1_pred_motifs$pfam_id %in% reqs_raw$pfam_id[which(reqs_raw$experimentalTechnique == "Prediction")]),]
  H1_pred_motifs = H1_pred_motifs[which(H1_pred_motifs$pfam_id %in% pfam_pdb$pfam_id),]

  after_struct_H1 = H1_pred_motifs %>% dplyr::select(Match, pfam_id)
  after_struct_H1_AF = H1_pred_motifs_AF %>% dplyr::select(Match, pfam_id)
  after_struct_H1 = plyr::match_df(before_struct_H1, after_struct_H1)
  after_struct_H1_AF = plyr::match_df(before_struct_H1, after_struct_H1_AF)


  pf_nr = dplyr::left_join(all_ints_pf_nr, pfam_doms_list, by = c("domain_protein" = "uniprot"))
  pf_nr = pf_nr %>% tidyr::separate(PFAM, into = c("pfam_id", "pfam_name"), sep = "--")
  pf_nr = pf_nr[!duplicated(pf_nr),]

  pred_motifs = new_hits_ints %>% dplyr::select(uniprot, Match, Org, Start_Pos, End_Pos)
  pred_motifs = pred_motifs[!duplicated(pred_motifs),]
  pred_motifs = pred_motifs[which(pred_motifs$uniprot %in% pf_nr$motif_protein),]

  pf_nr = dplyr::left_join(pf_nr, pred_motifs, by = c("motif_protein" = "uniprot"))

  before_struct_EM = pf_nr %>% dplyr::select(Match, pfam_id)
  before_struct_EM = before_struct_EM[!duplicated(before_struct_EM),]

  colnames(pf_nr)[which(names(pf_nr) == "domain_protein")] = "uniprot"
  #pf_nr = plyr::match_df(pf_nr, pfam_pdb) # remove domains with no available structures
  pf_nr_AF = pf_nr[which(pf_nr$pfam_id %in% reqs_raw$pfam_id[which(reqs_raw$experimentalTechnique == "Prediction")]),]
  pf_nr = pf_nr[which(pf_nr$pfam_id %in% pfam_pdb$pfam_id),]

  after_struct_EM = pf_nr %>% dplyr::select(Match, pfam_id)
  after_struct_EM = plyr::match_df(before_struct_EM, after_struct_EM)
  after_struct_EM_AF = pf_nr_AF %>% dplyr::select(Match, pfam_id)
  after_struct_EM_AF = plyr::match_df(before_struct_EM, after_struct_EM_AF)

  before_struct_all = rbind.data.frame(before_struct_H1, before_struct_EM)
  before_struct_all = before_struct_all[!duplicated(before_struct_all),]

  after_struct_all = rbind.data.frame(after_struct_H1, after_struct_EM)
  after_struct_all = after_struct_all[!duplicated(after_struct_all),]


  after_struct_all_AF = rbind.data.frame(after_struct_H1_AF, after_struct_EM_AF)
  after_struct_all_AF = after_struct_all_AF[!duplicated(after_struct_all_AF),]
  after_struct_all_AF = mismatch_df(after_struct_all_AF, after_struct_all)

  x1 = reqs_raw
  colnames(x1)[which(colnames(x1) == "motif_protein")] = "uniprot"
  x1 = plyr::match_df(x1, new_hits_ints)

  pepsite_reqs_count = dplyr::select(x1, peptide, pfam_id)
  pepsite_reqs_count = pepsite_reqs_count[!duplicated(pepsite_reqs_count),]
  colnames(pepsite_reqs_count)[1] = "Match"
  #pepsite_reqs_count = plyr::match_df(pepsite_reqs_count, after_struct_all)

  # colnames(pepsite_out)[c(4,5,6)] = c("Start_Pos", "End_Pos", "uniprot")
  # pepsite_out = plyr::match_df(pepsite_out, new_hits_ints, on = c("Start_Pos", "End_Pos", "uniprot"))
  #
  # colnames(pepsite_out)[c(4,5,6)] = c("start","end","motif_protein")
  # final_reqs_raw = reqs_raw
  # colnames(final_reqs_raw)[c(1,14,16,17)] = c("PDB_ID", "peptide_seq", "start", "end")
  # trial = plyr::match_df(final_reqs_raw, pepsite_out)
  #
  # motif_domain_ints = dplyr::left_join(trial, pepsite_out)
  #
  # before_domain_int = dplyr::select(motif_domain_ints, peptide_seq, pfam_id)
  # colnames(before_domain_int) = c("Match", "pfam_id")
  # before_domain_int = before_domain_int[!duplicated(before_domain_int),]
  # before_domain_int = plyr::match_df(before_domain_int, pepsite_reqs_count)

  after_domain_int = dplyr::select(all_motif_domain_ints, peptide_seq, pfam_id)
  colnames(after_domain_int) = c("Match", "pfam_id")
  after_domain_int = after_domain_int[!duplicated(after_domain_int),]
  after_domain_int = plyr::match_df(after_domain_int, pepsite_reqs_count)

  reported_numbers = c("Before_pdb_struct" = nrow(before_struct_all), "after_pdb_struct" = nrow(after_struct_all),
                       "after_AF_struct" = nrow(after_struct_all_AF),
                       "sent_to_pepsite" = nrow(pepsite_reqs_count),
                       "After_domain_int" = nrow(after_domain_int))
  return(reported_numbers)
}

classify_instance_type = function(new_hits_ints, compari_cutoff = 2){
  # FIXME fix this function based on new definition of novel and new instances
  qslim_insts = new_hits_ints %>% dplyr::select(uniprot, Start_Pos, End_Pos)
  qslim_insts = qslim_insts[!duplicated(qslim_insts),]

  known_insts = new_hits_ints[which(new_hits_ints$motif_distance == 0 & !is.na(new_hits_ints$Accession)),]
  known_insts = known_insts %>% dplyr::select(uniprot, Start_Pos, End_Pos)
  known_insts = known_insts[!duplicated(known_insts),]


  pb = txtProgressBar(min = 0, max = nrow(qslim_insts), initial = 0, style = 3)
  compiled = list()
  for (i in 1:nrow(qslim_insts)){
    filtered = dplyr::filter(new_hits_ints, uniprot == qslim_insts$uniprot[i], Start_Pos == qslim_insts$Start_Pos[i],
                             End_Pos == qslim_insts$End_Pos[i])
    filtered = filtered %>% dplyr::select(Dataset, uniprot,Start_Pos, End_Pos, Pattern,Id, Score)
    if (all(!is.na(filtered$Id))){
      filtered$Instance_type = ifelse(filtered$Score >= compari_cutoff, "New_instance", "Novel")
    }
    else{
      filtered$Instance_type = "Novel"
    }

    if (all(filtered$Instance_type == "Novel")){
      classified = qslim_insts[i,]
      classified$instance_type = "Novel"
    }
    else{
      classified = qslim_insts[i,]
      classified$instance_type = "New_Instance"
    }

    if(nrow(plyr::match_df(qslim_insts[i,], known_insts, on = c("uniprot", "Start_Pos", "End_Pos"))) != 0){
      classified$instance_type = "Known_Instance"
    }
    compiled[[i]] = classified
    setTxtProgressBar(pb = pb, i)
  }
  compiled = dplyr::bind_rows(compiled)
  return(compiled)
}

length_motif_regex = function(mot_regex){
  counter = 0
  non_wild_counter = 0
  splitted = unlist(strsplit(mot_regex, ""))
  if (splitted[1] == "^"){
    counter = counter + 1
    non_wild_counter = non_wild_counter + 1
    splitted = splitted[-1]
  }
  if (splitted[length(splitted)] == "$"){
    counter = counter + 1
    non_wild_counter = non_wild_counter + 1
    splitted = splitted[-length(splitted)]
  }

  if (any(splitted == "(")){
    to_remove_curved = c()
    curved_1 = which(splitted == "(")
    curved_2 = which(splitted == ")")
    curved_substr = c()
    for (i in 1:length(curved_1)){
      curved_substr = c(curved_substr, paste(splitted[c((curved_1[i] + 1):(curved_2[i] - 1))],collapse = ""))
      to_remove_curved = c(to_remove_curved, c(curved_1[i]:curved_2[i]))
    }

    for (i in curved_substr){
      counter = counter + length_motif_regex(i)[["total"]]
      non_wild_counter = non_wild_counter + length_motif_regex(i)[["non_wild"]]
    }
    splitted = splitted[-to_remove_curved]
  }
  if (any(splitted == "{")){
    to_remove_curly = c()
    curly_1 = which(splitted == "{")
    curly_2 = which(splitted == "}")
    for (i in 1:length(curly_1)){
      if (splitted[curly_1[i] - 1] == "]"){
        counter = counter + (as.numeric(splitted[curly_1[i] + 1]) -1)
        non_wild_counter = non_wild_counter + (as.numeric(splitted[curly_1[i] + 1]) -1)
        to_remove_curly = c(to_remove_curly,c(curly_1[i]:curly_2[i]))
      }
      else if (splitted[curly_1[i] - 1] == "."){
        counter = counter + as.numeric(splitted[curly_1[i] + 1])
        to_remove_curly = c(to_remove_curly,c((curly_1[i] - 1):curly_2[i]))
      }
      else{
        counter = counter + as.numeric(splitted[curly_1[i] + 1])
        non_wild_counter = non_wild_counter + as.numeric(splitted[curly_1[i] + 1])
        to_remove_curly = c(to_remove_curly,c((curly_1[i] - 1):curly_2[i]))
      }
    }
    splitted = splitted[-to_remove_curly]
  }

  if (any(splitted == "[")){
    to_remove_squared = c()
    squared_1 = which(splitted == "[")
    squared_2 = which(splitted == "]")
    for (i in 1:length(squared_1)){
      counter = counter + 1
      non_wild_counter = non_wild_counter + 1
      to_remove_squared = c(to_remove_squared, c(squared_1[i]:squared_2[i]))
    }
    splitted = splitted[-to_remove_squared]
  }
  counter = counter + length(splitted)
  non_wild_counter = non_wild_counter + length(splitted[which(splitted != ".")])
  return(list("total" = counter, "non_wild" = non_wild_counter))
}
