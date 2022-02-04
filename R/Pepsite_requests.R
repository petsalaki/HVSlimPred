Generate_pepsite_requests = function(new_hits_ints, EM_list, pfam_doms_list, pfam_pdb, human_human_int, host_viral = T){
  H1_pred_motifs = new_hits_ints %>% dplyr::select(Dataset,uniprot, Match, Org, Start_Pos, End_Pos)

  H1_pred_motifs = H1_pred_motifs[!duplicated(H1_pred_motifs),]
  if (host_viral){
    H1_pred_motifs = H1_pred_motifs %>% tidyr::separate(Dataset, into = c("V1", "H1"), sep = "_")
  }
  else{
    colnames(H1_pred_motifs)[1] = "H1"
  }
  H1_pred_motifs = H1_pred_motifs[which(H1_pred_motifs$H1 %in% pfam_pdb$uniprot),]
  H1_pred_motifs = dplyr::left_join(H1_pred_motifs, pfam_pdb, by = c("H1" = "uniprot"))
  H1_pred_motifs = H1_pred_motifs[!duplicated(H1_pred_motifs),]
  unique_H1_pfam = H1_pred_motifs %>% dplyr::select(H1, pfam_id)
  unique_H1_pfam = unique_H1_pfam[!duplicated(unique_H1_pfam),]
  unique_H1_pfam$is_H1 = "Yes"
  colnames(unique_H1_pfam)[1] = "uniprot"

  pf_nr = dplyr::bind_rows(EM_list$pfam_doms_nr)
  pf_nr = pf_nr %>% tidyr::separate(PFAM, into = c("pfam_id", "pfam_name"), sep = "--")
  pf_nr = pf_nr %>% dplyr::select(uniprot, pfam_id, pfam_name)

  pred_motifs = new_hits_ints %>% dplyr::select(uniprot, Match, Org, Start_Pos, End_Pos)
  pred_motifs = pred_motifs[!duplicated(pred_motifs),]
  pred_motifs = pred_motifs[which(pred_motifs$uniprot %in% pf_nr$uniprot),]

  pf_nr = dplyr::left_join(pf_nr, pred_motifs)
  pf_nr = pf_nr[which(pf_nr$pfam_id %in% pfam_pdb$pfam_id),] # remove domains with no available structures
  pred_motifs = pred_motifs[which(pred_motifs$uniprot %in% pf_nr$uniprot),]
  unique_nr = pf_nr[,c(1,2)]
  unique_nr = unique_nr[!duplicated(unique_nr),]
  unique_nr$is_H1 = "No"

  All_unique = rbind.data.frame(unique_H1_pfam, unique_nr)

  mapped_H1 = pfam2pdb_all(df = All_unique[which(All_unique$is_H1 == "Yes"),], pfam_pdb = pfam_pdb, pfam_doms_list = pfam_doms_list, H1_prot = T)
  mapped_not_H1 = pfam2pdb_all(df = All_unique[which(All_unique$is_H1 == "No"),], pfam_pdb = pfam_pdb, human_human_int = human_human_int, pfam_doms_list = pfam_doms_list, H1_prot = F)
  mapped_pdbs = rbind.data.frame(mapped_H1, mapped_not_H1)
  # mapped_pdbs = list()
  # for (i in 1:nrow(All_unique)){
  #   if (All_unique$is_H1[i] == "Yes"){
  #     result = pfam2pdb_per_mot(pfam_pdb = pfam_pdb, uniprot_id = All_unique$uniprot[i], pfam_ID = All_unique$pfam_id[i],H1_prot = T)
  #   }
  #   else{
  #     result = pfam2pdb_per_mot(pfam_pdb = pfam_pdb, uniprot_id = All_unique$uniprot[i], pfam_ID = All_unique$pfam_id[i], human_human_int = human_human_int, H1_prot = F)
  #   }
  #   mapped_pdbs[[i]] = result
  # }
  # mapped_pdbs = dplyr::bind_rows(mapped_pdbs)

  unique_pdbs = unique(mapped_pdbs$PDB)
  messed_up = unique_pdbs[which(nchar(unique_pdbs) > 4)]
  fixed_messed = c()
  for (i in 1:length(messed_up)){
    fixed_messed[i] = fix_messed_up_pdbIds(messed_up[i])
  }

  fixed_mapping = data.frame(old= messed_up, new = fixed_messed)
  mapped_pdbs = dplyr::left_join(mapped_pdbs, fixed_mapping, by = c("PDB" = "old"))
  mapped_pdbs$new[is.na(mapped_pdbs$new)] = mapped_pdbs$PDB[is.na(mapped_pdbs$new)]
  mapped_pdbs = mapped_pdbs %>% dplyr::select(-PDB)
  colnames(mapped_pdbs)[which(names(mapped_pdbs) == "new")] = "PDB"
  mapped_pdbs = mapped_pdbs[,c(ncol(mapped_pdbs), 1:(ncol(mapped_pdbs) - 1))]

  H1_pred_motifs = dplyr::left_join(H1_pred_motifs, fixed_mapping, by = c("PDB" = "old"))
  H1_pred_motifs$new[is.na(H1_pred_motifs$new)] = H1_pred_motifs$PDB[is.na(H1_pred_motifs$new)]
  H1_pred_motifs = H1_pred_motifs %>% dplyr::select(-PDB)
  colnames(H1_pred_motifs)[which(names(H1_pred_motifs) == "new")] = "PDB"

  mapped_pdbs = Annotate_pdb_API(mapped_df = mapped_pdbs)

  x = mapped_pdbs %>%  dplyr::select(domain_protein, pfam_id, motif_protein)
  x = x[!duplicated(x),]
  final_selected = list()
  for (i in 1:nrow(x)){
    if (is.na(x$motif_protein[i])){
      filtered = dplyr::filter(mapped_pdbs, domain_protein == x$domain_protein[i], pfam_id == x$pfam_id[i])
      filtered = filtered[is.na(filtered$motif_protein),]
    }
    else{
      filtered = dplyr::filter(mapped_pdbs, domain_protein == x$domain_protein[i], pfam_id == x$pfam_id[i], motif_protein == x$motif_protein[i])
    }
    final_selected[[i]] = Best_PDB(filtered, priortise_ligand = T)
  }
  final_selected = dplyr::bind_rows(final_selected)

  pdb_info = dplyr::select(final_selected, PDB, Chain, PDB_chain, ligandId, experimentalTechnique, resolution)
  pdb_info = pdb_info[!duplicated(pdb_info),]

  colnames(H1_pred_motifs)[which(names(H1_pred_motifs) == "H1")] = "domain_protein"
  H1_final_requests = plyr::match_df(H1_pred_motifs, final_selected)
  H1_final_requests = H1_final_requests[!is.na(H1_final_requests$Match),]
  colnames(H1_final_requests)[which(names(H1_final_requests) == "uniprot")] = "motif_protein"
  H1_final_requests = dplyr::left_join(H1_final_requests, pdb_info)
  if (host_viral){
    H1_final_requests = H1_final_requests %>% dplyr::select(-V1)
  }
  colnames(H1_final_requests)[which(names(H1_final_requests) == "Match")] = "peptide"

  x = H1_final_requests %>% dplyr::select(PDB, Chain, peptide,Start_Pos, End_Pos, motif_protein)

  final_requests = dplyr::left_join(final_selected[!is.na(final_selected$motif_protein),], pred_motifs, by = c("motif_protein" = "uniprot"))
  final_requests = final_requests[!is.na(final_requests$Match),]
  colnames(final_requests)[which(names(final_requests) == "Match")] = "peptide"

  y = final_requests %>% dplyr::select(PDB, Chain, peptide,Start_Pos, End_Pos, motif_protein)

  Total_final_requests = rbind.data.frame(y, x)
  Total_final_requests = Total_final_requests[!duplicated(Total_final_requests),]
  return(list("Total_requests" = Total_final_requests, "EM_requests" = final_requests, "H1_requests" = H1_final_requests))
}

Generate_pepsite_AF_reqs = function(new_hits_ints, EM_list, pfam_doms_list, pfam_human, pfam_pdb, human_human_int, host_viral = T,
                                    PLDDT = 70, plddt_overlap_thresh = 0.5){
  AF_pfam_human = pfam_human %>% dplyr::select(Uniprot.ID, envelope.start, envelope.end, hmm.acc)
  colnames(AF_pfam_human) = c("uniprot", "pdb_start", "pdb_end", "pfam_id")
  AF_pfam_human$uniprot_start = AF_pfam_human$pdb_start
  AF_pfam_human$uniprot_end = AF_pfam_human$pdb_end
  AF_pfam_human$PDB = paste0("AF_", AF_pfam_human$uniprot)
  AF_pfam_human$Chain = "A"
  AF_pfam_human = AF_pfam_human[,colnames(pfam_pdb)]
  AF_pfam_human = AF_pfam_human[which(AF_pfam_human$uniprot %in% human_AF_pdb_names$uniprot),]

  H1_pred_motifs = new_hits_ints %>% dplyr::select(Dataset,uniprot, Match, Org, Start_Pos, End_Pos)
  H1_pred_motifs = H1_pred_motifs[!duplicated(H1_pred_motifs),]
  if (host_viral){
    H1_pred_motifs = H1_pred_motifs %>% tidyr::separate(Dataset, into = c("V1", "H1"), sep = "_")
  }
  else{
    colnames(H1_pred_motifs)[1] = "H1"
  }
  H1_pred_motifs = H1_pred_motifs[which(H1_pred_motifs$H1 %nin% pfam_pdb$uniprot),]
  H1_pred_motifs = H1_pred_motifs[which(H1_pred_motifs$H1 %in% AF_pfam_human$uniprot),]

  H1_pred_motifs = dplyr::left_join(H1_pred_motifs, AF_pfam_human, by = c("H1" = "uniprot"))
  H1_pred_motifs = H1_pred_motifs[!duplicated(H1_pred_motifs),]
  unique_H1_pfam = H1_pred_motifs %>% dplyr::select(H1, pfam_id)
  unique_H1_pfam = unique_H1_pfam[!duplicated(unique_H1_pfam),]
  unique_H1_pfam$is_H1 = "Yes"
  colnames(unique_H1_pfam)[1] = "uniprot"

  pf_nr = dplyr::bind_rows(EM_list$pfam_doms_nr)
  pf_nr = pf_nr %>% tidyr::separate(PFAM, into = c("pfam_id", "pfam_name"), sep = "--")
  pf_nr = pf_nr %>% dplyr::select(uniprot, pfam_id, pfam_name)

  pred_motifs = new_hits_ints %>% dplyr::select(uniprot, Match, Org, Start_Pos, End_Pos)
  pred_motifs = pred_motifs[!duplicated(pred_motifs),]
  pred_motifs = pred_motifs[which(pred_motifs$uniprot %in% pf_nr$uniprot),]

  pf_nr = dplyr::left_join(pf_nr, pred_motifs)
  pf_nr = pf_nr[which(pf_nr$pfam_id %nin% pfam_pdb$pfam_id),] # remove domains with available structures
  pf_nr = pf_nr[which(pf_nr$pfam_id %in% AF_pfam_human$pfam_id),] # Include domains with alphafold available structures
  pred_motifs = pred_motifs[which(pred_motifs$uniprot %in% pf_nr$uniprot),]
  unique_nr = pf_nr[,c(1,2)]
  unique_nr = unique_nr[!duplicated(unique_nr),]
  unique_nr$is_H1 = "No"

  All_unique = rbind.data.frame(unique_H1_pfam, unique_nr)

  mapped_H1 = pfam2pdb_AF(df = All_unique[which(All_unique$is_H1 == "Yes"),], alphafold_pfam = AF_pfam_human, pfam_doms_list = pfam_doms_list, H1_prot = T)
  mapped_not_H1 = pfam2pdb_AF(df = All_unique[which(All_unique$is_H1 == "No"),], alphafold_pfam = AF_pfam_human, human_human_int = human_human_int, pfam_doms_list = pfam_doms_list, H1_prot = F)
  mapped_pdbs = rbind.data.frame(mapped_H1, mapped_not_H1)

  mapped_pdbs$PDB_chain = paste0(mapped_pdbs$PDB, ".", mapped_pdbs$Chain)
  mapped_pdbs$ligandId = NA
  mapped_pdbs$experimentalTechnique = "Prediction"
  mapped_pdbs$resolution = NA

  x = mapped_pdbs %>%  dplyr::select(domain_protein, pfam_id)
  x = x[!duplicated(x),]
  final_selected = list()
  message("Filtering domains with more than ", plddt_overlap_thresh * 100, " % overlap with high PLDDT ...")
  pb = txtProgressBar(min = 0, max = nrow(x), initial = 0, style = 3)
  for (i in 1:nrow(x)){
    filtered = dplyr::filter(mapped_pdbs, domain_protein == x$domain_protein[i], pfam_id == x$pfam_id[i])
    pass_plddt = AF_pLDDT_filter(df = filtered, PLDDT = PLDDT, plddt_overlap_thresh = plddt_overlap_thresh)
    final_selected[[i]] = pass_plddt
    setTxtProgressBar(pb = pb, i)
  }
  final_selected = dplyr::bind_rows(final_selected)

  pdb_info = dplyr::select(final_selected, PDB, Chain, PDB_chain, ligandId, experimentalTechnique, resolution, AF_PLDDT, AF_dom_overlap)
  pdb_info = pdb_info[!duplicated(pdb_info),]

  colnames(H1_pred_motifs)[which(names(H1_pred_motifs) == "H1")] = "domain_protein"
  H1_final_requests = plyr::match_df(H1_pred_motifs, final_selected)
  H1_final_requests = H1_final_requests[!is.na(H1_final_requests$Match),]
  colnames(H1_final_requests)[which(names(H1_final_requests) == "uniprot")] = "motif_protein"
  H1_final_requests = dplyr::left_join(H1_final_requests, pdb_info)
  if (host_viral){
    H1_final_requests = H1_final_requests %>% dplyr::select(-V1)
  }
  colnames(H1_final_requests)[which(names(H1_final_requests) == "Match")] = "peptide"

  x = H1_final_requests %>% dplyr::select(PDB, Chain, peptide,Start_Pos, End_Pos, motif_protein)

  final_requests = dplyr::left_join(final_selected[!is.na(final_selected$motif_protein),], pred_motifs, by = c("motif_protein" = "uniprot"))
  final_requests = final_requests[!is.na(final_requests$Match),]
  colnames(final_requests)[which(names(final_requests) == "Match")] = "peptide"

  y = final_requests %>% dplyr::select(PDB, Chain, peptide,Start_Pos, End_Pos, motif_protein)

  Total_final_requests = rbind.data.frame(y, x)
  Total_final_requests = Total_final_requests[!duplicated(Total_final_requests),]
  return(list("Total_requests" = Total_final_requests, "EM_requests" = final_requests, "H1_requests" = H1_final_requests))
}

Best_PDB = function(df,priortise_ligand = T){
  filtered = df
  if (priortise_ligand){
    if (any(is.na(filtered$ligandId))){
      filtered = filtered[is.na(filtered$ligandId),]
    }
    exp_type = unique(filtered$experimentalTechnique)
    if ("X-RAY DIFFRACTION" %in% exp_type || "ELECTRON MICROSCOPY" %in% exp_type || "ELECTRON CRYSTALLOGRAPHY" %in% exp_type){
      filtered = filtered[!is.na(filtered$resolution),]
      min_resol = min(filtered$resolution)
      filtered = filtered[which(filtered$resolution == min_resol),]
      return(filtered[1,])
    }
    else{
      return(filtered[1,])
    }
  }
  else{
    exp_type = unique(filtered$experimentalTechnique)
    if ("X-RAY DIFFRACTION" %in% exp_type || "ELECTRON MICROSCOPY" %in% exp_type){
      filtered = filtered[!is.na(filtered$resolution),]
      min_resol = min(filtered$resolution)
      filtered = filtered[which(filtered$resolution == min_resol),]

      if (any(is.na(filtered$ligandId))){
        filtered = filtered[is.na(filtered$ligandId),]
        return(filtered[1,])
      }
      else{
        return(filtered[1,])
      }
    }
    else{
      return(filtered[1,])
    }
  }
}
pfam2pdb_per_mot = function(pfam_pdb,uniprot_id, pfam_ID, human_human_int, H1_prot = F){
  alphabet = c(LETTERS, letters)
  if (H1_prot & missing(human_human_int)){
    pdbs = dplyr::filter(pfam_pdb, uniprot == uniprot_id, pfam_id == pfam_ID)
    if (length(intersect(pdbs$Chain, alphabet)) > 0){
      pdbs = pdbs[which(pdbs$Chain %in% alphabet),]
      pdbs = pdbs[order(match(pdbs$Chain, alphabet)),]
      pdbs = pdbs[!duplicated(pdbs$PDB),]
    }
    else{
      pdbs = pdbs[order(pdbs$Chain),]
      pdbs = pdbs[!duplicated(pdbs$PDB),]
    }
    colnames(pdbs)[3] = "domain_protein"
    pdbs$is_domain_prot_intr = "Yes"
    pdbs$motif_protein = NA
    return(pdbs)
  }
  else{
    ints_mots = get_human_interactors(human_human_intact = human_human_int, uniprot_id = uniprot_id, return_fasta = F)
    pfam_structures = pfam_pdb[which(pfam_pdb$pfam_id == pfam_ID),]
    if (nrow(pfam_structures) == 0){
      return(NULL)
    }

    if (length(intersect(pfam_structures$uniprot , ints_mots)) > 0){
      common_prots = intersect(pfam_structures$uniprot, ints_mots)
      pdbs = pfam_structures[which(pfam_structures$uniprot %in% common_prots),]
      if (length(intersect(pdbs$Chain, alphabet)) > 0){
        pdbs = pdbs[which(pdbs$Chain %in% alphabet),]
        pdbs = pdbs[order(match(pdbs$Chain, alphabet)),]
        pdbs = pdbs[!duplicated(pdbs$PDB),]
      }
      else{
        pdbs = pdbs[order(pdbs$Chain),]
        pdbs = pdbs[!duplicated(pdbs$PDB),]
      }
      colnames(pdbs)[3] = "domain_protein"
      pdbs$is_domain_prot_intr = "Yes"
      pdbs$motif_protein = uniprot_id


      return(pdbs)
    }
    else{
      pdbs = pfam_structures
      if (length(intersect(pdbs$Chain, alphabet)) > 0){
        pdbs = pdbs[which(pdbs$Chain %in% alphabet),]
        pdbs = pdbs[order(match(pdbs$Chain, alphabet)),]
        pdbs = pdbs[!duplicated(pdbs$PDB),]
      }
      else{
        pdbs = pdbs[order(pdbs$Chain),]
        pdbs = pdbs[!duplicated(pdbs$PDB),]
      }

      colnames(pdbs)[3] = "domain_protein"
      pdbs$is_domain_prot_intr = "No"
      pdbs$motif_protein = uniprot_id

      return(pdbs)
    }
  }
}
pfam2pdb_all = function(df, pfam_pdb, human_human_int, pfam_doms_list, H1_prot = F){
  alphabet = c(LETTERS, letters)
  pfam_doms_list = pfam_doms_list %>% tidyr::separate(PFAM, into = c("pfam_id", "pfam_name"), sep = "--") %>% dplyr::select(uniprot, pfam_id)
  if (H1_prot & missing(human_human_int)){
    pdbs = dplyr::filter(pfam_pdb, pfam_id %in% df$pfam_id, uniprot %in% df$uniprot)

    in_pfam_alphabet = pdbs[which(pdbs$Chain %in% alphabet),]
    in_pfam_alphabet = in_pfam_alphabet[order(match(in_pfam_alphabet$Chain, alphabet)),]
    in_pfam_alphabet = in_pfam_alphabet[!duplicated(in_pfam_alphabet[,c(which(names(in_pfam_alphabet) %in% c("PDB", "pfam_id", "uniprot")))]),]

    in_pfam_not_alphabet = pdbs[which(pdbs$Chain %nin% alphabet),]
    in_pfam_not_alphabet = in_pfam_not_alphabet[order(in_pfam_not_alphabet$Chain),]
    in_pfam_not_alphabet = in_pfam_not_alphabet[!duplicated(in_pfam_not_alphabet[,c(which(names(in_pfam_not_alphabet) %in% c("PDB", "pfam_id", "uniprot")))]),]
    colnames(pdbs)[3] = "domain_protein"
    pdbs$motif_protein = NA
    return(pdbs)
  }
  else{
    ints_mots = list()
    interactors = list()
    unique_prots = unique(df$uniprot)
    for (i in 1:length(unique_prots)){
      ints = get_human_interactors(human_human_intact = human_human_int, uniprot_id = unique_prots[i], return_fasta = F)
      interactors[[i]] = ints
      ints_mots[[i]] = data.frame(motif_protein = unique_prots[i], domain_protein = ints)
    }
    interactors = unlist(interactors)
    ints_mots = dplyr::bind_rows(ints_mots)
    ints_mots = dplyr::left_join(ints_mots, pfam_doms_list, by = c("domain_protein" = "uniprot"))
    colnames(ints_mots) = c("uniprot", "domain_protein", "pfam_id")
    ints_mots = plyr::match_df(ints_mots, df)
    colnames(ints_mots) = c("motif_protein", "uniprot", "pfam_id")
    ints_mots[!duplicated(ints_mots),]

    pfam_structures = pfam_pdb[which(pfam_pdb$pfam_id %in% df$pfam_id),]
    #pfam_structures = dplyr::left_join(pfam_structures, ints_mots)
    ints_mots_in_pfam = pfam_structures[which(pfam_structures$uniprot %in% interactors),]
    #ints_mots_not_in_pfam = pfam_structures[which(pfam_structures$uniprot %nin% interactors),]

    in_pfam_alphabet = ints_mots_in_pfam[which(ints_mots_in_pfam$Chain %in% alphabet),]
    in_pfam_alphabet = in_pfam_alphabet[order(match(in_pfam_alphabet$Chain, alphabet)),]
    in_pfam_alphabet = in_pfam_alphabet[!duplicated(in_pfam_alphabet[,c(which(names(in_pfam_alphabet) %in% c("PDB", "pfam_id", "uniprot")))]),]
    in_pfam_alphabet = dplyr::left_join(in_pfam_alphabet, ints_mots)


    in_pfam_not_alphabet = ints_mots_in_pfam[which(ints_mots_in_pfam$Chain %nin% alphabet),]
    in_pfam_not_alphabet = in_pfam_not_alphabet[order(in_pfam_not_alphabet$Chain),]
    in_pfam_not_alphabet = in_pfam_not_alphabet[!duplicated(in_pfam_not_alphabet[,c(which(names(in_pfam_not_alphabet) %in% c("PDB", "pfam_id", "uniprot")))]),]
    in_pfam_not_alphabet = dplyr::left_join(in_pfam_not_alphabet, ints_mots)

    # not_in_pfam_alphabet = ints_mots_not_in_pfam[which(ints_mots_not_in_pfam$Chain %in% alphabet),]
    # not_in_pfam_alphabet = not_in_pfam_alphabet[order(match(not_in_pfam_alphabet$Chain, alphabet)),]
    # not_in_pfam_alphabet = not_in_pfam_alphabet[!duplicated(not_in_pfam_alphabet[,c(which(names(not_in_pfam_alphabet) %in% c("PDB", "pfam_id", "uniprot")))]),]
    # #not_in_pfam_alphabet = dplyr::left_join(not_in_pfam_alphabet, ints_mots, by = c("uniprot" = "domain_protein"))
    # not_in_pfam_alphabet$motif_protein = NA
    #
    # not_in_pfam_not_alphabet = ints_mots_not_in_pfam[which(ints_mots_in_pfam$Chain %nin% alphabet),]
    # not_in_pfam_not_alphabet = not_in_pfam_not_alphabet[order(not_in_pfam_not_alphabet$Chain),]
    # not_in_pfam_not_alphabet = not_in_pfam_not_alphabet[!duplicated(not_in_pfam_not_alphabet[,c(which(names(not_in_pfam_not_alphabet) %in% c("PDB", "pfam_id", "uniprot")))]),]
    # #not_in_pfam_not_alphabet = dplyr::left_join(not_in_pfam_not_alphabet, ints_mots, by = c("uniprot" = "domain_protein"))
    # not_in_pfam_not_alphabet$motif_protein = NA

    pdbs = rbind.data.frame(in_pfam_alphabet, in_pfam_not_alphabet)
    colnames(pdbs)[3] = "domain_protein"
    #pdbs$is_domain_prot_intr = "Yes"
    #pdbs = dplyr::left_join(pdbs, ints_mots, by = c("domain_protein", "pfam_id"))
    pdbs = pdbs[!is.na(pdbs$motif_protein),]
    pdbs = pdbs[!duplicated(pdbs),]
    return(pdbs)
  }
}
pfam2pdb_AF = function(df, alphafold_pfam, human_human_int, pfam_doms_list, H1_prot = F){
  pfam_doms_list = pfam_doms_list %>% tidyr::separate(PFAM, into = c("pfam_id", "pfam_name"), sep = "--") %>% dplyr::select(uniprot, pfam_id)
  if (H1_prot & missing(human_human_int)){
    pdbs = dplyr::filter(alphafold_pfam, pfam_id %in% df$pfam_id, uniprot %in% df$uniprot)
    colnames(pdbs)[3] = "domain_protein"
    pdbs$motif_protein = NA
    return(pdbs)
  }
  else{
    ints_mots = list()
    interactors = list()
    unique_prots = unique(df$uniprot)
    for (i in 1:length(unique_prots)){
      ints = get_human_interactors(human_human_intact = human_human_int, uniprot_id = unique_prots[i], return_fasta = F)
      interactors[[i]] = ints
      ints_mots[[i]] = data.frame(motif_protein = unique_prots[i], domain_protein = ints)
    }
    interactors = unlist(interactors)
    ints_mots = dplyr::bind_rows(ints_mots)
    ints_mots = dplyr::left_join(ints_mots, pfam_doms_list, by = c("domain_protein" = "uniprot"))
    colnames(ints_mots) = c("uniprot", "domain_protein", "pfam_id")
    ints_mots = plyr::match_df(ints_mots, df)
    colnames(ints_mots) = c("motif_protein", "uniprot", "pfam_id")
    ints_mots[!duplicated(ints_mots),]

    pfam_structures = alphafold_pfam[which(alphafold_pfam$pfam_id %in% df$pfam_id),]
    #pfam_structures = dplyr::left_join(pfam_structures, ints_mots)
    ints_mots_in_pfam = pfam_structures[which(pfam_structures$uniprot %in% interactors),]
    #ints_mots_not_in_pfam = pfam_structures[which(pfam_structures$uniprot %nin% interactors),]


    pdbs = dplyr::left_join(ints_mots_in_pfam, ints_mots)
    colnames(pdbs)[3] = "domain_protein"
    #pdbs$is_domain_prot_intr = "Yes"
    #pdbs = dplyr::left_join(pdbs, ints_mots, by = c("domain_protein", "pfam_id"))
    pdbs = pdbs[!is.na(pdbs$motif_protein),]
    pdbs = pdbs[!duplicated(pdbs),]
    return(pdbs)
  }
}


fix_messed_up_pdbIds = function(pdbId){
  if (nchar(pdbId) == 5){
    return(toupper(stringr::str_remove(pdbId, "-")))
  }
  if (nchar(pdbId) == 8){
    trimmed = paste0(substr(pdbId, 1,1), "E", substr(pdbId, 7,8))
    return(trimmed)
  }
}
Annotate_pdb_API = function(mapped_df){
  mapped_df$PDB_chain = paste0(mapped_df$PDB, ".", mapped_df$Chain)

  unique_pdbs = unique(mapped_df$PDB)
  unique_pdbs_chains = unique(mapped_df$PDB_chain)

  qry <- ghql::Query$new()
  qry$query('pdb_info', paste0('{
  entries(entry_ids: ',paste0("[",paste(dQuote(unique_pdbs, q = F), collapse = ","),"]"),
                               ') {
    entry {
      id
    }
    exptl {
      method
    }
    pdbx_vrpt_summary {
      PDB_resolution
    }
  }
}'))
  qry$query('pdb_chain_info', paste0('{
  polymer_entity_instances(instance_ids: ',paste0("[",paste(dQuote(unique_pdbs_chains, q = F), collapse = ","),"]"),
                                     ') {
    rcsb_id
    rcsb_polymer_instance_feature {
      assignment_version
      description
      feature_id
      name
      provenance_source
      reference_scheme
      type
    }
  }
}'))

  con <- ghql::GraphqlClient$new('https://data.rcsb.org/graphql')
  pdb_info <- con$exec(qry$queries$pdb_info)
  pdb_info = jsonlite::fromJSON(pdb_info)
  pdb_info = as.data.frame(pdb_info$data$entries)
  pdb_info$entry = unlist(pdb_info$entry)
  pdb_info$exptl = unlist(lapply(pdb_info$exptl, function(x) x[1,]))
  pdb_info[,3] = unlist(pdb_info[,3])
  colnames(pdb_info) = c("structureId", "experimentalTechnique", "resolution")

  pdb_chain_info <- con$exec(qry$queries$pdb_chain_info)
  pdb_chain_info = jsonlite::fromJSON(pdb_chain_info)
  pdb_chain_info = clean_pdb_chain_info_json(pdb_chain_info)
  pdb_chain_info = pdb_chain_info %>% tidyr::separate(PDB_chain, into = c("structureId", "chainId"), sep = "[.]")

  pdb_annot = dplyr::left_join(pdb_chain_info, pdb_info)
  colnames(pdb_annot)[1:2] = c("PDB", "Chain")
  pdb_annot = plyr::match_df(pdb_annot, mapped_df, on = c("PDB", "Chain"))

  mapped_df = dplyr::left_join(mapped_df, pdb_annot)
  return(mapped_df)
}
clean_pdb_chain_info_json = function(res){
  final_result = list()
  for (i in 1:length(res$data$polymer_entity_instances$rcsb_id)){
    annots = res$data$polymer_entity_instances$rcsb_polymer_instance_feature[[i]]
    if ("BINDING_SITE" %in% annots[,7]){
      ligands = unique(annots[which(annots[,7] == "BINDING_SITE"),4])
      ligandIds = c()
      for (j in 1:length(ligands)){
        ligandIds[j] = strsplit(ligands[j], " ")[[1]][2]
      }
      ligandIds = paste(ligandIds, collapse = ",")
      final_result[[i]] = data.frame(PDB_chain = res$data$polymer_entity_instances$rcsb_id[i],
                                     ligandId = ligandIds)
    }
    else{
      final_result[[i]] = data.frame(PDB_chain = res$data$polymer_entity_instances$rcsb_id[i],
                                     ligandId = NA)
    }
  }
  final_result = dplyr::bind_rows(final_result)
  return(final_result)
}
Create_peptide_overlaps = function(requests){
  more_than_10 = requests[which(nchar(requests$peptide) > 10),]
  less_than_10 = mismatch_df(requests, more_than_10)

  overlaps = list()
  for (i in 1:nrow(more_than_10)){
    result = more_than_10[i,]
    result = result[rep(seq_len(nrow(result)), each = 2),]
    result$peptide[1] = substr(result$peptide[1], 1, 10)
    result$peptide[2] = substr(result$peptide[2], switch(as.character(nchar(result$peptide[2])),
                                                         "11" = 2, "12" = 3, "13" = 4), nchar(result$peptide[2]))
    overlaps[[i]] = result
  }
  overlaps = dplyr::bind_rows(overlaps)

  final_requests = rbind.data.frame(less_than_10, overlaps)
}

get_AF_pdb = function(uniprot_id){
  if(uniprot_id %in% human_AF_pdb_names$uniprot){
    pdb_file_path = human_AF_pdb_names$pdb_filename[which(human_AF_pdb_names$uniprot == uniprot_id)]
    if (length(pdb_file_path) > 1){
      return(NULL)
    }
    tmp = tempfile()
    y = paste0("https://alphafold.ebi.ac.uk/files/", pdb_file_path)
    data = download.file(y, destfile = tmp, quiet = T)
    pdb_res = bio3d::read.pdb(tmp)
    final = pdb_res$atom
    return(final)
  }
  else{
    return(NULL)
  }
}
AF_pLDDT_filter = function(df, PLDDT = 70, plddt_overlap_thresh = 0.5){
  upid = unique(df$domain_protein)
  pdb_atoms = get_AF_pdb(uniprot_id = upid)
  if (is.null(pdb_atoms)){
    return(NULL)
  }

  final_res = list()
  for (i in 1:nrow(df)){
    filt = pdb_atoms[which(pdb_atoms$resno %in% c(df$uniprot_start[i]:df$uniprot_end[i])),]
    conf_per_res = filt %>% dplyr::select(resno, b)
    conf_per_res = conf_per_res[!duplicated(conf_per_res),]
    dom_len = length(c(df$uniprot_start[i]:df$uniprot_end[i]))
    pass_conf = conf_per_res[which(conf_per_res$b >= PLDDT),]
    if((length(unique(pass_conf$resno)) / dom_len) >= plddt_overlap_thresh){
      res = df[i,]
      res$AF_PLDDT = mean(conf_per_res$b)
      res$AF_dom_overlap = length(unique(pass_conf$resno)) / dom_len
      final_res[[i]] = res
    }
    else{
      final_res[[i]] = NULL
    }
  }
  if (length(unlist(final_res)) == 0){
    return(NULL)
  }
  else{
    return(dplyr::bind_rows(final_res))
  }
}

#trial = readRDS("../Slimfinder_human_only/Package_analysis/Complete_Evaluation.rds")

#host_viral = readRDS("../New_Analysis_everything/Complete_Evaluation.rds")


