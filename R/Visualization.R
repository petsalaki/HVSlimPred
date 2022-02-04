Heatmap_all_evaluation = function(Eval_list, benchtype = c("ELM", "PRMdb")){
  heatmap_eval_ELM = data.frame(pfam_filt = pfam_dom_enrich_info$pfam_filter, OR_prot = 0,
                                RC_prot = 0, PR_prot = 0, F0.5_prot_dom = 0, RC_prot_dom = 0,
                                PR_prot_dom = 0,F0.5_HH = 0, RC_HH = 0, PR_HH = 0)

  colors = colorRampPalette(rev(RColorBrewer::brewer.pal(name = "RdYlBu", n = 5)))(100)
  for (i in 1:nrow(heatmap_eval_ELM)){
    heatmap_eval_ELM[i,-1] = map_all_eval_heat(Eval_list = Eval_list, heatmap_eval_ELM$pfam_filt[i], benchtype = benchtype)
  }
  #all_Eval_mat = heatmap_eval_ELM %>% dplyr::select(-RC_prot, -RC_HH, -RC_prot_dom)
  all_Eval_mat = heatmap_eval_ELM
  all_Eval_mat$pfam_filt = sub(".*PF_","",all_Eval_mat$pfam_filt)
  rownames(all_Eval_mat) = all_Eval_mat$pfam_filt
  all_Eval_mat = all_Eval_mat[,-1]
  for (i in 1:ncol(all_Eval_mat)){
    all_Eval_mat[,i] = all_Eval_mat[,i] / max(all_Eval_mat[,i])
  }
  #all_Eval_mat$Avg = rowMeans(all_Eval_mat[,c(1,2,3)])
  all_Eval_mat = as.matrix(all_Eval_mat)
  pheatmap::pheatmap(all_Eval_mat, cluster_rows = T, cluster_cols = F, cutree_rows = 10, color = colors, border_color = NA, main = paste0(benchtype, " Scores"))
}
map_all_eval_heat = function(Eval_list,pfam_filt, benchtype = c("ELM", "PRMdb")){
  final = rep(0,9)
  names(final) = c("OR_prot", "RC_prot", "PR_prot", "F0.5_prot_dom", "RC_prot_dom", "PR_prot_dom", "F0.5_HH", "RC_HH", "PR_HH")
  if (benchtype == "ELM"){
    HH = Eval_list[["ELM"]][["HH_motif"]][["All_F1_scores"]] %>% dplyr::filter(filter == pfam_filt, weight_site_resid == "residue site_wise")
    HH$F0.5 = (1.25 * HH$Rc * HH$Pr) / (0.25 * (HH$Rc + HH$Pr))
    prot = Eval_list[["ELM"]][["Prot_level"]] %>% dplyr::filter(filter == pfam_filt)
    prot_dom = Eval_list[["ELM"]][["Prot_dom"]] %>% dplyr::filter(filter == pfam_filt)
    prot_dom$F0.5 = (1.25 * prot_dom$recall * prot_dom$precision) / (0.25 * (prot_dom$recall + prot_dom$precision))
  }
  if (benchtype == "PRMdb"){
    HH = Eval_list[["PRMdb"]][["HH_motif"]][["All_F1_scores"]] %>% dplyr::filter(filter == pfam_filt, weight_site_resid == "residue site_wise")
    HH$F0.5 = (1.25 * HH$Rc * HH$Pr) / (0.25 * (HH$Rc + HH$Pr))
    prot = Eval_list[["PRMdb"]][["Prot_level"]] %>% dplyr::filter(filter == pfam_filt)
    prot_dom = Eval_list[["PRMdb"]][["Prot_dom"]] %>% dplyr::filter(filter == pfam_filt)
    prot_dom$F0.5 = (1.25 * prot_dom$recall * prot_dom$precision) / (0.25 * (prot_dom$recall + prot_dom$precision))
  }

  final["OR_prot"] = prot$odds_ratio
  final["RC_prot"] = prot$recall
  final["PR_prot"] = prot$precision
  final["F0.5_prot_dom"] = prot_dom$F0.5
  final["RC_prot_dom"] = prot_dom$recall
  final["PR_prot_dom"] = prot_dom$precision
  final["F0.5_HH"] = HH$F0.5
  final["RC_HH"] = HH$Rc
  final["PR_HH"] = HH$Pr
  return(as.numeric(final))
}

Heatmap_motif_metrics = function(metrics_per_mot, benchtype = c("ELM", "PRMdb"), return_non_zero = T){
  if (benchtype == "ELM"){
    HH_mat = metrics_per_mot
    rownames(HH_mat) = HH_mat$Id
    HH_mat = HH_mat[,-c(1,6)]
    HH_mat = as.matrix(HH_mat)
    x = pheatmap::pheatmap(HH_mat, cluster_rows = T, cluster_cols = T, height = 17, fontsize_row = 8, cutree_rows = 4, main = benchtype)
    if (return_non_zero){
      HH_mat = HH_mat[which(HH_mat[,2] != 0),]
      #y = cbind(HH_mat, cutree(x$tree_row, k = 4))
      #final = pheatmap::pheatmap(y[which(y[,8] != 1),][,-8], cluster_cols = T, main = benchtype)
      final = pheatmap::pheatmap(HH_mat, cluster_cols = T, cluster_rows = T, main = benchtype, treeheight_col = 5, treeheight_row = 5, cutree_rows = 4)
      return(final)
    }
    else{
      return(x)
    }
  }
  if (benchtype == "PRMdb"){
    HH_mat = metrics_per_mot
    HH_mat$combined_id = NA
    for(i in 1:nrow(HH_mat)){
      splitted = unlist(strsplit(HH_mat$Id[i], "_"))
      id_name = paste0(splitted[1], "_", splitted[2])
      filtered_info = PRM_info[which(PRM_info$Id == id_name),]
      if(length(splitted) > 2){
        combined = paste0(filtered_info$Protein_Name, "_", filtered_info$Domain_Group,"_", filtered_info$Domain_Number, ".", splitted[3])
      }
      else{
        combined = paste0(filtered_info$Protein_Name, "_", filtered_info$Domain_Group,"_", filtered_info$Domain_Number)
      }
      HH_mat$combined_id[i] = combined
    }
    HH_mat = HH_mat[!duplicated(HH_mat$combined_id),]

    rownames(HH_mat) = HH_mat$combined_id
    HH_mat = HH_mat[,-c(1,6,10)]
    HH_mat = as.matrix(HH_mat)
    x = pheatmap::pheatmap(HH_mat, cluster_rows = T, cluster_cols = T, height = 17, fontsize_row = 8, cutree_rows = 3, main = benchtype)
    if (return_non_zero){
      y = cbind(HH_mat, cutree(x$tree_row, k = 3))
      final = pheatmap::pheatmap(y[which(y[,8] != 1),][,-8], cluster_cols = T, main = benchtype)
      return(final)
    }
    else{
      return(x)
    }
  }

}

Sankey_ELM = function(new_hits_ints, fine_grained = T, compariMotif_cutoff = 0.6){
  data = new_hits_ints %>% dplyr::select(uniprot, Start_Pos, End_Pos, Id, new_score, Accession, motif_distance)
  colnames(data)[which(colnames(data) == "new_score")] = "Score"
  data = data[!duplicated(data),]
  data$ELM_class = substr(data$Id, 1,3)
  all_classes = unique(data$ELM_class)
  all_classes = all_classes[-7]
  if (!fine_grained){
    data$category = NA
    data[which(data$motif_distance == 0 & !is.na(data$Accession)),]$category = "Known Instances"
    data[which(data$Score >= compariMotif_cutoff & is.na(data$Accession)),]$category = "New Instances"
    data[which(data$Score < compariMotif_cutoff),]$category = "Novel Motif"

    links = data.frame(source = rep(all_classes,each = 3), target = rep(c("Known Instances", "New Instances", "Novel Motif"),6), value = 0)

    for (i in 1:nrow(links)){
      filtered = dplyr::filter(data, category == links$target[i], ELM_class == links$source[i])
      filtered = filtered[!duplicated(filtered[,c(1:3)]),]
      links$value[i] = nrow(filtered)
    }
    links = links[which(links$value != 0),]
    nodes <- data.frame(
      name=c(as.character(links$source),
             as.character(links$target)) %>% unique()
    )

    links$IDsource <- match(links$source, nodes$name)-1
    links$IDtarget <- match(links$target, nodes$name)-1

    # Make the Network
    p <- networkD3::sankeyNetwork(Links = links, Nodes = nodes,
                       Source = "IDsource", Target = "IDtarget",
                       Value = "value", NodeID = "name",
                       sinksRight= F, fontSize = 14)
    return(p)
  }
  else{
    data$category = NA
    data$conf = NA

    data[which(data$Score >= 0.67),]$conf = "Highest Similarity"
    data[which(data$Score >= 0.51 & data$Score < 0.67),]$conf = "High Similarity"
    data[which(data$Score >= 0.34 & data$Score < 0.51),]$conf = "Low Similarity"
    #data[which(data$Score >= 2 & data$Score < 2.5),]$conf = "Low Similarity"
    data[which(data$Score < 0.34),]$conf = "Lowest Similarity"

    data[which(data$Score >= 0.67),]$category = "New Instances"
    data[which(data$Score >= 0.34 & data$Score < 0.67),]$category = "Putative instances"
    data[which(data$Score < 0.34),]$category = "Novel Motif"
    data[which(data$motif_distance == 0 & !is.na(data$Accession)),]$category = "Known Instances"

    links_1 = data.frame(source = rep(all_classes,each = 4), target = rep(c("Highest Similarity","High Similarity", "Low Similarity", "Lowest Similarity"),6), value = 0)
    links_2 = data.frame(source = rep(c("Highest Similarity","High Similarity", "Low Similarity", "Lowest Similarity"),each = 4), target = rep(c("Known Instances", "New Instances", "Novel Motif", "Putative instances"),4), value = 0)

    for (i in 1:nrow(links_1)){
      filtered = dplyr::filter(data, conf == links_1$target[i], ELM_class == links_1$source[i])
      #filtered = filtered[!duplicated(filtered[,c(1:3)]),]
      links_1$value[i] = nrow(filtered)
    }

    for (i in 1:nrow(links_2)){
      filtered = dplyr::filter(data, category == links_2$target[i], conf == links_2$source[i])
      #filtered = filtered[!duplicated(filtered[,c(1:3)]),]
      links_2$value[i] = nrow(filtered)
    }
    links = rbind.data.frame(links_1, links_2)
    links = links[which(links$value != 0),]

    nodes <- data.frame(
      name=c(as.character(links$source),
             as.character(links$target)) %>% unique()
    )

    links$IDsource <- match(links$source, nodes$name)-1
    links$IDtarget <- match(links$target, nodes$name)-1

    # Make the Network
    p <- networkD3::sankeyNetwork(Links = links, Nodes = nodes,
                       Source = "IDsource", Target = "IDtarget",
                       Value = "value", NodeID = "name",
                       sinksRight= F, fontSize = 16)
    return(p)
  }
}

Sankey_H1_chembl = function(H1_clin_chembl, host_viral_ints, human_human_int, oxo_distance = 2){
  H1_clin_chembl = H1_clin_chembl[which(H1_clin_chembl$distance <= oxo_distance),]
  H1_clin_chembl$Clinvar_EFO_name = stringr::str_to_title(H1_clin_chembl$Clinvar_EFO_name)
  record_ids = unique(H1_clin_chembl$record_id)

  mol_names = list()
  for (i in 1:length(record_ids)){
    req = RCurl::getURL(paste0("https://www.ebi.ac.uk/chembl/api/data/compound_record/", record_ids[i], ".json"), async = T)
    req = jsonlite::fromJSON(req)
    mol_names[[i]] = data.frame(record_id = req$record_id, H1_Drug_name = stringr::str_to_title(req$compound_name))
  }
  mol_names = dplyr::bind_rows(mol_names)
  H1_clin_chembl = dplyr::left_join(H1_clin_chembl, mol_names)
  hv_ints = data.frame(V1 = c(host_viral_ints$uniprot_A, host_viral_ints$uniprot_B),
                       V1_sym = c(host_viral_ints$Gene_A, host_viral_ints$Gene_B),
                       V1_org = c(host_viral_ints$org_A, host_viral_ints$org_B))
  hv_ints = hv_ints[which(hv_ints$V1 %in% H1_clin_chembl$V1),]
  hv_ints = hv_ints[!duplicated(hv_ints),]
  hv_ints$V1_sym = paste0(hv_ints$V1_sym, "_", hv_ints$V1_org)
  H1_clin_chembl = dplyr::left_join(H1_clin_chembl, hv_ints[,c(1,2)])

  hh_ints = data.frame(H1 = c(human_human_int$uniprot_A, human_human_int$uniprot_B),
                       H1_sym = c(human_human_int$Gene_A, human_human_int$Gene_B))
  hh_ints = hh_ints[which(hh_ints$H1 %in% H1_clin_chembl$H1),]
  hh_ints = hh_ints[!duplicated(hh_ints),]
  H1_clin_chembl = dplyr::left_join(H1_clin_chembl, hh_ints)

  data = dplyr::select(H1_clin_chembl, H1_Drug_name, H1_sym, V1_sym, Clinvar_EFO_name)
  data = data[!duplicated(data),]

  data$V1_sym[which(data$V1_sym == "M2_A/New York/1682/2009(H1N1)")] = "M2_9INFA/H1N1"

  Drug_H1 = data %>% dplyr::select(H1_Drug_name, H1_sym)
  colnames(Drug_H1) = c("source", "target")
  Drug_H1$lk_grp = "Drug_H1"
  #Drug_H1$value = 1

  Drug_clin = data %>% dplyr::select(Clinvar_EFO_name, H1_Drug_name)
  colnames(Drug_clin) = c("source", "target")
  Drug_clin$lk_grp = "Drug_clin"
  #Drug_clin$value = 1

  H1_v1 = data %>% dplyr::select(H1_sym, V1_sym)
  colnames(H1_v1) = c("source", "target")
  H1_v1$lk_grp = "H1_v1"
  #H1_v1$value = 1

  links = rbind.data.frame(Drug_clin, Drug_H1, H1_v1)
  links = links %>% dplyr::group_by(source, target, lk_grp) %>% dplyr::summarise(value = dplyr::n())
  nodes <- data.frame(
    name=c(as.character(links$source),
           as.character(links$target)) %>% unique()
  )
  for (i in 1:nrow(nodes)){
    if (nodes$name[i] %in% data$H1_Drug_name){
      nodes$nd_grp[i] = "Drug"
    }
    if (nodes$name[i] %in% data$H1_sym){
      nodes$nd_grp[i] = "H1"
    }
    if (nodes$name[i] %in% data$V1_sym){
      nodes$nd_grp[i] = "V1"
    }
    if (nodes$name[i] %in% data$Clinvar_EFO_name){
      nodes$nd_grp[i] = "Disease"
    }
  }
  links$IDsource <- match(links$source, nodes$name)-1
  links$IDtarget <- match(links$target, nodes$name)-1

  colors = RColorBrewer::brewer.pal(name = "Set1", 4)

  my_color <- 'd3.scaleOrdinal() .domain(["Disease", "Drug", "H1", "V1"]) .range(["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"])'

  # Make the Network
  p <- networkD3::sankeyNetwork(Links = links, Nodes = nodes,
                                Source = "IDsource", Target = "IDtarget", NodeID = "name",Value = "value",
                                NodeGroup = "nd_grp",colourScale = my_color,
                                sinksRight= F, fontSize = 16, margin = list("left" = -50, "right" = -50))
  return(p)
}
