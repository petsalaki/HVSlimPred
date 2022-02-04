library(HVSlimPred)


# Prepare data for SLiMEnrich ---------------------------------------------
dir.create("data/input/SLiMEnrich")
dir.create("data/input/SLiMEnrich/Host_viral")
dir.create("data/input/SLiMEnrich/Human_only")
dir.create("data/input/SLiMEnrich/Host_viral/Default")
dir.create("data/input/SLiMEnrich/Human_only/Default")
dir.create("data/input/SLiMEnrich/Host_viral/Custom")
dir.create("data/input/SLiMEnrich/Human_only/Custom")

HV_final_hits = readRDS("data/output/HV_final_hits.rds")
HV_Interaction_file = HV_final_hits %>% tidyr::separate(Dataset, into = c("V1", "H1"), sep = "_") %>% dplyr::select(uniprot, H1)
HV_Interaction_file = HV_Interaction_file[!duplicated(HV_Interaction_file),]
colnames(HV_Interaction_file) = c("mProtein", "dProtein")
write.csv(HV_Interaction_file,"data/input/SLiMEnrich/Host_viral/Default/Interaction_file.csv", row.names = F)



HV_Motif_file = HV_final_hits[which(HV_final_hits$new_score > 0.5),]
HV_Interaction_file = HV_Motif_file %>% tidyr::separate(Dataset, into = c("V1", "H1"), sep = "_") %>% dplyr::select(uniprot, H1)
HV_Interaction_file = HV_Interaction_file[!duplicated(HV_Interaction_file),]
HV_Motif_file = HV_Motif_file %>% dplyr::select(Dataset, Pattern, uniprot, Org, Start_Pos,
                                                End_Pos, Match, Desc, Id)
HV_Motif_file = HV_Motif_file[!duplicated(HV_Motif_file),]

colnames(HV_Interaction_file) = c("mProtein", "dProtein")
write.csv(HV_Interaction_file,"data/input/SLiMEnrich/Host_viral/Custom/Interaction_file.csv", row.names = F)
write.csv(HV_Motif_file,"data/input/SLiMEnrich/Host_viral/Custom/Motif_file.csv", row.names = F)


human_final_hits = readRDS("data/output/human_final_hits.rds")
human_Interaction_file = human_final_hits %>% dplyr::select(uniprot, Dataset)
human_Interaction_file = human_Interaction_file[!duplicated(human_Interaction_file),]
colnames(human_Interaction_file) = c("mProtein", "dProtein")
write.csv(human_Interaction_file,"data/input/SLiMEnrich/Human_only/Default/Interaction_file.csv", row.names = F)


human_Motif_file = human_final_hits[which(human_final_hits$new_score > 0.5),]
human_Interaction_file = human_Motif_file %>% dplyr::select(uniprot, Dataset)
human_Interaction_file = human_Interaction_file[!duplicated(human_Interaction_file),]
human_Motif_file = human_Motif_file %>% dplyr::select(Dataset, Pattern, uniprot, Org, Start_Pos,
                               End_Pos, Match, Desc, Id)
human_Motif_file = human_Motif_file[!duplicated(human_Motif_file),]
colnames(human_Interaction_file) = c("mProtein", "dProtein")
write.csv(human_Interaction_file,"data/input/SLiMEnrich/Human_only/Custom/Interaction_file.csv", row.names = F)
write.csv(human_Motif_file,"data/input/SLiMEnrich/Human_only/Custom/Motif_file.csv", row.names = F)

# GOBP enrichment for H1 and motif_proteins -------------------------------
require(clusterProfiler)
require(org.Hs.eg.db)
require(enrichplot)

all_data = readRDS("data/output/HV_all_filters_data.rds")

HV_final_hits = readRDS("data/output/HV_final_hits.rds")

H1_prots = HV_final_hits %>% tidyr::separate(Dataset, into = c("V1", "H1"), "_")
H1_prots = unique(H1_prots$H1)

motif_prots = unique(HV_final_hits$uniprot)

H1_Interactome_entrez = AnnotationDbi::select(org.Hs.eg.db, keys = H1_prots, keytype = "UNIPROT",columns = "ENTREZID")
H1_Interactome_entrez = H1_Interactome_entrez[!is.na(H1_Interactome_entrez$ENTREZID),]

motif_Interactome_entrez = AnnotationDbi::select(org.Hs.eg.db, keys = motif_prots, keytype = "UNIPROT",columns = "ENTREZID")
motif_Interactome_entrez = motif_Interactome_entrez[!is.na(motif_Interactome_entrez$ENTREZID),]

prepare_GOBP_query = function(entrez_ids,req_data,filter_sym, prot_type = c("H1","motif")){
  df = req_data[[filter_sym]]

  matched = plyr::match_df(HV_final_hits, df, on = c("uniprot", "Start_Pos", "End_Pos"))
  if (nrow(matched) == 0){
    return(NULL)
  }
  if (prot_type == "H1"){
    matched = matched %>% tidyr::separate(Dataset, into = c("V1", "H1"), "_")
    q_prots = unique(matched$H1)
    if (stringr::str_detect(filter_sym, "C") & stringr::str_detect(filter_sym, "D")){
      q_prots = q_prots[which(q_prots %in% intersect(req_data$C$domain_protein, req_data$D$domain_protein))]
    }
    else if (stringr::str_detect(filter_sym, "C") & stringr::str_detect(filter_sym, "D", negate = T)){
      q_prots = q_prots[which(q_prots %in% req_data$C$domain_protein)]
    }
    else if (stringr::str_detect(filter_sym, "C", negate = T) & stringr::str_detect(filter_sym, "D")){
      q_prots = q_prots[which(q_prots %in% req_data$D$domain_protein)]
    }
  }
  else{
    q_prots = unique(matched$uniprot)
  }
  q_entrez = entrez_ids$ENTREZID[which(entrez_ids$UNIPROT %in% q_prots)]
  return(q_entrez)
}
enrich_GOBP = function(req_data, prot_type = c("H1","motif"), remove_iELM = F){
  # if (cons){
  #   req_data = req_data[["conserved_only"]]
  # }
  # else{
  #   req_data = req_data[["all_hits"]]
  # }
  req_cons = req_data[["conserved_only"]]
  req_all = req_data[["all_hits"]]

  if (prot_type == "H1"){
    entrez_universe = H1_Interactome_entrez
  }
  else{
    entrez_universe = motif_Interactome_entrez
  }

  gcSample = list()
  for (i in 1:nrow(req_all$Metadata)){
    res_all = prepare_GOBP_query(entrez_ids = entrez_universe, req_data = req_all,
                             filter_sym = req_all$Metadata$sym[i], prot_type = prot_type)
    res_cons = prepare_GOBP_query(entrez_ids = entrez_universe, req_data = req_cons,
                                  filter_sym = req_cons$Metadata$sym[i], prot_type = prot_type)
    if (is.null(res_all) | length(res_all) == 0){
      gcSample[[req_all$Metadata$Term[i]]] = NULL
    }
    else{
      gcSample[[req_all$Metadata$Term[i]]] = res_all
    }
    if (is.null(res_cons) | length(res_cons) == 0){
      gcSample[[paste0("Cons ",req_all$Metadata$Term[i])]] = NULL
    }
    else{
      gcSample[[paste0("Cons ",req_all$Metadata$Term[i])]] = res_cons
    }

  }
  #replace_filter_names = c("Enriched_Domain" = "E", "PepSite" = "P", "ClinVar_Path" = "C", "iELM_HMMs" = "I", "&" = "")
  replace_filter_names = c("Enriched_Domain" = "Enriched Domain", "PepSite" = "PepSite", "ClinVar_Path" = "ClinVar", "iELM_HMMs" = "iELM", "&" = ",")
  names(gcSample) = stringr::str_replace_all(names(gcSample), replace_filter_names)

  if (remove_iELM){
    gcSample = gcSample[stringr::str_detect(names(gcSample), "iELM", negate = T)]
  }
  ck <- compareCluster(geneCluster = gcSample, fun = "enrichGO", OrgDb = org.Hs.eg.db, ont = "BP", universe = entrez_universe$ENTREZID)
  return(ck)
}

pb = txtProgressBar(min = 0, max = 2, style = 3)
counter = 0
for (i in c("H1","motif")){
  counter = counter + 1
  final = enrich_GOBP(req_data = all_data, prot_type = i, remove_iELM = T)
  saveRDS(final, paste0("data/output/plot_objects/", "GOBP_enrichment_", i,".rds"))
  setTxtProgressBar(pb, counter)
}

