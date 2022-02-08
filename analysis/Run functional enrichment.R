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
require(GOSemSim)

all_data = readRDS("data/output/HV_all_filters_data.rds")

HV_final_hits = readRDS("data/output/HV_final_hits.rds")

H1_prots = HV_final_hits %>% tidyr::separate(Dataset, into = c("V1", "H1"), "_")
H1_prots = unique(H1_prots$H1)
motif_prots = unique(HV_final_hits$uniprot)

all_prots = unique(c(H1_prots, motif_prots))

prepare_GOBP_query = function(entrez_ids,req_data,filter_sym, prot_type = c("H1","motif", "Both")){
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
  else if (prot_type == "motif"){
    q_prots = unique(matched$uniprot)
  }
  else{
    H1_prots = matched %>% tidyr::separate(Dataset, into = c("V1", "H1"), "_")
    H1_prots = unique(H1_prots$H1)

    q_prots = unique(c(H1_prots, matched$uniprot))
  }
  q_entrez = entrez_ids$ENTREZID[which(entrez_ids$UNIPROT %in% q_prots)]
  return(q_entrez)
}
enrich_GOBP = function(req_data, prot_type = c("H1","motif", "Both"), remove_iELM = F, bg_slim_input = F){
  # if (cons){
  #   req_data = req_data[["conserved_only"]]
  # }
  # else{
  #   req_data = req_data[["all_hits"]]
  # }
  req_cons = req_data[["conserved_only"]]
  req_all = req_data[["all_hits"]]

  if (prot_type == "H1"){
    entrez_universe = AnnotationDbi::select(org.Hs.eg.db, keys = H1_prots, keytype = "UNIPROT",columns = "ENTREZID")
    entrez_universe = entrez_universe[!is.na(entrez_universe$ENTREZID),]
  }
  else if (prot_type == "motif"){
    entrez_universe = AnnotationDbi::select(org.Hs.eg.db, keys = motif_prots, keytype = "UNIPROT",columns = "ENTREZID")
    entrez_universe = entrez_universe[!is.na(entrez_universe$ENTREZID),]
  }
  else{
    if(bg_slim_input){
      entrez_universe = AnnotationDbi::select(org.Hs.eg.db, keys = unique(unlist(req_all[["all_input_prots"]])), keytype = "UNIPROT",columns = "ENTREZID")
      entrez_universe = entrez_universe[!is.na(entrez_universe$ENTREZID),]
    }
    else{
      entrez_universe = AnnotationDbi::select(org.Hs.eg.db, keys = all_prots, keytype = "UNIPROT",columns = "ENTREZID")
      entrez_universe = entrez_universe[!is.na(entrez_universe$ENTREZID),]
    }
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

final = enrich_GOBP(req_data = all_data, prot_type = "Both", remove_iELM = T, bg_slim_input = T)
saveRDS(final, paste0("data/output/plot_objects/", "GOBP_enrichment_all",".rds"))


closest_to_center <- function(x, centers) {
  # compute squared euclidean distance from each sample to each cluster center
  tmp <- sapply(seq_len(nrow(x)),
                function(i) apply(centers, 1,
                                  function(v) sum((x[i, ]-v)^2)))
  if (1 %in% dim(as.matrix(tmp))){
    names(tmp) = rownames(x)
    final = list()
    final[[1]] = names(tmp)[which(tmp == min(tmp))]
  }
  else{
    colnames(tmp) = rownames(x)
    tmp = t(tmp)
    final = list()
    for (j in 1:ncol(tmp)){
      final[[colnames(tmp)[j]]] = rownames(tmp)[which(tmp[,j] == min(tmp[,j]))]
    }
  }
  return(final)
}
dom_mat_kmeans = function(dom_mat, kmax = 10, bootstrap = 100, GoSim_thresh = 0.5, defined_K = F){
  counter = 0
  final_clusts = list()
  if (defined_K){
    km_clust = factoextra::eclust(dom_mat, "kmeans", k = kmax, nstart = 25, nboot = bootstrap)
  }
  else{
    for (l in c(10:2)){
      km_clust = try(factoextra::eclust(dom_mat, "kmeans", k.max = l, nstart = 25, nboot = bootstrap))
      if("try-error" %in% class(km_clust)){next()}
      else{break()}
    }
  }

  center_close = closest_to_center(km_clust$data, km_clust$centers)

  for (i in sort(unique(km_clust$cluster))){
    if (km_clust$size[i] == 1){
      counter = counter + 1
      final_clusts[[paste0("C_",counter)]][["members"]] = names(km_clust$cluster[which(km_clust$cluster == i)])
      final_clusts[[paste0("C_",counter)]][["min_dist_center"]] = names(km_clust$cluster[which(km_clust$cluster == i)])
      next()
    }
    x = km_clust$data
    x = x[,names(km_clust$cluster[which(km_clust$cluster == i)])]
    x = x[names(km_clust$cluster[which(km_clust$cluster == i)]),]

    if (min(x) < GoSim_thresh){
      if (nrow(x) == 2){
        for (k in 1:2){
          counter = counter + 1
          final_clusts[[paste0("C_",counter)]][["members"]] = rownames(x)[k]
          final_clusts[[paste0("C_",counter)]][["min_dist_center"]] = rownames(x)[k]
        }
        next()
      }
      res = dom_mat_kmeans(dom_mat = x, kmax = 2, defined_K = T)
      for (k in 1:length(res)){
        counter = counter + 1
        final_clusts[[paste0("C_", counter)]] = res[[k]]
      }
      # if (nrow(x) <= 10){
      # }
      # else{
      #   res = dom_mat_kmeans(dom_mat = x, defined_K = F)
      #   for (j in 1:length(res)){
      #     counter = counter + 1
      #     final_clusts[[paste0("C_", counter)]] = res[[j]]
      #   }
      # }
    }
    else{
      counter = counter + 1
      final_clusts[[paste0("C_",counter)]][["members"]] = names(km_clust$cluster[which(km_clust$cluster == i)])
      final_clusts[[paste0("C_",counter)]][["min_dist_center"]] = center_close[[as.character(i)]]
    }
  }
  return(final_clusts)

}
hsGO <- godata('org.Hs.eg.db', ont="BP", computeIC = F)

GO_ID = unique(final@compareClusterResult$ID)
GO_ID = GO_ID[which(GO_ID %in% hsGO@geneAnno$GO)]
GO_ID_pairs = expand.grid(GO_ID, GO_ID)
colnames(GO_ID_pairs) = c("GO1", "GO2")
GO_ID_pairs$GO1 = as.character(GO_ID_pairs$GO1)
GO_ID_pairs$GO2 = as.character(GO_ID_pairs$GO2)
pb = txtProgressBar(min = 0, max = nrow(GO_ID_pairs), style = 3)
GO_ID_pairs$SemSim = 0
for (i in 1:nrow(GO_ID_pairs)){
  GO_ID_pairs$SemSim[i] = goSim(GOID1 = GO_ID_pairs$GO1[i], GOID2 = GO_ID_pairs$GO2[i], semData=hsGO, measure="Wang")
  setTxtProgressBar(pb, i)
}

GO_ID_mat = tidyr::pivot_wider(GO_ID_pairs, names_from = GO1, values_from = SemSim)
rNames = GO_ID_mat$GO2
GO_ID_mat = GO_ID_mat[,-1]
rownames(GO_ID_mat) = rNames
GO_ID_mat = as.matrix(GO_ID_mat)
saveRDS(GO_ID_mat, "data/output/plot_objects/GOBP_Semsim_matrix.rds")

GO_id_clusts = dom_mat_kmeans(dom_mat = GO_ID_mat, GoSim_thresh = 0.5, bootstrap = 2)
all_min_dist_GO = list()
for (i in 1:length(GO_id_clusts)){
  all_min_dist_GO[[i]] = data.frame(Clust = names(GO_id_clusts)[i],
                                    GO_ID = GO_id_clusts[[i]]$min_dist_center)
}
all_min_dist_GO = dplyr::bind_rows(all_min_dist_GO)

GO_all = ontologyIndex::get_OBO("data/input/dcGOR/go-basic.obo")
GO_annot = data.frame(GO_all$name)
GO_annot$GO_ID = rownames(GO_annot)
colnames(GO_annot) = c("GO_name", "GO_ID")
all_min_dist_GO = dplyr::left_join(all_min_dist_GO, GO_annot)

final@compareClusterResult = final@compareClusterResult[which(final@compareClusterResult$ID %in% all_min_dist_GO$GO_ID),]
saveRDS(final, paste0("data/output/plot_objects/", "GOBP_enrichment_all",".rds"))

saveRDS(list("GOBP_clusters" = GO_id_clusts, "clust_GO_annotation" = all_min_dist_GO), "data/output/GOSemSim_GOBP_prot_clusts.rds")


