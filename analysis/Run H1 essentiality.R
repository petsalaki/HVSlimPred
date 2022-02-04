library(HVSlimPred)

human_human_int = read.csv("data/input/human_human_intact.csv.gz")
human_prots = unique(c(human_human_int$uniprot_A, human_human_int$uniprot_B))
x = human_human_int %>% dplyr::select(uniprot_A, Gene_A)
colnames(x) = c("uniprot", "symbol")
y = human_human_int %>% dplyr::select(uniprot_B, Gene_B)
colnames(y) = c("uniprot", "symbol")
human_sym_annot = rbind.data.frame(x, y)
human_sym_annot = human_sym_annot[!duplicated(human_sym_annot),]
human_sym_annot = human_sym_annot[!duplicated(human_sym_annot$uniprot),]

tmp = tempfile()
curl::curl_download(url = "https://www.ebi.ac.uk/biostudies/files/S-BSST479/EV_Table_1.xlsx", destfile = tmp)

Ev1_table = readxl::read_xlsx(tmp, sheet = 1, skip = 1, na = c("", NA))
colnames(Ev1_table)[1] = "Gene"
new_hits_ints = readRDS("data/output/HV_final_hits.rds")
H1_prots = unique(sub(".*_", "", new_hits_ints$Dataset))
H2_prots = unique(new_hits_ints$uniprot)
H1_mapping = human_sym_annot[which(human_sym_annot$uniprot %in% H1_prots),]
H2_mapping = human_sym_annot[which(human_sym_annot$uniprot %in% H2_prots),]

H1_missed = H1_mapping[which(H1_mapping$symbol %nin% Ev1_table$Gene),]
H2_missed = H2_mapping[which(H2_mapping$symbol %nin% Ev1_table$Gene),]

H2_prots = H2_prots[which(H2_prots %in% H2_mapping$uniprot)]

Get_centools_classif = function(upID, mapping_df){
  gene_sym = mapping_df$symbol[which(mapping_df$uniprot == upID)]
  if (gene_sym %in% Ev1_table$Gene){
    filtered = Ev1_table[which(Ev1_table$Gene == gene_sym),]
    if (is.na(filtered$Integrated_Cluster) & filtered$BAGEL_Training_set_INTEGRATED == "Not in training"){
      return(data.frame(H1 = upID, CEN_Clust = NA))
    }
    else{
      if (!is.na(filtered$Integrated_Cluster)){
        result = data.frame(H1 = upID, CEN_Clust = as.character(filtered$Integrated_Cluster))
        return(result)
      }
      else{
        result = data.frame(H1 = upID, CEN_Clust = as.character(filtered$BAGEL_Training_set_INTEGRATED))
        return(result)
      }
    }
  }
  else{
    return(NULL)
  }
}

H1_CEN = list()
for (i in 1:length(H1_prots)){
  H1_CEN[[i]] = Get_centools_classif(H1_prots[i], mapping_df = H1_mapping)
}
H1_CEN = dplyr::bind_rows(H1_CEN)
H1_CEN$CEN_Clust = stringr::str_replace_all(H1_CEN$CEN_Clust, c("Essentials" = "Essential", "Non Essential" = "Non_essential", "Rare_Context" = "Context"))
colnames(H1_CEN)[1] = "uniprot"
H2_CEN = list()
for (i in 1:length(H2_prots)){
  H2_CEN[[i]] = Get_centools_classif(H2_prots[i], mapping_df = H2_mapping)
}
H2_CEN = dplyr::bind_rows(H2_CEN)
H2_CEN$CEN_Clust = stringr::str_replace_all(H2_CEN$CEN_Clust, c("Essentials" = "Essential", "Non Essential" = "Non_essential", "Rare_Context" = "Context"))
colnames(H2_CEN)[1] = "uniprot"

Enrich_essentiality = function(new_hits_ints,filter_sym, bg_CEN, bg_type = c("H1", "H2")){
  df = all_data[[filter_sym]]
  main = plyr::match_df(new_hits_ints, df, on = c("uniprot", "Start_Pos", "End_Pos"))
  if (bg_type == "H1"){
    H1_prots = unique(sub(".*_", "", main$Dataset))
    H1_prots = H1_prots[which(H1_prots %in% bg_CEN$uniprot)]

    if (stringr::str_detect(filter_sym, "C") & stringr::str_detect(filter_sym, "D")){
      H1_prots = H1_prots[which(H1_prots %in% intersect(all_data$C$domain_protein, all_data$D$domain_protein))]
    }
    else if (stringr::str_detect(filter_sym, "C") & stringr::str_detect(filter_sym, "D", negate = T)){
      H1_prots = H1_prots[which(H1_prots %in% all_data$C$domain_protein)]
    }
    else if (stringr::str_detect(filter_sym, "C", negate = T) & stringr::str_detect(filter_sym, "D")){
      H1_prots = H1_prots[which(H1_prots %in% all_data$D$domain_protein)]
    }

    if (length(H1_prots) == 0){
      return(NULL)
    }

    final = list()
    for (i in c("Context", "Essential", "Non_essential")){
      CEN = bg_CEN[which(bg_CEN$CEN_Clust == i),]

      TP = length(intersect(H1_prots, CEN$uniprot))
      FN = length(unique(CEN$uniprot[which(CEN$uniprot %nin% H1_prots)]))
      FP = length(unique(H1_prots[which(H1_prots %nin% CEN$uniprot)]))
      TN = length(unique(bg_CEN$uniprot[which(bg_CEN$uniprot %nin% c(CEN$uniprot, H1_prots))]))

      numbs = c(TP, FN, FP, TN)
      fisher_results = fisher_and_F1(numbers = numbs, pval_method = "two.sided")
      fisher_results$CEN_clust = i
      fisher_results = fisher_results[,c(ncol(fisher_results), 1:(ncol(fisher_results) - 1))]
      final[[i]] = fisher_results
    }
    final = dplyr::bind_rows(final)
    return(final)
  }
  if (bg_type == "H2"){
    H2_prots = unique(main$uniprot)
    H2_prots = H2_prots[which(H2_prots %in% bg_CEN$uniprot)]
    if (length(H2_prots) == 0){
      return(NULL)
    }

    final = list()
    for (i in c("Context", "Essential", "Non_essential")){
      CEN = bg_CEN[which(bg_CEN$CEN_Clust == i),]

      TP = length(intersect(H2_prots, CEN$uniprot))
      FN = length(unique(CEN$uniprot[which(CEN$uniprot %nin% H2_prots)]))
      FP = length(unique(H2_prots[which(H2_prots %nin% CEN$uniprot)]))
      TN = length(unique(bg_CEN$uniprot[which(bg_CEN$uniprot %nin% c(CEN$uniprot, H2_prots))]))

      numbs = c(TP, FN, FP, TN)
      fisher_results = fisher_and_F1(numbers = numbs, pval_method = "greater")
      fisher_results$CEN_clust = i
      fisher_results = fisher_results[,c(ncol(fisher_results), 1:(ncol(fisher_results) - 1))]
      final[[i]] = fisher_results
    }
    final = dplyr::bind_rows(final)
    return(final)
  }

}

all_data = readRDS("data/output/HV_all_filters_data.rds")
conserved_all_data = all_data$conserved_only
all_data = all_data$all_hits


All_filters_essentiality = list()
for (i in 1:nrow(all_data$Metadata)){
  x = Enrich_essentiality(new_hits_ints = new_hits_ints, filter_sym = all_data$Metadata$sym[i], bg_CEN = H1_CEN, bg_type = "H1")
  y = Enrich_essentiality(new_hits_ints = new_hits_ints, filter_sym = all_data$Metadata$sym[i], bg_CEN = H2_CEN, bg_type = "H2")

  All_filters_essentiality[[all_data$Metadata$Term[i]]] = list("H1_Essentiality" = x, "H2_Essentiality" = y)
}
saveRDS(All_filters_essentiality, "data/output/plot_objects/All_filters_essentiality.rds")

all_data = conserved_all_data
All_filters_essentiality_conserved = list()
for (i in 1:nrow(all_data$Metadata)){
  x = Enrich_essentiality(new_hits_ints = new_hits_ints, filter_sym = all_data$Metadata$sym[i], bg_CEN = H1_CEN, bg_type = "H1")
  y = Enrich_essentiality(new_hits_ints = new_hits_ints, filter_sym = all_data$Metadata$sym[i], bg_CEN = H2_CEN, bg_type = "H2")

  All_filters_essentiality_conserved[[all_data$Metadata$Term[i]]] = list("H1_Essentiality" = x, "H2_Essentiality" = y)
}
saveRDS(All_filters_essentiality_conserved, "data/output/plot_objects/All_filters_essentiality_conserved.rds")
