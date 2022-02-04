require(HVSlimPred)
require(igraph)
require(dplyr)
require(RCy3)
require(pheatmap)
require(tidyverse)
require(enrichplot)
require(GOSemSim)
require(org.Hs.eg.db)
require(factoextra)
require(ComplexUpset)
require(UpSetR)
require(hrbrthemes)
require(viridis)
# Figure 2 and S1 ---------------------------------------------------------
HV_all_eval = readRDS("data/output/HV_all_filters_eval_OR.rds")

All_levels_bar_plot = function(req_data, prot = F, database, prot_comparison = T, letter_labels = T, remove_iELM = T){
  if (letter_labels){
    replace_filter_names = c("Enriched_Domain" = "E",
                             "PepSite" = "P", "ClinVar_Path" = "C", "iELM_HMMs" = "I",
                             "&" = "", "qslim_predicted" = "raw")
  }
  else{
    replace_filter_names = c("Enriched_Domain" = "Enriched Domain",
                             "PepSite" = "PepSite", "ClinVar_Path" = "ClinVar", "iELM_HMMs" = "iELM",
                             "&" = ",", "qslim_predicted" = "QSLiMFinder")
  }
  if (!prot){
    All_levels_result = list()
    filtered = req_data[["motif"]][[database]]
    filter_names = names(filtered[["all_hits"]])
    All_F1_scores = list()
    for (i in 1:length(filter_names)){
      if (filter_names[i] %in% names(filtered[["conserved_only"]])){
        x = filtered[["conserved_only"]][[filter_names[i]]]$All_F1_scores[which(filtered[["conserved_only"]][[filter_names[i]]]$All_F1_score$site_resid == "site_wise"),]
        x$Conserved = "Yes"
      }
      else{
        x = NULL
      }
      y = filtered[["all_hits"]][[filter_names[i]]]$All_F1_scores[which(filtered[["all_hits"]][[filter_names[i]]]$All_F1_score$site_resid == "site_wise"),]
      y$Conserved = "No"
      All_F1_scores[[i]] = rbind.data.frame(x,y)
    }
    All_F1_scores = dplyr::bind_rows(All_F1_scores)
    All_F1_scores$weight_site_resid = paste0(All_F1_scores$weight_type, "_", All_F1_scores$site_resid)
    All_F1_scores$F0.5 = (1.25 * All_F1_scores$Rc * All_F1_scores$Pr) / (0.25 * (All_F1_scores$Rc + All_F1_scores$Pr))
    All_F1_scores = All_F1_scores[!is.nan(All_F1_scores$F0.5),]
    All_F1_scores = All_F1_scores[!is.na(All_F1_scores$F0.5),]
    All_F1_scores$Reference = database

    All_levels_result = All_F1_scores
    All_levels_result = All_levels_result[which(All_levels_result$weight_site_resid == "residue_site_wise"),]

    #All_levels_result$title = stringr::str_replace_all(All_levels_result$title, "&", ",")
    #All_levels_result$title = stringr::str_replace_all(All_levels_result$title, "_", " ")
    All_levels_result$title = stringr::str_replace_all(All_levels_result$title, replace_filter_names)

    colnames(All_levels_result)[which(names(All_levels_result) == "odds_ratio")] = "Odds ratio"

    All_levels_result %>%
      dplyr::group_by(Conserved) %>%
      dplyr::arrange(Conserved, desc(F0.5), .by_group = T) %>%
      tidyr::unite("Conserved_filter", Conserved, title, sep = "_", remove = FALSE) %>%
      data.frame() %>%
      dplyr::mutate(Conserved_filter = factor(Conserved_filter, levels = Conserved_filter))

    if(remove_iELM){
      All_levels_result = All_levels_result[stringr::str_detect(All_levels_result$title, "I", negate = T),]
    }


    xc_abline = All_levels_result$F0.5[which(All_levels_result$title == "QSLiMFinder" & All_levels_result$Conserved == "No")]
    xc_abline_cons = All_levels_result$F0.5[which(All_levels_result$title == "QSLiMFinder" & All_levels_result$Conserved == "Yes")]

    p = ggplot2::ggplot(All_levels_result, ggplot2::aes(y=F0.5, x= reorder(title, -F0.5, max), fill = Conserved, label = NULL)) +
      ggplot2::geom_bar(stat="identity", alpha=0.7, width=0.8, position = ggplot2::position_dodge()) +
      ggplot2::theme_bw() +
      ggplot2::coord_flip() +
      ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 14), axis.text.x = ggplot2::element_text(size = 12)) +
      ggplot2::geom_hline(yintercept = xc_abline, color = "Red", linetype = "dashed") +
      ggplot2::geom_hline(yintercept = xc_abline_cons, color = "blue", linetype = "dashed") +
      ggplot2::ggtitle(paste0("F0.5", " ", database)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      # ggplot2::facet_wrap(~Conserved)

    return(p)
  }
  else{
    if (prot_comparison){
      cons = req_data[["prot"]][[database]][["conserved_only"]][["comparison_prot"]]
      cons$Conserved = "Yes"
      not_cons = req_data[["prot"]][[database]][["all_hits"]][["comparison_prot"]]
      not_cons$Conserved = "No"

      All_levels_result = rbind.data.frame(not_cons, cons)
      All_levels_result = All_levels_result[which(All_levels_result$odds_ratio != 0),]

      #All_levels_result$title = stringr::str_replace_all(All_levels_result$title, "&", ",")
      #All_levels_result$title = stringr::str_replace_all(All_levels_result$title, "_", " ")
      All_levels_result$title = stringr::str_replace_all(All_levels_result$title, replace_filter_names)

      colnames(All_levels_result)[which(names(All_levels_result) == "odds_ratio")] = "Odds ratio"

      if(remove_iELM){
        All_levels_result = All_levels_result[stringr::str_detect(All_levels_result$title, "I", negate = T),]
      }

      p = ggplot2::ggplot(All_levels_result, ggplot2::aes(y=log2(`Odds ratio`), x= reorder(title, -`Odds ratio`, max), fill = Conserved, color = "blue", label = NULL)) +
        ggplot2::geom_bar(stat="identity", alpha=0.7, width=0.8, position = ggplot2::position_dodge(), colour = "black") +
        ggplot2::theme_bw() +
        ggplot2::coord_flip() +
        ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 14), axis.text.x = ggplot2::element_text(size = 12)) +
        #ggplot2::geom_hline(yintercept = xc_abline, color = "Red", linetype = "dashed") +
        ggplot2::ggtitle(paste0(database)) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::ylab("log2(Odds ratio)")
      return(p)
    }
    else{
      cons = req_data[["prot"]][[database]][["conserved_only"]][["no_comparison_prot"]]
      cons$Conserved = "Yes"
      not_cons = req_data[["prot"]][[database]][["all_hits"]][["no_comparison_prot"]]
      not_cons$Conserved = "No"

      All_levels_result = rbind.data.frame(not_cons, cons)
      All_levels_result = All_levels_result[which(All_levels_result$odds_ratio != 0),]

      #All_levels_result$title = stringr::str_replace_all(All_levels_result$title, "&", ",")
      #All_levels_result$title = stringr::str_replace_all(All_levels_result$title, "_", " ")
      All_levels_result$title = stringr::str_replace_all(All_levels_result$title, replace_filter_names)

      colnames(All_levels_result)[which(names(All_levels_result) == "odds_ratio")] = "Odds ratio"

      if(remove_iELM){
        All_levels_result = All_levels_result[stringr::str_detect(All_levels_result$title, "I", negate = T),]
      }

      p = ggplot2::ggplot(All_levels_result, ggplot2::aes(y=log2(`Odds ratio`), x= reorder(title, -`Odds ratio`, max), fill = Conserved, color = "blue", label = NULL)) +
        ggplot2::geom_bar(stat="identity", alpha=0.7, width=0.8, position = ggplot2::position_dodge(), colour = "black") +
        ggplot2::theme_bw() +
        ggplot2::coord_flip() +
        ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 14), axis.text.x = ggplot2::element_text(size = 12)) +
        #ggplot2::geom_hline(yintercept = xc_abline, color = "Red", linetype = "dashed") +
        ggplot2::ggtitle(paste0(database)) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::ylab("log2(Odds ratio)")
      return(p)
    }
  }
}
upset_plot_motif = function(req_data, prot_only = F, database, prot_comparison = T, remove_iELM = T){
  if (!prot_only){
    upset_all_data = list()
    for (j in c("conserved_only", "all_hits")){
      req = req_data[["motif"]][[database]][[j]]
      upsetr_expression = c()
      for (i in 1:length(req)){
        total = c()
        metrics = req[[i]]$all_metrics_per_motif
        intersection = metrics$Id[which(metrics$Rc_res != 0)]
        upsetr_expression = c(upsetr_expression, length(intersection))
      }
      names(upsetr_expression) = names(req)
      names(upsetr_expression) = stringr::str_replace_all(names(upsetr_expression), c("Enriched_Domain" = "Enriched Domain",
                                                                                      "PepSite" = "PepSite", "ClinVar_Path" = "ClinVar", "iELM_HMMs" = "iELM",
                                                                                      "qslim_predicted" = "QSLiMFinder"))
      if (remove_iELM){
        upsetr_expression = upsetr_expression[stringr::str_detect(names(upsetr_expression),"iELM", negate = T)]
      }
      upset_all_data[[j]] = upsetr_expression
    }
    common_filters = dplyr::intersect(names(upset_all_data$conserved_only), names(upset_all_data$all_hits))
    for (i in common_filters){
      upset_all_data$all_hits[i] = upset_all_data$all_hits[i] - upset_all_data$conserved_only[i]
    }
    for (j in 1:length(upset_all_data)){
      upset_all_data[[j]] = fromExpression(upset_all_data[[j]])
      upset_all_data[[j]]$Conserved = ifelse(names(upset_all_data)[j] == "conserved_only", "Yes", "No")
    }
    upset_all_data = dplyr::bind_rows(upset_all_data)
    intersect_cols = colnames(upset_all_data)[which(colnames(upset_all_data) != "Conserved")]

    # p1 = ComplexUpset::upset(
    #   upset_all_data,
    #   intersect_cols,
    #   base_annotations = list(
    #     'Number of known motifs'= ComplexUpset::intersection_size()
    #   ),
    #   annotations = list(
    #     'Conserved Known motifs'=(
    #       ggplot2::ggplot(mapping=aes(fill=Conserved))
    #       + ggplot2::geom_bar(stat='count', position='fill')
    #       + ggplot2::scale_y_continuous(labels=scales::percent_format())
    #       + ggplot2::scale_fill_manual(values=c(
    #         'Yes'='#E41A1C', 'No'='#377EB8'))
    #       + ggplot2::ylab('Conserved Known motifs')
    #     )
    #   ),
    #   width_ratio=0.1
    # )

    p1 = ComplexUpset::upset(
      upset_all_data,
      intersect_cols,
      base_annotations = list(
        'Number of known motifs'= ComplexUpset::intersection_size(
          counts = T,
          mapping = aes(fill=Conserved)
        ) + ggplot2::scale_fill_manual(values=c(
          'Yes'='#E41A1C', 'No'='#377EB8'))
      ),
      width_ratio=0.1
    )


    #p1 = upset(x, order.by = "freq", point.size = 3,text.scale = c(1.25,1.25,1.25,1.25,1.25,1.25), mainbar.y.label = "Number of Known Motifs", line.size = 1)
    return(p1)
  }
  else{
    if (prot_comparison){
      upset_all_data = list()
      for (j in c("conserved_only", "all_hits")){
        req = req_data[["prot"]][[database]][[j]][["comparison_prot"]]
        upsetr_expression = req$TP
        names(upsetr_expression) = req$title
        names(upsetr_expression) = stringr::str_replace_all(names(upsetr_expression), c("Enriched_Domain" = "Enriched Domain",
                                                                                        "PepSite" = "PepSite", "ClinVar_Path" = "ClinVar", "iELM_HMMs" = "iELM",
                                                                                        "qslim_predicted" = "QSLiMFinder"))
        if (remove_iELM){
          upsetr_expression = upsetr_expression[stringr::str_detect(names(upsetr_expression),"iELM", negate = T)]
        }
        upset_all_data[[j]] = upsetr_expression
      }
      common_filters = dplyr::intersect(names(upset_all_data$conserved_only), names(upset_all_data$all_hits))
      for (i in common_filters){
        upset_all_data$all_hits[i] = upset_all_data$all_hits[i] - upset_all_data$conserved_only[i]
      }
      for (j in 1:length(upset_all_data)){
        upset_all_data[[j]] = fromExpression(upset_all_data[[j]])
        upset_all_data[[j]]$Conserved = ifelse(names(upset_all_data)[j] == "conserved_only", "Yes", "No")
      }
      upset_all_data = dplyr::bind_rows(upset_all_data)
      intersect_cols = colnames(upset_all_data)[which(colnames(upset_all_data) != "Conserved")]

      p1 = ComplexUpset::upset(
        upset_all_data,
        intersect_cols,
        base_annotations = list(
          'Number of known motif proteinss'= ComplexUpset::intersection_size(
            counts = T,
            mapping = aes(fill=Conserved)
          ) + ggplot2::scale_fill_manual(values=c(
            'Yes'='#E41A1C', 'No'='#377EB8'))
        ),
        width_ratio=0.1
      )


      #p1 = upset(x, order.by = "freq", point.size = 3,text.scale = c(1.25,1.25,1.25,1.25,1.25,1.25), mainbar.y.label = "Number of Known Motifs", line.size = 1)
      return(p1)
    }
    else{
      upset_all_data = list()
      for (j in c("conserved_only", "all_hits")){
        req = req_data[["prot"]][[database]][[j]][["no_comparison_prot"]]
        upsetr_expression = req$TP
        names(upsetr_expression) = req$title
        names(upsetr_expression) = stringr::str_replace_all(names(upsetr_expression), c("Enriched_Domain" = "Enriched Domain",
                                                                                        "PepSite" = "PepSite", "ClinVar_Path" = "ClinVar", "iELM_HMMs" = "iELM",
                                                                                        "qslim_predicted" = "QSLiMFinder"))

        if (remove_iELM){
          upsetr_expression = upsetr_expression[stringr::str_detect(names(upsetr_expression),"iELM", negate = T)]
        }
        upset_all_data[[j]] = upsetr_expression
      }
      common_filters = dplyr::intersect(names(upset_all_data$conserved_only), names(upset_all_data$all_hits))
      for (i in common_filters){
        upset_all_data$all_hits[i] = upset_all_data$all_hits[i] - upset_all_data$conserved_only[i]
      }
      for (j in 1:length(upset_all_data)){
        upset_all_data[[j]] = fromExpression(upset_all_data[[j]])
        upset_all_data[[j]]$Conserved = ifelse(names(upset_all_data)[j] == "conserved_only", "Yes", "No")
      }
      upset_all_data = dplyr::bind_rows(upset_all_data)
      intersect_cols = colnames(upset_all_data)[which(colnames(upset_all_data) != "Conserved")]

      p1 = ComplexUpset::upset(
        upset_all_data,
        intersect_cols,
        base_annotations = list(
          'Number of known motif proteins'= ComplexUpset::intersection_size(
            counts = T,
            mapping = aes(fill=Conserved)
          ) + ggplot2::scale_fill_manual(values=c(
            'Yes'='#E41A1C', 'No'='#377EB8'))
        ),
        width_ratio=0.1
      )
      #p1 = upset(x, order.by = "freq", point.size = 3,text.scale = c(1.25,1.25,1.25,1.25,1.25,1.25), mainbar.y.label = "Number of Known Motifs", line.size = 1)
      return(p1)
    }
  }
}

F2_A = upset_plot_motif(req_data = HV_all_eval, prot_only = T, prot_comparison = F, database = "ELM")
F2_B = upset_plot_motif(req_data = HV_all_eval, prot_only = T, prot_comparison = F, database = "PRMDB")
F2_C = All_levels_bar_plot(req_data = HV_all_eval, prot = T, prot_comparison = F, database = "ELM")
F2_D = All_levels_bar_plot(req_data = HV_all_eval, prot = T, prot_comparison = F, database = "PRMDB")


FS1_A = upset_plot_motif(req_data = HV_all_eval, prot_only = F, prot_comparison = F, database = "ELM")
FS1_B = upset_plot_motif(req_data = HV_all_eval, prot_only = F, prot_comparison = F, database = "PRMDB")
FS1_C = All_levels_bar_plot(req_data = HV_all_eval, prot = F, prot_comparison = F, database = "ELM")
FS1_D = All_levels_bar_plot(req_data = HV_all_eval, prot = F, prot_comparison = F, database = "PRMDB")
# Figure 3 ----------------------------------------------------------------
HV_final_hits = readRDS("data/output/HV_final_hits.rds")

hits_edges = HV_final_hits %>% dplyr::select(Pattern, Id, new_score)
hits_edges = hits_edges[!duplicated(hits_edges),]
hits_edges = hits_edges[!is.na(hits_edges$Id),]

recovered_edges = function(pattern, ELM_id){
    filtered = dplyr::filter(HV_final_hits, Pattern == pattern, Id == ELM_id)
    filtered = filtered[!duplicated(filtered),]

    if(!all(is.na(filtered$Accession))){
      return("Green")
    }
    else{
      return("Red")
    }
  }
generate_network_data_pattern_sim = function(remove_clv_trg = F){
  final_edges = hits_edges[which(hits_edges$new_score >= 0.67),]
  habal = final_edges %>% dplyr::select(Pattern, Id)
  habal = habal[!duplicated(habal),]

  for(i in 1:nrow(habal)){
    final_edges$edge_color[i] = recovered_edges(habal$Pattern[i], habal$Id[i])
  }

  final_edges$edge_weight = 0
  final_edges$edge_weight[which(final_edges$Score >= 0.5 & final_edges$Score < 0.6)] = 1
  final_edges$edge_weight[which(final_edges$Score >= 0.6 & final_edges$Score < 0.7)] = 2
  final_edges$edge_weight[which(final_edges$Score >= 0.8)] = 4
  final_edges$ELM_class =substr(final_edges$Id, 1,3)

  vertix_metadata = data.frame(V.name = unique(c(final_edges$Pattern, final_edges$Id)), node_shape = NA)
  vertix_metadata[which(vertix_metadata$V.name %in% final_edges$Pattern),]$node_shape = "Round Rectangle"
  vertix_metadata[which(vertix_metadata$V.name %in% final_edges$Id),]$node_shape = "Hexagon"

  if (remove_clv_trg){
    sample_network = final_edges[which(final_edges$ELM_class %nin% c("CLV", "TRG")),]
  }
  else{
    sample_network = final_edges
  }

  sample_v_metadata = vertix_metadata[which(vertix_metadata$V.name %in% c(sample_network$Pattern, sample_network$Id)),]

  g = graph_from_data_frame(d = sample_network, directed = F, vertices = sample_v_metadata)
  return(list(g, sample_network))
}

generate_heatmap_data_pattern_sim = function(compari_cutoff, remove_clv_trg = F){
  mat_data = hits_edges
  if (remove_clv_trg){
    mat_data$ELM_class =substr(mat_data$Id, 1,3)
    mat_data = mat_data[which(mat_data$ELM_class %nin% c("CLV", "TRG")),]
    mat_data = mat_data[,-4]
  }
  mat_data = mat_data[which(mat_data$new_score >= compari_cutoff),]
  mat_data = tidyr::spread(mat_data, key = Pattern, value = new_score)
  ids = mat_data$Id
  mat_data = mat_data[,-1]
  rownames(mat_data) = ids
  mat = as.matrix(mat_data)
  mat[is.na(mat)] = 0
  return(mat)
}

g = generate_network_data_pattern_sim(remove_clv_trg = T)
dd <- degree.distribution(g[[1]], cumulative=T, mode="all")
plot(dd, pch=19, cex=1, col="orange", xlab="Degree", ylab="Cumulative Frequency")

cytoscapePing()

createNetworkFromIgraph(g[[1]])

g_data = g[[2]]
table(g_data$edge_color)
length(unique(g_data$Id[which(g_data$edge_color == "Green")]))

x = g_data %>% group_by(Id) %>% summarise(n=n())
length(unique(g_data$Pattern[which(g_data$Id %in% c("LIG_PDZ_Class_1", "DOC_USP7_MATH_1", "DEG_SCF_TRCP1_1"))]))

# Figure 4 ----------------------------------------------------------------
tmp = tempfile()
curl::curl_download("https://www.ebi.ac.uk/biostudies/files/S-BSST668/Host_viral_Intact_PPI.xlsx", destfile = tmp)


host_viral_int = readxl::read_xlsx(tmp)
host_viral_symbols = data.frame(uniprot = c(host_viral_int$uniprot_A,host_viral_int$uniprot_B),
                                symbol = c(host_viral_int$Gene_A,host_viral_int$Gene_B),
                                Org = c(host_viral_int$org_A,host_viral_int$org_B))
host_viral_symbols = host_viral_symbols[!duplicated(host_viral_symbols),]

#pepsite_clinvar_filt = readxl::read_xlsx("../Evalutation results/updated_potential_candidates.xlsx", sheet = 2, skip = 1, na = c("", NA, "NA"))
pepsite_clinvar_filt = readRDS("data/output/HV_gold_insts.rds")
pepsite_clinvar_filt = pepsite_clinvar_filt %>% tidyr::separate(Dataset, into = c("V1","H1"), sep = "_")

generate_network_data_pepsite = function(compari_cutoff = 0.5, remove_not_provided_disease = T, pepsite_pval = 0.05,
                                         qslim.pval = 0.05, instance_type = c("Novel", "New_instance")){
  pepsite_net = pepsite_clinvar_filt
  pepsite_net = pepsite_net[which(pepsite_net$pepsite_p.val < pepsite_pval),]
  pepsite_net = pepsite_net[which(pepsite_net$qslim_pval < qslim.pval),]
  if (instance_type == "Novel"){
    pepsite_net = pepsite_net[which(pepsite_net$Instance_type == "Novel"),]
  }
  if (instance_type == "New_instance"){
    pepsite_net = pepsite_net[which(pepsite_net$Instance_type == "New_instance"),]
    pepsite_net = pepsite_net[which(pepsite_net$CompariMotif_norm_Score >= compari_cutoff),]
  }
  if (remove_not_provided_disease){
    pepsite_net = pepsite_net[which(pepsite_net$Disease_name %nin% "not provided"),]
  }

  pepsite_net = pepsite_net %>% dplyr::select(V1, H1, Pattern, motif_protein, GeneSymbol, ELM_Class, CompariMotif_norm_Score)
  pepsite_net = dplyr::left_join(pepsite_net, host_viral_symbols[,c(1,2)], by = c("H1" = "uniprot"))
  colnames(pepsite_net)[ncol(pepsite_net)] = "H1_sym"
  pepsite_net = dplyr::left_join(pepsite_net, host_viral_symbols[,c(1,2,3)], by = c("V1" = "uniprot"))
  colnames(pepsite_net)[c(9,10)] = c("V1_sym", "V1_org")
  pepsite_net$V1_sym = paste0(pepsite_net$V1_sym, "_", pepsite_net$V1_org)

  pattern_V1 = pepsite_net %>% dplyr::select(Pattern, V1_sym)
  pattern_V1 = pattern_V1[!duplicated(pattern_V1),]
  pattern_V1$edge_weight = 2
  pattern_V1$edge_line_type = "Solid"
  colnames(pattern_V1)[c(1,2)] = c("Int1","Int2")

  pattern_uniprot = pepsite_net %>% dplyr::select(Pattern, GeneSymbol)
  pattern_uniprot = pattern_uniprot[!duplicated(pattern_uniprot),]
  pattern_uniprot$edge_weight = 2
  pattern_uniprot$edge_line_type = "Dash"
  colnames(pattern_uniprot)[c(1,2)] = c("Int1","Int2")

  V1_H1 = pepsite_net %>% dplyr::select(V1_sym, H1_sym)
  V1_H1 = V1_H1[!duplicated(V1_H1),]
  V1_H1$edge_weight = 2
  V1_H1$edge_line_type = "Solid"
  colnames(V1_H1)[c(1,2)] = c("Int1","Int2")

  strain_info = dplyr::left_join(V1_H1, host_viral_int[,c(1,5,7)], by = c("Int1" = "uniprot_A"))
  strain_info = strain_info[!duplicated(strain_info),]

  H1_uniprot = pepsite_net %>% dplyr::select(H1_sym, GeneSymbol)
  H1_uniprot = H1_uniprot[!duplicated(H1_uniprot),]
  H1_uniprot$edge_weight = 2
  H1_uniprot$edge_line_type = "Solid"
  colnames(H1_uniprot)[c(1,2)] = c("Int1","Int2")

  if (instance_type == "New_instance"){
    Pattern_ELM = pepsite_net %>% dplyr::select(Pattern, ELM_Class, CompariMotif_norm_Score)
    Pattern_ELM = Pattern_ELM[order(Pattern_ELM$CompariMotif_norm_Score, decreasing = T),]
    Pattern_ELM = Pattern_ELM[!duplicated(Pattern_ELM),]
    colnames(Pattern_ELM) = c("Int1","Int2","edge_weight")
    Pattern_ELM$edge_weight = (Pattern_ELM$edge_weight * 5) + 2
    Pattern_ELM$edge_line_type = "Solid"
    final_edges = plyr::rbind.fill(pattern_V1, pattern_uniprot, V1_H1, H1_uniprot, Pattern_ELM)
  }
  else{
    final_edges = plyr::rbind.fill(pattern_V1, pattern_uniprot, V1_H1, H1_uniprot)
  }
  final_edges = final_edges[!duplicated(final_edges[,c(1,2)]),]
  vertix_metadata = data.frame(V.name = c(final_edges$Int1, final_edges$Int2), shape = NA, color = NA)

  vertix_metadata[which(vertix_metadata$V.name %in% pepsite_net$Pattern),]$shape = "Round Rectangle"
  vertix_metadata[which(vertix_metadata$V.name %in% c(pepsite_net$GeneSymbol,pepsite_net$H1_sym)),]$shape = "Ellipse"
  vertix_metadata[which(vertix_metadata$V.name %in% pepsite_net$V1_sym),]$shape = "Diamond"

  cols = RColorBrewer::brewer.pal(4,"Pastel1")
  vertix_metadata[which(vertix_metadata$V.name %in% pepsite_net$Pattern),]$color = cols[1]
  vertix_metadata[which(vertix_metadata$V.name %in% c(pepsite_net$GeneSymbol,pepsite_net$H1_sym)),]$color = cols[2]
  vertix_metadata[which(vertix_metadata$V.name %in% pepsite_net$V1_sym),]$color = cols[3]

  if (instance_type == "New_instance"){
    vertix_metadata[which(vertix_metadata$V.name %in% pepsite_net$ELM_Class),]$shape = "Octagon"
    vertix_metadata[which(vertix_metadata$V.name %in% pepsite_net$ELM_Class),]$color = cols[4]
  }

  vertix_metadata = vertix_metadata[!duplicated(vertix_metadata),]

  g = graph_from_data_frame(d = final_edges, directed = F, vertices = vertix_metadata)
  return(list(g, pepsite_net[!duplicated(pepsite_net),]))
}

g = generate_network_data_pepsite(instance_type = "Novel", compari_cutoff = 0.5, qslim.pval = 0.05, pepsite_pval = 0.05)
g_data = g[[2]]
length(unique(g_data$Pattern[which(g_data$GeneSymbol == "CFTR")]))

plot(g, layout = layout_nicely)

cytoscapePing()

createNetworkFromIgraph(g[[1]])

# Figure 5 ----------------------------------------------------------------
H1_GOBP = readRDS("data/output/plot_objects/GOBP_enrichment_H1.rds")
#H1_GOBP_cons = readRDS("data/output/plot_objects/GOBP_enrichment_H1_conserved.rds")
motif_GOBP = readRDS("data/output/plot_objects/GOBP_enrichment_motif.rds")
#motif_GOBP_cons = readRDS("data/output/plot_objects/GOBP_enrichment_motif_conserved.rds")
F5_H1 = enrichplot::dotplot(H1_GOBP) + theme(plot.margin = unit(c(-0.5,0,-0.3,-0.5), "cm")) + theme(axis.text.x = element_text(size = 10)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, vjust = 1, hjust = 1))

F5_motif = enrichplot::dotplot(motif_GOBP) + theme(plot.margin = unit(c(-0.5,0,-0.3,-0.5), "cm")) + theme(axis.text.x = element_text(size = 10)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, vjust = 1, hjust = 1))



# Figure 6  ---------------------------------------------------------
All_filters_dcGOR = readRDS("data/output/plot_objects/HV_all_filters_dcGOR.rds")
All_filters_dcGOR_cons = readRDS("data/output/plot_objects/All_filters_dcGOR_conserved.rds")
names(All_filters_dcGOR_cons) = paste0("Cons_", names(All_filters_dcGOR_cons))
All_filters_dcGOR_combined = list(All_filters_dcGOR, All_filters_dcGOR_cons)
All_filters_dcGOR_combined = unlist(All_filters_dcGOR_combined, recursive = F)

pfam_doms = readRDS("data/input/pfam_doms_human.rds")
pfam_doms = pfam_doms %>% tidyr::separate(PFAM, into = c("PFAM_ID", "PFAM_name"), sep = "--")

parse_dcGOR_eoutput = function(eout, adj_pval_cutoff = 0.05, join_pfam_name = F){
  final = eout@term_info
  for (i in 1:length(eout@overlap)){
    if (join_pfam_name){
      pfam_names = pfam_doms[which(pfam_doms$PFAM_ID %in% as.character(eout@overlap[[i]])),-1]
      pfam_names = pfam_names[!duplicated(pfam_names),]
      final$overlap[i] = paste(pfam_names$PFAM_name, collapse = "/")
    }
    else{
      final$overlap[i] = paste(eout@overlap[[i]], collapse = "/")
    }
  }
  final$zscore = as.numeric(eout@zscore)
  final$pvalue = as.numeric(eout@pvalue)
  final$adjp = as.numeric(eout@adjp)
  final = final[which(final$term_distance != 1),]
  return(final[which(final$adjp < adj_pval_cutoff),])
}


### Dot plot and Enrichment map plot
dcEnrich_to_compareClusRes = function(eout_filters, GO_algo = c("lea", "elim"), padjust_cutoff = 0.05, remove_iELM = F){
  replace_filter_names = c("Enriched_Domain" = "E", "PepSite" = "P", "ClinVar_Path" = "C", "iELM_HMMs" = "I", "&" = "", "Cons" = "CV")
  selected_eouts = list()
  gene_clusters = list()
  for (i in 1:length(eout_filters)){
    for (j in c("H1","Motif_prots")){
      if (!is.null(eout_filters[[i]][[j]])){
        if ("eoutput" %in% names(eout_filters[[i]][[j]])){
          if (class(eout_filters[[i]][[j]][["eoutput"]][[GO_algo]]) != "logical"){

            if (j == "H1"){
              selected_eouts[[paste0("H1_",names(eout_filters)[i])]] = eout_filters[[i]][[j]][["eoutput"]][[GO_algo]]
              gene_clusters[[paste0("H1_",names(eout_filters)[i])]] = eout_filters[[i]][[j]][["query"]]
            }
            else{
              selected_eouts[[paste0("DM_",names(eout_filters)[i])]] = eout_filters[[i]][[j]][["eoutput"]][[GO_algo]]
              gene_clusters[[paste0("DM_",names(eout_filters)[i])]] = eout_filters[[i]][[j]][["query"]]
            }
          }
        }
      }
    }
  }
  names(gene_clusters) = stringr::str_replace_all(names(gene_clusters), replace_filter_names)
  new_data = list()

  for (i in 1:length(selected_eouts)){
    parsed_eout = parse_dcGOR_eoutput(selected_eouts[[i]], adj_pval_cutoff = 10)
    parsed_eout = parsed_eout %>% dplyr::select(-term_namespace, -term_distance, -IC, -zscore)
    colnames(parsed_eout) = c("ID", "Description", "geneID", "pvalue", "p.adjust")
    count = as.numeric(lengths(selected_eouts[[i]]@overlap))
    Gene_ratio = paste0(count, "/", length(selected_eouts[[i]]@data))
    bg_ratio = paste0(as.numeric(lengths(selected_eouts[[i]]@anno)), "/", length(selected_eouts[[i]]@background))
    final = cbind(parsed_eout,data.frame(GeneRatio = Gene_ratio, BgRatio = bg_ratio, Count = count))
    final$Cluster = names(selected_eouts)[i]
    if (nrow(final) < 2 | max(final$pvalue) == 0){
      final$qvalue = 0
    }
    else{
      final$qvalue = qvalue::qvalue(p = final$pvalue, lambda = 0.05, pi0.method = "bootstrap")[["qvalues"]]
    }

    final = final %>% dplyr::select(Cluster, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count)
    new_data[[i]] = final
  }
  new_data = dplyr::bind_rows(new_data)
  new_data$Cluster = stringr::str_replace_all(new_data$Cluster, replace_filter_names)
  new_data = new_data[which(new_data$p.adjust < padjust_cutoff),]

  if (remove_iELM){
    new_data = new_data[stringr::str_detect(new_data$Cluster, "I", negate = T),]
    gene_clusters = gene_clusters[stringr::str_detect(names(gene_clusters), "I", negate = T)]
  }

  compare_cluster_obj = new("compareClusterResult", compareClusterResult = new_data, geneClusters = gene_clusters)
  return(compare_cluster_obj)
}
F6_A = dcEnrich_to_compareClusRes(All_filters_dcGOR_combined, GO_algo = "lea", remove_iELM = T)
F6_A = enrichplot::pairwise_termsim(F6_A)
F6_A = emapplot(F6_A, layout= "nicely", showCategory = 30) + theme(legend.position = "right") + scale_fill_manual(values =  RColorBrewer::brewer.pal(12, "Paired"))
F6_A

### PFAM SemSim
dcEnrich_PFAM_SemSim_GO = function(eout_filters, heatplot_title = NA,
                                   mot_prot_type = NULL, GO_algo = c("lea", "elim"),
                                   padjust_cutoff = 0.05, plot = T, remove_iELM = F, kmax = 20, bootstrap = 100){
  replace_filter_names = c("Enriched_Domain" = "E", "PepSite" = "P", "ClinVar_Path" = "C", "iELM_HMMs" = "I", "&" = "", "Cons" = "CV")
  selected_eouts = list()
  for (i in 1:length(eout_filters)){
    for (j in c("H1","Motif_prots")){
      if (!is.null(eout_filters[[i]][[j]])){
        if ("eoutput" %in% names(eout_filters[[i]][[j]])){
          if (class(eout_filters[[i]][[j]][["eoutput"]][[GO_algo]]) != "logical"){

            if (j == "H1"){
              selected_eouts[[paste0("H1_",names(eout_filters)[i])]] = eout_filters[[i]][[j]][["eoutput"]][[GO_algo]]
            }
            else{
              selected_eouts[[paste0("DM_",names(eout_filters)[i])]] = eout_filters[[i]][[j]][["eoutput"]][[GO_algo]]
            }
          }
        }
      }
    }
  }
  new_data = list()
  for (i in 1:length(selected_eouts)){
    parsed_eout = parse_dcGOR_eoutput(selected_eouts[[i]], adj_pval_cutoff = 10, join_pfam_name = F)
    parsed_eout = parsed_eout %>% dplyr::select(-term_namespace, -term_distance, -IC, -zscore)
    colnames(parsed_eout) = c("ID", "Description", "geneID", "pvalue", "p.adjust")
    count = as.numeric(lengths(selected_eouts[[i]]@overlap))
    Gene_ratio = paste0(count, "/", length(selected_eouts[[i]]@data))
    bg_ratio = paste0(as.numeric(lengths(selected_eouts[[i]]@anno)), "/", length(selected_eouts[[i]]@background))
    final = cbind(parsed_eout,data.frame(GeneRatio = Gene_ratio, BgRatio = bg_ratio, Count = count))
    final$Cluster = names(selected_eouts)[i]
    if (nrow(final) < 2 | max(final$pvalue) == 0){
      final$qvalue = 0
    }
    else{
      final$qvalue = qvalue::qvalue(p = final$pvalue, lambda = 0.05, pi0.method = "bootstrap")[["qvalues"]]
    }
    final = final %>% dplyr::select(Cluster, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count)
    new_data[[i]] = final
  }
  new_data = dplyr::bind_rows(new_data)
  new_data$Cluster = stringr::str_replace_all(new_data$Cluster, replace_filter_names)
  new_data = new_data[which(new_data$p.adjust < padjust_cutoff),]

  if (remove_iELM){
    new_data = new_data[stringr::str_detect(new_data$Cluster, "I", negate = T),]
  }

  if (!is.null(mot_prot_type)){
    if (mot_prot_type == "pepsite"){
      new_data = new_data[stringr::str_detect(new_data$Cluster, "P"),]
    }
    if (mot_prot_type == "EM"){
      new_data = new_data[stringr::str_detect(new_data$Cluster, "E"),]
    }
  }

  pfams = c()
  for (i in 1:nrow(new_data)){
    pfams = c(pfams, unlist(strsplit(new_data$geneID[i], split = "/")))
  }
  pfams = unique(pfams)
  pfams = pfams[which(pfams %in% pfam2GO$PFAM_ID)]
  pfam_annot = data.frame(ID = pfams, rm = NA)
  pfam_annot = dplyr::left_join(pfam_annot, pfam_doms[,-1], by = c("ID" = "PFAM_ID"))
  pfam_annot = pfam_annot[!duplicated(pfam_annot),]
  x = as.character(pfam_annot$PFAM_name)
  pfam_pairs = expand.grid(x, x)
  colnames(pfam_pairs) = c("PF1", "PF2")
  for (i in 1:nrow(pfam_pairs)){
    GO_1 = as.character(pfam2GO$ID[which(pfam2GO$PFAM_name == pfam_pairs$PF1[i])])
    GO_2 = as.character(pfam2GO$ID[which(pfam2GO$PFAM_name == pfam_pairs$PF2[i])])
    pfam_pairs$SemSim[i] = mgoSim(GO_1, GO_2, semData = hsGO, measure = "Wang", combine = "BMA")
  }
  pfam_mat = tidyr::pivot_wider(pfam_pairs, names_from = PF1, values_from = SemSim)
  rNames = pfam_mat$PF2
  pfam_mat = pfam_mat[,-1]
  rownames(pfam_mat) = rNames
  pfam_mat = as.matrix(pfam_mat)
  return(pfam_mat)
  # if (plot){
  #   pfam_heat = pheatmap::pheatmap(pfam_mat, main = heatplot_title, treeheight_row = 5, treeheight_col = 5, fontsize_row = 8)
  #   return(pfam_heat)
  # }
  # else{
  #   return(pfam_mat)
  # }
}
closest_to_center <- function(x, centers) {
  # compute squared euclidean distance from each sample to each cluster center
  tmp <- sapply(seq_len(nrow(x)),
                function(i) apply(centers, 1,
                                  function(v) sum((x[i, ]-v)^2)))
  colnames(tmp) = rownames(x)
  tmp = t(tmp)
  final = list()
  for (j in 1:ncol(tmp)){
    final[[colnames(tmp)[j]]] = rownames(tmp)[which(tmp[,j] == min(tmp[,j]))]
  }
  return(final)
}
dom_mat_kmeans = function(dom_mat, kmax = 10, bootstrap = 100, GoSim_thresh = 0.5, defined_K = F){
  counter = 0
  final_clusts = list()
  if (defined_K){
    km_clust = factoextra::eclust(dom_mat, "kmeans", k = max(kmax, nrow(dom_mat)), nstart = 25, nboot = bootstrap)
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
      if (nrow(x) <= 10){
        for (j in 2:nrow(x)){
          res = dom_mat_kmeans(dom_mat = x, kmax = j, defined_K = T)
          for (k in 1:length(res)){
            counter = counter + 1
            final_clusts[[paste0("C_", counter)]] = res[[k]]
          }
          break()
        }

      }
      else{
        res = dom_mat_kmeans(dom_mat = x, defined_K = F)
        for (j in 1:length(res)){
          counter = counter + 1
          final_clusts[[paste0("C_", counter)]] = res[[j]]
        }
      }
    }
    else{
      counter = counter + 1
      final_clusts[[paste0("C_",counter)]][["members"]] = names(km_clust$cluster[which(km_clust$cluster == i)])
      final_clusts[[paste0("C_",counter)]][["min_dist_center"]] = center_close[[as.character(i)]]
    }
  }
  return(final_clusts)

}


pfam2GO = read.table("data/input/dcGOR/pfam2GO_assoc.txt", header = T)
pfam2GO = dplyr::left_join(pfam2GO, pfam_doms[,-1])
pfam2GO = pfam2GO[!duplicated(pfam2GO),]
hsGO <- godata('org.Hs.eg.db', ont="BP")
pfam2GO = pfam2GO[which(pfam2GO$ID %in% names(hsGO@IC)),]
p2 = dcEnrich_PFAM_SemSim_GO(All_filters_dcGOR_combined, heatplot_title = "", GO_algo = "lea", bootstrap = 100, kmax = 10)
p2_clusts = dom_mat_kmeans(dom_mat = p2)

all_min_dist_GO = list()
for (i in 1:length(p2_clusts)){
  all_min_dist_GO[[i]] = data.frame(Clust = names(p2_clusts)[i],
                                 GO_ID = unique(pfam2GO$ID[which(pfam2GO$PFAM_name %in% p2_clusts[[i]]$min_dist_center)]))
}
all_min_dist_GO = dplyr::bind_rows(all_min_dist_GO)

GO_all = ontologyIndex::get_OBO("data/input/dcGOR/go-basic.obo")
GO_annot = data.frame(GO_all$name)
GO_annot$GO_ID = rownames(GO_annot)
colnames(GO_annot) = c("GO_name", "GO_ID")
all_min_dist_GO = dplyr::left_join(all_min_dist_GO, GO_annot)
all_min_dist_GO$label = 1

wide_mat = all_min_dist_GO %>% dplyr::select(Clust, GO_name, label)
wide_mat = tidyr::pivot_wider(wide_mat, names_from = Clust, values_from = label)
rNames = wide_mat$GO_name
wide_mat = as.matrix(wide_mat[,-1])
rownames(wide_mat) = rNames
wide_mat[is.na(wide_mat)] = 0
#colnames(wide_mat) = sub(".*_","",colnames(wide_mat))
F6_B = pheatmap::pheatmap(wide_mat, cluster_cols = F, cluster_rows = T,
                   border_color = "grey", cellwidth = 10,
                   color = c("#FFFFFF","#D73027"), treeheight_row = 2,
                   treeheight_col = 2, legend = F)


# fviz_dist(get_dist(p2$data, method = "spearman"))
# fviz_nbclust(p2$data, kmeans, method = "wss", k.max = 50, nstart = 25)
# gap_stat <- clusGap(p2$data, FUN = kmeans, nstart = 25,
#                     K.max = 50, B = 50)
#
# z = fviz_nbclust(p2$data, kmeans, method = "silhouette", k.max = 100, nstart = 25)
# res.km <- eclust(p2_data, "kmeans", k.max = 10, nstart = 25)
# res.hclust = eclust(p2$data, "hclust", nstart = 25, k.max = 20)
# fviz_gap_stat(p2$gap_stat)
# fviz_silhouette(p2)
# fviz_cluster(p2, repel = F,ellipse = T, ggtheme = theme_bw())
#
# fviz_cluster(res.hclust, geom = "point")

### Heatplots
dcEnrich_to_Enrichresult = function(eout_filters, GO_algo = c("lea", "elim"), padjust_cutoff = 0.05, remove_iELM = F){
  replace_filter_names = c("Enriched_Domain" = "E", "PepSite" = "P", "ClinVar_Path" = "C", "iELM_HMMs" = "I", "&" = "", "Cons" = "CV")
  selected_eouts = list()
  gene_sets = list()
  genes_2_sym = list()
  queries = list()
  bgs = list()
  for (i in 1:length(eout_filters)){
    for (j in c("H1", "Motif_prots")){
      if (!is.null(eout_filters[[i]][[j]])){
        if ("eoutput" %in% names(eout_filters[[i]][[j]])){
          if (class(eout_filters[[i]][[j]][["eoutput"]][[GO_algo]]) != "logical"){

            pfam_annot = data.frame(q = eout_filters[[i]][[j]][["query"]], rm = NA)
            pfam_annot = dplyr::left_join(pfam_annot, pfam_doms[,-1], by = c("q" = "PFAM_ID"))
            pfam_annot = pfam_annot[!duplicated(pfam_annot),]
            x = as.character(pfam_annot$PFAM_name)
            names(x) = as.character(pfam_annot$q)

            if (j == "H1"){
              selected_eouts[[paste0("H1_",names(eout_filters)[i])]] = eout_filters[[i]][[j]][["eoutput"]][[GO_algo]]
              gene_sets[[paste0("H1_",names(eout_filters)[i])]] = eout_filters[[i]][[j]][["eoutput"]][[GO_algo]]@anno
              genes_2_sym[[paste0("H1_",names(eout_filters)[i])]] = x

              queries[[paste0("H1_",names(eout_filters)[i])]] = eout_filters[[i]][[j]][["query"]]
              bgs[[paste0("H1_",names(eout_filters)[i])]] = eout_filters[[i]][[j]][["bg"]]
            }
            else{
              selected_eouts[[paste0("DM_",names(eout_filters)[i])]] = eout_filters[[i]][[j]][["eoutput"]][[GO_algo]]
              gene_sets[[paste0("DM_",names(eout_filters)[i])]] = eout_filters[[i]][[j]][["eoutput"]][[GO_algo]]@anno
              genes_2_sym[[paste0("DM_",names(eout_filters)[i])]] = x

              queries[[paste0("DM_",names(eout_filters)[i])]] = eout_filters[[i]][[j]][["query"]]
              bgs[[paste0("DM_",names(eout_filters)[i])]] = eout_filters[[i]][[j]][["bg"]]
            }
          }
        }
      }
    }
  }
  all_enrich_result = list()

  for (i in 1:length(selected_eouts)){
    parsed_eout = parse_dcGOR_eoutput(selected_eouts[[i]], adj_pval_cutoff = 10, join_pfam_name = T)
    parsed_eout = parsed_eout %>% dplyr::select(-term_namespace, -term_distance, -IC, -zscore)
    colnames(parsed_eout) = c("ID", "Description", "geneID", "pvalue", "p.adjust")
    count = as.numeric(lengths(selected_eouts[[i]]@overlap))
    Gene_ratio = paste0(count, "/", length(selected_eouts[[i]]@data))
    bg_ratio = paste0(as.numeric(lengths(selected_eouts[[i]]@anno)), "/", length(selected_eouts[[i]]@background))
    final = cbind(parsed_eout,data.frame(GeneRatio = Gene_ratio, BgRatio = bg_ratio, Count = count))
    final$Cluster = names(selected_eouts)[i]
    if (nrow(final) < 2 | max(final$pvalue) == 0){
      final$qvalue = 0
    }
    else{
      final$qvalue = qvalue::qvalue(p = final$pvalue, lambda = 0.05, pi0.method = "bootstrap")[["qvalues"]]
    }
    final = final %>% dplyr::select(Cluster, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count)
    final = final[which(final$p.adjust < padjust_cutoff),]

    if(nrow(final) != 0){
      en_res = new("enrichResult", result = final, pvalueCutoff = 0.05,
                   pAdjustMethod = "BH", qvalueCutoff = 0.2,
                   gene = queries[[i]], universe = bgs[[i]], geneSets = gene_sets[[i]],
                   organism = "Homo sapiens", keytype = "UNKNOWN", ontology = "GOBP", gene2Symbol = genes_2_sym[[i]],
                   readable = TRUE)
      all_enrich_result[[names(selected_eouts)[i]]] = en_res
    }
  }
  names(all_enrich_result) = stringr::str_replace_all(names(all_enrich_result), replace_filter_names)
  if (remove_iELM){
    all_enrich_result = all_enrich_result[stringr::str_detect(names(all_enrich_result), "I", negate = T)]
  }
  return(all_enrich_result)
}
heatplot_to_heatmap = function(enrich_res, cluster_info){
  mat_data = list()
  filt_dom_annot = list()
  for (i in 1:length(enrich_res)){
    res = enrich_res[[i]]@result
    per_filter = list()
    for (j in 1:nrow(res)){
      per_filter[[j]] = data.frame(GO_term = res$Description[j], domain = unlist(strsplit(res$geneID[j], "/")), filter = names(enrich_res)[i])
    }
    per_filter = dplyr::bind_rows(per_filter)
    filt_dom_annot[[i]] = data.frame(filter = names(enrich_res)[i], domain = unique(per_filter$domain))
    mat_data[[i]] = per_filter
  }
  mat_data = dplyr::bind_rows(mat_data)
  filt_dom_annot = dplyr::bind_rows(filt_dom_annot)
  dupl_dom_filt = unique(filt_dom_annot$domain[duplicated(filt_dom_annot$domain)])
  if (length(dupl_dom_filt) != 0){
    for (i in 1:length(dupl_dom_filt)){
      grouped_per_domain = mat_data[which(mat_data$domain == dupl_dom_filt[i]),]
      #grouped_per_domain = dplyr::left_join(filtered, mat_data)
      grouped_per_domain = grouped_per_domain %>% dplyr::group_by(filter) %>% dplyr::summarise(n=dplyr::n())
      remove_rows = which(filt_dom_annot$filter != grouped_per_domain$filter[which.max(grouped_per_domain$n)] & filt_dom_annot$domain == dupl_dom_filt[i])
      filt_dom_annot = filt_dom_annot[-remove_rows,]
    }
  }

  clust_annot = list()
  for (i in 1:length(cluster_info)){
    clust_annot[[i]] = data.frame(Cluster = names(cluster_info)[i], domain = unique(cluster_info[[i]]$members))
  }
  clust_annot = dplyr::bind_rows(clust_annot)

  filt_dom_annot = dplyr::left_join(filt_dom_annot,clust_annot)
  mat_data = dplyr::left_join(mat_data, clust_annot)

  filt_dom_annot = filt_dom_annot %>% dplyr::select(-domain)
  filt_dom_annot = filt_dom_annot[!duplicated(filt_dom_annot),]

  dupl_clust_filt = unique(filt_dom_annot$Cluster[duplicated(filt_dom_annot$Cluster)])
  if (length(dupl_clust_filt) != 0){
    for (i in 1:length(dupl_clust_filt)){
      grouped_per_clust = mat_data[which(mat_data$Cluster == dupl_clust_filt[i]),]
      #grouped_per_domain = dplyr::left_join(filtered, mat_data)
      grouped_per_clust = grouped_per_clust %>% dplyr::group_by(filter) %>% dplyr::summarise(n=dplyr::n())
      remove_rows = which(filt_dom_annot$filter != grouped_per_clust$filter[which.max(grouped_per_clust$n)] & filt_dom_annot$Cluster == dupl_clust_filt[i])
      filt_dom_annot = filt_dom_annot[-remove_rows,]
    }
  }

  #filt_dom_annot = filt_dom_annot[!duplicated(filt_dom_annot$domain),]

  rNames = filt_dom_annot$Cluster
  filt_dom_annot = as.data.frame(filt_dom_annot[,-2])
  rownames(filt_dom_annot) = as.character(rNames)
  colnames(filt_dom_annot) = "Filter"

  mat_data$label = 1
  mat_data = mat_data %>% dplyr::select(GO_term, Cluster, label)
  mat_data = mat_data[!duplicated(mat_data),]
  mat_data$Cluster = as.character(mat_data$Cluster)
  wide_mat = tidyr::pivot_wider(mat_data, names_from = Cluster, values_from = label)
  rNames = wide_mat$GO_term
  wide_mat = as.matrix(wide_mat[,-1])
  wide_mat[is.na(wide_mat)] = 0
  rownames(wide_mat) = rNames
  pheatmap(wide_mat, annotation_col = filt_dom_annot, border_color = "grey",
           annotation_names_col = F, legend = F,
           treeheight_col = 5, treeheight_row = 5, cluster_cols = F, show_colnames = T, color = c("#FFFFFF","#D73027"))

  # pheatmap(wide_mat, border_color = "grey",
  #          annotation_names_col = F, legend = F,
  #          treeheight_col = 5, treeheight_row = 5,fontsize_col = 8, cluster_cols = T, show_colnames = T)
}
debug(heatplot_to_heatmap)

p1 = dcEnrich_to_Enrichresult(All_filters_dcGOR_combined, GO_algo = "lea", remove_iELM = T)
F6_C = heatplot_to_heatmap(p1, cluster_info = p2_clusts)

# Figure 7 ----------------------------------------------------------------
new_H1_chembl = readRDS("data/output/HV_H1_chembl.rds")
#Sankey _drug
tmp = tempfile()
curl::curl_download("https://www.ebi.ac.uk/biostudies/files/S-BSST668/Host_viral_Intact_PPI.xlsx", destfile = tmp)
host_viral_int = readxl::read_xlsx(tmp)
human_human_int = read.csv("data/input/human_human_intact.csv.gz")
sankey_drug = HVSlimPred::Sankey_H1_chembl(new_H1_chembl$H1_and_Clinvar, host_viral_ints = host_viral_int, human_human_int = human_human_int, oxo_distance = 2)

out_dir = tempdir()
htmlwidgets::saveWidget(sankey_drug, file=file.path(out_dir,"Sankey_drug_2.html"))
webshot::webshot(file.path(out_dir,"Sankey_drug_2.html"), file.path(out_dir,"Sankey_drug.pdf"))

F7 = cowplot::ggdraw() +
  cowplot::draw_image(magick::image_read_pdf(file.path(out_dir,"Sankey_drug.pdf")))

# Figure S2 ---------------------------------------------------------------
shared_host_viral_eval = readRDS("data/output/shared_host_viral_eval.rds")
shared_human_eval = readRDS("data/output/shared_human_eval.rds")
host_viral_only = readRDS("data/output/HV_Complete_Evaluation.rds")
human_only = readRDS("data/output/human_Complete_Evaluation.rds")

clean_eval_scores = function(eval, dataset_type = c("ELM", "PRMdb")){
  if (dataset_type == "ELM"){
    prot = eval$ELM$Prot_level
    prot$F0.5 = (1.25 * prot$recall * prot$precision) / (0.25 * (prot$recall + prot$precision))
    motif = eval$ELM$HH_motif$All_F1_scores %>% dplyr::filter(weight_site_resid == "residue site_wise")
    motif$F0.5 = (1.25 * motif$Rc * motif$Pr) / (0.25 * (motif$Rc + motif$Pr))
    return(list("Prot" = prot, "Motif" = motif))
  }
  if (dataset_type == "PRMdb"){
    prot = eval$PRMdb$Prot_level
    prot$F0.5 = (1.25 * prot$recall * prot$precision) / (0.25 * (prot$recall + prot$precision))
    motif = eval$PRMdb$HH_motif$All_F1_scores %>% dplyr::filter(weight_site_resid == "residue site_wise")
    motif$F0.5 = (1.25 * motif$Rc * motif$Pr) / (0.25 * (motif$Rc + motif$Pr))
    return(list("Prot" = prot, "Motif" = motif))
  }
}
human_host_viral_boxplot = function(database){
  shared_hv_eval = clean_eval_scores(shared_host_viral_eval, dataset_type = database)
  shared_h_eval = clean_eval_scores(shared_human_eval, dataset_type = database)
  hv_only = clean_eval_scores(host_viral_only, dataset_type = database)
  h_only = clean_eval_scores(human_only, dataset_type = database)

  #Score = rep(c("F05_Mot", "F1_Prot", "PR_Mot", "PR_Prot", "RC_Mot", "RC_Prot"), each = 124)
  Score = rep(c("Motif F05", "Protein F1", "Motif Precision", "Protein Precision", "Motif Recall", "Protein Recall"), each = 124)
  Dataset = rep(c(rep("Host viral" , 62) , rep("Human only" , 62)) , 6)
  Value = c(shared_hv_eval$Motif$F0.5, shared_h_eval$Motif$F0.5, hv_only$Prot$F1, h_only$Prot$F1,
            shared_hv_eval$Motif$Pr, shared_h_eval$Motif$Pr, hv_only$Prot$precision, h_only$Prot$precision,
            shared_hv_eval$Motif$Rc, shared_h_eval$Motif$Rc, hv_only$Prot$recall, h_only$Prot$recall)

  data = data.frame(Score, Dataset ,  Value)
  new_order <- with(data, reorder(Score , Value, median , na.rm=T))

  par(mar=c(3,4,3,1))
  myplot <- boxplot(Value ~ Dataset*new_order , data=data  ,
                    boxwex=0.4 , ylab="Score",
                    main= paste0("Host viral Vs Human only", "-", database) ,
                    col=c("slateblue1" , "tomato"), xaxt="n")

  # To add the label of x axis
  my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
  my_names <- my_names[seq(1 , length(my_names) , 2)]
  # axis(1,
  #      at = seq(1.5 , 12 , 2),
  #      labels = my_names ,
  #      tick=FALSE , cex=0.3, cex.axis = 0.75)

  axis(1,labels = F, at = seq(1.5 , 12 , 2), tick = F)
  text(x = seq(1.5 , 12 , 2),
       y = par("usr")[3] - 0.013,
       labels = my_names,
       xpd = NA,
       ## Rotate the labels by 35 degrees.
       srt = 22.5,
       cex = 0.85, adj = 1)


  # Add the grey vertical lines
  for(i in seq(0.5 , 20 , 2)){
    abline(v=i,lty=1, col="grey")
  }

  # Add a legend only for ELM
  legend("topleft", legend = c("Host viral", "Human only"),
         col=c("slateblue1" , "tomato"),
         pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = c(0,0), text.width = 0.7, y.intersp = 0.5, x.intersp = 0.5)
}

p1 = cowplot::ggdraw() +
  cowplot::draw_image(magick::image_read("data/output/SLiMEnrich/Host_viral/Default/Histogram (1).png"))
p2 = cowplot::ggdraw() +
  cowplot::draw_image(magick::image_read("data/output/SLiMEnrich/Human_only/Default/Histogram (1).png"))


human_host_viral_boxplot("ELM")
p6 <- recordPlot()
plot.new() ## clean up device
p6 = cowplot::as_gtable(p6)
p6 = ggplotify::as.ggplot(p6)
human_host_viral_boxplot("PRMdb")
p7 <- recordPlot()
plot.new() ## clean up device
p7 = cowplot::as_gtable(p7)
p7 = ggplotify::as.ggplot(p7)

cowplot::plot_grid(p6,p1,p7,p2,nrow = 2, ncol = 2, rel_widths = c(1.2,1), labels = c("A", "C", "B", "D"))


# Figure S3 ---------------------------------------------------------------
HV_evaluation = readRDS("data/output/HV_Complete_Evaluation.rds")
HV_final_hits = readRDS("data/output/HV_final_hits.rds")

p1 = Heatmap_motif_metrics(HV_evaluation$ELM$HH_motif$all_metrics_per_motif, benchtype = "ELM")[[4]]

sankey_elm = Sankey_ELM(new_hits_ints = HV_final_hits)
htmlwidgets::saveWidget(sankey_elm, file="../Results/Plots_after_review/Sankey_ELM_2.html")
webshot::webshot("../Results/Plots_after_review/Sankey_ELM_2.html", "../Results/Plots_after_review/Sankey_ELM_2.pdf")

p2 = cowplot::ggdraw() +
  cowplot::draw_image(magick::image_read_pdf("../Results/Plots_after_review/Sankey_ELM_2.pdf"))
# cowplot::draw_plot_label("Host viral ELM", size = 12, hjust = -1.65, vjust = 3.5)
x = cowplot::plot_grid(p1,NULL, nrow = 2, rel_heights = c(1,0))
FS3 = cowplot::plot_grid(p2,x,ncol = 2, scale = c(1.25,1), labels = LETTERS[1:2], rel_widths = c(1.5,1))

# Figure S4 ---------------------------------------------------------------
#### Checking which new_score threshold works best
HV_final_hits = readRDS("data/output/HV_final_hits.rds")
y = c()
for (i in seq(0,1,0.01)){
  y = c(y,count_motif_instances(HV_final_hits[which(HV_final_hits$new_score >= i),]))
}
names(y) = seq(0,1,0.01)

score_insts = cbind(y,cutree(hclust(dist(y)), k = 5))
score_insts = as.data.frame(score_insts)
score_insts$x = rownames(score_insts)
rownames(score_insts) = NULL
score_insts = score_insts %>% dplyr::mutate_if(is.character, as.numeric)
score_insts$V2 = ifelse(score_insts$V2 == 3,2,score_insts$V2)
abline_idx = diff(score_insts$V2)
abline_idx = which(abline_idx != 0)
score_insts$V2 = as.character(score_insts$V2)
score_insts$V2 = stringr::str_replace_all(score_insts$V2, c("1" = "Lowest", "2" = "Low", "4" = "High", "5" = "Highest"))
colnames(score_insts)[which(colnames(score_insts) == "V2")] = "Similarity"
p1 = score_insts %>%
  ggplot2::ggplot(ggplot2::aes(x, y, color = Similarity)) +
  ggplot2::geom_point(size = 2.5) +
  # ggplot2::geom_vline(xintercept = score_insts$x[abline_idx]) +
  xlab("Normalized compariMotif") +
  ylab("Predicted Motif Instances") +
  theme_bw() +
  ggtitle("Normalized CompariMotif Score Cutoff") +
  theme(plot.title = element_text(hjust = 0.5))


### Breakdown of datasets between Slim and qslim
qslim_datasets_input = read.csv("data/input/qslimfinderoutperint.csv", na.strings = c("","NA"))
qslim_datasets_input = qslim_datasets_input[!duplicated(qslim_datasets_input),]
qslim_datasets_input = qslim_datasets_input[-3,]

HV_final_hits = readRDS("data/output/HV_final_hits.rds")
human_final_hits = readRDS("data/output/human_final_hits.rds")

slimfinder_datasets_input = readRDS("data/input/human_only_qslim_datasets.rds")

no_datasets = data.frame(Approach = rep(c("Human only", "Host viral"), each = 2),
                         Signif_Motifs = rep(c("Yes", "No"), 2),
                         No_datasets = c(length(unique(human_final_hits$Dataset)),
                                         length(slimfinder_datasets_input) - length(unique(human_final_hits$Dataset)),
                                         length(unique(HV_final_hits$Dataset)),
                                         length(unique(qslim_datasets_input$Dataset)) - length(unique(HV_final_hits$Dataset))))
totals <- no_datasets %>%
  dplyr::group_by(Approach) %>%
  dplyr::summarize(total = sum(No_datasets))

p2 = no_datasets %>%
  left_join(totals) %>%
  ggplot(aes(x = Approach, y = No_datasets, fill = Signif_Motifs)) +
  geom_bar(stat = "identity", width = 0.4) +
  geom_text(aes(Approach, total + 50, label = total, fill = NULL, vjust = 0), data = totals) +
  theme_ipsum() +
  theme(axis.title.x = element_text(size = 12, hjust = 0.5),
        axis.title.y = element_text(size = 12, hjust = 0.5),
        plot.title = element_text(size = 14, hjust = 0.5)) +
  ggtitle("No. Datasets per approach")



### Breakdown of number of motif instances across filters
#run all_data from prepapre_all_data for filters
all_data = readRDS("data/output/HV_all_filters_data.rds")
inst_per_filter = data.frame(Term = rep(all_data$all_hits$Metadata$Term, each = 2), No_insts = 0,
                             sym = rep(all_data$all_hits$Metadata$sym, each = 2))
inst_per_filter$Cons = rep(c("Yes","No"),15)

for (i in 1:nrow(inst_per_filter)){
  if (inst_per_filter$Cons[i] == "Yes"){
    inst_per_filter$No_insts[i] = count_motif_instances(all_data$conserved_only[[inst_per_filter$sym[i]]])
  }
  else{
    inst_per_filter$No_insts[i] = count_motif_instances(all_data$all_hits[[inst_per_filter$sym[i]]])
  }
}

inst_per_filter[31,] = c("QSLiMFinder", count_motif_instances(HV_final_hits), "-", "No")
inst_per_filter[32,] = c("QSLiMFinder", count_motif_instances(HV_final_hits[!is.na(HV_final_hits$Cons),]), "-", "Yes")

inst_per_filter = inst_per_filter[,-3]
inst_per_filter$Term = stringr::str_replace_all(inst_per_filter$Term, c("ClinVar_Path" = "ClinVar"))
inst_per_filter$Term = stringr::str_replace_all(inst_per_filter$Term, c("_" = " ", "&" = ","))
inst_per_filter$No_insts = as.numeric(inst_per_filter$No_insts)

# inst_per_filter$No_insts = (inst_per_filter$No_insts - min(inst_per_filter$No_insts)) /
#   (max(inst_per_filter$No_insts) - min(inst_per_filter$No_insts))


p3 = inst_per_filter %>%
  ggplot2::ggplot(ggplot2::aes(y=No_insts, x= reorder(Term, No_insts, max),fill = Cons, color = "blue", label = NULL)) +
  ggplot2::geom_bar(stat="identity", alpha=0.7, width=0.8, position = ggplot2::position_dodge(), colour = "black") +
  ggplot2::theme_bw() +
  ggplot2::coord_flip() +
  ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 14), axis.text.x = ggplot2::element_text(size = 12)) +
  ggplot2::ggtitle("Motif instances per filter") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))


### QC of top_ranked vs High IC
HV_main_out = read.csv("data/input/HV_allhitsperint.csv")
HV_occ = read.csv("data/input/HV_allhitsperint.occ.csv")

HV_cons_main_out = read.csv("data/input/HV_cons_allhitsperint.csv")
HV_cons_occ = read.csv("data/input/HV_cons_allhitsperint.occ.csv")

human_main_out = read.csv("data/input/human_allhitsint.csv")
human_occ = read.csv("data/input/human_allhitsint.occ.csv")
HV_cons_all_hits = preprocess_slimfinder_out(main_output = HV_cons_main_out, occ_out = HV_cons_occ)

HV_compari_output = readRDS("data/input/compari_motif_ELM_ids.RDS")
human_compari_output = read.table("data/input/human_allhits_slim_compari-compari_ELM_ids.compare.tdt", sep = "\t", header = T)

HV_final_hits_high_IC = parse_qslim(main_output = HV_main_out, occ_out = HV_occ, compari_output = HV_compari_output, preprocess_compari = T, select_top_rank = F)
HV_final_hits_top_rank = parse_qslim(main_output = HV_main_out, occ_out = HV_occ, compari_output = HV_compari_output, preprocess_compari = T, select_top_rank = T)
match_cols = c("Dataset", "Pattern", "uniprot", "Start_Pos", "End_Pos", "Match", "Cons")
x = HV_cons_all_hits$all_hitsint[,match_cols]
HV_final_hits = dplyr::left_join(HV_final_hits, x)

human_final_hits_high_IC = parse_qslim(main_output = human_main_out, occ_out = human_occ, compari_output = human_compari_output, preprocess_compari = T, select_top_rank = F)
human_final_hits_top_rank = parse_qslim(main_output = human_main_out, occ_out = human_occ, compari_output = human_compari_output, preprocess_compari = T, select_top_rank = T)

QC_redund_motifs = function(df_highIC, df_topRank, approach_title){
  common = plyr::match_df(df_highIC, df_topRank)
  common_motif_insts = count_motif_instances(common)

  common_known_insts = length(unique(common$Accession[!is.na(common$Accession)]))

  motif_insts = data.frame(Variant_type = c("Top_rank", "High_IC", "Both"),
                           Value = c(count_motif_instances(df_topRank) - common_motif_insts,
                                     count_motif_instances(df_highIC) - common_motif_insts,
                                     common_motif_insts))
  motif_insts$Approach = approach_title
  motif_insts$type = "Pred Motif Instances"
  known_insts = data.frame(Variant_type = c("Top_rank", "High_IC", "Both"),
                           Value = c(length(unique(df_topRank$Accession[!is.na(df_topRank$Accession)])) - common_known_insts,
                                     length(unique(df_highIC$Accession[!is.na(df_highIC$Accession)])) - common_known_insts,
                                     common_known_insts))
  known_insts$Approach = approach_title
  known_insts$type = "Known ELM Instances"

  df_highIC = df_highIC %>% dplyr::select(Pattern, uniprot, Start_Pos, End_Pos, Id, new_score)
  df_topRank = df_topRank %>% dplyr::select(Pattern, uniprot, Start_Pos, End_Pos, Id, new_score)

  df_highIC = df_highIC[!duplicated(df_highIC),]
  df_topRank = df_topRank[!duplicated(df_topRank),]

  df_highIC_scores = data.frame(Variant_type = "High_IC", new_score = df_highIC$new_score)
  df_topRank_scores = data.frame(Variant_type = "Top_rank", new_score = df_topRank$new_score)

  new_Score_distr = rbind.data.frame(df_highIC_scores, df_topRank_scores)
  new_Score_distr$Approach = approach_title

  return(list("motif_insts" = motif_insts, "known_insts" = known_insts, "score_distr" = new_Score_distr))
}

host_viral = QC_redund_motifs(df_highIC = HV_final_hits_high_IC, df_topRank = HV_final_hits_top_rank, approach_title = "Host Viral")
human_only = QC_redund_motifs(df_highIC = human_final_hits_high_IC, df_topRank = human_final_hits_top_rank, approach_title = "Human only")

counts_bar = rbind.data.frame(host_viral$motif_insts, host_viral$known_insts, human_only$motif_insts, human_only$known_insts)
score_distr = rbind.data.frame(host_viral$score_distr, human_only$score_distr)

p4 = score_distr %>%
  ggplot(aes(x=new_score, color=Variant_type, fill=Variant_type)) +
  geom_density(alpha=0.6) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum() +
  ylab("") +
  xlab("new_score") +
  facet_wrap(~Approach, scales = "free")

p5 = counts_bar %>%
  ggplot(aes(fill=Variant_type, y=Value, x=Approach)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~type, scales = "free") +
  ylab("Motif Instances")

### Boxplot evaluation comparable QSLiMFinder and SLiMFinder
shared_host_viral_eval = readRDS("data/output/shared_host_viral_eval_comparable_insts.rds")
shared_human_eval = readRDS("data/output/shared_human_eval_comparable_insts.rds")
host_viral_only = readRDS("data/output/HV_Complete_Evaluation.rds")
human_only = readRDS("data/output/human_Complete_Evaluation.rds")

p_qc = cowplot::plot_grid(p1,p2,p4,p5, labels = c("C","D","E","F"))

human_host_viral_boxplot("ELM")
p6 <- recordPlot()
plot.new() ## clean up device
p6 = cowplot::as_gtable(p6)
p6 = ggplotify::as.ggplot(p6)
human_host_viral_boxplot("PRMdb")
p7 <- recordPlot()
plot.new() ## clean up device
p7 = cowplot::as_gtable(p7)
p7 = ggplotify::as.ggplot(p7)

p8 = cowplot::plot_grid(p6,p7, nrow = 1, ncol = 2, labels = c("A","B"))
cowplot::plot_grid(p8, p_qc, rel_widths = c(0.5,1))

# Figure S5 ---------------------------------------------------------------
PS_info_all_filters_GO = readRDS("data/output/plot_objects/HV_PS_all_filters_GO.rds")
PS_info_all_filters_KEGG = readRDS("data/output/plot_objects/HV_PS_all_filters_KEGG.rds")

pfam_doms = readRDS("data/input/pfam_doms_human.rds")
pfam_doms = pfam_doms %>% tidyr::separate(PFAM, into = c("PFAM_ID", "PFAM_name"), sep = "--")


PSD_heatmap = function(orig,hm_title, dom_type = c("H1", "Motif_prots"), z_cutoff = 0){
  filt = orig[which(orig$dom_prot_type == dom_type),]
  filt$PS_z.score = scale(filt$PS_score)
  filt = filt[which(filt$PS_z.score > z_cutoff),]

  print(paste0("Number of possible filters is ", length(unique(filt$filter))))

  annot = pfam_doms[which(pfam_doms$PFAM_ID %in% filt$Accession),c(2,3)]
  annot = annot[!duplicated(annot),]

  long_score_mat = filt %>% dplyr::select(pathway, Accession, PS_z.score)
  long_score_mat = long_score_mat[!duplicated(long_score_mat),]
  long_score_mat = dplyr::left_join(long_score_mat, annot, by = c("Accession" = "PFAM_ID")) %>% dplyr::select(pathway, PFAM_name, PS_z.score)

  long_label_mat = filt %>% dplyr::select(pathway, Accession, filter)
  long_label_mat = long_label_mat[!duplicated(long_label_mat),]
  long_label_mat = long_label_mat %>% dplyr::group_by(pathway, Accession) %>% dplyr::summarise(count=dplyr::n())
  long_label_mat = as.data.frame(long_label_mat)
  long_label_mat = dplyr::left_join(long_label_mat, annot, by = c("Accession" = "PFAM_ID")) %>% dplyr::select(pathway, PFAM_name, count)

  wide_score_mat = tidyr::pivot_wider(long_score_mat, names_from = PFAM_name, values_from = PS_z.score)
  rNames = wide_score_mat$pathway
  wide_score_mat = wide_score_mat[,-1]
  rownames(wide_score_mat) = rNames
  wide_score_mat = as.matrix(wide_score_mat)
  wide_score_mat[is.na(wide_score_mat)] = 0

  wide_label_mat = tidyr::pivot_wider(long_label_mat, names_from = PFAM_name, values_from = count)
  rNames = wide_label_mat$pathway
  wide_label_mat = wide_label_mat[,-1]
  rownames(wide_label_mat) = rNames
  wide_label_mat = as.matrix(wide_label_mat)
  wide_label_mat[is.na(wide_label_mat)] = ""

  wide_label_mat = wide_label_mat[,colnames(wide_score_mat)]
  wide_label_mat = wide_label_mat[rownames(wide_score_mat),]

  pheatmap(wide_score_mat, show_colnames = T, cluster_cols = T, cluster_rows = T, display_numbers = wide_label_mat,
           number_color = "black", treeheight_row = 2, treeheight_col = 2, main = hm_title)
}

p1 = PSD_heatmap(orig = PS_info_all_filters_GO, dom_type = "H1", z_cutoff = 1, hm_title = "PSD H1 GO")
p2 = PSD_heatmap(orig = PS_info_all_filters_KEGG, dom_type = "H1", z_cutoff = 1, hm_title = "PSD H1 KEGG")

cowplot::plot_grid(p1[[4]], p2[[4]], label_size = 12, ncol = 1, nrow = 2, rel_heights = c(1.2,1), align = "hv", axis = "tblr")


# Figure S6 ---------------------------------------------------------------
All_filters_essentiality = readRDS("data/output/plot_objects/All_filters_essentiality.rds")
essen_plot_data = list()
counter = 0
for (i in 1:length(All_filters_essentiality)){
  for (j in 1:2){
    counter = counter + 1
    x = All_filters_essentiality[[i]][[j]]
    if (is.null(x)){
      next()
    }
    x$prot_type = sub("_.*", "", names(All_filters_essentiality[[i]])[j])
    x$filter = names(All_filters_essentiality)[i]
    x = x[,c(13, 12, 1:11)]
    essen_plot_data[[counter]] = x
  }
}
essen_plot_data = dplyr::bind_rows(essen_plot_data)
essen_plot_data$log_OR[which(essen_plot_data$TP == 0)] = 0

essen_plot_data$filter = stringr::str_replace_all(essen_plot_data$filter, "&", ",")
essen_plot_data$filter = stringr::str_replace_all(essen_plot_data$filter, "_", " ")
#replace_filter_names = c("Enriched Domain" = "E", "PepSite" = "P", "ClinVar Path" = "C", "iELM HMMs" = "I", "," = "")
replace_filter_names = c("ClinVar Path" = "ClinVar", "iELM HMMs" = "iELM")

essen_plot_data$filter = stringr::str_replace_all(essen_plot_data$filter, replace_filter_names)

# y_order = all_data$Metadata$Term
# y_order = stringr::str_replace_all(y_order, "&", ",")
# y_order = stringr::str_replace_all(y_order, "_", " ")
# y_order = stringr::str_replace_all(y_order, replace_filter_names)

plot_essentiality_bar = function(type = c("H1", "Motif_prots"), remove_iELM = F){
  if (type == "H1"){
    H1_essen_plot = essen_plot_data[which(essen_plot_data$prot_type == "H1"),]
    if (remove_iELM){
      H1_essen_plot = H1_essen_plot[stringr::str_detect(H1_essen_plot$filter, "iELM", negate = T),]
    }
    H1_essen_plot$signif_asterisk = ifelse(H1_essen_plot$p.val < 0.05, "*", "")
    H1_max_order = data.frame(filter = unique(H1_essen_plot$filter), max_OR = 0)
    for (i in 1:nrow(H1_max_order)){
      x = H1_essen_plot[which(H1_essen_plot$filter == H1_max_order$filter[i]),]
      H1_max_order$max_OR[i] = x$log_OR[which.max(x$log_OR)]
    }
    H1_max_order = H1_max_order$filter[order(H1_max_order$max_OR, decreasing = T)]
    H1_p = ggplot2::ggplot(H1_essen_plot, ggplot2::aes(y=log_OR, x=factor(filter, levels = rev(H1_max_order)), fill = CEN_clust)) +
      ggplot2::geom_bar(stat="identity", alpha=0.7, width=0.8, position = ggplot2::position_dodge(), colour = "black") +
      ggplot2::theme_bw() +
      ggplot2::coord_flip() +
      ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 14), axis.text.x = ggplot2::element_text(size = 14)) +
      #ggplot2::geom_hline(yintercept = xc_abline, color = "Red", linetype = "dashed") +
      ggplot2::ggtitle(paste0("H1 Essentiality")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::theme(legend.text=ggplot2::element_text(size=14)) +
      ggplot2::ylab("log2(odds ratio)")
    #ggfittext::geom_bar_text(position = "dodge")

    H1_p = H1_p + ggplot2::geom_text(data = H1_essen_plot,
                                     aes(x = factor(filter, levels = rev(H1_max_order)), group=CEN_clust, y = log_OR,
                                         label = format(signif_asterisk, nsmall = 0, digits=1, scientific = FALSE)),
                                     color="black", position=position_dodge(width = 0.9), hjust= 0, size = 8)
    H1_p = H1_p + ggplot2::geom_text(data = H1_essen_plot,
                                     aes(x = factor(filter, levels = rev(H1_max_order)), group=CEN_clust, y = log_OR,
                                         label = format(TP, nsmall = 0, digits=1, scientific = FALSE)),
                                     color="black", position=position_dodge(width = 0.8), hjust= 1)
    return(H1_p)
  }
  if (type == "Motif_prots"){
    H2_essen_plot = essen_plot_data[which(essen_plot_data$prot_type == "H2"),]
    if (remove_iELM){
      H2_essen_plot = H2_essen_plot[stringr::str_detect(H2_essen_plot$filter, "iELM", negate = T),]
    }
    H2_essen_plot$signif_asterisk = ifelse(H2_essen_plot$p.val < 0.05, "*", "")
    H2_max_order = data.frame(filter = unique(H2_essen_plot$filter), max_OR = 0)
    for (i in 1:nrow(H2_max_order)){
      x = H2_essen_plot[which(H2_essen_plot$filter == H2_max_order$filter[i]),]
      H2_max_order$max_OR[i] = x$log_OR[which.max(x$log_OR)]
    }
    H2_max_order = H2_max_order$filter[order(H2_max_order$max_OR, decreasing = T)]
    H2_p = ggplot2::ggplot(H2_essen_plot, ggplot2::aes(y=log_OR, x=factor(filter, levels = rev(H2_max_order)), fill = CEN_clust)) +
      ggplot2::geom_bar(stat="identity", alpha=0.7, width=0.8, position = ggplot2::position_dodge(), colour = "black") +
      ggplot2::theme_bw() +
      ggplot2::coord_flip() +
      ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 14), axis.text.x = ggplot2::element_text(size = 14)) +
      #ggplot2::geom_hline(yintercept = xc_abline, color = "Red", linetype = "dashed") +
      ggplot2::ggtitle(paste0("Motif Protein Essentiality")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::theme(legend.text=element_text(size=14)) +
      ggplot2::ylab("log2(odds ratio)")
    #ggfittext::geom_bar_text(position = "dodge")

    H2_p = H2_p + ggplot2::geom_text(data = H2_essen_plot,
                                     aes(x = factor(filter, levels = rev(H2_max_order)), group=CEN_clust, y = log_OR,
                                         label = format(signif_asterisk, nsmall = 0, digits=1, scientific = FALSE)),
                                     color="black", position=position_dodge(0.9), hjust= 0, size = 8)
    H2_p = H2_p + ggplot2::geom_text(data = H2_essen_plot,
                                     aes(x = factor(filter, levels = rev(H2_max_order)), group=CEN_clust, y = log_OR,
                                         label = format(TP, nsmall = 0, digits=1, scientific = FALSE)),
                                     color="black", position=position_dodge(0.8), hjust= 1)
    return(H2_p)
  }
}

p1 = plot_essentiality_bar("H1", remove_iELM = T)
p2 = plot_essentiality_bar("Motif_prots", remove_iELM = T)

x = cowplot::plot_grid(p1, NULL, labels = c("A", ""), label_size = 12, ncol = 2, nrow = 1, rel_widths = c(4,1))
y = cowplot::plot_grid(p2, NULL, labels = c("b", ""), label_size = 12, ncol = 2, nrow = 1, rel_widths = c(4,1))

cowplot::plot_grid(p1, p2, labels = c('', ''),ncol = 2, nrow = 1)

# Figure S7 ---------------------------------------------------------------
HV_evaluation = readRDS("data/output/HV_Complete_Evaluation.rds")
p1 = Heatmap_all_evaluation(HV_evaluation, benchtype = "ELM")[[4]]
p2 = Heatmap_all_evaluation(HV_evaluation, benchtype = "PRMdb")[[4]]
cowplot::plot_grid(p1,p2, ncol = 2, labels = LETTERS[1:2])
