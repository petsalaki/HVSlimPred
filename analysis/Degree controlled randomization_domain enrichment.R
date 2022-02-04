require(igraph)
require(snow)
require(snowfall)
require(VertexSort)

pfam_doms = readRDS("data/input/pfam_doms_human.rds")
human_fasta_seqs = readRDS("data/input/Human_prots_fasta.rds")
human_human_intact = read.csv("data/input/human_human_intact.csv.gz")
HV_cdhit_out = readRDS("data/input/HV_cdhit_output.rds")
human_cdhit_out = readRDS("data/input/human_only_cdhit_output.rds")
HV_final_hits = readRDS("data/output/HV_final_hits.rds")
human_final_hits = readRDS("data/output/human_final_hits.rds")

pfam_doms_igraph = graph_from_data_frame(pfam_doms)

randomized_pfam_doms = VertexSort::dpr(vgraph = pfam_doms_igraph, viteration_no = 1000, vparallel = F, vcpus = 1)

splitted = human_fasta %>% tidyr::separate(Header, into = c("rm_1", "uniprot", "name"), sep = "[|]") %>% dplyr::select(-rm_1) %>%
  tidyr::separate(name, into = c("uniprot_name", "rm_2"), sep = " ") %>% dplyr::select(-rm_2)
splitted$Header = human_fasta$Header
splitted$Desc = sub(".*HUMAN ","",splitted$Header)
splitted = splitted[which(splitted$uniprot %in% human_prots),]


get_human_interactors = function(human_human_intact, uniprot_id, header_mapping, return_fasta = T){
  ints_A = dplyr::filter(human_human_intact, uniprot_A == uniprot_id)
  ints_B = dplyr::filter(human_human_intact, uniprot_B == uniprot_id)
  if (nrow(ints_A) != 0){
    ints_A = unique(ints_A$uniprot_B)
  }
  else{
    ints_A = c()
  }

  if (nrow(ints_B) != 0){
    ints_B = unique(ints_B$uniprot_A)
  }
  else{
    ints_B = c()
  }
  all_ints = unique(c(ints_A, ints_B))

  if (return_fasta & !missing(header_mapping)){
    headers = header_mapping[which(header_mapping$uniprot %in% all_ints),]
    headers[["Header"]] = paste0(headers[["uniprot_name"]], "__", headers[["uniprot"]], " ", headers[["Desc"]])
    #final_fasta = humans_all_fasta[which(humans_all_fasta$Header %in% headers),]
    final_fasta = dplyr::select(headers, Header, Sequence)
    return(list("All_ints" = all_ints, "fasta_file" = final_fasta))
  }
  else{
    return(all_ints)
  }
}
Domains_in_H2_ints = function(H2_ints,pfam_doms_human = NULL){
  pfam_doms_list = pfam_doms_human
  domains = pfam_doms_list[which(pfam_doms_list$uniprot %in% H2_ints),]
  if (nrow(domains) != 0){
    return(domains)
  }
  else{
    return(NA)
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
calculate_empirical_pval = function(mot_carprots, cdhit_output, pfam_doms_list){
  final_res = list()
  for (i in 1:length(mot_carprots)){
    interactors = get_human_interactors(human_human_intact = human_human_int,mot_carprots[i], header_mapping = splitted)
    if (is.null(interactors$All_ints)){
      next()
    }
    else{
      cdhit_out_prot = cdhit_output[[mot_carprots[i]]]
      H2_int_prots = sub("[ ].*","",sub(".*[__]","",cdhit_out_prot$Header))

      nr_domains_H2 = Domains_in_H2_ints(H2_ints = H2_int_prots, pfam_doms_human = pfam_doms_list)
      if (is.na(nr_domains_H2)){
        next()
      }
      else{
        int_mot = interactors$All_ints

        domains = unique(nr_domains_H2$PFAM)
        final_TPs = list()
        for (j in 1:length(domains)){
          domlist_has_domain = pfam_doms_list$uniprot[which(pfam_doms_list$PFAM == domains[j])]
          int_domain = length(intersect(domlist_has_domain, int_mot))


          final_TPs[[j]] = data.frame(motif_protein = mot_carprots[i], PFAM = domains[j], TP = int_domain)
        }
        final_TPs = dplyr::bind_rows(final_TPs)

        final_res[[i]] = final_TPs
      }
    }
  }
  final_res = dplyr::bind_rows(final_res)
  return(final_res)
}


HV_motif_prots = unique(HV_final_hits$uniprot)
human_motif_prots = unique(human_final_hits$uniprot)

rand_pred_output_path = tempdir() # Replace with output path to store enrichment of each random network

for (i in 1:length(HV_motif_prots)){
  rand_network = igraph::as_data_frame(randomized_pfam_doms[[i]])
  colnames(rand_network) = c("uniprot","PFAM")
  rand_network = rand_network[!duplicated(rand_network),]

  tps = calculate_empirical_pval(mot_carprots = HV_motif_prots, cdhit_output = HV_cdhit_out,
                                 pfam_doms_list = rand_network)
  tps$network_idx = paste0("network_", i)
  saveRDS(tps, paste0(rand_pred_output_path,"rand_network_", i))
}

HV_enriched_doms = readRDS("data/output/HV_enriched_doms.rds")

observed_enrich = HV_enriched_doms$Enriched_domains
observed_enrich = mapply(cbind, observed_enrich, "uniprot"=names(observed_enrich), SIMPLIFY=F)
observed_enrich = dplyr::bind_rows(observed_enrich)
observed_enrich = observed_enrich %>% dplyr::select(uniprot, PFAM, TP)
colnames(observed_enrich) = c("motif_prot", "PFAM", "TP_obs")

for (i in 1:1000){
  rand_result = readRDS(paste0(rand_pred_output_path, "rand_network_", i))
  rand_result = rand_result[,-4]
  colnames(rand_result)[c(1,2)] = c("motif_prot", "PFAM")
  colnames(rand_result)[3] = paste0("TP_rand_", i)
  observed_enrich = dplyr::left_join(observed_enrich, rand_result, by = c("motif_prot", "PFAM"))
  observed_enrich[[paste0("TP_rand_", i)]][is.na(observed_enrich[[paste0("TP_rand_", i)]])] = 0
  observed_enrich[[paste0("TP_rand_", i)]] = ifelse(observed_enrich[[paste0("TP_rand_", i)]] >= observed_enrich[["TP_obs"]],1,0)
}

rand_emppval_res = observed_enrich
rand_emppval_res$emp_pval = rowSums(rand_emppval_res[,c(4:1003)])
rand_emppval_res$emp_pval = rand_emppval_res$emp_pval / 1000

final_emppval = rand_emppval_res %>% dplyr::select(motif_prot, PFAM, emp_pval)
colnames(final_emppval)[c(1,2)] = c("uniprot", "PFAM")

for (i in 1:length(HV_enriched_doms$Enriched_domains)){
  filtered = final_emppval[which(final_emppval$uniprot == names(HV_enriched_doms$Enriched_domains)[i]),c(2,3)]
  HV_enriched_doms$Enriched_domains[[i]] = dplyr::left_join(HV_enriched_doms$Enriched_domains[[i]], filtered)
}
saveRDS(HV_enriched_doms, "data/output/HV_enriched_doms_with_emppval.rds")
