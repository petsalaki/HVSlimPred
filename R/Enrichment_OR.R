#' Domain enrichment of motif-carrying proteins
#'
#' @author Bishoy Wadie, Evangelia Petsalaki
#' @return
#' A list of non_redundant domains and their enrichment scores for each motif_carrying protein.
#' @export
#' @importFrom utils read.delim read.table read.csv
Domain_Enrichment = function(human_fasta, human_human_int, mot_carprots, cdhit_output = NULL, cdhit_outdir = NULL,
                             pfam_doms_list = NULL, cdhit_path = NULL,
                             cdhit_ident_thresh = 0.7, cdhit_no_cores = 1, cdhit_wsl = T, keep_cdhit_out = T){
  #human_fasta = download_data(path = "Enrichment/Human_prots_fasta.rds", is.rds = T)
  #human_human_int = read.csv(download_data(path = "Enrichment/human_human_intact.csv.gz", is.rds = F))
  human_prots = unique(c(human_human_int$uniprot_A, human_human_int$uniprot_B))

  splitted = human_fasta %>% tidyr::separate(Header, into = c("rm_1", "uniprot", "name"), sep = "[|]") %>% dplyr::select(-rm_1) %>%
    tidyr::separate(name, into = c("uniprot_name", "rm_2"), sep = " ") %>% dplyr::select(-rm_2)
  splitted$Header = human_fasta$Header
  splitted$Desc = sub(".*HUMAN ","",splitted$Header)
  splitted = splitted[which(splitted$uniprot %in% human_prots),]

  pfam_indiv_count_nr = list()
  Enriched_domains = list()
  for (i in 1:length(mot_carprots)){
    interactors = get_human_interactors(human_human_intact = human_human_int,uniprot_id = mot_carprots[i], header_mapping = splitted)
    if (is.null(interactors$All_ints)){
      pfam_indiv_count_nr[[mot_carprots[i]]] = NULL
      Enriched_domains[[mot_carprots[i]]] = NULL
      next()
    }
    else{
      if (is.null(cdhit_path) & is.null(cdhit_output)){
        stop("Error : At least CD-hit path or CD-hit output should be provided")
      }
      else if (is.null(cdhit_path) & !is.null(cdhit_output)){
        cdhit_out_prot = cdhit_output[[mot_carprots[i]]]
        H2_int_prots = sub("[ ].*","",sub(".*[__]","",cdhit_out_prot$Header))
      }
      else{
        input_fasta = create_output_path(cdhit_outdir, file_name = paste0(mot_carprots[i], ".fa"))
        microseq::writeFasta(interactors$fasta_file, out.file = input_fasta)
        if (keep_cdhit_out & is.null(cdhit_outdir)){
          stop("Error : Please provide a working cdhit output directory path, otherwise set keep_cdhit_out = False")
        }
        else if (!keep_cdhit_out & is.null(cdhit_outdir)){
          cdhit_outdir = tempdir()
          Run_cdhit(path_to_cdhit = cdhit_path, outdir = cdhit_outdir, path_to_fasta = input_fasta,
                    identity_thresh = cdhit_ident_thresh, num_cores = cdhit_no_cores, use_wsl = use_wsl)

        }
        else if (keep_cdhit_out & !is.null(cdhit_outdir)){
          Run_cdhit(path_to_cdhit = cdhit_path, outdir = cdhit_outdir, path_to_fasta = input_fasta,
                    identity_thresh = cdhit_ident_thresh, num_cores = cdhit_no_cores, use_wsl = use_wsl)
        }
        output_fasta = create_output_path(cdhit_outdir, file_name = paste0(mot_carprots[i], "_nr.fa"))
        if (file.exists(output_fasta)){
          cdhit_out = microseq::readFasta(output_fasta)
          H2_int_prots = sub("[ ].*","",sub(".*[__]","",cdhit_out$Header))
        }
        else{
          stop(paste0("Error ", output_fasta, " doesn't exist", " please check file name"))
        }
      }
      nr_domains_H2 = Domains_in_H2_ints(H2_ints = H2_int_prots, pfam_doms_human = pfam_doms_list)
      if (is.na(nr_domains_H2)){
        pfam_indiv_count_nr[[mot_carprots[i]]] = NULL
        Enriched_domains[[mot_carprots[i]]] = NULL
        next()
      }
      else{
        pfam_indiv_count_nr[[mot_carprots[i]]] = data.frame(uniprot = mot_carprots[i], PFAM = unique(nr_domains_H2$PFAM),
                                                            total_domain_ints = length(unique(nr_domains_H2$uniprot)))

        Enriched_domains[[mot_carprots[i]]] = Enrich_domain_fisher(mot_car_prot = mot_carprots[i],
                                                                   pfam_domcount_nr = nr_domains_H2,
                                                                   H2_ints = interactors$All_ints, pfam_doms_human = pfam_doms_list)
      }
    }
  }
  # pfam_indiv_count_nr = dplyr::bind_rows(pfam_indiv_count_nr)
  # Enriched_domains = dplyr::bind_rows(Enriched_domains)
  return(list("pfam_doms_nr" = pfam_indiv_count_nr, "Enriched_domains" = Enriched_domains))
}

Run_cdhit = function(path_to_cdhit, path_to_fasta, outdir, identity_thresh = 0.7, num_cores = 1, use_wsl = T){
  if (use_wsl){
    outdir_system_path = system(paste("bash -c", shQuote(paste0("wslpath -u ", outdir))), intern = T)
    system_path_fasta = system(paste("bash -c", shQuote(paste0("wslpath -u ", path_to_fasta))), intern = T)
  }
  else{
    outdir_system_path = outdir
    system_path_fasta = path_to_fasta
  }

  output_filename = paste0(sub("[.].*","",sub(".*/","", path_to_fasta)), "_nr.fa")
  output_file = create_output_path(outdir_path = outdir_system_path, file_name = output_filename)

  command = paste0(path_to_cdhit, " -i ", system_path_fasta, " -o ", output_file, " -c ", identity_thresh, " -T ", num_cores)

  system(paste("bash -c", shQuote(command)), ignore.stdout = T)
}
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
  if (is.null(pfam_doms_human)){
    # TODO change path of input data
    pfam_doms_list = download_data(path = "Enrichment/pfam_doms_human.rds", is.rds = T)
  }
  else{
    pfam_doms_list = pfam_doms_human
  }
  domains = pfam_doms_list[which(pfam_doms_list$uniprot %in% H2_ints),]
  if (nrow(domains) != 0){
    return(domains)
  }
  else{
    return(NA)
  }
}
Enrich_domain_fisher = function(mot_car_prot, pfam_domcount_nr, H2_ints, pfam_doms_human = NULL){
  if (is.null(pfam_doms_human)){
    # TODO change the path of input data
    pfam_doms_list = download_data(path = "Enrichment/pfam_doms_human.rds", is.rds = T)
  }
  else{
    pfam_doms_list = pfam_doms_human
  }
  int_mot = H2_ints

  domains = unique(pfam_domcount_nr$PFAM)
  final = list()
  for (i in 1:length(domains)){
    domlist_has_domain = pfam_doms_list$uniprot[which(pfam_doms_list$PFAM == domains[i])]
    domlist_no_domain = pfam_doms_list$uniprot[which(pfam_doms_list$PFAM != domains[i])]

    int_domain = length(intersect(domlist_has_domain, int_mot))
    not_int_domain = length(unique(domlist_has_domain[which(domlist_has_domain %nin% int_mot)]))
    int_no_domain = length(intersect(domlist_no_domain, int_mot))
    not_int_no_domain = length(unique(domlist_no_domain[which(domlist_no_domain %nin% int_mot)]))

    fisher = stats::fisher.test(matrix(c(int_domain, not_int_domain, int_no_domain, not_int_no_domain), nrow = 2, ncol = 2,
                                byrow = T), alternative = "two.sided")


    final[[i]] = data.frame(PFAM = domains[i], TP = int_domain, FP = not_int_domain, FN = int_no_domain, TN = not_int_no_domain,
                            p.val = fisher$p.value, OR = as.numeric(fisher$estimate), conf_int_lb = fisher$conf.int[1],
                            conf_int_ub = fisher$conf.int[2])
  }
  final = dplyr::bind_rows(final)
  final$adj.p.val = stats::p.adjust(final$p.val, method = "BH")
  return(final)
}

#TODO Add function to calculate the empirical pval from randomization results



#' Motif carrying protein enrichment all filters
#'
#' @author Bishoy Wadie, Evangelia Petsalaki
#' @return
#' A dataframe with odds ratio of each filter.
#' @export
#' @importFrom utils read.delim read.table read.csv
All_prot_OR = function(EM_list, Clinvar_path, iELM_insts, new_hits_ints, pepsite_res, bench_data = NULL, min_dom = 5, adj_pval = 0.05, qslim_input = NULL, comparison = T, bg_possible = T, prot_only = T){
  required_data = prepare_required_data_OR(EM = EM_list, clin_path = Clinvar_path, iELM = iELM_insts, new_hits_ints = new_hits_ints, min_domain = min_dom,
                                           adj.pval = adj_pval, slim_input_prots = qslim_input, pepsite_res = pepsite_res, prot_only = prot_only)
  metadata = required_data[["Metadata"]]
  if (!comparison){
    B = all_levels_prot("all_pred_prots", required = required_data, bg = "qslim_input", consider_possible_only = bg_possible, bench_data = bench_data)
    B$title = "qslim_predicted"
    C = all_levels_prot("A", required = required_data ,bg = "qslim_input", consider_possible_only = bg_possible, bench_data = bench_data)
    C$title = metadata$Term[which(metadata$sym == "A")]
    D = all_levels_prot("AC", required = required_data, bg = "qslim_input", consider_possible_only = bg_possible, bench_data = bench_data)
    D$title = metadata$Term[which(metadata$sym == "AC")]
    E = all_levels_prot("ABC", required = required_data, bg = "qslim_input", consider_possible_only = bg_possible, bench_data = bench_data)
    E$title = metadata$Term[which(metadata$sym == "ABC")]

    G = all_levels_prot("ABCD", required = required_data, bg = "qslim_input", consider_possible_only = bg_possible, bench_data = bench_data)
    G$title = metadata$Term[which(metadata$sym == "ABCD")]

    All_levels_result = rbind.data.frame(B,C,D,E,G)
    return(All_levels_result)
  }
  else{
    All_levels_result = list()
    for (i in 1:nrow(metadata)){
      result = all_levels_prot(metadata$sym[i], required = required_data, bg = "qslim_pred", consider_possible_only = bg_possible, bench_data = bench_data)
      result$title = metadata$Term[which(metadata$sym == metadata$sym[i])]
      All_levels_result[[i]] = result
    }
    All_levels_result = dplyr::bind_rows(All_levels_result)
    for (i in 1:nrow(All_levels_result)){
      if (All_levels_result$TP[i] == 0){
        All_levels_result$odds_ratio[i] = 0
      }
      else{
        next()
      }
    }
    return(All_levels_result)
  }
}

prepare_required_data_OR = function(EM, clin_path, iELM, new_hits_ints,pepsite_res, min_domain = 5, adj.pval = 0.05, slim_input_prots = NULL, prot_only = T, cons = F){
  #ELM_PRM = download_data("Enrichment/ELM_PRM_bg_OR.xlsx", is.rds = F)
  #ELM = list("ELM_possible" = read_xlsx(ELM_PRM, sheet = "ELM_Human_Possible"), "ELM_Human" = read_xlsx(ELM_PRM, sheet = "ELM_Human"))
  #PRMdb = list("PRM_possible" = read_xlsx(ELM_PRM, sheet = "PRMdb_Human_Possible"), "PRM_Human" = read_xlsx(ELM_PRM, sheet = "PRMdb_Human"))
  #human_human_int = read.csv(download_data(path = "Enrichment/human_human_intact.csv.gz"))

  if (cons){
    new_hits_ints = new_hits_ints[!is.na(new_hits_ints$Cons),]
  }

  colnames(pepsite_res)[which(colnames(pepsite_res) %in% c("motif_protein", "start","end"))] = c("uniprot", "Start_Pos", "End_Pos")

  H1_prots = sub(".*_", "", unique(new_hits_ints$Dataset))

  all_pred_prots = unique(new_hits_ints$uniprot)
  qslim_instances = dplyr::select(new_hits_ints, uniprot, Start_Pos, End_Pos)
  qslim_instances = qslim_instances[!duplicated(qslim_instances),]

  if(is.null(slim_input_prots)){
    # TODO change input path
    qslim_input = download_data(path = "Evaluation_data/proteins_per_dataset_allhitsint.RDS", is.rds = T)
    qslim_input = unique(unlist(qslim_input))
  }
  else{
    qslim_input = slim_input_prots
  }

  A = EM$Enriched_domains
  A = mapply(cbind, A, "uniprot"=names(A), SIMPLIFY=F)
  A = dplyr::bind_rows(A)
  if ("emp_pval" %in% names(A)){
    A = A %>% dplyr::filter(TP >= min_domain, adj.p.val < adj.pval, emp_pval < 0.05)
  }
  else{
    A = A[which(A$TP >= min_domain & A$adj.p.val < adj.pval),]
  }
  A = A %>% tidyr::separate(PFAM, into = c("pfam_id", "pfam_name"), sep = "--")

  x = A %>% dplyr::select(pfam_id, uniprot)
  x = x[!duplicated(x),]
  A = dplyr::left_join(x,qslim_instances, by = "uniprot")
  A = A[!is.na(A$Start_Pos),]

  B = plyr::match_df(clin_path, qslim_instances)
  C = plyr::match_df(pepsite_res, qslim_instances)
  D = plyr::match_df(iELM, qslim_instances)

  # B = clin_path
  # C = pepsite_res
  # D = iELM
  colnames(D)[1] = "domain_protein"

  AB = plyr::match_df(B, A, on = c("uniprot", "Start_Pos", "End_Pos"))
  AC = plyr::match_df(C, A, on = c("uniprot", "Start_Pos", "End_Pos","pfam_id"))
  AD = plyr::match_df(D, A, on = c("uniprot", "Start_Pos", "End_Pos"))

  if (prot_only){
    BC = B[which(B$uniprot %in% C$uniprot),]
    BD = B[which(B$uniprot %in% D$uniprot),]
    CD = plyr::match_df(C, D, on = c("domain_protein", "uniprot"))
    ABC = BC[which(BC$uniprot %in% A$uniprot),]
    ABD = BD[which(BD$uniprot %in% A$uniprot),]
    ACD = plyr::match_df(AC, D, on = c("domain_protein", "uniprot"))
    BCD = CD[which(CD$uniprot %in% B$uniprot),]
    ABCD = ACD[which(ACD$uniprot %in% B$uniprot),]
  }
  else{
    BC = plyr::match_df(B,C, on = c("uniprot", "Start_Pos", "End_Pos"))
    BD = plyr::match_df(B, D, on = c("uniprot", "Start_Pos", "End_Pos"))
    CD = plyr::match_df(C, D, on = c("uniprot", "Start_Pos", "End_Pos", "domain_protein"))
    ABC = plyr::match_df(AC, B, on = c("uniprot", "Start_Pos", "End_Pos"))
    # ABD = BD[which(BD$uniprot %in% A$uniprot),]
    ABD = plyr::match_df(BD, A, on = c("uniprot", "Start_Pos", "End_Pos"))
    ACD = plyr::match_df(AC, D, on = c("uniprot", "Start_Pos", "End_Pos", "domain_protein"))
    BCD = plyr::match_df(CD, B, on = c("uniprot", "Start_Pos", "End_Pos"))
    ABCD = plyr::match_df(ACD, B, on = c("uniprot", "Start_Pos", "End_Pos"))
  }

  metadata = c("A" = "Enriched_Domain", "B" = "ClinVar_Path", "C" = "PepSite", "D" = "iELM_HMMs",
               "AB" = "Enriched_Domain&ClinVar_Path", "AC" = "Enriched_Domain&PepSite", "AD" = "Enriched_Domain&iELM_HMMs",
               "BC" = "ClinVar_Path&PepSite", "BD" = "ClinVar_Path&iELM_HMMs", "CD" = "PepSite&iELM_HMMs",
               "ABC" = "Enriched_Domain&ClinVar_Path&PepSite", "ABD" = "Enriched_Domain&ClinVar_Path&iELM_HMMs",
               "ACD" = "Enriched_Domain&PepSite&iELM_HMMs", "BCD" = "ClinVar_Path&PepSite&iELM_HMMs",
               "ABCD" = "Enriched_Domain&ClinVar_Path&PepSite&iELM_HMMs")
  metadata = as.data.frame(metadata)
  metadata$sym = rownames(metadata)
  colnames(metadata) = c("Term", "sym")
  rownames(metadata) = NULL

  final_data = list("A" = A, "B" = B, "C" = C, "D" = D, "AB" = AB, "AC" = AC, "AD" = AD, "BC" = BC, "BD" = BD,
                    "CD" = CD, "ABC" = ABC, "ABD" = ABD, "ACD" = ACD, "BCD" = BCD, "ABCD" = ABCD, "all_input_prots" = qslim_input, "all_pred_prots" =  all_pred_prots, "Metadata" = metadata)
  return(final_data)
}
all_levels_prot = function(Term,required, bench_data = NULL, bg = c("qslim_input", "qslim_pred"), consider_possible_only = T){
  if (is.null(bench_data)){
    bench_data = Get_ELM_all(ELM_data_type = "instances")
    bench_data = bench_data[stringr::str_detect(bench_data$Entry_name, "HUMAN"),]
    bench_data = bench_data[which(bench_data$Logic == "true positive"),]
    #bench_data = bench_data[which(bench_data$Class_type %nin% c("CLV", "TRG")),]
  }
  all_input_prots = unique(unlist(required$all_input_prots))
  if (consider_possible_only){
    bench_data = bench_data[which(bench_data$uniprot %in% all_input_prots),]

    #PRMdb_tp = required$PRMdb$PRM_possible
  }
  Total = required$all_pred_prots
  if (bg == "qslim_input"){
    if ("uniprot" %in% names(required[[Term]])){
      prots = unique(required[[Term]]$uniprot)
    }
    else{
      prots = unique(required[[Term]])
    }
    bench_prots = unique(bench_data$uniprot)
    #prmdb_prots = unique(PRMdb_tp$motif_protein)

    TP = length(intersect(prots, bench_prots))
    FP = length(prots[which(prots %nin% bench_prots)])
    FN = length(bench_prots[which(bench_prots %nin% prots)])
    TN = length(all_input_prots[which(all_input_prots %nin% c(prots, bench_prots))])

    metrics = c(TP, FN, FP, TN)
    result = fisher_and_F1(metrics, pval_method = "greater")
    #result$Reference = "ELM"

    # TP_prm = length(intersect(prots, prmdb_prots))
    # FP_prm = length(prots[which(prots %nin% prmdb_prots)])
    # FN_prm = length(prmdb_prots[which(prmdb_prots %nin% prots)])
    # TN_prm = length(all_input_prots[which(all_input_prots %nin% c(prots, prmdb_prots))])
    #
    # metrics_prm = c(TP_prm, FN_prm, FP_prm, TN_prm)
    # result_prm = fisher_and_F1(metrics_prm, pval_method = "greater")
    # result_prm$Reference = "PRMdb"

    #final = rbind.data.frame(result, result_prm)
    return(result)
  }
  if (bg == "qslim_pred"){
    if ("uniprot" %in% names(required[[Term]])){
      prots = unique(required[[Term]]$uniprot)
    }
    else{
      prots = unique(required[[Term]])
    }
    #ELM_prots = unique(ELM_tp[which(ELM_tp$uniprot %in% all_input_prots),]$uniprot)
    bench_prots = unique(bench_data$uniprot)
    #prmdb_prots = unique(PRMdb_tp$motif_protein)

    TP = length(intersect(prots, bench_prots))
    FP = length(prots[which(prots %nin% bench_prots)])
    FN = length(bench_prots[which(bench_prots %nin% prots)])
    TN = length(Total[which(Total %nin% c(prots, bench_prots))])

    metrics = c(TP, FN, FP, TN)
    result = fisher_and_F1(metrics, pval_method = "greater")
    #result$Reference = "ELM"

    # TP_prm = length(intersect(prots, prmdb_prots))
    # FP_prm = length(prots[which(prots %nin% prmdb_prots)])
    # FN_prm = length(prmdb_prots[which(prmdb_prots %nin% prots)])
    # TN_prm = length(Total[which(Total %nin% c(prots, prmdb_prots))])
    #
    # metrics_prm = c(TP_prm, FN_prm, FP_prm, TN_prm)
    # result_prm = fisher_and_F1(metrics_prm, pval_method = "greater")
    # result_prm$Reference = "PRMdb"
    #
    # final = rbind.data.frame(result, result_prm)
    return(result)
  }
}

#' Motif instances enrichment all filters
#'
#' @author Bishoy Wadie, Evangelia Petsalaki
#' @return
#' A dataframe with odds ratio of each filter.
#' @export
#' @importFrom utils read.delim read.table read.csv
All_motif_inst_OR = function(EM_list, Clinvar_path, iELM_insts, new_hits_ints, pepsite_res, prot_seqs, bench_data = NULL, min_dom = 5, adj_pval = 0.05, datasets_input = NULL, prot_only = F){

  All_Terms = prepare_required_data_OR(EM = EM_list, clin_path = Clinvar_path, iELM = iELM_insts, new_hits_ints = new_hits_ints, min_domain = min_dom,
                                       adj.pval = adj_pval, slim_input_prots = datasets_input, pepsite_res = pepsite_res, prot_only = prot_only)
  all_hits = new_hits_ints
  if (is.null(datasets_input)){
    # TODO change input path
    qslim_input = download_data(path = "Evaluation_data/proteins_per_dataset_allhitsint.RDS", is.rds = T)
  }
  else{
    qslim_input = datasets_input
  }

  if (is.null(bench_data)){
    bench_data = Get_ELM_all(ELM_data_type = "instances")
    ELM_all_instances = bench_data[stringr::str_detect(bench_data$Entry_name, "HUMAN"),]
    ELM_all_instances = ELM_all_instances[which(ELM_all_instances$Logic == "true positive"),]
    #ELM_all_instances = ELM_all_instances[which(ELM_all_instances$Class_type %nin% c("CLV", "TRG")),]
  }
  else{
    ELM_all_instances = bench_data
    # colnames(ELM_all_instances)[c(2,3,4)] = c("uniprot", "start", "end")
    # ELM_all_instances$Logic = "true positive"
  }

  # if (bench_data_type == "ELM"){
  #   ELM_all_instances = All_Terms$ELM$ELM_Human
  # }
  # else{
  #   ELM_all_instances = All_Terms$PRMdb$PRM_Human
  #   colnames(ELM_all_instances)[c(2,3,4)] = c("uniprot", "start", "end")
  #   ELM_all_instances$Logic = "true positive"
  # }
  if (is.null(prot_seqs)){
    # TODO change input path
    uniprot_seqs = read.csv(download_data(path = "Evaluation_data/uniprot_sequences.csv.gz"), header = T)
  }
  else{
    uniprot_seqs = prot_seqs
  }
  # uniprot_seqs[which(stringr::str_detect(uniprot_seqs$uniprot, "-") == TRUE),]$uniprot =
  #   substr(uniprot_seqs[which(stringr::str_detect(uniprot_seqs$uniprot, "-") == TRUE),]$uniprot, 1, nchar(uniprot_seqs[which(stringr::str_detect(uniprot_seqs$uniprot, "-") == TRUE),]$uniprot) -2 )

  required = list("pred_hits" = all_hits, "qslim_input" = qslim_input, "ELM_all" = ELM_all_instances,
                  "uniprot_seqs"  = uniprot_seqs)
  all_input_prots = unique(unlist(qslim_input))

  ELM_possible = ELM_all_instances[which(ELM_all_instances$uniprot %in% all_input_prots),]
  ELM_possible = ELM_possible[which(ELM_possible$Logic == "true positive"),]
  qslim_instances = dplyr::select(all_hits, uniprot, Start_Pos, End_Pos)
  qslim_instances = qslim_instances[!duplicated(qslim_instances),]

  possible_motifs = unique(ELM_possible$Id)

  all_metrics_per_motif = list()
  weights_per_motif = list()
  for (i in 1:length(possible_motifs)){
    result = calc_metrics_per_motif(df = qslim_instances, motif = possible_motifs[i], min_annot = 3, required_data = required, match_insts = T)
    all_metrics_per_motif[[i]] = result[["metrics"]]
    weights_per_motif[[i]] = result[["weights"]]
  }
  all_metrics_per_motif = dplyr::bind_rows(all_metrics_per_motif)
  all_metrics_per_motif$OR = ((all_metrics_per_motif$Rc_res * -1) / (all_metrics_per_motif$Rc_res -1)) /
    ((all_metrics_per_motif$SP_res -1) / (all_metrics_per_motif$SP_res * -1))
  weights_per_motif = dplyr::bind_rows(weights_per_motif)
  qslim_F1 = calc_F1_score_HH(weights = weights_per_motif, metrics = all_metrics_per_motif)
  qslim_F1$OR = 0
  for (i in 1:nrow(qslim_F1)){
    if (qslim_F1$site_resid[i] != "residue_wise"){
      qslim_F1$OR[i] = NA
    }
    else{
      qslim_F1$OR[i] = ((qslim_F1$Rc[i] * -1) / (qslim_F1$Rc[i] -1)) /
        ((qslim_F1$SP[i] -1) / (qslim_F1$SP[i] * -1))
    }
  }
  qslim_F1$title = "qslim_predicted"
  qslim_result = list("All_F1_scores" = qslim_F1, "all_metrics_per_motif" = all_metrics_per_motif,
                      "weights_per_motif" = weights_per_motif)

  metadata = All_Terms[["Metadata"]]

  All_F1_scores = list()
  All_F1_scores[["qslim_predicted"]] = qslim_result
  for (i in 1:nrow(metadata)){
    if (nrow(All_Terms[[metadata$sym[i]]]) == 0){
      next()
    }
    filtered = all_levels_motif_inst(All_Terms = All_Terms , metadata$sym[i], prereqs = required, qslim_predicted_bg = qslim_result$all_metrics_per_motif)
    All_F1_scores[[metadata$Term[which(metadata$sym == metadata$sym[i])]]] = filtered
  }
  return(All_F1_scores)
}

all_levels_motif_inst = function(All_Terms,sym, prereqs, qslim_predicted_bg){
  metadata = All_Terms[["Metadata"]]
  ELM_all_instances = prereqs[["ELM_all"]]
  qslim_input = prereqs[["qslim_input"]]
  all_input_prots = unique(unlist(qslim_input))

  possible_motifs = unique(qslim_predicted_bg$Id)

  all_metrics_per_motif = list()
  weights_per_motif = list()
  for (i in 1:length(possible_motifs)){
    result = calc_metrics_per_motif(All_Terms[[sym]], possible_motifs[i], min_annot = 3, required_data = prereqs, match_insts = T)
    all_metrics_per_motif[[i]] = result[["metrics"]]
    weights_per_motif[[i]] = result[["weights"]]
  }
  all_metrics_per_motif = dplyr::bind_rows(all_metrics_per_motif)
  all_metrics_per_motif$OR = ((all_metrics_per_motif$Rc_res * -1) / (all_metrics_per_motif$Rc_res -1)) /
    ((all_metrics_per_motif$SP_res -1) / (all_metrics_per_motif$SP_res * -1))
  weights_per_motif = dplyr::bind_rows(weights_per_motif)

  final_result = calc_F1_score_HH(weights = weights_per_motif, metrics = all_metrics_per_motif)

  final_result$OR = 0
  for (i in 1:nrow(final_result)){
    if (final_result$site_resid[i] != "residue_wise"){
      final_result$OR[i] = NA
    }
    else{
      final_result$OR[i] = ((final_result$Rc[i] * -1) / (final_result$Rc[i] -1)) /
        ((final_result$SP[i] -1) / (final_result$SP[i] * -1))
    }
  }

  final_result$title = metadata$Term[which(metadata$sym == sym)]
  return(list("All_F1_scores" = final_result, "all_metrics_per_motif" = all_metrics_per_motif,
              "weights_per_motif" = weights_per_motif))
}

#' @importFrom utils download.file
pepsite_clinvar_OR = function(slim_hits, clinvar_path){

  all_inst = dplyr::select(slim_hits, uniprot, Start_Pos, End_Pos)
  all_inst = all_inst[!duplicated(all_inst),]
  all_inst$inst_id = 1:nrow(all_inst)

  clinvar_path = clinvar_path %>% dplyr::select(uniprot, Start_Pos, End_Pos)
  clinvar_path = clinvar_path[!duplicated(clinvar_path),]
  clinvar_path = dplyr::left_join(clinvar_path, all_inst)

  # TODO change input path to Biostudies

  final_motif_domain_ints =  readRDS("data/output/HV_pepsite_final_with_AF.rds")

  pepsite_inst = final_motif_domain_ints %>% dplyr::select(motif_protein, start, end)
  pepsite_inst = pepsite_inst[!duplicated(pepsite_inst),]
  colnames(pepsite_inst) = colnames(clinvar_path)[c(1:3)]
  pepsite_inst = dplyr::left_join(pepsite_inst, all_inst)

  TP = length(intersect(pepsite_inst$inst_id, clinvar_path$inst_id))
  FP = length(unique(pepsite_inst$inst_id[which(pepsite_inst$inst_id %nin% clinvar_path$inst_id)]))
  FN = length(unique(clinvar_path$inst_id[which(clinvar_path$inst_id %nin% pepsite_inst$inst_id)]))
  TN = length(unique(all_inst$inst_id[which(all_inst$inst_id %nin% unique(c(clinvar_path$inst_id, pepsite_inst$inst_id)))]))

  final = c(TP, FN, FP, TN)
  return(fisher_and_F1(final, pval_method = "greater"))
}

#' @importFrom utils download.file
Novel_pepsite_clinvar_Chisq = function(pep_clin_EM,main_output = NULL, occ_out = NULL, compari_out = NULL, parse_occ = F, preprocess_compari = T, clinvar_path, compari_score_cutoff = 2){
  if (parse_occ & is.null(occ_out)){
    all_hits = parse_qslim(main_output = main_output, occ_out = occ_out, compari_output = compari_out, preprocess_compari = preprocess_compari)
  }
  else{
    all_hits = occ_out
  }

  all_inst = dplyr::select(all_hits, uniprot, Start_Pos, End_Pos)
  all_inst = all_inst[!duplicated(all_inst),]
  all_inst$inst_id = 1:nrow(all_inst)

  # clinvar_path = clinvar_path %>% dplyr::select(uniprot, Start_Pos, End_Pos)
  # clinvar_path = clinvar_path[!duplicated(clinvar_path),]
  # clinvar_path = dplyr::left_join(clinvar_path, all_inst)
  #
  # final_motif_domain_ints =  pepsit_hits #Pepsite and Enriched_Domains
  #
  # pepsite_inst = final_motif_domain_ints %>% dplyr::select(motif_protein, start, end)
  # pepsite_inst = pepsite_inst[!duplicated(pepsite_inst),]
  # colnames(pepsite_inst) = colnames(clinvar_path)[c(1:3)]
  # pepsite_inst = dplyr::left_join(pepsite_inst, all_inst)

  pepsite_inst = dplyr::left_join(pep_clin_EM, all_inst)

  #pepsite_clinvar_insts = intersect(pepsite_inst$inst_id, clinvar_path$inst_id)
  pepsite_clinvar_insts = unique(pepsite_inst$inst_id)


  x = all_hits %>% dplyr::select(uniprot, Start_Pos, End_Pos, Id, Score)
  x = x[!duplicated(x),]
  x = dplyr::left_join(x, all_inst)

  x$high_score = ifelse(x$Score > compari_score_cutoff, "Yes", "No")
  y = x %>% dplyr::group_by(inst_id) %>% dplyr::summarise(n = dplyr::n())
  colnames(y)[2] = "Total"
  z = x %>% dplyr::group_by(inst_id, high_score) %>% dplyr::summarise(n = dplyr::n())
  z = dplyr::left_join(z, y)
  z = z[which(z$high_score == "No"),]
  z$No_pct = z$n / z$Total
  z = z[which(z$No_pct == 1),]
  Novel_inst_ids = unique(z$inst_id)

  all_inst_ids = unique(all_inst$inst_id)
  Novel_pepclin = length(intersect(Novel_inst_ids, pepsite_clinvar_insts))
  Novel_not_pepclin = length(unique(Novel_inst_ids[which(Novel_inst_ids %nin% pepsite_clinvar_insts)]))
  pepclin_not_novel = length(unique(pepsite_clinvar_insts[which(pepsite_clinvar_insts %nin% Novel_inst_ids)]))
  not_novel_notpepclin = length(unique(all_inst_ids[which(all_inst_ids %nin% c(Novel_inst_ids, pepsite_clinvar_insts))]))

  chisq_matrix = matrix(c(Novel_pepclin, Novel_not_pepclin, pepclin_not_novel, not_novel_notpepclin), nrow = 2, ncol = 2, byrow = T)
  return(chisq.test(chisq_matrix))
}
