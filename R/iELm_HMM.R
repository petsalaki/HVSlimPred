#' iELM HMM analysis on H1 proteins
#'
#' Preprocess and evaluate identified motif-binding domains in H1 proteins using both iELM and PFAM HMMs.
#'
#' @return
#' Combined output from HMMER against both iELM and pfam HMMs for all H1 proteins.
#' @export
#' @author Bishoy Wadie, Evangelia Petsalaki
#' @importFrom utils read.delim read.table read.csv
#' @references
#' Weatheritt, Robert J., et al. "iELM—a web server to explore short linear motif-mediated interactions." \emph{Nucleic acids research} 40.W1 (2012): W364-W369.
#' @examples
#' HMM_results = iELM_HMM()
iELM_HMM = function(path_HMM_pfam_domtblouts = NULL, path_HMM_iELM_domtblouts = NULL, pct_coverage = 0.8, ELM_renamed_classes = NULL){
  HMM_final_results = read_iELM_HMMER_domtblouts(pfam_domtblouts = path_HMM_pfam_domtblouts, iELM_domtblouts = path_HMM_iELM_domtblouts)
  HMM_final_results = HMM_final_results %>% tidyr::separate(domain_name, into = c("remove_col", "H1", "H1_sym"), sep = "[|]") %>% dplyr::select(-remove_col)
  HMM_final_results = HMM_final_results[which(HMM_final_results$sequence_evalue < 0.01),]
  HMM_final_results$pct_coverage = ((HMM_final_results$hmm_to - HMM_final_results$hmm_from) + 1) / HMM_final_results$qlen
  HMM_final_results = HMM_final_results[which(HMM_final_results$pct_coverage >= pct_coverage),]

  HMM_pfam_names = sub("\\.*.txt","", pfam_HMMs)
  HMM_pfam_names = data.frame(HMM_ELM = HMM_pfam_names, HMM_ELM_type = "pfam")
  HMM_domain_all_names = sub("\\.*.txt","", iELM_HMMs)
  HMM_domain_all_names = data.frame(HMM_ELM = HMM_domain_all_names, HMM_ELM_type = "iELM")
  HMM_types_names = rbind.data.frame(HMM_pfam_names, HMM_domain_all_names)
  HMM_final_results = dplyr::left_join(HMM_final_results, HMM_types_names)

  if(is.null(ELM_renamed_classes)){
    # TODO change input path
    ELM_renamed = read.delim(download_data(path = "iELM_HMM/ELM_renamed.tsv", is.rds = F))
  }
  else{
    ELM_renamed = ELM_renamed_classes
  }
  HMM_final_results = as.data.frame(HMM_final_results)
  HMM_final_results[which(HMM_final_results$HMM_ELM_type == "pfam"),25] = stringr::str_replace_all(HMM_final_results[which(HMM_final_results$HMM_ELM_type == "pfam"),25], "[+]", "-")

  HMM_final_results = HMM_final_results %>% tidyr::separate(HMM_ELM, into = c("ELM_Class", "known_pdb"), sep = "[+]", extra = "merge")

  for (i in 1:nrow(ELM_renamed)){
    HMM_final_results$ELM_Class = stringr::str_replace_all(HMM_final_results$ELM_Class, ELM_renamed$Renamed_From[i], ELM_renamed$Renamed_To[i])
  }
  #the duplication of the for loop below is on purpose
  for (i in 1:nrow(ELM_renamed)){
    HMM_final_results$ELM_Class = stringr::str_replace_all(HMM_final_results$ELM_Class, ELM_renamed$Renamed_From[i], ELM_renamed$Renamed_To[i])
  }
  return(HMM_final_results)
}

read_iELM_HMMER_domtblouts = function(pfam_domtblouts = NULL, iELM_domtblouts = NULL){
  if (is.null(pfam_domtblouts) & is.null(iELM_domtblouts)){
    x = list()
    y = list()
    for (i in 1:length(pfam_HMMs)){
      # TODO change input path
      x[[i]] = parse_Hmmer_output(pfam_HMMs[i], download_data(path = file.path("iELM_HMM", "pfam_H1_domtblouts", pfam_HMMs[i]), is.rds = F))
    }
    x = dplyr::bind_rows(x)

    for (i in 1:length(iELM_HMMs)){
      # TODO change input path
      y[[i]] = parse_Hmmer_output(iELM_HMMs[i], download_data(path = file.path("iELM_HMM", "domain_all_H1_domtblouts", iELM_HMMs[i]), is.rds = F))
    }
    y = dplyr::bind_rows(y)
    HMM_final_results = rbind.data.frame(x,y)
  }
  else{
    HMM_pfam_files = list.files(pfam_domtblouts)
    HMM_domain_all_files = list.files(iELM_domtblouts)

    x = list()
    y = list()
    for (i in 1:length(HMM_pfam_files)){
      x[[i]] = parse_Hmmer_output(HMM_pfam_files[i], file.path(pfam_domtblouts, HMM_pfam_files[i]))
    }
    x = dplyr::bind_rows(x)

    for (i in 1:length(HMM_domain_all_files)){
      y[[i]] = parse_Hmmer_output(HMM_domain_all_files[i], file.path(iELM_domtblouts, HMM_domain_all_files[i]))
    }
    y = dplyr::bind_rows(y)

    HMM_final_results = rbind.data.frame(x,y)
  }
  return(HMM_final_results)
}

# TODO Add a function or some note to prepare input and run HMMER

parse_Hmmer_output = function(file_name, file_path){
  ELM_name = strsplit(file_name, "[.]")[[1]][1]
  output = rhmmer::read_domtblout(file_path)
  if (nrow(output) == 0){
    return(NULL)
  }
  else{
    output$HMM_ELM = ELM_name
    return(output)
  }
}
