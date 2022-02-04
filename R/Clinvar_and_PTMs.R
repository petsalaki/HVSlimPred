#' Map predicted motifs to PTMs
#'
#' Maps the predicted motifs to PTMs retrieved from uniprot.
#'
#' @author Bishoy Wadie, Evangelia Petsalaki
#' @return
#' A data frame containing PTMs mapped to all predicted motifs.
#' @export
#' @importFrom utils read.delim read.table read.csv
#' @examples
#' PTM_mapping = Map_PTM()
Map_PTM = function(slim_hits, all_pfam_dom = NULL, uniprot_mod_res = NULL){
  new_hits_int = slim_hits
  # if (parse_occ & is.null(occ_out)){
  #   new_hits_int = parse_qslim(main_output = main_output ,occ_out = occ_out, compari_output = compari_out, preprocess_compari = preprocess_compari)
  # }
  # else{
  #   new_hits_int = occ_out
  # }
  # TODO Change this function after Get_ELM_all is changed
  ELM_all_instances = Get_ELM_all(ELM_data_type = "instances")
  if (is.null(uniprot_mod_res)){
    # TODO change input path
    uniprot_PTM = read.table(download_data(path = "PTM/UP000005640_9606_mod_res.bed", is.rds = F), sep = "\t", quote = "")
  }
  else{
    uniprot_PTM = uniprot_mod_res
  }

  Mapped_qslim = new_hits_int %>% dplyr::select(-Accession, -start, -end, -Logic)
  Mapped_qslim = Mapped_qslim[!duplicated(Mapped_qslim),]


  # FIXME change index of ELM_all_instances to dplyr select
  Mapped_PL = dplyr::left_join(Mapped_qslim, ELM_all_instances[,c(5,7,8)], by = c("uniprot"))
  Mapped_PL = qslim_ELM_motif_distance(Mapped_PL)
  Mapped_PL = Mapped_PL[which(Mapped_PL$motif_distance == 0),]

  Mapped_qslim = new_hits_int

  # FIXME change index of uniprot_PTM to dplyr:select
  uniprot_PTM = uniprot_PTM[,c(4,14)]
  uniprot_PTM = uniprot_PTM[!duplicated(uniprot_PTM),]
  colnames(uniprot_PTM) = c("uniprot", "Notes")
  uniprot_PTM = uniprot_PTM[which(uniprot_PTM$uniprot %in% Mapped_qslim$uniprot),]
  uniprot_PTM = uniprot_PTM %>% tidyr::separate(Notes, into = c("PTM_res", "PTM_type", "PTM_comments"), sep = ";")
  uniprot_PTM$res_num = substr(uniprot_PTM$PTM_res, 2, nchar(uniprot_PTM$PTM_res))

  qslim_insts = new_hits_int %>% dplyr::select(uniprot, Start_Pos, End_Pos)
  qslim_insts = qslim_insts[which(qslim_insts$uniprot %in% uniprot_PTM$uniprot),]
  qslim_insts = qslim_insts[!duplicated(qslim_insts),]

  mapped_PTM_qslim = list()
  for (i in 1:nrow(qslim_insts)){
    mapped_PTM_qslim[[i]] = map_uniprot_PTMs(PTMs = uniprot_PTM, qslim_insts$uniprot[i], qslim_insts$Start_Pos[i],
                                             qslim_insts$End_Pos[i])
  }
  mapped_PTM_qslim = dplyr::bind_rows(mapped_PTM_qslim)
  mapped_PTM_qslim = mapped_PTM_qslim[!duplicated(mapped_PTM_qslim),]

  Mapped_qslim_uniprot = dplyr::left_join(Mapped_qslim, mapped_PTM_qslim)
  Mapped_qslim_uniprot = Label_ELM_instance(ELM = ELM_all_instances,df = Mapped_qslim_uniprot, PL = Mapped_PL)

  pfam_filters_uniprot = Motifcarr_to_pfam_dom_enrich(mot_car_prots = unique(Mapped_qslim_uniprot$uniprot), all_pfam_dom_list = all_pfam_dom)
  Mapped_qslim_uniprot = dplyr::left_join(Mapped_qslim_uniprot, pfam_filters_uniprot)
  return(Mapped_qslim_uniprot)

}

map_uniprot_PTMs = function(PTMs,uniprot_id, start, end){
  if (uniprot_id %nin% PTMs$uniprot){
    return(NULL)
  }
  final_result = data.frame(uniprot = uniprot_id, Start_Pos = start, End_Pos = end)
  filtered_PTMs = PTMs %>% dplyr::filter(uniprot == uniprot_id)
  matched = list()
  for (i in 1:nrow(filtered_PTMs)){
    if(filtered_PTMs$res_num[i] >= start && filtered_PTMs$res_num[i] <= end){
      matched[[i]] = filtered_PTMs[i,]
    }
    else{
      matched[[i]] = NULL
    }
  }
  matched = dplyr::bind_rows(matched)
  if (nrow(matched) != 0){
    final_result = dplyr::left_join(final_result, matched[,c(1:4)], by = "uniprot")
    return(final_result)
  }
  else{
    return(NULL)
  }
}
Label_ELM_instance = function(ELM, PL, df){
  df$Right_P_ELM = NA
  df$entry_id = 1:nrow(df)
  df = df[,c(length(df),1:(length(df) -1))]
  df$Right_M_ELM = NA
  df$Right_L_ELM = NA


  df[which(df$uniprot %in% ELM$uniprot),]$Right_P_ELM = "Yes"
  df[which(df$uniprot %nin% ELM$uniprot),]$Right_P_ELM = "No"


  right_mot = df[!is.na(df$Accession),]
  df[which(df$entry_id %in% right_mot$entry_id),]$Right_M_ELM = "Yes"
  df[which(df$Right_P_ELM == "Yes" & df$entry_id %nin% right_mot$entry_id),]$Right_M_ELM = "No"

  right_loc = plyr::match_df(df, PL, on = c("Dataset","uniprot", "Start_Pos", "End_Pos", "Id"))

  df[which(df$entry_id %in% right_loc$entry_id),]$Right_L_ELM = "Yes"
  df[which(df$Right_P_ELM == "Yes" & df$entry_id %nin% right_loc$entry_id),]$Right_L_ELM = "No"
  df[which(df$Right_M_ELM == "Yes" & df$motif_distance != 0),]$Right_L_ELM = "No"

  return(df[,-1])

}

#' Map predicted motifs to ClinVar variants
#'
#' Maps the predicted motifs to disease variants retrieved from ClinVar.
#'
#' @author Bishoy Wadie, Evangelia Petsalaki
#' @return
#' A data frame containing ClinVar variants mapped to all predicted motifs.
#' @export
#' @importFrom utils read.delim read.table read.csv
#' @examples
#' PTM_mapping = Map_PTM()
Map_ClinVar = function(slim_hits, Clinvar_variant_summary=NULL, pathogenic_only = T, qslim_noAA_gncoords = NULL){
  new_hits_int = slim_hits
  qslim_human_prots = new_hits_int[which(new_hits_int$Org == "HUMAN"),]$uniprot
  qslim_human_prots = qslim_human_prots[!duplicated(qslim_human_prots)]
  uniprot_to_entrez = AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = qslim_human_prots, columns = "ENTREZID", keytype = "UNIPROT")

  if (is.null(Clinvar_variant_summary)){
    # FIXME change reading from download data after download_data has been fixed
    clinvar = read.delim(download_data(path = "ClinVar/variant_summary_2021-03.txt.gz", is.rds = F), header = T, sep = "\t")
    message("ClinVar Variant_summary was downloaded on 04-03-2021\n", "Please run use_latest_clinvar() to download the latest version\n")
  }
  else{
    clinvar = Clinvar_variant_summary
    rm(Clinvar_variant_summary)
  }

  clinvar = clinvar[which(clinvar$GeneID %in% uniprot_to_entrez$ENTREZID),]
  clinvar = clinvar[,c(1:length(clinvar_colnames))]
  names(clinvar) = clinvar_colnames
  uniprot_entrez_clinvar = dplyr::filter(uniprot_to_entrez, ENTREZID %in% clinvar$GeneID)

  clinvar_AA = clinvar[stringr::str_detect(clinvar$Name, "p.", negate = F),]
  clinvar_AA$GeneID = as.character(clinvar_AA$GeneID)
  clinvar_AA$AA_change = sub(".*p[.]", "", clinvar_AA$Name)
  clinvar_AA = dplyr::filter(clinvar_AA, Type %in% c("Deletion", "single nucleotide variant", "deletion"))
  clinvar_AA = clinvar_AA %>% tidyr::separate(AA_change, into = c("fromAA", "toAA"), sep = "_")
  clinvar_AA[is.na(clinvar_AA$toAA),]$toAA = clinvar_AA[is.na(clinvar_AA$toAA),]$fromAA
  clinvar_AA$fromAA = gsub(".*?([0-9]+).*", '\\1', clinvar_AA$fromAA)
  clinvar_AA$toAA = gsub(".*?([0-9]+).*", '\\1', clinvar_AA$toAA)
  clinvar_AA$AA_change = sub(".*p[.]", "", clinvar_AA$Name)
  clinvar_AA$fromAA = as.integer(clinvar_AA$fromAA)
  clinvar_AA$toAA = as.integer(clinvar_AA$toAA)


  qslim_insts = new_hits_int %>% dplyr::select(uniprot, Start_Pos, End_Pos)
  qslim_insts = qslim_insts[!duplicated(qslim_insts),]
  qslim_insts = dplyr::left_join(qslim_insts, uniprot_entrez_clinvar, by = c("uniprot" = "UNIPROT"))
  qslim_insts_AA = qslim_insts[which(qslim_insts$ENTREZID %in% clinvar_AA$GeneID),]
  qslim_insts_AA = qslim_insts_AA[!duplicated(qslim_insts_AA),]

  qslim_clinvar_AA = list()
  for (i in 1:nrow(qslim_insts)){
    qslim_clinvar_AA[[i]] = compare_clinvar_qslim(clin_bg = clinvar_AA ,uniprot_id = qslim_insts_AA$uniprot[i], entrez = qslim_insts_AA$ENTREZID[i],
                                                     start_inst = qslim_insts_AA$Start_Pos[i], end_inst = qslim_insts_AA$End_Pos[i], noAA = F)
  }
  qslim_clinvar_AA = dplyr::bind_rows(qslim_clinvar_AA)

  clinvar_no_AA = clinvar[stringr::str_detect(clinvar$Name, "p.", negate = T),]
  clinvar_no_AA$GeneID = as.character(clinvar_no_AA$GeneID)
  clinvar_no_AA = dplyr::filter(clinvar_no_AA, Type %in% c("Deletion", "single nucleotide variant", "deletion"))

  qslim_insts_noAA = qslim_insts[which(qslim_insts$ENTREZID %in% clinvar_no_AA$GeneID),]
  qslim_insts_noAA = qslim_insts_noAA[!duplicated(qslim_insts_noAA),]

  qslim_insts_clinvar_no_AA = qslim_clin_noAA(clin_entrez = uniprot_entrez_clinvar, clin_noAA = clinvar_no_AA, slim_insts = qslim_insts_noAA, ready_reqs = qslim_noAA_gncoords)
  qslim_clinvar_no_AA = list()
  for (i in 1:nrow(qslim_insts_clinvar_no_AA)){
    qslim_clinvar_no_AA[[i]] = compare_clinvar_qslim(clin_bg = clinvar_no_AA, qslim_inst_noAA = qslim_insts_clinvar_no_AA,uniprot_id = qslim_insts_clinvar_no_AA$uniprot[i], entrez = qslim_insts_clinvar_no_AA$ENTREZID[i],
                                                    start_inst = qslim_insts_clinvar_no_AA$DNA_start[i], end_inst = qslim_insts_clinvar_no_AA$DNA_end[i], noAA = T)
  }
  qslim_clinvar_no_AA = dplyr::bind_rows(qslim_clinvar_no_AA)
  qslim_clinvar_no_AA = qslim_clinvar_no_AA[!duplicated(qslim_clinvar_no_AA),]

  # TODO try to standardize the selection of columns to apply for any downloaded ClinVar data

  qslim_clinvar_no_AA = qslim_clinvar_no_AA %>% dplyr::select(-AlleleID, -HGNC_ID, -LastEvaluated,
                                                       -nsv.esv_.dbVar.,-Origin, -Assembly, -ChromosomeAccession,
                                                       -Start, -Stop, -ReferenceAllele, -AlternateAllele,
                                                       -Cytogenetic, -ReviewStatus, -NumberSubmitters, -Guidelines,
                                                       -TestedInGTR, -OtherIDs, -SubmitterCategories, -VariationID)

  qslim_clinvar_AA = qslim_clinvar_AA %>% dplyr::select(-AlleleID, -HGNC_ID, -LastEvaluated,
                                                 -nsv.esv_.dbVar.,-Origin, -Assembly, -ChromosomeAccession,
                                                 -Start, -Stop, -ReferenceAllele, -AlternateAllele,
                                                 -Cytogenetic, -ReviewStatus, -NumberSubmitters, -Guidelines,
                                                 -TestedInGTR, -OtherIDs, -SubmitterCategories, -VariationID)
  Final_clinvar_mapping = plyr::rbind.fill(qslim_clinvar_AA, qslim_clinvar_no_AA)
  Final_clinvar_mapping = Final_clinvar_mapping[!duplicated(Final_clinvar_mapping),]

  new_hits_short = new_hits_int %>% dplyr::select(-Id, -Sim_rank, -Score)
  Final_clinvar_mapping$Start_Pos = as.integer(Final_clinvar_mapping$Start_Pos)
  Final_clinvar_mapping$End_Pos = as.integer(Final_clinvar_mapping$End_Pos)
  Final_clinvar_mapping = dplyr::left_join(new_hits_short, Final_clinvar_mapping, by = c("uniprot", "Start_Pos", "End_Pos"))
  Final_clinvar_mapping = Final_clinvar_mapping[!duplicated(Final_clinvar_mapping),]

  if (pathogenic_only){
    clin_sig = c("Likely pathogenic","Likely pathogenic, drug response",
                 "Likely pathogenic, other","Pathogenic","Pathogenic, Affects","Pathogenic, drug response",
                 "Pathogenic, drug response, other","Pathogenic, other","Pathogenic, protective",
                 "Pathogenic, risk factor","Pathogenic/Likely pathogenic","Pathogenic/Likely pathogenic, drug response",
                 "Pathogenic/Likely pathogenic, other")
    final_clinvar_mapping_path = Final_clinvar_mapping[which(Final_clinvar_mapping$ClinicalSignificance %in% clin_sig),]
    final_clinvar_mapping_path = final_clinvar_mapping_path[!duplicated(final_clinvar_mapping_path),]
    return(list("ClinVar_all" = Final_clinvar_mapping, "ClinVar_path" = final_clinvar_mapping_path))
  }
  else{
    return(Final_clinvar_mapping)
  }
}

#' @importFrom utils download.file
use_latest_clinvar = function(){
  latest_clinvar_name = paste0("variant_summary_", format(Sys.Date(), "%Y"), "-", format(Sys.Date(), "%m"), ".txt.gz")
  if (nchar(as.numeric(format(Sys.Date(), "%m"))) == 1){
    previous_month = paste0("0", as.character((as.numeric(format(Sys.Date(), "%m")) - 1)))
    previous_clinvar = paste0("variant_summary_", format(Sys.Date(), "%Y"), "-", previous_month ,".txt.gz")
  }
  else{
    previous_month = as.character((as.numeric(format(Sys.Date(), "%m")) - 1))
    previous_clinvar = paste0("variant_summary_", format(Sys.Date(), "%Y"), "-", previous_month ,".txt.gz")
  }

  ftp_base <- "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/"
  list_files <- curl::new_handle()
  curl::handle_setopt(list_files, ftp_use_epsv = TRUE, dirlistonly = FALSE)

  con <- curl::curl(url = ftp_base, "r", handle = list_files)
  files <- readLines(con)
  if (any(stringr::str_detect(files, latest_clinvar_name))){
    tmp = tempfile()
    clinvar = download.file(paste0(ftp_base, latest_clinvar_name), destfile = tmp)
    close(con)
    return(tmp)
  }
  else{
    tmp = tempfile()
    clinvar = download.file(paste0(ftp_base, previous_clinvar), destfile = tmp)
    close(con)
    return(tmp)
  }
}

uniprot_gncoords = function(clin_entrez, clin_noAA,slim_insts){
  all_uniprot_gncoords = list()
  no_AA_prots = clin_entrez[which(clin_entrez$ENTREZID %in% clin_noAA$GeneID),]$UNIPROT
  slim_insts = slim_insts[which(slim_insts$uniprot %in% no_AA_prots),]
  for (i in 1:nrow(slim_insts)){
    requestURL <- paste0("https://www.ebi.ac.uk/proteins/api/coordinates/location/", slim_insts$uniprot[i],
                         ":", slim_insts$Start_Pos[i], "-", slim_insts$End_Pos[i])
    r <- httr::GET(requestURL, httr::accept("application/json"))
    if (r$status_code != 200){
      all_uniprot_gncoords[[i]] = NULL
      next()
    }
    json <- jsonlite::toJSON(httr::content(r))
    json = jsonlite::fromJSON(json)
    all_uniprot_gncoords[[i]] = json[["locations"]]
  }
  all_uniprot_gncoords = dplyr::bind_rows(all_uniprot_gncoords)
  all_uniprot_gncoords = all_uniprot_gncoords[,-4]
  all_uniprot_gncoords = all_uniprot_gncoords[!duplicated(all_uniprot_gncoords),]
  all_uniprot_gncoords = dplyr::mutate_all(all_uniprot_gncoords, as.character)
  all_uniprot_gncoords = all_uniprot_gncoords[which(all_uniprot_gncoords$chromosome %in%
                                                      c(as.character(c(1:22)), "X", "Y")),]
  all_uniprot_gncoords = dplyr::left_join(all_uniprot_gncoords, clin_entrez, by = c("accession" = "UNIPROT"))
  all_uniprot_gncoords = all_uniprot_gncoords[,which(colnames(all_uniprot_gncoords) %in% c("accession","taxid","chromosome", "proteinStart",
                                                                                           "geneStart", "proteinEnd", "geneEnd", "ENTREZID"))]
  colnames(all_uniprot_gncoords) = c("uniprot", "taxid", "chromosome", "Start_Pos", "DNA_start", "End_Pos", "DNA_end", "ENTREZID")
  return(all_uniprot_gncoords)
}
compare_clinvar_qslim = function(clin_bg, qslim_inst_noAA,uniprot_id,entrez,start_inst, end_inst, noAA = F){
  if (!noAA & missing(qslim_inst_noAA)){
    clinvar_filtered = dplyr::filter(clin_bg, GeneID == entrez)
    clinvar_filtered = clinvar_filtered[which(clinvar_filtered$fromAA >= start_inst & clinvar_filtered$toAA <= end_inst),]
    if(nrow(clinvar_filtered) != 0){
      final = data.frame(uniprot = uniprot_id, ENTREZID = entrez, Start_Pos = start_inst, End_Pos = end_inst)
      final = dplyr::left_join(final, clinvar_filtered, by = c("ENTREZID" = "GeneID"))
      return(final)
    }
    else{
      return(NULL)
    }
  }
  else{
    clinvar_filtered = dplyr::filter(clin_bg, GeneID == entrez)
    clinvar_filtered = clinvar_filtered[which(clinvar_filtered$Start >= start_inst & clinvar_filtered$Stop <= end_inst),]
    if(nrow(clinvar_filtered) != 0){
      final = data.frame(uniprot = uniprot_id, ENTREZID = entrez, DNA_start = start_inst, DNA_end = end_inst)
      final = plyr::match_df(qslim_inst_noAA, final)
      final = final %>% dplyr::select(uniprot, Start_Pos, End_Pos, ENTREZID)
      final = dplyr::left_join(final, clinvar_filtered, by = c("ENTREZID" = "GeneID"))
      return(final)
    }
    else{
      return(NULL)
    }
  }
}
qslim_clin_noAA = function(clin_entrez, clin_noAA,slim_insts, ready_reqs){
  if (is.null(ready_reqs)){
    # TODO change input path
    qslim_insts_clinvar_no_AA = download_data(path = "ClinVar/qslim_insts_clinvar_no_AA.rds", is.rds = T)
  }
  else{
    qslim_insts_clinvar_no_AA = ready_reqs
  }
  qslim_insts_clinvar_no_AA = dplyr::left_join(qslim_insts_clinvar_no_AA, clin_entrez, by = c("accession" = "UNIPROT"))
  colnames(qslim_insts_clinvar_no_AA) = c("uniprot", "taxid", "chromosome", "Start_Pos", "DNA_start", "End_Pos", "DNA_end", "ENTREZID")
  send_new = mismatch_df(slim_insts, qslim_insts_clinvar_no_AA)
  if (nrow(send_new) == 0){
    return(qslim_insts_clinvar_no_AA)
  }
  else{
    new_insts_clinvar_no_AA = uniprot_gncoords(clin_entrez = clin_entrez, clin_noAA = clin_noAA, slim_insts = send_new)
    final = rbind.data.frame(qslim_insts_clinvar_no_AA, new_insts_clinvar_no_AA)
    return(final)
  }
}
