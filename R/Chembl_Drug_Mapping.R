#' Maps Chembl Drugs to H1 proteins
#'
#' Maps ChEMBL small molecule compounds to H1 proteins and connects the drug indication with the ClinVar disease associated with the predicted motifs in H2 proteins.
#'
#' @author Bishoy Wadie, Evangelia Petsalaki
#' @return
#' A list containing all Drugs targeting H1 proteins from ChEMBL database and the ones that share indication with ClinVar disease.
#' @export
#' @importFrom utils read.delim read.table read.csv
#' @examples
#' new_hits = parse_qslim()
#' Clinvar_mapped = Map_ClinVar(new_hits_int = new_hits)
#' H1_chembl = H1_to_ChEMBL(new_hits_int = new_hits, clinvar_df = Clinvar_mapped)
H1_to_ChEMBL = function(slim_hits, clinvar_df ,zooma_output = NULL, use_zooma = F, host_viral = T, chemble_source = NULL, Drug_inds = NULL, indirect_mapping_dist = 2){
  new_hits_int = slim_hits
  # TODO change input path
  if (is.null(chemble_source)){
    chembl_data = read.csv(download_data(path = "Drug_mapping/filtered_chembl_map.csv.gz", is.rds = F), header = T)
  }
  else{
    chembl_data = chemble_source
  }

  if (is.null(Drug_indication)){
    Drug_indication = read.csv(download_data(path = "Drug_mapping/Drug_indication_chembl_27.csv", is.rds = F), header = T)
  }
  else{
    Drug_indication = Drug_inds
  }
  tmp = tempfile()
  utils::download.file("https://www.ebi.ac.uk/efo/efo.obo", destfile = tmp)
  efo_obo = ontologyIndex::get_ontology(tmp)

  if(host_viral){
    H1_v1 = new_hits_int %>% tidyr::separate(Dataset, into = c("V1","H1"), sep = "_")
    H1_v1 = H1_v1 %>% dplyr::select(V1, H1, uniprot)
  }
  else{
    H1_v1 = new_hits_int
    colnames(H1_v1)[which(names(H1_v1) == "Dataset")] = "H1"
    H1_v1 = H1_v1 %>% dplyr::select(H1, uniprot)
  }
  H1_v1 = H1_v1[!duplicated(H1_v1),]
  H1_entrez = AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = unique(H1_v1$H1), keytype = "UNIPROT", columns = "ENTREZID")
  H1_entrez = H1_entrez[!duplicated(H1_entrez$UNIPROT),]
  H1_v1 = dplyr::left_join(H1_v1, H1_entrez, by = c("H1" = "UNIPROT"))
  colnames(H1_v1)[which(names(H1_v1) == "ENTREZID")] = "H1_entrez"
  chembl_uniprot = chembl_data %>% dplyr::select(compound_chembl_id, accession)
  chembl_uniprot = chembl_uniprot[!duplicated(chembl_uniprot),]
  colnames(chembl_uniprot) = c("H1_Drug", "H1")
  H1_v1 = dplyr::left_join(H1_v1, chembl_uniprot, by = "H1")

  h1v1_molregno = chembl_data[which(chembl_data$compound_chembl_id %in% H1_v1$H1_Drug),c(1,2)]
  h1v1_molregno = h1v1_molregno[!duplicated(h1v1_molregno),]
  h1v1_molregno = dplyr::left_join(h1v1_molregno, Drug_indication)

  if (use_zooma){
    if(is.null(zooma_output)){
      # TODO change input path
      Zooma_mapped_EFO = read.csv(download_data(path = "Drug_mapping/Zooma_mapped_EFO_host_viral.csv", is.rds = F), header = T)
    }
    else{
      cleaned_phenotypes = send_to_zooma(clinvar_mapped = clinvar_df, output_all_disease_info = T)
      Zooma_mapped_EFO = zooma_output
      Zooma_mapped_EFO$ONTOLOGY.TERM.S. = stringr::str_replace(Zooma_mapped_EFO$ONTOLOGY.TERM.S., "_", ":")
      Zooma_mapped_EFO = dplyr::left_join(Zooma_mapped_EFO, cleaned_phenotypes[,c(5,2)], by = c("PROPERTY.VALUE" = "Disease"))
      Zooma_mapped_EFO = Zooma_mapped_EFO[!duplicated(Zooma_mapped_EFO),]
      colnames(Zooma_mapped_EFO)[5] = "EFO_id"
      colnames(Zooma_mapped_EFO)[9] = "PhenotypeList"
    }
    clinvar_EFO = dplyr::left_join(clinvar_df, Zooma_mapped_EFO[,c(9,5)])
  }
  else{
    clinvar_EFO = EFO_mapping(clinpath = clinvar_df, efo = efo_obo)
  }
  common_EFOs = intersect(clinvar_EFO[["clinvar_EFO"]]$EFO_id, h1v1_molregno$efo_id)
  overlap = h1v1_molregno[which(h1v1_molregno$efo_id %in% common_EFOs),]
  overlap = overlap[,c(2,4,5,8)]
  overlap = overlap[!duplicated(overlap),]
  colnames(overlap) = c("H1_Drug", "record_id", "max_phase_for_ind", "H1_Drug_EFO")
  H1_v1 = dplyr::left_join(H1_v1, overlap)

  if (host_viral){
    H1_v1$Dataset = paste0(H1_v1$V1,"_", H1_v1$H1)
  }
  else{
    H1_v1$Dataset = H1_v1$H1
  }
  x = clinvar_EFO[["clinvar_EFO"]] %>% dplyr::select(Dataset, uniprot, EFO_id, distance)
  x = x[!duplicated(x),]
  colnames(x)[3] = "Clinvar_EFO_id"
  H1_v1 = dplyr::left_join(H1_v1, x, by = c("Dataset", "uniprot"))

  Drug_H1_clinvar = Filter_H1_clinvar(df = H1_v1, efo = efo_obo)
  Drug_H1_clinvar = Drug_H1_clinvar[which(Drug_H1_clinvar$distance <= indirect_mapping_dist),]
  return(list("All_H1_mapped" = H1_v1, "H1_and_Clinvar" = Drug_H1_clinvar, "Mapped_phenoIDs" = clinvar_EFO[["Mapped_phenoIDs"]]))
}

cleanup_phenotypeONTs_clinvar = function(phenoIDs, phenoList){
  diseases = strsplit(phenoList, ";")[[1]]
  IDs = strsplit(phenoIDs, ";")[[1]]
  different = data.frame(phenoID = IDs, Disease = diseases)
  IDs_per_disease = list()
  counter = 0
  for (i in 1:nrow(different)){
    diff_Ids = strsplit(different$phenoID[i], ",")[[1]]

    for (j in 1:length(diff_Ids)){
      counter = counter + 1
      ont = strsplit(diff_Ids[j], ":")[[1]]
      ont_type = ont[1]
      ont_ID = paste(ont[-1], collapse = ":")
      result = data.frame(pheno_type = ont_type, pheno_ID = ont_ID, Disease = different$Disease[i])
      IDs_per_disease[[counter]] = result
    }
  }
  final = dplyr::bind_rows(IDs_per_disease)
  return(final)
}

#' Prepare Input to Zooma
#'
#' Send the output to Zooma as query
#'
#' @author Bishoy Wadie, Evangelia Petsalaki
#' @return
#' A character vector of Disease names
#' @export
send_to_zooma = function(clinvar_mapped, output_all_disease_info = F){
  clinvar_pheno = clinvar_mapped %>% dplyr::select(PhenotypeIDS, PhenotypeList)
  clinvar_pheno = clinvar_pheno[!duplicated(clinvar_pheno),]

  cleaned_phenotypes = list()
  for (i in 1:nrow(clinvar_pheno)){
    x = cleanup_phenotypeONTs_clinvar(clinvar_pheno$PhenotypeIDS[i], clinvar_pheno$PhenotypeList[i])
    y = data.frame(PhenotypeIDS = clinvar_pheno$PhenotypeIDS[i], phenotypeList = clinvar_pheno$PhenotypeList[i])
    cleaned_phenotypes[[i]] = cbind.data.frame(y,x)
  }
  cleaned_phenotypes = dplyr::bind_rows(cleaned_phenotypes)
  if (output_all_disease_info){
    return(cleaned_phenotypes)
  }
  else{
    return(unique(cleaned_phenotypes$Disease))
  }
}
Filter_H1_clinvar = function(df, efo){
  H1_v1_efos = df[!is.na(df$H1_Drug_EFO),]
  H1_v1_efos = H1_v1_efos %>% dplyr::select(H1_Drug_EFO, Clinvar_EFO_id)
  H1_v1_efos = H1_v1_efos[!is.na(H1_v1_efos$Clinvar_EFO_id),]
  H1_v1_efos = H1_v1_efos[!duplicated(H1_v1_efos),]

  selected = list()
  for (i in 1:nrow(H1_v1_efos)){
    if (H1_v1_efos[i,1] %in% efo[["children"]][[H1_v1_efos[i,2]]] ||
        H1_v1_efos[i,2] %in% efo[["children"]][[H1_v1_efos[i,1]]] ||
        H1_v1_efos[i,1] == H1_v1_efos[i,2]){
      selected[[i]] = H1_v1_efos[i,]
    }
    else{
      selected[[i]] = NULL
    }
  }
  selected = dplyr::bind_rows(selected)
  final = plyr::match_df(df, selected)

  efo_names = data.frame(efo$name)
  efo_names$Clinvar_EFO_id = rownames(efo_names)
  colnames(efo_names) = c("Clinvar_EFO_name", "Clinvar_EFO_id")
  rownames(efo_names) = NULL
  efo_names = efo_names[which(efo_names$Clinvar_EFO_id %in% final$Clinvar_EFO_id),]
  final = dplyr::left_join(final, efo_names)

  return(final)
}

EFO_mapping = function(clinpath, efo){
  phenoIDs = unique(clinpath$PhenotypeIDS)
  split_phenoIDs = list()
  for (i in 1:length(phenoIDs)){
    result = unlist(strsplit(phenoIDs[i], "[,|;]"))
    x = data.frame(PhenotypeIDS = phenoIDs[i], ontology = result)
    split_phenoIDs[[i]] = x
  }
  split_phenoIDs = dplyr::bind_rows(split_phenoIDs)
  split_phenoIDs = split_phenoIDs %>% tidyr::separate(ontology, into = c("ontology", "ID"), sep = ":", extra = "merge")
  split_phenoIDs = split_phenoIDs[!is.na(split_phenoIDs$ID),]

  for (i in 1:nrow(split_phenoIDs)){
    if (stringr::str_detect(split_phenoIDs$ID[i], ":")){
      split_phenoIDs$correct_id[i] = split_phenoIDs$ID[i]
    }
    else{
      split_phenoIDs$correct_id[i] = paste0(split_phenoIDs$ontology[i], ":", split_phenoIDs$ID[i])
    }
  }
  split_phenoIDs = split_phenoIDs[which(split_phenoIDs$ontology %nin% c("Gene", "MedGen")),]
  # EFO_onts_mapping is saved as internal data to package
  split_phenoIDs = dplyr::left_join(split_phenoIDs, EFO_onts_mapping[,c(1,3,7)], by = c("correct_id" = "mapped_curie"))

  x = split_phenoIDs[!is.na(split_phenoIDs$curie_id),]
  y = split_phenoIDs[is.na(split_phenoIDs$curie_id),]
  y[which(y$correct_id %in% efo$id),]$curie_id = y[which(y$correct_id %in% efo$id),]$correct_id
  y = y[!is.na(y$curie_id),]
  y$distance = 1

  final = rbind.data.frame(x, y)
  colnames(final)[5] = "EFO_id"

  clinvar_EFO = dplyr::left_join(clinpath, final[,c(1,5,6)])
  return(list("clinvar_EFO" = clinvar_EFO, "Mapped_phenoIDs" = final))
}
